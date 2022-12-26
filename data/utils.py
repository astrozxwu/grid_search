import numpy as np
import VBBinaryLensing
import scipy.optimize as op
from astropy import time, coordinates as coord, units as u
VBBL = VBBinaryLensing.VBBinaryLensing()
VBBL.RelTol = 1e-03

def load(fn,subs=True):
    df = np.loadtxt(fn,usecols=(0,1,2)).T
    if subs:
        if df[0][0] > 2450000:
            df[0] = df[0] - 2450000
    return df

def hjdcor(ra,dec,t,fmt='jd'):
    cod = coord.SkyCoord(ra,dec,unit=(u.deg,u.deg),frame='icrs')
    offset = 0 if t[0] > 2450000 else 2450000
    ts = time.Time(t+offset,format=fmt,scale='utc')

    greenwich = coord.EarthLocation.of_site('greenwich')

    ltt_helio = ts.light_travel_time(cod,'heliocentric',location=greenwich)
    return (t + ltt_helio).value


def mag2flux(data):
    refmag = 18
    mag = data[1]
    err = data[2]
    flux = 3631*10**(0.4*(refmag-mag))
    ferr = flux*err*0.4*np.log(10)
    return np.array([data[0],flux,ferr])

def flux2mag(data):
    refmag = 18
    if len(data) == 1:
        return -2.5*np.log10(data[0]/3631)+refmag
    flux = data[1]
    ferr = data[2]
    mag = -2.5*np.log10(flux/3631)+refmag
    magerr = ferr/flux*2.5/np.log(10)
    return np.array([data[0],mag,magerr])

def format_param(param):
    # format to t0,u0,tE,s,q,alpha,rho,pi1,pi2
    keys = param.keys()
    if not 'tE' in keys:
        param['tE'] = param['teff'] / abs(param['u0'])
    if (not 's' in keys) or (not 'q' in keys):
        param['s'] = 10 ** param['logs']
        param['q'] = 10 ** param['logq']
    if not 'rho' in keys:
        param['rho'] = 10 ** param['log_rho']

    return param

def sample_gen(par1, par2, n, method='gauss'):
    if method == 'gauss':
        return par1+par2*np.random.randn(n)
    if method == 'uniform':
        return np.random.uniform(low=par1, high=par2, size=n)
    if method == 'log-uniform':
        return np.exp(np.random.uniform(low=np.log(par1), high=np.log(par2), size=n))

def Binary_model(_t, param, plx=False, a1=0):
    '''
    using dict param : 't0','u0','tE','s','q','alpha','rho'
    '''
    _t0,_u0,_tE = [param[i] for i in ['t0','u0','tE']]
    _s,_q,_alpha,_rho = [param[i] for i in ['s','q','alpha','rho']]
    if a1 : VBBL.a1 = a1
    y1 = np.zeros(len(_t))
    y2 = np.zeros(len(_t))
    params = [np.log(_s), -np.log(_q), _u0, 2*np.pi - _alpha, np.log(_rho), np.log(_tE), _t0]
    if plx:
        pi1,pi2 = [param[i] for i in ['pi1','pi2']]
        params.append(pi1)
        params.append(pi2)
        return np.array(VBBL.BinaryLightCurveParallax(params, _t, y1, y2))
    return np.array(VBBL.BinaryLightCurve(params, _t, y1, y2))

def getchi2_single(dat, model, blending=False, ap=False):
    '''
    dat like [fluxes,errors]
    model like [magnifications]
    '''
    y = dat[1]/dat[2]
    if ap and blending:
        A = np.vstack([model/dat[2], 1/dat[2], dat[3]/dat[2]]).T
    elif blending:
        A = np.vstack([model/dat[2], 1/dat[2]]).T
    elif ap:
        A = np.vstack([model/dat[2],dat[3]/dat[2]]).T
    else:
        A = np.vstack([model/dat[2]]).T
    res = np.linalg.lstsq(A, y, rcond=None)

    f = np.append(res[0],[0,0])
    if ap and not blending:
        f[2] = f[1]
        f[1] = 0
    return res[1][0], f[:3]

def usual2VBBL(param):
    s = 10**param[5]
    q = 10**param[6]
    a = param[7]
    t0 = param[1]
    u0 = param[2]
    tE = np.abs(param[3]/param[2])
    u0_ = u0+q/(1+q)*(s-1/s)*np.sin(a)
    teff_ = tE*u0_
    t0_ = t0+q/(1+q)*(s-1/s)*np.cos(a)*tE
    param[2] = u0_
    param[3] = np.abs(teff_)
    param[1] = t0_
    return param

def VBBL2mcmc_planet(param):
    '''
    INPUT : t0 u0 tE rho s q alpha piEN piEE
    OUTPUT : new_t0 new_u0/w new_teff t_star logw logq alpha piEN piEE
    '''
    _param = {i:param[i] for i in ['alpha']}
    t0 = param['t0']
    u0 = param['u0']
    tE = param['tE']
    rho = param['rho']
    s = param['s']
    q = param['q']
    a = param['alpha']
    qf = qfac(s)
    w = qf * q/(1+q)**2
    t_star = rho*tE
    if s > 1:
        new_u0 = u0+q/(1+q)*(s-1/s)*np.sin(a)
        new_t0 = t0-q/(1+q)*(s-1/s)*np.cos(a)*tE
        new_teff = abs(tE*new_u0)
    else:
        new_u0 = u0
        new_t0 = t0
        new_teff = abs(tE*new_u0)
    _param['t0'] = new_t0
    _param['u0/w'] = new_u0/w
    _param['teff'] = new_teff
    _param['t_star'] = t_star
    _param['logw'] = np.log10(w)
    _param['logq'] = np.log10(q)
    _param['flag_wide'] = 1*(param['s']>1)
    return _param

def mcmc2VBBL_planet(param):
    '''
    INPUT : new_t0 new_u0/w new_teff t_star logw logq alpha piEN piEE
    OUTPUT: t0 u0 teff rhos logs logq alpha piEN piEE
    '''
    _param = {i:param[i] for i in ['alpha']}
    t0 = param['t0']
    u0_w = param['u0/w']
    teff = param['teff']
    t_star = param['t_star']
    w = 10**param['logw']
    q = 10**param['logq']
    a = param['alpha']
    qf = w/q *(1+q)**2
    s = invqfac(qf)
    u0 = u0_w*w
    tE = abs(teff/u0)
    rho = t_star/tE
    if param['flag_wide']:
        s = 1/s
        u0 -= q/(1+q)*(s-1/s)*np.sin(a)
        t0 += q/(1+q)*(s-1/s)*np.cos(a)*tE
    _param['t0'] = t0
    _param['u0'] = u0
    _param['tE'] = tE
    _param['rho'] = rho
    _param['s'] = s
    _param['q'] = q
    return _param

def qfac(s):
    st = (s+1/s)
    cphi = 3/4*st*(1-(1-32/9/st**2)**0.5)
    sphi = (1-cphi**2)**0.5
    qfac = 4*sphi**3/(st-2*cphi)**2
    return qfac

def invqfac(q):
    def f(x):
        return qfac(x)-q
    return op.brentq(f,1e-5,1)

def VBBL2mcmc_stellar(param):
    '''
    INPUT : t0 u0 tE rho s q alpha piEN piEE
    OUTPUT : new_t0 new_u0/w new_teff t_star logw logq alpha piEN piEE
    '''
    _param = {i:param[i] for i in ['alpha','rho']}
    t0 = param['t0']
    u0 = param['u0']
    tE = param['tE']
    s = param['s']
    q = param['q']
    a = param['alpha']
    if s > 1:
        u0 += q/(1+q)*(s-1/s)*np.sin(a)
        t0 -= q/(1+q)*(s-1/s)*np.cos(a)*tE
    teff = abs(tE*u0)
    _param['t0'] = t0
    _param['u0'] = u0
    _param['teff'] = teff
    _param['logs'] = np.log10(s)
    _param['logq'] = np.log10(q)
    return _param

def mcmc2VBBL_stellar(param):
    '''
    INPUT : new_t0 new_u0/w new_teff t_star logw logq alpha piEN piEE
    OUTPUT: t0 u0 teff rhos logs logq alpha piEN piEE
    '''
    _param = {i:param[i] for i in ['alpha','rho']}
    t0 = param['t0']
    u0 = param['u0']
    teff = param['teff']
    s = 10**param['logs']
    q = 10**param['logq']
    a = param['alpha']
    tE = abs(teff/u0)
    if s > 1:
        u0 -= q/(1+q)*(s-1/s)*np.sin(a)
        t0 += q/(1+q)*(s-1/s)*np.cos(a)*tE
    _param['t0'] = t0
    _param['u0'] = u0
    _param['tE'] = tE
    _param['s'] = s
    _param['q'] = q
    return _param

def c2w(param):
    qc = param['q']
    sc = param['s']
    qw = qc*(1-qc)**-2
    sw = (sc**-1)*(1+qc)*(qc**2-qc+1)**-0.5
    _param = param.copy()
    _param['s'] = sw
    _param['q'] = qw
    return _param

if __name__ == '__main__':
    param = {'s':2.1,'q':0.4,'alpha':1,'t0':9600,'tE':100,'u0':-0.05,'rho':0.15}
    print(param)
    print(VBBL2mcmc_stellar(param))
    print((mcmc2VBBL_stellar(VBBL2mcmc_stellar(param))))
    import matplotlib.pyplot as plt

    logs = np.linspace(0.05,2,1000)
    logw = np.log10(qfac(10**logs))
    plt.scatter(logs,logw)
    plt.show()
