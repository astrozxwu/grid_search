import numpy as np
from scipy.optimize import minimize
import os
from config import *
from utils import *
from astropy.table import Table

init = np.loadtxt(os.path.join(sPATH,'init.dat'))
def chi2(theta,param_fix):
    _t0,_u0,_tE,_tstar = theta
    param = {'t0':_t0,'u0':_u0,'tE':_tE,'tstar':_tstar} 
    param = mcmc2VBBL_stellar({**param, **param_fix})
    _chi2 = 0
    for i,_a,_bl in zip(flux,a1,bl):
        model = Binary_model(i[0], param, plx, a1=_a)
        _chi2 += getchi2_single(i, model, blending=_bl)[0]
    return _chi2 
def grid(param_fix):
    res = minimize(chi2, x0=init, method='nelder-mead', args=(param_fix), bounds=[[9400,9800],[-2,2],[10,500],[0.5,30]], options={'xatol':5e-3,'fatol':2,'disp':False})
    _t0,_u0,_tE,_tstar = res.x
    param = {'chi2':res.fun,'covergence':1*res.success,'t0':_t0,'u0':_u0,'tE':_tE,'tstar':_tstar} 
    return {**param,**param_fix}
# test consistency
# {'chi2': 970.799908, 't0': 9558.128689355906, 'u0': -0.48598393424984426, 'teff': 27.989189881928354, 'rho': 0.063489, 'log_s': 0.4888744596374558, 'log_q': -0.30720131713610244, 'alpha': 3.754068, 'tE': 57.59282953485272, 's': 3.0822968300285827, 'q': 0.4929452462967452}

grids = np.loadtxt('./grids')
all_res = Table(names=['chi2','convergence','t0','u0','tE','tstar','logs','logq','alpha'],dtype=['f4','i4','f4','f4','f4','f4','f4','f4','f4'])
formatter = {'chi2':'.3f','t0':'.4f','u0':'.4f','tE':'.4f','tstar':'.4f','logs':'.5f','logq':'.5f','alpha':'.5f'}
def format_table(tab,formatter):
    for i in formatter.keys():
        tab[i].info.format = formatter[i]
format_table(all_res,formatter)
for i in range(len(grids[:])):
    param_fix = {'logs':grids[i][0],'logq':grids[i][1],'alpha':grids[i][2]}
    res = grid(param_fix)
    all_res.add_row(res.values())
    if i%5 == 4:
        all_res.write(os.path.join(oPATH,'best.dat'),overwrite=True,format='ascii.tab')
all_res.write(os.path.join(oPATH,'best.dat'),overwrite=True,format='ascii.tab')
