import numpy as np
#from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import jenkspy

def load(filename):
    df = np.loadtxt(filename,usecols=(0,1,2))
    if df[0][0] > 2450000:
        df = df.T
        df[0] = df[0]-2450000
        df = df.T
    bad = (df.T[2] > 0)*1
    return np.vstack((df.T[:3],bad)).T
def getbad(data,limit):
    bad = True 
    while bad:
        bad = False
        effn = np.sum(data[3])
        if effn < 3:
            bad = False
            continue
        data = data.T[data[3]==1].T
        w_tot = np.sum(1/data[2]**2*data[3])
        dat_tot = np.sum(data[1]/data[2]**2*data[3])
        avg = dat_tot/w_tot
        chi2list = (data[1]-avg)**2/data[2]**2*data[3]
        ind = np.argmax(chi2list)
        chi2lim = np.sum(chi2list)/(effn-1)*limit**2
        if chi2list[ind] > chi2lim:
            bad = True
            data[3][ind] = 0
    return data



def bin(data,fac=1,limit=3):
    nclass = len(set([ int(i) for i in data[0]]))/fac
    breaks = jenkspy.jenks_breaks(data[0], nb_class=int(nclass))
    # print('days=',len(set([ int(i) for i in data[0]])),'bins=',int(nclass))
    output = []
    for i in range(len(breaks)-1):
        lower = breaks[i]-(i==0)
        upper = breaks[i+1]
        index = np.logical_and(data[0]>lower,data[0]<=upper)
#        lower = breaks[i]
#        upper = breaks[i+1]
#        index = np.logical_and(data[0]>=lower,data[0]<upper)
        tmp_data = data.T[index].T
        tmp_data = getbad(tmp_data,limit=limit)
        w_tot = np.sum(1/(tmp_data[2]**2)*tmp_data[3])
        mag_tot = np.sum(tmp_data[1]/(tmp_data[2]**2)*tmp_data[3])
        date_tot = np.sum(tmp_data[0]/(tmp_data[2]**2)*tmp_data[3])
        mag_avg = mag_tot/w_tot
        date_avg = date_tot/w_tot
        err = w_tot**(-0.5)
        output.append([date_avg,mag_avg,err,1])
    return np.array(output).T

def mag2flux(data):
    _data = data.copy()
    refmag = 15
    mag = data[1]
    err = data[2]
    flux = 10**(0.4*(refmag-mag))
    ferr = flux*err*0.4*np.log(10)
    _data[1] = flux
    _data[2] = ferr
    return _data

def flux2mag(data,add=0):
    _data = data.copy()
    refmag = 15
    flux = data[1]+add
    ferr = data[2]
    mag = -2.5*np.log10(flux)+refmag
    magerr = ferr/flux*2.5/np.log(10)
    _data[1] = mag
    _data[2] = magerr
    return _data

def plot_bin(data,fac=1):
    # plt.scatter(data[0],data[1],facecolors='none',marker='s',edgecolors='blue')
    plt.errorbar(data[0],data[1],yerr=data[2], fmt='o', markersize=8,
                      fillstyle='none', markeredgewidth=1.5)
    n,ap = file.split('.')
    # data = mag2flux(data)
    bin_data = bin(data,fac=fac)
    # bin_data = flux2mag(bin_data,add=0)
    np.savetxt(f'{n}_bin2.{ap}',bin_data[:3].T,fmt="%12.6f %12.6f %12.6f")
    plt.errorbar(bin_data[0],bin_data[1],yerr=bin_data[2],fmt='s')
    # plt.gca().invert_yaxis()
    return bin_data


file = 'PEST_ip_sub.dia'
data = load(file).T
# data = data.T[data[0] > 9608].T
# data = plot_bin(data,fac=1)
data = data.T[data[0] < 9608].T
data = plot_bin(data,fac=0.3)

plt.show()
