import os
from utils import *
dPATH = '../data/'
oPATH = './output/'
sPATH = './startups/'

filelist = ['LCOGT-g_bin.hjdcor.pho', 'LCOGT-r_bin.hjdcor.pho', 'LCOGT-i_bin.hjdcor.pho', 'AO_g_bin.hjdcor.pho', 'AO_r_bin.hjdcor.pho',
            'AO_i_bin.hjdcor.pho', 'PEST_gp_sub_bin.hjdcor.dia', 'PEST_rp_sub_bin.hjdcor.dia', 'PEST_ip_sub_bin.hjdcor.dia', 'bE_ap_g_bin.pho', 'bi_ap_g_bin.pho', 'bm_ap_g_bin.pho']
VBBL.SetObjectCoordinates(os.path.join(dPATH,'coord.txt'),'.')
labels = [i.split('.')[0].replace('_', ' ').replace('-', ' ').replace('bin',' ')
          for i in filelist]
clist = ['r', 'black', 'darkgreen', 'blue', 'orange', 'magenta',
         'lime', 'olivedrab', 'darkslategray', 'cyan', 'gold', 'darkviolet']
parameters = ['t0', 'u0', 'tE', 'tstar',
              'log_s', 'log_q', 'alpha', 'pi1', 'pi2']
n_params = len(parameters)
headers = ['chi2'] + parameters
nobs = len(filelist)
#errnorm = [1, 0.7, 0.6, 1., 0.5, 0.8, 4., 9., 7., 1.6, 1.5, 1.3]
errnorm = [0.9, 0.6, 0.6, 1., 0.5, 0.55, 2.4, 3.0, 1.8, 1.2, 1.4, 1.1]
erradd = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
bl = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
ra, dec = (15 + 1./60 + 00.65/3600)*15, -(54 + 24./60 + 00.5/3600)
gamma_lld = [0.87, 0.76, 0.64, 0.87, 0.76,
             0.64, 0.87, 0.76, 0.64, 0.87, 0.87, 0.87]
a1 = [3*u/(2+u) for u in gamma_lld]  # 3*gamma_lld/(2+gamma_lld)
plx = False 
if not plx:
    headers = headers[:-2]

data = [(load(dPATH+i)) for i in filelist]
nobs = len(filelist)
for i in range(nobs):
    data[i][2] = (data[i][2]**2 + erradd[i]**2)**0.5
    data[i][2] = data[i][2]*errnorm[i]
    if i in [6, 7, 8]:
        data[i][2] = data[i][2]/10000
        data[i][1] = data[i][1]/10000
flux = []
for i in range(nobs):
    if i in [6, 7, 8]:
        flux.append(data[i])
    else:
        flux.append(mag2flux(data[i]))
