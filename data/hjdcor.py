from utils import *


ra, dec = (15 + 1./60 + 00.65/3600)*15, -(54 + 24./60 + 00.5/3600)
fn = 'PEST_gp_sub_bin.dia'
df = load(fn)
df[0] = hjdcor(ra,dec,df[0])
n,ap = fn.split('.')


np.savetxt(f'{n}.hjdcor.{ap}',df.T,fmt='%12.6f  %12.6f  %12.6f')
