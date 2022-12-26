import numpy as np
import os
import subprocess

logs = np.linspace(0,0.8,17)
logq = np.linspace(-2,2,21)
alpha = np.linspace(0,2*np.pi,31)
ntot = int(len(logs)*len(logq)*len(alpha))
n_split = 22
n_single = int(ntot/n_split)

grids = []
for i in alpha:
    for j in logq:
        for k in logs:
            grids.append([k,j,i])

for i in range(n_split):
    os.system('mkdir t_'+str(i))
    os.system('cp -r ./source_code/* '+'./t_'+str(i)+'/')
    iend = (i+1)*n_single if i<n_split-1 else ntot
    np.savetxt('./t_'+str(i)+'/grids',grids[i*n_single:iend],fmt='%10.5f')

my_env = os.environ.copy()
my_env["OMP_NUM_THREADS"] = "1"

for i in range(n_split):
    #os.chdir('t_'+str(i))
    #os.system('python3 run.py &')
    subprocess.Popen(['python3','run.py'],cwd='t_'+str(i),shell=False,env=my_env)
