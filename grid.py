import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.table import Table,vstack
import pandas as pd

def plot(ax,data,cols,seeds=[],truths=[]):
    ax.scatter(data[cols[0]],data[cols[1]],marker='s',s=20,c=norm(data[0]-m),cmap=cmap)
    # se = data.T[data[0] < 1500].T
    # ax.scatter(se[cols[0]],se[cols[1]],marker='o',s=60,facecolors='none',edgecolors='r')
    for (i,j) in zip(seeds,truths):
        c=norm(i[0]-m)
        c=cmap(c)
        ax.scatter(i[cols[0]],i[cols[1]],marker='o',s=80,facecolors='none',edgecolors=c,lw=1)
        c=norm(j[0]-m)
        c=cmap(c)
        ax.scatter(j[cols[0]],j[cols[1]],marker='s',s=80,facecolors='none',edgecolors=c,lw=2)

        ax.arrow(i[cols[0]],i[cols[1]],j[cols[0]]-i[cols[0]],j[cols[1]]-i[cols[1]],fc='k', ec='k',lw=1)
    # print('min chi2:',m)

flist = ['all.tab']
df = Table.read(flist[0],format='ascii.tab')



m = min(df['chi2'])
print('min chi2:',m)
columns = df.colnames
df = df.to_pandas()
#df = df[df.alpha < 5.6]
#df = df[df.alpha > 5.0]


_s_q = df.sort_values('chi2').groupby(by=['logs','logq'],as_index=False).first()
_s_q = _s_q[['chi2','logs','logq','alpha']].to_numpy().T
_s_alpha = df.sort_values('chi2').groupby(by=['logs','alpha'],as_index=False).first()
_s_alpha = _s_alpha[['chi2','logs','logq','alpha']].to_numpy().T
_q_alpha = df.sort_values('chi2').groupby(by=['logq','alpha'],as_index=False).first()
_q_alpha = _q_alpha[['chi2','logs','logq','alpha']].to_numpy().T

cmap = mpl.cm.Set1
cmap = mpl.colors.ListedColormap(['red','gold','limegreen','blue','gray'])
delta = 30

bounds = [0, delta, delta*2, delta*3, delta*4, delta*5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
fig = plt.figure(figsize=(11,10))#,dpi=400)
gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
axs = gs.subplots(sharex=False, sharey=False)
parmlist = ['logs','logq','alpha']
axs[0][1].set_visible(False)

ax1 = axs[0][0]
ax2 = axs[1][0]
ax3 = axs[1][1]

# test = (df.sort_values('chi2').groupby(by=['logs','logq'],as_index=False).first()[['logs','logq','chi2']])
# test['chi2'] = test['chi2'] - m
# print(test)
# test = df.reset_index().pivot('logs','logq','chi2')
# import seaborn as sns
# ax1 = sns.heatmap(test,cmap=cmap)


ax2.get_shared_x_axes().join(ax1, ax2)
ax2.get_shared_y_axes().join(ax2, ax3)
# ax1.set_ylim(ax3.get_xlim())
# ax1.set_yticks(ax3.get_xticks())
import numpy as np

def pp(ax,keys=['logs','logq']):
    mesh = np.meshgrid(np.unique(df[keys[0]]),np.unique(df[keys[1]]))
    grids = np.array(mesh).T.reshape(-1,2)
    grids = np.concatenate([grids.T,[len(grids)*[np.nan]]]).T
    x = pd.DataFrame(grids,columns=[keys[0],keys[1],'chi2']).groupby(by=keys).first()
    _s_q_ = df.sort_values('chi2').groupby(by=keys).first()
    x.update(_s_q_)
    x = x.reset_index()
    x = x.pivot(index=keys[1], columns=keys[0], values='chi2')
    ax.pcolormesh(mesh[0],mesh[1],norm(x.values-m),cmap=cmap)


pp(ax1,['logs','logq'])
pp(ax2,['logs','alpha'])
pp(ax3,['logq','alpha'])
# plot(ax1,_s_q,cols=[1,2])
# plot(ax2,_s_alpha,cols=[1,3])
# plot(ax3,_q_alpha,cols=[2,3])
ax1.get_xaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)

plt.subplots_adjust(wspace=0, hspace=0)
ax1.set_ylabel('log_q')
ax2.set_ylabel('alpha')
ax2.set_xlabel('log_s')
ax3.set_xlabel('log_q')
cbar_ax = fig.add_axes([0.9,0.11,0.03,0.77])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),cax = cbar_ax)
#plt.savefig(f'alpha{prefix}_chi2_{m:.0f}.pdf')
plt.suptitle(r'min $\chi^2$: ' + f'{m:.2f}')
plt.show()

