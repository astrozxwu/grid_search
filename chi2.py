import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MaxNLocator
nparm = 7
parmlist = ['t0','u0','tE','tstar','log_s','log_q','alpha']
#parmlist = [parmlist[i] for i in [3,4,5]]
#   9397.2606521    0.1426160     13.4144779     -0.500000000      0.600000000      0.802444444

np.set_printoptions(suppress = True)

def draw(sample,delta=1,nsigma=4,clip=True):
    chi2 = min(sample[0])
    if clip : sample = sample.T[sample[0] < 2000].T#chi2+((nsigma+3)**2)*delta].T
    sample = np.array(sorted(sample.T,key = lambda x:-x[0])).T
    print(chi2)
    # data in [parm1,parm2....,parmn,chi2]
    '''
    plot in 
    0
    1
    2
    3
    4
    5
     0 1 2 3 4 5 
    '''
    '''
    x:down
    y:right
    '''
    import matplotlib.cm as cm
    import matplotlib as mpl
    cmap = mpl.cm.viridis_r
    cmap = mpl.colors.ListedColormap(['red','gold','limegreen','blue','gray','lightblue'][:nsigma+1])
    bounds = delta*np.arange(nsigma+1)**2
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(11,10))#,dpi=400)
    gs = fig.add_gridspec(nparm, nparm, hspace=0, wspace=0)
    try:
        axs = gs.subplots(sharex=False, sharey=False)
    except:
        axs = [[plt.subplot(gs[i,j]) for j in range(nparm)] for i in range(nparm)]
    #fig, axs = plt.subplots(nparm, nparm)


    for i in range(nparm):
        for j in range(nparm):
            if not i == nparm-1:
                axs[i][j].get_xaxis().set_visible(False)
            if not j == 0:
                axs[i][j].get_yaxis().set_visible(False)
    for x in range(nparm):
        for y in range(nparm):
            if y > x :
                axs[x][y].set_visible(False)
                continue
            if y == x :
                _,bins,_ = axs[x][x].hist(sample[x+1],bins=50,density=True)
                #mu, sigma = scipy.stats.norm.fit(sample[x+1])
                #best_fit_line = scipy.stats.norm.pdf(bins, mu, sigma)
                #axs[x][x].plot(bins, best_fit_line,c='black')
                continue
            #(x=1,y=0)
            axs[x][y].scatter(sample[y+1],sample[x+1],s=1,c=norm(sample[0]-chi2),cmap=cmap)#,alpha=0.5)
            #bin_x = sample[y+1]
            #_y = sample[x+1]
            #_z = sample[0]-chi2
            #axs[x][y].hexbin(_x,_y,C=_z,norm=norm,gridsize=50,cmap=cmap,bins=None)
    for x in range(nparm):
            axs[x][0].set_ylabel(parmlist[x],fontsize=12)
            axs[nparm-1][x].set_xlabel(parmlist[x],fontsize=12)
            axs[x][0].tick_params(axis='both', which='major', labelsize=9)
            axs[nparm-1][x].tick_params(axis='both', which='major', labelsize=9)
            #axs[nparm-1][x].xaxis.set_ticks(np.linspace(start, end, 5))
            axs[x][0].yaxis.set_major_locator(MaxNLocator(nbins=4))
            axs[nparm-1][x].xaxis.set_major_locator(MaxNLocator(nbins=4))
            axs[nparm-1][x].set_xticks(axs[nparm-1][x].get_xticks()[1:-1])
            axs[x][0].set_yticks(axs[x][0].get_yticks()[1:-1])

            plt.setp( axs[nparm-1][x].xaxis.get_majorticklabels(), rotation=45)
            plt.setp( axs[x][0].yaxis.get_majorticklabels(), rotation=45)
#            ax.tick_params(axis='both', which='minor', labelsize=8)

    fig.tight_layout()
    fig.subplots_adjust(right=0.9)
    l = 0.92
    b = 0.09
    w = 0.03
    h = 1.07 - 2*b
    rect = [l,b,w,h]
    cbar_ax = fig.add_axes(rect)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),cax = cbar_ax)
    cbar.ax.tick_params(labelsize=14)
    plt.suptitle(r'min $\chi^2$: ' + f'{chi2:.2f}')
    plt.show()



sample = np.loadtxt('all.tab').T[[0,2,3,4,5,6,7,8]]
sample[3] = sample[3]/sample[2]
sample[4] = sample[3]*sample[4]
#sample[0] = sample[0]*(-2)
draw(sample,delta=20,nsigma=4,clip=True)
