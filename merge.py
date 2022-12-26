from astropy.table import Table, vstack

n_split = 20


tab = Table(names=['chi2','convergence','t0','u0','teff','rho','logs','logq','alpha'],dtype=['f4','i4','f4','f4','f4','f4','f4','f4','f4'])
formatter = {'chi2':'.3f','t0':'.4f','u0':'.4f','teff':'.4f','rho':'.4f','logs':'.5f','logq':'.5f','alpha':'.5f'}
def format_table(tab,formatter):
    for i in formatter.keys():
        tab[i].info.format = formatter[i]
format_table(tab,formatter)

for i in range(n_split):
    f = f't_{i}/output/best.dat'
    df = Table.read(f,format='ascii.tab')
    tab = vstack([tab,df])

tab.sort('chi2')
print(tab)


print(tab[tab['convergence'] == 0])
tab.write('all.tab',overwrite=True,format='ascii.tab')
