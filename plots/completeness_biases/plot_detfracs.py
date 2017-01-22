from __future__ import division
import sys
sys.path.insert(0,'/Users/yotam/research/SAM-CANDELS_pipeline/plots/')
from includes import *
init_plotting_custom()
from numpy.random import random as rand

def montecarlo(sdf,param,bins,N):
    allcs,bs = np.histogram(sdf[param],bins=bins)
    Ncs = []
    Nsdf = sdf.shape[0]
    dps = sdf.dp.values
    for i in range(N):
        # assign new random numbers and do MC
        rands = rand(Nsdf)
        det = np.full(Nsdf,np.nan)
        det[rands <= dps] = 1
        det[~(rands <= dps)] = 0
        """
        sdf['rand'] = rand(sdf.shape[0])
        sdf['detected?'] == np.nan
        sdf['detected?'][sdf.rand <= sdf.dp] = 1
        sdf['detected?'][~(sdf.rand <= sdf.dp)] = 0
        sdf['detected?'] = sdf['detected?'].astype(int)
        """
        # histogram and get counts
        ###cs,bs = np.histogram(sdf[param][sdf['detected?']==1],bins=bins)
        cs,bs = np.histogram(sdf[param][det==1],bins=bins)
        Ncs.append(cs)
    # debug
    Ncs = np.array(Ncs)
    print np.mean(Ncs,axis=0)
    print np.std(Ncs/allcs,axis=0)
    # get average and stds
    Ncs = Ncs/allcs
    #for item in Ncs: print item
    mean_cs = np.mean(Ncs,axis=0)
    stds = np.std(Ncs,axis=0)
    return mean_cs,bs,stds,allcs

df = read_hdf('/Volumes/G-RAID/yotam/SAM-CANDELS_pipeline/catalogs/SC/12.16/goodsn.h5')
fieldArea = (32.*32.)*60.*60. # arcsec^2

#=== plot ========================================================
# redshift bins
#zbins = [0.0,0.5,1.0,2.0,3.0]
#zbins = [0.1,0.5,1.2,1.7,2.2]
zbins = [0.1,0.5,1.24,1.72,2.15,2.57,3.0]
df['zbin'] = np.digitize(df.redshift,zbins)

# magnitude cuts
mcuts = [26,25.5,25,24.5]
#mcuts = [25.5,25,24.5]
H = df['wfc3f160w_dust_mag']

# which parameter
df['U-V'] = df['U_rest'] - df['V_rest']
df['lreff'] = np.log10(df['reff'])
df['lmstar'] = np.log10(df['mstar']) + 10
param = 'lmstar'
paramlabel = r'$\log{[M_{\star}/M_{\odot}]}$'
#param = 'lreff'
#paramlabel = r'$\log{[R_{\rm eff}/{\rm kpc}]}$'
#param = 'n_sersic'
#paramlabel = 'n_sersic'
#param = 'U-V'
#paramlabel = 'U-V'
nbins = 15
bins = np.linspace(np.min(df[param]),np.max(df[param]),nbins)
NMC = 100

# setup figure
scale = 1.25
fig,axes = plt.subplots(len(zbins)-1,len(mcuts),
            figsize=(6*scale,6*scale),
            sharey=False,sharex=True)
ax = axes.flatten()

# loop over redshift bins
for i,zbin in enumerate(np.unique(df.zbin)[1:-1]):
    inzbin = (df.zbin == zbin)
    # loop over magnitude cuts
    for j,mcut in enumerate(mcuts):
        # DO n_s < 2.5
        gdf = df[(H < mcut)&(inzbin)]
        gdf = gdf[gdf['n_sersic'] < 2.5]
        cs,bs,es,allcs = montecarlo(gdf[(H < mcut)&(inzbin)],param,bins,NMC)
        # plot
        xs = (bs[1:] + bs[:-1])/2.
        ys = cs
        yerrs = es
        # get rid of bins with only one object
        thresh = 10
        goodbins = np.where(allcs > thresh)[0]
        xs,ys,yerrs = xs[goodbins],ys[goodbins],yerrs[goodbins]
        axes[i][j].errorbar(xs,ys,yerr=yerrs,c='b',
                            marker='o',ms=2,
                            ls='None')
        # DO n_s < 2.5
        gdf = df[(H < mcut)&(inzbin)]
        gdf = gdf[gdf['n_sersic'] > 2.5]
        cs,bs,es,allcs = montecarlo(gdf[(H < mcut)&(inzbin)],param,bins,NMC)
        # plot
        xs = (bs[1:] + bs[:-1])/2.
        ys = cs
        yerrs = es
        # get rid of bins with only one object
        thresh = 10
        goodbins = np.where(allcs > thresh)[0]
        xs,ys,yerrs = xs[goodbins],ys[goodbins],yerrs[goodbins]
        axes[i][j].errorbar(xs,ys,yerr=yerrs,c='r',
                            marker='o',ms=2,
                            ls='None')

        print 'done with %d,%d' %(i,j)

# formatting
for i in range(len(mcuts)):
    axes[0][i].set_title(r'$\mathrm{H160W} < \, %s$' %mcuts[i])
for i in range(len(zbins)-1):
    zstring = r'$%s < z < %s$' %(zbins[i],zbins[i+1])
    axes[i][len(mcuts)-1].text(1.05,0.5,zstring,rotation=-90,
                ha='center',va='center',
                transform=axes[i][len(mcuts)-1].transAxes)

ax[20].set_xlabel(paramlabel)
#ax[20].set_ylabel(r'$N_{\rm detected}/N_{\rm true}$')
ax[20].set_ylabel('detected fraction')

for i in range(len(ax)):
    ax[i].set_ylim(ymax = 1.)
for i in range(len(ax)):
    ticks = ax[i].get_yticks()
    if len(ticks) > 6:
        ax[i].set_yticks(ticks[::-2])


fig.subplots_adjust(wspace=0.1,hspace=0.15)

# save
plt.tight_layout()
fig.subplots_adjust(right=0.95)
#plt.savefig('testing_detfrac_%s.png' %param,dpi=300)
plt.savefig('testing_detfrac_%s.pdf' %param)

#+======================================================================


