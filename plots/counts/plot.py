import sys
sys.path.insert(0,'/Users/yotam/research/SAM-CANDELS_pipeline/plots/')
from includes import *
from numpy.random import random as rand
init_plotting_custom()

# define parameters
field = 'goodsn'
candels_fieldsize = 147.
mocks_fieldsize = 32.*32.

# grab catalogs
cdf = read_hdf('/Users/yotam/research/mocks/pipeline/lib/CANDELS/%s.h5' %field)
mdf = read_hdf('/Volumes/G-RAID/yotam/catalogs/SC/12.16/%s.h5' %field)

# selection on CANDELS data
cdf = cdf[cdf.PhotFlag == 0]
cdf = cdf[cdf.CLASS_STAR < 0.8]
# fluxes to magnitudes for CANDELS
for col in cdf.columns:
    if col.endswith('FLUX'):
        magcol = col.rstrip('FLUX')+'MAG'
        cdf[magcol] = -2.5*np.log10(cdf[col]) + 23.9


#-------------------------
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
    print np.std(Ncs,axis=0)
    mean_cs = np.mean(Ncs,axis=0)
    stds = np.std(Ncs,axis=0)
    return mean_cs,bs,stds,allcs

def returnxy(df,band,bins):
    # digitize by bins
    mbin = np.digitize(df[band],bins)
    # histogram
    counts,bins = np.histogram(df[band],bins)
    bincenters = 0.5*(bins[1:] + bins[:-1])
    # return 
    x = bincenters
    y = counts
    return x,y

# redshift bins
zbins = [(0.1,0.5),(0.5,1.24),(1.24,1.72),
         (1.72,2.15),(2.15,2.57),(2.57,3.0),
         (0.0,2.0),(2.0,4.0),(0.0,4.0)]


# set up figure
fig,axes = plt.subplots(3,3,figsize=(6,6),sharex=True)
ax = axes.flatten()

# set magnitude bins
dbin = 1.0
bins = np.arange(18,31,dbin)

# loop over redshifts
mocksz = 'redshift'
#mocksz = 'z_est'
candelsz = 'zbest'
mdf_all = mdf
cdf_all = cdf
del mdf,cdf

for i,(zmin,zmax) in enumerate(zbins):
    mdf = mdf_all[(mdf_all[mocksz] > zmin)&(mdf_all[mocksz] < zmax)]
    cdf = cdf_all[(cdf_all[candelsz] > zmin)&(cdf_all[candelsz] < zmax)] 

    # mocks
    mdf_bands = ['wfc3f160w_dust_mag','wfc3f160w_dust_mag_noisy']
    colors = ['red','orange','blue','green']
    labels=['raw','raw + noise','CC','CC + noise']
    for j,mdf_band in enumerate(mdf_bands):
        ### RAW
        # no baddies
        sel = np.isfinite(mdf[mdf_band])
        # grab counts and bincenters
        x,y = returnxy(mdf[sel],mdf_band,bins)
        # normalize to field size and binsize
        y = y/dbin/mocks_fieldsize
        # plot mocks
        ax[i].plot(x,np.log10(y),c=colors[j],label=mdf_band)
        ### COMPLETENESS CORRECTED
        det = (mdf['detected?'] == 1)
        # grab counts and bincenters
        x,y = returnxy(mdf[sel&det],mdf_band,bins)
        NMC = 10
        mean_cs,bs,stds,allcs = montecarlo(mdf[sel],mdf_band,bins,NMC)    
	# normalize to field size and binsize
        y = mean_cs
        y = y/dbin/mocks_fieldsize
        yerr = stds
        yerr = yerr/dbin/mocks_fieldsize
        # plot mocks
        if (j+2 == 3):
            ax[i].errorbar(x,np.log10(y),yerr=yerr,c=colors[j+2]) 
        else:
            ax[i].plot(x,np.log10(y),c=colors[j+2],label='detected %s' %mdf_band)
        

    # obs. data
    cdf_band = 'WFC3_F160W_MAG'
    # no baddies
    sel = np.isfinite(cdf[cdf_band])
    # counts,bincenters
    x,y = returnxy(cdf[sel],cdf_band,bins)
    # normalize
    y = y/dbin/candels_fieldsize
    # plot
    ax[i].plot(x,np.log10(y),'ko',ms=4,label='CANDELS')

    # formatting
    ax[i].set_title(r'$%s < z < %s$' %(zmin,zmax))

# formatting
#plt.legend(loc='best')
ax[7].set_xlabel('WFC3F160W (AB mag)')
ax[3].set_ylabel(r'$N \; \mathrm{arcmin^{-1}} \; \mathrm{mag^{-1}}$')
# save
plt.tight_layout()
plt.savefig('WFC3F160W_counts_truez.png',dpi=300)





