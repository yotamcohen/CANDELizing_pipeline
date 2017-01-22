import sys
sys.path.insert(0,'/Users/yotam/research/SAM-CANDELS_pipeline/plots/')
from includes import *
init_plotting()

cols = np.loadtxt('filters.txt',unpack=True,dtype=str)
cdf = read_hdf('/Users/yotam/research/mocks/pipeline/lib/CANDELS/goodsn.h5')
cdf = cdf[cdf.PhotFlag == 0]

# set up figure
fig,axes = plt.subplots(3,3,figsize=(6,6),sharex=True,sharey=True)
ax = axes.flatten()
# loop over filterbands
for i,col in enumerate(cols):
    # take out baddies
    df = cdf[cdf[col+'_FLUX'] > -99]
    # grab fluxes and magnitudes for CANDELS
    cflux = df[col+'_FLUX']
    cfluxerr = df[col+'_FLUXERR']
    cmag = -2.5*np.log10(cflux) + 23.9
    # define quantities to plot
    x = cmag
    ##x = -2.5*np.log10(df['WFC3_F160W_FLUX']) + 23.9
    y = np.log10(cfluxerr/cflux)
    # debugging
    print col,np.min(x),np.max(x)
    
    # scatter plot
    ax[i].scatter(x,y,s=1,lw=0,c='lightgray')
    # bin up y to show averages
    mbins = np.linspace(np.min(x),np.max(x),10)
    df['mbin'] = np.digitize(x,mbins)
    xavg,yavg = [],[]
    for mbin in np.unique(df.mbin):
        # subset in this bin
        inbin = (df.mbin == mbin)
        # average x value in this bin
        xavg.append(np.mean(x[inbin]))
        # average y value in this bin
        yavg.append(np.mean(y[inbin]))
    # plot averages over scatter points
    ax[i].plot(xavg,yavg,'b.-',ms=4)

    # set x label
    ax[i].set_xlabel('%s (magAB)' %col,fontsize=8) 
    ##ax[7].set_xlabel('WFC3F160W (magAB)',fontsize=12) 
    # set y label
    ax[3].set_ylabel(r'$\log{[\, \delta F / F \,]}$',fontsize=12)
    


plt.tight_layout()
#plt.savefig('deltaFbyF_vs_mag_GOODSN.png',dpi=300)
plt.savefig('deltaFbyF_vs_mag_GOODSN.pdf')

