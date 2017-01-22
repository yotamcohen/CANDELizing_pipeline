import sys
sys.path.insert(0,'/Users/yotam/research/SAM-CANDELS_pipeline/plots/')
from includes import *
init_plotting_custom()

def main(mocks,md_disk,md_devauc,title):

    #============== lookup table =======================
    Xdisk,Xdevauc = md_disk['X'],md_devauc['X']
    Ydisk,Ydevauc = md_disk['Y'],md_devauc['Y']
    Hdisk,Hdevauc = md_disk['H'],md_devauc['H']
    Xdisk += vega_to_AB
    Xdevauc += vega_to_AB

    # REPLACE -0 WITH NAN
    Hdisk[np.isneginf(1./Hdisk)] = np.nan
    Hdevauc[np.isneginf(1./Hdevauc)] = np.nan

    Hdisk_unq,Hdevauc_unq = np.unique(Hdisk),np.unique(Hdevauc)
    """
    stride_disk = np.round(len(np.unique(Hdisk))/5.)
    stride_devauc = np.round(len(np.unique(Hdevauc))/5.)
    levels_disk = Hdisk_unq[1::stride_disk]
    levels_devauc = Hdevauc_unq[1::stride_devauc]
    """

    Xmin,Xmax = np.min([Xdisk,Xdevauc]),np.max([Xdisk,Xdevauc])
    Ymin,Ymax = np.min([Ydisk,Ydevauc]),np.max([Ydisk,Ydevauc])
    md_extent = [Xmin,Xmax,Ymin,Ymax]
    #print 'meshdict_extent:', md_extent
    levels = np.array([0.25,0.5,0.75,0.9])
    levels_disk = levels
    levels_devauc = levels

    #================== mocks ========================
    # first get global range
    lrpix,Hmag = np.log10(mocks['r_pixels']),mocks['wfc3f160w_dust']
    mXmin,mXmax = np.min(Hmag),np.max(Hmag)
    mYmin,mYmax = np.min(lrpix),np.max(lrpix)
    mocks_extent = [[mXmin,mXmax],[mYmin,mYmax]]
    
    #print Xmin,Xmax,Ymin,Ymax
    #print mXmin,mXmax,mYmin,mYmax

    #================== figure ========================
    nrow,ncol = len(zbins),2
    scale = 1.0
    fig,axes = plt.subplots(nrows=nrow,ncols=ncol,figsize=(5*scale,8*scale),
                            sharex=True,sharey=True)
    ax = axes.flatten()
    for i,(zmin,zmax) in enumerate(zbins):
           
        ### CANDELS data
        cmap,alpha = 'gnuplot',0.2
        lw = 1.0
        vmin,vmax = 0.0,1.0
        #vmin,vmax = np.min(levels),np.max(levels)
        # lookup table disk
        im_ctrs = axes[i][0].contour(Hdisk,extent=md_extent,levels=levels_disk,
                            cmap=cmap,lw=lw,vmin=vmin,vmax=vmax)
        # lookup table devauc
        im_ctrs = axes[i][1].contour(Hdevauc,extent=md_extent,levels=levels_devauc,
                            cmap=cmap,lw=lw,vmin=vmin,vmax=vmax)
        # box outline
        lw = 0.75
        sty = 'k--'
        for j in [0,1]:
            axes[i][j].plot([Xmin,Xmax],[Ymin,Ymin],sty,lw=lw)
            axes[i][j].plot([Xmax,Xmax],[Ymin,Ymax],sty,lw=lw)
            axes[i][j].plot([Xmin,Xmax],[Ymax,Ymax],sty,lw=lw)
            axes[i][j].plot([Xmin,Xmin],[Ymin,Ymax],sty,lw=lw)
        
        ### MOCKS
        # redshift slice
        zsel = (mocks['redshift'] >= zmin)&(mocks['redshift'] <= zmax)
        zdf = mocks[zsel]
        # n_sersic split
        n_sersic_cut = 2.5
        disks = mocks['n_sersic'] <= n_sersic_cut
        devaucs = mocks['n_sersic'] > n_sersic_cut
        
        ### HIST2d
        nbins = 100
        # disks 
        axes[i][0].hist2d(Hmag[zsel][disks],lrpix[zsel][disks],bins=nbins,
                            range=mocks_extent,norm=LogNorm(),cmap='binary')
        # devaucs
        axes[i][1].hist2d(Hmag[zsel][devaucs],lrpix[zsel][devaucs],bins=nbins,
                            range=mocks_extent,norm=LogNorm(),cmap='binary')
        ### format
        aXmin = np.min([Xmin,mXmin])
        aXmax = np.max([Xmax,mXmax])
        aYmin = np.min([Ymin,mYmin])
        aYmax = np.max([Ymax,mYmax])
        xticks = [20,22,24,26,28]
        yticks = [-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5]
        for j in [0,1]:
            axes[i][j].set_xlim(20,30)
            axes[i][j].set_ylim(-1.5,3)
            axes[i][j].set_xticks(xticks)
            axes[i][j].set_yticks(yticks)

    for i,(zmin,zmax) in enumerate(zbins):
        zstring = '%s < z < %s' %(zmin,zmax)
        axes[i][1].text(1.05,0.5,zstring,rotation=-90,
                        ha='center',va='center',transform=axes[i][1].transAxes)

    axes[0][0].set_title('expdisk (n_sersic < 2.5)')
    axes[0][1].set_title('devauc (n_sersic > 2.5)')
    
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    #fig.text(0.5,0.05,'observed H-band magnitude',ha='center',va='center')
    fig.text(0.5,0.05,r'WFC3F160W (magAB)',fontsize=12,
             ha='center',va='center')
    fig.text(0.04,0.5,r'$\log{[R_{\rm eff}]}$ [pixels]',
             rotation='vertical',ha='center',va='center')
    
    # save and finish
    #figname = '%s__completeness_table_withmocks.png' %(title.replace(' ','_'))
    figname = '%s__completeness_table_withmocks.pdf' %(title.replace(' ','_'))
    #figname = '%s__completeness_table_withmocks_highz.png' %(title.replace(' ','_'))
    #plt.savefig(figname,dpi=600)
    plt.savefig(figname)
    print 'saved figure!'

mdf = read_hdf('/Volumes/G-RAID/yotam/catalogs/SC/12.16/goodsn.h5')
diskf = 'hist2d_n_f160w_4_expdisk_XYH.npz'
meshdict_disk = np.load('/Users/yotam/research/mocks/pipeline/lib/completeness_data/%s' %diskf)
devaucf = 'hist2d_n_f160w_4_devauc_XYH.npz'
meshdict_devauc = np.load('/Users/yotam/research/mocks/pipeline/lib/completeness_data/%s' %devaucf)

#zbins = [0.1,0.5,1.24,1.72,2.15,2.57,3.0]
#zbins = [(0.0,0.5),(0.5,1.0),(1.0,2.0),(2.0,3.0)]
zbins = [(0.1,0.5),(0.5,1.24),(1.24,1.72),(1.72,2.15),(2.15,2.57),(2.57,3.0)]
vega_to_AB = 1.2514

main(mdf,meshdict_disk,meshdict_devauc,'GOODSN')



