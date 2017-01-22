import sys
sys.path.insert(0,'/Users/yotam/research/SAM-CANDELS_pipeline/plots/')
from includes import *
init_plotting()
from matplotlib.colors import LogNorm

diskf = 'hist2d_n_f160w_4_expdisk_XYH.npz'
md_disk = np.load('/Users/yotam/research/mocks/pipeline/lib/completeness_data/%s' %diskf)
devaucf = 'hist2d_n_f160w_4_devauc_XYH.npz'
md_devauc = np.load('/Users/yotam/research/mocks/pipeline/lib/completeness_data/%s' %devaucf)

vega_to_AB = 1.2514


#============== lookup table =======================
Xdisk,Xdevauc = md_disk['X'],md_devauc['X']
Ydisk,Ydevauc = md_disk['Y'],md_devauc['Y']
Hdisk,Hdevauc = md_disk['H'],md_devauc['H']
# REPLACE -0 WITH NAN
Hdisk[np.isneginf(1/Hdisk)] = np.nan
Hdevauc[np.isneginf(1/Hdevauc)] = np.nan

Xdisk += vega_to_AB
Xdevauc += vega_to_AB

Hdisk_unq,Hdevauc_unq = np.unique(Hdisk),np.unique(Hdevauc)
levels_disk,levels_devauc = Hdisk_unq[1::20],Hdevauc_unq[1::20]

Xmin,Xmax = np.min([Xdisk,Xdevauc]),np.max([Xdisk,Xdevauc])
Ymin,Ymax = np.min([Ydisk,Ydevauc]),np.max([Ydisk,Ydevauc])
md_extent = [Xmin,Xmax,Ymin,Ymax]
levels = np.array([0.0,0.25,0.5,0.75,0.9])

#================== CANDELS ========================
morph = 'disk'

df = read_hdf('/Users/yotam/research/mocks/pipeline/lib/CANDELS/goodsn.h5')
df = df[df.PhotFlag == 0]
if morph == 'disk':
    df = df[df.gf_n < 2.5]
elif morph == 'ellip':
    df = df[df.gf_n > 2.5]

# first get global range
lrpix,Hmag = np.log10(df['FLUX_RADIUS_2_F160W'].values),df['Hmag'].values
mXmin,mXmax = np.min(Hmag),np.max(Hmag)
mYmin,mYmax = np.min(lrpix),np.max(lrpix)
extent = [[mXmin,mXmax],[mYmin,mYmax]]
    
print Xmin,Xmax,Ymin,Ymax
print mXmin,mXmax,mYmin,mYmax

#================== figure ========================
scale = 0.5
plt.figure(figsize=(8*scale,5*scale))
cmap,alpha = 'gnuplot',1.0
vmin,vmax = 0.0,1.0
if morph == 'disk':
    Hdat = Hdisk
elif morph == 'ellip':
    Hdat = Hdevauc
plt.contour(Hdat,extent=md_extent,levels=levels,
            cmap=cmap,lw=1.0,alpha=alpha,vmin=vmin,vmax=vmax)

#im_ctrs = axes[i][1].contour(Hdevauc,extent=md_extent,levels=levels,
#            cmap=cmap,lw=1.0,alpha=alpha,vmin=vmin,vmax=vmax)
# box outline
plt.plot([Xmin,Xmax],[Ymin,Ymin],'k-',lw=1)
plt.plot([Xmax,Xmax],[Ymin,Ymax],'k-',lw=1)
plt.plot([Xmin,Xmax],[Ymax,Ymax],'k-',lw=1)
plt.plot([Xmin,Xmin],[Ymin,Ymax],'k-',lw=1)

# CANDELS
"""
plt.hist2d(Hmag,lrpix,bins=100,
           range=extent,
           cmap='binary',
           norm=LogNorm())
"""
plt.scatter(Hmag,lrpix,s=1,c='b',lw=0)
plt.xlabel('HF160W')
plt.ylabel('log[R_eff (pixels)]')
if morph == 'disk':
    plt.title('n < 2.5')
elif morph == 'ellip':
    plt.title('n > 2.5')
plt.tight_layout()
if morph == 'disk':
    plt.savefig('CANDELScompleteness_disk.png',dpi=300)
elif morph == 'ellip':
    plt.savefig('CANDELScompleteness_ellip.png',dpi=300)


