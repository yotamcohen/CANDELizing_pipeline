# Last update: 09/19/2016

from pandas import *
import os
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import time
from matplotlib import rcParams
params = {
    'axes.labelsize':10,
    'axes.titlesize':10,
    'xtick.labelsize':8,
    'ytick.labelsize':8,
    'font.size':10,
    'font.family':'serif',
    'legend.fontsize':8,
    'legend.handlelength':3,
    'text.usetex':False,
    'savefig.dpi':1000 }
rcParams.update(params)
sky = 4*np.pi * (180*60/np.pi)**2 # arcmin^2

def get_catalogs():
    # get individual catalog data and concatenate

    maindf = DataFrame()
    fields = ['cosmos','egs','goodsn','goodss','uds']
    goodwords = ['flux','fluxerr']
    for field in fields:
        store = HDFStore('lib/CANDELS/%s.h5' %field)
        df = store.dat
        store.close()
        #///////////////////////////////////////////////////////////////////
        # selection cuts
        goodphot = df['PhotFlag'] == 0
        notstar = df['CLASS_STAR'] < 0.8
        notAGN = df['AGNFlag'] == 0
        df = df[goodphot&notstar&notAGN]
        #///////////////////////////////////////////////////////////////////
        cols = df.columns
        newcols = map(lambda x: x.lower(), cols)
        df.columns = newcols
        ### take only photometry columns
        goodcols = []
        for col in df.columns:
            for w in goodwords:
                if col.endswith(w):
                    goodcols.append(col)
        goodcols = np.unique(goodcols)
        df = df[goodcols]
        
        maindf = concat([maindf,df],axis=0,ignore_index=True)
        print '\n\n\n %s shape: ' %field,maindf.shape,'\n\n\n'
    
    maindf.replace(to_replace=-99.0,value=np.nan,inplace=True)
    newdf = DataFrame()
    for col in maindf.columns:
        if col.endswith('flux'):
            mag = -2.5*np.log10(maindf[col]) + 23.9
            newdf[col] = mag
        elif col.endswith('fluxerr'):
            flux = maindf[col.rstrip('err')]
            fluxerr = maindf[col]
            magerr = np.abs((-2.5/np.log(10))*(fluxerr/flux))
            newdf[col] = magerr
    magcols = map(lambda x: x.replace('fluxerr','magerr').replace('flux','mag'),
                  newdf.columns)
    newdf.columns = magcols

    #n_orig = maindf.shape[0]
    n_orig = newdf.shape[0]
    #print 'Original CANDELS catalog size: ', n_orig
    print 'New CANDELS catalog size: ', n_orig

    #return maindf
    return newdf

def fit_errors():
    # plot magnitude and corresponding errors

    ofile = open('lib/CANDELS_errors_fit_parameters.txt','w+')
    ofile.write('filter slope intcpt slope_error intcpt_error \n')
    hcols = ['WFC','ACS','IRAC']
    for col in candels.columns:
        if col.endswith('FLUX') and any(c in col for c in hcols):
            flux = candels[col]
            fluxerr = candels['%sERR' %col]
            ## calculate AB magnitude and error
            ABzp = 23.9
            mag = -2.5*np.log10(candels[col]) + ABzp
            magerr = np.abs((-2.5/np.log(10))*(fluxerr/flux))            
            ## define selection
            print 'before',candels.shape[0]
            selection = (~np.isnan(mag))&(~np.isnan(magerr))
            selection = selection & (~np.isinf(mag))&(~np.isinf(magerr))
            print 'after',candels[selection].shape[0]
            ## throw away nulls
            mag = mag[selection]
            magerr = magerr[selection]
            ## send off to fitting 
            popt,perr = fit(mag,np.log10(magerr))
            print popt,perr
            m,b = popt[0],popt[1]
            sigma_m,sigma_b = perr[0],perr[1]
            ## write to file
            ofile.write('%s %f %f %f %f \n' %(col.rstrip('_FLUX'),m,b,sigma_m,sigma_b))
            ## plot data
            plt.figure(figsize=(5,4))
            
            #/////////////////////////////////////////////////////////////////
            ### errors in magnitude units
            plt.yscale('log')
            #plt.scatter(mag,magerr/mag,s=2,alpha=0.75,lw=0,
            plt.scatter(mag,magerr,s=2,alpha=0.75,lw=0,
                        c=candels['CLASS_STAR'][selection],cmap='seismic')
            #/////////////////////////////////////////////////////////////////
            
            cbar = plt.colorbar()
            cbar.set_label('stellarity index',labelpad=10,fontsize=8,rotation=270)
            ## plot fit
            
            #/////////////////////////////////////////////////////////////////
            ### when doing errors in magnitude units
            xmin,xmax = np.min(mag),np.max(mag)
            _ = np.linspace(xmin,xmax,100)
            plt.plot(_,10**((m*_)+b),'k-',lw=1.5)
            #/////////////////////////////////////////////////////////////////

            ## formatting and save
            magname = col.rstrip('_FLUX')
            magerrname = magname+'_ERR'
            plt.xlabel('%s (magnitude)' %magname)
            plt.ylabel('%s (magnitude)' %magerrname)
            #plt.ylabel('fractional error (magnitude)')
            plt.savefig('plots/new_CANDELS_errors/%s.png' %magerrname,dpi=500)
            #plt.savefig('plots/new_CANDELS_errors/%s_fracerror.png' %magerrname,dpi=500)
            plt.close()
            print 'saved figure %s' %magerrname

    ofile.close()

def linfit(x,m,b):
    return m*x + b

def fit(x,y):
    popt,pcov = curve_fit(linfit,x,y)
    perr = np.sqrt(np.diag(pcov))
    return popt,perr


if __name__ == '__main__':
    """
    t1 = time.time()
    df = get_catalogs()
    t2 = time.time()
    print 'Grabbed. Took %f seconds' %(t2-t1)
    
    store = HDFStore('lib/CANDELS/photometry.h5')
    store['dat'] = df
    store.close()
    """

    store = HDFStore('lib/CANDELS/photometry.h5')
    df = store.dat
    store.close()
    ploterrs(df)





