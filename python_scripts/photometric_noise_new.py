from __future__ import division
from pandas import *
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time

def main(h5file):

    # grab mocks
    store = HDFStore(h5file)
    main_mdf = store.dat
    store.close()
    print 'grabbed lightcone'

    # grab CANDELS
    store = HDFStore('./lib/CANDELS/photometry_fluxes_magnitudes.h5')
    main_cdf = store.dat
    store.close()
    print 'grabbed photometry from CANDELS'

    with open('lib/MOCKS_noise_bands.txt') as f:
        mcols = f.readlines()
    mcols = [x.rstrip('\n') for x in mcols]
    with open('lib/CANDELS_noise_bands.txt') as f:
	ccols = f.readlines()
    ccols = [x.rstrip('\n') for x in ccols]

    # add errors
    mdf = main_mdf
    for ccol,mcol in zip(ccols,mcols):
        
        # column names
        cflux,cfluxerr = ccol+'_flux',ccol+'_fluxerr'
        cmag,cmagerr = ccol+'_mag',ccol+'_magerr'
        mflux,mfluxerr = mcol+'_flux',mcol+'_fluxerr'
        mmag,mmagerr = mcol+'_mag',mcol+'_magerr'
        mfluxdust,mfluxdusterr = mcol+'_dust_flux',mcol+'dust_fluxerr'
        mmagdust,mmagdusterr = mcol+'_dust_mag',mcol+'_dust_magerr'

        # create flux and magnitude columns in df
        mdf[mmag] = mdf[mcol]
        mdf[mflux] = 10**((mdf[mcol]-23.9)/-2.5)
        mdf[mmagerr],mdf[mfluxerr] = np.nan,np.nan
        mdf[mmagdust] = mdf[mcol+'_dust']
        mdf[mfluxdust] = 10**((mdf[mcol+'_dust']-23.9)/-2.5)
        mdf[mmagdusterr],mdf[mfluxdusterr] = np.nan,np.nan

        # create fitting function from CANDELS errors
        # first throw away baddies
        notnans = (~np.isnan(main_cdf[cmag]))&(~np.isnan(main_cdf[cmagerr]))
        notinfs = (~np.isinf(main_cdf[cmag]))&(~np.isinf(main_cdf[cmagerr]))
        cdf = main_cdf[(notnans) & (notinfs)]
        print cmag,cdf.shape
        # fit fractional flux error
        popt,perr = fit(cdf[cmag],np.log10(cdf[cfluxerr]/cdf[cflux]))
        m_mean,b_mean = popt[0],popt[1]
        sigma_m_mean,sigma_b_mean = perr[0],perr[1]
        # now deal with mocks
        # TODO: do it for dust
        ## step 1: get the value of the fitting function
        fitvalue = 10**(m_mean*mdf[mmag]+b_mean)
        fitvaluedust = 10**(m_mean*mdf[mmagdust] + b_mean)
        ## step 2: multiply by true flux to get delta(flux)
        deltaf = fitvalue*mdf[mflux]
        deltafdust = fitvaluedust*mdf[mfluxdust]
        ## step 3: choose a gaussian deviate with mean 0 and spread delta(flux)
        gaussnoise = np.random.normal(0,deltaf)
        gaussnoisedust = np.random.normal(0,deltafdust)
        ## step 4: that is the random noise
        mdf[mfluxerr] = gaussnoise
        mdf[mfluxdusterr] = gaussnoisedust
        ## step 5: add the random noise to the true flux to get noisy flux
        mfluxnoisy = mflux+'_noisy'
        mdf[mfluxnoisy] = mdf[mflux] + mdf[mfluxerr]
        mfluxdustnoisy = mfluxdust+'_noisy'
        mdf[mfluxdustnoisy] = mdf[mfluxdust] + mdf[mfluxdusterr]
        ## step 6: get a noisy magnitude too
        mmagnoisy = mmag+'_noisy'
        mdf[mmagnoisy] = -2.5*np.log10(mdf[mfluxnoisy]) + 23.9
        mdf[mmagerr] = mdf[mmagnoisy] - mdf[mmag]
        mmagdustnoisy = mmagdust+'_noisy'
        mdf[mmagdustnoisy] = -2.5*np.log10(mdf[mfluxdustnoisy]) + 23.9
        mdf[mmagdusterr] = mdf[mmagdustnoisy] - mdf[mmagdust]
        ## should be done at this point
    
        # drop extraneous columns
        extracols = [mfluxerr,mmagerr,mfluxdusterr,mmagdusterr]
        for col in extracols:
            mdf.drop(col,axis=1,inplace=True)

    # save
    store = HDFStore(h5file)
    store['dat'] = mdf
    store.close()
    print 'dataframe now includes photometric noise'
    print 'saved to %s' %h5file

    return mdf


def linfit(x,m,b):
    return m*x + b

def fit(x,y):
    popt,pcov = curve_fit(linfit,x,y)
    perr = np.sqrt(np.diag(pcov))
    return popt,perr


if __name__ == '__main__':
    h5file = sys.argv[1]
    df = main(h5file)        







