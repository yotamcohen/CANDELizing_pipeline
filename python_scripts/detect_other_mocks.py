# Last update: 09/12/2016

from __future__ import division
from pandas import *
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import time
import sys
import os


#==================================================================
# MC error estimation function
#==================================================================
def detfracErrs(df,N):
    """ this should basically do the MC detection method N many times
    in order to get an estimate on the random uncertainty in the
    detected fraction """
    detfracs = []
    for n in range(N):
        df['rand'] = random(df.shape[0])
        df['detected?'] = np.nan
        df['detected?'][df.rand <= df.dp] = 1
        df['detected?'][~(df.rand <= df.dp)] = 0
        df['detected?'] = df['detected?'].astype(int)
        detectedfrac = (df[df['detected?']==1].shape[0])/(df.shape[0])
        detfracs.append(detectedfrac)
    mean,std = np.mean(detfracs),np.std(detfracs)
    #print 'mean,std: %f,%f' %(mean,std)
    return mean,std


def main(field,h5file,mainparam='wfc3f160w_dust_noisy'):
    
    #==================================================================
    # grab SC mock
    #==================================================================
    #zbins = [(0.0,1.0),(1.0,2.0),(2.0,3.0),(3.0,4.0)]
    
    # grab columns to take from reference file
    with open('lib/detfracs_takecols.txt') as f:
        goodcols = f.readlines()
    goodcols = [x.rstrip('\n') for x in goodcols]
    # grab mocks catalogs
    datdir = '/Volumes/G-RAID/yotam/SAM-CANDELS_pipeline/catalogs/SC'
    store = HDFStore('%s/%s.h5' %(datdir,field))
    maindf = store.dat
    store.close()
    maindf['lmstar'] = np.log10(maindf['mstar']*1E10)
    maindf = maindf[maindf[mainparam] < 99] ### cut for observed magnitudes non-null valued
    
    # find endpoints of range
    Mmin_global = np.min(maindf[mainparam])
    Mmax_global = np.max(maindf[mainparam])
    n_Ms = 21
    Ms = np.linspace(Mmin_global,Mmax_global,n_Ms)
    Mbins = np.array([(l,r) for l,r in zip(Ms[:-1],Ms[1:])])
    Mbincenters = np.array((Ms[:-1]+Ms[1:])/2.)

    # loop over bins
    detfracs,errs = [],[]
    for i,(Mmin,Mmax) in enumerate(Mbins):
        Msel = (maindf[mainparam] >= Mmin) & (maindf[mainparam] <= Mmax)
        Mdf = maindf[Msel]
        
        if Mdf.shape[0] == 0: ### empty magnitude bin
            detfrac,std = np.nan,np.nan
        else:
            # send Mdf off to detfracErrs
            detfrac,std = detfracErrs(Mdf,5)
            # return from detfracErrs
            
        # append arrays to arrays for this magnitude bin
        detfracs.append(detfrac)
        errs.append(std)

    #==================================================================
    # plot
    #==================================================================
    fig = plt.figure(figsize=(6,4))
    plt.plot(Mbincenters,detfracs,'ro-',ms=3,lw=1.0)
    plt.errorbar(Mbincenters,detfracs,yerr=errs,fmt='none')
    plt.ylim(-0.1,1.1)
    
    fig.text(0.5,0.04,mainparam,ha='center',va='center')
    fig.text(0.04,0.5,'detected fraction',ha='center',va='center',rotation='vertical')
    fig.savefig('detfracs_SC_%s_%s.png' %(field,mainparam),dpi=500)


    print 'DONE WITH %s' %mainparam


    """ use information from SC detected fractions to attempt to
    CANDELize the other groups' mocks """
    # grab mock catalog
    print '\n\n working on %s... \n\n' %h5file
    store = HDFStore(h5file)
    odf = store.dat
    store.close()
    odf['lmstar'] = np.log10(odf['mstar']*1E10)

    # add detected? column
    odf['detected?'] = np.nan
    
    # loop over bins of main parameter
    for j,(Mmin,Mmax) in enumerate(Mbins):
        oMsel = (odf[mainparam] >= Mmin) & (odf[mainparam] <= Mmax)
        oMdf = odf[oMsel]
        print '%s < M < %s, size = %d' %(Mmin,Mmax,oMdf.shape[0])
        
        if oMdf.shape[0] == 0:
            continue
        
        oMdfinds = oMdf.index
        detfrac = detfracs[j]
        print 'detected fraction = %f' %detfrac
        if np.isnan(detfrac):
            continue

        det = oMdf.sample(frac=detfrac)
        detinds = det.index
        notdetinds = np.array(filter(lambda x: x not in detinds, oMdfinds))
        assert len(detinds) + len(notdetinds) == len(oMdfinds)

        odf['detected?'].ix[detinds] = 1
        odf['detected?'].ix[notdetinds] = 0
        
    
    print '\n\n writing to %s... \n\n' %h5file
    store = HDFStore(h5file)
    store['dat'] = odf
    store.close()

    print 'ALL DONE!'



if __name__ == '__main__':
    fields = ['cosmos','egs','goodsn','goodss','uds']
    datdir = '/Volumes/G-RAID/yotam/SAM-CANDELS_pipeline/catalogs'
    for field in fields:
        #main(field,'data/SAGE/%s.h5' %field,mainparam='wfc3f160w_dust_noisy')
        #main(field,'data/LU/%s.h5' %field,mainparam='wfc3f160w_dust_noisy')
        #main(field,'%s/SAGE/%s.h5' %(datdir,field),mainparam='lmstar')
        #main(field,'%s/LU/%s.h5' %(datdir,field),mainparam='lmstar')
        main(field,'%s/SAGE/%s.h5' %(datdir,field),mainparam='wfc3f160w_dust_noisy')
        main(field,'%s/LU/%s.h5' %(datdir,field),mainparam='wfc3f160w_dust_noisy')
