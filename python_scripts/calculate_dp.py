# Last update: 01/11/2017


from __future__ import division
import sys
import os
from pandas import *
import pandas as pd
pd.set_option('chained_assignment',None)
import numpy as np
import datetime
vega_to_AB = 1.2514

def main(h5file,disk_fname,ellip_fname):
    "Calculate the detection probability (dp) of galaxies in the catalog."

    print 'working on %s' %h5file

    #=====================================================
    # get completeness data
    #=====================================================
    meshdict_disk = np.load(disk_fname)
    Xd = meshdict_disk['X'] + vega_to_AB
    Yd = meshdict_disk['Y']
    Hd = meshdict_disk['H']
    meshdict_ellip = np.load(ellip_fname)
    Xe = meshdict_ellip['X'] + vega_to_AB
    Ye = meshdict_ellip['Y']
    He = meshdict_ellip['H']

    # REPLACE -0 WITH NAN
    Hd[np.isneginf(1/Hd)] = np.nan
    He[np.isneginf(1/He)] = np.nan

    # get master lists
    Xds = []
    for block in Xd:
        _ = (block[1:]+block[:-1])/2
        Xds.append(_)
    Yds = []
    for block in Yd:
        _ = (block[1:]+block[:-1])/2
        Yds.append(_)
    Xds,Yds = np.array(Xds),np.array(Yds)
    Xds,Yds = Xds.flatten(),Yds.flatten()
    Hds = Hd.flatten()
    completeness_d = np.array(zip(Xds,Yds,Hds)) # master list
    Xes = []
    for block in Xe:
        _ = (block[1:]+block[:-1])/2
        Xes.append(_)
    Yes = []
    for block in Ye:
        _ = (block[1:]+block[:-1])/2
        Yes.append(_)
    Xes,Yes = np.array(Xes),np.array(Yes)
    Xes,Yes = Xes.flatten(),Yes.flatten()
    Hes = He.flatten()
    completeness_e = np.array(zip(Xes,Yes,Hes)) # master list

    #print np.min(Xd),np.min(Xe),np.max(Xd),np.max(Xe)
    #print np.min(Yd),np.min(Ye),np.max(Yd),np.max(Ye)

    # get bin extents in both variables
    uXds = Xd[0]
    uYds = Yd[:,0]
    uXes = Xe[0]
    uYes = Ye[:,0]
    dXd = np.mean([(uXds[i+1] - uXds[i]) for i in range(len(uXds)-1)])
    dYd = np.mean([(uYds[i+1] - uYds[i]) for i in range(len(uYds)-1)])
    dXe = np.mean([(uXes[i+1] - uXes[i]) for i in range(len(uXes)-1)])
    dYe = np.mean([(uYes[i+1] - uYes[i]) for i in range(len(uYes)-1)])

    #========================================================
    # grab input catalog
    #========================================================
    store = HDFStore(h5file)
    catalog = store.dat
    store.close()
    #catalog['dp'] = np.nan # add detection probability column
    catalog['dp'] = np.nan # add detection probability column
    ## parameters in completeness space
    lReff = np.log10(catalog['r_pixels'])
    #Hmag = catalog['wfc3f160w']
    #Hmag = catalog['wfc3f160w_dust']
    #Hmag = catalog['wfc3f160w_dust_noisy'] # new way, july2016
    Hmag = catalog['wfc3f160w_dust_mag_noisy'] # jan2017
    ## create selections based on morphology
    catalog_d = catalog['n_sersic'] <= 2.5
    catalog_e = catalog['n_sersic'] > 2.5

    ## deal with objects outside the completeness data parameter space (baddies)
    Hmin, Hmax = 20+vega_to_AB,28+vega_to_AB
    lRmin,lRmax = -1.0,2.5
    catalog['dp'][(Hmag < Hmin)] = 1.0
    catalog['dp'][(lReff < 0) & (Hmag < Hmax)] = 1.0 # point sources
    catalog['dp'][(Hmag > Hmax) | (lReff > lRmax)] = 0.0 # dim or huge or both
    
    ## calculate dp for all the rest
    for tup in completeness_d:
        Xsel = (Hmag > (tup[0] - dXd)) & (Hmag < (tup[0] + dXd))
        Ysel = (lReff > (tup[1] - dYd)) & (lReff < (tup[1] + dYd))
        catalog_sel = catalog[(Xsel & Ysel) & catalog_d]
        #print tup[0],tup[1],tup[2],catalog_sel.shape[0]
        catalog['dp'][(Xsel & Ysel) & catalog_d] = tup[2]
    for tup in completeness_e:
        Xsel = (Hmag > (tup[0] - dXe)) & (Hmag < (tup[0] + dXe))
        Ysel = (lReff > (tup[1] - dYe)) & (lReff < (tup[1] + dYe))
        catalog_sel = catalog[(Xsel & Ysel) & catalog_e]
        #print tup[0],tup[1],tup[2],catalog_sel.shape[0]
        catalog['dp'][(Xsel & Ysel) & catalog_e] = tup[2]

    # handle the really edge case of nan dp cells, but good Hband
    nandp = np.isnan(catalog['dp'])
    nanH = np.isnan(Hmag)
    catalog['dp'][nandp & ~nanH] = 0.0

    """
    # change -0 to 0
    catalog['dp'].replace(-0,0,inplace=True)
    """

    store = HDFStore(h5file)
    store['dat'] = catalog
    store.close()

    #========================================================
    # get quick stats regarding the detection probabilities
    #========================================================
    ps = np.arange(0,1,0.1)
    for p in ps:
	pcent = catalog[catalog.dp > p].shape[0] / catalog.shape[0]
	#print p, np.round(pcent,4)
    perfectfrac = catalog[catalog.dp == 1.0].shape[0] / catalog.shape[0]
    #print '1.0',np.round(perfectfrac,4),'\n'

    # return catalog dataframe to play with (optional)
    print 'dataframe now includes detection probabilities'
    return catalog


if __name__ == '__main__':
    h5file = sys.argv[1]
    diskfile = sys.argv[2]
    ellipfile = sys.argv[3]
    df = main(h5file,diskfile,ellipfile)



