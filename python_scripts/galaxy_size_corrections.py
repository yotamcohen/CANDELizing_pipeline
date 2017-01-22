# Last Update: 01/16/2016

from pandas import *
import numpy as np
import sys
import os
import warnings
warnings.filterwarnings("ignore")

def main(h5file):
    "make size corrections as prescribed by rss"

    df = read_hdf(h5file)
    print 'grabbed dataframe. working...'
    
    chi_gas = 1.7 # converts from baryonic (stars+gas) disk to stellar disk scale radius
    fproj_disk = 1.0 # from 3D to 2D projected radius for disk
    fproj_sph = 2.0*0.340 # from 3D to 2D projected radius for spheroid
    rstar_to_rv = 0.76 # convert from stellar mass radius to V-band radius
    # the extra 1.678 below converts from scale radius to reff
    rd = 1.678*fproj_disk*(df['r_disk']/chi_gas)/rstar_to_rv
    rb = fproj_sph*df['rbulge']/rstar_to_rv    
    # also include stellar mass projected sizes in final dataframe [kpc]
    df['r_disk_proj_star'] = 1.678*fproj_disk*(df['r_disk']/chi_gas)
    df['r_bulge_proj_star'] = fproj_sph*df['rbulge']
   
    # van der Wel paper
    wavdepd=-0.35+(0.12*df['redshift'])-0.25*np.log10(df['mstar'])
    wavdepb=-0.25
    zp=np.zeros(len(df['redshift']))
    zp[df['redshift']<=1.5]=1.5
    zp[df['redshift']>1.5]=2.2

    # h-band radius
    rd=rd/(((1+df['redshift'])/(1+zp))**wavdepd)
    rb=rb/(((1+df['redshift'])/(1+zp))**wavdepb)
    
    # include H-band projected sizes in final dataframe [kpc]
    df['r_disk_proj_H'] = rd
    df['r_bulge_proj_H'] = rb

    # B/T using light
    dustcorrection160=10**((df['wfc3f160w']-df['wfc3f160w_dust'])/2.5)
    wfc160dustbulge=-2.5*np.log10(dustcorrection160)+df['wfc3f160w_bulge']
    btl=10**((df['wfc3f160w_dust']-wfc160dustbulge)/2.5)
    df['BT_light'] = btl
    df['BT_mass'] = df['mbulge']/df['mstar']

    # write to file
    df.to_hdf(h5file,'dat')
    print 'saved dataframe to %s' %h5file

    return df

if __name__ == '__main__':
    h5file = sys.argv[1]
    df = main(h5file)
