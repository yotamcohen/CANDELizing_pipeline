# Last update: 01/16/2017

from __future__ import division
import sys
import os
import time
from pandas import *
import numpy as np
from scipy.io import readsav as readsav
from scipy.interpolate import RectBivariateSpline

# set cosmological parameters
import cosmolopy.distance as dist
cosmo = {'omega_M_0':0.308, 'omega_lambda_0':0.692, 'h':0.678}
cosmo = dist.set_omega_k_0(cosmo)

def convertBT(bt, Rbd):
    "get light profile parameters for galaxies in mocks"
    lookuptable = './lib/lookup_toy_BD.w_area.save'
    wuyts=readsav(lookuptable, python_dict=True)
    wuytsbt, wuytsreb = wuyts['bt_arr'], wuyts['reb_arr']
    wuytsret, wuytsn = wuyts['ret_grid'], wuyts['n_grid']
    nfunc=RectBivariateSpline(wuytsbt, wuytsreb, wuytsn)
    refunc=RectBivariateSpline(wuytsbt, wuytsreb, wuytsret)
    #print nfunc(bt, Rbd), refunc(bt, Rbd)
    return nfunc(bt, Rbd), refunc(bt, Rbd)
convertBT=np.vectorize(convertBT)

def main(h5file):
    "add relevant profile parameters to dataframe and save to specified h5 file"
    # grab df from input h5 file
    df = read_hdf(h5file)
    print 'grabbed dataframe from %s' %h5file
    # add BT ratio and BD ratio parameters
    #df['BT_mass'] = df['mbulge']/df['mstar']
    #df['RBD'] = df['rbulge']/df['r_disk'] # old version
    df['RBD'] = df['r_bulge_proj_H']/df['r_disk_proj_H'] # new version 04/24/16
    print 'dataframe now includes bulge/disk ratio from projected sizes'
    # add sersic n and effective radius
    n,re = convertBT(df['BT_light'],df['RBD'])
    #df['n'],df['re'] = n,re*df['r_disk'] # old version
    #df['re'] = re
    df['n_sersic'],df['reff'] = n,re*df['r_disk_proj_H'] # new version 03/31/16
    # ^^^ reff here should be in kpc
    print 'dataframe now includes sersic profile parameters'
    # add angular diameter distance in kpc
    df['d_a'] = dist.angular_diameter_distance(df['redshift'],z0=0,**cosmo)*10**3
    # add luminosity distance in kpc
    df['d_L'] = dist.luminosity_distance(df['redshift'],**cosmo)*10**3
    # add angular radius from small angle approximation
    df['angular_radius'] = df['reff']/df['d_a'] # in radians
    df['angular_radius'] *= 180/np.pi * 3600 # convert to arcsec
    # define resolution of wfc3f160w in arcsec/pixel
    res = 0.06
    df['r_pixels'] = df['angular_radius'] / res
    print 'dataframe now includes galaxy radii in PIXELS as seen by HSTWFC3F160W'
    # write new df to output h5 file (can of course be same as input)
    df.to_hdf(h5file,'dat')
    # return dataframe to play with (optional)
    print 'saved dataframe to %s' %h5file
    return df

if __name__ == '__main__':
    h5file = sys.argv[1]
    df = main(h5file)


