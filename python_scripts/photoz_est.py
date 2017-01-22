# Last update: 01/16/2017

from __future__ import division
from pandas import *
import pandas as pd
pd.set_option('chained_assignment',None)
import numpy as np
import numpy.random as random
import sys


def main(h5file):
    """
    do rough photo-z estimation based on the old fashioned
    'ruler method' of Figure 10 of Dahlen et al. 2013
    """
    
    # grab dataframe
    df = read_hdf(h5file)

    # set up bins and sigmas
    Hbins = np.array([23.5,24.5,25.5,26.5])
    sigmas = np.array([0.04,0.045,0.055,0.07])
    # decide which H-band to use
    Hmag = 'wfc3f160w_dust'
    whichHbin = np.digitize(df[Hmag],Hbins)
    # set sigma based on bin
    sigma_O = np.ones((df.shape[0]))
    for i in range(len(Hbins)):
        sigma_O[whichHbin == i] = sigmas[i]
    sigma_O[whichHbin == len(Hbins)] = sigmas[-1]
    # get list of gaussian deviates with spreads sigma_O
    deviates = random.normal(0.0,sigma_O,sigma_O.shape[0])
    # calculate z_est
    z_true = df['redshift']
    z_est = z_true - deviates*(1.+z_true)
    df['z_est'] = z_est.values
    # set baddies to 0
    df['z_est'][df.z_est < 0.] = 0.

    # write
    df.to_hdf(h5file,'dat')
    return df

if __name__ == '__main__':
    h5file = sys.argv[1]
    df = main(h5file)
