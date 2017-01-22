# Last update: 10/19/2016

from __future__ import division
from pandas import *
import pandas as pd
pd.set_option('chained_assignment',None)
import numpy as np
from numpy.random import random as rand
import sys

def main(h5file,keyname='dat'):
    "simple MC detection scheme using detection probabilities"

    # grab file
    store = HDFStore(h5file)
    df = store[keyname]
    store.close()

    # assign a random number from a uniform distribution to each galaxy
    df['rand'] = rand(df.shape[0])
    # create 'detected' column
    df['detected?'] = np.nan
    # use simple MC to determine if detected
    df['detected?'][df.rand <= df.dp] = 1
    df['detected?'][~(df.rand <= df.dp)] = 0
    df['detected?'] = df['detected?'].astype(int)
    # print detected fraction
    detectedfrac = df[df['detected?'] == 1].shape[0]/df.shape[0]
    print 'detected fraction of sample in %s:' %h5file, detectedfrac
    # save to output h5 file
    store = HDFStore(h5file)
    store[keyname] = df
    store.close()
    return df # optionally play with dataframe


if __name__ == '__main__':
    h5file = sys.argv[1]
    df = main(h5file)
