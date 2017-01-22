# Last update: 01/19/2017

from __future__ import division
from pandas import *
import pandas as pd
pd.set_option('chained_assignment',None)
import numpy as np
import sys

def main(h5file,datfilename):
    # read in library file containing final columns to write
    colsfile = './lib/CANDELized_mocks_finalcols.txt'
    with open(colsfile) as f:
        foo = f.readlines()
    # grab intrinsic column names and units (if applicable)
    keepcols = []
    colnames = []
    for line in foo:
        l = line.split()
        keepcols.append(l[0])
        colnames.append(line.rstrip('\n'))
    
    # grab dataframe
    df = read_hdf(h5file)
    df.rename(columns={'reff':'r_eff_proj_H'},inplace=True)

    # grab intersection of keepcols and df.columns
    finalcols = []
    finalcolnames = []
    for keepcol,colname in zip(keepcols,colnames):
        if keepcol in df.columns:
            finalcols.append(keepcol)
            finalcolnames.append(colname)

    # write the column names into the header of the .dat file
    with open(datfilename,'w+') as f:
        for i,line in enumerate(finalcolnames):
            headerline = '# %i %s\n' %(i,line)
            f.write(headerline)
    # write the data
    df = df[finalcols]
    df.to_csv(datfilename,sep=' ',mode='a',
              na_rep=str(np.nan),
              header=False,index=False)

    print 'done!'
    return df

if __name__ == '__main__':

    h5file = sys.argv[1]
    datfilename = sys.argv[2]
    df = main(h5file,datfilename)


