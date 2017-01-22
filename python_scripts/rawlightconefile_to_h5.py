# last update: 10/19/2016

from pandas import *
import sys
import os


def main(headerfile,datfile,h5file):
    "simply take raw ASCII lightcone.dat files and turn into h5 files"

    # get column names from header file
    # rows should be in the following format:
    # <number> <column name (no spaces, only underscores)> <optional info like units>
    columns = []
    with open(headerfile,'r') as f:
        for line in f:
            l = line.split(' ')
            c = l[2].rstrip('\n')
            columns.append(c)
    
    # now, use pandas to read in the raw lightcone file
    df = read_table(datfile,comment='#',header=None,names=columns,
                        sep='\s*',engine='python')
    print 'grabbed lightcone of shape', df.shape
    # grab filename and write dataframe to h5 file at specified path
    store = HDFStore(h5file)
    store['dat'] = df
    store.close()
    print 'writing to h5 file...'
    print 'saved lightcone data to %s' %h5file
    print 'all units are same as in original lightcone data file'
    # return the dataframe if you want to play with it
    return df

if __name__ == '__main__':
    headerfile = sys.argv[1]
    datfile = sys.argv[2]
    h5file = sys.argv[3]
    df = main(headerfile,datfile,h5file)    
