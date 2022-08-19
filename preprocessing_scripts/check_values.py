#!/usr/bin/env python3
import xarray
from subprocess import call
def main(inputfile):

    if True:
        # This would be nice but takes forever on big datasets
        with xarray.open_dataset(inputfile) as ds:

            print('{:<25} {:<15} {:<15} {:<15}'.format('var', 'minval', 'maxval', 'numnan'))
            print(''.join(['-' for i in range(25+15+15+15)]))
            for v in ds.variables:
                try:
                    minval = ds[v].min().values
                    maxval = ds[v].max().values
                    numnan = ds[v].isnull().sum().values
                    print(f'{v:<25} {minval:<15.5} {maxval:<15.5} {numnan:<15}')
                except:
                    print(f'Failed to process {v}, continuing...')
                    continue

    else:
        with xarray.open_dataset(inputfile) as ds:
            varnames = ds.variables.keys()

        def ncmax(inputfile, v, tmpfile='foo.nc'):
            call(f'ncwa -y max -O -v {v} {inputfile} {tmpfile}'.split(' '))
            call(f'ncks --trd -H -C -v {v} {tmpfile}'.split(' ')) # | cut -f 3- -d \' \''.split(' '))

        tmpfile = '/global/cscratch1/sd/bhillma/foo.nc'
        for v in varnames:
            print(v)
            ncmax(inputfile, v, tmpfile=tmpfile)

if __name__ == '__main__':
    import plac; plac.call(main)
