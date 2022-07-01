#!/usr/bin/env python3

import xarray

def main(file1, file2):
    ds1 = xarray.open_dataset(file1)
    ds2 = xarray.open_dataset(file2)

    for v in ds1.variables.keys():
        #print(f'{v:10} file1: {ds1[v].min().values:10} - {ds1[v].max().values:10}; file2: {ds2[v].min().values:10} - {ds2[v].max().values:10}')
        print(f'{v} file1: {ds1[v].min().values} - {ds1[v].max().values}; file2: {ds2[v].min().values} - {ds2[v].max().values}')
    ds1.close()
    ds2.close()

if __name__ == '__main__':
    import plac; plac.call(main)
