#!/usr/bin/env python3
import xarray, numpy

def main(inputfile, topofile, gasfile, outputfile):

    # Stuff we need to get from cami file
    screami_from_cami = {
        'time': 'time',
        'time_bnds': 'time_bnds',
        'lat' : 'lat',
        'lon' : 'lon',
        'hyam': 'hyam',
        'hybm': 'hybm',
        'hyai': 'hyai',
        'hybi': 'hybi',
        'P0'  : 'P0',
        'qv'  : 'Q',
        'ps'  : 'PS',
        'T_mid': 'T',
        #'horiz_winds': ('U', 'V'),
        'pref_mid': 'hybm', #'lev',
        'qc': 'CLDLIQ',
        'qi': 'CLDICE',
        'qr': 'RAINQM',
        'nc': 'NUMLIQ',
        'ni': 'NUMICE',
        'nr': 'NUMRAI',
    }

    # Other things we need:
    #   gas concentrations for trace gas constituents
    #   surface geopotential
    gas_names = ('ch4', 'co', 'co2', 'n2', 'n2o', 'o2', 'o3')

    with xarray.open_dataset(inputfile) as ds_in:

        ds_out = xarray.Dataset()

        ntime = ds_in.sizes['time']
        ncol = ds_in.sizes['ncol']
        nlev = ds_in.sizes['lev']

        # Fields we can copy directly
        for scream_name, eam_name in screami_from_cami.items():
            print(scream_name)
            ds_out[scream_name] = ds_in[eam_name]

        # Handle U,V to horiz_winds(time,ncol,dim2,lev)
        print('Map U,V to horiz_winds')
        ds_out['horiz_winds'] = xarray.DataArray(
            numpy.zeros([ntime,nlev,2,ncol]),
            dims=('time','lev','dim2','ncol'),
        )
        ds_out['horiz_winds'][:,:,0,:] = ds_in['U']
        ds_out['horiz_winds'][:,:,1,:] = ds_in['V']

        # Grab phis from topo file
        with xarray.open_dataset(topofile) as ds_topo:
            if 'PHIS_d' in ds_topo.variables.keys():
                ds_out['phis'] = ds_topo['PHIS_d'].rename({'ncol_d': 'ncol'})
            else:
                ds_out['phis'] = ds_topo['PHIS']
        # Broadcast phis to include time dim
        ds_out['phis'], *__ = xarray.broadcast(ds_out['phis'], ds_out['ps'])

        # Grab trace gas concentrations from file
        with xarray.open_dataset(gasfile).isel(time=0).drop('time') as ds_gas:
            for v in gas_names:
                ds_out[v], *__ = xarray.broadcast(ds_gas[v], ds_out['ps'])

        # Permute dimensions
        ds_out = ds_out.transpose('time','ncol','dim2','lev','ilev','nbnd')

        ds_out.to_netcdf(outputfile)

if __name__ == '__main__':
    import plac; plac.call(main)
