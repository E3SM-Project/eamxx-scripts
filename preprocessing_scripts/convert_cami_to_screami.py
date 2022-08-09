#!/usr/bin/env python3
import xarray, numpy, dask

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
        'pref_mid': 'lev', #hybm
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

    with xarray.open_dataset(inputfile, decode_cf=False) as ds_in:

        ds_out = xarray.Dataset()

        ntime = ds_in.sizes['time']
        ncol = ds_in.sizes['ncol']
        nlev = ds_in.sizes['lev']

        # Fields we can copy directly
        print('Copy fields from cami')
        for scream_name, eam_name in screami_from_cami.items():
            print(f'  {eam_name} -> {scream_name}')
            if eam_name in ds_in.variables.keys():
                ds_out[scream_name] = ds_in[eam_name]
            else:
                print(f'WARNING: {eam_name} not found in inputfile; hope you did not need that one...')

        # Make sure ps field contains a time dimension
        if 'time' not in ds_out['ps'].dims:
            ds_out['ps'] = ds_out['ps'].expand_dims('time', axis=0)

        # Handle U,V to horiz_winds(time,ncol,dim2,lev)
        print('Map U,V to horiz_winds')
        ds_out['horiz_winds'] = xarray.concat([ds_in['U'], ds_in['V']], 'dim2')

        # Grab phis from topo file
        print('Grab phis from topo file')
        with xarray.open_dataset(topofile, decode_cf=False) as ds_topo:
            if 'PHIS_d' in ds_topo.variables.keys():
                print('  use PHIS_d')
                ds_out['phis'] = ds_topo['PHIS_d'].rename({'ncol_d': 'ncol'})
            else:
                print('  use PHIS')
                ds_out['phis'] = ds_topo['PHIS']
        # Broadcast phis to include time dim
        ds_out['phis'], *__ = xarray.broadcast(ds_out['phis'], ds_out['ps'])

        # Grab trace gas concentrations from file
        print('Get trace gases')
        with xarray.open_dataset(gasfile, decode_cf=False).isel(time=0).drop('time') as ds_gas:
            for v in gas_names:
                ds_out[v], *__ = xarray.broadcast(ds_gas[v], ds_out['ps'])

        # Permute dimensions
        print('Permute dimensions')
        try:
            if 'nbnd' in ds_out.dims:
                ds_out = ds_out.transpose('time','ncol','dim2','lev','ilev','nbnd')
            else:
                ds_out = ds_out.transpose('time','ncol','dim2','lev','ilev')
        except:
            print('  permute dimensions failed, but continuing anyways.')

        print('Write to netcdf')
        ds_out.to_netcdf(outputfile)

        print('Done.')

if __name__ == '__main__':
    import plac; plac.call(main)
