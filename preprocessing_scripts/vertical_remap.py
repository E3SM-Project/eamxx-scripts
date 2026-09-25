#!/usr/bin/env python3

import xarray, numpy, argparse

def main():
    parser = argparse.ArgumentParser(description="Do vertical remap from one sigma hybrid grid to another")
    parser.add_argument("input_file"                , help="Input data file")
    parser.add_argument("output_vertical_grid_file" , help="File with vertical coordinate information for target grid")
    parser.add_argument("output_file"               , help="Output data file with remapped fields")
    parser.add_argument("--input_vertical_grid_file", help="File with vertical coordinate information for source grid, if not contained in input_file", default=None)
    parser.add_argument("--vertical_dimension"      , help="Name of vertical dimension", default='lev')
    parser.add_argument("--use_log_surface_pressure", help="Use lnsp for hybrid calculation (ECMWF IFS)", default=True)
    parser.add_argument("--fields"                  , help="Names of fields in input_file to remap", nargs="+", type=str, default=None)
    args = parser.parse_args()
    _main(args.input_file, args.input_vertical_grid_file,
          args.output_vertical_grid_file, args.output_file,
          args.vertical_dimension, args.fields, args.use_log_surface_pressure)


# Do linear interpolation with extrapolation at the endpoints
def linear_interp(x, xp, fp):

    # Start with linear interpolation within domain of xp
    y = numpy.interp(x, xp, fp)

    # Extrapolate left
    mask_left = x < xp[0]
    if numpy.any(mask_left):
        slope_left = (fp[1] - fp[0]) / (xp[1] - xp[0])
        y[mask_left] = fp[0] + (x[mask_left] - xp[0]) * slope_left

    # Extrapolate right
    mask_right = x > xp[-1]
    if numpy.any(mask_right):
        slope_right = (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
        y[mask_right] = fp[-1] + (x[mask_right] - xp[-1]) * slope_right

    return y


def _main(input_file, input_vertical_grid_file, output_vertical_grid_file,
          output_file, vertical_dimension, fields, use_log_surface_pressure):

    # Open input dataset
    # TODO: figure out how to abstract out the chunking for large datasets, or
    # allow passing this on the command line
    # Make sure we do NOT chunk on vertical dimension
    ds_in = xarray.open_dataset(input_file, chunks={vertical_dimension: -1}) #chunks={'ncells': 10000, vertical_dimension: -1, 'time': 1})

    # Open input_vertical_grid_file to get grid information, if separate from input file
    if input_vertical_grid_file is not None:
        with xarray.open_dataset(input_vertical_grid_file) as ds_vert_in:
            for v in ('hyam', 'hybm', 'hyai', 'hybi'):
                ds_in[v] = ds_vert_in[v]

    # Get list of fields to remap
    if fields is None: fields = [f for f in ds_in.variables.keys() if vert_dim in ds_in[f].dims]

    # Compute pressure on input vertical grid
    if use_log_surface_pressure and all([v in ds_in.variables.keys() for v in ('lnsp', 'hyam', 'hybm', 'hyai', 'hybi')]):
        # ECMWF IFS vertical grid
        ps = numpy.exp(ds_in['lnsp'].isel({vertical_dimension: 0}))
        p_in = (ds_in['hyam'] + ds_in['hybm'] * ps).rename({'nhym': vertical_dimension})
    elif all([v in ds_in.variables.keys() for v in ('P0', 'PS', 'hyam', 'hybm', 'hyai', 'hybi')]):
        # EAM/EAMxx vertical grid
        ps = ds_in['PS']
        p_in = ds_in['hyam'] * ds_in['P0'] + ds_in['hybm'] * ps
    else:
        raise RuntimeError('Cannot infer vertical grid type of input_file.')

    # Compute pressure on output vertical grid and figure out name of output vertical coordinate dimension
    # Also copy this information to output dataset
    with xarray.open_dataset(output_vertical_grid_file) as ds_vert:
        output_vertical_dimension = ds_vert['hyam'].dims[0]
        p_out = ds_vert['P0'] * ds_vert['hyam'] + ps * ds_vert['hybm']
        ds_out = xarray.Dataset({v: ds_vert[v] for v in ('hyam', 'hybm', 'hyai', 'hybi')})
        ds_out['PS'] = ps  # Make sure surface pressure is copied as well
        ds_out['P0'] = ds_vert['P0']

    # Remap fields to target vertical grid, linear in pressure
    for f in fields:
        # Sanity check: make sure input does not contain missing data
        if ds_in[f].isnull().any():
            raise RuntimeError('Cannot handle missing data.')

        # Do vertical interpolation vectorized over non-level dimensions
        print(p_out.dims, p_in.dims, ds_in[f].dims)
        ds_out[f] = xarray.apply_ufunc(
                numpy.interp,
                p_out, p_in, ds_in[f],
                input_core_dims=[[output_vertical_dimension], [vertical_dimension], [vertical_dimension]],
                output_core_dims=[[output_vertical_dimension]],
                vectorize=True, dask='parallelized',
                output_dtypes=[float],
            )

        # Copy attributes from input data
        ds_out[f].attrs = ds_in[f].attrs

    print('Write to file...', end='')
    ds_out.to_netcdf(output_file)
    print('done.')
    ds_in.close()
    ds_out.close()


if __name__ == '__main__':
    main()
