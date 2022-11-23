# This python routine is meant to create a regional output
# remap file for use in SCREAM.
#
# Author: Aaron S. Donahue (donahue5@llnl.gov)
# Date: 2022-11-23
#
# -----------------------------------------------------------
#
# Inputs:
#   1. grid_file: 
#      A netCDF file that has grid information for the simulation.
#      Must have the variables lat and lon.
#   2. regional_sites:
#      A CSV file that contains a list of lat/lon boxes with the format
#              site_name, lon_min, lon_max, lat_min, lat_max
#                string     Real    Real      Real     Real
#      e.g.
#        site 1, 10.0, 20.0, -10.0, 10.0
#        site 2, -20.0, -10.0, -10.0, 10.0
#        ...
#   3. case_name:
#      A string that will be used to for the file names of the output
#
# Outputs:
#   1. remap_file:
#     A netCDF file that contains a set of remap triplets for the regional
#     mapping.  The output will contain the variables:
#       col: The source grid global column id (1-based)
#       row: The target grid column id (1-based)
#       S:   The associated weight for the mapping of col->row
#       lat: The latitude in the source grid for a specific target column
#       lon: The longitude in the source grid for a specific target column
#   2. remap_key:
#     A CSV file that contains the column ids in the remap file for each site.
#     The file will have the format,
#        site_name, col_id_start, col_id_len
#          string      int           int
#     where the site_name matches the regional_sites input name,
#           the col_id_start is the first column in the remapped list corresponding
#               to this site 
#           the col_id_len is the number of remap columns used for this site.
# -------------------------------------------------------------

import argparse, textwrap, csv
import numpy as np
import netCDF4

class Grid:
  lat      = np.array([])
  lon      = np.array([])

  def __init__(self,filename):
    f = netCDF4.Dataset(filename,"r")
    self.lat = f["lat"][:].data
    self.lon = f["lon"][:].data
    f.close()

class Site:
  name     = ""
  lon_bnds = np.zeros(2)
  lat_bnds = np.zeros(2)
  gids     = np.array([])
  lons     = np.array([])
  lats     = np.array([])
  offset   = 0
  length   = 0
  my_id    = -1

  def __init__(self,name_,bnds):
    self.name = name_
    self.lon_bnds = bnds[:2] #np.array([bnds[0],bnds[1]])
    self.lat_bnds = bnds[2:] #np.array([bnds[2],bnds[3]])

class Sites:
  m_sites = []

  def __init__(self,filename):
    with open(filename, newline='') as csvfile:
      csv_in = csv.reader(csvfile, delimiter=',')
      for row in csv_in:
        bnds = np.array([float(x) for x in row[1:]])
        l_site = Site(row[0],bnds)
        self.m_sites.append(l_site)

def construct_remap(casename,sites,grid):
  # Determine the full set of remap triplets
  col  = np.array([])
  row  = np.array([])
  S    = np.array([])
  lats = np.array([])
  lons = np.array([])
  s_ids = np.array([])
  offset = 0
  idx    = 1
  for s in sites:
    n    = len(s.gids)
    col  = np.append(col,s.gids+1)
    row  = np.append(row,np.arange(n)+offset+1)
    S    = np.append(S,np.ones(n))
    lats = np.append(lats,s.lats)
    lons = np.append(lons,s.lons)
    s_ids = np.append(s_ids,np.ones(n)*idx)
    s.offset = offset
    s.length = n
    s.my_id = idx
    offset += n
    idx += 1

  # Create remap netCDF file
  f = netCDF4.Dataset(casename+"_map.nc","w")
  f.createDimension('n_a',len(grid.lat))
  f.createDimension('n_b',offset)
  f.createDimension('n_s',offset)

  v_col = f.createVariable('col', np.int32, ('n_s',))
  v_row = f.createVariable('row', np.int32, ('n_s',))
  v_S   = f.createVariable('S',   np.float32, ('n_s',))
  v_lat = f.createVariable("lat", np.float32, ('n_s',))
  v_lon = f.createVariable("lon", np.float32, ('n_s',))
  v_ids = f.createVariable("site_ID", np.int32, ('n_s',))

  v_col[:] = col
  v_row[:] = row
  v_S[:]   = S  
  v_lat[:] = lats
  v_lon[:] = lons
  v_ids[:] = s_ids

  f.close() 

  # Create key to check sites
  with open(casename+"_key.csv","w") as csvfile:
    csvwriter = csv.writer(csvfile)
    for s in sites:
      csvwriter.writerow([s.my_id,s.name,str(s.offset+1),str(s.length)])
  

def main():
  help = textwrap.dedent("""
  HELP
  """)
  p = argparse.ArgumentParser(description=help, formatter_class=argparse.RawTextHelpFormatter)
  p.add_argument('grid_file', type=str,
                  help='Path to a netCDF file containing the source grid information for remapping.')
  p.add_argument('site_info', type=str,
                  help='Path to a csv format file that contains the lat/lon bounds for regional remapping.')
  p.add_argument('casename', type=str,
                  help='Case name for regional remapping, this will used for the file names of output.')
  p.add_argument('neg_lon', default=False,
                  help='Are longitude values able to be negative, i.e. bound by [-180,180] rather than [0,360]')
  m_args = p.parse_args(); 

  # Load the data from args
  src_grid    = Grid(m_args.grid_file)
  remap_sites = Sites(m_args.site_info) 

  # Determine which gids in source grid match each site
  lon_adj = 0.0
  if m_args.neg_lon:
    lon_adj = 180.0
  for idx, isite in enumerate(remap_sites.m_sites):
    isite.lon_bnds = isite.lon_bnds + lon_adj
    print("Finding gids for site: %s\n    with bounds [%f,%f] x [%f,%f]" %(isite.name,isite.lon_bnds[0],isite.lon_bnds[1],isite.lat_bnds[0],isite.lat_bnds[1]))
    lids = np.where( (src_grid.lon >= isite.lon_bnds[0]) & 
                     (src_grid.lon <= isite.lon_bnds[1]) & 
                     (src_grid.lat >= isite.lat_bnds[0]) & 
                     (src_grid.lat <= isite.lat_bnds[1]))
    isite.gids = lids[0]
    isite.lons = src_grid.lon[isite.gids]
    isite.lats = src_grid.lat[isite.gids]

  # Construct output files
  construct_remap(m_args.casename,remap_sites.m_sites,src_grid)

if __name__ == '__main__':
    main()
