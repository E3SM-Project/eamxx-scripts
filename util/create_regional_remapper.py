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
#      If instead of a box the user desires a circular boundary around
#      a center point the CSV file can have the following format,
#              site_name, lon, lat, radius
#      e.g.
#        site 1, 35.0, 20.0, 5.0
#      Finally, it is okay to have a CSV file with a mix of these two
#      formats.  Any row with 3 real values will be interpreted as a
#      circular domain.  4 real values will result in a box.
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
import sys

class Grid:
  lat       = np.array([])
  lon       = np.array([])
  chunk_idx = 0
  size      = 0
  filename  = ""
  lat_name  = ""  # Name for lat coordinate in file
  lon_name  = ""  # Name for lon coordinate in file
  grid_type = ""  # Type of grid file
  # Additional SCRIP variables
  area      = np.array([])
  corner_lat = np.array([])
  corner_lon = np.array([])
  imask     = np.array([])
  dims      = np.array([])
  num_corners = 0

  def __init__(self,filename_,grid_filetype_):
    self.filename = filename_
    self.grid_type = grid_filetype_
    f = netCDF4.Dataset(self.filename,"r")
    if grid_filetype_ == "scrip":
      self.size     = f.dimensions["grid_size"].size
      self.lat_name = "grid_center_lat"
      self.lon_name = "grid_center_lon"
      if "grid_corners" in f.dimensions:
        self.num_corners = f.dimensions["grid_corners"].size
    else:
      self.size     = f.dimensions["ncol"].size
      self.lat_name = "lat"
      self.lon_name = "lon"
    f.close()

  def grab_chunk(self,chunk,load_full_scrip=False):
    f = netCDF4.Dataset(self.filename,"r")
    chunk_end = np.min([self.size,self.chunk_idx+chunk]);
    self.lat  = f[self.lat_name][self.chunk_idx:chunk_end].data
    self.lon  = f[self.lon_name][self.chunk_idx:chunk_end].data
    
    # Load additional SCRIP variables if requested
    if load_full_scrip and self.grid_type == "scrip":
      self.area = f["grid_area"][self.chunk_idx:chunk_end].data
      if "grid_corner_lat" in f.variables:
        self.corner_lat = f["grid_corner_lat"][self.chunk_idx:chunk_end,:].data
        self.corner_lon = f["grid_corner_lon"][self.chunk_idx:chunk_end,:].data
      if "grid_imask" in f.variables:
        self.imask = f["grid_imask"][self.chunk_idx:chunk_end].data
      if self.chunk_idx == 0 and "grid_dims" in f.variables:
        self.dims = f["grid_dims"][:].data
    
    self.chunk_idx = chunk_end
    f.close()

class Site:
  name     = ""
  inv_lon_bnds = False
  lon_bnds = np.zeros(2)
  lat_bnds = np.zeros(2)
  gids     = np.array([])
  lon_gids = np.array([])
  lat_gids = np.array([])
  lons     = np.array([])
  lats     = np.array([])
  offset   = 0
  length   = 0
  my_id    = -1
  # Additional SCRIP grid data for regional output (separate for lon and lat filtering)
  lon_areas    = np.array([])
  lon_imasks   = np.array([])
  lon_corner_lats = np.array([])
  lon_corner_lons = np.array([])
  lat_areas    = np.array([])
  lat_imasks   = np.array([])
  lat_corner_lats = np.array([])
  lat_corner_lons = np.array([])
  # Final filtered SCRIP data
  areas    = np.array([])
  corner_lats = np.array([])
  corner_lons = np.array([])
  imasks   = np.array([])

  def __init__(self,name_,bnds,bnds_type):
    self.name = name_
    if bnds_type == "box":
      self.lat_bnds = bnds[2:]
      self.lon_bnds = bnds[:2]
    elif bnds_type == "radius":
      slon = bnds[0]
      slat = bnds[1]
      srad = bnds[2]
      self.lat_bnds = slat + np.array([-1,1])*srad
      self.lon_bnds = slon + np.array([-1,1])*srad

  def filter_gids(self,filter_scrip=False):
    mask_lon = np.isin(self.lon_gids,self.lat_gids)
    mask_lat = np.isin(self.lat_gids,self.lon_gids)
    self.gids = self.lon_gids[mask_lon]
    self.lons = self.lons[mask_lon]
    self.lats = self.lats[mask_lat]
    if filter_scrip:
      self.areas = self.lon_areas[mask_lon]
      self.imasks = self.lon_imasks[mask_lon]
      if len(self.lon_corner_lats) > 0:
        self.corner_lats = self.lon_corner_lats[mask_lon,:]
        self.corner_lons = self.lon_corner_lons[mask_lon,:]
    
    # Sort by global IDs to preserve original grid order
    sort_idx = np.argsort(self.gids)
    self.gids = self.gids[sort_idx]
    self.lons = self.lons[sort_idx]
    self.lats = self.lats[sort_idx]
    if filter_scrip:
      self.areas = self.areas[sort_idx]
      self.imasks = self.imasks[sort_idx]
      if len(self.corner_lats) > 0:
        self.corner_lats = self.corner_lats[sort_idx,:]
        self.corner_lons = self.corner_lons[sort_idx,:]

class Sites:
  m_sites = []

  def __init__(self,filename):
    with open(filename, newline='') as csvfile:
      csv_in = csv.reader(csvfile, delimiter=',')
      for row in csv_in:
        bnds = np.array([float(x) for x in row[1:]])
        if len(row[1:])==4:
            # Then a full box has been defined
            l_site = Site(row[0],bnds,"box")
        elif len(row[1:])==3:
            # Then a site location plus radius is defined
            l_site = Site(row[0],bnds,"radius")
        else:
            raise
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
  f = netCDF4.Dataset(casename+"_map.nc","w",format="NETCDF3_CLASSIC")
  f.createDimension('n_a',grid.size)
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
      csvwriter.writerow([s.my_id,s.name,str(s.offset+1),str(s.offset+s.length)])

def construct_regional_grid(casename,sites,grid):
  """Create a SCRIP format grid file containing only the regional subset."""
  # Combine all sites into a single regional grid
  total_size = sum(s.length for s in sites)
  
  # Allocate arrays
  lats = np.zeros(total_size)
  lons = np.zeros(total_size)
  areas = np.zeros(total_size)
  imasks = np.zeros(total_size, dtype=np.int32)
  has_corners = any(len(s.corner_lats) > 0 for s in sites)
  num_corners = sites[0].corner_lats.shape[1] if has_corners else 0
  corner_lats = np.zeros((total_size, num_corners)) if has_corners else None
  corner_lons = np.zeros((total_size, num_corners)) if has_corners else None
  
  # Fill arrays from sites
  offset = 0
  for s in sites:
    n = s.length
    lats[offset:offset+n] = s.lats
    lons[offset:offset+n] = s.lons
    areas[offset:offset+n] = s.areas
    imasks[offset:offset+n] = s.imasks
    if has_corners:
      corner_lats[offset:offset+n,:] = s.corner_lats
      corner_lons[offset:offset+n,:] = s.corner_lons
    offset += n
  
  # Create regional grid netCDF file in SCRIP format
  f = netCDF4.Dataset(casename+"_regional_grid.nc","w",format="NETCDF3_CLASSIC")
  f.createDimension('grid_size', total_size)
  if has_corners:
    f.createDimension('grid_corners', num_corners)
  f.createDimension('grid_rank', 1)
  
  # Create variables
  v_area = f.createVariable('grid_area', np.float64, ('grid_size',))
  v_area.units = "radians^2"
  
  v_clat = f.createVariable('grid_center_lat', np.float64, ('grid_size',), fill_value=9.96920996838687e+36)
  v_clat.units = "degrees"
  
  v_clon = f.createVariable('grid_center_lon', np.float64, ('grid_size',), fill_value=9.96920996838687e+36)
  v_clon.units = "degrees"
  
  if has_corners:
    v_cornlat = f.createVariable('grid_corner_lat', np.float64, ('grid_size', 'grid_corners'), fill_value=9.96920996838687e+36)
    v_cornlat.units = "degrees"
    
    v_cornlon = f.createVariable('grid_corner_lon', np.float64, ('grid_size', 'grid_corners'), fill_value=9.96920996838687e+36)
    v_cornlon.units = "degrees"
  
  v_imask = f.createVariable('grid_imask', np.int32, ('grid_size',))
  
  v_dims = f.createVariable('grid_dims', np.int32, ('grid_rank',))
  
  # Write data
  v_area[:] = areas
  v_clat[:] = lats
  v_clon[:] = lons
  if has_corners:
    v_cornlat[:,:] = corner_lats
    v_cornlon[:,:] = corner_lons
  v_imask[:] = imasks
  v_dims[:] = [total_size]
  
  # Add global attributes
  f.api_version = 5.0
  f.version = 5.0
  f.floating_point_word_size = 8
  f.file_size = 0
  
  f.close()
  print(f"Regional grid file created: {casename}_regional_grid.nc")
  print(f"  Total grid cells: {total_size}")
  

def main():
  help = textwrap.dedent("""
  HELP
  """)
  p = argparse.ArgumentParser(description=help, formatter_class=argparse.RawTextHelpFormatter)
  p.add_argument('grid_file', type=str,
                  help='Path to a netCDF file containing the source grid information for remapping.')
  p.add_argument('grid_filetype', type=str, default='scrip', choices=['scrip','init'],
                  help='Type of grid file containing grid information.  This helps the tool recognize what the variable names are in the file')
  p.add_argument('site_info', type=str,
                  help='Path to a csv format file that contains the lat/lon bounds for regional remapping.')
  p.add_argument('casename', type=str,
                  help='Case name for regional remapping, this will used for the file names of output.')
  p.add_argument('chunksize', type=int, default=48602,
                  help='Option to read grid data in chunks, useful for large simulation grids')
  p.add_argument('--output-regional-grid', action='store_true',
                  help='Output a SCRIP format grid file containing only the regional subset')
  m_args = p.parse_args(); 

  # Load the data from args
  src_grid    = Grid(m_args.grid_file,m_args.grid_filetype)
  remap_sites = Sites(m_args.site_info) 

  # Determine which gids in source grid match each site
  load_scrip = m_args.output_regional_grid and src_grid.grid_type == "scrip"
  print("Finding global dof's that could fit site dimensions...")
  while(src_grid.chunk_idx < src_grid.size):
    chnk_start = src_grid.chunk_idx+1
    chnk_end   = np.amin([src_grid.chunk_idx+m_args.chunksize,src_grid.size])
    chnk_perc  = np.round(chnk_end/src_grid.size*100,decimals=2)
    print("  Searching cols (%4.2f%%): %12d - %12d, out of %12d total" %(chnk_perc,chnk_start,chnk_end,src_grid.size))
    chunk_idx = src_grid.chunk_idx
    src_grid.grab_chunk(m_args.chunksize, load_full_scrip=load_scrip)
    for idx, isite in enumerate(remap_sites.m_sites):
      isite.lon_bnds[0]=isite.lon_bnds[0] % 360
      isite.lon_bnds[1]=isite.lon_bnds[1] % 360
      if (isite.lon_bnds[0]>isite.lon_bnds[1]):
        lids = np.where( (src_grid.lon >= isite.lon_bnds[0]) | 
                         (src_grid.lon <= isite.lon_bnds[1])) 
      else:
        lids = np.where( (src_grid.lon >= isite.lon_bnds[0]) & 
                         (src_grid.lon <= isite.lon_bnds[1])) 
      isite.lons     = np.append(isite.lons, src_grid.lon[lids[0]])
      isite.lon_gids = np.append(isite.lon_gids, lids[0]+chunk_idx)
      
      # Store SCRIP data for lon matches
      if load_scrip:
        isite.lon_areas = np.append(isite.lon_areas, src_grid.area[lids[0]])
        isite.lon_imasks = np.append(isite.lon_imasks, src_grid.imask[lids[0]])
        if len(src_grid.corner_lat) > 0:
          if len(isite.lon_corner_lats) == 0:
            isite.lon_corner_lats = src_grid.corner_lat[lids[0],:]
            isite.lon_corner_lons = src_grid.corner_lon[lids[0],:]
          else:
            isite.lon_corner_lats = np.vstack([isite.lon_corner_lats, src_grid.corner_lat[lids[0],:]])
            isite.lon_corner_lons = np.vstack([isite.lon_corner_lons, src_grid.corner_lon[lids[0],:]])

      lids = np.where( (src_grid.lat >= isite.lat_bnds[0]) & 
                       (src_grid.lat <= isite.lat_bnds[1]))
      isite.lats = np.append(isite.lats, src_grid.lat[lids[0]])
      isite.lat_gids = np.append(isite.lat_gids, lids[0]+chunk_idx)
      
      # Store SCRIP data for lat matches
      if load_scrip:
        isite.lat_areas = np.append(isite.lat_areas, src_grid.area[lids[0]])
        isite.lat_imasks = np.append(isite.lat_imasks, src_grid.imask[lids[0]])
        if len(src_grid.corner_lat) > 0:
          if len(isite.lat_corner_lats) == 0:
            isite.lat_corner_lats = src_grid.corner_lat[lids[0],:]
            isite.lat_corner_lons = src_grid.corner_lon[lids[0],:]
          else:
            isite.lat_corner_lats = np.vstack([isite.lat_corner_lats, src_grid.corner_lat[lids[0],:]])
            isite.lat_corner_lons = np.vstack([isite.lat_corner_lons, src_grid.corner_lon[lids[0],:]])
            
  # Consilidate the gids for each site
  for idx, isite in enumerate(remap_sites.m_sites):
      print("Setting gids for site: %s\n    with bounds [%f,%f] x [%f,%f]" %(isite.name,isite.lon_bnds[0],isite.lon_bnds[1],isite.lat_bnds[0],isite.lat_bnds[1]))
      isite.filter_gids(filter_scrip=load_scrip)
  # Construct output files
  construct_remap(m_args.casename,remap_sites.m_sites,src_grid)
  
  # Construct regional grid file if requested
  if m_args.output_regional_grid:
    if src_grid.grid_type == "scrip":
      construct_regional_grid(m_args.casename,remap_sites.m_sites,src_grid)
    else:
      print("Warning: --output-regional-grid only supported for SCRIP format grid files. Skipping.")

if __name__ == '__main__':
    main()
