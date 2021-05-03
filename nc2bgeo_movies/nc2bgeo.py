#!/usr/bin/env python
"""
Module containing helper functions to convert E3SM
 netcdf output from CESM/E3SM/SCREAMv0 into 
bgeo format for ingestion by Houdini. It is based on
code from Kalina Borkiewicz (kalina@illinois.edu) and
modified for E3SM by Peter Caldwell (caldwell19@llnl.gov)
"""

#IMPORT STUFF:
#==============
import struct
import numpy as np
import pylab as pl

#DEFINE CONSTANTS:
#==============
earthRadius = 6371. #in km
Rd=287.15 #gas constant for dry air
grav=9.8  #gravity

def make_3d_pres(ps,hya,hyb):
    """
    Create time-and-space varying 3d pressure from 
    hybrid weights hya and hyb (both with only vertical 
    dimension) and surface pressure ps (with just lat and
    lon dims).
    """

    #COPY ALL ARRAYS TO HAVE 3 DIMENSIONS
    # Below, expand hya and hyb so that they have 1+# of dimensions of ps
    for i in np.arange(len(ps.shape)):
        hya=np.expand_dims(hya,axis=i+1)
        hyb=np.expand_dims(hyb,axis=i+1)
        
    # Below, using broadcasting feature in numpy to create HYA,HYB, and PS
    # which have dimensions [len(hya),ps.shape[0],ps.shape[1]] if PS is
    # regridded or [len(hya),ps.shape[0]] if PS is native (1D)
    # https://numpy.org/doc/stable/user/basics.broadcasting.html
    HYA=(hya*np.ones([1]+list(ps.shape),np.float32))
    HYB=(hyb*np.ones([1]+list(ps.shape),np.float32))
    PS= (np.ones(hya.shape,np.float32)*np.reshape(ps,[1]+list(ps.shape)))

    P=HYA*100000. + HYB*PS

    return P

def resize_latlon(lat,lon):
    """
    For 2d variables (no height coordinate), convert 1d lat and lon arrays into
    2d arrays with same dims as var. This is just a reminder function since meshgrid
    does exactly what we want already.
    """
    LAT,LON = pl.meshgrid(lat,lon) #make lat and lon have values for each point in var
    return LAT,LON

#------------------------------------------------------------
def resize_latlonlev(lat,lon,lev):
    """
    For 3d variables (ie including a height coordinate), convert 1d lat and lon
    arrays into 3d arrays with the same dims as var. pl.meshgrid does this, but with
    inscrutable options that this function just keeps track of. Note that lev should
    already be 3d since E3SM uses pressure coordinates so heights change with lat and lon.
    """
    if np.isscalar(lev) or len(lev.shape)!=3:
        raise Exception('resize_latlonlev only works with 3 dim lev array')
    
    junk,LAT,LON = pl.meshgrid(lev[:,0,0],lat,lon,indexing='ij') #make lat and lon same size as

    return LAT,LON

#------------------------------------------------------------
def resize_ncollev(lat,lon,lev):
    """
    For 3d variables (ie including a height coordinate) with a ncol dim instead of separate
    lat/lon dims (i.e. for native rather than regridded data), this function makes lat and lon
    be copied over the vertical dimension so they have the same size as var.
    """
    if np.isscalar(lev) or len(lev.shape)!=2:
        raise Exception('resize_ncollev only works with 2 dim lev array')
    
    junk,LAT = pl.meshgrid(lev[:,0],lat)
    junk,LON = pl.meshgrid(lev[:,0],lon)

    return LAT,LON
   
#------------------------------------------------------------

def latlonht2xyz(LAT,LON,HT):
    """
    Given scalars or n-dimensional arrays for latitude, longitude, and heights above
    the earth's surface, return same-dimensioned X, Y, and Z arrays giving 
    the x,y,z position of each aray element (in meters). Works on arrays for
    efficiency.
    """
    Radius = earthRadius*1000. + HT #in meters
    X = Radius * np.cos( LAT*np.pi/180. ) * np.cos(LON*np.pi/180.)
    Y = Radius * np.cos( LAT*np.pi/180. ) * np.sin(LON*np.pi/180.)
    Z = Radius * np.sin( LAT*np.pi/180.)
    return X,Y,Z

#------------------------------------------------------------

def subset_latlonht2xyz(lat,lon,hts,latslice,lonslice):
    """
    Compute x,y,z,var points from (lat,lon,hts) points that fall in the
    selected latslice x lonslice bounds (in degrees). 
    NOTE: this fn returns 1d (flattened) x,y,z outputs.
    """

    #FIRST FLATTEN ALL INPUTS TO AVOID WORRYING ABOUT DIMS WHEN LOOPING
    lat=lat.flatten()
    lon=lon.flatten()
    Var=var.flatten()

    #Making scalar hts same dim as lat and lon simplifies code later
    if np.isscalar(hts):
        hts=hts*np.ones(len(lat))
    else:
        hts=hts.flatten()

    mask=np.logical_and(
        np.logical_and(np.greater_equal(lat,latslice[0]),
                       np.less(lat,latslice[1])),
        np.logical_and(np.greater_equal(lon,lonslice[0]),
                       np.less(lon,lonslice[1]))
    )

    LAT=lat[mask]
    LON=lon[mask]
    HT=hts[mask]
    V=Var[mask]

    X,Y,Z=latlonht2xyz(LAT,LON,HT)

    return X,Y,Z,V

#------------------------------------------------------------

def subset_latlonht(lat,lon,hts,var,latslice,lonslice):
    """
    Return (lat,lon,hts,var) points that fall in the
    selected latslice x lonslice bounds (in degrees). 
    NOTE: this fn returns 1d (flattened) x,y,z outputs.
    It is basically just a simpler form of subset_latlonht2xyz.
    """

    #FIRST FLATTEN ALL INPUTS TO AVOID WORRYING ABOUT DIMS WHEN LOOPING
    lat=lat.flatten()
    lon=lon.flatten()
    Var=var.flatten()

    #Making scalar hts same dim as lat and lon simplifies code later
    if np.isscalar(hts):
        hts=hts*np.ones(len(lat))
    else:
        hts=hts.flatten()

    mask=np.logical_and(
        np.logical_and(np.greater_equal(lat,latslice[0]),
                       np.less(lat,latslice[1])),
        np.logical_and(np.greater_equal(lon,lonslice[0]),
                       np.less(lon,lonslice[1]))
    )

    X=lat[mask]
    Y=lon[mask]
    Z=hts[mask]
    V=Var[mask]

    return X,Y,Z,V

#------------------------------------------------------------

def defineBgeoAttrib(outfile, name, valueType, size, defaultValue):
    """Declares an attribute in a .bgeo file."""
    typeInt = 0

    if valueType == 'i':
        typeInt = 1  # int32
    elif valueType == 'x':
	    typeInt = 4  # index
    elif valueType == 'f':
        typeInt = 5   # vector of float32s
        #typeInt = 0  # a single float32. Unless we define 'f' more precisely, can't access this option.

    else:
        raise Exception("Invalid Attribute Type. Please use f (vector of floats), i (int), or x (index)")
        
    outfile.write(struct.pack('>' + 'h', len(name))) # length of name 
    outfile.write(bytearray(name, "utf-8")) # name

    # Index type has its own different rules
    if typeInt == 4:
        outfile.write(struct.pack('>' + 'h', 1))
        outfile.write(struct.pack('>' + 'I', typeInt))
        outfile.write(struct.pack('>' + 'i', size))
        for v in defaultValue:
            outfile.write(struct.pack('>' + 'h', len(v))) # length of default value
            outfile.write(bytearray(list(v)))

    # For all other types
    else:
        outfile.write(struct.pack('>' + 'h', size)) # size
        outfile.write(struct.pack('>' + 'I', typeInt)) # type
        if(size > 1):
            for v in defaultValue:
                outfile.write(struct.pack('>' + valueType, v)) # default value
        else:
            outfile.write(struct.pack('>' + valueType, defaultValue)) # default value

#-------------------------------------------------------------------------------------
            
def writeBgeoFile(FileName,varname,var,X,Y,Z,year=0,dataTime=0):
    """
    Open and write a Bgeo file (using python's text file writing utility)
    .
    INPUTS:
    Filename=string containing name of Bgeo file to open/write to (should end w/ .bgeo)
    varname=string containing name to give this var in the bgeo file.
    var=numpy array containing data to write 
    X=x positions to write (can be latitude or Cartesian x value)
    Y=y positions to write (can be longitude or Cartesian y value)
    Z=z positions to write (height from surface of planet or Cartesian z value)
    dayofyear=day of year (as a float; fractional values are fractions of day). 
        Set to 0 if omitted.
    .
    OUTPUTS:
    none
    .
    NOTES: 
    1). var, X, Y, and Z must all have the same dimensionality
    2). Maybe eventually var and varname can be lists of variables?
    """

    #CHECK FOR ERRORS:
    #=======================
    if var.size!=X.size:
        raise Exception('X should have same dims as var')
    if X.size!=Y.size:
        raise Exception('X and Y should have same dims')
    if X.size!=Z.size:
        raise Exception('X and Z should have same dims')
        
    #GET TOTAL NUMBER OF POINTS TO WRITE
    #=====================
    #This is trivial once we've computed all the (x,y,z) locations
    npoints=np.product(X.size) 
        
    # Write header
    #=====================
    outfile = open(FileName, "wb")
    bgeov = "BgeoV".encode("ascii")
    outfile.write(bgeov)
    outfile.write(struct.pack('>' + 'i', 5)) # version
    nprims = 0
    npointgroups = 0
    nprimgroups = 0
    npointattrib = 1
    nvertexattrib = 0
    nprimattrib = 0
    nattrib = 2
    counts =  [npoints, nprims, npointgroups, nprimgroups, npointattrib, nvertexattrib, nprimattrib, nattrib]
    outfile.write(struct.pack('>'+'i'*len(counts), *counts))
    
    # Define point attributes
    #=======================
    defineBgeoAttrib(outfile, varname, 'f', 1, 0) #open bgeo file, var name, valueType, size, defaultValue
    
    #Write out each point:
    #=======================
    #Stitch x,y,z,1,var for each point together into a zip object that can be iterated over fast
    #It may be more efficient for latlonht2xyz to always return flattened arrays rather than doing here.
    XX=zip(X.flatten(),Y.flatten(),Z.flatten(),np.ones(Z.shape,np.int).flatten(),var.flatten())
    for point in XX:
        outfile.write(struct.pack('>'+'f'*len(point), *point))
                    
    # Write detail attributes
    #=======================
    defineBgeoAttrib(outfile, "year", 'i', 1, 0)
    defineBgeoAttrib(outfile, "day", 'f', 1, 0)  
    outfile.write(struct.pack('>'+'i', year))
    outfile.write(struct.pack('>'+'f', dataTime))

    # EOF
    #=======================
    outfile.write(b"\x00") #beginExtra
    outfile.write(b"\xff") #endExtra

    print('Wrote file')

