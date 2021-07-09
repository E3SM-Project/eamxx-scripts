#!/usr/bin/env python
"""
Copy of horizontal_interpolate.F90 by Jack Chen from 
CESM/E3SM code ported to python and simplified.

Thoughts:
1. This code computes cell edge lat and lon as ave between consecutive cell centers. This
   obviously won't work if lats and lons aren't monotonically increasing (which physics data 
   aren't). I'm confused how this code works at all. I think it operates on dycore grid variables?
2. It would make more sense to just map the net contribution from a given src cell
   to the corresponding dest cell rather than sweeping in the lat then lon directions 
   as done here. The ESMF/Tempest mapping files already do this, so we should probably
   switch to them. 
2. xy_interp relies on the destination grid being unstructured so the number of 
   dest_lats == the number of dest_lons. The src grid must be lat x lon (with potentially
   different numbers of lats and lons). This is why count_lats and count_lons are potentially 
   different. 
3. weights2indices is essentially a data sparsification routine: it allows us to only loop
   over cells that have nonzero contribution to a given destination cell. While it reduces 
   looping, it *doesn't* reduce the size of the data we need to tote around.

"""

import numpy as np

###############################################################
def xy_interp_init(dest_lats,dest_lons,src_lats,src_lons):
    """
    PURPOSE: 
    Compute the weights to apply to each src data cell to convert it
    into the destination grid. Note that destination lats and lons 
    should be just the ones needed for a given MPI rank, while source 
    lats and lons need to include the global sets of lats and lons 
    because we don't know a priori which source cells will be involved
    in the interpolation.
    
    INPUT VARIABLE DESCRIPTION:
    dest = destination grid (the model in this case)
    src  = source grid (data read in from file in this case)
    lats = cell-center latitudes *in degrees*
    lons = cell-center longitudes *in degrees*

    OUTPUTS:
    weight_lats(dest index,src_index) = weights for interp from dest to src in latitude direction
    weight_lons(dest index,src_index) = weights for interp from dest to src in longitude direction
    """

    #INIT WEIGHTS TO ZERO
    #=====================
    #we will only change weights that are nonzero, so this step is important
    weight_lats=np.zeros([len(dest_lats),len(src_lats)])
    weight_lons=np.zeros([len(dest_lons),len(src_lons)])

    #CHECK THAT CELL-CENTER LAT AND LON FOLLOW INTENDED CONVENTION:
    #=====================
    if np.max(dest_lats)>90:
        raise Exception('Dest lats should be <=90 deg')
    if np.min(dest_lats)<-90.:
        raise Exception('Dest lats should be >=-90 deg')
    if np.max(src_lats)>90:
        raise Exception('Src lats should be <=90 deg')
    if np.min(src_lats)<-90.:
        raise Exception('Src lats should be >=-90 deg')

    if np.max(dest_lons)>360:
        raise Exception('Dest lons should be <=360 deg')
    if np.min(dest_lons)<0.:
        raise Exception('Dest lons should be >=0 deg')
    if np.max(src_lons)>360:
        raise Exception('Src lons should be <=360 deg')
    if np.min(src_lons)<0.:
        raise Exception('Src lons should be >=0 deg')

    #LATS AND LONS NEED TO BE MONOTONICALLY INCREASING AND LINEARLY INCREASING
    #IN ORDER TO BE ABLE TO GET EDGES FROM CELL CENTERS BELOW!
    dest_lats_diff=dest_lats[1:]-dest_lats[0:-1]
    if np.min(dest_lats_diff)<0:
        raise Exception('Dest lats should be monotonically increasing')
    src_lats_diff=src_lats[1:]-src_lats[0:-1]
    if np.min(src_lats_diff)<0:
        raise Exception('Src lats should be monotonically increasing')
    dest_lons_diff=dest_lons[1:]-dest_lons[0:-1]
    if np.min(dest_lons_diff)<0:
        raise Exception('Dest lons should be monotonically increasing')
    src_lons_diff=src_lons[1:]-src_lons[0:-1]
    if np.min(src_lons_diff)<0:
        raise Exception('Src lons should be monotonically increasing')
    
    #CONVERT CELL-CENTER LAT AND LON INTO "CELL-EDGE" LAT AND LON:
    #=====================
    src_lat_edges=np.zeros(len(src_lats)+1)
    #for first and last edge, just extrapolate halfway past the first
    #or last point
    src_lat_edges[0]=src_lats[0] - (src_lats[1]-src_lats[0])/2.
    src_lat_edges[-1]=src_lats[-1] + (src_lats[-1]-src_lats[-2])/2.
    #cell edges in interior of dataset are just the average of adjacent midpts
    for i in range(1,len(src_lats)):
        src_lat_edges[i]=(src_lats[i-1]+src_lats[i])/2.

    src_lon_edges=np.zeros(len(src_lons)+1)
    #for first and last edge, just extrapolate halfway past the first
    #or last point
    src_lon_edges[0]=src_lons[0] - (src_lons[1]-src_lons[0])/2.
    src_lon_edges[-1]=src_lons[-1] + (src_lons[-1]-src_lons[-2])/2.
    #cell edges in interior of dataset are just the average of adjacent midpts
    for i in range(1,len(src_lons)):
        src_lon_edges[i]=(src_lons[i-1]+src_lons[i])/2.

    dest_lat_edges=np.zeros(len(dest_lats)+1)
    #for first and last edge, just extrapolate halfway past the first
    #or last point
    dest_lat_edges[0]=src_lats[0] - (dest_lats[1]-dest_lats[0])/2.
    dest_lat_edges[-1]=dest_lats[-1] + (dest_lats[-1]-dest_lats[-2])/2.
    #cell edges in interior of dataset are just the average of adjacent midpts
    for i in range(1,len(dest_lats)):
        dest_lat_edges[i]=(dest_lats[i-1]+dest_lats[i])/2.
        
    dest_lon_edges=np.zeros(len(dest_lons)+1)
    #for first and last edge, just extrapolate halfway past the first
    #or last point
    dest_lon_edges[0]=src_lons[0] - (dest_lons[1]-dest_lons[0])/2.
    dest_lon_edges[-1]=dest_lons[-1] + (dest_lons[-1]-dest_lons[-2])/2.
    #cell edges in interior of dataset are just the average of adjacent midpts
    for i in range(1,len(dest_lons)):
        dest_lon_edges[i]=(dest_lons[i-1]+dest_lons[i])/2.

    #CESM VERSION ADDS 360 TO ALL LONGITUDES HERE!!!

    #GET WEIGHTS IN LONGITUDE DIRECTION:
    #======================
    #loop over all dest cells and find src cells that contribute to it.
    for dest_i in range(len(dest_lons)): #note dest_i indices are local to MPI rank
        dest_west=dest_lon_edges[dest_i]
        dest_east=dest_lon_edges[dest_i+1]

        for src_i in range(len(src_lons)): #src_i indices span the global src dataset
            src_west=src_lon_edges[src_i]
            src_east=src_lon_edges[src_i+1]

            #For each dest cell, check whether this src cell contributes. If no overlap
            #weighting should be zero (which weights_lon was already initialized to).
            if (src_west>=dest_west) and (src_east<=dest_east):
                # case 1: src contained in dest
                #        src_west             src_east
                #          |-------------------|
                #    |---------------------------------|
                #  dest_west                           dest_east
                weight_lons[dest_i,src_i] =  (src_east-src_west)/(dest_east-dest_west)

            elif (src_west>=dest_west) and (src_west<dest_east):
                # case 2: src extends off end of dest
                #        src_west                          src_east
                #          |--------------------------------|
                #    |---------------------------------|
                #  dest_west                           dest_east
                weight_lons[dest_i,src_i] = (dest_east-src_west)/(dest_east-dest_west)

            elif  (src_east>dest_west) and (src_east<=dest_east):
                # case 3: src extends beyond start of dest
                #       src_west                          src_east
                # |--------------------------------|
                #        |---------------------------------|
                #      dest_west                           dest_east
                weight_lons[dest_i,src_i] = (src_east-dest_west)/(dest_east-dest_west)

    # End points may need a contribution from the other side of the dataset. 
    if src_lon_edges[-1]>dest_lon_edges[-1]:
        # case 1: src wraps back to beginning
        #   src_lon_edges[-2]      src_lon_edges[-1] <--- end point
        #      |-------------------------|
        #    |----------------|......................|
        #dest_lon_edges[-2]  dest_lon_edges[-1]  dest_lon_edges[1] <---using 0 indexing here!
        
        #Trickery below: 1). weight_lons needs to be added to because above loop could have already
        #assigned some weight, 2). (src_lon_edges[-1]-dest_lon_edges[-1]) is a difference of the
        #biggest lon values and (dest_lon_edges[1]-dest_lon_edges[0]) is a difference of the smallest
        #lons. We only avoid the need for modulo 360 arithmetic because taking the difference removes
        #the absolute magnitude of the terms included.
        weight_lons[0,-1]= weight_lons[0,-1]\
            +(src_lon_edges[-1]-dest_lon_edges[-1])/(dest_lon_edges[1]-dest_lon_edges[0])

    if src_lon_edges[-1]<dest_lon_edges[-1]:
        # case 2: dest requires contribution from beginning of src dataset.
        #   src_lon_edges[-2]   src_lon_edges[-1]    src_lon_edges[1]
        #      |-------------------------|.................|
        #           |-------------------------------|
        #    dest_lon_edges[-2]           dest_lon_edges(-1] <--- end point
        
        #Note [-1] indices above are equal to [0] indices modulo 360. Also, the next line
        #is *not* identical to horizontal_interpolate.F90, which I think was wrong.
        weight_lons[-1,0] = weight_lons[-1,0]\
            +(dest_lon_edges[-1]-src_lon_edges[-1])/(dest_lon_edges[-1]-dest_lon_edges[-2]) 

    #GET WEIGHTS IN LATITUDE DIRECTION
    #======================================
    #NOTE: I STILL DON'T UNDERSTAND WHY GAUSSIAN WEIGHTS ARE NEEDED OR WHY DENOMINATORS
    #OF PARTIAL WEIGHTINGS USE SRC RATHER THAN DEST.
    for dest_i in range(len(dest_lats)): #note dest_i indices are local to MPI rank
        dest_south=dest_lat_edges[dest_i]
        dest_north=dest_lat_edges[dest_i+1]
        dest_gw=np.sin(dest_north*np.pi/180.) - np.sin(dest_south*np.pi/180.)
        
        for src_i in range(len(src_lats)): #src_i indices span the global src dataset
            src_south=src_lat_edges[src_i]
            src_north=src_lat_edges[src_i+1]
            src_gw=np.sin(src_north*np.pi/180.) - np.sin(src_south*np.pi/180.)

            #Check if this dest cell gets contribution from the given src cell. If
            #not, do nothing since weight_lats is already initialized to zero.
            #Include Gaussian weighting (gw) in latitude direction.
            if (src_south>=dest_south) and (src_north<=dest_north):
                # case 1: 
                #                src_south             src_north
                #                  |-------------------|
                #            |---------------------------------|
                #          dest_south                           dest_north
                weight_lats[dest_i,src_i] =  src_gw/dest_gw

            elif (src_south>=dest_south) and (src_south<dest_north):
                # case 2: 
                #                src_south                          src_north
                #                  |--------------------------------|
                #            |---------------------------------|
                #          dest_south                           dest_north
                weight_lats[dest_i,src_i] = (dest_north-src_south)/(src_north-src_south)*src_gw/dest_gw

            elif (src_north>dest_south) and (src_north<=dest_north):
                # case 3: 
                #       src_south                          src_north
                #         |--------------------------------|
                #                |---------------------------------|
                #              dest_south                           dest_north
                weight_lats[dest_i,src_i] = (src_north-dest_south)/(src_north-src_south)*src_gw/dest_gw

    #SANITY CHECK RESULTING WEIGHTS:
    #===============================
    if np.max(weight_lons)>1:
        raise Exception('Lon map weights should be <1')
    if np.min(weight_lons)<0:
        raise Exception('Lon map weights should be >0')
    if np.max(weight_lats)>1:
        print('dest_ind,src_ind=%i, %i'%(dest_ind,src_ind))
        raise Exception('Lat map weights should be <1 but max is %f'%(np.max(weight_lats)))
    if np.min(weight_lats)<0:
        raise Exception('Lat map weights should be >0')

    #each dest grid cell should have contributions from src cells which sum to 1
    for dest_i in range(len(dest_lons)):
        if np.abs(np.sum(weight_lons[dest_i,:])-1.)>1e-10:
            raise Exception("Lon weights don't sum to 1 for dest index %i"%(dest_i) )
    for dest_i in range(len(dest_lats)):
        if np.abs(np.sum(weight_lats[dest_i,:])-1.)>1e-10:
            raise Exception("Lat weights don't sum to 1 for dest index %i"%(dest_i) )

    #EXIT, RETURNING NEWLY COMPUTED WEIGHTS
    #===================================
    return weight_lats,weight_lons

###############################################################

def weights2indices(weight_lats,weight_lons):
    """
    Given remapping weights, construct arrays which just 
    tell you the indices of the src domain you need to remap to 
    each destination cell.

    INPUTS:
    weight_lats[dest index,src lat index): how much of a given src lat goes into a given dest cell.
    weight_lons(dest index,src lon index): how much of a given src lon goes into a given dest cell.

    OUTPUTS:
    count_{lats,lons}(dest index): number of src lats or lons with nonzero weight for a given dest cell.
    index_{lats,lons}(dest_index,count_{lats,lons}[dest index]): gives the index for the lat or lon in the 
       original src data with a nonzero weight. This is kind of a matrix sparsifier.

    """

    count_lats=np.zeros(weight_lats.shape[0],np.int)
    count_lons=np.zeros(weight_lons.shape[0],np.int)

    #only the first few src indices will be nonzero, but
    #regridding very fine src to very coarse dest could
    #conceivably require all src cells for a single dest,
    #so making index 2nd dim len(src_lats) just in case.
    index_lats=np.nan*np.ones(weight_lats.shape,np.int) 
    index_lons=np.nan*np.ones(weight_lons.shape,np.int)

    print('weight_lats.shape=',weight_lats.shape)
    print('weight_lons.shpae=',weight_lons.shape)
    
    for dest_i in range(weight_lons.shape[0]):

        #GET LON INDICES
        for src_i in range(weight_lons.shape[1]):
            if weight_lons[dest_i,src_i]>0:
                count_lons[dest_i]+=1
                index_lons[dest_i,count_lons[dest_i]] = src_i

        #GET LAT INDICES
        for src_i in range(weight_lats.shape[1]):
            if weight_lats[dest_i,src_i]>0:
                count_lats[dest_i]+=1
                index_lats[dest_i,count_lats[dest_i]] = src_i

    return count_lats,index_lats,count_lons,index_lons

###############################################################

def xy_interp(weight_lats,count_lats,index_lats,weight_lons,count_lons,index_lons,src_data):
    """
    Remap src_data onto the destination grid using precomputed weights.
    """

    #Next line is a python hack: know number of horiz cells handled here is
    #weight_lons.shape[0]==weight_lats.shape[0] and len(vertical layers) is
    #src_data.shape[-1]
    dest_data=np.zeros([weight_lons.shape[0],src_data.shape[-1]])

    #Loop over dest indices. Uses fact that model data is unstructured rather than lat x lon
    for dest_i in range(weight_lons.shape[0]): 

        for k in range(src_data.shape[-1]): #assuming vertical dim is last
        
            #Loop over the src indices known to have nonzero weight for this dest cell
            for lat_i in range(count_lats[dest_i]):

                #First construct the longitude-interpolated value
                temp_sum=0.
                for lon_i in range(count_lons[dest_i]): 
                    temp_sum += src_data[index_lons[dest_i,lon_i],index_lats[dest_i,lat_i],k]\
                        *weight_lons[dest_i,index_lons[dest_i,lon_i]]

                #Now sum up the contributions in the latitude direction.
                dest_data[dest_i,k] += temp_sum*weight_lats[targ_i,index_lats[targ_i,lat_i]]
    
    return dest_data

###############################################################
