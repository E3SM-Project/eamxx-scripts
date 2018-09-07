#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools for plotting atmospheric output from the E3SM model, 
which is on a particular unstructured grid.
"""

#IMPORT STUFF
#===================
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import os
from tqdm import tqdm
import PIL

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=
class E3SM_Map(object):
    def __init__(self,lat,lon):
        """
        This object sets up the (lat,lon) map to use for all 
        E3SM_Map.egg(data) calls.
        """
        self.lat = lat
        self.wrapwidth = 90.0*np.sqrt(6.0/(len(self.lat)-2))*3
        self.wraplat = [np.argwhere(self.lat>=90.0-self.wrapwidth)]
        self.wraplat.append(np.argwhere(self.lat<=-90+self.wrapwidth))
        # Longitude needs to be converted for basemap
        self.lon = lon
        tmp = np.where(self.lon>180.0)[0]
        self.lon[tmp] = self.lon[tmp]-360.0
        self.wraplon = [np.argwhere(self.lon>=180-self.wrapwidth)]
        self.wraplon.append(np.argwhere(self.lon<=-180+self.wrapwidth))
        #
        # Setup default values
        self.defaults = { 'outdir' : os.getcwd() \
                         ,'figsize' : [20,14] \
                         ,'projection' : 'cyl' \
                         ,'fontsize' : 14 \
                         ,'draw_parallels' : True \
                         ,'parallel_spacing' : 15 \
                         ,'draw_meridians' : True \
                         ,'meridian_spacing' : 30 \
                         ,'draw_coasts' : True \
                         ,'map_border' : [-180,180,-90,90] \
                         ,'savefig' : True \
                         ,'figname' : 'basemap'\
                         ,'cmap' : cm.seismic \
                         ,'ncontours' : 11 \
                         ,'dpi' : 200 \
                         ,'showplot' : True \
                         ,'clabel' : '' \
                         ,'title' : '' \
                         }
        
#------------------------------------------------------------
    def list_defaults(self):
        print "\nCurrent Default Values:\n"
        for it in self.defaults:
            mystr = "%20s:\t" %(it)
            print mystr, self.defaults[it]

#------------------------------------------------------------            
    def set_default(self, **kwargs):
        for it in kwargs:
            if not self.defaults.get(it) == None:
                self.defaults[it] = kwargs.get(it)
            else:
                print "Error: ARG = %s not included in default list" %(it)
                print "Issue 'list_defaults' to see a list of valid options"
                
#------------------------------------------------------------
    def list_kwargs(self,func):
        myd = { 'valid options' : ['global_plot'] }
        if func == 'global_plot':
           myd = { \
                 'savefig' : 'LOGICAL: whether or not to save the figure as png' \
                 ,'border' : 'REAL(dim=4): lon/lat border of plot (left,right,bottom,top)' \
                 ,'figname' : 'STR: output figure name' \
                 ,'figsize' : 'REAL(dim=2): over ride figure size (X,Y)' \
                 ,'cmap' : 'matplotlib cmap to use for color scheme' \
                 ,'ncon' : 'INT: number of contours to use in plot' \
                 ,'dpi' : 'INT: dots per inch in saved figure' \
                 ,'data' : 'List of REAL Arrays: Use this data instead of loading from self' \
                 ,'units' : 'STR: use these units for plot labels' \
                 ,'lattick' : 'REAL: spacing for latitude ticks' \
                 ,'lontick' : 'REAL: spacing for longitude ticks' \
                 ,'showplot' : 'LOGICAL: True to display figure' \
                 ,'clabel' : 'STR: label for colorbar axis' \
                 ,'title' : 'STR: title for figure' \
                 ,'vmin'  : 'REAL: minimum contour value' \
                 ,'vmax'  : 'REAL: maximum contour value' \
                 }
        for it in myd:
            mystr = "%20s:\t" %(it)
            print mystr, myd[it]
      

#------------------------------------------------------------
    def egg(self,data,**kwargs):
        fs = self.defaults.get('fontsize')
        border = kwargs.get('border',self.defaults.get('map_border'))
        savefig = kwargs.get('savefig',self.defaults.get('savefig'))
        showplot = kwargs.get('showplot',self.defaults.get('showplot'))
        mlat = self.lat
        mlon = self.lon
        mfigsz = kwargs.get('figsize',self.defaults.get('figsize'))
        cmap = kwargs.get('cmap',self.defaults.get('cmap'))
        ncon = kwargs.get('ncontours',self.defaults.get('ncontours'))
        dpi = kwargs.get('dpi',self.defaults.get('dpi'))
        clabel = kwargs.get('clabel',self.defaults.get('clabel'))
        title  = kwargs.get('title',self.defaults.get('title'))
        stations = kwargs.get('stations',[])
        stnlabel = kwargs.get('station_labels',[])
        cbar_format = kwargs.get('cbar_format','default')

        ## Set up wrapped data for full map plot
        mlat = np.append(mlat,mlat[self.wraplon[0]])
        mlon = np.append(mlon,mlon[self.wraplon[0]]-360.0)
        wrap = self.wraplon[0]
        mlat = np.append(mlat,mlat[self.wraplon[1]])
        mlon = np.append(mlon,mlon[self.wraplon[1]]+360.0)
        wrap = np.append(wrap,self.wraplon[1])
        mlat = np.append(mlat,mlat[self.wraplat[0]]-180.0)
        mlon = np.append(mlon,mlon[self.wraplat[0]])
        wrap = np.append(wrap,self.wraplat[0])
        mlat = np.append(mlat,mlat[self.wraplat[1]]+180.0)
        mlon = np.append(mlon,mlon[self.wraplat[1]])
        wrap = np.append(wrap,self.wraplat[1])
        
        # Setup figure
        fig = pl.figure(figsize=mfigsz)
        mymap = Basemap(projection=self.defaults.get('projection')\
                        ,resolution='c' \
                        ,lat_0=np.amax([np.amin(mlat),-90.0])\
                        ,lon_0=0\
                        ,llcrnrlon=border[0]\
                        ,urcrnrlon=border[1]\
                        ,llcrnrlat=border[2]\
                        ,urcrnrlat=border[3]\
                    )
        xx,yy = mymap(mlon,mlat)
        if self.defaults.get('draw_parallels'):
            latspace = kwargs.get('lattick',self.defaults.get('parallel_spacing'))
            mymap.drawparallels(np.arange(np.floor(border[2]),np.ceil(border[3])+1,latspace),\
                                labels=[1,0,0,0],fontsize=fs)
        if self.defaults.get('draw_meridians'):
            lonspace = kwargs.get('lontick',self.defaults.get('meridian_spacing'))
            mymap.drawmeridians(np.arange(np.floor(border[0]),np.ceil(border[1])+1,lonspace),\
                                labels=[0,0,0,1],fontsize=fs)
        if self.defaults.get('draw_coasts'):
            mymap.drawcoastlines(linewidth=0.25)
        ## Draw Contour
        if len(data)>0: 
        	vmin = kwargs.get('vmin',np.amin(data))
	        vmax = kwargs.get('vmax',np.amax(data))
        
	        im1 = mymap.contourf(lon,lat,data,np.linspace(vmin,vmax,ncon)\
                             ,tri=True\
                             ,cmap=cmap\
                             )
                if cbar_format == 'default':
		        cb1 = fig.colorbar(im1,orientation='horizontal',shrink=0.5) #,cax = cbaxes
		else:
		        cb1 = fig.colorbar(im1,orientation='horizontal',shrink=0.5,format=cbar_format) #,cax = cbaxes
        	cb1.set_label(clabel,fontsize=fs) 
	else:
		mymap.bluemarble()
	if len(stations)>0:
		mymap.plot(stations[0],stations[1],'o',color='red',linewidth=6,markerfacecolor='red')
		for si, ss in enumerate(stnlabel):
			plt.text(stations[0][si]-2,stations[1][si]+3,ss,color='red',fontsize=14)


                    
        pl.title(title,fontsize=fs)
            
        if savefig:
            figfile=kwargs.get('figname','egg.png')
            pl.savefig(figfile,dpi=dpi,bbox_inches='tight')
        if showplot:
            pl.show()
        pl.close(fig) 

        
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=
#STUFF FOR TESTING SCRIPT:
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=

from netCDF4 import Dataset

#f=Dataset('/global/cscratch1/sd/petercal/ACME_simulations/theta.20180216.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG/run/theta.20180216.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG.cam.h0.0001-12.nc')

f = []
f.append(Dataset("/p/cscratch/acme/donahue5/E3SM_Runs/parallel_split/20180305/ftype30/output/FQcase1/case_scripts.cam.h1.0001-01-01-00000.nc"))
f.append(Dataset("/p/cscratch/acme/donahue5/E3SM_Runs/parallel_split/20180305/ftype30/output/ftype2/ftype2_split_ps_ss.cam.h1.0001-01-01-00000.nc"))

lat=f[0].variables['lat'][:]
lon=f[0].variables['lon'][:]

a=E3SM_Map(lat,lon)


