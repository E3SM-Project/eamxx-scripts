#!/usr/bin/env python
"""
Tools for doing model analysis.
"""

#import cdms2 as cdms
from netCDF4 import Dataset
from netCDF4 import MFDataset
import numpy as np
import os.path
import time

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def get_fis(dir,run_name,start_yr,end_yr):

    fis=[]
    times=[]
    i=-1
    for yr in range(start_yr,end_yr+1): #include endpoint
        i+=1
        #print 'yr = ',yr
        j=-1
        for mon in range(1,12+1):
            j+=1
            #print '  mon = ',mon
            fi=dir+'/'+run_name+'.cam.h0.'+leftpad(str(yr),4)\
                    +'-'+leftpad(str(mon),2)+'.nc'
            if os.path.isfile(fi):
                fis.append(fi)
                times.append(yr+(mon-1)/12.)
            else:
                print "  Couldn't find "+fi
                
                #if first year, allow partial year when time_res='mo' (checked elsewhere)
                if i==0: #allow partial year if 1st year:
                    if j==0:
                        return fis,times

                #quit looping once found a file that doesn't exist... 
                #assume files after that don't exist either. 
                #This is useful when a simulation hasn't finished 
                #a complete year
                elif j!=0: #if only part of a year is available
                    for k in range(j):
                        nm=fis.pop()
                        times.pop()
                        #print 'k=%i, removing %s'%(k,nm) #for debugging
                    return fis,times 
                else: #if j==0, which means this is the 1st month of the year.
                    return fis,times #don't need to strip other months, but do need to exit.

    return fis,times

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def global_ave_timeseries(run_name,vars,dir,start_yr,end_yr,time_res='mon'):
    """
    Open history files in dir/run_name.cam.h0.* with
    times between start_yr and end_yr (including all 
    months in both start_yr and end_yr) and compute the global-average
    of each variable in the list vars (where each entry in vars 
    is a string containing the name of a 2d variable in h0 files).
    time_res can be mon (for monthly output; default) or ann (for annual output).
    """

    #GET LIST OF FILES TO LOAD:
    #==========================
    ttime=time.time()
    btime=time.time()

    fis,times=get_fis(dir,run_name,start_yr,end_yr)

    if fis==[]:
        raise Exception('No usable files found')

    if time_res=='ann' and len(fis)%12!=0:
        raise Exception('Incomplete year was returned. time_res="ann" requires full year')

    print '  Got file list. Took %f sec.'%(time.time() - btime)

    #OPEN ALL FILES AS IF THEY WERE A SINGLE FILE:
    #===========================
    btime=time.time()
    try:
        f = MFDataset(fis)
    except:
        print 'fis = ',fis
        raise Exception("Files exist but Couldn't be opened by MFDataset.")

    print '  Opened files. Took %f sec.'%(time.time() - btime)

    #GET AREA WEIGHTS/MAKE APPROPRIATE SIZE
    #============================
    btime=time.time()
    area=f.variables['area'][:].squeeze() #get area weights for glob ave

    #Next line relies on area being 1-D (e.g. for SE dycore) not 2-D (e.g. FV dycore)
    AREA=np.dot(np.ones([len(fis),1]),np.reshape(area,[1,len(area)]))

    print '  Got weights. Took %f sec.'%(time.time() - btime)

    #GET NINO3.4 WEIGHTS
    #============================
    btime=time.time()

    lats=f.variables['lat'][:].squeeze() #pos for deg N, neg for deg S
    lons=f.variables['lon'][:].squeeze() #deg E, ranges btwn 0 and 360

    nino34_area=np.where( np.logical_and( lats>=-5.0,lats<=5.0),area,0.) #zero weight outside of 10 deg band around equator.
    nino34_area=np.where( np.logical_and( lons>=210.0,lons<=240.0),nino34_area,0.) #zero weight outside Pacific longitudes.

    NINO34_AREA=np.dot(np.ones([len(fis),1]),np.reshape(nino34_area,[1,len(area)]))

    print '  Got Nino3.4 weights. Took %f sec.'%(time.time() - btime)

    #LOOP OVER VARIABLES AND GET GLOB AVE FOR EACH TIME:
    #==========================================
    btime=time.time()

    y=np.zeros([len(vars),len(fis)])
    k=-1
    for var in vars:
        k+=1

        print '    var = ',var

        if var.lower()=='nino3.4':
            x=f.variables['TS'][:].squeeze()
            y[k,:]=np.sum(x*NINO34_AREA,axis=1)/np.sum(nino34_area)

        else:
            x=f.variables[var][:].squeeze()

            #the next line relies on x and area having dims=[time,lat/lon index] and will
            #probably fail with a message about dimensions if this
            #isn't true.

            y[k,:]=np.sum(x*AREA,axis=1)/np.sum(area)

    f.close()

    print '  Computed global aves. Took %f sec.'%(time.time() - btime)

    #MAKE ANNUAL AVERAGES IF REQUESTED:
    #==================================
    if time_res=='mon':
        OUT=y
    elif time_res=='ann':
        #GET FRACTION OF A YEAR FOR EACH MONTH (assuming noleap calendar):
        dpm=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        frac_yr=dpm/365.

        OUT=y[:,0::12]*frac_yr[0] #first month of each year
        times=times[0::12] #just take 1st month of year
        for i in range(1,12):
            OUT+=y[:,i::12]*frac_yr[i]
    else:
        raise Exception('time_res can only be mon or ann')

    print '  total time in global_ave_timeseries = %f'%(time.time()-ttime)

    return OUT,times

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def leftpad(str,num=2):
    """
    Takes a string as input, and if it has length < num, left-pad w/ zeros until
    its size is correct. For example, leftpad('2',4) returns '0002'.
    """
    if len(str)<num:
        str='0'+str
        str=leftpad(str,num) #recursion here
    return str
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def linreg(x,y):
    """
    [m,b]=linreg(x,y) gives the least-squares linear regression of x
    onto y in the form Y=m*X+b.Converted from cdat/MV version 5/19/17.

    Only works for 1-d x and y. Behavior with missing values or nans unclear.
    """
    xp=x-np.average(x)
    yp=y-np.average(y)

    m=np.average(np.multiply(xp,yp))/np.average(xp**2)
    b=np.average(y)-m*np.average(x)

    return m,b
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def rmse(var,run_name,dir,start_yr,end_yr):
    """
    return spatial RMSE of the temporal average of "var" for run "run_name" with
    h0 files in "dir".
    """

    from scipy.interpolate import griddata
    from pylab import meshgrid


    #FIRST, COMPUTE TEMPORAL AVERAGE:
    #======================================

    #GET LIST OF FILES TO LOAD:
    #--------------------------
    ttime=time.time()
    btime=time.time()

    fis,times=get_fis(dir,run_name,start_yr,end_yr)

    if fis==[]:
        raise Exception('No usable files found')

    print '  Got file list. Took %f sec.'%(time.time() - btime)

    #OPEN ALL FILES AS IF THEY WERE A SINGLE FILE:
    #--------------------------
    btime=time.time()
    try:
        f = MFDataset(fis)
    except:
        print 'fis = ',fis
        raise Exception("Files exist but Couldn't be opened by MFDataset.")

    if var=='PRECT':
        x=f.variables['PRECL'][:].squeeze()+f.variables['PRECC'][:].squeeze()
        #x=x*86400.*1000. #convert m/s to mm/day for comparison against GPCP below
        x=x*2.5e6*1000. #converts m/s into W/m2 (need to do this for obs below too!)

    else:
        x=f.variables[var][:].squeeze()

    lats=f.variables['lat'][:].squeeze()
    lons=f.variables['lon'][:].squeeze()

    print '  Opened files, got var. Took %f sec.'%(time.time() - btime)

    #COMPUTE TEMPORAL AVERAGE:
    #--------------------------
    btime=time.time()

    #GET FRACTION OF A YEAR FOR EACH MONTH (assuming noleap calendar):
    dpm=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    frac_yr=dpm/365.

    #Note: dims=(time,latlonind)
    y=x[0::12,:]*frac_yr[0] #first month of each year
    for i in range(1,12):
        y+=x[i::12,:]*frac_yr[i]

    z=np.average(y,axis=0).squeeze() #if more than 1 yrs, compute ave over years.

    print '  Computed temporal ave. Took %f sec.'%(time.time() - btime)


    #NEXT, ACQUIRE OBS DATA:
    #======================================
    btime=time.time()

    if var=='PRECT':
        fo=MFDataset('/global/cscratch1/sd/petercal/ACME_diags/obs_data_20140804/GPCP_ANN_climo.nc')
        z_obs=fo.variables['PRECT'][:].squeeze()
        lats_obs=fo.variables['lat'][:].squeeze()
        lons_obs=fo.variables['lon'][:].squeeze()
        fo.close()

        z_obs = z_obs/86400./1000.*2.5e6*1000. #convert mm/day->m/sec, then-> W/m2

    elif var in ['SWCF','LWCF']:
        fo=MFDataset('/global/cscratch1/sd/petercal/ACME_diags/obs_data_20140804/CERES-EBAF_ANN_climo.nc')
        z_obs=fo.variables[var][:].squeeze()
        lats_obs=fo.variables['lat'][:].squeeze()
        lons_obs=fo.variables['lon'][:].squeeze()
        fo.close()
    else:
        raise Exception('var can only be PRECT, SWCF, or LWCF!')

    print '  Loaded obs data. Took %f sec.'%(time.time() - btime)


    #NOW REGRID MODEL DATA ONTO OBS GRID:
    #======================================
    #following https://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
    #checked this works by pcolor plotting obs prect vs regridded model prect and 
    #noting that they look suitably similar.

    btime=time.time()

    orig_pts=np.array([lons,lats]).T
    LONS_OBS,LATS_OBS=meshgrid(lons_obs,lats_obs)

    z_rgr = griddata( orig_pts, z, (LONS_OBS,LATS_OBS), method='linear')

    if np.sum(np.isnan(z_rgr).astype(np.int))>0:
        print '  ** have %i nans in time-ave model output after regridding'\
            %(np.sum(np.isnan(z_rgr).astype(np.int)))
        print '  *** switching to nearest neighbor.'

        z_rgr = griddata( orig_pts, z, (LONS_OBS,LATS_OBS), method='nearest')


    print '  Regridded model data. Took %f sec.'%(time.time() - btime)

    """FOR DEBUGGING
    import pylab as pl
    pl.pcolor(LONS_OBS,LATS_OBS,z_rgr)
    pl.colorbar()
    pl.show()
    """

    #NOW COMPUTE RMSE:
    #======================================
    #note: using latitude weighting. Checked yrs 11-20 of 
    #20170512.beta1_05-v0.4atm-60Locn.A_WCYCL1850S.ne30_oECv3_ICG.edison 
    #against AMWG. I had 1.02 and AMWG had 0.998 for PRECT, 13.79 vs 13.428 
    #for SWCF, and 5.90 vs 5.80 for LWCF.

    btime=time.time()

    rmse = np.sqrt(np.sum( np.cos(LATS_OBS*np.pi/180.)*(z_rgr - z_obs)**2. )/np.sum(np.cos(LATS_OBS*np.pi/180.)))

    print '  total time in rmse = %f sec'%(time.time()-ttime)

    return rmse
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def timing_data(f):
    """
    Return the data loaded from a timing file (f) in a structure suitable for
    organization into a tree.
    """
    data = None
    sta = "CPL:INIT"      # First line of timing tree
    fin = "Overhead sum"  # Last line of timing tree
    with open(f,'r') as a:
        found = False
        for line in a:
            if '*' in line:  # remove * from line but retain proper level
                line = line[1:]
                line = ' ' + line
            if found:
                #if fin in line or len(line) <= 1: return data
                if line.find('"') == lstrip or len(line) <= 1: return data
                data.append(line[lstrip:])
            else:
                if sta in line:
                    found = True
                    data = list()
                    lstrip = line.find('"')
                    data.append(line[lstrip:])
    return data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

class Timing_Node(object):
    """
    This class stores each line of timing information as a node which can
    then be used with the Timing_Tree class.
    """
    def __init__(self):
        self.name = None
        self.parent = None
        self.data = None
        self.children = []
        self.level = None
    
    def __str__(self):
        print self.name
    
    def printNode(self):
        if self.name:
            tmp = ''
            for i in range(self.level):
                tmp = tmp + ' '
            tmp = tmp + self.name
            print tmp
        if not self.children:
            return
        
    def printChildren(self):
        if self.children:
            for child in self.children:
                child.printNode()
        return
    
    def getdata(self,unit):
        if unit == 'wallclock':
            return self.wallclock()
        elif unit == 'minclock':
            return self.minclock()
        elif unit == 'maxclock':
            return self.maxclock()
        elif unit == 'numcalls':
            return self.numcalls()
        elif unit == 'meanclock':
            return self.wallclock()/self.numcalls()
        else:
            return None
    
    def wallclock(self):
        return self.data[1]
    
    def minclock(self):
        return self.data[3]
    
    def maxclock(self):
        return self.data[2]
    
    def numcalls(self):
        return self.data[0]
    
    def checkchildren(self,datatype,a,b,chk):
        a = self.getdata(datatype)
        b = 0.0
        chk = True
        for child in self.children:
            b = b+child.getdata(datatype)
        if b>a:
            chk = False
        return a,b,chk
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

class Timing_Tree(object):
    """
    This class stores all timing nodes into a tree stucture.
    """
    def __init__(self):
        self.root = None
        self.treename = None
        self.maxlevel = 0
        self.npes = None
        self.dpes = None
        self.ppes = None
        self.ncol = None
        self.nele = None

    def istree(self):
        return True
    
    def addNode(self, _node, _current):
        self.maxlevel = max(self.maxlevel,_node.level)
        if _node.level == 0:
            self.root = _node
            return
        if _node.level == _current.level + 1:
            _current.children.append(_node)
            _node.parent = _current
            return
        if _node.level == _current.level:
            _current.parent.children.append(_node)
            _node.parent = _current.parent
            return
        if _node.level < _current.level:
            _node.parent = self.findParent(_node, _current)
            _node.parent.children.append(_node)
            return
        
    def findParent(self, _node, _current):
        if _node.level == _current.level:
            return _current.parent
        else:
            return self.findParent(_node, _current.parent)
    
    def printTree(self, current):
        current.printNode()
        if not current.children:
            return
        for child in current.children:
            self.printTree(child)
        
    def cumulsum(self,current,cs,datatype):
        if not current.parent:
            return cs
        for child in current.parent.children:
            if child.name == current.name:
                break
            cs = cs + child.getdata(datatype)
        cs = self.cumulsum(current.parent,cs,datatype)
        return cs
        
    def findNodesAtLevel(self, level, current, nodeList):
        if current.level == level:
            nodeList.append(current)
            return nodeList
        if level > current.level:
            for n in current.children:
                self.findNodesAtLevel(level, n, nodeList)
            return
        if level < current.level:
            return nodeList
        
    def findNodesByName(self, name, current, nodeList):
        if name in current.name:
            nodeList.append(current)
        for n in current.children:
            self.findNodesByName(name, n, nodeList)
        return nodeList
        
    def findNode(self, _node, _name,**kwargs):
        if 'exact' in kwargs:
            exact = kwargs.get('exact')
        else:
            exact = False
        if _name in _node.name:
            if exact and _name == _node.name:
                return _node
            elif not exact:
                return _node
        if not _node.children:
            return None
        for child in _node.children:
            target = self.findNode(child, _name,**kwargs)
            if target is not None:
                break
        return target

    def checkTree(self, current,datatype):
        a,b,chk = current.checkchildren(datatype,0.0,0.0,True)
        if not chk:
            print "%s FAILS with (parent,children) = (%f,%f)" %(current.name,a,b)
        if not current.children:
            return
        for child in current.children:
            self.checkTree(child,datatype)
