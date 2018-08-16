"""
This module just contains all the variables needed to run p3_main in
a configuration where both liquid and ice processes are exercised.

The "np.float32,order='F'" stuff is needed to make the variables of the
type fortran expects.
"""

#IMPORT NEEDED STUFF:
#======================
import numpy as np

# --- INPUT: DIMENSIONS ---
it= 1  #how many steps the model has taken. If first step, behaves differently to avoid uninitialized vals.
its=1  #lower index for column chunk. always 1
ite=1  #upper index for column chunk. use single column for tests.
kts=1  #1st vertical layer (model top)
kte=72 #top vertical layer (near surface)
ncat=1 #number of ice categories (will always be 1)
n_diag_2d = 1 #number of user-specified 2d diagnostic output fields (we aren't using any, just returns zeros)
n_diag_3d = 1 #number of user-specified 3d diagnostic output fields (we aren't using any, just returns zeros)
log_predictnc = True      # true = 2 moment, 1 = prescribed number conc
typediags_on = True       # compute in-depth precip types for output
model = 'GEM'             # model is used to decide if index=1 is surf. GEM has index=1 top like E3SM.

# --- INOUT VARIABLES ---
qc = np.zeros([ite,kte],np.float32,order='F')       #cloud liquid water mixing ratio    kg/kg
qc[:,-20:-5]=1e-4*(1.-np.arange(15)/14.)            #max cld at ~700mb, decreasing to 0 at 900mb.
nc = np.ones([ite,kte],np.float32,order='F')*1e6    #cloud liquid drop number           #/kg
qr = np.zeros([ite,kte],np.float32,order='F')       #rain water mixing ratio            kg/kg
qr[:,-20:]=1e-5*(1.-np.arange(20)/19.)              #max rain at 700mb, decreasing to zero at surf.
nr = np.ones([ite,kte],np.float32,order='F')*1e6    #rain drop number                   #/kg

#note: indexing for making some levels nonzero assumes there's only 1 column and 1 category.
qitot = np.zeros([ite,kte,ncat],np.float32,order='F')    #total ice mass mixing ratio   kg/kg
qitot[0,-20:-5,0]=1e-4 #*(1.-np.arange(15)/14.)  #max ice at ~700mb, decreasing to 0 at 900mb.
nitot = np.ones([ite,kte,ncat],np.float32,order='F')*1e6 #total ice number              #/kg
qirim = np.zeros([ite,kte,ncat],np.float32,order='F')    #rime ice mass mixing ratio    kg/kg
qirim[0,-20:-5,0]=1e-4*(1.-np.arange(15)/14.)  #max rim at ~700mb, decreasing to 0 at 900mb.
birim = np.zeros([ite,kte,ncat],np.float32,order='F')    #rime ice volume mixing ratio  m3/kg
birim[0,-20:-5,0]=1e-2    #guess at reasonable value based on: m3/kg is 1/density and liquid water has a density of 1000 kg/m3
ssat = np.zeros([ite,kte],np.float32,order='F')          #supersaturation (qv - qs)     kg/kg

#qv goes to zero halfway through profile (to avoid condensate near model top)
tmp=-5e-4+np.reshape(1e-3/float(kte)\
                     *np.arange(kte),[1,kte])
qv = np.where(tmp>0.,tmp,0.,)                            # water vapor mixing ratio  kg kg-1
qv[:,-20:-5]=5e-3                                        #make layer with qc saturated
qv=np.array(qv,np.float32,order='F')

#pres is actually an input variable, but needed here to compute theta.
pres = 100.+np.reshape(1e5/float(kte)\
                  *np.arange(kte),[1,kte])  # pressure (note: must be >0 to avoid nans!)  Pa
pres=np.array(pres,np.float32,order='F')

#To get potential temperature, start by making absolute temperature
#vary between 150K at top of atmos and 300k at surface, then convert to potential temp.
T = (150.+np.reshape(150./float(kte)\
                *np.arange(kte),[1,kte]))
T=np.array(T,np.float32,order='F')
Rd=287.    #Gas constant for dry air                  J/kg/K
cp=1004.   #heat constant of air at constant pressure J/kg
p0=100000. #reference pressure                        Pa
th=T*(p0/pres)**(Rd/cp)                                        # potential temperature           K

#THE NEXT SECTION MODIFIES INOUT VARIABLES TO SATISFY WEIRD CONDITIONS NEEDED FOR CODE COVERAGE
#make qsmall<qitot<1e-8 where T>273.15 in 1 cell to test code instantaneously melting ice.
qitot[0,-1,0]=1.e-9
qv[0,-1]=5.e-2 #also needs to be supersaturated to avoid getting set to 0 earlier.

#make lowest-level qc and qr>0 to trigger surface rain and drizzle calculation
qr[0,-1]=1.e-6
qc[0,-1]=1.e-6

#make qitot>1e-8 where qr=0 to test rain collection conditional
qitot[0,-25,0]=5.e-8

#make qc>0 and qr>0 where T<233.15 to test homogeneous freezing
qc[0,35]=1.e-7
qv[0,35]=1.e-6

#deposition/condensation-freezing needs t<258.15 and >5% supersat
qv[0,33]=1e-4

#ONLY NOW THAT QV IS FINALIZED CAN WE DEFINE OLD VALUES OF IT
th_old = np.array(th,np.float32,order='F')                     # theta @ beginning of timestep   K
qv_old = np.array(qv,np.float32,order='F')                     # qv at beginning of timestep    kg kg-1


# --- INPUT VARIABLES ---
dt =np.float32(1800.)                                          # model time step                s
uzpl = np.zeros([ite,kte],np.float32,order='F')                # vertical air velocity          m s-1

#compute vertical grid spacing dzq (in m) from pres and theta:
g=9.8                                                          #gravity       m/s2
edge_pres=np.zeros([ite,kte+1],np.float32,order='F')           #pressure at cell edges Pa
edge_pres[:,1:-1]=(pres[:,0:-1]+pres[:,1:])/2.
edge_pres[:,0]=max(0.,pres[:,0] - 0.5*(pres[:,1]-pres[:,0])/(1.-0.))
edge_pres[:,-1]=pres[:,-1] + 0.5*(pres[:,-1]-pres[:,-2])/1.0
dpres=edge_pres[:,1:]-edge_pres[:,0:-1]                        #diff in pres btwn cell bottom and top Pa

#note dpres is >0 since index=1 is top of atmos. dz should also be pos, 
#so dropped minus sign in dp/dz=-rho*g
dzq  = Rd*T/(g*pres)*dpres                                     #diff in geometric ht btwn bottom and top m


# --- OUTPUT VARIABLES (included here for reference) ---
#prt_liq   #liquid surface precip          m/s
#prt_sol     #frozen surface precip          m/s
#diag_ze     #equiv reflectivity             dBZ
#diag_effc   #liq effective radius           m
#diag_effi   #ice effective radius           m
#diag_vmi    #mass-weighted ice fall speed   m/s
#diag_di     #mean ice diameter              m
#diag_rhoi   #bulk density of ice            kg/m            
#diag_2d     #user-specified 2d outputs      N/A
#diag_3d     #user-specified 3d outputs      N/A
