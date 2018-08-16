"""
This module just contains all the variables needed to run p3_main in
a configuration where no microphysical processes should be active. 
In other words, all condensate species should be zero at the end of this 
call and theta and qv should be unchanged.
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
qc = np.zeros([ite,kte],np.float32)       #cloud liquid water mixing ratio    kg/kg
nc = np.ones([ite,kte],np.float32)*1e6    #cloud liquid drop number           #/kg
qr = np.zeros([ite,kte],np.float32)       #rain water mixing ratio            kg/kg
nr = np.ones([ite,kte],np.float32)*1e6    #rain drop number                   #/kg

qitot = np.zeros([ite,kte,ncat],np.float32)    #total ice mass mixing ratio   kg/kg
nitot = np.ones([ite,kte,ncat],np.float32)*1e6 #total ice number              #/kg
qirim = np.zeros([ite,kte,ncat],np.float32)    #rime ice mass mixing ratio    kg/kg
birim = np.zeros([ite,kte,ncat],np.float32)    #rime ice volume mixing ratio  m3/kg
ssat = np.zeros([ite,kte],np.float32)          #supersaturation (qv - qs)     kg/kg

#qv goes to zero halfway through profile (to avoid condensate near model top)
tmp=-5e-4+np.reshape(1e-3/float(kte)\
                     *np.arange(kte),[1,kte]).astype(np.float32)
qv = np.where(tmp>0.,tmp,0.,)                                 # water vapor mixing ratio  kg kg-1
#th = np.ones([ite,kte],np.float32)*300.                      #potential temperature          K

#potential temperature should increase with height. This simple profile is missing a stratosphere.
th = (500.-np.reshape(200./float(kte)\
                *np.arange(kte),[1,kte])).astype(np.float32) # potential temperature           K
th_old = th.copy()                                           # theta @ beginning of timestep   K
qv_old = qv.copy()                                           # qv at beginning of timestep    kg kg-1

# --- INPUT VARIABLES ---
dt =np.float32(1800.)                                          # model time step                s
uzpl = np.zeros([ite,kte],np.float32)                          # vertical air velocity          m s-1
pres = 1.+np.reshape(1e5/float(kte)\
                  *np.arange(kte),[1,kte]).astype(np.float32)  # pressure (note: must be >0 to avoid nans!)  Pa

#compute vertical grid spacing dzq (in m) from pres and theta:
Rd=287.    #Gas constant for dry air                  J/kg/K
cp=1004.   #heat constant of air at constant pressure J/kg
p0=100000. #reference pressure                        Pa
g=9.8      #gravity                                   m/s2
edge_pres=np.zeros([ite,kte+1],np.float32)
edge_pres[:,1:-1]=(pres[:,0:-1]+pres[:,1:])/2.
edge_pres[:,0]=max(0.,pres[:,0] - 0.5*(pres[:,1]-pres[:,0])/(1.-0.))
edge_pres[:,-1]=pres[:,-1] + 0.5*(pres[:,-1]-pres[:,-2])/1.0
dpres=edge_pres[:,1:]-edge_pres[:,0:-1]

#convert potential temp to temp.
T=th*(pres/p0)**(Rd/cp)

#note dpres is >0 since index=1 is top of atmos. dz should also be pos, so dropped minus sign in dp/dz=-rho*g
dzq  = Rd*T/(g*pres)*dpres  


# --- OUTPUT VARIABLES (included for reference) ---
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
