#!/usr/bin/env python
"""
Implement qv_sat formula from Flatau et al (1992): 
https://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281992%29031%3C1507%3APFTSVP%3E2.0.CO%3B2)
in python to make sure it matches our result in C++.

[BSINGH: Extended the original code to add new test cases along with Murphy and Koop specific
print statements. Test cases are obtained from Murphy and Koop (2005) paper (accessed 07/05/2020):
https://pdfs.semanticscholar.org/519c/bb52a54c5abadff7263bf98678fccce5b594.pdf

Added control variables "print_flatau" and "print_murphyKoop" to decide whether to
print Flatau specific or MurphyKoop specific (or both) print statements]


#Notes: 
1. Units of output from Flatau are never given. They end up being hPa, which can be computed by 
looking at output when T=T0 => es=a1 = 6.12... which is 100x smaller than the 612 Pa we know is the
correct answer.

2. F90 uses MWH2O/MWdry (ratio of molecular weights) as the "0.622" scaling factor for converting
pressure to mixing ratio. Wallace+Hobbes suggest this should be Rd/Rv. Note that Rd=R_star/MWdry and 
Rv=R_star/MWH20 (Wallace and Hobbes eq 2.14), explaining the equivalency. Rounding differences between 
the two versions result in Rd/Rv agreeing with MWH2O/MWdry to ~5 decimal places but disagreeing after. 
I'm sticking with MWH20/MWdry because it is more fundamental and because it is what's already done.

3. Flatau is ambiguous about whether ice saturation should be set equal to liquid saturation at 273.15 K, 
or whether to use the ice saturation equation at that value. In other words, is the upper bound of the 
ice saturation bounds inclusive or exclusive of its endpoint. It would be nice if it was inclusive, because
then we could test that both ice and liquid matched the table 4 (T-T0) intercept values at 273.15 K. Alas, 
F90 and C++ both use liquid values at 273.15 K, so this code does as well.

4. My initial implementation of the 8th degree polynomial used sum(const_i*(T-T0)**i) for i=0:8 while
the F90 and C++ impl used const0+dt*(const1+dt*(...)). This discrepancy resulted in 1e13 errors at 
double precision. Tol is set to 1e-14 so that was unacceptable. I switched to the F90 impl to pass the
test, but I'm grumpy that the methods give such different answers.

5. Getting single-precision to pass is hard because polysvp1 returns a value which is the result of O(20) 
operations, which provides substantial potential for differences relative to the double precision value. 

"""

#IMPORT STUFF:
#=============================
import numpy as np
import pylab as pl

#Setting print_* variables to "True" would print values of that scheme
#currently implemented schemes are "Flatau" and "Murphy and Koop"
print_flatau     = False
print_murphyKoop = True


#DATA USED IN CALCULATIONS:
#=============================
T0=273.15 #"conversion between Celsius an Kelvins" from Flateau. Note this was what was wrong in initial P3 impl.
print('*************************************')
print('********* T0 = '+str(T0)+' **********')
print('*************************************')

Rd=287.042
Rv=461.505
MWH2O= 18.016;
MWdry= 28.966;
ep_2=MWH2O/MWdry
Lv = 2.5e6
Lf = 3.34e5
rho_liq=1000.0
rho_ice=917.0
tboil = 373.16    #from wv_saturation.F90, comments say this val is too high, but most often used.
h2otrip = 273.16  #from shr_const_mod.F90
tmelt = 273.15    #also from shr_const_mod.F90, but a bit confusing. Just seems right.

#values I transcribed directly from the RHS of table 4 in Flatau:
a_liq=np.array(
    [6.11239921,0.443987641,0.142986287e-1,0.26484743e-3,0.302950461e-5,\
     0.206739458e-7,0.640689451e-10,-0.952447341e-13,-0.976195544e-15],np.float64)

a_ice=np.array(
    [6.11147274,0.503160820,0.188439774e-1,0.420895665e-3,0.615021634e-5,\
     0.602588177e-7,0.385852041e-9,0.146898966e-11,0.252751365e-14],np.float64)

#values from original F90 impl (which I've confirmed are identical):
"""
a_liq_F90=np.array(
    [6.11239921,      0.443987641,     0.142986287e-1,\
     0.264847430e-3,  0.302950461e-5,  0.206739458e-7,\
     0.640689451e-10,-0.952447341e-13,-0.976195544e-15],np.float64)

a_ice_F90=np.array(
    [6.11147274,     0.503160820,     0.188439774e-1,\
     0.420895665e-3, 0.615021634e-5,  0.602588177e-7,\
     0.385852041e-9, 0.146898966e-11, 0.252751365e-14],np.float64)

#the next 2 lines check that I copied the coefficients correctly. I have.
print('Diff btwn my liq and F90 ver: ',a_liq - a_liq_F90)
print('Diff btwn my ice and F90 ver: ',a_ice - a_ice_F90)
"""

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def polysvp_liq(T):
    """
    8th order polynomial fit to Wexler's formula for vapor saturation,
    taken from Flatau et al 1992, RHS of table 4. Note that at T=273.15, 
    only the first term is non-zero. This term is ~6.11 but units are unclear.
    #Looking on the internet, sat water vapor at 0 C is 611 Pa... so output
    #from this calculation is in hPa. I convert to Pa before outputting.
    """
    
    """
    esat=a_liq[0]\
          +a_liq[1]*(T-T0)\
          +a_liq[2]*(T-T0)**2.\
          +a_liq[3]*(T-T0)**3.\
          +a_liq[4]*(T-T0)**4.\
          +a_liq[5]*(T-T0)**5.\
          +a_liq[6]*(T-T0)**6.\
          +a_liq[7]*(T-T0)**7.\
          +a_liq[8]*(T-T0)**8.
    """
    dt=T-T0
    esat=a_liq[0]\
          +dt*(a_liq[1]\
               +dt*(a_liq[2]\
                    +dt*(a_liq[3]\
                         +dt*(a_liq[4]
                              +dt*(a_liq[5]\
                                   +dt*(a_liq[6]
                                        +dt*(a_liq[7]
                                             +dt*a_liq[8])))))))
    
    
    return esat*100.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def polysvp_ice_naive(T):
    """
    8th order polynomial fit to Wexler's formula for ice saturation,
    taken from Flatau et al 1992, RHS of table 4. As noted for _liq, 
    output units from the polynomial calc are hPa and I convert to 
    Pa before outputting. I call this version "naive" because it 
    implements Flatau directly, while F90 and C++ use clever grouping
    to avoid using power functions. This difference ends up causing 
    bigger errors, so I use the F90 impl for constructing my test #s. 
    """

    if T>=T0: #at freezing, F90 and C++ versions set ice to liq value.
        return polysvp_liq(T)
    else:

        dt=(T-T0)

        esat=a_ice[0]\
              +a_ice[1]*dt\
              +a_ice[2]*dt**2.\
              +a_ice[3]*dt**3.\
              +a_ice[4]*dt**4.\
              +a_ice[5]*dt**5.\
              +a_ice[6]*dt**6.\
              +a_ice[7]*dt**7.\
              +a_ice[8]*dt**8.
           
        return esat*100.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def polysvp_ice(T):
    """
    Copy F90 formulation of polysvp for ice because my own
    is giving different answers.
    """
    
    if T>=T0: #C++ and F90 assume ice takes liq value at T0.
        return polysvp_liq(T)
    else:
        dt=T-T0
        esat=a_ice[0]\
              +dt*(a_ice[1]\
                   +dt*(a_ice[2]\
                        +dt*(a_ice[3]\
                             +dt*(a_ice[4]
                                  +dt*(a_ice[5]\
                                       +dt*(a_ice[6]
                                            +dt*(a_ice[7]
                                                 +dt*a_ice[8])))))))
    


        return esat*100.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def vap2mix(vap_pres,pres):
    """
    Convert vapor pressure to mixing ratio following eq 2.64 
    from Wallace and Hobbes edition 1. Needs input in Pa,
    which are not the natural output units from polysvp. This
    calculation is exact.
    """

    #scal_fact= Rd/Rv #This seems right to Peter based on Wallace and Hobbes eq 2.64
    scal_fact=MWH2O/MWdry #this is what's actually used in F90
    mix_ratio = scal_fact*vap_pres/(pres - vap_pres)

    return mix_ratio

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def es_clausius_liq(T):
    """
    test my derivation of Clausius-Clapeyron for liquid saturation vapor pressure.
    Ultimately this was just used as a sanity check that my implementation in C++ was correct.
    """

    es=np.zeros(len(T))
    for i in range(len(T)):
        es[i]=polysvp_liq(T[i])

    des_dT=Lv/T * es/(Rv*T - es/rho_liq)

    return des_dT

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def es_clausius_ice(T):
    """
    test my derivation of Clausius-Clapeyron for ice saturation vapor pressure.
    As above, this function doesn't end up providing any numbers for C++, but just
    confirmed that the formulation I use in C++ testing is right.
    """

    es=np.zeros(len(T))
    for i in range(len(T)):
        es[i]=polysvp_ice(T[i])

    des_dT=(Lv+Lf)/T * es/(Rv*T - es/rho_ice)

    return des_dT

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def qs_exact(t_atm,p_atm,is_ice):
    """
    Exact analytical formulation of qs for testing expected 
    error relative to Wexler approximation.
    """


    if (t_atm < T0 and is_ice):
        es_T0=611.147274 
        L=Lv+Lf
    else:
        es_T0=611.239921
        L=Lv
    
    qs_T0 = ep_2 * es_T0 / max(p_atm-es_T0, 1.e-3)
    result = qs_T0*np.exp( -L/Rv*(1/t_atm - 1/T0) );

    return result

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def dqs_dT_exact(T,qs,is_ice):
    """
    Just computes dqs_dT exactly using Clausius Clapeyron
    for use in qs_numeric.
    """

    if (T < T0 and is_ice):
        L=Lv+Lf
    else:
        L=Lv

    dqs_dT = L*qs/(Rv*T**2.)

    return dqs_dT

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def qs_numeric(T_targ,p,is_ice):
    """
    Numerical integration of qs to confirm exact formulation
    is correct. T_targ and p must be scalars. is_ice is a logical.
    """
    dT = 1e-5 #how finely to divide T.

    if (T_targ < T0 and is_ice):
        es_T0=611.147274 
        L=Lv+Lf
    else:
        es_T0=611.239921
        L=Lv
    
    #GET INITIAL CONDITION
    qs = ep_2 * es_T0 / max(p-es_T0, 1.e-3) #always start with qs at T0=273.15
    T = T0

    if T_targ==T0:
        return qs

    elif T_targ < T0:
        while T>T_targ:
            qs -= dT*dqs_dT_exact(T,qs,is_ice)
            T  -= dT
    else:
        while T<T_targ:
            qs += dT*dqs_dT_exact(T,qs,is_ice)
            T  += dT

    return qs

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def GoffGratch_esl(t): 
    """
    For comparison, test Goff Gratch 1946 liquid  used by CAM/E3SM since 2012.
    copied from components/cam/src/physics/cam/wv_sat_methods.F90
    Uncertain below -70 C

    input:  t  = Temperature in Kelvin
    output: es = SVP in Pa
    """

    es = 10.0**(-7.90298*(tboil/t-1.0)+ 
                5.02808*np.log10(tboil/t)- 
                1.3816e-7*(10.**(11.344*(1.-t/tboil))-1.)+ 
                8.1328e-3*(10.**(-3.49149*(tboil/t-1.))-1.)+ 
                np.log10(1013.246))*100.0
    
    return es

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def GoffGratch_esi(t):
    """
    As GoffGratch_esl, but for ice.
    good down to -100 C
    """
    es = 10.**(-9.09718*(h2otrip/t-1.)-3.56654* 
               np.log10(h2otrip/t)+0.876793*(1.-t/h2otrip)+ 
               np.log10(6.1071))*100.
    
    return es

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def MurphyKoop_esl(t):
    """
    As in GoffGratch above, but using Murphy+Koop 2005 formulation.
    good for 123 < T < 332 K
    """

    es = np.exp(54.842763 - (6763.22 / t) - (4.210 * np.log(t)) +
                (0.000367 * t) + (np.tanh(0.0415 * (t - 218.8)) *
                                  (53.878 - (1331.22 / t) - (9.44523 * np.log(t)) +
                                   0.014025 * t)))
    
    return es
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def MurphyKoop_esi(t):
    """
    As above, but for Murphy+Koop 2005 ice. 
    good down to 110 K
    """

    es = np.exp(9.550426 - (5723.265 / t) + (3.53068 * np.log(t))
                - (0.00728332 * t))

    return es

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# END OF FUNCTIONS

#TEST THAT WHEN T=T0=273.15, ESAT = THE CONSTANT TERM FROM THE POLYNOMIAL FIT
#=============================================================================================
#because polynomial expansion is relative to T-T0, all terms but first cancel for T=T0, 
#allowing us to read the expected value directly off the table. Multiplying by 100 b/c 
#table values are in hPa.

print('When T=T0, polysvp should exactly match 1st Flatau coef. Diff = '+str(polysvp_liq(T0) - a_liq[0]*100.))
print('\n')

#unfortunately, polysvp uses liquid saturation at exactly T=T0, so this test is expected to fail.
#print('polysvp_ice(T0) - a_ice[0] = ',polysvp_ice(T0) - a_ice[0]*100.)

#COMPARE MY FLATAU VALUES AGAINST OLD VALUES FROM KYLE AND PRODUCE NEW VALUES FOR C++ TESTS.
#============================================================================================
#info in the following lists was grabbed from p3_unit_tests.cpp before starting my rewrite... so these are Kyle's
#old values when he thought T0=273.16. Description of the terms in the arrays below are:
#[T, pres, old expected sat ice pres, old expected sat liq pres, old expected qsat_ice, old expected qsat_liq]

# Test values @ the melting point of H20 @ 1e5 Pa
tmelt_data=[273.15, 1e5, 610.7960763188032, 610.7960763188032,\
            0.003822318507864685,  0.003822318507864685]

#Test values @ 243.15K @ 1e5 Pa
T243_data=[243.15, 1e5, 37.98530141245404, 50.98455924912173,\
           0.00023634717905493638,  0.0003172707211143376]

#Test values @ 303.15 @ 1e5 Pa
T303_data=[303.15, 1e5, 4242.757341329608, 4242.757341329608,\
           0.0275579183092878, 0.0275579183092878]

#NEW: Test values @ 243.15K @ 500 mb
T243P500_data=[243.15, 5e4, 37.98530141245404, 50.98455924912173,\
           0.00023634717905493638,  0.0003172707211143376]

#Following Test values are directly coming from Murphy and Koop (2005) paper
#Table C1, titled "VALUES RECOMMENDED FOR CHECKING COMPUTER CODES"

#DEFINITIONS:ice vapor pressure(ivp); liquid vapor pressure(lvp)
#"expected values" are from the paper

#Test values @150  @ 1e5 Pa (expected values-> ivp: 6.106e-6; lvp:1.562e-5)
T150_data    = [150,    1e5]

#Test values @180  @ 1e5 Pa (expected values-> ivp: 0.0053975; lvp: 0.011239)
T180_data    = [180,    1e5]

#Test values @210  @ 1e5 Pa (expected values-> ivp: 0.70202; lvp: 1.2335)
T210_data    = [210,    1e5]

#Test values @240  @ 1e5 Pa (expected values-> ivp: 27.272; lvp: 37.667)
T240_data    = [240,    1e5]

#Test values @273.16  @ 1e5 Pa (expected values-> ivp: 611.657; lvp: 611.657)
T273_16_data = [273.16, 1e5]

#Test values @300  @ 1e5 Pa (expected values-> lvp:3536.8, no ivp mentioned as temp is 300K)
T300_data    = [300,    1e5]


#The last line in this loop prints stuff to copy/paste into C++ tests.
#Tell users what that data consists of:
print('Below, "data to copy" is: [sat ice pres, sat liq pres, qsat_ice, qsat_liq]' )


ind=0
for data in [tmelt_data,T243_data,T303_data, T243P500_data, T150_data, T180_data,
             T210_data, T240_data, T273_16_data, T300_data]:
    ind+=1
    print("Case %i: T=%5.2f p=%5.1f mb"%(ind,data[0],data[1]/100.))

    #SATURATION VAPOR PRESSURE CHECKS:
    #---------------------------------
    e_ice_f=polysvp_ice_naive(data[0])
    e_ice=polysvp_ice(data[0])
    e_liq=polysvp_liq(data[0])
    e_liq_GG=GoffGratch_esl(data[0])
    e_ice_GG=GoffGratch_esi(data[0])
    e_liq_MK=MurphyKoop_esl(data[0])
    e_ice_MK=MurphyKoop_esi(data[0])

    #SATURATION MIXING RATIO CHECKS:
    #---------------------------------
    qsat_l=vap2mix(e_liq,data[1])
    qsat_l_GG=vap2mix(e_liq_GG,data[1])
    qsat_l_MK=vap2mix(e_liq_MK,data[1])
    qsat_exact_l=qs_exact(data[0],data[1],is_ice=False)
    qsat_num_l = qs_numeric(data[0],data[1],is_ice=False)
    #if above freezing, don't bother with expensive repeat calculations
    if data[0]<T0:
        qsat_i=vap2mix(e_ice,data[1])
        qsat_i_GG=vap2mix(e_ice_GG,data[1])
        qsat_i_MK=vap2mix(e_ice_MK,data[1])
        qsat_exact_i=qs_exact(data[0],data[1],is_ice=True)
        qsat_num_i = qs_numeric(data[0],data[1],is_ice=True)
    else:
        qsat_i    = qsat_l
        qsat_i_MK = qsat_l_MK

    if(print_flatau):
        print('----------------------------------------------------------')
        print('--                Flatau                                --')
        print('----------------------------------------------------------')
        #print('  Naive minus old expected e_ice='+str(e_ice_f - data[2]))
        #print('  New minus old expected e_ice='+str(e_ice - data[2]))
        print('  F90 ice p vs naive ice p='+str(e_ice - e_ice_f))
        print('  Flatau vs GoffGratch ice p='+str(e_ice_GG - e_ice))
        print('  Flatau vs MurphyKoop ice p='+str(e_liq_MK - e_ice))

        print('  Flatau vs GoffGratch liq p='+str(e_liq_GG - e_liq))
        print('  Flatau vs MurphyKoop liq p='+str(e_liq_MK - e_liq))
        #print('  actual minus old expected qsat_liq='+str(qsat_l - data[5]))
        print('  exact minus numerical qsat_liq='+str(qsat_exact_l - qsat_num_l))
        print('  exact minus Flatau qsat_liq='+str(qsat_exact_l - qsat_l))
        print('  Percent error in Flatau liq = %e'%((qsat_exact_l - qsat_l)/qsat_exact_l*100.) )
        print('  Percent error in GoffGratch liq = %e'%((qsat_exact_l - qsat_l_GG)/qsat_exact_l*100.) )
        print('  Percent error in MurphyKoop liq = %e'%((qsat_exact_l - qsat_l_MK)/qsat_exact_l*100.) )
        #print('  actual minus old expected e_liq='+str(e_liq - data[3]))
        if data[0]<T0:
            #print('  actual minus old expected qsat_ice='+str(qsat_i - data[4]))
            print('  exact minus numerical qsat_ice='+str(qsat_exact_i - qsat_num_i))
            print('  exact minus Flatau qsat_ice='+str(qsat_exact_i - qsat_i))
            print('  Percent error in Flatau ice = %e'%((qsat_exact_i - qsat_i)/qsat_exact_i*100.) )
            print('  Percent error in GoffGratch ice = %e'%((qsat_exact_i - qsat_i_GG)/qsat_exact_i*100.) )
            print('  Percent error in MurphyKoop ice = %e'%((qsat_exact_i - qsat_i_MK)/qsat_exact_i*100.) )
    #PRINT DATA TO COPY
    #---------------------------------
        print('  Data to copy: ',e_ice,e_liq,qsat_i,qsat_l)
    if(print_murphyKoop):
        print('----------------------------------------------------------')
        print('--             Murphy and Koop                          --')
        print('----------------------------------------------------------')
        #print('  F90 ice p vs naive ice p='+str(e_ice - e_ice_f))
        print('  MurphyKoop vs GoffGratch ice p='+str(e_ice_GG - e_ice_MK))
        print('  MurphyKoop vs Flatau     ice p='+str(e_ice    - e_ice_MK))

        print('  MurphyKoop vs GoffGratch liq p='+str(e_liq_GG - e_liq_MK))
        print('  MurphyKoop vs Flatau     liq p='+str(e_liq    - e_liq_MK))
        print('  exact minus numerical qsat_liq='+str(qsat_exact_l - qsat_num_l))
        print('  exact minus MurphyKoop qsat_liq='+str(qsat_exact_l - qsat_l_MK))
        print('  Percent error in MurphyKoop liq = %e'%((qsat_exact_l - qsat_l_MK)/qsat_exact_l*100.) )
        print('  Percent error in GoffGratch liq = %e'%((qsat_exact_l - qsat_l_GG)/qsat_exact_l*100.) )
        print('  Percent error in Flatau     liq = %e'%((qsat_exact_l - qsat_l)/qsat_exact_l*100.) )
        if data[0]<T0:
            #print('  actual minus old expected qsat_ice='+str(qsat_i - data[4]))
            print('  exact minus numerical qsat_ice='+str(qsat_exact_i - qsat_num_i))
            print('  exact minus MurphyKoop qsat_ice='+str(qsat_exact_i - qsat_i_MK))
            print('  Percent error in MurphyKoop ice = %e'%((qsat_exact_i - qsat_i_MK)/qsat_exact_i*100.) )
            print('  Percent error in GoffGratch ice = %e'%((qsat_exact_i - qsat_i_GG)/qsat_exact_i*100.) )
            print('  Percent error in Flatau ice = %e'%((qsat_exact_i - qsat_i)/qsat_exact_i*100.) )
    #PRINT DATA TO COPY
    #---------------------------------
        print('  Data to copy: ',e_ice_MK,e_liq_MK,qsat_i_MK,qsat_l_MK)

#CHECK THAT Pressure CLAUSIUS-CLAPEYRON FORMULATION IS RIGHT
#=============================
#postscript: confirmed my des/dT formulae are right!
"""
dT=10.
Ts=np.arange(-85.,70.,dT)+273.15

esl=np.zeros(len(Ts))
for i in range(len(Ts)):
    esl[i]=polysvp_liq(Ts[i])

desl_dT = (esl[1:]-esl[0:-1])/dT

mid_Ts=(Ts[1:]+Ts[0:-1])/2.
cc_desl_dT = es_clausius_liq(mid_Ts)

pl.figure(1)
pl.plot(mid_Ts,desl_dT,'b-')
pl.plot(mid_Ts,cc_desl_dT,'r--')
pl.title('Liq des/dT')

dT=10.
Ts=np.arange(-90,0.,dT)+273.15

esi=np.zeros(len(Ts))
for i in range(len(Ts)):
    esi[i]=polysvp_ice(Ts[i])

desi_dT = (esi[1:]-esi[0:-1])/dT

mid_Ts=(Ts[1:]+Ts[0:-1])/2.
cc_desi_dT = es_clausius_ice(mid_Ts)

pl.figure(2)
pl.plot(mid_Ts,desi_dT,'b-')
pl.plot(mid_Ts,cc_desi_dT,'r--')
pl.title('Ice des/dT')

pl.show()
"""
