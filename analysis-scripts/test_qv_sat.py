#!/usr/bin/env python
"""
Implement qv_sat formula from Flatau et al (1992): 
https://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281992%29031%3C1507%3APFTSVP%3E2.0.CO%3B2)
in python to make sure it matches our result in C++.

#Notes: 
1. Units of output from Flatau are never given. They end up being hPa, which can be computed by 
looking at output when T=T0 => es=a1 = 6.12... which is 100x smaller than the 612 Pa we know is the
correct answer.

2. F90 uses MWH2O/MWdry (ratio of molecular weights) as the "0.622" scaling factor for converting
pressure to mixing ratio. Wallace+Hobbes suggest this should be Rd/Rv. Note that Rd=R_star/MWdry and 
Rv=R_star/MWH20 (Wallace and Hobbes eq 2.14), explaining the equivalency. Rounding differences between 
the two versions result in Rd/Rv agreeing with MWH2O/MWdry to ~5 decimal places but disagreeing after. 
I'm sticking with MWH20/MWdry because it is more fundamental and because it is what's already done.

3. F90 sets ice saturation vapor pressure equal to liquid's value above (but not *at* freezing). 
You can see from table 4 of Flatau that liquid and ice saturation exactly at freezing are not quite equal. 
Interestingly, Kyle's value with T0=273.16 used the liq val at freezing rather than the ice 
value. It seems like this would have induced relatively large errors in his unit test at 0 C??? Flatau 
is ambiguous on whether the end of the esat_ice curve includes or excludes its 0 degree C endpoint or 
not. I'm going to assume it does because it allows us to test that we match their table values
and because that's what the F90 code already does. 

"""

#IMPORT STUFF:
#=============================
import numpy as np
import pylab as pl

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

#values I transcribed directly from the RHS of table 4 in Flatau:
a_liq=np.array(
    [6.11239921,0.443987641,0.142986287e-1,0.26484743e-3,0.302950461e-5,\
     0.206739458e-7,0.640689451e-10,-0.952447341e-13,-0.976195544e-15])

a_ice=np.array(
    [6.11147274,0.503160820,0.188439774e-1,0.420895665e-3,0.615021634e-5,\
     0.602588177e-7,0.385852041e-9,0.146898966e-11,0.252751365e-14])

#values from original F90 impl:
a_liq_F90=np.array(
    [6.11239921,      0.443987641,     0.142986287e-1,\
     0.264847430e-3,  0.302950461e-5,  0.206739458e-7,\
     0.640689451e-10,-0.952447341e-13,-0.976195544e-15])

a_ice_F90=np.array(
    [6.11147274,     0.503160820,     0.188439774e-1,\
     0.420895665e-3, 0.615021634e-5,  0.602588177e-7,\
     0.385852041e-9, 0.146898966e-11, 0.252751365e-14])

#the next 2 lines check that I copied the coefficients correctly. I have.
#print('Diff btwn my liq and F90 ver: ',a_liq - a_liq_F90)
#print('Diff btwn my ice and F90 ver: ',a_ice - a_ice_F90)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def polysvp_liq(T):
    """
    8th order polynomial fit to Wexler's formula for vapor saturation,
    taken from Flatau et al 1992, RHS of table 4. Note that at T=273.15, 
    only the first term is non-zero. This term is ~6.11 but units are unclear.
    #Looking on the internet, sat water vapor at 0 C is 611 Pa... so output
    #from this calculation is in hPa. I convert to Pa before outputting.
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

    return esat*100.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def polysvp_ice(T):
    """
    8th order polynomial fit to Wexler's formula for ice saturation,
    taken from Flatau et al 1992, RHS of table 4. As noted for _liq, 
    output units from the polynomial calc are hPa and I convert to 
    Pa before outputting.
    """

    if T>T0: #using strictly greater than is consistent w/ F90 ver.
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
def polysvp_F90_ice(T):
    """
    Copy F90 formulation of polysvp for ice because my own
    is giving different answers.
    """
    
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
    which are not the natural output units from polysvp.
    """

    #scal_fact= Rd/Rv #This seems right to Peter based on Wallace and Hobbes eq 2.64
    scal_fact=MWH2O/MWdry #this is what's actually used in F90
    mix_ratio = scal_fact*vap_pres/(pres - vap_pres)

    return mix_ratio

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#info in the following lists grabbed from p3_unit_tests.cpp:
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

#The last line in this loop prints stuff to copy/paste into C++ tests.
#Tell users what that data consists of:
print('Last lin in loop gives: expected sat ice pres, expected sat liq pres, expected qsat_ice, expected qsat_liq')


ind=0
for data in [tmelt_data,T243_data,T303_data]:
    ind+=1
    print('Case '+str(ind))

    #e_ice_f=polysvp_F90_ice(data[0])
    #print('  F90 actual minus old expected e_ice='+str(e_ice_f - data[2]))

    e_ice=polysvp_ice(data[0])
    print('  my actual minus old expected e_ice='+str(e_ice - data[2]))

    e_liq=polysvp_liq(data[0])
    print('  actual minus old expected e_liq='+str(e_liq - data[3]))

    qsat_i=vap2mix(e_ice,data[1])
    print('  actual minus old expected qsat_ice='+str(qsat_i - data[4]))

    qsat_l=vap2mix(e_liq,data[1])
    print('  actual minus old expected qsat_liq='+str(qsat_l - data[5]))

    print('  Data to copy: ',e_ice,e_liq,qsat_i,qsat_l)

#===============================
#WHEN T=T0, ESAT SHOULD EQUAL THE FIRST TERM OF ITS POLYNOMIAL FIT
#BECAUSE ALL OTHER TERMS ARE MULTIPLIED BY (T-T0)

print('polysvp_liq(T0) - a_liq[0] = ',polysvp_liq(T0) - a_liq[0]*100.)
print('polysvp_ice(T0) - a_ice[0] = ',polysvp_ice(T0) - a_ice[0]*100.)
#print('polysvp_F90_ice(T0) - a_ice[0] = ',polysvp_F90_ice(T0) - a_ice[0]*100.)
