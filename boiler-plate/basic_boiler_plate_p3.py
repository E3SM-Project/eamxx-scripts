#!/usr/bin/python

#---------------------------------------------------------------------------------------
# Purpose: A very basic script to output boiler plate for C++ porting. This script
#          doesn't print out the entire boiler plate but it helps with many repetitive
#          parts especially where we have huge argument lists.
#
# Input: All intent ins, in-outs and outs from the fortran code
#
# Author: Balwinder Singh
#---------------------------------------------------------------------------------------

from __future__ import print_function

# Form groups of vars
#---------------------------------------------------------------------------------------
# Break reals into multiple groups, so that we can maintain order in which they appear
#---------------------------------------------------------------------------------------
intent_in_real_1 =  ["t","pres","rho","xxlv","xxls","qvs","qvi"]

intent_in_real_2 = []

intent_in_bool   = []

intent_in_sclr   = []#["dt"]

#For bfb unit test, name of the "data" variable
data_var = 'gtspvd'


#intent-outs
intent_out_real  = ["mu","dv","sc","dqsdt","dqsidt","ab","abi","kap","eii"]


#----------- USER INPUT ENDS (Minor changes may be needed below for special cases)------------------------


# Hopefully all vars will be in the order in which they appear in Fortran
intent_in_real = intent_in_real_1 + intent_in_real_2

ordered_intent_in_vars = intent_in_real_1 + intent_in_bool + intent_in_real_2 + intent_in_sclr
ordered_intent_out_vars = intent_out_real

ordered_vars = intent_in_real_1 + intent_in_bool + intent_in_real_2 + intent_in_sclr + intent_out_real


print (" === p3_functions_f90.cpp ====")
for i, var in enumerate(intent_in_real):
    if i:
        print(', ',end='')
    print('Real '+var+'_',end='')
print('')
for i, var in enumerate(intent_in_bool):
    if i:
        print(', ',end='')
    print('bool '+var+'_',end='')
print('')

for i, var in enumerate(intent_out_real):
    if i:
        print(', ',end='')
    print('Real* '+var+'_',end='')
print('')
print('')
for i, var in enumerate(intent_out_real):
    print('Real local_'+var+' = *'+var+'_;' )


print("typename P3F::Spack ", end='')
for i, var in enumerate(intent_in_real):
    if i:
        print(', ',end='')
    print(var+'('+var+'_)',end='')
print(';')

print()
print()

print("bool ",end='')
for i, var in enumerate(intent_in_bool):
    if i:
        print(', ',end='')
    print(var+'('+var+'_)',end='')
print(';')

print()
print()
print('typename P3F::Scalar',end='')
for i, var in enumerate(intent_in_sclr):
    if i:
        print(', ',end='')
    print(var+'('+var+'_)',end='')
print(';')

print()
print()
print('typename P3F::Spack ',end='')
for i, var in enumerate(intent_out_real):
    if i:
        print(', ',end='')
    print(var+'(local_'+var+')',end='')
print(';')

for i, var in enumerate(intent_out_real):
    print('t_d('+str(i)+') = '+var+'[0];')
print('')

for i, var in enumerate(intent_out_real):
    print('*'+var+'_ = t_h('+str(i)+');')
print('')

for i, var in enumerate(ordered_intent_in_vars):
    if i:
        print(', ',end='')
    print('d.'+var,end='')
print('')

for i, var in enumerate(ordered_intent_out_vars):
    if i:
        print(', ',end='')
    print('&d.'+var,end='')
print('')
print('')

for i, var in enumerate(intent_in_real):
    if i:
        print(', ',end='')
    print('Real '+var,end='')
print('')
for i, var in enumerate(intent_in_bool):
    if i:
        print(', ',end='')
    print('bool '+var,end='')
print('')

for i, var in enumerate(intent_out_real):
    if i:
        print(', ',end='')
    print('Real* '+var,end='')
print('')


print (" === p3_functions_f90.cpp ENDS====")
print("======================================================================")
print (" === p3_functions_f90.hpp ===")

print('Real ',end='')
for i, var in enumerate(intent_in_real):
    if i:
        print(', ',end='')
    print(var,end='')
print('')

print('bool ')
for i, var in enumerate(intent_in_bool):
    if i:
        print(', ',end='')
    print('bool '+var,end='')
print('')

print('//Outs')
print('Real ',end='')
for i, var in enumerate(intent_out_real):
    if i:
        print(', ',end='')
    print(var,end='')
print('')


print (" === p3_functions_f90.hpp ENDS===")


#--------------------------------------------------
# For p3_functions.hpp file subroutine args
#--------------------------------------------------
print()
print()
print('------------------NOT IN ORDER!!!---------------------------------------')
print()
print()

for i, var in enumerate(intent_in_real):
    if i:
        print(', ',end='')
    print('const Spack& '+var,end='')
print(',')
for i, var in enumerate(intent_in_bool):
    if i:
        print(', ',end='')
    print('const bool '+var,end='')
print(',')
for i, var in enumerate(intent_in_sclr):
    if i:
        print(', ',end='')
    print('const Scalar '+var,end='')
print(',')
for i, var in enumerate(intent_out_real):
    if i:
        print(', ',end='')
    print('Spack& '+var,end='')
print(';')


#--------------------------------------------------
# For priting var vals for BFB unit test
#--------------------------------------------------



print('\'{\',', end='')
for i, var in enumerate(ordered_vars):
    if i:
        print(',\',\',',end='')
    print(var,end='')
print(',\'},\'')

#--------------------------------------------------
# For p3_unit_tests.cpp
#--------------------------------------------------
print()
print()
print('Spack ',end='')
for i, var in enumerate(intent_in_real + intent_out_real):
    if i:
        print(', ',end='')
    print(var,end='')
print(';')
print('bool ',end='')
for i, var in enumerate(intent_in_bool):
    if i:
        print(', ',end='')
    print(var,end='')
print(';')

print('Scalar ',end='')
for i, var in enumerate(intent_in_sclr):
    if i:
        print(', ',end='')
    print(var,end='')
print(';')


print()
print()

#single val assignment
for i, var in enumerate(intent_in_bool + intent_in_sclr):
    print(var+ ' = '+data_var+'_device(0).'+var+';')
print()
print()

#intent-ins packs
for i, var in enumerate(intent_in_real + intent_out_real):
     print(var+'[s] = '+data_var+'_device(s).'+var+';')

#----------------------------------------------------------
#copy back results, so opposite of the above
#----------------------------------------------------------
print()
print()

for i, var in enumerate(intent_in_sclr +intent_in_bool):
    print(data_var+'_device(0).'+var+' = '+var+';')
print()
print()

#intent-ins packs (reverse)
for i, var in enumerate(intent_in_real+intent_out_real):
     print(data_var+'_device(s).'+var+' = '+var+'[s];')
#----------------------------------------------------------

#sync back to host to verify output
print()
print()
for i, var in enumerate(intent_out_real):
     print('REQUIRE('+data_var+'[s].'+var+' == '+data_var+'_host(s).'+var+');')
