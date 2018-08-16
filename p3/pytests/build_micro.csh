#!/usr/bin/env csh

#This script uses f2py to compile the p3 microphysics code and 
#sets up the necessary linkage to be able to call it from python.
# for details about f2py, see https://acme-climate.atlassian.net/
#wiki/spaces/NGDNA/pages/764739624/Making+Standalone+Component+
#Models+using+F2PY

#GCOV notes: implementation here follows: https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/765263889/gcov+for+code+coverage. gcov is already defined on nersc so step 1 on that page not needed. 

#BASIC SET UP:
#=================
# NM will be the name of the module. In python, subroutines will be 
# accessible under <NM>.<filename>.<fortran subroutine name>.

set NM='p3'

#This is the source code directory.
#ON NERSC:
#set DIR='/global/u1/p/petercal/E3SM_code/scream/components/cam/src/physics/cam'
#ON LC:
set DIR='/g/g11/caldwep/gitwork/scream/components/cam/src/physics/cam'

#These are the files to load:
set FILES="${DIR}/micro_p3.F90"

#NOTE: The line below defines all the functions/subroutines we want python 
#access to. Each of these functions must be defined as 'public' in the 
#fortran code we are parsing. Note also that the functions/subroutines below 
#are separated by spaces: the code still compiles if commas are used instead,
#but it crashes when you try to load it into python. 

set PUB_FNS='p3_init p3_main'

#CREATE .pyf FILE WITH INTERFACES TO PYTHON:
#=================
#The .pyf file is a fortran file which just contains the interfaces needed
#to link the PUB_FNS to python. Note that the first line below is needed f2py stuff and the 2nd line is gcov stuff.

f2py ${FILES} -m ${NM} --overwrite-signature -h ${NM}.pyf only: $PUB_FNS 

#FIX PROBLEMS WITH AUTO-GENERATED .pyf FILE :
#=================
#sed -ie 's/kind=r8/kind=8/g' ${NM}.pyf

#NOW ACTUALLY COMPILE THE CODE
#=================
f2py -c ${NM}.pyf ${FILES} --fcompiler=gfortran --f90flags="-fprofile-arcs -ftest-coverage" -lgcov --build-dir .
