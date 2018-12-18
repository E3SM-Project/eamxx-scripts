#ifndef INCLUDE_TYPES_F
#define INCLUDE_TYPES_F

! Use c_real in Fortran code to abstract single/double precision.

#ifdef DOUBLE_PRECISION
#define c_real c_double
#else
#define c_real c_float
#endif

#endif
