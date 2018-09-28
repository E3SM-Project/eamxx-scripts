#include "types.f.h"

module cpp_bridge
  interface
     function array_io_file_exists(filename) result(exists) bind(c)
       use iso_c_binding
       character(kind=c_char), intent(in) :: filename(*)
       logical(kind=c_bool) :: exists
     end function array_io_file_exists

     function array_io_write(filename, a, n) result(ok) bind(c)
       use iso_c_binding
       character(kind=c_char), intent(in) :: filename(*)
       type(c_ptr), intent(in) :: a
       integer(kind=c_int), intent(in), value :: n
       logical(kind=c_bool) :: ok
     end function array_io_write

     function array_io_read(filename, a, n) result(ok) bind(c)
       use iso_c_binding
       character(kind=c_char), intent(in) :: filename(*)
       type(c_ptr) :: a
       integer(kind=c_int), intent(in), value :: n
       logical(kind=c_bool) :: ok
     end function array_io_read

     function dump_all(filename, qr, nr, th, dzq, pres, prt_liq, ni, nk, dt, ts) result(ok) bind(c)
       use iso_c_binding
       character(kind=c_char), intent(in) :: filename(*)
       type(c_ptr), intent(in) :: qr, nr, th, dzq, pres, prt_liq
       integer(kind=c_int), intent(in), value :: ni, nk, ts
       real(kind=c_real), value, intent(in) :: dt
       logical(kind=c_bool) :: ok
     end function dump_all

     subroutine populate_input_from_fortran(nk, kdir, qr, nr, th, dzq, pres) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in) :: qr, nr, th, dzq, pres
       integer(kind=c_int), intent(in), value :: nk, kdir
     end subroutine populate_input_from_fortran

     subroutine dump_arch_f90() bind(c)
       use iso_c_binding
     end subroutine dump_arch_f90

  end interface
end module cpp_bridge
