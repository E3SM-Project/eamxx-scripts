#include "types.f.h"

module li_cpp_bridge
  interface

     function dump_all_li(filename, y2, ncol, km1, km2, minthresh) result(ok) bind(c)
       use iso_c_binding
       character(kind=c_char), intent(in) :: filename(*)
       type(c_ptr), intent(in) :: y2
       integer(kind=c_int), intent(in), value :: ncol, km1, km2
       real(kind=c_real), intent(in), value :: minthresh
       logical(kind=c_bool) :: ok
     end function dump_all_li

     subroutine populate_li_input_from_fortran(km1, km2, x1_i, y1_i, x2_i) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in) :: x1_i, y1_i, x2_i
       integer(kind=c_int), intent(in), value :: km1, km2
     end subroutine populate_li_input_from_fortran

  end interface
end module li_cpp_bridge
