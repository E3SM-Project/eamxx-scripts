#include "types.f.h"

module shoc

  implicit none

contains

  !=============================================================================!
  subroutine populate_li_input(km1, km2, x1_i, y1_i, x2_i)
    use li_cpp_bridge
    use iso_c_binding

    implicit none

    integer, intent(in) :: km1, km2

    real, dimension(1:km1), intent(inout), target :: x1_i, y1_i
    real, dimension(1:km2), intent(inout), target :: x2_i

    call populate_li_input_from_fortran(km1, km2, c_loc(x1_i), c_loc(y1_i), c_loc(x2_i))

  end subroutine populate_li_input

  !=============================================================================!
  subroutine linear_interp_wrap(km1,km2,ncol,minthresh,repeat) bind(c)
    use li_cpp_bridge
    use iso_c_binding
    use omp_lib

    implicit none

    integer(kind=c_int), intent(in) :: km1, km2, ncol, repeat
    real(kind=c_real), intent(in)   :: minthresh

    real, dimension(:,:), allocatable, target :: x1, y1, x2, y2

    real :: x1_i(km1), y1_i(km1), x2_i(km2), y2_i(km2) ! single col

    integer :: i, r, k
    logical :: ok
    character (kind=c_char, len=*), parameter :: filename = c_char_"fortran"//char(0)
    real(8) :: start, finish

    allocate(x1(ncol, km1))
    allocate(y1(ncol, km1))
    allocate(x2(ncol, km2))
    allocate(y2(ncol, km2))

    ! call dump_arch_f90()
    print '("Running with ncol=",I0," km1=",I0," km2=",I0)', ncol, km1, km2

    call populate_li_input(km1, km2, x1_i, y1_i, x2_i)

    do i = 1, ncol
       do k = 1, km1
          x1(i, k) = x1_i(k)
          y1(i, k) = y1_i(k)
       end do
       do k = 1, km2
          x2(i, k) = x2_i(k)
       end do
    end do

    y2_i(:) = 0

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, k)
    do r = 1, repeat+1

       !$OMP DO
       do i = 1, ncol
          do k = 1, km2
             y2(i, k) = y2_i(k)
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do i = 1, ncol
          call linear_interp(x1(i,:), x2(i,:), y1(i,:), y2(i,:), km1, km2, 1, minthresh)
       end do
       !$OMP END DO

       if (r.eq.1) then
          start = omp_get_wtime()
       endif

    end do
    !$OMP END PARALLEL

    finish = omp_get_wtime()

    print '("Time = ",E20.3," seconds.")', (finish - start) / repeat

    ok = dump_all_li(filename, c_loc(y2), ncol, km1, km2, minthresh)

    deallocate(x1, y1, x2, y2)

  end subroutine linear_interp_wrap

  subroutine linear_interp_c(x1,x2,y1,y2,km1,km2,ncol,minthresh) bind(c)

    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: km1, km2, ncol
    real(kind=c_real), value, intent(in) :: minthresh
    real(kind=c_real), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(kind=c_real), intent(in) :: x2(ncol,km2)
    real(kind=c_real), intent(out) :: y2(ncol,km2)

    call linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)

  end subroutine linear_interp_c

  !==============================================================
  ! Linear interpolation to get values on various grids

  subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)
    implicit none

    integer, intent(in) :: km1, km2
    integer, intent(in) :: ncol
    real, intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real, intent(in) :: x2(ncol,km2)
    real, intent(in) :: minthresh
    real, intent(out) :: y2(ncol,km2)

    integer :: k1, k2, i

    do i=1,ncol
       do k2=1,km2
          if( x2(i,k2) <= x1(i,1) ) then
             y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
          elseif( x2(i,k2) >= x1(i,km1) ) then
             y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
          else
             do k1 = 2,km1
                if( (x2(i,k2)>=x1(i,k1-1)).and.(x2(i,k2)<x1(i,k1)) ) then
                   y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
                endif
             enddo ! end k1 loop
          endif

          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif

       enddo ! end k2 loop
    enddo ! i loop

    return

end subroutine linear_interp

end module shoc
