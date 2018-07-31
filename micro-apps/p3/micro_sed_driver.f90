program micro_sed

  use micro_sed_mod

  implicit none

  ! kts, kte, ni, nk, its, ite, dt
  integer, dimension(8) :: args

  integer :: i
  character(len=32) :: arg
  real :: dt

  if (iargc().ne.8) then
     write (*,*) 'Usage: micro_sed kts kte ni nk its ite dt ts'
     call exit(1)
  end if

  do i = 1, 8
     call getarg(i, arg)
     read (arg, *) args(i)
  end do

  call p3_init()

  dt = args(7)
  call micro_sed_func_wrap(args(1), args(2), 1, args(3), args(4), args(5), args(6), dt, args(8))

end program micro_sed
