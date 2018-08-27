program micro_sed

  use micro_sed_mod

  implicit none

  integer :: ni, nk, ts, kdir
  character(len=32) :: arg
  real :: dt

  if (iargc().ne.5) then
     write (*,*) 'Usage: micro_sed ni nk time_step_len num_steps kdir'
     call exit(1)
  end if

  call getarg(1, arg)
  read (arg, *) ni

  call getarg(2, arg)
  read (arg, *) nk

  call getarg(3, arg)
  read (arg, *) dt

  call getarg(4, arg)
  read (arg, *) ts

  call getarg(5, arg)
  read (arg, *) kdir

  if (kdir.ne.-1 .and. kdir.ne.1) then
     write (*,*) 'kdir must be -1 or 1'
     call exit(1)
  end if

  call p3_init()

  call micro_sed_func_wrap(ni, nk, dt, ts, kdir)

end program micro_sed
