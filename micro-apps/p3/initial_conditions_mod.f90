module initial_conditions_mod
  interface
     subroutine fully_populate_input_data(ni, nk, qr, nr, th, dzq, pres) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in) :: qr, nr, th, dzq, pres
       integer(kind=c_int), intent(in), value :: ni, nk
     end subroutine fully_populate_input_data
  end interface
end module initial_conditions_mod
