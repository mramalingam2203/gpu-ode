program main
use anchor

implicit none


type (params) :: par_1,par_2

par_1%r = 1.0001
par_1%delta = -1.0991
par_1%v = 0.1919

par_2%r = 1.0212
par_2%delta = -1.210
par_2%v = 0.32122


call sub_1(par_1)
call sub_2(par_2)

! write(*,*)par_1
! write(*,*)par_2



end program main
