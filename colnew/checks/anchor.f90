module anchor
implicit none

type params
      real::r, delta, v
end type

type (params)::p1,p2

contains 


subroutine sub_1(par_1L)
	type (params) :: par_1L
	par_1L%r = 2.0*par_1L%r
	par_1L%delta = -(par_1L%delta)**2.0
	par_1L%v = cos(par_1L%v)
	p1%r = 20
	p1%delta = 20
	p1%v = 20

	p2%r = -2
	p2%delta = -3
	p2%v = -5

	write(*,*)p1,p2

end subroutine sub_1

subroutine sub_2(par_2L)
	type (params) :: par_2L
	par_2L%r = 2.0*par_2L%r
	par_2L%delta = -(par_2L%delta)**3.0
	par_2L%v = cos(par_2L%v)
end subroutine sub_2


end module anchor