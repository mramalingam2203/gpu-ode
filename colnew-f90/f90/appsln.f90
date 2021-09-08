!*==APPSLN.spg  processed by SPAG 6.72Dc at 05:08 on  8 Sep 2021
!
!----------------------------------------------------------------------
!                             p a r t  4
!               polynomial and service routines
!----------------------------------------------------------------------
!
      SUBROUTINE APPSLN(X,Z,Fspace,Ispace)
!
!*****************************************************************
!
!c APPSLN sets up a call to approximate the solution.
!
!     purpose
!
!           set up a standard call to  approx  to evaluate the
!           approximate solution  z = z( u(x) )  at a point x
!           (it has been computed by a call to  colnew ).
!           the parameters needed for  approx  are retrieved
!           from the work arrays  ispace  and  fspace .
!
!*****************************************************************
!
      IMPLICIT NONE
!*--APPSLN25
!*** Start of declarations inserted by SPAG
      REAL*8 a , dummy , Fspace , X , Z
      INTEGER i , is4 , is5 , is6 , Ispace
!*** End of declarations inserted by SPAG
      DIMENSION Z(1) , Fspace(1) , Ispace(1) , a(28) , dummy(1)
      is6 = Ispace(6)
      is5 = Ispace(1) + 2
      is4 = is5 + Ispace(4)*(Ispace(1)+1)
      i = 1
      CALL APPROX(i,X,Z,a,Fspace(is6),Fspace(1),Ispace(1),Fspace(is5),  &
     &            Fspace(is4),Ispace(2),Ispace(3),Ispace(5),Ispace(8),  &
     &            Ispace(4),2,dummy,0)
      END
