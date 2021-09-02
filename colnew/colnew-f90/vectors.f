      subroutine dgefa(a,lda,n,ipvt,info)

c*********************************************************************72
c
cc dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72
c
cc dgesl solves a linear system factored by dgefa.
c
c  dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
         end do
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
         end do
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc daxpy computes constant times a vector plus a vector.
c
c  discussion:
c
c    this routine uses double precision real arithmetic.
c
c    this routine uses unrolled loops for increments equal to one.
c
c  modified:
c
c    18 december 2008
c
c  author:
c
c    jack dongarra
c
c  reference:
c
c    jack dongarra, jim bunch, cleve moler, pete stewart,
c    linpack user's guide,
c    siam, 1979,
c    isbn13: 978-0-898711-72-1,
c    lc: qa214.l56.
c
c    charles lawson, richard hanson, david kincaid, fred krogh,
c    basic linear algebra subprograms for fortran usage,
c    acm transactions on mathematical software,
c    volume 5, number 3, pages 308-323, 1979.
c
c  parameters:
c
c    input, integer n, the number of elements in dx and dy.
c
c    input, double precision da, the multiplier of dx.
c
c    input, double precision dx(*), the first vector.
c
c    input, integer incx, the increment between successive entries of dx.
c
c    input/output, double precision dy(*), the second vector.
c    on output, dy(*) has been replaced by dy(*) + da * dx(*).
c
c    input, integer incy, the increment between successive entries of dy.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        ix = 1
        iy = 1
        if ( incx .lt. 0 ) ix = (-n+1)*incx + 1
        if ( incy .lt. 0 ) iy = (-n+1)*incy + 1

        do i = 1, n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      else

        m = mod(n,4)

        do i = 1, m
          dy(i) = dy(i) + da*dx(i)
        end do

        do i = m + 1, n, 4
          dy(i) = dy(i) + da*dx(i)
          dy(i + 1) = dy(i + 1) + da*dx(i + 1)
          dy(i + 2) = dy(i + 2) + da*dx(i + 2)
          dy(i + 3) = dy(i + 3) + da*dx(i + 3)
        end do

      end if

      return
      end
      function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc ddot forms the dot product of two vectors.
c
c  discussion:
c
c    this routine uses double precision real arithmetic.
c
c    this routine uses unrolled loops for increments equal to one.
c
c  modified:
c
c    07 july 2007
c
c  author:
c
c    jack dongarra
c
c  reference:
c
c    jack dongarra, jim bunch, cleve moler, pete stewart,
c    linpack user's guide,
c    siam, 1979,
c    isbn13: 978-0-898711-72-1,
c    lc: qa214.l56.
c
c    charles lawson, richard hanson, david kincaid, fred krogh,
c    basic linear algebra subprograms for fortran usage,
c    acm transactions on mathematical software,
c    volume 5, number 3, pages 308-323, 1979.
c
c  parameters:
c
c    input, integer n, the number of entries in the vectors.
c
c    input, double precision dx(*), the first vector.
c
c    input, integer incx, the increment between successive entries in dx.
c
c    input, double precision dy(*), the second vector.
c
c    input, integer incy, the increment between successive entries in dy.
c
c    output, double precision ddot, the sum of the product of the 
c    corresponding entries of dx and dy.
c
      implicit none

      double precision ddot
      double precision dx(*)
      double precision dy(*)
      double precision dtemp
      integer i,incx,incy,ix,iy,m,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      end do
      if( n .lt. 5 ) go to 60
   40 continue
      do i = m+1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      end do

   60 ddot = dtemp

      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc dscal scales a vector by a constant.
c
c  discussion:
c
c    this routine uses double precision real arithmetic.
c
c  modified:
c
c    07 july 2007
c
c  author:
c
c    jack dongarra
c
c  reference:
c
c    jack dongarra, jim bunch, cleve moler, pete stewart,
c    linpack user's guide,
c    siam, 1979,
c    isbn13: 978-0-898711-72-1,
c    lc: qa214.l56.
c
c    charles lawson, richard hanson, david kincaid, fred krogh,
c    basic linear algebra subprograms for fortran usage,
c    acm transactions on mathematical software,
c    volume 5, number 3, pages 308-323, 1979.
c
c  parameters:
c
c    input, integer n, the number of entries in the vector.
c
c    input, double precision sa, the multiplier.
c
c    input/output, double precision x(*), the vector to be scaled.
c
c    input, integer incx, the increment between successive entries of x.
c
      implicit none

      double precision da
      double precision dx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      end do
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      end do
      if( n .lt. 5 ) return
   40 continue
      do i = m+1, n, 5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      end do

      return
      end
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc idamax finds the index of element having maximum absolute value.
c
c  discussion:
c
c    this routine uses double precision real arithmetic.
c
c  modified:
c
c    07 july 2007
c
c  author:
c
c    jack dongarra
c
c  reference:
c
c    jack dongarra, jim bunch, cleve moler, pete stewart,
c    linpack user's guide,
c    siam, 1979,
c    isbn13: 978-0-898711-72-1,
c    lc: qa214.l56.
c
c    charles lawson, richard hanson, david kincaid, fred krogh,
c    basic linear algebra subprograms for fortran usage,
c    acm transactions on mathematical software,
c    volume 5, number 3, pages 308-323, 1979.
c
c  parameters:
c
c    input, integer n, the number of entries in the vector.
c
c    input, double precision x(*), the vector to be examined.
c
c    input, integer incx, the increment between successive entries of sx.
c
c    output, integer idamax, the index of the element of sx of maximum
c    absolute value.
c
      implicit none

      double precision dx(*),dmax
      integer idamax
      integer i,incx,ix,n

      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do  i = 2,n
         if ( dmax .lt. dabs ( dx(ix) ) ) then
           idamax = i
           dmax = dabs(dx(ix))
         end if
         ix = ix + incx
      end do
      return
c
c  code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do i = 2,n
        if( dmax .lt. dabs(dx(i)) ) then
          idamax = i
          dmax = dabs(dx(i))
        end if
      end do

      return
      end
