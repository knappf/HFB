      MODULE geom

      public :: v_mosh,cleb,sixj_int,sixj_real,racah,
     &                                   ninej_int,ninej_real,calc_fact

      double precision, save :: fac(0:170),gam(0:170),bb(200,200)
      logical, save :: first_fact = .true.

! Strongly recommented to use the subroutine "cleb" for Clebsch-Gordan coefficients
! Do not use cg_int and cg_real - these subroutines gives strange values!!!!

      CONTAINS

      double precision function cleb(j1,m1,j2,m2,j,m)
       implicit none
! calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
! arguments are integer and twice the true value.

       double precision :: factor,sum
       integer :: j1,m1,j2,m2,j,m,par,z,zmin,zmax

       if(first_fact) then
        call calc_fact
        first_fact=.false.
       endif

! some checks for validity (let's just return zero for bogus arguments)

       if(2*(j1/2)-int(2*(j1/2.0)).ne.2*(abs(m1)/2)-int(2*(abs(m1)/2.0))
     &.or.2*(j2/2)-int(2*(j2/2.0)).ne.2*(abs(m2)/2)-int(2*(abs(m2)/2.0))
     &.or.2*(j/2)-int(2*(j/2.0)).ne.2*(abs(m)/2)-int(2*(abs(m)/2.0)) 
     &.or.j1.lt.0.or.j2.lt.0.or.j.lt.0.or.abs(m1).gt.j1.or.abs(m2).gt.j2
     &.or.abs(m).gt.j.or.j1+j2.lt.j.or.abs(j1-j2).gt.j.or.m1+m2.ne.m)
     & then
        cleb= 0.0
       else

        factor=0.0
        factor=bb(1+j1,1+(j1+j2-j)/2)/bb(1+(j1+j2+j+2)/2,1+(j1+j2-j)/2)
        factor=factor * bb(1+j2,1+(j1+j2-j)/2) / bb(1+j1,1+(j1-m1)/2)
        factor=factor / bb(1+j2,1+(j2-m2)/2) / bb(1+j,1+(j-m)/2)
        factor=dsqrt(factor)

        zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
        zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)

        sum=0.0
        do z = zmin,zmax
          par=1
          if(2*(z/2)-int(2*(z/2.0)) .ne. 0) par=-1
          sum=sum+par*bb(1+(j1+j2-j)/2,1+z)
     7            *bb(1+(j1-j2+j)/2,1+(j1-m1)/2-z)*
     &             bb(1+(-j1+j2+j)/2,1+(j2+m2)/2-z)
        end do

        cleb = factor*sum
       end if

       return
      end function cleb

      double precision function v_mosh(lam,n,l,bn,bl,n1,l1,n2,l2) 

! Gives the Moshinsky transformation coefficient <N L n l; lam|n1 l1 n2 l2; lam>.
! for a special case of l=0
!
! Input parameters:
!
! lam      :  Total orbital angular momentum
! n        :  Number of radial nodes in the relative H.O. radial wave function
! l        :  Relative motion orbital angular momentum quantum number
! bn       :  Number of radial nodes in the center-of-mass H.O. radial wave function
! bl       :  center-of-mass motion orbital angular momentum quantum number
! n1, n2   :  Laboratory system wave functions 1 and 2, number of radial nodes
! l1, l2   :  Laboratory system wave functions 1 and 2, orbital angular momentum quantum numbers
!
! Based on J. Phys. A. 13, 1977-1989, (1980)
!
! Petr Vesely 2010

       include 'define.inc'
       include 'parameters.inc'

       INTEGER :: lam, n, l, bn, bl, n1, l1, n2, l2 
       INTEGER :: t1,t2
       double precision :: sum1,sum2,q1,q2,q3
!       double precision :: fac(0:100),gam(0:100)

       if(first_fact) then 
        call calc_fact
        first_fact=.false.
       endif

       if(l.ne.0) then
        write(*,*) 'Wrong input. (works only for l=0)'
        stop
       endif

       if(lam.ne.bl) then
        write(*,*) 'Wrong input. (lambda=L)'
        stop
       endif

       if(2*n1+l1+2*n2+l2.ne.2*n+2*bn+bl) then
!        write(*,*) 'Wrong input. Selection rule should be valid.'
!        write(*,*) 'n1=',n1,'l1=',l1,'n2=',n2,'l2=',l2,'n=',n,'l=',0,
!     &                                               'bn=',bn,'bl=',bl
        v_mosh=0.d0
        return
!        stop
       endif

       sum1 = 0.5d0 * sqrt(pi) * (-1)**bn * cleb(2*l1,0,2*l2,0,2*bl,0)
     &       * SQRT((dble(fac(n1)*fac(n2)*fac(n)*(2*l1+1)*(2*l2+1))
     &       *gam(2*n1+2*l1+3)*gam(2*n2+2*l2+3))/(dble(fac(bn)
     &       *(2*bl+1))*gam(2*n+3)*gam(2*bn+2*bl+3)))

       sum2 = 0.d0
       DO t1=0,n1
        DO t2=0,n2
         if(t1+t2+(l1+l2-bl)/2-bn.ge.0) then
          q1=gam(2*t1+2*t2+l1+l2+bl+3)/(gam(2*t1+2*l1+3)
     &                                             *gam(2*t2+2*l2+3))
          q2=dble(fac(t1+t2+(l1+l2-bl)/2))/dble(fac(t1+t2
     &                                             +(l1+l2-bl)/2-bn))
          q3=1.d0/DBLE(fac(t1)*fac(t2)*fac(n1-t1)*fac(n2-t2))
          sum2=sum2+(-1)**(t1+t2)*q1*q2*q3/2.d0**(dble(t1+t2)
     &                                             +dble(l1+l2)/2.d0)
         endif
        ENDDO
       ENDDO

       v_mosh = sum1 * sum2

       return
      end function v_mosh

      double precision function sixj_int(ia,ib,ic,id,ie,iff)

!     6j coefficient, integer arguments

       INTEGER :: ia,ib,ic,id,ie,iff
       REAL(8) :: a,b,e,d,c,f

       a = 0.5d0*dble(ia)
       b = 0.5d0*dble(ib)
       c = 0.5d0*dble(ic)
       d = 0.5d0*dble(id)
       e = 0.5d0*dble(ie)
       f = 0.5d0*dble(iff)

       sixj_int = racah(a,b,e,d,c,f)*(-1d0)**((ia+ib+id+ie)/2)

      end function sixj_int

      double precision function sixj_real(a,b,c,d,e,f) 

! 6j coefficient, double precision arguments

       double precision :: a,b,e,d,c,f

       sixj_real = racah(a,b,e,d,c,f)*(-1d0)**(NINT(a+b+d+e))

      end function sixj_real

      double precision function ninej_int(ia,ib,ic,id,ie,iff,ig,ih,io) 

! 9-j coefficient

       include 'define.inc'

       integer :: ia,ib,ic,id,ie,iff,ig,ih,io
       double precision :: a,b,c,d,e,f,g,h,o

       a = 0.5d0*dble(ia)
       b = 0.5d0*dble(ib)
       c = 0.5d0*dble(ic)
       d = 0.5d0*dble(id)
       e = 0.5d0*dble(ie)
       f = 0.5d0*dble(iff)
       g = 0.5d0*dble(ig)
       h = 0.5d0*dble(ih)
       o = 0.5d0*dble(io)
       ninej_int = ninej_real(a,b,c,d,e,f,g,h,o)

       return
      end function ninej_int

      double precision function ninej_real(a,b,c,d,e,f,g,h,o)

       include 'define.inc'

       integer :: bbdim
       double precision :: s,fll,fm1
!       double precision :: bb(200,200)

       if(first_fact) then
        call calc_fact
        first_fact=.false.
       endif

       ninej_real = 0.d0
       x1 = a+b-c
       x2 = a+d-g
       x3 = d+e-f
       x4 = b+h-e
       x5 = h+o-g
       x6 = f+o-c
       IF (MIN(x1,x2,x3,x4,x5,x6,c-ABS(a-b),g-ABS(a-d),f-ABS(d-e),
     &    e-ABS(b-h),g-ABS(h-o),c-ABS(f-o)) .lt. -0.1d0) RETURN
       i1 = x1+1.01d0
       i2 = x2+1.01d0
       i3 = x3+1.01d0
       i4 = x4+1.01d0
       i5 = x5+1.01d0
       i6 = x6+1.01d0
       l1 = b+c-a+1.01d0
       l2 = d+g-a+1.01d0
       l3 = d+f-e+1.01d0
       l4 = e+h-b+1.01d0
       l5 = g+h-o+1.01d0
       l6 = c+f-o+1.01d0
       j1 = a+d+h+o+2.01d0
       j2 = b+d+f+h+2.01d0
       j3 = a+b+f+o+2.01d0
       k3 = d+e+f+1.01d0
       k5 = g+h+o+1.01d0
       k6 = c+f+o+1.01d0
       ja = a+a+1.01d0
       jb = b+b+1.01d0
       jd = d+d+1.01d0
       jf = f+f+1.01d0
       jh = h+h+1.01d0
       fkn = MAX(ABS(a-o),ABS(d-h),ABS(b-f))
       fkx = MIN(a+o,d+h,b+f)
       kn = fkx-fkn+1.01d0
       sum = 0.d0
       DO kk = 1, kn
        fk = fkn+kk-1
        n1 = a+o-fk+1.01d0
        n2 = b+f-fk+1.01d0
        n3 = d+h-fk+1.01d0
        m1 = d+g-o+fk+1.01d0
        m2 = e+h-f+fk+1.01d0
        m3 = b+c-o+fk+1.01d0
        ii1 = a-o+fk+1.01d0
        ii2 = b-f+fk+1.01d0
        jj1 = a+o+fk+1.01d0
        jj2 = b+f+fk+1.01d0
        ff = (2*fk+1)/(jj1*bb(jj1,ja))*SQRT(bb (j1, n1) 
     &       *bb (m1, l2)*bb (j2, n2)*bb (m2, l4)*bb (j3, n1)        
     &       *bb (m3, l1) / (jj2*bb (jh, n3)*bb (m1, l5)*bb (jj2, jb)
     &       *bb (jd, n3)*bb (m2, i3)*bb (jf, n2)*bb (m3, l6) ) )
        ixn = MAX(0.d0, a+h-fk-g, d+o-fk-g)+1.01d0
        ixx = MIN(n1, n3, i2, i5)
        sx = 0.d0
        zx = (-1) **ixn
        DO ix = ixn, ixx
          zx =-zx
          mm1 = n3-ix+1
          mm2 = i2-ix+1
          sx = sx+zx*bb (n1, ix)*bb (l2, mm1)*bb (ii1, mm2)      
     &          *bb (i5, ix) / bb (j1, ix)
        END DO
        iyn = MAX(0.d0, b+d-fk-e, f+h-fk-e)+1.01d0
        iyx = MIN(n2, n3, i4, l3)
        sy = 0.d0
        zy = (-1) **iyn
        DO iy = iyn, iyx
          zy =-zy
          mm1 = n3-iy+1
          mm2 = i4-iy+1
          sy = sy+zy*bb (n2, iy)*bb (l4, mm1)*bb (ii2, mm2)      
     &          *bb (l3, iy) / bb (j2, iy)
        END DO
        izn = MAX(0.d0, a+f-fk-c, b+o-fk-c)+1.01d0
        izx = MIN(n1, n2, i1, i6)
        sz = 0.d0
        zz = (-1) **izn
        DO  iz = izn, izx
         zz =-zz
         mm1 = n2-iz+1
         mm2 = i1-iz+1
         sz = sz+zz*bb(n1,iz)*bb(l1,mm1)*bb(ii1,mm2)*bb(i6,iz)/bb(j3,iz)
        END DO
        sum = sum+ff*sx*sy*sz
       END DO
       ff = bb(j1,i5)*bb(j2,l3)*bb(j3,i6) / (k3*k5*k6*bb 
     &    (ja, i2)*bb(k5,jh)*bb(jb, i4)*bb(k3,jd)*bb(ja,i1)  
     &    *bb (k6, jf) )
       ninej_real = sum*SQRT (ff)

       return
      end function ninej_real

      double precision function racah(a,b,c,d,e,f)

! racah coefficient

       IMPLICIT REAL (8)(a-h, o-z)

       integer :: bbdim
       double precision :: s,fll,fm1
!       double precision :: bb(200,200)

       if(first_fact) then
        call calc_fact
        first_fact=.false.
       endif

       racah = 0.0d0
       IF (MIN(e-ABS(a-b),e-ABS(c-d),f-ABS(a-c),f-ABS(b-d)).lt.-0.1) 
     &                    RETURN
       a1 = a+1.01d0
       d1 = d+1.01d0
       f1 = f+1.01d0
       i1 = a1+b-e
       i2 = c+d1-e
       i3 = a1+c-f
       i4 = b+d1-f
       IF (MIN(i1,i2,i3,i4) .lt. 1) RETURN
       izn = MAX(0.0d0,a+d-e-f,b+c-e-f)+1.01d0
       izx = MIN(i1,i2,i3,i4)
       IF (izn .gt. izx) RETURN
       i0 = a1+b+c+d1
       j1 = a1+e-b
       j3 = c+f1-a
       j4 = d1+f-b
       k1 = a1+b+e
       k4 = b+d1+f
       l1 = j1+j3-1
       ia = a1+a
       id = d1+d
       sum = 0.0d0
       sn = (-1.0d0)**izn
       DO 10 iz = izn, izx
        sn =-sn
        izm = iz-1
        m1 = i2-izm
        m2 = i3-izm
        sum = sum+sn*bb(i1,iz)*bb(j3,m1)*bb(j1,m2) 
     &        *bb(i4,iz)/bb(i0,iz)
10     END DO
       fac1 = bb(i0,i1)*bb(i0,i4)*bb(l1,j1)/(k1*k4*bb(k1,ia)*bb(ia,i3) 
     &      *bb(k4,id)*bb(id,i2)*bb(l1,j4))
       racah = SQRT(fac1)*sum

       return
      end function racah

      subroutine calc_fact
       IMPLICIT REAL (8)(a-h, o-z)

       fac(0)=1
       DO i = 1,170
        fac(i) = i*fac(i-1)
       END DO

       gam(0) = 1.d0
       gam(1) = SQRT(pi)
       gam(2) = 1.d0
       DO i = 3,170
        gam(i) = (dble(i)/2.d0-1.d0)*gam(i-2)
       END DO

       bbdim=200
       DO i = 1,bbdim
        DO j = 1,bbdim
          bb(i,j) = 0d0
        END DO
       END DO
       bb(1,1) = 1d0
       DO l = 2, bbdim
        bb(l,1) = 1d0
        bb(l,l) = 1d0
        IF (l .eq. 2) CYCLE
         lm = (l+1) / 2
         s = 1d0
         DO m = 2,lm
          ll = l-m+1
          fll = dble(ll)
          fm1 = dble(m-1)
          s = s*fll/fm1
          bb(l,m) = s
          bb(l,ll) = s
         END DO
       END DO

       return
      end subroutine calc_fact

      END MODULE geom
