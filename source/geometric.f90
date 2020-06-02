MODULE geometric

  ! Sisaltaa kaikki impulssimomenttikytkentaan liittyvat
  ! aliohjelmat.

  PRIVATE

  PUBLIC :: CGCOEF, sixj, ninej, fact, first_fact, dfact, first_dfact, factor, dfactor, rme, moshinsky2, mosh, racah, ninej_int, sixj_int

  INTERFACE CGCOEF
     MODULE PROCEDURE CGCOEF
  END INTERFACE


  INTERFACE sixj
     MODULE PROCEDURE sixj_real,sixj_int
  END INTERFACE

  INTERFACE ninej
     MODULE PROCEDURE ninej_real,ninej_int
  END INTERFACE

  INTEGER, PARAMETER :: bbdim = 1021
  INTEGER, PARAMETER :: factdim = 170
  INTEGER, PARAMETER :: dfactdim = 300
  INTEGER, PARAMETER :: gamdim = 300

  ! CG-kertoimet tallettava vektori.

  REAL(8), ALLOCATABLE, SAVE :: cgv(:,:,:),cgv12(:)

  ! kertoma n!, tuplakertoma n!!, gammafunktio ja binomikerroin.

  REAL(8), SAVE :: fact(0:factdim), dfact(0:dfactdim), gam(0:gamdim), bb(bbdim, bbdim)

  ! Loogiset muuttujat jotka kertovat, onko em. vektorit taytetty

  LOGICAL, SAVE :: first_fact = .TRUE., first_dfact = .TRUE.
  LOGICAL, SAVE :: first_gamma = .TRUE., first_binom = .TRUE.

CONTAINS

  FUNCTION sixj_int(ia,ib,ic,id,ie,IF) RESULT(res)

    ! 6j coefficient, integer arguments

    INTEGER,INTENT(in) :: ia,ib,ic,id,ie,IF
    REAL(8) :: a,b,e,d,c,f,res

    a = 0.5d0*real(ia,kind=8)
    b = 0.5d0*real(ib,kind=8)
    c = 0.5d0*real(ic,kind=8)
    d = 0.5d0*real(id,kind=8)
    e = 0.5d0*real(ie,kind=8)
    f = 0.5d0*real(IF,kind=8)

    res = racah(a,b,e,d,c,f)*(-1d0)**((ia+ib+id+ie)/2)

  END FUNCTION sixj_int

  FUNCTION sixj_real(a,b,c,d,e,f) RESULT(res)

    ! 6j coefficient, real(8) arguments

    REAL(8), INTENT(in) :: a,b,e,d,c,f
    REAL(8) :: res,ai,bi,ei,di,ci,fi

    ai = a; bi = b; ci = c; di = d; ei = e; fi = f
    res = racah(ai,bi,ei,di,ci,fi)*(-1d0)**(NINT(a+b+d+e))

  END FUNCTION sixj_real

  FUNCTION racah(a,b,c,d,e,f)

    ! racah coefficient                                                 

    IMPLICIT REAL (8)(a-h, o-z) 

    IF (first_binom) THEN
       CALL binom
    END IF
    racah = 0.0d0 
    IF (MIN(e-ABS(a-b),e-ABS(c-d),f-ABS(a-c),f-ABS(b-d)) < -0.1) RETURN
    IF (a < 0.0) RETURN
    IF (b < 0.0) RETURN
    IF (c < 0.0) RETURN
    IF (d < 0.0) RETURN
    IF (e < 0.0) RETURN
    IF (f < 0.0) RETURN
    a1 = a+1.01d0 
    d1 = d+1.01d0 
    f1 = f+1.01d0 
    i1 = a1+b-e 
    i2 = c+d1-e 
    i3 = a1+c-f 
    i4 = b+d1-f 
    IF (MIN(i1,i2,i3,i4) < 1) RETURN 
    izn = MAX(0.0d0,a+d-e-f,b+c-e-f)+1.01d0 
    izx = MIN(i1,i2,i3,i4) 
    IF (izn > izx) RETURN 
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
       sum = sum+sn*bb(i1,iz)*bb(j3,m1)*bb(j1,m2) &
            *bb(i4,iz)/bb(i0,iz)                                    
10  END DO
    fac = bb(i0,i1)*bb(i0,i4)*bb(l1,j1)/(k1*k4*bb(k1,ia)*bb(ia,i3) &
         *bb(k4,id)*bb(id,i2)*bb(l1,j4))     
    racah = SQRT(fac)*sum 

    RETURN 

  END FUNCTION racah

  FUNCTION ninej_int(ia,ib,ic,id,ie,if,ig,ih,io) RESULT(coex)

    ! 9-j coefficient                                                   

    integer, intent(in) :: ia,ib,ic,id,ie,if,ig,ih,io
    real(8) :: a,b,c,d,e,f,g,h,o,coex

    a = 0.5d0*real(ia,kind=8)
    b = 0.5d0*real(ib,kind=8)
    c = 0.5d0*real(ic,kind=8)
    d = 0.5d0*real(id,kind=8)
    e = 0.5d0*real(ie,kind=8)
    f = 0.5d0*real(if,kind=8)
    g = 0.5d0*real(ig,kind=8)
    h = 0.5d0*real(ih,kind=8)
    o = 0.5d0*real(io,kind=8)
    coex = ninej_real(a,b,c,d,e,f,g,h,o)

  END FUNCTION ninej_int

  FUNCTION ninej_real(a,b,c,d,e,f,g,h,o) RESULT(coex)

    ! 9-j coefficient                                                   

    IMPLICIT REAL(8) (a-h, o-z)

    IF (first_binom) THEN
       CALL binom
    END IF

    coex = 0.d0 
    x1 = a+b-c 
    x2 = a+d-g 
    x3 = d+e-f 
    x4 = b+h-e 
    x5 = h+o-g 
    x6 = f+o-c 
    IF (MIN(x1,x2,x3,x4,x5,x6,c-ABS(a-b),g-ABS(a-d),f-ABS(d-e),&
         e-ABS(b-h),g-ABS(h-o),c-ABS(f-o)) < -0.1d0) RETURN
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
    sum = 0. 
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
       ff = (2*fk+1)/(jj1*bb(jj1,ja))*SQRT(bb (j1, n1) &
            *bb (m1, l2)*bb (j2, n2)*bb (m2, l4)*bb (j3, n1)        &
            *bb (m3, l1) / (jj2*bb (jh, n3)*bb (m1, l5)*bb (jj2, jb)&
            *bb (jd, n3)*bb (m2, i3)*bb (jf, n2)*bb (m3, l6) ) )    
       ixn = MAX(0.d0, a+h-fk-g, d+o-fk-g)+1.01d0 
       ixx = MIN(n1, n3, i2, i5) 
       sx = 0.d0 
       zx = (-1) **ixn 
       DO ix = ixn, ixx 
          zx =-zx 
          mm1 = n3-ix+1 
          mm2 = i2-ix+1 
          sx = sx+zx*bb (n1, ix)*bb (l2, mm1)*bb (ii1, mm2)      &
               *bb (i5, ix) / bb (j1, ix)  
       END DO
       iyn = MAX(0.d0, b+d-fk-e, f+h-fk-e)+1.01d0 
       iyx = MIN(n2, n3, i4, l3) 
       sy = 0.d0 
       zy = (-1) **iyn 
       DO iy = iyn, iyx 
          zy =-zy 
          mm1 = n3-iy+1 
          mm2 = i4-iy+1 
          sy = sy+zy*bb (n2, iy)*bb (l4, mm1)*bb (ii2, mm2)      &
               *bb (l3, iy) / bb (j2, iy)
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
    ff = bb(j1,i5)*bb(j2,l3)*bb(j3,i6) / (k3*k5*k6*bb &
         (ja, i2)*bb(k5,jh)*bb(jb, i4)*bb(k3,jd)*bb(ja,i1)  &
         *bb (k6, jf) )
    coex = sum*SQRT (ff) 

    RETURN 

  END FUNCTION ninej_real

  !C=======================================================================
  REAL(8) FUNCTION CGCOEF(J1,M1,J2,M2,J,M)

    !C
    !C     CALCULATION OF CLEBSCH-GORDON COEFFICIENTS
    !C     ALSO KNOWN AS VECTOR COUPLING COEFFICIENTS
    !C     DEFINITION OF THE ARGUMENTS ARE STANDARD
    !C     ARGUMENTS ARE DOUBLE THEIR ACTUAL VALUES TO AVOID HALF-INTEGRALS
    !C
    !C   ADAPTED FROM THE CERN LIBRARY ROUTINE CLEBSG (U110),
    !C   WHICH IN TURN WAS ADAPTED FROM THE HARWELL LIBRARY, VERSION 08/01/74
    !C
    !C
    IMPLICIT REAL(8) (A-H,O-Z)
    PARAMETER (MV25=170)
    DIMENSION H(MV25),L(MV25)
    INTEGER, SAVE :: LSIZ=MV25
    !C
    DATA H(1) / 1.D0 /
    DATA L(1) / 0 /
    DATA LAST / 1 /
    DATA SCALE / 8.D0 /
    !C
    CGCOEF=0.0
    IF (J .LT.0)                                 RETURN
    IF (J1.LT.0)                                 RETURN
    IF (J2.LT.0)                                 RETURN
    IF (IABS(M ).GT.J )                          RETURN
    IF (IABS(M1).GT.J1)                          RETURN
    IF (IABS(M2).GT.J2)                          RETURN
    IF (M1+M2.NE.M)                              RETURN
    IF (MOD(J1+J2+J,2).NE.0)                     RETURN
    I0=(J1+J2+J)/2+1
    IF (J1+J2-J.LT.0)                            RETURN
    I1=(J1+J2-J)/2
    IF (J.LT.IABS(J1-J2))                        RETURN
    I2=(J-(J1-J2))/2
    I3=(J+(J1-J2))/2
    IF (MOD(J -M ,2).NE.0)                       RETURN
    I8=(J +M )/2
    I9=(J -M )/2
    IF (MOD(J1-M1,2).NE.0)                       RETURN
    I4=(J1+M1)/2
    I5=(J1-M1)/2
    IF (MOD(J2-M2,2).NE.0)                       RETURN
    I6=(J2+M2)/2
    I7=(J2-M2)/2
    N2=J2-J-M1
    N3=J1-J+M2
    N4=MIN(I1,I5,I6)
    N5=MAX(0,N2,N3)
    N1=N5/2
    IF (N1.GT.N4)                                RETURN
    MM1=MAX(I0,I1,I2,I3,I4,I5,I6,I7,I8,I9,N4-N2/2,N4-N3/2)+1
    IF (MM1.LE.LAST)  GO TO 500
    IF (MM1.LE.LSIZ)  GO TO 100
666 FORMAT("CGCOEF CALLED WITH ARGUMENTS",6I5&
         &,/,8X,"THIS REQUIRES CALCULATIONS WITH A FACTORIAL OF",I5&
         &,/,8X,"THIS VERSION ASSUMED THE MAXIMUM FACTORIAL EVER CA&
         &LCULATE D TO BE OF",I5,/,8X,"THE SIZE OF INTERNAL ARRAYS& 
         &H AND L MUST BE INCREASED")
    WRITE(6, '("CGCOEF CALLED WITH ARGUMENTS",6I5&
         &,/,8X,"THIS REQUIRES CALCULATIONS WITH A FACTORIAL OF",I5&
         &,/,8X,"THIS VERSION ASSUMED THE MAXIMUM FACTORIAL EVER CA&
         &LCULATE D TO BE OF",I5,/,8X,"THE SIZE OF INTERNAL ARRAYS& 
         &H AND L MUST BE INCREASED")') J,M,J1,M1,J2,M2, MM1, LSIZ
    STOP
100 CONTINUE
    X=LAST
    MM2=LAST+1
    LAST=MM1
    DO 400  MM3=MM2,MM1
       H(MM3)=H(MM3-1)*X
       L(MM3)=L(MM3-1)
200    IF (H(MM3).LT.SCALE)  GO TO 400
       H(MM3)=H(MM3)/SCALE**2
       L(MM3)=L(MM3)+2
       GO TO 200
400    X=X+1.D0
500    CONTINUE
       IY =  (L(I1+1)+L(I2+1)+L(I3+1)+L(I4+1)+L(I5+1)+L(I6+1)+L(I7+1)+&
            L(I8+1)+L(I9+1)-L(I0+1))/2
       Y=SQRT(H(I1+1)*H(I2+1)*H(I3+1)*H(I4+1)*H(I5+1)*H(I6+1)*H(I7+1)*&
            H(I8+1)*H(I9+1)/H(I0+1)*(J+1))
       IF ((N5/4)*4-N5+1.LT.0)  Y=-Y
       Z=0.
       DO 510  N5=N1,N4
          MM1=I1-N5
          MM2=I5-N5
          MM3=I6-N5
          MM4=N5-N2/2
          MM5=N5-N3/2
          X=Y/H(MM1+1)/H(MM2+1)/H(MM3+1)/H(MM4+1)/H(MM5+1)/H(N5+1)
          IX= L(MM1+1)+L(MM2+1)+L(MM3+1)+L(MM4+1)+L(MM5+1)+L(N5+1)
800       IF (IY-IX) 900,210,110
900       X=X/SCALE
          IX=IX-1
          GO TO 800
110       X=X*SCALE
          IX=IX+1
          GO TO 800
210       Z=Z+X
          Y=-Y
510    END DO
       CGCOEF=Z
       !write(6,*) "      CGCOEF=",Z
       RETURN
     END FUNCTION CGCOEF

  FUNCTION moshinsky2(lam,n,l,bn,bl,n1,l1,n2,l2,Dfrac) RESULT(sum)

    ! Gives the Moshinsky transformation coefficient <N L n l; lam|n1 l1 n2 l2; lam>.  
    ! This routine uses angular momentum selection rules to give exact 
    ! zeroes when necessary.
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
    ! Dfrac    :  m2/m1 - fraction of the masses of the first and second particle

    IMPLICIT NONE

    INTEGER, INTENT(in) :: lam, n, l, bn, bl, n1, l1, n2, l2
    REAL(8), INTENT(in) :: Dfrac

    INTEGER :: ia,ib,ic,id,k1,k2,p1,p2,t1,t2,phas
    REAL(8) :: fl,fbl,fl1,fl2,flam,sum1,sum2,sum3,sum4,sum
    REAL(8) :: z1,z2

    sum = 0d0
    ia = 2*n1+l1
    ib = 2*n2+l2
    ic = 2*n+l
    id = 2*bn+bl
    flam = REAL(lam,kind=8)
    fl1 = REAL(l1,kind=8)
    fl2 = REAL(l2,kind=8)
    fl = REAL(l,kind=8)
    fbl = REAL(bl,kind=8)

    ! Energy and angular momentum tests

    IF (ia+ib /= ic+id) RETURN
    IF ((lam < ABS (bl-l) ) .OR. (lam > bl+l)) RETURN
    IF ((lam < ABS (l1-l2) ) .OR. (lam > l1+l2)) RETURN
    if (n1 < 0) stop " moshinsky: n1<0 !!!"
    if (n2 < 0) stop " moshinsky: n2<0 !!!"
    if (n < 0) stop " moshinsky: n<0 !!!"
    if (bn < 0) stop " moshinsky: bn<0 !!!"

    ! Initialise factorials etc. if necessary

    IF (first_fact) CALL factor
    IF (first_dfact) CALL dfactor
    IF (first_gamma) CALL gamma

    ! Calculate...

    sum = 0.d0

    sum1=0.5d0*SQRT(ACOS(-1d0))*dble((-1)**(n1+n2+n+lam))*dble(2*l+1)
    sum1=sum1/(1.d0+Dfrac)**(dble(l)/2.d0)
    sum1=sum1*dsqrt(fact(n)*fact(n1)*fact(n2)/fact(bn))
    sum1=sum1*dsqrt(gam(2*n1+2*l1+3)/gam(2*n+2*l+3))
    sum1=sum1*dsqrt(gam(2*n2+2*l2+3)/gam(2*bn+2*bl+3))

    sum3=0.d0
    do k1=0,l
     k2=l-k1
     do p1=abs(l1-k1),l1+k1
      do p2=abs(l2-k2),l2+k2
       sum2=ninej_int(2*p1,2*p2,2*bl,2*k1,2*k2,2*l,2*l1,2*l2,2*lam)
       sum2=sum2*CGCOEF(2*p1,0,2*k1,0,2*l1,0)*CGCOEF(2*p2,0,2*k2,0,2*l2,0)*CGCOEF(2*p1,0,2*p2,0,2*bl,0)
       sum2=sum2*dble((2*p1+1)*(2*p2+1))*dsqrt(fact(2*l)/(fact(2*k1)*fact(2*k2)))
       do t1=0,n1-(k1+p1-l1)/2
        do t2=0,n2-(k2+p2-l2)/2
         sum4=dble((-1)**(k1+t1+t2+(p1+p2+bl)/2))*Dfrac**(dble(t2)+dble(k1+p2)/2.d0)/(1.d0+Dfrac)**(dble(t1+t2)+dble(p1+p2)/2.d0)
         sum4=sum4*gam(p1+p2+bl+2*t1+2*t2+3)/(gam(2*p1+2*t1+3)*gam(2*p2+2*t2+3))
         sum4=sum4/(fact(t1)*fact(t2)*fact(n1-t1-(k1+p1-l1)/2)*fact(n2-t2-(k2+p2-l2)/2))
         if (t1+t2+(p1+p2-bl)/2.lt.0) then
          z1=0.d0
         else
          z1=dble(fact(t1+t2+(p1+p2-bl)/2))
         endif
         if (t1+t2+(p1+p2-bl)/2-bn.lt.0) then
          z2=0.d0
         else
          z2=1.d0/dble(fact(t1+t2+(p1+p2-bl)/2-bn))
         endif
         sum3=sum3+z1*z2*sum2*sum4
        enddo
       enddo
      enddo
     enddo
    enddo

    sum = sum1 * sum3

    RETURN

  END FUNCTION moshinsky2
               
  FUNCTION mosh(lam,n,l,bn,bl,n1,l1,n2,l2) RESULT(sum)

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

    IMPLICIT NONE

    INTEGER(4), INTENT(in) :: lam, n, l, bn, bl, n1, l1, n2, l2 ! To compile with g95/Gillis
!    INTEGER(8), INTENT(in) :: lam, n, l, bn, bl, n1, l1, n2, l2

    INTEGER(4) :: t1,t2
    REAL(8) :: sum1,sum2,sum,q1,q2,q3

    REAL(8), PARAMETER :: pi = 3.141592653589793d0

    IF (first_fact) THEN
       CALL factor
    END IF

    IF (first_gamma) THEN
       CALL gamma
    END IF

    if(l.ne.0) then
     write(*,*) 'Wrong input. (works only for l=0)'
     stop
    endif

    if(lam.ne.bl) then
     write(*,*) 'Wrong input. (lambda=L)'
     stop
    endif

!    if(l1.ne.l2) then
!     write(*,*) 'Wrong input. Angular momentum should be conserved.'
!     stop
!    endif

    if(2*n1+l1+2*n2+l2.ne.2*n+2*bn+bl) then
     write(*,*) 'Geometric: mosh: Wrong input. Selection rule should be valid.'
     write(*,*) 'n1=',n1,'l1=',l1,'n2=',n2,'l2=',l2,'n=',n,'l=',0,'bn=',bn,'bl=',bl
     stop
    endif

    sum1 = 0.5d0 * sqrt(pi) * (-1)**bn * CGCOEF(2*l1,0,2*l2,0,2*bl,0) & 
           * SQRT((dble(fact(n1)*fact(n2)*fact(n)*(2*l1+1)*(2*l2+1))*gam(2*n1+2*l1+3) &
           *gam(2*n2+2*l2+3))/(dble(fact(bn)*(2*bl+1))*gam(2*n+3)*gam(2*bn+2*bl+3)))

    sum2 = 0.d0
     DO t1=0,n1
      DO t2=0,n2
       if(t1+t2+(l1+l2-bl)/2-bn.lt.0) cycle
        q1 = gam(2*t1+2*t2+l1+l2+bl+3)/(gam(2*t1+2*l1+3)*gam(2*t2+2*l2+3))
        q2 = dble(fact(t1+t2+(l1+l2-bl)/2))/dble(fact(t1+t2+(l1+l2-bl)/2-bn))
        q3 = 1.d0 / DBLE(fact(t1)*fact(t2)*fact(n1-t1)*fact(n2-t2))
        sum2 = sum2 + (-1)**(t1+t2) * q1 * q2 * q3 / 2.d0**(dble(t1+t2)+dble(l1+l2)/2.d0)
      ENDDO
     ENDDO

    sum = sum1 * sum2

  RETURN
  END FUNCTION mosh

  FUNCTION rme(n1,l1,n2,l2,m,oscpar) RESULT(res)

    ! Laskee r**m -operaattorin matriisielementin harmonisen 
    ! oskillaattorin tilojen v�lill�. oscpar on oskillaattori-
    ! parametri.

    IMPLICIT NONE

    REAL(8), PARAMETER :: pi = 3.141592653589793d0

    REAL(8), INTENT(in) :: oscpar
    INTEGER, INTENT(in) :: n1,l1,n2,l2,m

    REAL(8) :: xnorm, x, res
    INTEGER :: i,k1,k2
    LOGICAL, SAVE :: first = .TRUE.

    ! Alustetaan vektorit fact, dfact ja gam ensimmaisella
    ! kutsukerralla

    IF (first_fact) THEN
       CALL factor
    END IF
    IF (first_dfact) THEN
       CALL dfactor
    END IF
    IF (first_gamma) THEN
       CALL gamma
    END IF

    ! Matriisielementin arvo lasketaan tassa.

    x = 0d0
    xnorm = SQRT(2d0**(l1+l2-n1-n2+4)*dfact(2*n1+2*l1+1) &
         * dfact(2*n2+2*l2+1)/(fact(n1)*fact(n2)*pi)) &
         / (dfact(2*l1+1)*dfact(2*l2+1))

    DO k1 = 0,n1
       DO k2 = 0,n2
          x = x+gg(n1,k1,l1)*gg(n2,k2,l2)*gam(2*k1+2*k2+l1+l2+m+3)
       END DO
    END DO
    res = 0.5d0*xnorm*x*oscpar**m

  CONTAINS
      
    ! Apufunktio

    FUNCTION gg(n,k,l) RESULT(gres)
      INTEGER :: n,k,l
      REAL(8) :: gres
      gres = (-2.0d0)**k*fact(n)*dfact(2*l+1) &
           /(fact(k)*fact(n-k)*dfact(2*l+2*k+1))
    END FUNCTION gg

  END FUNCTION rme

  SUBROUTINE binom
    ! Binomikerroin.
    REAL(8) :: s,fll,fm1
    DO i = 1,bbdim 
       DO j = 1,bbdim 
          bb(i,j) = 0d0 
       END DO
    END DO
    bb(1,1) = 1d0 
    DO l = 2, bbdim
       bb(l,1) = 1d0 
       bb(l,l) = 1d0 
       IF (l == 2) CYCLE
       lm = (l+1) / 2 
       s = 1d0
       DO m = 2,lm 
          ll = l-m+1 
          fll = REAL(ll,kind=8)
          fm1 = REAL(m-1,kind=8) 
          s = s*fll/fm1 
          bb(l,m) = s 
          bb(l,ll) = s 
       END DO
    END DO
    first_binom = .FALSE.
    RETURN 
  END SUBROUTINE binom

  SUBROUTINE factor
    ! Kertoma n!
    fact(0) = 1d0 
    DO i = 1,factdim
       fact(i) = i*fact(i-1) 
    END DO
    first_fact = .FALSE. 
    RETURN 
  END SUBROUTINE factor

  SUBROUTINE dfactor
    ! Kaksoiskertoma n!!
    dfact(0) = 1d0 
    dfact(1) = 1d0 
    DO i = 2, dfactdim
       dfact(i) = i*dfact(i-2) 
    END DO
    first_dfact = .FALSE. 
    RETURN 
  END SUBROUTINE dfactor
       
  SUBROUTINE gamma
    ! Gammafunktio.
    gam(0) = 1d0 
    gam(1) = SQRT(ACOS(-1d0)) 
    gam(2) = 1d0 
    DO i = 3,gamdim
       gam(i) = (i/2d0-1d0)*gam(i-2) 
    END DO
    first_gamma = .FALSE. 
    RETURN 
  END SUBROUTINE gamma

END MODULE geometric
