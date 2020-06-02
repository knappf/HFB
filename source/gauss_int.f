      double precision function gauss_int(lam,ni,li,nj,lj,b1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: fc1(0:170),gm1(0:170)
       double precision :: fac2(0:170),gi(0:170)

       fc1(0)=1
       do i = 1,170
        fc1(i) = i*fc1(i-1)
       enddo

       gm1(0) = 1.d0
       gm1(1) = dsqrt(pi)
       gm1(2) = 1.d0
       do i = 3,170
        gm1(i) = (dble(i)/2.d0-1.d0)*gm1(i-2)
       enddo

       fac2(0)=1.d0
       fac2(1)=1.d0
       do i=2,170
        fac2(i)=dble(i)*fac2(i-2)
       enddo

       gi(0)=dsqrt(pi)/2.d0
       gi(1)=0.5d0
       do i=2,170
        if(mod(i,2).eq.0) gi(i)=dsqrt(pi)*fac2(i-1)/dble(2**(1+i/2))
        if(mod(i,2).eq.1) gi(i)=fac2(i-1)/dble(2**(1+(i-1)/2))
       enddo

       gauss_int=0.d0

       dd1=4.d0/(dsqrt(pi)*b1**lam)
       dd2=dsqrt(fc1(ni)*fc1(nj)*dble(2**(ni+nj+li+lj))
     &                     /(fac2(2*ni+2*li+1)*fac2(2*nj+2*lj+1)))
       ss1=0.d0
       do mi=0,ni
        do mj=0,nj
         ss1=ss1+dble((-1)**(mi+mj))*gm1(2*ni+2*li+3)*gm1(2*nj+2*lj+3)
     &            *gi(lam+2+li+lj+2*mi+2*mj)/(fc1(mi)*fc1(mj)*fc1(ni-mi)
     &                    *fc1(nj-mj)*gm1(2*li+2*mi+3)*gm1(2*lj+2*mj+3))
        enddo
       enddo

       gauss_int=dd1*dd2*ss1

       return
      end
