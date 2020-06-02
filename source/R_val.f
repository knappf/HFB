      function R_val(n,l,j2,bo,r1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: fc1(0:170),fc2(0:170)

       fc1(0)=1
       DO i = 1,170
        fc1(i) = i*fc1(i-1)
       END DO

       fc2(0)=1
       fc2(1)=1
       DO i = 2,170
        fc2(i) = i*fc2(i-2)
       END DO

       dnorm=dsqrt(2.d0**(dble(n+l+3)))*bo**1.5d0/dsqrt(dsqrt(4.d0*pi))
       dnorm=dnorm*dsqrt(fc1(n)/fc2(2*n+2*l+1))

       i_max=n

       val_L=0.d0
       do i=0,i_max
        if(i.eq.0) then
         a1=1.d0
         a2=1.d0
        endif
        if(i.ne.0) then
         a1=a1*dble(-n+i-1)
         a2=a2*(dble(l+i-1)+1.5d0)
        endif
        val_L=val_L+(bo*r1)**(2.d0*dble(i))*a1/(a2*dble(fc1(i)))
       enddo

       do i=0,n
        if(i.eq.0) a3=1.d0
        if(i.ne.0) a3=a3*(dble(l+i-1)+1.5d0)
       enddo
 
       R_val=dnorm*((bo*r1)**dble(l))*a3*val_L/fc1(n)
       R_val=R_val*dexp(-0.5d0*(bo*r1)**2.d0)

       return
      end
