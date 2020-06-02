      function B_val(l,r1,p1) ! the spherical Bessel function  sqrt(2/pi) * j_l(k*r)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision,allocatable :: bb(:)

       allocate(bb(0:l))
       bb=0.d0

       dnorm=dsqrt(2.d0/pi)
 
       val=0.d0

       z=max(1.d-7,r1*p1)

       if(l.le.1) then
       if(l.eq.0) then
        val=dsin(z)/z
       endif

       if(l.eq.1) then
        val=dsin(z)/z**2.d0-dcos(z)/z
       endif
       endif

       if(l.gt.1) then                                                  !we use the recursive formula j_{l-1}(z)+j_{l+1}(z)=((2l+1)/z)*j_{l}(z)
        bb(0)=dsin(z)/z
        bb(1)=dsin(z)/z**2.d0-dcos(z)/z
        do m=2,l
         vval=0.d0
         vval=-bb(m-2)+bb(m-1)*dble(2*m-1)/z
         bb(m)=vval
        enddo
        val=bb(l)
       endif

       B_val=dnorm*val

       deallocate(bb)

       return
      end
