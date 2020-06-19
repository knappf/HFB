      double precision function VNNtz(i,j,k,l,Jp)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       VNNtz=0.d0

       if(levtz(i)%tz.eq.-1.and.levtz(j)%tz.eq.-1) then
        if(levtz(k)%tz.eq.-1.and.levtz(l)%tz.eq.-1) then
         VNNtz=Vpp(levtz(i)%point,levtz(j)%point,
     &                            levtz(k)%point,levtz(l)%point,Jp)
        endif
       endif

       if(levtz(i)%tz.eq.1.and.levtz(j)%tz.eq.1) then
        if(levtz(k)%tz.eq.1.and.levtz(l)%tz.eq.1) then
         VNNtz=Vnn(levtz(i)%point,levtz(j)%point,
     &                            levtz(k)%point,levtz(l)%point,Jp)
        endif
       endif

       if(levtz(i)%tz.eq.-1.and.levtz(j)%tz.eq.1) then
        if(levtz(k)%tz.eq.-1.and.levtz(l)%tz.eq.1) then
         VNNtz=Vpn(levtz(i)%point,levtz(j)%point,
     &                            levtz(k)%point,levtz(l)%point,Jp)
        endif
       endif
       if(levtz(i)%tz.eq.1.and.levtz(j)%tz.eq.-1) then
        if(levtz(k)%tz.eq.1.and.levtz(l)%tz.eq.-1) then
         VNNtz=Vpn(levtz(j)%point,levtz(i)%point,
     &                            levtz(l)%point,levtz(k)%point,Jp)
     &        *dble((-1)**((levtz(i)%j2+levtz(j)%j2+levtz(k)%j2+
     &                                             levtz(l)%j2)/2))
        endif
       endif
       if(levtz(i)%tz.eq.-1.and.levtz(j)%tz.eq.1) then
        if(levtz(k)%tz.eq.1.and.levtz(l)%tz.eq.-1) then
         VNNtz=-Vpn(levtz(i)%point,levtz(j)%point,
     &                            levtz(l)%point,levtz(k)%point,Jp)
     &        *dble((-1)**(Jp+(levtz(k)%j2+levtz(l)%j2)/2))
        endif
       endif
       if(levtz(i)%tz.eq.1.and.levtz(j)%tz.eq.-1) then
        if(levtz(k)%tz.eq.-1.and.levtz(l)%tz.eq.1) then
         VNNtz=-Vpn(levtz(j)%point,levtz(i)%point,
     &                            levtz(k)%point,levtz(l)%point,Jp)
     &        *dble((-1)**(Jp+(levtz(i)%j2+levtz(j)%j2)/2))
        endif
       endif

       return
      end
