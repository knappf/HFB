      double precision function VNLtz(i,j,k,l,Jp)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       VNLtz=0.d0

       if(levtz(i)%tz.eq.-1.and.levtz(k)%tz.eq.-1) then
        VNLtz=VpY(levtz(i)%point,j,levtz(k)%point,l,Jp)
       endif
       if(levtz(i)%tz.eq.1.and.levtz(k)%tz.eq.1) then
        VNLtz=VnY(levtz(i)%point,j,levtz(k)%point,l,Jp)
       endif

       return
      end
