      subroutine fermi_energy_p

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       if(.not.ifp_hfb) then
        e_filled=0.d0
        e_empty=0.d0
        do i=1,id
         if(dabs(lhfp(i)%vi-dsqrt(dble(lhfp(i)%j2+1))).lt.precis) 
     &                                               e_filled=lhfp(i)%ei
         if(dabs(lhfp(id-i+1)%ui-dsqrt(dble(lhfp(id-i+1)%j2+1))).lt.
     &                                    precis)e_empty=lhfp(id-i+1)%ei
        enddo
        ferp=0.5d0*(e_empty+e_filled)
       endif

       if(ifp_hfb) then
        ferp=ferp+dferp
       endif

       return
      end
