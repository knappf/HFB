      subroutine fermi_energy_n

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       if(.not.ifn_hfb) then
        e_filled=0.d0
        e_empty=0.d0
        do i=1,id
         if(dabs(lhfn(i)%vi-dsqrt(dble(lhfn(i)%j2+1))).lt.precis) 
     &                                               e_filled=lhfn(i)%ei
         if(dabs(lhfn(id-i+1)%ui-dsqrt(dble(lhfn(id-i+1)%j2+1))).lt.
     &                                    precis)e_empty=lhfn(id-i+1)%ei
        enddo
        fern=0.5d0*(e_empty+e_filled)
       endif

       if(ifn_hfb) then
        fern=fern+dfern
       endif

       return
      end
