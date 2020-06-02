      subroutine spur_vectz(spuv,i2,JJ,ph_s)

       USE technical
       USE math
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: spuv(i2)
       TYPE(twoquas_type) :: ph_s(i2)

       spuv=0.d0

       if(JJ.eq.1) then
        do k=1,i2
         l1=ph_s(k)%q1
         l1tz=levtz(l1)%tz
         l1p=levtz(l1)%point
         l2=ph_s(k)%q2
         l2tz=levtz(l2)%tz
         l2p=levtz(l2)%point
         if(l1tz.eq.l2tz.and.l1tz.eq.1) then
          spuv(k)=dsqrt(4.d0*pi/9.d0)*trE1_n(l2p,l1p)/dble(AZ+AN)
         endif
         if(l1tz.eq.l2tz.and.l1tz.eq.-1) then
          spuv(k)=dsqrt(4.d0*pi/9.d0)*trE1_p(l2p,l1p)/dble(AZ+AN)
         endif
        enddo
       endif

       dnor=0.d0
       do k=1,i2
        dnor=dnor+spuv(k)**2.d0
       enddo
       do k=1,i2
        spuv(k)=spuv(k)/dsqrt(dnor)
       enddo

       return
      end
