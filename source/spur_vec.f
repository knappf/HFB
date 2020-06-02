      subroutine spur_vec(ii1,spuv,i2,JJ,ph_s)

       USE technical
       USE math
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: spuv(i2)
       TYPE(twoquas_type) :: ph_s(i2)

       spuv=0.d0

       if(JJ.eq.0) then      !  the N spuriosity
        if(ii1.eq.1) then
         do k=1,i2
          if(ph_s(k)%tz.eq.1) then
           if(ph_s(k)%q1.eq.ph_s(k)%q2) then
            ll=ph_s(k)%q1
            spuv(k)=dsqrt(dble(2*(lhfn(ll)%j2+1)))
     &                       *lhfn(ll)%ui*lhfn(ll)%vi
           endif
          endif
          if(ph_s(k)%tz.eq.-1) then
           if(ph_s(k)%q1.eq.ph_s(k)%q2) then
            ll=ph_s(k)%q1
            spuv(k)=dsqrt(dble(2*(lhfp(ll)%j2+1)))
     &                       *lhfp(ll)%ui*lhfp(ll)%vi
           endif
          endif
         enddo
        endif
       endif

       if(JJ.eq.1) then      !  the CM spuriosity
        do k=1,i2
         if(ph_s(k)%tz.eq.1) then
          l1=ph_s(k)%q1
          l2=ph_s(k)%q2
          if(ii1.eq.1)
     &     spuv(k)=dsqrt(4.d0*pi/9.d0)*trE1_n(l2,l1)/dble(AZ+AN)
          if(ii1.eq.0) 
     &     spuv(k)=dsqrt(4.d0*pi/9.d0)*trE1_n(l1,l2)/dble(AZ+AN)
!          if(ii1.eq.1) then
           spuv(k)=spuv(k)
     &       *(lhfn(l1)%ui*lhfn(l2)%vi-lhfn(l1)%vi*lhfn(l2)%ui)
!          endif
          if(ii1.eq.0) then
           spuv(k)=spuv(k)*(-1)**((lhfn(l1)%j2-lhfn(l2)%j2)/2)
          endif
         endif
         if(ph_s(k)%tz.eq.-1) then
          l1=ph_s(k)%q1
          l2=ph_s(k)%q2
          if(ii1.eq.1)
     &     spuv(k)=dsqrt(4.d0*pi/9.d0)*trE1_p(l2,l1)/dble(AZ+AN)
          if(ii1.eq.0)
     &     spuv(k)=dsqrt(4.d0*pi/9.d0)*trE1_p(l1,l2)/dble(AZ+AN)
!          if(ii1.eq.1) then
           spuv(k)=spuv(k)
     &       *(lhfp(l1)%ui*lhfp(l2)%vi-lhfp(l1)%vi*lhfp(l2)%ui)
!          endif
          if(ii1.eq.0) then
           spuv(k)=spuv(k)*(-1)**((lhfp(l1)%j2-lhfp(l2)%j2)/2)
          endif
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
