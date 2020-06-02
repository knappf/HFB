      subroutine occupations

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       logical :: ocp,ocn

       allocate(Up_HFB(id,id),Un_HFB(id,id),Vp_HFB(id,id),Vn_HFB(id,id))
       Up_HFB=0.d0
       Un_HFB=0.d0
       Vp_HFB=0.d0
       Vn_HFB=0.d0

       allocate(Ap_HFB(id,id),An_HFB(id,id),Bp_HFB(id,id),Bn_HFB(id,id))
       Ap_HFB=0.d0
       An_HFB=0.d0
       Bp_HFB=0.d0
       Bn_HFB=0.d0

!*************************************************************************
!   Here the initial occupations are determined.                         *
       ipp=0
       inn=0
       ocp=.true.
       ocn=.true.
       do i=1,id
        ipp=ipp+levp(i)%j2+1
        inn=inn+levn(i)%j2+1
        if(.not.ocp) then
         levp(i)%vi=0.d0
         levp(i)%ui=1.d0 !dsqrt(dble(levp(i)%j2+1))
         Vp_HFB(i,i)=0.d0
         Up_HFB(i,i)=1.d0 !dsqrt(dble(levp(i)%j2+1))
        endif
        if(ocp) then
         if(ipp.le.AZ) then
          levp(i)%vi=1.d0 !dsqrt(dble(levp(i)%j2+1))
          levp(i)%ui=0.d0
          Vp_HFB(i,i)=1.d0 !dsqrt(dble(levp(i)%j2+1))
          Up_HFB(i,i)=0.d0
         endif
         if(ipp.gt.AZ) then
          ocp=.false.
          levp(i)%vi=
     &           (dble(levp(i)%j2+1-ipp+AZ)/dble(levp(i)%j2+1))**0.5d0 !dsqrt(dble(levp(i)%j2+1-ipp+AZ))
          levp(i)%ui=dsqrt(dble(ipp-AZ)/dble(levp(i)%j2+1)) !dsqrt(dble(ipp-AZ))
          Vp_HFB(i,i)=levp(i)%vi
          Up_HFB(i,i)=levp(i)%ui
         endif
        endif
        if(.not.ocn) then
         levn(i)%vi=0.d0
         levn(i)%ui=1.d0 !dsqrt(dble(levn(i)%j2+1))
         Vn_HFB(i,i)=0.d0
         Un_HFB(i,i)=1.d0 !dsqrt(dble(levn(i)%j2+1))
        endif
        if(ocn) then
         if(inn.le.AN) then
          levn(i)%vi=1.d0 !dsqrt(dble(levn(i)%j2+1))
          levn(i)%ui=0.d0
          Vn_HFB(i,i)=1.d0 !dsqrt(dble(levn(i)%j2+1))
          Un_HFB(i,i)=0.d0
         endif
         if(inn.gt.AN) then
          ocn=.false.
          levn(i)%vi=
     &           (dble(levn(i)%j2+1-inn+AN)/dble(levn(i)%j2+1))**0.5d0 !dsqrt(dble(levn(i)%j2+1-inn+AN))
          levn(i)%ui=dsqrt(dble(inn-AN)/dble(levn(i)%j2+1)) !dsqrt(dble(inn-AN))
          Vn_HFB(i,i)=levn(i)%vi
          Un_HFB(i,i)=levn(i)%ui
         endif
        endif
       enddo
!*************************************************************************

       open(1,file='HO_scheme.out',status='unknown',form='formatted')
        write(1,*) 'List of HO levels:'
        write(1,*)'i,    N,    nn,    l,    2*j,    ener,   vi'
        do ii = 1, id
         write(1,'(1x,5(i4,1x),2(f12.5,1x))') 
     &    levp(ii)%index,levp(ii)%N,levp(ii)%nn,levp(ii)%l,
     &       levp(ii)%j2,levp(ii)%spenrg,levp(ii)%vi
        end do
       close(1)

       lhfp=levp
       lhfn=levn

       allocate(rhop_HFB(id,id), rhon_HFB(id,id))
       allocate(kapp_HFB(id,id), kapn_HFB(id,id))
       rhop_HFB=0.d0
       rhon_HFB=0.d0
       kapp_HFB=0.d0
       kapn_HFB=0.d0

       call make_densities

       if(.not.ifp_hfb) then
        e_filled=0.d0
        e_empty=0.d0
        do i=1,id
         if(dabs(lhfp(i)%vi-1.d0).lt.precis) e_filled=lhfp(i)%ei
         if(dabs(lhfp(id-i+1)%ui-1.d0).lt.precis)e_empty=lhfp(id-i+1)%ei
        enddo
        ferp=0.5d0*(e_empty+e_filled)
       endif

       if(.not.ifn_hfb) then
        e_filled=0.d0
        e_empty=0.d0
        do i=1,id
         if(dabs(lhfn(i)%vi-1.d0).lt.precis) e_filled=lhfn(i)%ei
         if(dabs(lhfn(id-i+1)%ui-1.d0).lt.precis)e_empty=lhfn(id-i+1)%ei
        enddo
        fern=0.5d0*(e_empty+e_filled)
       endif

       if(ifp_hfb) then
        do i=1,id
         if(dabs(lhfp(i)%vi).gt.precis) then
          ferp=lhfp(i)%ei
         endif
        enddo
        dferp=0.05d0
       else
        dferp=0.d0
       endif

       if(ifn_hfb) then
        do i=1,id
         if(dabs(lhfn(i)%vi).gt.precis) then
          fern=lhfn(i)%ei
         endif
        enddo
        dfern=0.05d0
       else
        dfern=0.d0
       endif

!*******************************************************************************
!    Hypernuclear part

       allocate(UY_HFB(id,id),VY_HFB(id,id),AY_HFB(id,id),BY_HFB(id,id))
       UY_HFB=0.d0
       VY_HFB=0.d0
       AY_HFB=0.d0
       BY_HFB=0.d0

       allocate(rhoY_HFB(id,id), kapY_HFB(id,id))
       rhoY_HFB=0.d0
       kapY_HFB=0.d0

       lhfY=levY

       if(if_self.eq.1) then
        rhoY_HFB(n_lam_occup,n_lam_occup)=1.d0
        lhfY(n_lam_occup)%vi=1.d0
        lhfY(n_lam_occup)%ui=0.d0
        levY(n_lam_occup)%vi=1.d0
        levY(n_lam_occup)%ui=0.d0
       endif

       return
      end
