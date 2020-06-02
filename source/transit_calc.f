      subroutine transit_calc(iipar,Jc,phh,TD,ww,idi)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       TYPE(twoquas_type) :: phh(idi)
       double precision :: TD(idi,idi),ww(idi)

       if(iipar.eq.+1.and.Jc.eq.0) then
        open(11,file='E0_is.out',form='formatted',status='unknown')
        open(12,file='E0_iv.out',form='formatted',status='unknown')
        open(13,file='N_transit.out',form='formatted',status='unknown')
        open(14,file='E0_ph.out',form='formatted',status='unknown')
        write(11,*) 'isE0: E_i      B_i'
        write(12,*) 'ivE0: E_i      B_i'
        write(13,*) '<TDA|N|0>: E_i      B_i'
        write(14,*) 'physE0: E_i      B_i'
        do i=1,idi
         E0_is=0.d0
         E0_iv=0.d0
         EN=0.d0
         E0_ph=0.d0
         do j=1,idi
          if(if_QTDA.eq.0) ip1=phh(j)%q2
          if(if_QTDA.eq.0) ih1=phh(j)%q1
          if(if_QTDA.eq.1) ip1=phh(j)%q1
          if(if_QTDA.eq.1) ih1=phh(j)%q2
          itz=phh(j)%tz
          if(itz.eq.-1) ff1=(lhfp(ip1)%ui*lhfp(ih1)%vi
     &                           +lhfp(ip1)%vi*lhfp(ih1)%ui)
          if(itz.eq.1) ff1=(lhfn(ip1)%ui*lhfn(ih1)%vi
     &                           +lhfn(ip1)%vi*lhfn(ih1)%ui)
          if(if_QTDA.eq.0.and.itz.eq.-1) ff1=ff1*
     &                  (-1)**((lhfp(ip1)%j2-lhfp(ih1)%j2)/2)
          if(if_QTDA.eq.0.and.itz.eq.+1) ff1=ff1*
     &                  (-1)**((lhfn(ip1)%j2-lhfn(ih1)%j2)/2)
          if(ip1.eq.ih1) ff1=ff1/dsqrt(2.d0)
          if(itz.eq.-1) E0_is=0.5d0*trE0_p(ih1,ip1)*ff1*TD(j,i)+E0_is
          if(itz.eq.+1) E0_is=0.5d0*trE0_n(ih1,ip1)*ff1*TD(j,i)+E0_is
          if(itz.eq.-1) E0_iv=-0.5d0*trE0_p(ih1,ip1)*ff1*TD(j,i)+E0_iv
          if(itz.eq.+1) E0_iv=0.5d0*trE0_n(ih1,ip1)*ff1*TD(j,i)+E0_iv
          if(itz.eq.-1) EN=trEN_p(ih1,ip1)*ff1*TD(j,i)+EN
          if(itz.eq.+1) EN=trEN_n(ih1,ip1)*ff1*TD(j,i)+EN
          if(itz.eq.-1) E0_ph=1.0d0*trE0_p(ih1,ip1)*ff1*TD(j,i)+E0_ph
         enddo
         write(11,*) ww(i),E0_is**2.d0
         write(12,*) ww(i),E0_iv**2.d0
         write(13,*) ww(i),EN**2.d0
         write(14,*) ww(i),E0_ph**2.d0
        enddo
        close(11)
        close(12)
        close(13)
        close(14)
       endif

       if(iipar.eq.-1.and.Jc.eq.1) then
        open(11,file='E1_is.out',form='formatted',status='unknown')
        open(12,file='E1_iv.out',form='formatted',status='unknown')
        open(13,file='Squeeze.out',form='formatted',status='unknown')
        open(14,file='E1_ph.out',form='formatted',status='unknown')
        write(11,*) 'isE1: E_i      B_i'
        write(12,*) 'ivE1: E_i      B_i'
        write(13,*) 'Squeeze: E_i      B_i'
        write(14,*) 'physE1: E_i      B_i'
        do i=1,idi
         E1_is=0.d0
         E1_iv=0.d0
         S1=0.d0
         E1_ph=0.d0
         do j=1,idi
          if(if_QTDA.eq.0) ip1=phh(j)%q2
          if(if_QTDA.eq.0) ih1=phh(j)%q1
          if(if_QTDA.eq.1) ip1=phh(j)%q1
          if(if_QTDA.eq.1) ih1=phh(j)%q2
          itz=phh(j)%tz
          if(itz.eq.-1) ff1=(lhfp(ip1)%ui*lhfp(ih1)%vi
     &                           -lhfp(ip1)%vi*lhfp(ih1)%ui)
          if(itz.eq.1) ff1=(lhfn(ip1)%ui*lhfn(ih1)%vi
     &                           -lhfn(ip1)%vi*lhfn(ih1)%ui)
          if(if_QTDA.eq.0.and.itz.eq.-1) ff1=ff1*
     &                  (-1)**((lhfp(ip1)%j2-lhfp(ih1)%j2)/2)
          if(if_QTDA.eq.0.and.itz.eq.+1) ff1=ff1*
     &                  (-1)**((lhfn(ip1)%j2-lhfn(ih1)%j2)/2)
          if(ip1.eq.ih1) ff1=ff1/dsqrt(2.d0)
          if(itz.eq.-1) E1_is=0.5d0*trE1_p(ih1,ip1)*ff1*TD(j,i)+E1_is
          if(itz.eq.+1) E1_is=0.5d0*trE1_n(ih1,ip1)*ff1*TD(j,i)+E1_is
          if(itz.eq.-1) E1_iv=-0.5d0*trE1_p(ih1,ip1)*ff1*TD(j,i)+E1_iv
          if(itz.eq.+1) E1_iv=0.5d0*trE1_n(ih1,ip1)*ff1*TD(j,i)+E1_iv
          if(itz.eq.-1) S1=trS1_p(ih1,ip1)*ff1*TD(j,i)+S1
          if(itz.eq.+1) S1=trS1_n(ih1,ip1)*ff1*TD(j,i)+S1
          if(itz.eq.-1) E1_ph=1.0d0*trE1_p(ih1,ip1)*ff1*TD(j,i)+E1_ph
         enddo
         write(11,*) ww(i),E1_is**2.d0
         write(12,*) ww(i),E1_iv**2.d0
         write(13,*) ww(i),S1**2.d0
         write(14,*) ww(i),E1_ph**2.d0
        enddo
        close(11)
        close(12)
        close(13)
        close(14)
       endif

       if(iipar.eq.+1.and.Jc.eq.2) then
        open(11,file='E2_is.out',form='formatted',status='unknown')
        open(12,file='E2_iv.out',form='formatted',status='unknown')
        open(13,file='E2_ph.out',form='formatted',status='unknown')
        write(11,*) 'isE2: E_i      B_i'
        write(12,*) 'ivE2: E_i      B_i'
        write(13,*) 'physE2: E_i      B_i'
        do i=1,idi
         E2_is=0.d0
         E2_iv=0.d0
         E2_ph=0.d0
         do j=1,idi
          if(if_QTDA.eq.0) ip1=phh(j)%q2
          if(if_QTDA.eq.0) ih1=phh(j)%q1
          if(if_QTDA.eq.1) ip1=phh(j)%q1
          if(if_QTDA.eq.1) ih1=phh(j)%q2
          itz=phh(j)%tz
          if(itz.eq.-1) ff1=(lhfp(ip1)%ui*lhfp(ih1)%vi
     &                           +lhfp(ip1)%vi*lhfp(ih1)%ui)
          if(itz.eq.1) ff1=(lhfn(ip1)%ui*lhfn(ih1)%vi
     &                           +lhfn(ip1)%vi*lhfn(ih1)%ui)
          if(if_QTDA.eq.0.and.itz.eq.-1) ff1=ff1*
     &                  (-1)**((lhfp(ip1)%j2-lhfp(ih1)%j2)/2)
          if(if_QTDA.eq.0.and.itz.eq.+1) ff1=ff1*
     &                  (-1)**((lhfn(ip1)%j2-lhfn(ih1)%j2)/2)
          if(ip1.eq.ih1) ff1=ff1/dsqrt(2.d0)
          if(itz.eq.-1) E2_is=0.5d0*trE2_p(ih1,ip1)*ff1*TD(j,i)+E2_is
          if(itz.eq.+1) E2_is=0.5d0*trE2_n(ih1,ip1)*ff1*TD(j,i)+E2_is
          if(itz.eq.-1) E2_iv=-0.5d0*trE2_p(ih1,ip1)*ff1*TD(j,i)+E2_iv
          if(itz.eq.+1) E2_iv=0.5d0*trE2_n(ih1,ip1)*ff1*TD(j,i)+E2_iv
          if(itz.eq.-1) E2_ph=1.0d0*trE2_p(ih1,ip1)*ff1*TD(j,i)+E2_ph
         enddo
         write(11,*) ww(i),E2_is**2.d0
         write(12,*) ww(i),E2_iv**2.d0
         write(13,*) ww(i),E2_ph**2.d0
        enddo
        close(11)
        close(12)
        close(13)
       endif

       if(iipar.eq.-1.and.Jc.eq.3) then
        open(11,file='E3_is.out',form='formatted',status='unknown')
        open(12,file='E3_iv.out',form='formatted',status='unknown')
        open(13,file='E3_ph.out',form='formatted',status='unknown')
        write(11,*) 'isE3: E_i      B_i'
        write(12,*) 'ivE3: E_i      B_i'
        write(13,*) 'physE3: E_i      B_i'
        do i=1,idi
         E3_is=0.d0
         E3_iv=0.d0
         E3_ph=0.d0
         do j=1,idi
          if(if_QTDA.eq.0) ip1=phh(j)%q2
          if(if_QTDA.eq.0) ih1=phh(j)%q1
          if(if_QTDA.eq.1) ip1=phh(j)%q1
          if(if_QTDA.eq.1) ih1=phh(j)%q2
          itz=phh(j)%tz
          if(itz.eq.-1) ff1=(lhfp(ip1)%ui*lhfp(ih1)%vi
     &                           -lhfp(ip1)%vi*lhfp(ih1)%ui)
          if(itz.eq.1) ff1=(lhfn(ip1)%ui*lhfn(ih1)%vi
     &                           -lhfn(ip1)%vi*lhfn(ih1)%ui)
          if(if_QTDA.eq.0.and.itz.eq.-1) ff1=ff1*
     &                  (-1)**((lhfp(ip1)%j2-lhfp(ih1)%j2)/2)
          if(if_QTDA.eq.0.and.itz.eq.+1) ff1=ff1*
     &                  (-1)**((lhfn(ip1)%j2-lhfn(ih1)%j2)/2)
          if(ip1.eq.ih1) ff1=ff1/dsqrt(2.d0)
          if(itz.eq.-1) E3_is=0.5d0*trE3_p(ih1,ip1)*ff1*TD(j,i)+E3_is
          if(itz.eq.+1) E3_is=0.5d0*trE3_n(ih1,ip1)*ff1*TD(j,i)+E3_is
          if(itz.eq.-1) E3_iv=-0.5d0*trE3_p(ih1,ip1)*ff1*TD(j,i)+E3_iv
          if(itz.eq.+1) E3_iv=0.5d0*trE3_n(ih1,ip1)*ff1*TD(j,i)+E3_iv
          if(itz.eq.-1) E3_ph=1.0d0*trE3_p(ih1,ip1)*ff1*TD(j,i)+E3_ph
         enddo
         write(11,*) ww(i),E3_is**2.d0
         write(12,*) ww(i),E3_iv**2.d0
         write(13,*) ww(i),E3_ph**2.d0
        enddo
        close(11)
        close(12)
        close(13)
       endif

       if(iipar.eq.1.and.Jc.eq.1) then
        qprot=1.d0 !0.68d0
        qneut=1.d0 !0.64d0
        gsprot=5.59d0
        gsneut=-3.83d0
        gyrmu=1.d0 !0.105155d0
        open(11,file='M1spin.out',form='formatted',status='unknown')
        open(12,file='M1orbi.out',form='formatted',status='unknown')
        open(13,file='M1.out',form='formatted',status='unknown')
        write(11,*) 'spin M1: E_i      B_i'
        write(12,*) 'orbital M1: E_i      B_i'
        write(13,*) 'phys M1: E_i      B_i'
        do i=1,idi
         V1s=0.d0
         V1l=0.d0
         V1ph=0.d0
         do j=1,idi
          if(if_QTDA.eq.0) ip1=phh(j)%q2
          if(if_QTDA.eq.0) ih1=phh(j)%q1
          if(if_QTDA.eq.1) ip1=phh(j)%q1
          if(if_QTDA.eq.1) ih1=phh(j)%q2
          itz=phh(j)%tz
          if(itz.eq.-1) ff1=(lhfp(ip1)%ui*lhfp(ih1)%vi
     &                           -lhfp(ip1)%vi*lhfp(ih1)%ui)
          if(itz.eq.1) ff1=(lhfn(ip1)%ui*lhfn(ih1)%vi
     &                           -lhfn(ip1)%vi*lhfn(ih1)%ui)
          if(if_QTDA.eq.0.and.itz.eq.-1) ff1=ff1*
     &                  (-1)**((lhfp(ip1)%j2-lhfp(ih1)%j2)/2)
          if(if_QTDA.eq.0.and.itz.eq.+1) ff1=ff1*
     &                  (-1)**((lhfn(ip1)%j2-lhfn(ih1)%j2)/2)
          if(ip1.eq.ih1) ff1=ff1/dsqrt(2.d0)
          if(itz.eq.-1) V1s=qprot*gsprot*gyrmu*
     &                            trM1s_p(ih1,ip1)*ff1*TD(j,i)+V1s
          if(itz.eq.+1) V1s=qneut*gsneut*gyrmu*
     &                            trM1s_n(ih1,ip1)*ff1*TD(j,i)+V1s
          if(itz.eq.-1) V1l=1.d0*trM1l_p(ih1,ip1)*ff1*TD(j,i)+V1l
          if(itz.eq.-1) V1ph=qprot*gsprot*gyrmu*
     &                            trM1s_p(ih1,ip1)*ff1*TD(j,i)
     &                      +1.d0*trM1l_p(ih1,ip1)*ff1*TD(j,i)+V1ph
          if(itz.eq.+1) V1ph=qneut*gsneut*gyrmu*
     &                            trM1s_n(ih1,ip1)*ff1*TD(j,i)+V1ph
         enddo
         write(11,*) ww(i),V1s**2.d0
         write(12,*) ww(i),V1l**2.d0
         write(13,*) ww(i),V1ph**2.d0
        enddo
        close(11)
        close(12)
        close(13)
       endif

       return
      end
