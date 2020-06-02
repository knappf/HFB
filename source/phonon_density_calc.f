      subroutine phonon_density_calc(iipar,Jc,phh,TD,ww,idi)

       USE technical
       USE math
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       TYPE(twoquas_type) :: phh(idi)
       double precision :: TD(idi,idi),ww(idi)

       if(.not.(if_ort.eq.1.and.((iipar.eq.-1.and.Jc.eq.1).or.
     &        ((iipar.eq.1.and.Jc.eq.0).and.if_QTDA.eq.1)))) then

       if(iipar.eq.-1.and.Jc.eq.1) then
        open(11,file='Tran_dens.out',form='formatted',status='unknown')
        write(11,*) 'E1 transition densities'
        do i=1,idi
         write(11,*) 'TDA phonon=',ww(i)
         do m=0,igrid2
          radi=dble(m)*sizebox/dble(igrid2)
          E1_densp=0.d0
          E1_densn=0.d0
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
           if(itz.eq.-1) E1_densp=trE1_p_dens(ih1,ip1,m)*ff1*TD(j,i)
     &                                                        +E1_densp
           if(itz.eq.+1) E1_densn=trE1_n_dens(ih1,ip1,m)*ff1*TD(j,i)
     &                                                        +E1_densn
          enddo
          write(11,*) radi,E1_densp,E1_densn
         enddo
        enddo
        close(11)
       endif

       elseif(iipar.eq.-1.and.Jc.eq.1) then

       if(iipar.eq.-1.and.Jc.eq.1) then
        open(11,file='Tran_dens.out',form='formatted',status='unknown')
        write(11,*) 'E1 transition densities'
        do i=2,idi
         write(11,*) 'TDA phonon=',ww(i)
         do m=0,igrid2
          radi=dble(m)*sizebox/dble(igrid2)
          E1_densp=0.d0
          E1_densn=0.d0
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
           if(itz.eq.-1) E1_densp=trE1_p_dens(ih1,ip1,m)*ff1*TD(j,i)
     &                                                        +E1_densp
           if(itz.eq.+1) E1_densn=trE1_n_dens(ih1,ip1,m)*ff1*TD(j,i)
     &                                                        +E1_densn
          enddo
          write(11,*) radi,E1_densp,E1_densn
         enddo
        enddo
        close(11)
       endif

       endif

       return
      end
