      subroutine MBPT_energy(ene2,ene2Y)

       USE technical
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision, allocatable, save :: ccj(:,:,:,:,:)

       ene2=0.d0
       ene2Y=0.d0

       allocate(ccj(jmax,-jmax:jmax,jmax,-jmax:jmax,0:jmax))
       ccj=0.d0

       do ii1=1,jmax,2
        do ii2=-jmax,jmax,2
         do ii3=1,jmax,2
          do ii4=-jmax,jmax,2
           do ii5=0,jmax
            ccj(ii1,ii2,ii3,ii4,ii5)=cleb(ii1,ii2,ii3,ii4,2*ii5,ii2+ii4)
           enddo
          enddo
         enddo
        enddo
       enddo

       do ip1=1,id
        if(dabs(lhfp(ip1)%ui-1.d0).lt.1.d-7) then
        do ip2=ip1,id
         if(dabs(lhfp(ip2)%ui-1.d0).lt.1.d-7) then
         do ih1=1,id
          if(dabs(lhfp(ih1)%vi-1.d0).lt.1.d-7) then
          do ih2=ih1,id
           if(dabs(lhfp(ih2)%vi-1.d0).lt.1.d-7) then

       do mp1=-lhfp(ip1)%j2,lhfp(ip1)%j2,2
        if(ip1.eq.ip2) then
         do mp2=mp1+2,lhfp(ip2)%j2,2

          do mh1=-lhfp(ih1)%j2,lhfp(ih1)%j2,2
           if(ih1.eq.ih2) then
            do mh2=mh1+2,lhfp(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vpp(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfp(ip1)%j2,mp1,lhfp(ip2)%j2,mp2,jjc)
     &         *ccj(lhfp(ih1)%j2,mh1,lhfp(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei-lhfp(ip2)%ei)
              endif
             enddo
            enddo
           else
            do mh2=-lhfp(ih2)%j2,lhfp(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vpp(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfp(ip1)%j2,mp1,lhfp(ip2)%j2,mp2,jjc)
     &         *ccj(lhfp(ih1)%j2,mh1,lhfp(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei-lhfp(ip2)%ei)
              endif
             enddo
            enddo
           endif
          enddo

         enddo
        else
         do mp2=-lhfp(ip2)%j2,lhfp(ip2)%j2,2

          do mh1=-lhfp(ih1)%j2,lhfp(ih1)%j2,2
           if(ih1.eq.ih2) then
            do mh2=mh1+2,lhfp(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vpp(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfp(ip1)%j2,mp1,lhfp(ip2)%j2,mp2,jjc)
     &         *ccj(lhfp(ih1)%j2,mh1,lhfp(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei-lhfp(ip2)%ei)
              endif
             enddo
            enddo
           else
            do mh2=-lhfp(ih2)%j2,lhfp(ih2)%j2,2
             do jjc=0,jmax
               if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vpp(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfp(ip1)%j2,mp1,lhfp(ip2)%j2,mp2,jjc)
     &         *ccj(lhfp(ih1)%j2,mh1,lhfp(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei-lhfp(ip2)%ei)
              endif
             enddo
            enddo
           endif
          enddo

         enddo
        endif
       enddo

           endif
          enddo
          endif
         enddo
         endif
        enddo
        endif
       enddo

       do ip1=1,id
        if(dabs(lhfn(ip1)%ui-1.d0).lt.1.d-7) then
        do ip2=ip1,id
         if(dabs(lhfn(ip2)%ui-1.d0).lt.1.d-7) then
         do ih1=1,id
          if(dabs(lhfn(ih1)%vi-1.d0).lt.1.d-7) then
          do ih2=ih1,id
           if(dabs(lhfn(ih2)%vi-1.d0).lt.1.d-7) then

       do mp1=-lhfn(ip1)%j2,lhfn(ip1)%j2,2
        if(ip1.eq.ip2) then
         do mp2=mp1+2,lhfn(ip2)%j2,2

          do mh1=-lhfn(ih1)%j2,lhfn(ih1)%j2,2
           if(ih1.eq.ih2) then
            do mh2=mh1+2,lhfn(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vnn(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfn(ip1)%j2,mp1,lhfn(ip2)%j2,mp2,jjc)
     &         *ccj(lhfn(ih1)%j2,mh1,lhfn(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei-lhfn(ip2)%ei)
              endif
             enddo
            enddo
           else
            do mh2=-lhfn(ih2)%j2,lhfn(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vnn(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfn(ip1)%j2,mp1,lhfn(ip2)%j2,mp2,jjc)
     &         *ccj(lhfn(ih1)%j2,mh1,lhfn(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei-lhfn(ip2)%ei)
              endif
             enddo
            enddo
           endif
          enddo

         enddo
        else
         do mp2=-lhfn(ip2)%j2,lhfn(ip2)%j2,2

          do mh1=-lhfn(ih1)%j2,lhfn(ih1)%j2,2
           if(ih1.eq.ih2) then
            do mh2=mh1+2,lhfn(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vnn(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfn(ip1)%j2,mp1,lhfn(ip2)%j2,mp2,jjc)
     &         *ccj(lhfn(ih1)%j2,mh1,lhfn(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei-lhfn(ip2)%ei)
              endif
             enddo
            enddo
           else
            do mh2=-lhfn(ih2)%j2,lhfn(ih2)%j2,2
             do jjc=0,jmax
              if(mp1+mp2.eq.mh1+mh2) then
               ene2=ene2+(Vnn(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfn(ip1)%j2,mp1,lhfn(ip2)%j2,mp2,jjc)
     &         *ccj(lhfn(ih1)%j2,mh1,lhfn(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei-lhfn(ip2)%ei)
              endif
             enddo
            enddo
           endif
          enddo

         enddo
        endif
       enddo

           endif
          enddo
          endif
         enddo
         endif
        enddo
        endif
       enddo

       do ip1=1,id
        if(dabs(lhfp(ip1)%ui-1.d0).lt.1.d-7) then
        do ip2=1,id
         if(dabs(lhfn(ip2)%ui-1.d0).lt.1.d-7) then
         do ih1=1,id
          if(dabs(lhfp(ih1)%vi-1.d0).lt.1.d-7) then
          do ih2=1,id
           if(dabs(lhfn(ih2)%vi-1.d0).lt.1.d-7) then

       do mp1=-lhfp(ip1)%j2,lhfp(ip1)%j2,2
        do mp2=-lhfn(ip2)%j2,lhfn(ip2)%j2,2
         do mh1=-lhfp(ih1)%j2,lhfp(ih1)%j2,2
          do mh2=-lhfn(ih2)%j2,lhfn(ih2)%j2,2

           do jjc=0,jmax
            if(mp1+mp2.eq.mh1+mh2) then
             ene2=ene2+(Vpn(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfp(ip1)%j2,mp1,lhfn(ip2)%j2,mp2,jjc)
     &         *ccj(lhfp(ih1)%j2,mh1,lhfn(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfp(ih1)%ei+lhfn(ih2)%ei-lhfp(ip1)%ei-lhfn(ip2)%ei)
            endif
           enddo

          enddo
         enddo
        enddo
       enddo

           endif
          enddo
          endif
         enddo
         endif
        enddo
        endif
       enddo

       do ip1=1,id
        if(dabs(lhfp(ip1)%ui-1.d0).lt.1.d-7) then
        do ip2=2,id
!         if(dabs(lhfY(ip2)%ui-1.d0).lt.1.d-7) then
         do ih1=1,id
          if(dabs(lhfp(ih1)%vi-1.d0).lt.1.d-7) then
          do ih2=1,1!id
!           if(dabs(lhfY(ih2)%vi-1.d0).lt.1.d-7) then

       do mp1=-lhfp(ip1)%j2,lhfp(ip1)%j2,2
        do mp2=-lhfY(ip2)%j2,lhfY(ip2)%j2,2
         do mh1=-lhfp(ih1)%j2,lhfp(ih1)%j2,2
          do mh2=-lhfY(ih2)%j2,lhfY(ih2)%j2,2

           do jjc=0,jmax
            if(mp1+mp2.eq.mh1+mh2) then
             ene2Y=ene2Y+(VpY(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfp(ip1)%j2,mp1,lhfY(ip2)%j2,mp2,jjc)
     &         *ccj(lhfp(ih1)%j2,mh1,lhfY(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfp(ih1)%ei+lhfY(ih2)%ei-lhfp(ip1)%ei-lhfY(ip2)%ei)
            endif
           enddo

          enddo
         enddo
        enddo
       enddo

!           endif
          enddo
          endif
         enddo
!         endif
        enddo
        endif
       enddo

       do ip1=1,id
        if(dabs(lhfn(ip1)%ui-1.d0).lt.1.d-7) then
        do ip2=2,id
!         if(dabs(lhfY(ip2)%ui-1.d0).lt.1.d-7) then
         do ih1=1,id
          if(dabs(lhfn(ih1)%vi-1.d0).lt.1.d-7) then
          do ih2=1,1!id
!           if(dabs(lhfY(ih2)%vi-1.d0).lt.1.d-7) then

       do mp1=-lhfn(ip1)%j2,lhfn(ip1)%j2,2
        do mp2=-lhfY(ip2)%j2,lhfY(ip2)%j2,2
         do mh1=-lhfn(ih1)%j2,lhfn(ih1)%j2,2
          do mh2=-lhfY(ih2)%j2,lhfY(ih2)%j2,2

           do jjc=0,jmax
            if(mp1+mp2.eq.mh1+mh2) then
             ene2Y=ene2Y+(VnY(ip1,ip2,ih1,ih2,jjc)
     &         *ccj(lhfn(ip1)%j2,mp1,lhfY(ip2)%j2,mp2,jjc)
     &         *ccj(lhfn(ih1)%j2,mh1,lhfY(ih2)%j2,mh2,jjc))**2.d0
     &         /(lhfn(ih1)%ei+lhfY(ih2)%ei-lhfn(ip1)%ei-lhfY(ip2)%ei)
            endif
           enddo

          enddo
         enddo
        enddo
       enddo

!           endif
          enddo
          endif
         enddo
!         endif
        enddo
        endif
       enddo

       ene2Y=ene2Y/dble(lhfY(1)%j2+1)

       deallocate(ccj)

       return
      end
