      subroutine hfb_iteration

       USE technical
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision, allocatable, save :: argu2(:)
       double precision, allocatable, save :: argup(:)
       double precision, allocatable, save :: argun(:)
       double precision, allocatable, save :: rad_denp(:)
       double precision, allocatable, save :: rad_denn(:)
       double precision, allocatable, save :: rad_den2(:)

       double precision :: hp(id,id),Dp(id,id)
       double precision :: hn(id,id),Dn(id,id)
       double precision :: wp(id),wn(id)

       double precision :: hY(id,id),wY(id)
       integer :: jlY(id)

       double precision :: rp(id,id),rn(id,id),rY(id,id)
       double precision :: kp(id,id),kn(id,id)

!       double precision :: h2p(id,id),h2n(id,id)
       double precision :: w2p(id),w2n(id)

!       double precision :: D2p(id,id),D2n(id,id)

       double precision :: E_diff, E_orig
       integer :: i_count1,i_count2,count_max

       alpha=0.90d0

       count_max=2000
       E_diff=1.d9

       call HFB_energy

       call hp_field(hp)
       call hn_field(hn)
       call Dp_field(Dp)
       call Dn_field(Dn)

       if(if_self.eq.1) call hY_field(hY)

       allocate(H11p(id,id),H11n(id,id))
       H11p=Dp
       H11n=Dn

       write(*,*) 'E_HFB,       lambda_p,         lambda_n'
       write(*,*) E_HFB,ferp,fern

       i_count1=0
       E_diff1=1.d9

       do while(((ifp_hfb.or.ifn_hfb).and.(dabs(E_diff1).gt.precis
     &         .and.i_count1.le.2)).or.((.not.(ifp_hfb.or.ifn_hfb)).and.
     &         (dabs(E_diff1).gt.precis.and.i_count1.le.count_max)))
        i_count1=i_count1+1
        E_orig1=E_HFB

        call diagonalization(hp,wp,id)
        call diagonalization(hn,wn,id)

        if(if_self.eq.1) call diagonalization(hY,wY,id)

        call make_sp_levels(hp,hn,wp,wn,i_count1)
        tran_p=hp
        tran_n=hn

        if(if_self.eq.1) then
         tran_Y=hY

         do i=1,id
          d_max=0.d0
          do j=1,id
           if(dabs(hY(j,i)).gt.d_max) then
            d_max=dabs(hY(j,i))
            j_pointer=j
           endif
          enddo
          jlY(i)=1000*levY(j_pointer)%l+levY(j_pointer)%j2
         enddo
         do i=1,id
          llll=jlY(i)/1000
          jjjj=mod(jlY(i),1000)
          lhfY(i)%ipar=(-1)**llll
          lhfY(i)%l=llll
          lhfY(i)%j2=jjjj
          lhfY(i)%ei=wY(i)
          lhfY(i)%qei=dabs(wY(i)-wY(1)+0.5d0)
         enddo
         do ll=0,(jmax+1)/2
          do jj=1,jmax,2
           nn=-1
           do i=1,id
            if(lhfY(i)%j2.eq.jj.and.lhfY(i)%l.eq.ll) then
             nn=nn+1
             lhfY(i)%nn=nn
            endif
           enddo
          enddo
         enddo
         do i=1,id
          do j=1,id
           if(levY(i)%nn.eq.lhfY(j)%nn.and.(levY(i)%l.eq.lhfY(j)%l
     &                           .and.levY(i)%j2.eq.lhfY(j)%j2)) then
            lhfY(j)%ui=levY(i)%ui
            lhfY(j)%vi=levY(i)%vi
           endif
          enddo
         enddo

         VY_HFB=0.d0
         UY_HFB=0.d0
         do i=1,id
          UY_HFB(i,i)=lhfY(i)%ui
          VY_HFB(i,i)=lhfY(i)%vi!*(-1)**((lhfY(i)%j2-1)/2)
         enddo
         call dgemm('N','N',id,id,id,1.d0,hY,max(1,id),VY_HFB,
     &                        max(1,id),0.d0,BY_HFB,max(1,id))
         VY_HFB=BY_HFB

         call dgemm('N','T',id,id,id,1.d0,VY_HFB,max(1,id),VY_HFB,
     &                          max(1,id),0.d0,rhoY_HFB,max(1,id))

         open(1,file='HF_Y.out',status='unknown',form='formatted')
          write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
          do ii = 1, id
           write(1,*) lhfY(ii)%index,lhfY(ii)%l,lhfY(ii)%j2,
     &      lhfY(ii)%ei,lhfY(ii)%qei,lhfY(ii)%vi
          end do
         close(1)
        endif

        call make_densities
        call check_densities

        bos1=dsqrt(0.5d0*(zmp+zmn)*hbarom/(hbarc**2.d0))
        bos2=dsqrt(zmY*hbarom/(hbarc**2.d0))
        rad_den1=0.d0
        do i=1,igrid
         val=0.d0
         radi=zcross(i)  !dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_den1(i)=val/(4.d0*pi)
        enddo
!********************************************************************************
       r4tot=0.d0
       r2tot=0.d0
       r4Z=0.d0
       r2Z=0.d0
       r4N=0.d0
       r2N=0.d0
       allocate(argu2(0:igrid2),argup(0:igrid2),argun(0:igrid2))
       argu2=0.d0
       argup=0.d0
       argun=0.d0
       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),
     &                                    rad_den2(0:igrid2))
       rad_denp=0.d0 
       rad_denn=0.d0
       rad_den2=0.d0
        do i=0,igrid2
         val1=0.d0
         val2=0.d0
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val1=val1+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val2=val2+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_denp(i)=val1/(4.d0*pi)
         rad_denn(i)=val2/(4.d0*pi)
         rad_den2(i)=(val1+val2)/(4.d0*pi)
        enddo
       do i=0,igrid2
        radi=dble(i)*sizebox/dble(igrid2)
        argu2(i)=radi**2.d0*rad_den2(i)
        argup(i)=radi**2.d0*rad_denp(i)
        argun(i)=radi**2.d0*rad_denn(i)
       enddo
       dx=sizebox/dble(igrid2)
       do i=0,igrid2-1
        radi=(dble(i)+0.5d0)*sizebox/dble(igrid2)
        r4tot=r4tot+dx*0.5d0*(argu2(i)+argu2(i+1))*radi**2.d0
        r2tot=r2tot+dx*0.5d0*(argu2(i)+argu2(i+1))
        r4Z=r4Z+dx*0.5d0*(argup(i)+argup(i+1))*radi**2.d0
        r2Z=r2Z+dx*0.5d0*(argup(i)+argup(i+1))
        r4N=r4N+dx*0.5d0*(argun(i)+argun(i+1))*radi**2.d0
        r2N=r2N+dx*0.5d0*(argun(i)+argun(i+1))
       enddo

       open(1,file='radial.out',status='unknown',form='formatted')
        write(1,*) 'sqrt(<r^2>_tot)=',dsqrt(r4tot/r2tot),'fm'
        write(1,*) 'sqrt(<r^2>_p)=',dsqrt(r4Z/r2Z),'fm'
        write(1,*) 'sqrt(<r^2>_n)=',dsqrt(r4N/r2N),'fm'
       close(1)
       deallocate(argu2,argup,argun)
       deallocate(rad_denp,rad_denn,rad_den2)
!********************************************************************************

        call HFB_energy

        E_diff1=E_HFB-E_orig1

        write(*,*) E_HFB,ferp,fern

        call hp_field(hp)
        call hn_field(hn)
        call Dp_field(Dp)
        call Dn_field(Dn)

        if(if_self.eq.1) call hY_field(hY) 

        H11p=Dp
        H11n=Dn
       enddo

       if(ifp_hfb.or.ifn_hfb) then
        i_count2=0
        E_diff2=1.d9

        do while(dabs(E_diff2).gt.precis.and.i_count2.le.count_max)
         i_count2=i_count2+1
         E_orig2=E_HFB

         rp=rhop_HFB
         rn=rhon_HFB
         kp=kapp_HFB
         kn=kapn_HFB

         if(if_self.eq.1) rY=rhoY_HFB

         if(.not.ifp_hfb) then
          call diagonalization(hp,wp,id)

          call make_sp_levels(hp,hn,wp,wn,i_count2+2)
          tran_p=hp
         endif

         if(.not.ifn_hfb) then
          call diagonalization(hn,wn,id)

          call make_sp_levels(hp,hn,wp,wn,i_count2+2)
          tran_n=hn
         endif

         if(ifp_hfb) then
          if(i_count2.eq.1) then
           do i=1,id
            do j=1,id
             if(levp(i)%j2.eq.levp(j)%j2.and.levp(i)%l.eq.levp(j)%l)
     &                                                            then
              Dp(i,j)=Dp(i,j)+0.1d0
             endif
            enddo
           enddo
          endif

          do i=1,id
           hp(i,i)=hp(i,i)-ferp
          enddo

!          call dgemm('N','N',id,id,id,1.d0,hp,max(1,id),hp,
!     &                        max(1,id),0.d0,h2p,max(1,id))
!          call dgemm('N','T',id,id,id,1.d0,Dp,max(1,id),Dp,
!     &                        max(1,id),0.d0,D2p,max(1,id))
!          do i=1,id
!           do j=1,id
!            h2p(i,j)=h2p(i,j)+D2p(i,j)
!           enddo
!          enddo

!          call diagonalization(h2p,w2p,id)
          call diagonalization(hp,wp,id)

          do i=1,id
           w2p(i)=Dp(i,i)
           w2n(i)=Dn(i,i)
          enddo

          call make_qsp_levels(wp,wn,hp,hn,w2p,w2n)
          tran_p=hp
         endif

         if(ifn_hfb) then
          if(i_count2.eq.1) then
           do i=1,id
            do j=1,id
             if(levn(i)%j2.eq.levn(j)%j2.and.levn(i)%l.eq.levn(j)%l)
     &                                                            then
              Dn(i,j)=Dn(i,j)+0.1d0
             endif
            enddo
           enddo
          endif

          do i=1,id
           hn(i,i)=hn(i,i)-fern
          enddo

!          call dgemm('N','N',id,id,id,1.d0,hn,max(1,id),hn,
!     &                        max(1,id),0.d0,h2n,max(1,id))
!          call dgemm('N','T',id,id,id,1.d0,Dn,max(1,id),Dn,
!     &                        max(1,id),0.d0,D2n,max(1,id))
!          do i=1,id
!           do j=1,id
!            h2n(i,j)=h2n(i,j)+D2n(i,j)
!           enddo
!          enddo

!          call diagonalization(h2n,w2n,id)
          call diagonalization(hn,wn,id)

          do i=1,id
           w2p(i)=Dp(i,i)
           w2n(i)=Dn(i,i)
          enddo

          call make_qsp_levels(wp,wn,hp,hn,w2p,w2n)
          tran_n=hn
         endif

         if(if_self.eq.1) call diagonalization(hY,wY,id)

         if(if_self.eq.1) then
          tran_Y=hY

          do i=1,id
           d_max=0.d0
           do j=1,id
            if(dabs(hY(j,i)).gt.d_max) then
             d_max=dabs(hY(j,i))
             j_pointer=j
            endif
           enddo
           jlY(i)=1000*levY(j_pointer)%l+levY(j_pointer)%j2
          enddo
          do i=1,id
           llll=jlY(i)/1000
           jjjj=mod(jlY(i),1000)
           lhfY(i)%ipar=(-1)**llll
           lhfY(i)%l=llll
           lhfY(i)%j2=jjjj
           lhfY(i)%ei=wY(i)
           lhfY(i)%qei=dabs(wY(i)-wY(1)+0.5d0)
          enddo
          do ll=0,(jmax+1)/2
           do jj=1,jmax,2
            nn=-1
            do i=1,id
             if(lhfY(i)%j2.eq.jj.and.lhfY(i)%l.eq.ll) then
              nn=nn+1
              lhfY(i)%nn=nn
             endif
            enddo
           enddo
          enddo
          do i=1,id
           do j=1,id
            if(levY(i)%nn.eq.lhfY(j)%nn.and.(levY(i)%l.eq.lhfY(j)%l
     &                           .and.levY(i)%j2.eq.lhfY(j)%j2)) then
             lhfY(j)%ui=levY(i)%ui
             lhfY(j)%vi=levY(i)%vi
            endif
           enddo
          enddo

          VY_HFB=0.d0
          UY_HFB=0.d0
          do i=1,id
           UY_HFB(i,i)=lhfY(i)%ui
           VY_HFB(i,i)=lhfY(i)%vi!*(-1)**((lhfp(i)%j2-1)/2)
          enddo
          call dgemm('N','N',id,id,id,1.d0,hY,max(1,id),VY_HFB,
     &                         max(1,id),0.d0,BY_HFB,max(1,id))
          VY_HFB=BY_HFB

          call dgemm('N','T',id,id,id,1.d0,VY_HFB,max(1,id),VY_HFB,
     &                           max(1,id),0.d0,rhoY_HFB,max(1,id))

          open(1,file='HF_Y.out',status='unknown',form='formatted')
           write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
           do ii = 1, id
            write(1,*) lhfY(ii)%index,lhfY(ii)%l,lhfY(ii)%j2,
     &       lhfY(ii)%ei,lhfY(ii)%qei,lhfY(ii)%vi
           end do
          close(1)
         endif

         call make_densities
         call check_densities

         bos1=dsqrt(0.5d0*(zmp+zmn)*hbarom/(hbarc**2.d0))
         bos2=dsqrt(zmY*hbarom/(hbarc**2.d0))
         rad_den1=0.d0
         do i=1,igrid
          val=0.d0
          radi=zcross(i)    !dble(i)*sizebox/dble(igrid)
          do j=1,id
           val_p=0.d0
           val_n=0.d0
           do k=1,id
            val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
            val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
           enddo
           val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
           val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
          enddo
          rad_den1(i)=val/(4.d0*pi)
         enddo
!********************************************************************************
       r4tot=0.d0
       r2tot=0.d0
       r4Z=0.d0
       r2Z=0.d0
       r4N=0.d0
       r2N=0.d0
       allocate(argu2(0:igrid2),argup(0:igrid2),argun(0:igrid2))
       argu2=0.d0
       argup=0.d0
       argun=0.d0
       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),
     &                                    rad_den2(0:igrid2))
       rad_denp=0.d0
       rad_denn=0.d0
       rad_den2=0.d0
        do i=0,igrid2
         val1=0.d0
         val2=0.d0
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val1=val1+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val2=val2+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_denp(i)=val1/(4.d0*pi)
         rad_denn(i)=val2/(4.d0*pi)
         rad_den2(i)=(val1+val2)/(4.d0*pi)
        enddo
       do i=0,igrid2
        radi=dble(i)*sizebox/dble(igrid2)
        argu2(i)=radi**2.d0*rad_den2(i)
        argup(i)=radi**2.d0*rad_denp(i)
        argun(i)=radi**2.d0*rad_denn(i)
       enddo
       dx=sizebox/dble(igrid2)
       do i=0,igrid2-1
        radi=(dble(i)+0.5d0)*sizebox/dble(igrid2)
        r4tot=r4tot+dx*0.5d0*(argu2(i)+argu2(i+1))*radi**2.d0
        r2tot=r2tot+dx*0.5d0*(argu2(i)+argu2(i+1))
        r4Z=r4Z+dx*0.5d0*(argup(i)+argup(i+1))*radi**2.d0
        r2Z=r2Z+dx*0.5d0*(argup(i)+argup(i+1))
        r4N=r4N+dx*0.5d0*(argun(i)+argun(i+1))*radi**2.d0
        r2N=r2N+dx*0.5d0*(argun(i)+argun(i+1))
       enddo

       open(1,file='radial.out',status='unknown',form='formatted')
        write(1,*) 'sqrt(<r^2>_tot)=',dsqrt(r4tot/r2tot),'fm'
        write(1,*) 'sqrt(<r^2>_p)=',dsqrt(r4Z/r2Z),'fm'
        write(1,*) 'sqrt(<r^2>_n)=',dsqrt(r4N/r2N),'fm'
       close(1)
       deallocate(argu2,argup,argun)
       deallocate(rad_denp,rad_denn,rad_den2)
!********************************************************************************

         rhop_HFB=alpha*rp+(1.d0-alpha)*rhop_HFB
         rhon_HFB=alpha*rn+(1.d0-alpha)*rhon_HFB
         kapp_HFB=alpha*kp+(1.d0-alpha)*kapp_HFB
         kapn_HFB=alpha*kn+(1.d0-alpha)*kapn_HFB

         if(if_self.eq.1) rhoY_HFB=alpha*rY+(1.d0-alpha)*rhoY_HFB

         call HFB_energy

         E_diff2=E_HFB-E_orig2

         write(*,*) E_HFB,ferp,fern

         call hp_field(hp)
         call hn_field(hn)
         call Dp_field(Dp)
         call Dn_field(Dn)

         if(if_self.eq.1) call hY_field(hY)

         H11p=Dp
         H11n=Dn
        enddo
       endif

       open(1,file='Energy.out',status='unknown',form='formatted')
        write(1,*) E_HFB,E_pair
       close(1)

       open(1,file='HF_basis_p.out',status='unknown',form='formatted')
       open(2,file='HF_basis_n.out',status='unknown',form='formatted')
        write(1,*) tran_p 
        write(2,*) tran_n
       close(1)
       close(2)

       do i=1,id
!        lhfp(i)%vi=lhfp(i)%vi*(-1)**((lhfp(i)%j2-1)/2)
!        lhfn(i)%vi=lhfn(i)%vi*(-1)**((lhfn(i)%j2-1)/2)
       enddo
       open(1,file='HF_p.out',status='unknown',form='formatted')
        write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
        do ii = 1, id
         write(1,*) lhfp(ii)%index,lhfp(ii)%l,lhfp(ii)%j2,
     &    lhfp(ii)%ei,lhfp(ii)%qei,lhfp(ii)%vi
        end do
       close(1)

       open(1,file='HF_n.out',status='unknown',form='formatted')
        write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
        do ii = 1, id
         write(1,*) lhfn(ii)%index,lhfn(ii)%l,lhfn(ii)%j2,
     &    lhfn(ii)%ei,lhfn(ii)%qei,lhfn(ii)%vi
        end do
       close(1)

       s1=0.d0
       s2=0.d0
       do i=1,id
        s1=s1+rhop_HFB(i,i)*dble(levp(i)%j2+1)
        s2=s2+rhon_HFB(i,i)*dble(levn(i)%j2+1)
       enddo

       if(dabs(s1-dble(AZ)).gt.dsqrt(precis)) then
        write(*,*) 'Error: wrong proton number AZ!!!'
        write(*,*) s1,'=/=',AZ
        stop
       endif
       if(dabs(s2-dble(AN)).gt.dsqrt(precis)) then
        write(*,*) 'Error: wrong proton number AN!!!'
        write(*,*) s2,'=/=',AN
        stop
       endif

       call check_densities

       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),
     &                                    rad_den2(0:igrid2))
       rad_denp=0.d0
       rad_denn=0.d0
       rad_den2=0.d0
        do i=0,igrid2
         val1=0.d0
         val2=0.d0
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val1=val1+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val2=val2+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_denp(i)=val1/(4.d0*pi)
         rad_denn(i)=val2/(4.d0*pi)
         rad_den2(i)=(val1+val2)/(4.d0*pi)
        enddo

       open(47,file='rad_density.out',status='unknown',form='formatted')
       rad_sum=0.d0
       dx=sizebox/dble(igrid2)
       do i=0,igrid2
        radi=dble(i)*sizebox/dble(igrid2)
        rad_sum=rad_sum+rad_den2(i)*dx*radi**2.d0
        write(47,*) radi,rad_den2(i)
       enddo
       write(47,*) 'Int rho = ',rad_sum*4.d0*pi
       close(47)

       deallocate(rad_denp,rad_denn,rad_den2)

!*********************************************************************************
!     Calculation of the Lambda s.p. basis                                       *

       if(if_self.eq.0) then

       hY=0.d0
       call hY_field(hY)

       call diagonalization(hY,wY,id)
       tran_Y=hY

       do i=1,id
        d_max=0.d0
        do j=1,id
         if(dabs(hY(j,i)).gt.d_max) then
          d_max=dabs(hY(j,i))
          j_pointer=j
         endif
        enddo
        jlY(i)=1000*levY(j_pointer)%l+levY(j_pointer)%j2
       enddo
       do i=1,id
        llll=jlY(i)/1000
        jjjj=mod(jlY(i),1000)
        lhfY(i)%ipar=(-1)**llll
        lhfY(i)%l=llll
        lhfY(i)%j2=jjjj
        lhfY(i)%ei=wY(i)
        lhfY(i)%qei=dabs(wY(i)-wY(1)+0.5d0)
        lhfY(i)%ui=1.d0
        lhfY(i)%vi=0.d0
       enddo

       open(1,file='HF_Y.out',status='unknown',form='formatted')
        write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
        do ii = 1, id
         write(1,*) lhfY(ii)%index,lhfY(ii)%l,lhfY(ii)%j2,
     &    lhfY(ii)%ei,lhfY(ii)%qei,lhfY(ii)%vi
        end do
       close(1)

       endif

!*********************************************************************************

       return
      end
