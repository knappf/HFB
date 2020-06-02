      Subroutine hY_field(h1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: h1(id,id)
       double precision :: h2(id,id)

       integer :: mpp(6),mtp(6)

       h1=0.d0
       h2=0.d0
       do i=1,id
        do j=1,id
         val=0.d0
         do k=1,id
          do l=1,id
           if(levY(i)%j2.eq.levY(j)%j2
     &                .and.levY(k)%j2.eq.levY(l)%j2) then
            do Jp=0,jmax
            val=val+VpY(k,i,l,j,Jp)
     &              *rhop_HFB(l,k)*dble(2*Jp+1)/dble(levY(i)%j2+1)
     &             +VnY(k,i,l,j,Jp)
     &              *rhon_HFB(l,k)*dble(2*Jp+1)/dble(levY(i)%j2+1)
            enddo
           endif
          enddo
         enddo
         h1(i,j)=kin_Y(i,j)+val
        enddo
       enddo

       do i=1,idm
        Ni=lev1pnm(i)%N
        li=lev1pnm(i)%l
        ji=lev1pnm(i)%j2
        mi=lev1pnm(i)%m2
        ii=lev1pnm(i)%jsch
        do j=1,idm
         Nj=lev1pnm(j)%N
         lj=lev1pnm(j)%l
         jj=lev1pnm(j)%j2
         mj=lev1pnm(j)%m2
         ij=lev1pnm(j)%jsch
         if(li.eq.lj.and.(ji.eq.jj.and.mi.eq.mj)) then
         if(mi.eq.1.and.mj.eq.1) then
          val=0.d0

          do k=1,idm
           Nk=lev1pnm(k)%N
           lk=lev1pnm(k)%l
           jk=lev1pnm(k)%j2
           mk=lev1pnm(k)%m2
           ik=lev1pnm(k)%jsch
           do m=1,idm 
            Nm=lev1pnm(m)%N
            lm=lev1pnm(m)%l
            jm=lev1pnm(m)%j2
            mm=lev1pnm(m)%m2
            im=lev1pnm(m)%jsch
            if(lk.eq.lm.and.(jk.eq.jm.and.mk.eq.mm)) then 

          do l=1,idm
           Nl=lev1pnm(l)%N
           ll=lev1pnm(l)%l
           jl=lev1pnm(l)%j2
           ml=lev1pnm(l)%m2
           il=lev1pnm(l)%jsch
           do n=1,idm
            Nn=lev1pnm(n)%N
            ln=lev1pnm(n)%l
            jn=lev1pnm(n)%j2
            mn=lev1pnm(n)%m2
            inn=lev1pnm(n)%jsch
            if(ll.eq.ln.and.(jl.eq.jn.and.ml.eq.mn)) then

            if(Ni+Nk.le.noscmax12.and.Nj+Nm.le.noscmax12) then
            if(Ni+Nk+Nl.le.noscmax123.and.Nj+Nm+Nn.le.noscmax123) then

            r1=rhop_HFB(lp1(im),lp1(ik))*rhop_HFB(lp1(inn),lp1(il))
            r2=rhon_HFB(lp1(im),lp1(ik))*rhon_HFB(lp1(inn),lp1(il))
            r3=rhop_HFB(lp1(im),lp1(ik))*rhon_HFB(lp1(inn),lp1(il))

            Vppl=0.d0
            Vnnl=0.d0
            Vpnl=0.d0

            val=val+0.5d0*Vppl*r1+0.5d0*Vnnl*r2+Vpnl*r3

            endif
            endif

            endif
           enddo
          enddo

            endif
           enddo
          enddo

          h2(ii,ij)=val

         endif
         endif
        enddo
       enddo

       do i=1,id
        do j=1,id
         h1(i,j)=h1(i,j)+h2(lp2(i),lp2(j))
        enddo
       enddo

       return
      end
