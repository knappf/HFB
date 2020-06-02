      Subroutine hn_field(h1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
   
       double precision :: h1(id,id)
       double precision :: h2(id,id)
       type(twoquas_type), allocatable, save :: prho(:)
       integer :: iprho
       integer :: mpp(6),mtp(6)

       h1=0.d0
       h2=0.d0
       do i=1,id
        do j=1,id
         val=0.d0
         do k=1,id
          do l=1,id
           if(levn(i)%j2.eq.levn(j)%j2
     &                .and.levn(k)%j2.eq.levn(l)%j2) then
            do Jp=0,jmax
            val=val+Vnn(i,k,j,l,Jp)
     &              *rhon_HFB(l,k)*dble(2*Jp+1)/dble(levn(i)%j2+1)
     &             +Vpn(k,i,l,j,Jp)
     &              *rhop_HFB(l,k)*dble(2*Jp+1)/dble(levn(i)%j2+1)
            if(if_self.eq.1) then
             val=val+VnY(i,k,j,l,Jp)
     & *rhoY_HFB(l,k)*dble(2*Jp+1)/(dble((levn(i)%j2+1)*(levY(k)%j2+1)))
            endif
            enddo
           endif
          enddo
         enddo
         h1(i,j)=kin_n(i,j)+val
        enddo
       enddo

       iprho=0
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
         iprho=iprho+1
         !prho(iprho)%q1=i
         !prho(iprho)%q2=j
         end if
       end do
       end do
       if (allocated(prho)) deallocate(prho)
       allocate(prho(iprho))
       iprho=0
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
         iprho=iprho+1
         prho(iprho)%q1=i
         prho(iprho)%q2=j
         end if
       end do
       end do


       do indij=1,iprho 
        i=prho(indij)%q1
        j=prho(indij)%q2

        Ni=lev1pnm(i)%N
        li=lev1pnm(i)%l
        ji=lev1pnm(i)%j2
        mi=lev1pnm(i)%m2
        ii=lev1pnm(i)%jsch

         Nj=lev1pnm(j)%N
         lj=lev1pnm(j)%l
         jj=lev1pnm(j)%j2
         mj=lev1pnm(j)%m2
         ij=lev1pnm(j)%jsch

         if(mi.eq.1.and.mj.eq.1) then
         val=0.d0

       do indkm=1,iprho
        k=prho(indkm)%q1
        m=prho(indkm)%q2

           Nk=lev1pnm(k)%N
           lk=lev1pnm(k)%l
           jk=lev1pnm(k)%j2
           mk=lev1pnm(k)%m2
           ik=lev1pnm(k)%jsch

            Nm=lev1pnm(m)%N
            lm=lev1pnm(m)%l
            jm=lev1pnm(m)%j2
            mm=lev1pnm(m)%m2
            im=lev1pnm(m)%jsch

       if(Ni+Nk.le.noscmax12.and.Nj+Nm.le.noscmax12) then

       do indln=1,iprho
        l=prho(indln)%q1
        n=prho(indln)%q2

           Nl=lev1pnm(l)%N
           ll=lev1pnm(l)%l
           jl=lev1pnm(l)%j2
           ml=lev1pnm(l)%m2
           il=lev1pnm(l)%jsch

            Nn=lev1pnm(n)%N
            ln=lev1pnm(n)%l
            jn=lev1pnm(n)%j2
            mn=lev1pnm(n)%m2
            inn=lev1pnm(n)%jsch

            if(Ni+Nk+Nl.le.noscmax123.and.Nj+Nm+Nn.le.noscmax123) then

            r1=rhon_HFB(lp1(im),lp1(ik))*rhon_HFB(lp1(inn),lp1(il))
            r2=rhop_HFB(lp1(im),lp1(ik))*rhop_HFB(lp1(inn),lp1(il))
            r3=rhop_HFB(lp1(im),lp1(ik))*rhon_HFB(lp1(inn),lp1(il))
            r4=rhon_HFB(lp1(im),lp1(ik))*rhoY_HFB(lp1(inn),lp1(il))
            r5=rhop_HFB(lp1(im),lp1(ik))*rhoY_HFB(lp1(inn),lp1(il))

            larrow1%a=ii
            larrow1%b=ik
            larrow1%c=il
            larrow1%d=ij
            larrow1%e=im
            larrow1%f=inn
            iphase=1

            mpp(1)=mi
            mpp(2)=mk
            mpp(3)=ml
            mpp(4)=mj
            mpp(5)=mm
            mpp(6)=mn

            call order

            mtp(1)=-1
            mtp(2)=-1
            mtp(3)=-1
            mtp(4)=-1
            mtp(5)=-1
            mtp(6)=-1

            Vnnn=0.d0
            Vnpp=0.d0
            Vnpn=0.d0
            Vnnl=0.d0
            Vpnl=0.d0

            iip=larrow2%a
            jip=lev1pn(larrow2%a)%j2
            mip=mpp(larrowm%a)
            ikp=larrow2%b
            jkp=lev1pn(larrow2%b)%j2
            mkp=mpp(larrowm%b)
            ilp=larrow2%c
            jlp=lev1pn(larrow2%c)%j2
            mlp=mpp(larrowm%c)
            ijp=larrow2%d
            jjp=lev1pn(larrow2%d)%j2
            mjp=mpp(larrowm%d)
            imp=larrow2%e
            jmp=lev1pn(larrow2%e)%j2
            mmp=mpp(larrowm%e)
            inp=larrow2%f
            jnp=lev1pn(larrow2%f)%j2
            mnp=mpp(larrowm%f)

            do J12p=abs(jip-jkp)/2,(jip+jkp)/2
               c1=cg3(jip,mip,jkp,mkp,2*J12p)
             do J45p=abs(jjp-jmp)/2,(jjp+jmp)/2
               c2=cg3(jjp,mjp,jmp,mmp,2*J45p)
              do jtot2=max(abs(2*J12p-jlp),
     &           abs(2*J45p-jnp)),min(2*J12p+jlp,2*J45p+jnp),2
               c3=cg3(2*J12p,mip+mkp,jlp,mlp,jtot2)
               c4=cg3(2*J45p,mjp+mmp,jnp,mnp,jtot2)
!******************************

            mtp(1)=-1
            mtp(2)=-1
            mtp(3)=-1
            mtp(4)=-1
            mtp(5)=-1
            mtp(6)=-1

            mtip=mtp(larrowm%a)
            mtkp=mtp(larrowm%b)
            mtlp=mtp(larrowm%c)
            mtjp=mtp(larrowm%d)
            mtmp=mtp(larrowm%e)
            mtnp=mtp(larrowm%f)
            do itabp=0,1
               c5=cg3(1,mtip,1,mtkp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtjp,1,mtmp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtkp,1,mtlp,itot2)
               c8=cg3(2*itdep,mtjp+mtmp,1,mtnp,itot2)

               ipo1=lpoint(iip,ikp,ilp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ijp,imp,inp,J45p,jtot2,itdep,itot2)
               Vnnn=Vnnn+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo
!****************************
               mtp(1)=-1
               mtp(2)=1
               mtp(3)=1
               mtp(4)=-1
               mtp(5)=1
               mtp(6)=1

               mtip=mtp(larrowm%a)
               mtkp=mtp(larrowm%b)
               mtlp=mtp(larrowm%c)
               mtjp=mtp(larrowm%d)
               mtmp=mtp(larrowm%e)
               mtnp=mtp(larrowm%f)

       do itabp=0,1
         c5=cg3(1,mtip,1,mtkp,2*itabp)
        do itdep=0,1
         c6=cg3(1,mtjp,1,mtmp,2*itdep)
         do itot2=max(abs(2*itabp-1),
     &            abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
         c7=cg3(2*itabp,mtip+mtkp,1,mtlp,itot2)
         c8=cg3(2*itdep,mtjp+mtmp,1,mtnp,itot2)
         ipo1=lpoint(iip,ikp,ilp,J12p,jtot2,itabp,itot2)
         ipo2=lpoint(ijp,imp,inp,J45p,jtot2,itdep,itot2)
         Vnpp=Vnpp+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
          enddo
         enddo
        enddo
!****************************

        mtp(1)=-1
        mtp(2)=1
        mtp(3)=-1
        mtp(4)=-1
        mtp(5)=1
        mtp(6)=-1

        mtip=mtp(larrowm%a)
        mtkp=mtp(larrowm%b)
        mtlp=mtp(larrowm%c)
        mtjp=mtp(larrowm%d)
        mtmp=mtp(larrowm%e)
        mtnp=mtp(larrowm%f)

        do itabp=0,1
        c5=cg3(1,mtip,1,mtkp,2*itabp)
         do itdep=0,1
        c6=cg3(1,mtjp,1,mtmp,2*itdep)
          do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
        c7=cg3(2*itabp,mtip+mtkp,1,mtlp,itot2)
        c8=cg3(2*itdep,mtjp+mtmp,1,mtnp,itot2)
        ipo1=lpoint(iip,ikp,ilp,J12p,jtot2,itabp,itot2)
        ipo2=lpoint(ijp,imp,inp,J45p,jtot2,itdep,itot2)
        Vnpn=Vnpn+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
          enddo
         enddo
        enddo
!****************************
              enddo
             enddo
            enddo

            val=val+0.5d0*Vnnn*r1+0.5d0*Vnpp*r2+Vnpn*r3+Vnnl*r4+Vpnl*r5

             endif
           enddo

             end if
           enddo

          h2(ii,ij)=val

         endif
       enddo

       do i=1,id
        do j=1,id
         h1(i,j)=h1(i,j)+h2(lp2(i),lp2(j))
        enddo
       enddo

       return
      end
