      subroutine transf_interaction
       
       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpn_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: VpY_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: VnY_HFB(:,:,:,:,:)

       double precision :: Tp(id,id),Tn(id,id),TY(id,id)
       double precision :: bp(id,id),bn(id,id),bY(id,id)
       double precision :: timef,timein

       type(twoquas_type), allocatable, save :: prho(:)

       integer :: point_p(id),point_n(id)
       integer :: pinv_p(id),pinv_n(id)

       integer :: mpp(6),mtp(6)

!       Vpp=Vpp+3.d0*Vpp_DD
!       Vnn=Vnn+3.d0*Vnn_DD
!       Vpn=Vpn+3.d0*Vpn_DD
!       VpY=VpY+3.d0*VpY_DD
!       VnY=VnY+3.d0*VnY_DD

       iprho=0
        do im=1,idm
         Nm=lev1pnm(im)%N
         lm=lev1pnm(im)%l
         jm=lev1pnm(im)%j2
         mm=lev1pnm(im)%m2
         m=lev1pnm(im)%jsch

         do inn=1,idm
          Nn=lev1pnm(inn)%N
          ln=lev1pnm(inn)%l
          jn=lev1pnm(inn)%j2
          mn=lev1pnm(inn)%m2
          n=lev1pnm(inn)%jsch
          if(ln.eq.lm.and.(jn.eq.jm.and.mn.eq.mm)) then
           iprho=iprho+1
          endif
         enddo
        enddo
       allocate(prho(iprho))
       iprho=0
        do im=1,idm
         Nm=lev1pnm(im)%N
         lm=lev1pnm(im)%l
         jm=lev1pnm(im)%j2
         mm=lev1pnm(im)%m2
         m=lev1pnm(im)%jsch

         do inn=1,idm
          Nn=lev1pnm(inn)%N
          ln=lev1pnm(inn)%l
          jn=lev1pnm(inn)%j2
          mn=lev1pnm(inn)%m2
          n=lev1pnm(inn)%jsch
          if(ln.eq.lm.and.(jn.eq.jm.and.mn.eq.mm)) then
           iprho=iprho+1
           prho(iprho)%q1=im
           prho(iprho)%q2=inn
          endif
         enddo
        enddo

       call cpu_time(timein)

       do i=1,id
        write(*,*) 'Index i=',i,'from',id
        Ni=lev1pn(i)%N
        li=lev1pn(i)%l
        ji=lev1pn(i)%j2
        do j=1,id
         Nj=lev1pn(j)%N
         lj=lev1pn(j)%l
         jj=lev1pn(j)%j2
         if(Ni+Nj.le.noscmax12) then
         do k=1,id
          Nk=lev1pn(k)%N
          lk=lev1pn(k)%l
          jk=lev1pn(k)%j2
          do l=1,id
           Nl=lev1pn(l)%N
           ll=lev1pn(l)%l
           jl=lev1pn(l)%j2
           if(Nk+Nl.le.noscmax12) then
           do Jp=max(abs(ji-jj)/2,abs(jk-jl)/2),min((ji+jj)/2,(jk+jl)/2)

            vv1=0.d0
            vv2=0.d0
            vv3=0.d0
            vv4=0.d0
            vv5=0.d0

            do mi=-min(ji,jj),min(ji,jj),2
             mj=-mi
             do mk=-min(jk,jl),min(jk,jl),2
              ml=-mk

              cv1=cg3(ji,mi,jj,mj,2*Jp)
              cv2=cg3(jk,mk,jl,ml,2*Jp)

            do mprho=1,iprho
             im=prho(mprho)%q1
             inn=prho(mprho)%q2

!            do im=1,idm 
             Nm=lev1pnm(im)%N
             lm=lev1pnm(im)%l
             jm=lev1pnm(im)%j2
             mm=lev1pnm(im)%m2
             m=lev1pnm(im)%jsch

!            do inn=1,idm
             Nn=lev1pnm(inn)%N
             ln=lev1pnm(inn)%l
             jn=lev1pnm(inn)%j2
             mn=lev1pnm(inn)%m2
             n=lev1pnm(inn)%jsch

             if(Ni+Nj+Nm.le.noscmax123.and.Nk+Nl+Nn.le.noscmax123) then

             rp=rhop_HFB(lp1(n),lp1(m))
             rn=rhon_HFB(lp1(n),lp1(m))
             rY=rhoY_HFB(lp1(n),lp1(m))

             larrow1%a=i
             larrow1%b=j
             larrow1%c=m
             larrow1%d=k
             larrow1%e=l
             larrow1%f=n
             iphase=1

             mpp(1)=mi
             mpp(2)=mj
             mpp(3)=mm
             mpp(4)=mk
             mpp(5)=ml
             mpp(6)=mn

             call order

             mtp(1)=1
             mtp(2)=1
             mtp(3)=1
             mtp(4)=1
             mtp(5)=1
             mtp(6)=1

             Vppp=0.d0
             Vnnp=0.d0
             Vppn=0.d0
             Vnnn=0.d0
             Vpnp=0.d0
             Vpnn=0.d0

             Vppl=0.d0
             Vnnl=0.d0
             Vpnl=0.d0
             Vplp=0.d0
             Vpln=0.d0
             Vnln=0.d0
             Vnlp=0.d0

            iip=larrow2%a
            jip=lev1pn(larrow2%a)%j2
            mip=mpp(larrowm%a)
            ijp=larrow2%b
            jjp=lev1pn(larrow2%b)%j2
            mjp=mpp(larrowm%b)
            imp=larrow2%c
            jmp=lev1pn(larrow2%c)%j2
            mmp=mpp(larrowm%c)
            ikp=larrow2%d
            jkp=lev1pn(larrow2%d)%j2
            mkp=mpp(larrowm%d)
            ilp=larrow2%e
            jlp=lev1pn(larrow2%e)%j2
            mlp=mpp(larrowm%e)
            inp=larrow2%f
            jnp=lev1pn(larrow2%f)%j2
            mnp=mpp(larrowm%f)

            do J12p=abs(jip-jjp)/2,(jip+jjp)/2
               c1=cg3(jip,mip,jjp,mjp,2*J12p)
             do J45p=abs(jkp-jlp)/2,(jkp+jlp)/2
               c2=cg3(jkp,mkp,jlp,mlp,2*J45p)
              do jtot2=max(abs(2*J12p-jmp),
     &           abs(2*J45p-jnp)),min(2*J12p+jmp,2*J45p+jnp),2
               c3=cg3(2*J12p,mip+mjp,jmp,mmp,jtot2)
               c4=cg3(2*J45p,mkp+mlp,jnp,mnp,jtot2)
!***********************
             mtp(1)=1
             mtp(2)=1
             mtp(3)=1
             mtp(4)=1
             mtp(5)=1
             mtp(6)=1

             mtip=mtp(larrowm%a)
             mtjp=mtp(larrowm%b)
             mtmp=mtp(larrowm%c)
             mtkp=mtp(larrowm%d)
             mtlp=mtp(larrowm%e)
             mtnp=mtp(larrowm%f)

            do itabp=0,1
               c5=cg3(1,mtip,1,mtjp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtkp,1,mtlp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtjp,1,mtmp,itot2)
               c8=cg3(2*itdep,mtkp+mtlp,1,mtnp,itot2)
               ipo1=lpoint(iip,ijp,imp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ikp,ilp,inp,J45p,jtot2,itdep,itot2)
               Vppp=Vppp+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo

             mtp(1)=1
             mtp(2)=1
             mtp(3)=-1
             mtp(4)=1
             mtp(5)=1
             mtp(6)=-1

            mtip=mtp(larrowm%a)
            mtjp=mtp(larrowm%b)
            mtmp=mtp(larrowm%c)
            mtkp=mtp(larrowm%d)
            mtlp=mtp(larrowm%e)
            mtnp=mtp(larrowm%f)

            do itabp=0,1
               c5=cg3(1,mtip,1,mtjp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtkp,1,mtlp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtjp,1,mtmp,itot2)
               c8=cg3(2*itdep,mtkp+mtlp,1,mtnp,itot2)
               ipo1=lpoint(iip,ijp,imp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ikp,ilp,inp,J45p,jtot2,itdep,itot2)
               Vppn=Vppn+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo

             mtp(1)=-1
             mtp(2)=-1
             mtp(3)=-1
             mtp(4)=-1
             mtp(5)=-1
             mtp(6)=-1

            mtip=mtp(larrowm%a)
            mtjp=mtp(larrowm%b)
            mtmp=mtp(larrowm%c)
            mtkp=mtp(larrowm%d)
            mtlp=mtp(larrowm%e)
            mtnp=mtp(larrowm%f)

            do itabp=0,1
               c5=cg3(1,mtip,1,mtjp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtkp,1,mtlp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtjp,1,mtmp,itot2)
               c8=cg3(2*itdep,mtkp+mtlp,1,mtnp,itot2)
               ipo1=lpoint(iip,ijp,imp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ikp,ilp,inp,J45p,jtot2,itdep,itot2)
               Vnnn=Vnnn+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo

             mtp(1)=-1
             mtp(2)=-1
             mtp(3)=1
             mtp(4)=-1
             mtp(5)=-1
             mtp(6)=1

            mtip=mtp(larrowm%a)
            mtjp=mtp(larrowm%b)
            mtmp=mtp(larrowm%c)
            mtkp=mtp(larrowm%d)
            mtlp=mtp(larrowm%e)
            mtnp=mtp(larrowm%f)

            do itabp=0,1
               c5=cg3(1,mtip,1,mtjp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtkp,1,mtlp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtjp,1,mtmp,itot2)
               c8=cg3(2*itdep,mtkp+mtlp,1,mtnp,itot2)
               ipo1=lpoint(iip,ijp,imp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ikp,ilp,inp,J45p,jtot2,itdep,itot2)
               Vnnp=Vnnp+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo

             mtp(1)=1
             mtp(2)=-1
             mtp(3)=1
             mtp(4)=1
             mtp(5)=-1
             mtp(6)=1

            mtip=mtp(larrowm%a)
            mtjp=mtp(larrowm%b)
            mtmp=mtp(larrowm%c)
            mtkp=mtp(larrowm%d)
            mtlp=mtp(larrowm%e)
            mtnp=mtp(larrowm%f)

            do itabp=0,1
               c5=cg3(1,mtip,1,mtjp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtkp,1,mtlp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtjp,1,mtmp,itot2)
               c8=cg3(2*itdep,mtkp+mtlp,1,mtnp,itot2)
               ipo1=lpoint(iip,ijp,imp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ikp,ilp,inp,J45p,jtot2,itdep,itot2)
               Vpnp=Vpnp+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo

             mtp(1)=1
             mtp(2)=-1
             mtp(3)=-1
             mtp(4)=1
             mtp(5)=-1
             mtp(6)=-1

            mtip=mtp(larrowm%a)
            mtjp=mtp(larrowm%b)
            mtmp=mtp(larrowm%c)
            mtkp=mtp(larrowm%d)
            mtlp=mtp(larrowm%e)
            mtnp=mtp(larrowm%f)

            do itabp=0,1
               c5=cg3(1,mtip,1,mtjp,2*itabp)
             do itdep=0,1
               c6=cg3(1,mtkp,1,mtlp,2*itdep)
              do itot2=max(abs(2*itabp-1),
     &              abs(2*itdep-1)),min(2*itabp+1,2*itdep+1),2
               c7=cg3(2*itabp,mtip+mtjp,1,mtmp,itot2)
               c8=cg3(2*itdep,mtkp+mtlp,1,mtnp,itot2)
               ipo1=lpoint(iip,ijp,imp,J12p,jtot2,itabp,itot2)
               ipo2=lpoint(ikp,ilp,inp,J45p,jtot2,itdep,itot2)
               Vpnn=Vpnn+c1*c2*c3*c4*c5*c6*c7*c8*V3B(ipo1,ipo2)
     &                                             *dble(iphase)
              enddo
             enddo
            enddo
!***********************
              enddo
             enddo
            enddo

             vv1=vv1+cv1*cv2*(Vppp*rp+Vppn*rn+Vppl*rY)
             vv2=vv2+cv1*cv2*(Vnnn*rn+Vnnp*rp+Vnnl*rY)
             vv3=vv3+cv1*cv2*(Vpnp*rp+Vpnn*rn+Vpnl*rY)
             vv4=vv4+cv1*cv2*(Vplp*rp+Vpln*rn)
             vv5=vv5+cv1*cv2*(Vnln*rn+Vnlp*rp)

             endif!N123

             enddo!mprho

             enddo!mk
            enddo!mi

           Vpp(lp1(i),lp1(j),lp1(k),lp1(l),Jp)=
     &                        Vpp(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv1
           Vnn(lp1(i),lp1(j),lp1(k),lp1(l),Jp)=
     &                        Vnn(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv2
           Vpn(lp1(i),lp1(j),lp1(k),lp1(l),Jp)=
     &                        Vpn(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv3

           VpY(lp1(i),lp1(j),lp1(k),lp1(l),Jp)=
     &                        VpY(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv4
           VnY(lp1(i),lp1(j),lp1(k),lp1(l),Jp)=
     &                        VnY(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv5

           enddo!Jp
          endif!kl
          enddo!l
         enddo!k
        endif!ij
        enddo!j
       enddo!i

       call cpu_time(timef)
       write(190,*)'time vers 9=',timef-timein

       deallocate(prho)

       Tp=tran_p
       Tn=tran_n
       TY=tran_Y

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE0_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE0_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trE0_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trE0_p=bp
       trE0_n=bn
       trE0_Y=bY

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE1_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE1_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trE1_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trE1_p=bp
       trE1_n=bn
       trE1_Y=bY

!************************
!       open(61,file='E1p.out',status='unknown',form='formatted')
!        do i=1,id
!         write(61,*) (trE1_p(i,j),j=1,id)
!        enddo
!       close(61)
!       open(62,file='E1n.out',status='unknown',form='formatted')
!        do i=1,id
!         write(62,*) (trE1_n(i,j),j=1,id)
!        enddo
!       close(62)
!************************

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE2_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE2_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trE2_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trE2_p=bp
       trE2_n=bn
       trE2_Y=bY

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE3_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE3_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trE3_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trE3_p=bp
       trE3_n=bn
       trE3_Y=bY

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trEN_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trEN_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trEN_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trEN_p=bp
       trEN_n=bn
       trEN_Y=bY

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trS1_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trS1_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trS1_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trS1_p=bp
       trS1_n=bn
       trS1_Y=bY

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trM1s_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trM1s_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trM1s_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trM1s_p=bp
       trM1s_n=bn
       trM1s_Y=bY

       bp=0.d0
       bn=0.d0
       bY=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         val3=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trM1l_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trM1l_n(k,l)*Tn(k,i)*Tn(l,j)
           val3=val3+trM1l_Y(k,l)*TY(k,i)*TY(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
         bY(i,j)=val3
        enddo
       enddo
       trM1l_p=bp
       trM1l_n=bn
       trM1l_Y=bY

       do m=0,igrid2
        bp=0.d0
        bn=0.d0
        do i=1,id
         do j=1,id
          val1=0.d0
          val2=0.d0
          do k=1,id
           do l=1,id
            val1=val1+trE1_p_dens(k,l,m)*Tp(k,i)*Tp(l,j)
            val2=val2+trE1_n_dens(k,l,m)*Tn(k,i)*Tn(l,j)
           enddo
          enddo
          bp(i,j)=val1
          bn(i,j)=val2
         enddo
        enddo
        do i=1,id
         do j=1,id
          trE1_p_dens(i,j,m)=bp(i,j)
          trE1_n_dens(i,j,m)=bn(i,j)
         enddo
        enddo
       enddo

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+H11p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+H11n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       H11p=bp
       H11n=bn

       do i=1,id
        do j=1,id
         H11p(i,j)=-(lhfp(i)%ui*lhfp(j)%vi+lhfp(i)%vi*lhfp(j)%ui)
     &                                                  *H11p(i,j)
         H11n(i,j)=-(lhfn(i)%ui*lhfn(j)%vi+lhfn(i)%vi*lhfn(j)%ui)
     &                                                  *H11n(i,j)
        enddo
       enddo
       do i=1,id
        H11p(i,i)=H11p(i,i)+(lhfp(i)%ui**2.d0-lhfp(i)%vi**2.d0)
     &                                        *(lhfp(i)%ei-ferp)
        H11n(i,i)=H11n(i,i)+(lhfn(i)%ui**2.d0-lhfn(i)%vi**2.d0)
     &                                        *(lhfn(i)%ei-fern)
       enddo

       open(4,file='H11p.dat',status='unknown',form='formatted')
       write(4,*) ferp
        do i1=1,id
         write(4,*) (H11p(i1,j1),j1=1,id) 
        enddo
       close(4)
       open(4,file='H11n.dat',status='unknown',form='formatted')
       write(4,*) fern
        do i1=1,id
         write(4,*) (H11n(i1,j1),j1=1,id) 
        enddo
       close(4)

       allocate(Vpp_HFB(id,id,id,id,0:jmax))
       allocate(Vnn_HFB(id,id,id,id,0:jmax))
       allocate(Vpn_HFB(id,id,id,id,0:jmax))
       allocate(VpY_HFB(id,id,id,id,0:jmax))
       allocate(VnY_HFB(id,id,id,id,0:jmax))
       Vpp_HFB=0.d0
       Vnn_HFB=0.d0
       Vpn_HFB=0.d0
       VpY_HFB=0.d0
       VnY_HFB=0.d0

       write(*,*) 'Transf. interaction in 1st index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            valpY=0.d0
            valnY=0.d0
            do i2=1,id
             valp=valp+Vpp(i2,j1,k1,l1,Jp)*Tp(i2,i1)
             valn=valn+Vnn(i2,j1,k1,l1,Jp)*Tn(i2,i1)
             val0=val0+Vpn(i2,j1,k1,l1,Jp)*Tp(i2,i1)
             valpY=valpY+VpY(i2,j1,k1,l1,Jp)*Tp(i2,i1)
             valnY=valnY+VnY(i2,j1,k1,l1,Jp)*Tn(i2,i1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
            VpY_HFB(i1,j1,k1,l1,Jp)=valpY
            VnY_HFB(i1,j1,k1,l1,Jp)=valnY
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB
       VpY=VpY_HFB
       VnY=VnY_HFB

       write(*,*) 'Transf. interaction in 2nd index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            valpY=0.d0
            valnY=0.d0
            do j2=1,id
             valp=valp+Vpp(i1,j2,k1,l1,Jp)*Tp(j2,j1)
             valn=valn+Vnn(i1,j2,k1,l1,Jp)*Tn(j2,j1)
             val0=val0+Vpn(i1,j2,k1,l1,Jp)*Tn(j2,j1)
             valpY=valpY+VpY(i1,j2,k1,l1,Jp)*TY(j2,j1)
             valnY=valnY+VnY(i1,j2,k1,l1,Jp)*TY(j2,j1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
            VpY_HFB(i1,j1,k1,l1,Jp)=valpY
            VnY_HFB(i1,j1,k1,l1,Jp)=valnY
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB
       VpY=VpY_HFB
       VnY=VnY_HFB

       write(*,*) 'Transf. interaction in 3rd index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            valpY=0.d0
            valnY=0.d0
            do k2=1,id
             valp=valp+Vpp(i1,j1,k2,l1,Jp)*Tp(k2,k1)
             valn=valn+Vnn(i1,j1,k2,l1,Jp)*Tn(k2,k1)
             val0=val0+Vpn(i1,j1,k2,l1,Jp)*Tp(k2,k1)
             valpY=valpY+VpY(i1,j1,k2,l1,Jp)*Tp(k2,k1)
             valnY=valnY+VnY(i1,j1,k2,l1,Jp)*Tn(k2,k1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
            VpY_HFB(i1,j1,k1,l1,Jp)=valpY
            VnY_HFB(i1,j1,k1,l1,Jp)=valnY
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB
       VpY=VpY_HFB
       VnY=VnY_HFB

       write(*,*) 'Transf. interaction in 4th index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            valpY=0.d0
            valnY=0.d0
            do l2=1,id
             valp=valp+Vpp(i1,j1,k1,l2,Jp)*Tp(l2,l1)
             valn=valn+Vnn(i1,j1,k1,l2,Jp)*Tn(l2,l1)
             val0=val0+Vpn(i1,j1,k1,l2,Jp)*Tn(l2,l1)
             valpY=valpY+VpY(i1,j1,k1,l2,Jp)*TY(l2,l1)
             valnY=valnY+VnY(i1,j1,k1,l2,Jp)*TY(l2,l1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
            VpY_HFB(i1,j1,k1,l1,Jp)=valpY
            VnY_HFB(i1,j1,k1,l1,Jp)=valnY
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB
       VpY=VpY_HFB
       VnY=VnY_HFB

       open(1,file='vlk_hfb.dat',status='unknown',form='formatted')
       write(1,*) '  Tz Par  2J   a   b   c   d          <ab|V|cd> '
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1-1
            j=2*i2-1
            k=2*i3-1
            l=2*i4-1
            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              factab=1.d0
              factcd=1.d0
              if(i.eq.j) factab=dsqrt(2.d0)
              if(k.eq.l) factcd=dsqrt(2.d0)
              xnorm=factab*factcd
              if(dabs(Vpp_HFB(i1,i2,i3,i4,Jp)).gt.precis) write(1,*)
     &        -1,0,2*Jp,i,j,k,l,Vpp_HFB(i1,i2,i3,i4,Jp)/xnorm
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1-1
            j=2*i2
            k=2*i3-1
            l=2*i4
!            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              if(dabs(Vpn_HFB(i1,i2,i3,i4,Jp)).gt.precis)
     &        write(1,*) 0,0,2*Jp,i,j,k,l,Vpn_HFB(i1,i2,i3,i4,Jp)
             endif
!            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1
            j=2*i2
            k=2*i3
            l=2*i4
            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              factab=1.d0
              factcd=1.d0
              if(i.eq.j) factab=dsqrt(2.d0)
              if(k.eq.l) factcd=dsqrt(2.d0)
              xnorm=factab*factcd
              if(dabs(Vnn_HFB(i1,i2,i3,i4,Jp)).gt.precis) write(1,*)
     &        1,0,2*Jp,i,j,k,l,Vnn_HFB(i1,i2,i3,i4,Jp)/xnorm
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       close(1)

       open(1,file='vlk_hfb_NY.dat',status='unknown',form='formatted')
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            if(dabs(VpY_HFB(i1,i2,i3,i4,Jp)).gt.precis) write(1,*) 
     &  'pY',i1,i2,i3,i4,Jp,VpY_HFB(i1,i2,i3,i4,Jp)
            if(dabs(VnY_HFB(i1,i2,i3,i4,Jp)).gt.precis) write(1,*) 
     &  'nY',i1,i2,i3,i4,Jp,VnY_HFB(i1,i2,i3,i4,Jp)
           enddo
          enddo
         enddo
        enddo
       enddo
       close(1)

!****************************************************************************
!      The many-body perturbation theory E(2) energy is calculated here     *
       ener2pt=0.d0
       ener2ptY=0.d0
       if((.not.ifp_hfb).and.(.not.ifn_hfb)) then
        Vpp=Vpp_HFB
        Vnn=Vnn_HFB
        Vpn=Vpn_HFB
        VpY=VpY_HFB
        VnY=VnY_HFB
        call MBPT_energy(ener2pt,ener2ptY)

        open(1,file='Energy_2.out',status='unknown',form='formatted')
         write(1,*) 'E^(2)=',ener2pt,'MeV    for nuclear core'
         write(1,*) 'E^(2)=',ener2ptY,'MeV    for Lambda'
        close(1)
       endif
!****************************************************************************

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=-1
         do i=1,id
          if(lhfp(i)%j2.eq.jj.and.lhfp(i)%l.eq.ll) then
           nn=nn+1
           lhfp(i)%nn=nn
          endif
         enddo
        enddo
       enddo

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=-1
         do i=1,id
          if(lhfn(i)%j2.eq.jj.and.lhfn(i)%l.eq.ll) then
           nn=nn+1
           lhfn(i)%nn=nn
          endif
         enddo
        enddo
       enddo

       i1p=0
       i1n=0
       do i=1,id
        if(lhfp(i)%vi**2.d0.gt.0.98) then
         i1p=i1p+1
         point_p(i)=id+1-i1p
         pinv_p(id+1-i1p)=i
        else
         point_p(i)=i-i1p
         pinv_p(i-i1p)=i
        endif
        if(lhfn(i)%vi**2.d0.gt.0.98) then
         i1n=i1n+1
         point_n(i)=id+1-i1n
         pinv_n(id+1-i1n)=i
        else
         point_n(i)=i-i1n
         pinv_n(i-i1n)=i
        endif
       enddo

       open(2,file='proton_HF.dat',status='unknown',form='formatted')
        write(2,'(1x,a61)')'n,     l,    2*j,      Tz,     qei,
     &  ui,      vi'
        itz=-1
        do ii = 1, id
         jj=ii !pinv_p(ii)
         write(2,'(1x,4(i4,1x),3(f12.5,1x))') lhfp(jj)%nn,lhfp(jj)%l,lhf
     &p(jj)%j2,itz,lhfp(jj)%qei,DSQRT(1.D0-lhfp(jj)%vi**2.d0),
     &lhfp(jj)%vi
        end do
       close(2)
       open(2,file='neutron_HF.dat',status='unknown',form='formatted')
        write(2,'(1x,a61)')'n,     l,    2*j,      Tz,     qei,
     &  ui,      vi'
        itz=1
        do ii = 1, id
         jj=ii !pinv_n(ii)
         write(2,'(1x,4(i4,1x),3(f12.5,1x))') lhfn(jj)%nn,lhfn(jj)%l,lhf
     &n(jj)%j2,itz,lhfn(jj)%qei,DSQRT(1.D0-lhfn(jj)%vi**2.d0),
     &lhfn(jj)%vi
        end do
       close(2)

!       open(3,file='gmat_HF_pp.dat',status='unknown',form='formatted')
!       itz=-2
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_p(i1)
!             j=i2 !pinv_p(i2)
!             k=i3 !pinv_p(i3)
!             l=i4 !pinv_p(i4)
!             if(i.le.j.and.k.le.l) then
!              if(1000*i+j.le.1000*k+l) then
!               if(dabs(Vpp_HFB(i,j,k,l,Jp)).gt.precis) then
!                factab=1.d0
!                factcd=1.d0
!                if(i.eq.j) factab=dsqrt(2.d0)
!                if(k.eq.l) factcd=dsqrt(2.d0)
!                xnorm=factab*factcd
!                write(3,'(1x,14(i4,1x),1(f12.5,1x))') 
!     &            lhfp(i)%nn,lhfp(i)%l,lhfp(i)%j2,
!     &            lhfp(j)%nn,lhfp(j)%l,lhfp(j)%j2,
!     &            lhfp(k)%nn,lhfp(k)%l,lhfp(k)%j2,
!     &            lhfp(l)%nn,lhfp(l)%l,lhfp(l)%j2,
!     &            itz,2*Jp,Vpp_HFB(i,j,k,l,Jp)/xnorm
!               endif
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

!       open(3,file='gmat_HF_nn.dat',status='unknown',form='formatted')
!       itz=2
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_n(i1)
!             j=i2 !pinv_n(i2)
!             k=i3 !pinv_n(i3)
!             l=i4 !pinv_n(i4)
!             if(i.le.j.and.k.le.l) then
!              if(1000*i+j.le.1000*k+l) then
!               if(dabs(Vnn_HFB(i,j,k,l,Jp)).gt.precis) then
!                factab=1.d0
!                factcd=1.d0
!                if(i.eq.j) factab=dsqrt(2.d0)
!                if(k.eq.l) factcd=dsqrt(2.d0)
!                xnorm=factab*factcd
!                write(3,'(1x,14(i4,1x),1(f12.5,1x))')
!     &            lhfn(i)%nn,lhfn(i)%l,lhfn(i)%j2,
!     &            lhfn(j)%nn,lhfn(j)%l,lhfn(j)%j2,
!     &            lhfn(k)%nn,lhfn(k)%l,lhfn(k)%j2,
!     &            lhfn(l)%nn,lhfn(l)%l,lhfn(l)%j2,
!     &            itz,2*Jp,Vnn_HFB(i,j,k,l,Jp)/xnorm
!               endif
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

!       open(3,file='gmat_HF_pn.dat',status='unknown',form='formatted')
!       itz=0
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_p(i1)
!             j=i2 !pinv_n(i2)
!             k=i3 !pinv_p(i3)
!             l=i4 !pinv_n(i4)
!             if(1000*i+j.le.1000*k+l) then
!              if(dabs(Vpn_HFB(i,j,k,l,Jp)).gt.precis) then
!               write(3,'(1x,14(i4,1x),1(f12.5,1x))')
!     &           lhfp(i)%nn,lhfp(i)%l,lhfp(i)%j2,
!     &           lhfn(j)%nn,lhfn(j)%l,lhfn(j)%j2,
!     &           lhfp(k)%nn,lhfp(k)%l,lhfp(k)%j2,
!     &           lhfn(l)%nn,lhfn(l)%l,lhfn(l)%j2,
!     &           itz,2*Jp,Vpn_HFB(i,j,k,l,Jp)
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

       open(4,file='unitar_p.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (Tp(i1,j1),j1=1,id) !(Tp(i1,pinv_p(j1)),j1=1,id)
        enddo
       close(4)

       open(4,file='unitar_n.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (Tn(i1,j1),j1=1,id) !(Tn(i1,pinv_n(j1)),j1=1,id)
        enddo
       close(4)

       open(4,file='unitar_Y.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (TY(i1,j1),j1=1,id) 
        enddo
       close(4)

       deallocate(Vpp_HFB,Vnn_HFB,Vpn_HFB,VpY_HFB,VnY_HFB)

       return
      end
