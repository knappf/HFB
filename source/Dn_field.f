      subroutine Dn_field(D1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: D1(id,id)

       type(twoquas_type), allocatable, save :: prho(:)
       integer :: iprho
       integer :: mpp(6),mtp(6)

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_gen(:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_gen(:,:,:,:)

       allocate(Vpp_gen(id,id,id,id),Vnn_gen(id,id,id,id))
       Vpp_gen=0.d0
       Vnn_gen=0.d0

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
       if (allocated(prho)) deallocate(prho)       
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

       do i=1,id
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
           if(Nk+Nl.le.noscmax12.and.(ji.eq.jj.and.jk.eq.jl)) then
           Jp=0 !max(abs(ji-jj)/2,abs(jk-jl)/2),min((ji+jj)/2,(jk+jl)/2)

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
!***********************
              enddo
             enddo
            enddo

             vv1=vv1+cv1*cv2*(Vppp*rp+Vppn*rn)
             vv2=vv2+cv1*cv2*(Vnnn*rn+Vnnp*rp)

             endif!N123

             enddo!mprho

             enddo!mk
            enddo!mi

           Vpp_gen(lp1(i),lp1(j),lp1(k),lp1(l))=
     &                        Vpp(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv1
           Vnn_gen(lp1(i),lp1(j),lp1(k),lp1(l))=
     &                        Vnn(lp1(i),lp1(j),lp1(k),lp1(l),Jp)+vv2


!           enddo!Jp
          endif!kl
          enddo!l
         enddo!k
        endif!ij
        enddo!j
       enddo!i

       D1=0.d0
       do i=1,id
        do j=1,id
         if(levn(i)%j2.eq.levn(j)%j2.and.levn(i)%l.eq.levn(j)%l) then
          val=0.d0
          do k=1,id
           do l=1,id
            val=val+Vnn_gen(k,l,i,j)
     &                *kapn_HFB(l,k)*dsqrt(dble(levn(k)%j2+1))
           enddo
          enddo
          D1(i,j)=0.5d0*val/dsqrt(dble(levn(i)%j2+1))
         endif
        enddo
       enddo

       deallocate(Vpp_gen,Vnn_gen)

       return
      end
