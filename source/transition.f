      subroutine transition 

       USE technical
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       bos=dsqrt(0.5d0*(zmp+zmn)*hbarom/(hbarc**2.d0))

       allocate(trE0_p(id,id),trE0_n(id,id),trE0_Y(id,id))
       allocate(trE1_p(id,id),trE1_n(id,id),trE1_Y(id,id))
       allocate(trE2_p(id,id),trE2_n(id,id),trE2_Y(id,id))
       allocate(trE3_p(id,id),trE3_n(id,id),trE3_Y(id,id))
       allocate(trEN_p(id,id),trEN_n(id,id),trEN_Y(id,id))
       allocate(trS1_p(id,id),trS1_n(id,id),trS1_Y(id,id))
       allocate(trM1s_p(id,id),trM1s_n(id,id),trM1s_Y(id,id))
       allocate(trM1l_p(id,id),trM1l_n(id,id),trM1l_Y(id,id))
       trE0_p=0.d0
       trE0_n=0.d0
       trE0_Y=0.d0
       trE1_p=0.d0
       trE1_n=0.d0
       trE1_Y=0.d0
       trE2_p=0.d0
       trE2_n=0.d0
       trE2_Y=0.d0
       trE3_p=0.d0
       trE3_n=0.d0
       trE3_Y=0.d0
       trEN_p=0.d0
       trEN_n=0.d0
       trEN_Y=0.d0
       trS1_p=0.d0
       trS1_n=0.d0
       trS1_Y=0.d0
       trM1s_p=0.d0
       trM1s_n=0.d0
       trM1s_Y=0.d0
       trM1l_p=0.d0
       trM1l_n=0.d0
       trM1l_Y=0.d0

       do i=1,id
        do j=1,id
         lam=0
         radial=0.d0
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         radial=gauss_int(2,ni,li,nj,lj,bos)
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         trE0_p(i,j)=angular*radial
         trE0_n(i,j)=angular*radial
         trE0_Y(i,j)=angular*radial
        enddo
       enddo

       do i=1,id
        do j=1,id
         lam=1
         radial=0.d0
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         radial=gauss_int(1,ni,li,nj,lj,bos)
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         trE1_p(i,j)=angular*radial
         trE1_n(i,j)=angular*radial
         trE1_Y(i,j)=angular*radial
        enddo
       enddo

       do i=1,id
        do j=1,id
         lam=2
         radial=0.d0
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         radial=gauss_int(2,ni,li,nj,lj,bos)
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         trE2_p(i,j)=angular*radial
         trE2_n(i,j)=angular*radial
         trE2_Y(i,j)=angular*radial
        enddo
       enddo

       do i=1,id
        do j=1,id
         lam=3
         radial=0.d0
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         radial=gauss_int(3,ni,li,nj,lj,bos)
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         trE3_p(i,j)=angular*radial
         trE3_n(i,j)=angular*radial
         trE3_Y(i,j)=angular*radial
        enddo
       enddo

       do i=1,id
        do j=1,id
         lam=0
         radial=0.d0
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         if(ni.eq.nj.and.li.eq.lj) radial=1.d0
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         trEN_p(i,j)=angular*radial
         trEN_n(i,j)=angular*radial
         trEN_Y(i,j)=angular*radial
        enddo
       enddo

       do i=1,id
        do j=1,id
         lam=1
         radial=0.d0
         radial2=0.d0
         r2=1.44d0*(dble(AZ+AN)**(2.d0/3.d0))
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         radial=gauss_int(3,ni,li,nj,lj,bos)
         radial2=gauss_int(1,ni,li,nj,lj,bos)
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         trS1_p(i,j)=angular*(radial-r2*radial2)
         trS1_n(i,j)=angular*(radial-r2*radial2)
         trS1_Y(i,j)=angular*(radial-r2*radial2)
        enddo
       enddo

       do i=1,id
        do j=1,id
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         radial=0.d0
         if(ni.eq.nj.and.li.eq.lj) radial=1.d0
         angulars=0.d0
         angularl=0.d0
         cl1=cleb(2*li,0,1,1,ji,1)
         cl2=cleb(2*lj,0,1,1,jj,1)
         cl3=cleb(2*li,2,1,-1,ji,1)
         cl4=cleb(2*lj,2,1,-1,jj,1)
         cl5=cleb(2,0,jj,1,ji,1)
         if(dabs(cl5).lt.precis) cl5=1.d24
         phase=dble((-1)**((ji-jj)/2+1))
         angulars=0.5d0*(cl1*cl2-cl3*cl4)
         angularl=cl3*cl4
         trM1s_p(i,j)=dsqrt(3.d0)*phase*dsqrt(dble(ji+1))
     &                       *radial*angulars/(dsqrt(4.d0*pi)*cl5)
         trM1s_n(i,j)=dsqrt(3.d0)*phase*dsqrt(dble(ji+1))
     &                       *radial*angulars/(dsqrt(4.d0*pi)*cl5)
         trM1s_Y(i,j)=dsqrt(3.d0)*phase*dsqrt(dble(ji+1))
     &                       *radial*angulars/(dsqrt(4.d0*pi)*cl5)
         trM1l_p(i,j)=dsqrt(3.d0)*phase*dsqrt(dble(ji+1))
     &                       *radial*angularl/(dsqrt(4.d0*pi)*cl5)
         trM1l_n(i,j)=dsqrt(3.d0)*phase*dsqrt(dble(ji+1))
     &                       *radial*angularl/(dsqrt(4.d0*pi)*cl5)
         trM1l_Y(i,j)=dsqrt(3.d0)*phase*dsqrt(dble(ji+1))
     &                       *radial*angularl/(dsqrt(4.d0*pi)*cl5)
        enddo
       enddo

       allocate(trE1_p_dens(id,id,0:igrid2),trE1_n_dens(id,id,0:igrid2))
       trE1_p_dens=0.d0
       trE1_n_dens=0.d0

       do i=1,id
        do j=1,id
         lam=1
         li=levp(i)%l
         lj=levp(j)%l
         ji=levp(i)%j2
         jj=levp(j)%j2
         ni=levp(i)%nn
         nj=levp(j)%nn
         angular=0.d0
         if(mod(lam+li+lj,2).eq.0) then
          phase=dble((-1)**(li+lj+lam))
          ff1=dsqrt(dble((2*lam+1)*(jj+1)))/dsqrt(4.d0*pi)
          ff2=cleb(2*lam,0,jj,1,ji,1)
          angular=phase*ff1*ff2
         endif
         do k=0,igrid2
          radial=0.d0
          radi=dble(k)*sizebox/dble(igrid2)
          radial=radi**3.d0*R_val(ni,li,ji,bos,radi)
     &                           *R_val(nj,lj,jj,bos,radi)
          trE1_p_dens(i,j,k)=angular*radial
          trE1_n_dens(i,j,k)=angular*radial
         enddo
        enddo
       enddo

       return
      end
