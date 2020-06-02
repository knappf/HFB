      subroutine initialize_basis

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       allocate(levp(id),levn(id),levY(id))
       allocate(lhfp(id),lhfn(id),lhfY(id))
       allocate(lev1pn(id),lp1(id),lp2(id))

       im = 0

       do N = 0,noscmax
        do l = N,0,-1
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
           do jj = 2*l+1, 2*l-1, -2
            if(jj.gt.0) then
             im = im + 1
             levp(im)%index = im
             levn(im)%index = im
             levY(im)%index = im
             levp(im)%ipar = (-1)**l
             levn(im)%ipar = (-1)**l
             levY(im)%ipar = (-1)**l
             levp(im)%N = N
             levn(im)%N = N
             levY(im)%N = N
             levp(im)%nn = nn
             levn(im)%nn = nn
             levY(im)%nn = nn
             levp(im)%l = l
             levn(im)%l = l
             levY(im)%l = l
             levp(im)%j2 = jj
             levn(im)%j2 = jj
             levY(im)%j2 = jj
             levp(im)%spenrg = hbarom * (dble(N)+1.5d0)
             levn(im)%spenrg = hbarom * (dble(N)+1.5d0)
             levY(im)%spenrg = hbarom * (dble(N)+1.5d0)
             levp(im)%ei = levp(im)%spenrg
             levn(im)%ei = levn(im)%spenrg
             levY(im)%ei = levY(im)%spenrg
             levp(im)%qei = levp(im)%spenrg
             levn(im)%qei = levn(im)%spenrg
             levY(im)%qei = levY(im)%spenrg
             levY(im)%ui = 1.d0
             levY(im)%vi = 0.d0
            endif
           end do
         endif
        end do
       end do

       im = 0

       do N = 0,noscmax
        do l = 0,N
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
           do jj = 2*l-1, 2*l+1, 2
            if(jj.gt.0) then
             im = im + 1
             lev1pn(im)%index = im
             lev1pn(im)%ipar = (-1)**l
             lev1pn(im)%N = N
             lev1pn(im)%nn = nn
             lev1pn(im)%l = l
             lev1pn(im)%j2 = jj
             lev1pn(im)%spenrg = hbarom * (dble(N)+1.5d0)
             lev1pn(im)%ei = lev1pn(im)%spenrg
             lev1pn(im)%qei = lev1pn(im)%spenrg
            endif
           end do
         endif
        end do
       end do

       im = 0
       imm= 0

       do N = 0,noscmax
        do l = 0,N
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
           do jj = 2*l-1, 2*l+1, 2
            if(jj.gt.0) then
             im = im + 1
             do mm = -jj,jj,2
              imm = imm + 1 
             enddo
            endif
           end do
         endif
        end do
       end do

       idm=imm
       allocate(lev1pnm(idm))

       im = 0
       imm= 0

       do N = 0,noscmax
        do l = 0,N
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
           do jj = 2*l-1, 2*l+1, 2
            if(jj.gt.0) then
             im = im + 1
             do mm = -jj,jj,2
              imm = imm + 1
              lev1pnm(imm)%index = imm
              lev1pnm(imm)%ipar = (-1)**l
              lev1pnm(imm)%N = N
              lev1pnm(imm)%nn = nn
              lev1pnm(imm)%l = l
              lev1pnm(imm)%j2 = jj
              lev1pnm(imm)%m2 = mm
              lev1pnm(imm)%jsch = im
              lev1pnm(imm)%spenrg = hbarom * (dble(N)+1.5d0)
              lev1pnm(imm)%ei = lev1pn(im)%spenrg
              lev1pnm(imm)%qei = lev1pn(im)%spenrg
             enddo
            endif
           end do
         endif
        end do
       end do

       lp1=0
       lp2=0

       do i=1,id
        do j=1,id
         if(levp(i)%nn.eq.lev1pn(j)%nn.and.(levp(i)%l.eq.lev1pn(j)%l
     &                         .and.levp(i)%j2.eq.lev1pn(j)%j2)) then
          lp1(j)=i
          lp2(i)=j
         endif
        enddo
       enddo

       jmax=0
       do i=1,id
        if(levp(i)%j2.gt.jmax) jmax=levp(i)%j2
       enddo

       i_lnl=0
       do N = 0,noscmax
        do l = N,0,-1
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
          i_lnl=i_lnl+1
         endif
        enddo
       enddo

       id_lnl=i_lnl

       i_lnl=0
       do N = 0,2*noscmax
        do l = N,0,-1
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
          i_lnl=i_lnl+1
         endif
        enddo
       enddo

       allocate(lnl(i_lnl))

       i_lnl=0
       do N = 0,2*noscmax
        do l = N,0,-1
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
          i_lnl=i_lnl+1
          lnl(i_lnl)%index = i_lnl
          lnl(i_lnl)%ipar = (-1)**l
          lnl(i_lnl)%N = N
          lnl(i_lnl)%nn = nn
          lnl(i_lnl)%l = l
          lnl(i_lnl)%spenrg = hbarom * (dble(N)+1.5d0)
         endif
        enddo
       enddo

       id_lnl2=i_lnl

       id3=0

       do i=1,id
        Ni=lev1pn(i)%N
        if(Ni.le.noscmax) then
        do j=1,i !id
         Nj=lev1pn(j)%N
         if(Ni+Nj.le.noscmax12) then
         do k=1,j !id
          Nk=lev1pn(k)%N
          if(Ni+Nj+Nk.le.noscmax123) then

           do J12=abs(lev1pn(i)%j2-lev1pn(j)%j2)/2,
     &                         (lev1pn(i)%j2+lev1pn(j)%j2)/2
            do Jtot=abs(2*J12-lev1pn(k)%j2),2*J12+lev1pn(k)%j2,2
           do itab=0,1
            do itot=abs(2*itab-1),2*itab+1,2 !1,3,2
             id3=id3+1
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

       allocate(lev3(id3))
       allocate(lpoint(id,id,id,0:jmax,1:3*jmax,0:1,3))
       lpoint=0

       id3=0

       do i=1,id
        Ni=lev1pn(i)%N
        if(Ni.le.noscmax) then
        do j=1,i !id
         Nj=lev1pn(j)%N
         if(Ni+Nj.le.noscmax12) then
         do k=1,j !id
          Nk=lev1pn(k)%N
          if(Ni+Nj+Nk.le.noscmax123) then

           do J12=abs(lev1pn(i)%j2-lev1pn(j)%j2)/2,
     &                         (lev1pn(i)%j2+lev1pn(j)%j2)/2
            do Jtot=abs(2*J12-lev1pn(k)%j2),2*J12+lev1pn(k)%j2,2
           do itab=0,1
            do itot=abs(2*itab-1),2*itab+1,2 !1,3,2
             id3=id3+1
             lev3(id3)%index=id3
             lev3(id3)%i=i
             lev3(id3)%j=j
             lev3(id3)%k=k
             lev3(id3)%Jab=J12
             lev3(id3)%JJ=Jtot
             lev3(id3)%Tab=itab
             lev3(id3)%TT=itot

             lpoint(i,j,k,J12,Jtot,itab,itot)=id3
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

       open(1,file='3b_basis.out',status='unknown',form='formatted')
        write(1,*) 'Total dimension = ',id3
        write(1,*)'n       i       j        k      Jab       Jtot
     &       Tab         Ttot'
        do ii = 1, id3
         write(1,*) lev3(ii)%index,lev3(ii)%i,lev3(ii)%j,
     &    lev3(ii)%k,lev3(ii)%Jab,lev3(ii)%JJ,lev3(ii)%Tab,lev3(ii)%TT
        end do
       close(1)

       if(n_lam_occup.lt.1.or.n_lam_occup.gt.id) then
        write(*,*) 'Error!!!'
        write(*,*) 'Value of n_lam_occup must be in the interval 1-id '
        stop
       endif

      return
      end
