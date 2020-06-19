      subroutine TDAtzpp

       USE technical
       USE math
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       type(twoquas_type), allocatable, save :: ph(:)
!       type(phonon_type), allocatable, save :: list_phon(:)
       double precision, allocatable, save :: TDA_matrix(:,:),w(:)
       double precision, allocatable, save :: spur_CM(:)
!       double precision, allocatable, save :: sixj1(:,:,:,:,:,:)
       integer, allocatable, save :: ityp1(:),ityp2(:)
       integer :: p1,h1,tz1,p2,h2,tz2

       idtz = max_p-min_p+1+max_n-min_n+1
!       allocate(levtz(idtz))
       nn=0
       do i=min_p,max_p !1,id
        nn=nn+1
        levtz(nn)%index=nn
        levtz(nn)%ipar=lhfp(i)%ipar
        levtz(nn)%N=lhfp(i)%N
        levtz(nn)%nn=lhfp(i)%nn
        levtz(nn)%l=lhfp(i)%l
        levtz(nn)%j2=lhfp(i)%j2
        levtz(nn)%tz=-1  
        levtz(nn)%ph=0
        levtz(nn)%point=i
        levtz(nn)%spenrg=lhfp(i)%spenrg
        levtz(nn)%ei=lhfp(i)%ei
        levtz(nn)%qei=lhfp(i)%qei
        levtz(nn)%ui=lhfp(i)%ui
        levtz(nn)%vi=lhfp(i)%vi
        if(dabs(levtz(nn)%ui).gt.0.5d0) levtz(nn)%ph=1
       enddo
       do i=min_n,max_n !1,id
        nn=nn+1
        levtz(nn)%index=nn
        levtz(nn)%ipar=lhfn(i)%ipar
        levtz(nn)%N=lhfn(i)%N
        levtz(nn)%nn=lhfn(i)%nn
        levtz(nn)%l=lhfn(i)%l
        levtz(nn)%j2=lhfn(i)%j2
        levtz(nn)%tz=1  
        levtz(nn)%ph=0
        levtz(nn)%point=i
        levtz(nn)%spenrg=lhfn(i)%spenrg
        levtz(nn)%ei=lhfn(i)%ei
        levtz(nn)%qei=lhfn(i)%qei
        levtz(nn)%ui=lhfn(i)%ui
        levtz(nn)%vi=lhfn(i)%vi
        if(dabs(levtz(nn)%ui).gt.0.5d0) levtz(nn)%ph=1
       enddo

       open(1,file='levtzpp.out',status='unknown',form='formatted')
        write(1,*)'i,     l,    2*j,      ei,      qei,     vi
     &   tz     p/h   pointer'
        do ii = 1, idtz
         write(1,*) levtz(ii)%index,levtz(ii)%l,levtz(ii)%j2,
     &    levtz(ii)%ei,levtz(ii)%qei,levtz(ii)%vi,levtz(ii)%tz,
     &    levtz(ii)%ph,levtz(ii)%point
        end do
       close(1)

       open(1,file='LevIsopp.out',status='unknown',form='unformatted')
        write(1) idtz
        write(1) levtz
       close(1)

       jmax=0
       do i=min_p,max_p !1,id
        if(lhfp(i)%j2.gt.jmax) jmax=lhfp(i)%j2
       enddo
       do i=min_n,max_n !1,id
        if(lhfn(i)%j2.gt.jmax) jmax=lhfn(i)%j2
       enddo
!       do i=min_Y,max_Y !1,id
!        if(lhfY(i)%j2.gt.jmax) jmax=lhfY(i)%j2
!       enddo

!       allocate(sixj1(jmax,jmax,0:jmax,jmax,jmax,0:jmax))
!       sixj1=0.d0

!       do j1=1,jmax,2
!        do j2=1,jmax,2
!         do j3=0,jmax
!          do j4=1,jmax,2
!           do j5=1,jmax,2
!            do j6=0,jmax
!             sixj1(j1,j2,j3,j4,j5,j6)=sixj_int(j1,j2,2*j3,j4,j5,2*j6)
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       enddo

       open(10,file='VNNtzpp.out',form='formatted',status='unknown')
       do i=1,idtz
        do j=1,idtz
         do k=1,idtz
          do l=1,idtz
           do Jp=0,jmax
            if(dabs(VNNtz(i,j,k,l,Jp)).gt.1.d-7)then
             write(10,*) i,j,k,l,Jp,VNNtz(i,j,k,l,Jp)
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       close(10)

       open(1,file='TDAtz_enerpp.out',status='unknown',form='formatted')
        write(1,*)
       open(2,file='TDAtz_dimenspp.out',status='unknown',
     &                                                 form='formatted')
        write(2,*)
       open(3,file='phontz_strucpp.out',status='unknown',
     &                                                 form='formatted')
        write(3,*)
       open(4,file='phz_storepp.out',status='unknown',
     &                                               form='unformatted')
       open(5,file='phz_store2pp.out',status='unknown',
     &                                               form='unformatted')

       write(4) idtz,min_p,max_p,min_n,max_n,min_Y,max_Y

       do ipar=-1,1,2
        do Jp=0,jmax

        write(*,*) 'TDA calculation for pi=',ipar,'J=',Jp

       if(if_QTDA.eq.0) then
        i1=0
        do i=1,idtz
         if(levtz(i)%ph.eq.1) then !particle state
          do j=i,idtz
           if(levtz(j)%ph.eq.1) then !particle state
            if(levtz(i)%ipar*levtz(j)%ipar.eq.ipar) then
             if((levtz(i)%j2+levtz(j)%j2)/2.ge.Jp) then
              if(abs(levtz(i)%j2-levtz(j)%j2)/2.le.Jp) then
               i1=i1+1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo
        
        if(i1.gt.0) allocate(ph(i1),ityp1(i1))
        if(i1.gt.0) ityp1=0
!        if(i1.gt.0) allocate(list_phon(i1))

        i1=0
        do i=1,idtz
         if(levtz(i)%ph.eq.1) then !particle state
          do j=i,idtz
           if(levtz(j)%ph.eq.1) then !particle state
            if(levtz(i)%ipar*levtz(j)%ipar.eq.ipar) then
             if((levtz(i)%j2+levtz(j)%j2)/2.ge.Jp) then
              if(abs(levtz(i)%j2-levtz(j)%j2)/2.le.Jp) then
               i1=i1+1
               ph(i1)%q1=i  !hole
               ph(i1)%q2=j  !particle
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo      
       endif

       i1_spur=0
!       if(if_ort.eq.1.and.(ipar.eq.-1.and.Jp.eq.1)) i1_spur=1

       write(2,*) ' pi=',ipar,'J=',Jp
       write(2,*) 'dimension=',i1-i1_spur

       if(i1.gt.0) allocate(TDA_matrix(i1,i1))
       if(i1.gt.0) allocate(w(i1))
       if(i1.gt.0) TDA_matrix=0.d0

!       if(if_ort.eq.1) then
!       if(ipar.eq.-1.and.Jp.eq.1) allocate(spur_CM(i1))
!       if(ipar.eq.-1.and.Jp.eq.1) 
!     &  call spur_vectz(spur_CM,i1,Jp,ph)   
!       endif

       if(i1.gt.0) then

        if(if_QTDA.eq.0) then
         do i=1,i1
          p1=ph(i)%q2
          h1=ph(i)%q1
          do j=1,i1
           p2=ph(j)%q2
           h2=ph(j)%q1
            if(i.eq.j) TDA_matrix(i,j)=levtz(p1)%ei+levtz(h1)%ei
            TDA_matrix(i,j)=TDA_matrix(i,j)+VNNtz(p2,h2,p1,h1,Jp)
          enddo
         enddo
        endif

!        if(ipar.eq.-1.and.Jp.eq.1) then
!         write(5) i1
!         write(5) TDA_matrix
!        endif

!        if(.not.(if_ort.eq.1.and.(ipar.eq.-1.and.Jp.eq.1))) then

        call diagonalization(TDA_matrix,w,i1)

         write(1,*) ' pi=',ipar,'J=',Jp
         write(1,*) w(1:i1) 

        write(3,*) ' pi=',ipar,'J=',Jp
        write(3,*)
        do i=1,i1
         write(3,*) 'phonon  E_i=',w(i)
         write(3,*) 
         do j=1,i1
          m1=ph(j)%q1
          m2=ph(j)%q2
          if(TDA_matrix(j,i)**2.d0.gt.1.d-2) then
            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))') 
     &             'T',levtz(m1)%l,levtz(m1)%j2,levtz(m1)%ei,
     &                 levtz(m2)%l,levtz(m2)%j2,levtz(m2)%ei,
     &                 TDA_matrix(j,i)**2.d0
          endif
         enddo
         write(3,*) 
        enddo

!        elseif(ipar.eq.-1.and.Jp.eq.1) then
!         call diag_subspace(TDA_matrix,w,spur_CM,i1)
!         write(1,*) ' pi=',ipar,'J=',Jp
!         write(1,*) w(2:i1)
!         write(3,*) ' pi=',ipar,'J=',Jp
!         write(3,*)
!         do i=2,i1
!          write(3,*) 'phonon  E_i=',w(i)
!          write(3,*)
!          do j=1,i1
!           m1=ph(j)%q1
!           m2=ph(j)%q2
!           if(TDA_matrix(j,i)**2.d0.gt.1.d-2) then
!            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
!     &             'T',levtz(m1)%l,levtz(m1)%j2,levtz(m1)%ei,
!     &                 levtz(m2)%l,levtz(m2)%j2,levtz(m2)%ei,
!     &                 TDA_matrix(j,i)**2.d0
!           endif
!          enddo
!          write(3,*)
!         enddo
!        endif

       write(4) ipar,Jp,i1
       write(4) (ph(i),i=1,i1)
       write(4) (w(i),i=1,i1)

       do i=1,i1
        write(4) (TDA_matrix(j,i),j=1,i1)   !list_phon(i)
       enddo

       do i=1,i1
        valmax=0.d0
        do j=1,i1
         if(dabs(TDA_matrix(j,i)).gt.valmax) then 
           valmax=dabs(TDA_matrix(j,i))
           k=j
         endif
        enddo
        ityp1(i)=levtz(ph(k)%q2)%tz+levtz(ph(k)%q1)%tz
       enddo
       write(4) ityp1

       endif

        if(i1.gt.0) deallocate(ph,TDA_matrix,w,ityp1)  !,list_phon)

!        if(if_ort.eq.1) then
!        if(ipar.eq.-1.and.Jp.eq.1) deallocate(spur_CM)
!        endif

        enddo
       enddo

       close(1)
       close(2)
       close(3)
       close(4)
       close(5)
       write(*,*) 'TDA calculation has finished.'

       return
      end
