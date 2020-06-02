      subroutine pnTDAph

       USE technical
       USE math
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       type(twoquas_type), allocatable, save :: ph(:)
!       type(phonon_type), allocatable, save :: list_phon(:)
       double precision, allocatable, save :: TDA_matrix(:,:),w(:)
       double precision, allocatable, save :: sixj1(:,:,:,:,:,:)
       integer :: p1,h1,tz1,p2,h2,tz2

       if(if_QTDA.eq.0) write(*,*) 'TDA calculation starts'
       if(if_QTDA.eq.1) write(*,*) 'QTDA calculation starts'

       allocate(sixj1(jmax,jmax,0:jmax,jmax,jmax,0:jmax))
       sixj1=0.d0

       do j1=1,jmax,2
        do j2=1,jmax,2
         do j3=0,jmax
          do j4=1,jmax,2
           do j5=1,jmax,2
            do j6=0,jmax
             sixj1(j1,j2,j3,j4,j5,j6)=sixj_int(j1,j2,2*j3,j4,j5,2*j6)
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo

       open(1,file='pnTDAph_en.out',status='unknown',form='formatted')
        write(1,*)
       open(2,file='pnTDAph_dim.out',status='unknown',form='formatted')
        write(2,*)
       open(3,file='pnphonph_str.out',status='unknown',form='formatted')
        write(3,*)
       open(4,file='pnTDAph_store.out',status='unknown',
     &                                               form='unformatted')

!       call dimm(id,nosc_TDA)

       jmax=0
       do i=min_p,max_p !1,id
        if(lhfp(i)%j2.gt.jmax) jmax=lhfp(i)%j2
       enddo
       do i=min_n,max_n !1,id
        if(lhfn(i)%j2.gt.jmax) jmax=lhfn(i)%j2
       enddo
       do i=min_Y,max_Y !1,id
        if(lhfY(i)%j2.gt.jmax) jmax=lhfY(i)%j2
       enddo

       write(4) min_p,max_p,min_n,max_n,min_Y,max_Y

      do ipar=-1,1,2
        do Jp=0,jmax

        write(*,*) 'pnTDAph calculation for pi=',ipar,'J=',Jp

       if(if_QTDA.eq.0) then
        i1=0
        do i=min_p,max_p !1,id
         if(dabs(lhfp(i)%ui**2.d0-1.d0).lt.1.d-2) then
          do j=min_n,max_n !1,id
           if(dabs(lhfn(j)%vi**2.d0-1.d0).lt.1.d-2) then
            if(lhfp(i)%ipar*lhfn(j)%ipar.eq.ipar) then
             if((lhfp(i)%j2+lhfn(j)%j2)/2.ge.Jp) then
              if(abs(lhfp(i)%j2-lhfn(j)%j2)/2.le.Jp) then
               i1=i1+1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo

        if(i1.gt.0) allocate(ph(i1))
!        if(i1.gt.0) allocate(list_phon(i1))

        i1=0
        do i=min_p,max_p !1,id
         if(dabs(lhfp(i)%ui**2.d0-1.d0).lt.1.d-2) then
          do j=min_n,max_n !1,id
           if(dabs(lhfn(j)%vi**2.d0-1.d0).lt.1.d-2) then
            if(lhfp(i)%ipar*lhfn(j)%ipar.eq.ipar) then
             if((lhfp(i)%j2+lhfn(j)%j2)/2.ge.Jp) then
              if(abs(lhfp(i)%j2-lhfn(j)%j2)/2.le.Jp) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=-1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo

       endif

       i1_spur=0

       write(2,*) ' pi=',ipar,'J=',Jp
       write(2,*) 'dimension=',i1-i1_spur

       if(i1.gt.0) allocate(TDA_matrix(i1,i1))
       if(i1.gt.0) allocate(w(i1))
       if(i1.gt.0) TDA_matrix=0.d0

       if(i1.gt.0) then

        if(if_QTDA.eq.0) then
         do i=1,i1
          p1=ph(i)%q2  !neutron index
          h1=ph(i)%q1  !proton index
          tz1=ph(i)%tz
          do j=1,i1
           p2=ph(j)%q2  !neutron index
           h2=ph(j)%q1  !proton index
           tz2=ph(j)%tz
           if(tz1.eq.-1.and.tz2.eq.-1) then
            phase=dble(-1)
            if(i.eq.j) TDA_matrix(i,j)=lhfp(h1)%ei-lhfn(p1)%ei
             do Jpp=0,jmax
             TDA_matrix(i,j)=TDA_matrix(i,j)+Vpn(h2,p1,h1,p2,Jpp)
     &        *phase*dble(2*Jpp+1)*sixj1(lhfp(h2)%j2,lhfn(p1)%j2,Jpp,
     &                                    lhfp(h1)%j2,lhfn(p2)%j2,Jp)
             enddo
           endif
          enddo
         enddo
        endif

        if(if_QTDA.eq.1) then
         write(*,*) 'Stop: pnTDAph not coded in the qp version.'
         stop
        endif

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
          m3=ph(j)%tz
          if(TDA_matrix(j,i)**2.d0.gt.1.d-2) then
           if(m3.eq.-1) then
            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'pn',lhfp(m1)%l,lhfp(m1)%j2,lhfp(m1)%ei,lhfn(m2)%l,
     &                 lhfn(m2)%j2,lhfn(m2)%ei,TDA_matrix(j,i)**2.d0
           endif
          endif
         enddo
         write(3,*)
        enddo

       write(4) ipar,Jp,i1
       write(4) (ph(i),i=1,i1)
       write(4) (w(i),i=1,i1)

       do i=1,i1
        write(4) (TDA_matrix(j,i),j=1,i1)  !list_phon(i)
       enddo

       endif

        if(i1.gt.0) deallocate(ph,TDA_matrix,w)  !,list_phon)

        enddo
      enddo

       close(1)
       close(2)
       close(3)
       close(4)

       deallocate(sixj1)

       write(*,*) 'pnTDAph calculation has finished.'

       return
      end
