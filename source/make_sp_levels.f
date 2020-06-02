      subroutine make_sp_levels(a1,a2,w1,w2,loop)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: w1(id),w2(id)
       double precision :: a1(id,id),a2(id,id)
       integer :: jlp(id),jln(id)

       do i=1,id
        d_max=0.d0
        do j=1,id
         if(dabs(a1(j,i)).gt.d_max) then
          d_max=dabs(a1(j,i))
          j_pointer=j
         endif
        enddo
        jlp(i)=1000*levp(j_pointer)%l+levp(j_pointer)%j2
       enddo
       do i=1,id
        llll=jlp(i)/1000
        jjjj=mod(jlp(i),1000)
        lhfp(i)%ipar=(-1)**llll
        lhfp(i)%l=llll
        lhfp(i)%j2=jjjj
        lhfp(i)%ei=w1(i)
        lhfp(i)%qei=dabs(w1(i)-ferp)
       enddo
!       if(loop.eq.1) then
!       i_min=0
!       i_max=0
!        do i=1,id
!         i_max=i_min+lhfp(i)%j2+1
!         if(AZ.ge.i_max) lhfp(i)%ui=0.d0
!         if(AZ.ge.i_max) lhfp(i)%vi=1.d0
!         if(AZ.le.i_min) lhfp(i)%ui=1.d0
!         if(AZ.le.i_min) lhfp(i)%vi=0.d0
!         if(AZ.lt.i_max.and.AZ.gt.i_min) lhfp(i)%ui=
!     &              dsqrt(1.d0-dble(AZ-i_min)/dble(i_max-i_min))
!         if(AZ.lt.i_max.and.AZ.gt.i_min) lhfp(i)%vi=
!     &              dsqrt(dble(AZ-i_min)/dble(i_max-i_min))
!         i_min=i_min+lhfp(i)%j2+1
!        enddo
!       endif
       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=0
         do i=1,id
          if(levp(i)%j2.eq.jj.and.levp(i)%l.eq.ll) then
           nn=nn+1
           mm=0
           do j=1,id
            if(lhfp(j)%j2.eq.jj.and.lhfp(j)%l.eq.ll) then
             mm=mm+1
             if(mm.eq.nn) lhfp(j)%ui=levp(i)%ui
             if(mm.eq.nn) lhfp(j)%vi=levp(i)%vi
            endif
           enddo
          endif
         enddo
        enddo
       enddo
       Vp_HFB=0.d0
       Up_HFB=0.d0
       do i=1,id
        Up_HFB(i,i)=lhfp(i)%ui
        Vp_HFB(i,i)=lhfp(i)%vi!*(-1)**((lhfp(i)%j2-1)/2)
       enddo
       call dgemm('N','N',id,id,id,1.d0,a1,max(1,id),Vp_HFB,
     &                        max(1,id),0.d0,Bp_HFB,max(1,id))
       Vp_HFB=Bp_HFB
       call dgemm('N','N',id,id,id,1.d0,a1,max(1,id),Up_HFB,
     &                        max(1,id),0.d0,Ap_HFB,max(1,id))
       Up_HFB=Ap_HFB
       if(.not.ifp_hfb) then
        e_filled=0.d0
        e_empty=0.d0
        do i=1,id
         if(dabs(lhfp(i)%vi-1.d0).lt.precis) e_filled=lhfp(i)%ei
         if(dabs(lhfp(id-i+1)%ui-1.d0).lt.precis)e_empty=lhfp(id-i+1)%ei
        enddo
        ferp=0.5d0*(e_empty+e_filled)
       endif
       if(ifp_hfb.and.loop.le.2) then
        do i=1,id
         if(dabs(lhfp(i)%vi).gt.precis) ferp=lhfp(i)%ei
        enddo
       endif

       do i=1,id
        d_max=0.d0
        do j=1,id
         if(dabs(a2(j,i)).gt.d_max) then
          d_max=dabs(a2(j,i))
          j_pointer=j
         endif
        enddo
        jln(i)=1000*levn(j_pointer)%l+levn(j_pointer)%j2
       enddo
       do i=1,id
        llll=jln(i)/1000
        jjjj=mod(jln(i),1000)
        lhfn(i)%ipar=(-1)**llll
        lhfn(i)%l=llll
        lhfn(i)%j2=jjjj
        lhfn(i)%ei=w2(i)
        lhfn(i)%qei=dabs(w2(i)-fern)
       enddo
!       if(loop.eq.1) then
!        i_min=0
!        i_max=0
!        do i=1,id
!         i_max=i_min+lhfn(i)%j2+1
!         if(AN.ge.i_max) lhfn(i)%ui=0.d0
!         if(AN.ge.i_max) lhfn(i)%vi=1.d0
!         if(AN.le.i_min) lhfn(i)%ui=1.d0
!         if(AN.le.i_min) lhfn(i)%vi=0.d0
!         if(AN.lt.i_max.and.AN.gt.i_min) lhfn(i)%ui=
!     &               dsqrt(1.d0-dble(AN-i_min)/dble(i_max-i_min))
!         if(AN.lt.i_max.and.AN.gt.i_min) lhfn(i)%vi=
!     &               dsqrt(dble(AN-i_min)/dble(i_max-i_min))
!         i_min=i_min+lhfn(i)%j2+1
!        enddo
!       endif
       do ll=0,(jmax-1)/2
        do jj=1,jmax,2
         nn=0
         do i=1,id
          if(levn(i)%j2.eq.jj.and.levn(i)%l.eq.ll) then
           nn=nn+1
           mm=0
           do j=1,id
            if(lhfn(j)%j2.eq.jj.and.lhfn(j)%l.eq.ll) then
             mm=mm+1
             if(mm.eq.nn) lhfn(j)%ui=levn(i)%ui
             if(mm.eq.nn) lhfn(j)%vi=levn(i)%vi
            endif
           enddo
          endif
         enddo
        enddo
       enddo
       Vn_HFB=0.d0
       Un_HFB=0.d0
       do i=1,id
        Un_HFB(i,i)=lhfn(i)%ui
        Vn_HFB(i,i)=lhfn(i)%vi!*(-1)**((lhfn(i)%j2-1)/2)
       enddo
       call dgemm('N','N',id,id,id,1.d0,a2,max(1,id),Vn_HFB,
     &                           max(1,id),0.d0,Bn_HFB,max(1,id))
       Vn_HFB=Bn_HFB
       call dgemm('N','N',id,id,id,1.d0,a2,max(1,id),Un_HFB,
     &                           max(1,id),0.d0,An_HFB,max(1,id))
       Un_HFB=An_HFB
       if(.not.ifn_hfb) then
        e_filled=0.d0
        e_empty=0.d0
        do i=1,id
         if(dabs(lhfn(i)%vi-1.d0).lt.precis) e_filled=lhfn(i)%ei
         if(dabs(lhfn(id-i+1)%ui-1.d0).lt.precis)e_empty=lhfn(id-i+1)%ei
        enddo
        fern=0.5d0*(e_empty+e_filled)
       endif
       if(ifn_hfb.and.loop.le.2) then
        do i=1,id
         if(dabs(lhfn(i)%vi).gt.precis) fern=lhfn(i)%ei
        enddo
       endif

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

       return
      end
