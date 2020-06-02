      subroutine make_qsp_levels(ww1,ww2,h1,h2,w1,w2)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: w1(id),w2(id)
       double precision :: ww1(id),ww2(id)
!       double precision :: a1(id,id),a2(id,id)
       double precision :: h1(id,id),h2(id,id)
       double precision :: Z(id,id)
       integer :: jlp(id),jln(id)

       if(ifp_hfb) then
       i_it_pair=0
       vvsum=1.d9
       do while(dabs(vvsum-dble(AZ)).ge.precis.and.i_it_pair.le.100000)
        i_it_pair=i_it_pair+1
!        call dgemm('T','N',id,id,id,1.d0,a1,max(1,id),h1,
!     &                          max(1,id),0.d0,Z,max(1,id))
        do i=1,id
         d_max=0.d0
         do j=1,id
          if(dabs(h1(j,i)).gt.d_max) then
           d_max=dabs(h1(j,i))
           j_pointer=j
          endif
         enddo
         jlp(i)=1000*levp(j_pointer)%l+levp(j_pointer)%j2
        enddo
        do i=1,id
         llll=jlp(i)/1000
         jjjj=mod(jlp(i),1000)
!         d_max=0.d0
!         do k1=1,id
!          if(dabs(Z(k1,i)).gt.d_max) then
!           d_max=dabs(Z(k1,i))
!           nklm=k1
!          endif
!         enddo
         h_elem=ww1(i)
         denom=dsqrt(h_elem**2.d0+w1(i)**2.d0)
         uusum=0.5d0*(1.d0+h_elem/denom) !0.5d0*(1.d0+h_elem/dsqrt(w1(nklm)))
         vvsum=0.5d0*(1.d0-h_elem/denom) !0.5d0*(1.d0-h_elem/dsqrt(w1(nklm)))
         lhfp(i)%ipar=(-1)**llll
         lhfp(i)%l=llll
         lhfp(i)%j2=jjjj
         lhfp(i)%qei=denom !dsqrt(w1(nklm))
         lhfp(i)%ei=h_elem+ferp
         lhfp(i)%ui=dsqrt(uusum)
         lhfp(i)%vi=dsqrt(vvsum)
        enddo
        Vp_HFB=0.d0
        Up_HFB=0.d0
        do i=1,id
         Up_HFB(i,i)=lhfp(i)%ui
         Vp_HFB(i,i)=lhfp(i)%vi!*(-1)**((lhfp(i)%j2-1)/2)
        enddo
        call dgemm('N','N',id,id,id,1.d0,h1,max(1,id),Vp_HFB,
     &                              max(1,id),0.d0,Bp_HFB,max(1,id))
        Vp_HFB=Bp_HFB
        call dgemm('N','N',id,id,id,1.d0,h1,max(1,id),Up_HFB,
     &                              max(1,id),0.d0,Ap_HFB,max(1,id))
        Up_HFB=Ap_HFB
        vvsum=0.d0
        do j=1,id
         vvsum=vvsum+lhfp(j)%vi**2.d0*dble(lhfp(j)%j2+1)
        enddo
        dferp=-0.001d0*(vvsum-dble(AZ))
        ferp=ferp+dferp
         do j=1,id
          ww1(j)=ww1(j)-dferp
         enddo
        enddo
       endif

       if(ifn_hfb) then
       i_it_pair=0
       vvsum=1.d9
       do while(dabs(vvsum-dble(AN)).ge.precis.and.i_it_pair.le.100000)
        i_it_pair=i_it_pair+1
!        call dgemm('T','N',id,id,id,1.d0,a2,max(1,id),h2,
!     &                          max(1,id),0.d0,Z,max(1,id))
        do i=1,id
         d_max=0.d0
         do j=1,id
          if(dabs(h2(j,i)).gt.d_max) then
           d_max=dabs(h2(j,i))
           j_pointer=j
          endif
         enddo
         jln(i)=1000*levn(j_pointer)%l+levn(j_pointer)%j2
        enddo
        do i=1,id
         llll=jln(i)/1000
         jjjj=mod(jln(i),1000)
!         d_max=0.d0
!         do k1=1,id
!          if(dabs(Z(k1,i)).gt.d_max) then
!           d_max=dabs(Z(k1,i))
!           nklm=k1
!          endif
!         enddo
         h_elem=ww2(i)
         denom=dsqrt(h_elem**2.d0+w2(i)**2.d0)
         uusum=0.5d0*(1.d0+h_elem/denom) !0.5d0*(1.d0+h_elem/dsqrt(w2(nklm)))
         vvsum=0.5d0*(1.d0-h_elem/denom) !0.5d0*(1.d0-h_elem/dsqrt(w2(nklm)))
         lhfn(i)%ipar=(-1)**llll
         lhfn(i)%l=llll
         lhfn(i)%j2=jjjj
         lhfn(i)%qei=denom !dsqrt(w2(nklm))
         lhfn(i)%ei=h_elem+fern
         lhfn(i)%ui=dsqrt(uusum)
         lhfn(i)%vi=dsqrt(vvsum)
        enddo
        Vn_HFB=0.d0
        Un_HFB=0.d0
        do i=1,id
         Un_HFB(i,i)=lhfn(i)%ui
         Vn_HFB(i,i)=lhfn(i)%vi!*(-1)**((lhfn(i)%j2-1)/2)
        enddo
        call dgemm('N','N',id,id,id,1.d0,h2,max(1,id),Vn_HFB,
     &                               max(1,id),0.d0,Bn_HFB,max(1,id))
        Vn_HFB=Bn_HFB
        call dgemm('N','N',id,id,id,1.d0,h2,max(1,id),Un_HFB,
     &                               max(1,id),0.d0,An_HFB,max(1,id))
        Un_HFB=An_HFB
        vvsum=0.d0
        do j=1,id
         vvsum=vvsum+lhfn(j)%vi**2.d0*dble(lhfn(j)%j2+1)
        enddo
        dfern=-0.001d0*(vvsum-dble(AN))
        fern=fern+dfern
         do j=1,id
          ww2(j)=ww2(j)-dfern
         enddo
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
