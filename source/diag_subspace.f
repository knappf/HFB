      subroutine diag_subspace(TDA_m,w1,spu,iii)

       USE technical
       USE math
       USE geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: TDA_m(iii,iii),w1(iii),spu(iii)
       double precision :: vectors(iii,iii)

       double precision :: TDA_sub(iii-1,iii-1),w2(iii-1)
       double precision :: transform(iii,iii-1),pp(iii,iii-1)

       w1=0.d0
       vectors=0.d0
       do k=1,iii
        vectors(k,k)=1.d0
       enddo

!*******************************************************************
!      Here we project basis vectors to the subspace which is 
!      the supplement to the spurious vector...

       do k=1,iii
        angle=0.d0
        do l=1,iii
         angle=angle+spu(l)*vectors(l,k)
        enddo
        do l=1,iii
         vectors(l,k)=vectors(l,k)-angle*spu(l)
        enddo
       enddo
!*******************************************************************

!*******************************************************************
!      Here we make Gram-Schmidt orthogonalization of the 
!      projected basis vectors
      
       do k=2,iii
        do l=1,k-1
         angle=0.d0
         amp=0.d0
         do l2=1,iii
          w1(l2)=vectors(l2,k)
         enddo
         do l2=1,iii
          angle=angle+w1(l2)*vectors(l2,l)
          amp=amp+vectors(l2,l)**2.d0
         enddo
         if(dabs(amp).gt.1.d-6) then
          do l2=1,iii
           vectors(l2,k)=vectors(l2,k)-vectors(l2,l)*angle/amp
          enddo
         endif
        enddo
       enddo
!*******************************************************************

       do k=1,iii
        amp=0.d0
        do l=1,iii
         amp=amp+vectors(l,k)**2.d0
        enddo
        if(amp.gt.1.d-6) then
         do l=1,iii
          vectors(l,k)=vectors(l,k)/dsqrt(amp)
         enddo
        endif
       enddo

       do k=1,iii
        h_max=-1.d0
        do l=1,iii
         if(dabs(vectors(l,k)).gt.h_max) h_max=dabs(vectors(l,k))
        enddo
        if(h_max.lt.1.d-7) l_pointer=k
       enddo

       do k=1,l_pointer-1
        do l=1,iii
         transform(l,k)=vectors(l,k)
        enddo
       enddo
       do k=l_pointer+1,iii
        do l=1,iii
         transform(l,k-1)=vectors(l,k)
        enddo
       enddo

       pp=0.d0
       TDA_sub=0.d0
       do k=1,iii 
        do l=1,iii-1
         do m=1,iii
          pp(k,l)=pp(k,l)+TDA_m(k,m)*transform(m,l)
         enddo
        enddo
       enddo
       do k=1,iii-1
        do l=1,iii-1
         do m=1,iii
          TDA_sub(k,l)=TDA_sub(k,l)+transform(m,k)*pp(m,l)
         enddo
        enddo
       enddo

       call diagonalization(TDA_sub,w2,iii-1)

       pp=0.d0
       do k=1,iii
        do l=1,iii-1
         do m=1,iii-1
          pp(k,l)=pp(k,l)+transform(k,m)*TDA_sub(m,l)
         enddo
        enddo
       enddo

       w1(1)=0.d0
       do k=2,iii
        w1(k)=w2(k-1)
       enddo

       do l=1,iii
        TDA_m(l,1)=spu(l)
       enddo

       do k=2,iii
        do l=1,iii
         TDA_m(l,k)=pp(l,k-1)
        enddo
       enddo

       return
      end
