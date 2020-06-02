      subroutine diagonalization(A,e,i)

       include 'define.inc'

       double precision :: A(i,i),e(i)
       double precision :: B(i,i),C(i,i),Y(i,i)
       double precision :: work(max(1,6*i))
       logical :: ortn

       ortn=.false.

       icou1=max(1,6*i)

       call dsyev('V','U',i,A,max(1,i),e,work,icou1,info)

       if(info.ne.0) then
        write(*,*) 'The diagonalization failed!!!'
        stop
       endif

       B=A
       C=0.d0
       Y=0.d0

       call dgemm('T','N',i,i,i,1.d0,B,max(1,i),B,
     &                        max(1,i),0.d0,C,max(1,i))

       do i1=1,i                                              ! This is a check that 
        do i2=1,i                                             ! the eigen-vectors really
         if(i1.ne.i2) then                                    ! form the orthonormal
          if(dabs(C(i1,i2)).gt.precis) then                   ! set of vectors
           ortn=.true.                                        ! 
          endif                                               ! Usually in the case of 
         endif                                                ! degenerate spectrum 
         if(i1.eq.i2) then                                    ! the eigen-vectors are 
          if(dabs(C(i1,i2)-1.d0).gt.precis) then              ! not orthonormal on the 
           ortn=.true.                                        ! output of dsyev routine
          endif
         endif
        enddo
       enddo

       if(ortn) then                   ! Here the Gram-Schmidt orthonormalization follows
        do i1=1,i
         n2=i1
         do i2=1,i
          Y(i2,i1)=A(i2,i1)
         enddo
         do j1=1,i1-1
          ccv=0.d0
          do j2=1,i
           ccv=ccv+A(j2,i1)*Y(j2,j1)
          enddo
          do j2=1,i
           Y(j2,i1)=Y(j2,i1)-ccv*Y(j2,j1)
          enddo
         enddo
         anorm=0.d0
         do m=1,i
          anorm=anorm+Y(m,i1)**2.d0
         enddo
         if(dabs(anorm).lt.precis) then
          write(*,*) 'Error in the diagonalization: wrong norm!!!'
          stop
         endif
         do n=1,i
          Y(n,i1)=Y(n,i1)/dsqrt(anorm)
         enddo
        enddo
       endif

       return
      end
