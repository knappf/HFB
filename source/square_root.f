      subroutine square_root(i,A)

       include 'define.inc'

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: w(:),work(:)

       double precision :: A(i,i)
       double precision :: B(i,i),C(i,i),D(i,i)

       B=A
       C=0.d0
       D=A

       icou1=max(1,6*i) !originally max(1,3*i-1)
       allocate(w(i))
       allocate(work(icou1))

       call dsyev('V','U',i,B,max(1,i),w,work,icou1,info)    !diagonalization routine

       if(info.ne.0) then
        write(*,*) 'The diagonalization of matrix failed!'
        stop
       endif

       do k=1,i
        if(w(k).lt.0.d0) then
         write(*,*) 'Square root of matrix failed!'
         stop
        endif
       enddo

       do k=1,i
        do l=1,i
         C(k,l)=B(k,l)*dsqrt(w(l))
        enddo
       enddo

       call dgemm('N','T',i,i,i,1.d0,C,max(1,i),B,
     &                        max(1,i),0.d0,D,max(1,i))

       A=D

       deallocate(w,work)

       return
      end
