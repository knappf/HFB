      subroutine power2(ch,i,A)

       include 'define.inc'

       double precision :: A(i,i)
       double precision :: B(i,i),C(i,i),D(i,i)

       character*1 :: ch

       B=A

       if(ch.eq.'T') then
        call dgemm('N','T',i,i,i,1.d0,A,max(1,i),B,
     &                        max(1,i),0.d0,D,max(1,i))
       endif

       if(ch.eq.'N') then
        call dgemm('N','N',i,i,i,1.d0,A,max(1,i),B,
     &                        max(1,i),0.d0,D,max(1,i))
       endif

       A=D

       return
      end
