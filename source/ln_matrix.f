      subroutine ln_matrix(i,A)

       include 'define.inc'

       double precision :: A(i,i)
       double precision :: B(i,i),C(i,i),D(i,i)

       B=A
       C=0.d0
       D=A

       kk=1
       value_min=1.d10

       do while(value_min.gt.precis)

        kk=kk+1

        call dgemm('N','N',i,i,i,-1.d0,A,max(1,i),B,max(1,i),
     &       0.d0,C,max(1,i))
        B=C
        D=D+B/dble(kk)

        value_min=0.d0
        do ii1=1,i
         do ii2=1,i
          if(dabs(C(ii1,ii2)/dble(kk)).gt.value_min) 
     &                               value_min=dabs(C(ii1,ii2)/dble(kk))
         enddo
        enddo

        if(kk.gt.1000) then
         write(*,*) 'Error!! ln(A) did not converge'
         stop
        endif

       enddo

!       do ii2=1,i
!        D(ii2,ii2)=D(ii2,ii2)+1.d0
!       enddo

       A=D

       return
      end
