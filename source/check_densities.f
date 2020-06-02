      subroutine check_densities

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision, allocatable, save :: A1(:,:)

       precis2=dsqrt(precis)

       allocate(A1(id,id))
       A1=0.d0

       s1=0.d0
       s2=0.d0
       do i=1,id
        s1=s1+rhop_HFB(i,i)*dble(levp(i)%j2+1)
        s2=s2+rhon_HFB(i,i)*dble(levn(i)%j2+1)
       enddo

       if(dabs(s1-dble(AZ)).gt.precis2) then
        if(.not.ifp_hfb) then
         write(*,*) 'Error: wrong proton number AZ!!!'
         write(*,*) s1,'=/=',AZ
         stop
        endif
       endif
       if(dabs(s2-dble(AN)).gt.precis2) then
        if(.not.ifn_hfb) then
         write(*,*) 'Error: wrong neutron number AN!!!'
         write(*,*) s2,'=/=',AN
         stop
        endif
       endif

       call dgemm('N','T',id,id,id,1.d0,Vp_HFB,max(1,id),Vp_HFB,
     &                        max(1,id),0.d0,A1,max(1,id))
       call dgemm('N','T',id,id,id,1.d0,Up_HFB,max(1,id),Up_HFB,
     &                        max(1,id),1.d0,A1,max(1,id))

       do i=1,id
        do j=1,id
         if(i.eq.j) then
          if(dabs(A1(i,j)-1.d0).gt.precis2) then
           write(*,*) 'Error: Not unitarity on U,V matrices!!!'
           stop
          endif
         else
          if(dabs(A1(i,j)).gt.precis2) then
           write(*,*) 'Error: Not unitarity on U,V matrices!!!'
           stop
          endif
         endif
        enddo
       enddo

       A1=0.d0
       call dgemm('N','T',id,id,id,1.d0,Vn_HFB,max(1,id),Vn_HFB,
     &                        max(1,id),0.d0,A1,max(1,id))
       call dgemm('N','T',id,id,id,1.d0,Un_HFB,max(1,id),Un_HFB,
     &                        max(1,id),1.d0,A1,max(1,id))

       do i=1,id
        do j=1,id
         if(i.eq.j) then
          if(dabs(A1(i,j)-1.d0).gt.precis2) then
           write(*,*) 'Error: Not unitarity on U,V matrices!!!'
           stop
          endif
         else
          if(dabs(A1(i,j)).gt.precis2) then
           write(*,*) 'Error: Not unitarity on U,V matrices!!!'
           stop
          endif
         endif
        enddo
       enddo

       A1=0.d0
       call dgemm('N','T',id,id,id,1.d0,Up_HFB,max(1,id),Vp_HFB,
     &                        max(1,id),0.d0,A1,max(1,id))
       call dgemm('N','T',id,id,id,1.d0,Vp_HFB,max(1,id),Up_HFB,
     &                        max(1,id),1.d0,A1,max(1,id))

       do i=1,id
        do j=1,id
         if(dabs(A1(i,j)).gt.precis2) then
!          write(*,*) 'Error: Not unitarity on U,V matrices!!!'
!          stop
         endif
        enddo
       enddo

       A1=0.d0
       call dgemm('N','T',id,id,id,1.d0,Un_HFB,max(1,id),Vn_HFB,
     &                        max(1,id),0.d0,A1,max(1,id))
       call dgemm('N','T',id,id,id,1.d0,Vn_HFB,max(1,id),Un_HFB,
     &                        max(1,id),1.d0,A1,max(1,id))

       do i=1,id
        do j=1,id
         if(dabs(A1(i,j)).gt.precis2) then
!          write(*,*) 'Error: Not unitarity on U,V matrices!!!'
!          stop
         endif
        enddo
       enddo

       deallocate(A1)

       return
      end
