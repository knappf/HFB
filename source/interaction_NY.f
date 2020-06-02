      subroutine interaction_NY

       USE technical
       USE geometric

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       open(unit=1,file='VNL.out',status='old',form='formatted')
        do while(.not.eof(1))
         read(1,*) i1,i2,i3,i4,ij,val1

         if((i1.le.id.and.i2.le.id).and.(i3.le.id.and.i4.le.id)) then
          if((if_2b.eq.1.and.(levY(i1)%N+levY(i2)%N.le.noscmax12.and.
     &                  levY(i3)%N+levY(i4)%N.le.noscmax12)).or.
     &      (if_2b.eq.0.and.(levY(i1)%N+levY(i2)%N.le.2*noscmax.and.
     &                  levY(i3)%N+levY(i4)%N.le.2*noscmax))) then
          VpY(i1,i2,i3,i4,ij)=VpY(i1,i2,i3,i4,ij)+val1
          VnY(i1,i2,i3,i4,ij)=VnY(i1,i2,i3,i4,ij)+val1
          endif
         endif
        end do
       close(1)

!       open(1,file='Vpy.out',status='unknown',form='formatted')
!       open(2,file='Vny.out',status='unknown',form='formatted')

!       do i1=1,id
!        do i2=1,id
!         do i3=1,id
!          do i4=1,id
!           do ij=0,2*noscmax+1

!            if(dabs(VpY(i1,i2,i3,i4,ij)).gt.1.d-7) 
!     &               write(1,*) i1,i2,i3,i4,ij,VpY(i1,i2,i3,i4,ij)
!            if(dabs(VnY(i1,i2,i3,i4,ij)).gt.1.d-7) 
!     &               write(2,*) i1,i2,i3,i4,ij,VnY(i1,i2,i3,i4,ij)

!           enddo
!          enddo
!         enddo
!        enddo
!       enddo

!       close(1)
!       close(2)

       return
      end
