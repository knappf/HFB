      subroutine dimm(ii,n)           ! = dimension in m-scheme

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       i = 0

!******************************************************
!      This calculates the total dimension in m-scheme*
!       do j=0,n
!        i = i + (j+1)*(j+2)
!       end do
!        ii = i
!******************************************************

!******************************************************
!      This calculates the total dimension in j-scheme*
       do j=0,n
        i = i + (j+1)
       end do
        ii = i
!******************************************************

      return
      end subroutine
