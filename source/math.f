      module math

       include 'define.inc'

        PUBLIC :: delta1,delta2,delta3,delta4

        CONTAINS

        function  delta1(i1,j1)
         integer :: i1,j1
         if(i1.eq.j1) then
          delta1=1.d0
         else
          delta1=0.d0
         endif
        end function

        function  delta2(i1,i2,j1,j2)
         integer :: i1,i2,j1,j2
         if(i1.eq.j1.and.i2.eq.j2) then
          delta2=1.d0
         elseif(i1.eq.j2.and.i2.eq.j1) then
          delta2=-1.d0
         else
          delta2=0.d0
         endif
        end function

        function  delta3(i1,i2,i3,j1,j2,j3)
         integer :: i1,i2,i3,j1,j2,j3
         delta3=0.d0
        end function

        function  delta4(i1,i2,i3,i4,j1,j2,j3,j4)
         integer :: i1,i2,i3,i4,j1,j2,j3,j4
         delta4=0.d0
        end function

      end module math
