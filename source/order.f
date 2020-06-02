      subroutine order

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       larrowm%a=1
       larrowm%b=2
       larrowm%c=3
       larrowm%d=4
       larrowm%e=5
       larrowm%f=6

       if(larrow1%a.lt.larrow1%b) then
        ii=larrow1%a
        larrow1%a=larrow1%b
        larrow1%b=ii
        ii=larrowm%a
        larrowm%a=larrowm%b
        larrowm%b=ii
        iphase=iphase*(-1)
       endif 

       if(larrow1%b.lt.larrow1%c) then
        ii=larrow1%b
        larrow1%b=larrow1%c
        larrow1%c=ii
        ii=larrowm%b
        larrowm%b=larrowm%c
        larrowm%c=ii
        iphase=iphase*(-1)
       endif

       if(larrow1%a.lt.larrow1%b) then
        ii=larrow1%a
        larrow1%a=larrow1%b
        larrow1%b=ii
        ii=larrowm%a
        larrowm%a=larrowm%b
        larrowm%b=ii
        iphase=iphase*(-1)
       endif

       if(larrow1%b.lt.larrow1%c) then
        ii=larrow1%b
        larrow1%b=larrow1%c
        larrow1%c=ii
        ii=larrowm%b
        larrowm%b=larrowm%c
        larrowm%c=ii
        iphase=iphase*(-1)
       endif

       if(larrow1%d.lt.larrow1%e) then
        ii=larrow1%d
        larrow1%d=larrow1%e
        larrow1%e=ii
        ii=larrowm%d
        larrowm%d=larrowm%e
        larrowm%e=ii
        iphase=iphase*(-1)
       endif

       if(larrow1%e.lt.larrow1%f) then
        ii=larrow1%e
        larrow1%e=larrow1%f
        larrow1%f=ii
        ii=larrowm%e
        larrowm%e=larrowm%f
        larrowm%f=ii
        iphase=iphase*(-1)
       endif

       if(larrow1%d.lt.larrow1%e) then
        ii=larrow1%d
        larrow1%d=larrow1%e
        larrow1%e=ii
        ii=larrowm%d
        larrowm%d=larrowm%e
        larrowm%e=ii
        iphase=iphase*(-1)
       endif

       if(larrow1%e.lt.larrow1%f) then
        ii=larrow1%e
        larrow1%e=larrow1%f
        larrow1%f=ii
        ii=larrowm%e
        larrowm%e=larrowm%f
        larrowm%f=ii
        iphase=iphase*(-1)
       endif

       larrow2=larrow1

       if(larrow1%a*1000000+larrow1%b*1000+larrow1%c.lt.
     &      larrow1%d*1000000+larrow1%e*1000+larrow1%f) then
        larrow2%a=larrow1%d
        larrow2%b=larrow1%e
        larrow2%c=larrow1%f
        larrow2%d=larrow1%a
        larrow2%e=larrow1%b
        larrow2%f=larrow1%c

        ii=larrowm%a
        jj=larrowm%b
        kk=larrowm%c
        larrowm%a=larrowm%d
        larrowm%b=larrowm%e
        larrowm%c=larrowm%f
        larrowm%d=ii
        larrowm%e=jj
        larrowm%f=kk
       endif

       return
      end
