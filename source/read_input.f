      subroutine read_input

       include 'define.inc'
       include 'parameters.inc' 

       open(unit=1, file='input.dat', status= 'old', form= 'formatted')
        read(1,*) i,j
        read(1,*) k,k12,k123
        read(1,*) if_2b
        read(1,*) k1,k2,k3,k4,k5,k6
        read(1,*) if_q
        read(1,*) d
        read(1,*) e
        read(1,*) qp,qn,qpn,F0
        read(1,*) v1
        read(1,*) v1LN
        read(1,*) if_o
        read(1,*) if_s
        read(1,*) n_lam_occup
       close(1)

       AZ = j
       AN = i-j
       noscmax = k
       noscmax12 = k12
       noscmax123 = k123
       min_p = k1
       max_p = k2
       min_n = k3
       max_n = k4
       min_Y = k5
       max_Y = k6
       if_QTDA = if_q
       precis = d
       hbarom = e
       quenp = qp
       quenn = qn
       quenpn = qpn
       V3b_NNN = v1
       V3b_LNN = v1LN
       F_0=F0
       if_ort = if_o
       if_self = if_s

       write(*,*) 'HFB reads input'
       write(*,*) 'A=',AZ+AN,'Z=',AZ
       write(*,*) 'noscmax=',noscmax
!       write(*,*) 'nosc_TDA=',nosc_TDA
       write(*,*) 'precis=',precis
       write(*,*) 'hbarom=',hbarom,' MeV'
       write(*,*) 'quenching factors =',quenp,quenn,quenpn,F_0
       write(*,*) 'calculation in j-scheme'

       if(mod(AZ,2).eq.1) then
        write(*,*) 'Error! Wrong input AZ - just even nuclei allowed.'
        stop
       endif

       if(mod(AN,2).eq.1) then
        write(*,*) 'Error! Wrong input AN - just even nuclei allowed.'
        stop
       endif

       if(AZ.eq.2.or.(AZ.eq.8.or.(AZ.eq.20.or.(AZ.eq.28.or.
     &                   (AZ.eq.50.or.(AZ.eq.82.or.AZ.eq.126)))))) then
        ifp_hfb = .false.
       else
        ifp_hfb = .true.
       endif

       if(AN.eq.2.or.(AN.eq.8.or.(AN.eq.20.or.(AN.eq.28.or.
     &                   (AN.eq.50.or.(AN.eq.82.or.AN.eq.126)))))) then
        ifn_hfb = .false.
       else
        ifn_hfb = .true.
       endif

       if(ifp_hfb.and.ifn_hfb) then
        write(*,*) 
        write(*,*) 'This option (open-shell nucleus in both protons and 
     &neutrons) '
        write(*,*) 'was not tested and may be it will not run properly.'
        write(*,*) 
!        stop
       endif

       if(if_QTDA.ne.0.and.if_QTDA.ne.1) then
        write(*,*) 'The if_QTDA must be =0 or =1!!!'
        stop
       endif

       if(if_self.ne.0.and.if_self.ne.1) then
        write(*,*) 'The if_self must be =0 or =1!!!'
        stop
       endif

       if(if_QTDA.eq.0.and.(ifp_hfb.or.ifn_hfb)) then
        write(*,*) 'The particle-hole TDA cannot be calculated if 
     &non-zero pairing is active!!!'
        stop
       endif

       if(.not.(if_2b.eq.0.or.if_2b.eq.1)) then
        write(*,*) 'if_2b must be = 0 or = 1 !! '
        stop
       endif

      return
      end
