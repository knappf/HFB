
       parameter (pi=3.141592653589793d0)
       parameter (hbarc=197.3269680000d0)        ! hbar*c [MeV*fm]
       parameter (zmp=938.272310000000d0)        ! proton mass
       parameter (zmn=939.565630000000d0)        ! neutron mass
       parameter (zmY=1115.68300000000d0)        ! Lambda mass
       parameter (charge2=1.4398840000d0)        ! electric charge to the power of 2

       INTEGER :: AZ, AN        ! proton and neutron number
       INTEGER :: noscmax       ! the maximal oscillator shell (defines the configuration space: #shells = noscmax + 1) for HFB
       INTEGER :: noscmax12     ! the maximum of N1 + N2 <= noscmax12
       INTEGER :: noscmax123    ! the maximum of N1 + N2 + N3 <= noscmax123
       INTEGER :: nosc_TDA       ! the maximal oscillator shell (defines the configuration space: #shells = noscmax + 1) for TDA
       INTEGER :: if_QTDA       ! = 0 for particle-hole TDA, = 1 for QTDA
       INTEGER :: if_ort        ! = 0 without TDA orthogonalization, = 1 with TDA orthogonalization
       INTEGER :: if_self       ! = 0 Lambda non treated self-consistently, = 1 Lambda treated self-consistently
       INTEGER :: if_2b         ! = 0 noscmax12 is not applied on V^{2B}, = 1 noscmax12 is applied also on V^{2B}
       INTEGER :: igrid,igrid2
       INTEGER :: min_p,max_p,min_n,max_n,min_Y,max_Y
       INTEGER :: n_lam_occup

       DOUBLE PRECISION :: precis  ! parameter to stop the iteration step
       DOUBLE PRECISION :: hbarom  ! hbar*omega [MeV]
       DOUBLE PRECISION :: quenp,quenn,quenpn,F_0
       DOUBLE PRECISION :: V3b_NNN          ! strength of the NNN interaction
       DOUBLE PRECISION :: V3b_LNN          ! strength of the LambdaNN interaction
       DOUBLE PRECISION :: sizebox

       logical :: ifp_hfb, ifn_hfb

       common /input/ AZ, AN, noscmax, noscmax12, noscmax123, nosc_TDA, 
     &                                  if_QTDA, if_ort, if_self, if_2b
       common /spacetda/ min_p,max_p,min_n,max_n,min_Y,max_Y
       common /n_lam_occupation/ n_lam_occup
       common /hb/ hbarom
       common /pre/ precis
       common /logi/ ifp_hfb, ifn_hfb
       common /qu/ quenp,quenn,quenpn,F_0
       common /DD/ V3b_NNN,V3b_LNN,sizebox,igrid,igrid2

!       b = sqrt(m*omega/hbar) = sqrt(m*c^2*hbar*omega)/(hbar*c) = sqrt(938.5*41*A^(-1/3) MeV^2)/(197.327 MeV*fm)
