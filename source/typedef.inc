!   Used types

      TYPE level_type
c      SEQUENCE
        INTEGER           :: index        ! the number of level
        INTEGER           :: ipar         ! parity
	INTEGER           :: N            ! the quantum number N=2*n+l
	INTEGER           :: nn           ! the quantum number n
	INTEGER           :: l            ! the quantum number l
        INTEGER           :: j2           ! the quantum number 2*j
        DOUBLE PRECISION  :: spenrg       ! s.p.-energy e_i = hbar*omega*(N+3/2)
	DOUBLE PRECISION  :: ei           ! s.p. HF energy
        DOUBLE PRECISION  :: qei          ! quasi-s.p. HFB energy
        DOUBLE PRECISION  :: ui           ! Bogolyubov coefficient u_i 
        DOUBLE PRECISION  :: vi           ! Bogolyubov coefficient v_i 
      END TYPE level_type

      TYPE leveltz_type
c      SEQUENCE
        INTEGER           :: index        ! the number of level
        INTEGER           :: ipar         ! parity
        INTEGER           :: N            ! the quantum number N=2*n+l
        INTEGER           :: nn           ! the quantum number n
        INTEGER           :: l            ! the quantum number l
        INTEGER           :: j2           ! the quantum number 2*j
        INTEGER           :: tz           ! projection of isospin tz= -1/2 for proton +1/2 for neutrons
        INTEGER           :: ph           ! = 1 for particle, = 0 for hole
        INTEGER           :: point        ! pointer to lhfp or lhfn
        DOUBLE PRECISION  :: spenrg       ! s.p.-energy e_i = hbar*omega*(N+3/2)
        DOUBLE PRECISION  :: ei           ! s.p. HF energy
        DOUBLE PRECISION  :: qei          ! quasi-s.p. HFB energy
        DOUBLE PRECISION  :: ui           ! Bogolyubov coefficient u_i
        DOUBLE PRECISION  :: vi           ! Bogolyubov coefficient v_i
      END TYPE leveltz_type

      TYPE levelm_type
c      SEQUENCE
        INTEGER           :: index        ! the number of level
        INTEGER           :: ipar         ! parity
        INTEGER           :: N            ! the quantum number N=2*n+l
        INTEGER           :: nn           ! the quantum number n
        INTEGER           :: l            ! the quantum number l
        INTEGER           :: j2           ! the quantum number 2*j
        INTEGER           :: m2           ! the quantum number 2*m
        INTEGER           :: jsch         ! the pointer to j-scheme
        DOUBLE PRECISION  :: spenrg       ! s.p.-energy e_i = hbar*omega*(N+3/2)
        DOUBLE PRECISION  :: ei           ! s.p. HF energy
        DOUBLE PRECISION  :: qei          ! quasi-s.p. HFB energy
        DOUBLE PRECISION  :: ui           ! Bogolyubov coefficient u_i
        DOUBLE PRECISION  :: vi           ! Bogolyubov coefficient v_i
      END TYPE levelm_type

      TYPE level3_type
c      SEQUENCE
        INTEGER           :: index        ! the number of level
        INTEGER           :: i            ! the i-th state
        INTEGER           :: j            ! the j-th state
        INTEGER           :: k            ! the k-th state
        INTEGER           :: Jab          ! the angular momentum [ji x jj]Jab
        INTEGER           :: JJ           ! 2-times the total angular momentum 
        INTEGER           :: Tab          ! the isospin [ti x tj]Tab
        INTEGER           :: TT           ! 2-times the total isospin
      END TYPE level3_type

      TYPE phonon_type
        INTEGER           :: noscmax      ! show the space in which the TDA phonon was calculated
	INTEGER           :: ipar         ! parity of the phonon
	INTEGER           :: jproj        ! j-projection of the phonon
	INTEGER           :: dimp         ! the dimension of the fields cph_p
	INTEGER           :: dimn         ! the dimension of the fields cph_n
	DOUBLE PRECISION  :: ener         ! energy of the phonon
	DOUBLE PRECISION  :: cph_p(112,112)   ! the amplitudes c_ph for protons
	DOUBLE PRECISION  :: cph_n(112,112)   ! the amplitudes c_ph for neutrons
      END TYPE phonon_type

      TYPE phonontz_type
        INTEGER           :: noscmax      ! show the space in which the TDA phonon was calculated
        INTEGER           :: ipar         ! parity of the phonon
        INTEGER           :: jproj        ! j-projection of the phonon
        INTEGER           :: dim          ! the dimension of the fields cph
        INTEGER           :: tz           ! = 2, = 0, = -2
        DOUBLE PRECISION  :: ener         ! energy of the phonon
        DOUBLE PRECISION  :: cph(250,250)   ! the amplitudes c_ph 
      END TYPE phonontz_type

      TYPE twoquas_type
        INTEGER           :: q1           ! the first index 
        INTEGER           :: q2           ! the second index
        INTEGER           :: tz           ! -1 for protons, +1 for neutrons
      END TYPE twoquas_type

      TYPE elem_type
c      SEQUENCE
        INTEGER           :: a            ! the 1-th state
        INTEGER           :: b            ! the 2-th state
        INTEGER           :: c            ! the 3-th state
        INTEGER           :: d            ! the 4-th state
        INTEGER           :: e            ! the 5-th state
        INTEGER           :: f            ! the 6-th state
      END TYPE elem_type
