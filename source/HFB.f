      ! test Git
      program Phonon_EDF

       USE technical
       USE math

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

        call read_input

        call dimm(id,noscmax)  !calculates the dimension 'id' for no. of shells 'noscmax'

        call initialize_basis
        call occupations

        call kinetic
        call interaction
 
        call transition 

        call hfb_iteration

        call transf_interaction

        if(if_self.eq.0) call TDA

        if(if_self.eq.0) call NYTDA

        if(if_self.eq.0) call pnTDApp

        if(if_self.eq.0) call pnTDAhh

        if(if_self.eq.0) call pnTDAph

        if(if_self.eq.0) call pnTDAhp

        if(if_self.eq.0) call ppTDApp

        if(if_self.eq.0) call nnTDApp

        if(if_self.eq.0) call TDAtz

        if(if_self.eq.0) call TDAtzpp

        call deallocate_all

       stop
      end

