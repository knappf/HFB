      subroutine kinetic

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       allocate(kin_p(id,id),kin_n(id,id),kin_Y(id,id))
       kin_p=0.d0
       kin_n=0.d0
       kin_Y=0.d0

       do j=1,id
        kin_p(j,j) = 0.5d0*hbarom*(dble(levp(j)%N)+1.5d0)
        kin_n(j,j) = 0.5d0*hbarom*(dble(levn(j)%N)+1.5d0)
        kin_Y(j,j) = 0.5d0*hbarom*(dble(levY(j)%N)+1.5d0)
        do k=1,id
         if(levp(j)%l.eq.levp(k)%l) then
          if(levp(j)%j2.eq.levp(k)%j2) then
           if(levp(j)%nn.eq.levp(k)%nn+1) then
            kin_p(j,k) = 0.5d0*hbarom*dsqrt(dble(levp(j)%nn)*     
     &     (dble(levp(j)%nn+levp(j)%l)+0.5d0))
           endif
           if(levp(j)%nn+1.eq.levp(k)%nn) then
            kin_p(j,k) = 0.5d0*hbarom*dsqrt(dble(levp(k)%nn)*
     &     (dble(levp(k)%nn+levp(k)%l)+0.5d0))
           endif
          endif
         endif
         if(levn(j)%l.eq.levn(k)%l) then
          if(levn(j)%j2.eq.levn(k)%j2) then
           if(levn(j)%nn.eq.levn(k)%nn+1) then
            kin_n(j,k) = 0.5d0*hbarom*dsqrt(dble(levn(j)%nn)*
     &     (dble(levn(j)%nn+levn(j)%l)+0.5d0))
           endif
           if(levn(j)%nn+1.eq.levn(k)%nn) then
            kin_n(j,k) = 0.5d0*hbarom*dsqrt(dble(levn(k)%nn)*
     &     (dble(levn(k)%nn+levn(k)%l)+0.5d0))
           endif
          endif
         endif
         if(levY(j)%l.eq.levY(k)%l) then
          if(levY(j)%j2.eq.levY(k)%j2) then
           if(levY(j)%nn.eq.levY(k)%nn+1) then
            kin_Y(j,k) = 0.5d0*hbarom*dsqrt(dble(levY(j)%nn)*
     &     (dble(levY(j)%nn+levY(j)%l)+0.5d0))
           endif
           if(levY(j)%nn+1.eq.levY(k)%nn) then
            kin_Y(j,k) = 0.5d0*hbarom*dsqrt(dble(levY(k)%nn)*
     &     (dble(levY(k)%nn+levY(k)%l)+0.5d0))
           endif
          endif
         endif
        enddo
       enddo

       if(if_self.eq.1) ACM = dble(AZ+AN)+1.19d0
       if(if_self.eq.0) ACM = dble(AZ+AN)

       do i=1,id
        do j=1,id
         kin_p(i,j)=kin_p(i,j)*(dble(1)-dble(1)/ACM)
         kin_n(i,j)=kin_n(i,j)*(dble(1)-dble(1)/ACM)
         kin_Y(i,j)=kin_Y(i,j)*(dble(1)-dble(1)/(ACM*0.84d0))
        enddo
       enddo

       return
      end
