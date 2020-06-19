      double precision function FNNtz(i,j,k,l,Jp)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       FNNtz=0.d0

       val=0.d0
        do Jpp=0,jmax
         phase=(-1)**((levtz(j)%j2+levtz(k)%j2)/2-Jp-Jpp)
         tri=dble(2*Jpp+1)*sixj1(levtz(i)%j2,levtz(j)%j2,Jp,
     &                          levtz(l)%j2,levtz(k)%j2,Jpp)
         val=val+phase*tri*VNNtz(i,k,j,l,Jpp)
        enddo
       FNNtz=val

       return
      end
