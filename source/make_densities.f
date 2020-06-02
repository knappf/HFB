      subroutine make_densities

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       call dgemm('N','T',id,id,id,1.d0,Vp_HFB,max(1,id),Vp_HFB,
     &                        max(1,id),0.d0,rhop_HFB,max(1,id))

       call dgemm('N','T',id,id,id,1.d0,Vp_HFB,max(1,id),Up_HFB,
     &                        max(1,id),0.d0,kapp_HFB,max(1,id))

       call dgemm('N','T',id,id,id,1.d0,Vn_HFB,max(1,id),Vn_HFB,
     &                        max(1,id),0.d0,rhon_HFB,max(1,id))

       call dgemm('N','T',id,id,id,1.d0,Vn_HFB,max(1,id),Un_HFB,
     &                        max(1,id),0.d0,kapn_HFB,max(1,id))

       return
      end
