      subroutine deallocate_all
       
       USE technical

       deallocate(levp,levn,levY,lhfp,lhfn,lhfY,lnl,lev3,lpoint)
       deallocate(lev1pn,lev1pnm,levtz)
       deallocate(kin_p,kin_n,kin_Y)
       deallocate(Vpp,Vnn,Vpn)
       deallocate(VpY,VnY)
       deallocate(FNNtz,VNLtz)
       deallocate(V3B)
       deallocate(Up_HFB,Un_HFB,Vp_HFB,Vn_HFB)
       deallocate(Ap_HFB,An_HFB,Bp_HFB,Bn_HFB)
       deallocate(UY_HFB,VY_HFB,AY_HFB,BY_HFB)
       deallocate(rhop_HFB,rhon_HFB,kapp_HFB,kapn_HFB)
       deallocate(rhoY_HFB,kapY_HFB)
       deallocate(tran_p,tran_n,tran_Y)
       deallocate(trE0_p,trE0_n,trE1_p,trE1_n)
       deallocate(trE2_p,trE2_n,trE3_p,trE3_n)
       deallocate(trEN_p,trEN_n,trS1_p,trS1_n)
       deallocate(trM1s_p,trM1s_n,trM1l_p,trM1l_n)
       deallocate(trE0_Y,trE1_Y,trE2_Y,trE3_Y)
       deallocate(trEN_Y,trS1_Y,trM1s_Y,trM1l_Y)
       deallocate(trE1_p_dens,trE1_n_dens)
       deallocate(rad_den1)
       deallocate(zcross,weight)
       deallocate(H11p,H11n)
       deallocate(lp1,lp2,cg3)

       return
      end
