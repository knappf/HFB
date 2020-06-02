      MODULE technical

       include 'typedef.inc'

        PUBLIC :: kin_p, kin_n, kin_Y, levp, levn, lhfp, lhfn,
     &            levY,lhfY, UY_HFB, VY_HFB, AY_HFB, BY_HFB, 
     &            rhoY_HFB, kapY_HFB,
     &            Up_HFB, Vp_HFB, Un_HFB, Vn_HFB,
     &            Ap_HFB, Bp_HFB, An_HFB, Bn_HFB,
     &            rhop_HFB, rhon_HFB, kapp_HFB, kapn_HFB,
     &            Vpp,Vnn,Vpn,Fpp,Fnn,Fpn,tran_p,tran_n,tran_Y,
     &            VpY,VnY,V3B,VNNtz,FNNtz,VNLtz,
     &            trE0_p,trE0_n,trE1_p,trE1_n,trE2_p,trE2_n,
     &            trE3_p,trE3_n,trEN_p,trEN_n,trS1_p,trS1_n,
     &            trM1s_p,trM1s_n,trM1l_p,trM1l_n,
     &            trE0_Y,trE1_Y,trE2_Y,trE3_Y,trEN_Y,trS1_Y,
     &            trM1s_Y,trM1l_Y,
     &            trE1_p_dens,trE1_n_dens,
     &            rad_den1,H11p,H11n,zcross,weight,
     &            lnl, lev3, lpoint, 
     &            lev1pn,lp1,lp2,cg3,larrow1,larrow2,larrowm,iphase,
     &            lev1pnm,levtz

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kin_p(:,:), kin_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kin_Y(:,:)

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(level_type),dimension(:),allocatable, save :: lhfp,lhfn

        type(level_type),dimension(:),allocatable, save :: levY,lhfY

        type(level_type),dimension(:),allocatable, save :: lnl

        type(level_type),dimension(:),allocatable, save :: lev1pn
        type(levelm_type),dimension(:),allocatable, save :: lev1pnm

        type(level3_type),dimension(:),allocatable,save :: lev3

        type(leveltz_type),dimension(:),allocatable,save :: levtz

        INTEGER, ALLOCATABLE, SAVE :: lpoint(:,:,:,:,:,:,:)
        INTEGER, ALLOCATABLE, SAVE :: lp1(:),lp2(:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Up_HFB(:,:), Vp_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Un_HFB(:,:), Vn_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: UY_HFB(:,:), VY_HFB(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Ap_HFB(:,:), Bp_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: An_HFB(:,:), Bn_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: AY_HFB(:,:), BY_HFB(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rhop_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rhon_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rhoY_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kapp_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kapn_HFB(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: kapY_HFB(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: VpY(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: VnY(:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: VNNtz(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: FNNtz(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: VNLtz(:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: V3B(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpp(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fnn(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpn(:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: tran_p(:,:),tran_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: tran_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE0_p(:,:),trE0_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE0_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_p(:,:),trE1_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE2_p(:,:),trE2_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE2_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE3_p(:,:),trE3_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE3_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trEN_p(:,:),trEN_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trEN_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trS1_p(:,:),trS1_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trS1_Y(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trM1s_p(:,:),trM1s_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trM1s_Y(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trM1l_p(:,:),trM1l_n(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trM1l_Y(:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_p_dens(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: trE1_n_dens(:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: rad_den1(:)

        double precision, allocatable, save :: H11p(:,:),H11n(:,:)

        double precision, allocatable, save :: zcross(:),weight(:)

        double precision, allocatable, save :: cg3(:,:,:,:,:)

        type(elem_type) :: larrow1,larrow2,larrowm
        INTEGER :: iphase

      END MODULE technical
