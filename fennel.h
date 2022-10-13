      SUBROUTINE biology (ng,tile)
!
!svn $Id: redoxh.h UCPH-IGN: 2017 Daryabor & Bjerrum$
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                               Katja Fennel   !
!****************************************** Alexander F. Shchepetkin ***
!  Biogeochemical model based insight into the coupling of Nitrogen    !
!  and Sulfur cycles Azhar et al. (2014). This routine is developed    !
!  on the Fennel et al. (2006) ecosystem model based on the sulfur     !
!  cycling component to compute the biological source and sinks. The   !
!  detailed equations of the sulfur cycle are given in Azhar et al.    ! 
!  (2014).                                                             !
!                                                                      !
!  By activation of the "BIO_REDOXH" option a biogeochemical model    ! 
!  coupled with the nitrogen and sulfur cycles will be implemented.    !
!                                                                      !
!  It is recommended to activate always the  "BIO_SEDIMENT" option     !
!  to ensure conservation of mass by converting the organic matter     !
!  that is sinking out of the bottom most grid cell into inorganic     !
!  nutrients (i.e.,  instantanaous remineralization  at the water-     !
!  sediment interface). Additionally, the "DENITRIFICATION" option     !
!  can be activated.  Hence, a fraction of the instantenous bottom     !
!  remineralization is  assumed to  occur  through  the  anearobic     !
!  (denitrification)  pathway  and  thus  lost  from the  pool  of     !
!  biologically availalbe fixed nitrogen. See Fennel et al. (2006)     !
!  for details.                                                        !
!                                                                      !
!  Additional  options can be  activated to  enable  simulation of     !
!  inorganic carbon and dissolved oxygen.  Accounting of inorganic     !
!  carbon is activated by the "CARBON" option,  and results in two     !
!  additional  biological  tracer  variables:  DIC and alkalinity.     !
!  See Fennel et al. (2008) for details.                               !
!                                                                      !
!  If the "pCO2_RZ" options is activated, in addition to "CARBON",     !
!  the carbonate system  routines by Zeebe and Wolf-Gladrow (2001)     !
!  are used,  while the  OCMIP  standard routines are the default.     !
!  There are two different ways of treating alkalinity.  It can be     !
!  treated diagnostically (default),  in this case alkalinity acts     !
!  like a passive tracer  that is  not affected  by changes in the     !
!  concentration of  nitrate or ammonium.  However,  if the option     !
!  "TALK_NONCONSERV" is used,  the alkalinity  will be affected by     !
!  sources and sinks in nitrate. See Fennel et al. (2008) for more     !
!  details.                                                            !
!                                                                      !
!  If the "OXYGEN" option is activated,  one additional biological     !
!  tracer variable for dissolved oxygen. "OXYGEN" can be activated     !
!  independently of the  "CARBON"  option. If "OCMIP_OXYGEN_SC" is     !
!  used, in addition to "OXYGEN",  the Schmidt number of oxygen in     !
!  seawater will be  computed  using the  formulation  proposed by     !
!  Keeling et al. (1998, Global Biogeochem. Cycles,  12, 141-163).     !
!  Otherwise, the Wanninkhof's (1992) formula will be used.            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Azhar, M. A., Canfield, D. E., Fennel, K., Thamdrup, B., &        !
!      Bjerrum, C. J. (2014). A model‐based insight into the coupling  !
!      of nitrogen and sulfur cycles in a coastal upwelling system.    !
!      Journal of Geophysical Research: Biogeosciences, 119(3),        !
!      264-285.                                                        !
!                                                                      !
!    Fennel, K., Wilkin, J., Levin, J., Moisan, J., O Reilly, J.,      !
!      Haidvogel, D., 2006: Nitrogen cycling in the Mid Atlantic       !
!      Bight and implications for the North Atlantic nitrogen          !
!      budget: Results from a three-dimensional model.  Global         !
!      Biogeochemical Cycles 20, GB3007, doi:10.1029/2005GB002456.     !
!                                                                      !
!    Fennel, K., Wilkin, J., Previdi, M., Najjar, R. 2008:             !
!      Denitrification effects on air-sea CO2 flux in the coastal      !
!      ocean: Simulations for the Northwest North Atlantic.            !
!      Geophys. Res. Letters 35, L24608, doi:10.1029/2008GL036147.     !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping

!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                   GRID(ng) % rmask_full,                         &
# endif
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
#if defined CARBON || (defined OXYGEN || defined H_SULF) 
# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
# else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
# endif
#endif
#ifdef CARBON
     &                   OCEAN(ng) % pH,                                &
#endif
#if defined USECOS_BURIAL
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
#endif
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif

      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                         rmask_full,                              &
# endif
#endif
     &                         Hz, z_r, z_w, srflx,                     &
#if defined CARBON || (defined OXYGEN || defined H_SULF)
# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
# else
     &                         sustr, svstr,                            &
# endif
#endif
#ifdef CARBON
     &                         pH,                                      &
#endif
#if defined USECOS_BURIAL 
                               bustr, bvstr,                            &
#endif
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                      &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE dateclock_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# if defined CARBON || (defined OXYGEN || defined H_SULF)
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif
# endif
# if defined USECOS_BURIAL 
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# if defined CARBON || (defined OXYGEN || defined H_SULF)
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# if defined USECOS_BURIAL
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
# endif      
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
#ifdef CARBON
      integer, parameter :: Nsink = 9
#else
      integer, parameter :: Nsink = 7
#endif

      integer :: Iter, i, ibio, isink, itrc, ivar, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

#if defined CARBON || (defined OXYGEN || defined H_SULF)
      real(r8) :: u10squ
#endif
!
#if defined OXYGEN || defined H_SULF 
      real(r8), parameter :: OA0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: OA1 = 3.22014_r8       ! saturation
      real(r8), parameter :: OA2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: OA3 = 4.94457_r8
      real(r8), parameter :: OA4 =-0.256847_r8
      real(r8), parameter :: OA5 = 3.88767_r8
      real(r8), parameter :: OB0 =-0.00624523_r8
      real(r8), parameter :: OB1 =-0.00737614_r8
      real(r8), parameter :: OB2 =-0.0103410_r8
      real(r8), parameter :: OB3 =-0.00817083_r8
      real(r8), parameter :: OC0 =-0.000000488682_r8
!      real(r8), parameter :: rOxNO3= 8.625_r8       ! 138/16
!      real(r8), parameter :: rOxNH4= 6.625_r8       ! 106/16
      real(r8), parameter :: rOxNO3= 9.625_r8       ! 138/16
      real(r8), parameter :: rOxNH4= 7.625_r8       ! 106/16
      real(r8), parameter :: rOxPO4= 106.0_r8       ! 106/1
      real(r8), parameter :: rNxP_D= 45.0_r8        ! 45 
      real(r8), parameter :: rNxP_P= 16.0_r8        ! 16
      real(r8), parameter :: NxP_D= 0.02222222_r8   ! 45
      real(r8), parameter :: NxP_P= 0.0625_r8       ! 16
      !real(r8), parameter :: kinhO2dn = 0.1_r8      ! 7.0  (10-15)
      real(r8), parameter :: kinhO2dn = 1.0_r8      ! 7.0  (10-15)
      real(r8), parameter :: kinhO2an = 0.1_r8      ! 0.5  (1-5)
      real(r8), parameter :: kinhNO3an= 4.0_r8      ! 5.0  (1-5)
      real(r8), parameter :: kO2 = 0.3_r8           ! 5.0  (3)
      real(r8), parameter :: kNO3= 15.0_r8          !  7.0 (30)
!      real(r8), parameter :: kNO2= 30.0_r8          !  1.0  (?)
      real(r8), parameter :: kNO2= 1.0_r8          !  1.0  (?)	
      real(r8), parameter :: kSO2= 1.0_r8           ! Bo
      real(r8), parameter :: kSNO3= 3.0_r8          ! Bo
      real(r8), parameter :: kSNO2= 6.0_r8          ! Bo (?)
      real(r8), parameter :: KH2SO= 0.93_r8         ! Bo
      real(r8), parameter :: KH2SN1= 0.93_r8        ! Bo
      real(r8), parameter :: KH2SN2= 0.33_r8        ! Bo
      real(r8), parameter :: KANMX= 0.07_r8         ! Yakhusev 0.03
      real(r8), parameter :: LNO3= 1.0_r8           ! Anderson, 1982
      real(r8), parameter :: LNO2= 1.0_r8           ! 1-LNO3
      real(r8), parameter :: eNO2= 1.0_r8           ! electron NO2 to NO3
!     real(r8), parameter :: gNO2= 0.055_r8         ! rate nitrifi NO2 (d-1)
      real(r8), parameter :: I_thNO2= 0.0364_r8     ! Olson, 1981
      real(r8), parameter :: D_p5NO2= 0.074_r8      ! Olson, 1981 0.074
!      real(r8), parameter :: NitriR2= 0.005_r8      ! 0.005
      real(r8), parameter :: NitriR2= 0.05_r8      ! 0.05


      real(r8) :: l2mol = 1000.0_r8/22.3916_r8      ! liter to mol
#endif
#ifdef CARBON
      integer :: iday, month, year

      integer, parameter :: DoNewton = 0            ! pCO2 solver

      real(r8), parameter :: Acoef = 2073.1_r8      ! Schmidt
      real(r8), parameter :: Bcoef = 125.62_r8      ! number
      real(r8), parameter :: Ccoef = 3.6276_r8      ! transfer
      real(r8), parameter :: Dcoef = 0.043219_r8    ! coefficients

      real(r8), parameter :: A1 = -60.2409_r8       ! surface
      real(r8), parameter :: A2 = 93.4517_r8        ! CO2
      real(r8), parameter :: A3 = 23.3585_r8        ! solubility
      real(r8), parameter :: B1 = 0.023517_r8       ! coefficients
      real(r8), parameter :: B2 = -0.023656_r8
      real(r8), parameter :: B3 = 0.0047036_r8

      real(r8) :: pmonth                            ! months since Jan 1951
      real(r8) :: pCO2air_secular
      real(r8) :: yday, hour

      real(r8), parameter :: pi2 = 6.2831853071796_r8

      real(r8), parameter :: D0 = 282.6_r8          ! coefficients
      real(r8), parameter :: D1 = 0.125_r8          ! to calculate
      real(r8), parameter :: D2 =-7.18_r8           ! secular trend in
      real(r8), parameter :: D3 = 0.86_r8           ! atmospheric pCO2
      real(r8), parameter :: D4 =-0.99_r8
      real(r8), parameter :: D5 = 0.28_r8
      real(r8), parameter :: D6 =-0.80_r8
      real(r8), parameter :: D7 = 0.06_r8
#endif

      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR, K_PO41, K_PO42
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp, Epp2, Vp2
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4, t_PPmax2

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, cff7
      real(r8) :: cff22, cff222, cff55, cff41, cff71, cff51, ranmx
      real(r8) :: cff61, cff11, cffi, fac33, cff33, fanmx, fac22

      real(r8) :: fac1, fac2, fac3, fac12, fac4, cff8, cff9
      real(r8) :: cffL, cffR, cu, dltL, dltR, fac5, cff66, cff10

      real(r8) :: total_N, T_N, K_f, K_p, PhyIS2, L_PO4
      real(r8) :: ws, wscrit, Tcrit, amp, dlim2, fd2, cff32
      real(r8) :: oxlim, dlim1, srlim, tolim, fox, fd1, fsr
      real(r8) :: N_Flux_Lim, N_Flux_N, N_Flux_P, N_Flux_PmortalD
      real(r8) :: N_Flux_RemineSc, N_Flux_RemineLc, N_Flux_RemineDc
      real(r8) :: NO3_scratch, NH4_scratch, PO4_scratch
      real(r8) :: PNfrac1, PNfrac2, N_Flux_NewProd_NO3_lim
      real(r8) :: N_Flux_RegProd_lim, N_Flux_NewProd_P_lim

      real(r8) :: TSS

#ifdef DIAGNOSTICS_BIO
      real(r8) :: fiter
#endif

#if defined OXYGEN || defined H_SULF
      real(r8) :: SchmidtN_Ox, O2satu, O2_Flux
      real(r8) :: TS, AA
#endif

#ifdef CARBON
      real(r8) :: CO2_Flux, CO2_sol, SchmidtN, TempK, cff1C, cff2C
      real(r8) :: C_Flux_CoagP, C_Flux_CoagD, C_Flux_Remine
#endif

      real(r8) :: N_Flux_Assim, P_Flux_Assim
      real(r8) :: N_Assim_Excess
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: N_Flux_EgestP, N_Flux_EgestD
      real(r8) :: P_Flux_EgestP, P_Flux_EgestD
      real(r8) :: N_Flux_NewProd_NO3, N_Flux_RegProd
      real(r8) :: N_Flux_Nitrifi, N_Flux_NewProd_NFix
      real(r8) :: N_Flux_Nitrifi2, N_Flux_NewProd_P
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_Remine, P_Flux_Remine      
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo
      real(r8) :: cff3P, fac12P, fac1P, fac33P, fac34
      real(r8) :: fac34P, P_Flux_Zexcret
      real(r8) :: cff1P, cff2P, P_Flux_CoagP, P_Flux_CoagD
      real(r8) :: N_Flux_Zexcret_Diaz, N_Flux_Zexcret_Phyt
      real(r8) :: P_Flux_Zexcret_Diaz, P_Flux_Zexcret_Phyt

      real(r8) :: N_Flux_RemineD
      real(r8) :: P_Flux_RemineD
      real(r8) :: rDON
      real(r8) :: Cbe, CNbur, ReSuspR, Ustarb
      real(r8) :: cff12, cff13, cff14, cff15, cff16, cff17, cff18
      real(r8) :: N_Flux_SolubS, N_Flux_SolubL, P_Flux_SolubS, P_Flux_SolubL 
      real(r8) :: N_Flux_RemineSa, N_Flux_RemineLa, N_Flux_RemineDa
      real(r8) :: P_Flux_RemineSa, P_Flux_RemineLa, P_Flux_RemineDa
      real(r8) :: f_NTR, f_DNF,cff3c,cff6c,cff9c
      real(r8) :: cff3a,cff3b,cff13a,cff13b,cff16a,cff16b

      real(r8), dimension(Nsink) :: Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur
#ifdef CARBON
      real(r8), dimension(IminS:ImaxS) :: pCO2
#endif

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

#include "set_bounds.h"
#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                            &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
#ifdef DIAGNOSTICS_BIO
!
!  A factor to account for the number of iterations in accumulating
!  diagnostic rate variables.
!
      fiter=1.0_r8/REAL(BioIter(ng),r8)
#endif
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPhyt
      idsink(2)=iChlo
      idsink(3)=iSDeN
      idsink(4)=iLDeN
      idsink(5)=iSDeP
      idsink(6)=iLDeP
      idsink(7)=iDiaz
#ifdef CARBON
      idsink(8)=iSDeC
      idsink(9)=iLDeC
#endif
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! phytoplankton
      Wbio(2)=wPhy(ng)                ! chlorophyll
      Wbio(3)=wSDet(ng)               ! small Nitrogen-detritus
      Wbio(4)=wLDet(ng)               ! large Nitrogen-detritus
      Wbio(5)=wSDet(ng)               ! small P-detritus
      Wbio(6)=wLDet(ng)               ! large P-detritus
      Wbio(7)=wPhy(ng)                ! Diazotroph
#ifdef CARBON
      Wbio(8)=wSDet(ng)               ! small Carbon-detritus
      Wbio(9)=wLDet(ng)               ! large Carbon-detritus
#endif
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
        END DO
#ifdef CARBON
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio_old(i,k,iTIC_)=MIN(Bio_old(i,k,iTIC_),3000.0_r8)
            Bio_old(i,k,iTIC_)=MAX(Bio_old(i,k,iTIC_),400.0_r8)
            Bio(i,k,iTIC_)=Bio_old(i,k,iTIC_)
          END DO
        END DO
#endif
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            AttFac=0.0_r8
            IF (PARsur(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1

!ADD CHESROMS_BGC_ATT from original chesroms_ECB

                TSS = (Bio(i,k,iZoop)*ZooCN(ng) +                       &
     &                 Bio(i,k,iPhyt)*PhyCN(ng) +                       &
     &                 Bio(i,k,iSDeC)           +                       &
     &                 Bio(i,k,iLDeC))*12.0_r8/1000.0_r8
!               Convert TSS from g-C/m3 to g/m3:
                TSS = TSS * 2.9_r8 ! Cerco et al., 2017.
!# ifdef ISS_2_SIZE_CLASSES
!                TSS = TSS + Bio(i,k,iISS1) + Bio(i,k,iISS2)
!# endif
                Att=(rkd1(ng)+rkdChl1(ng)*Bio(i,k,iChlo)+               &
     &                rkdTSS1(ng)*TSS+rkdS1(ng)*Bio(i,k,isalt))
                Att=max( Att, 0.6_r8 ) ! The 0.6 is from Fei Da, 201705.
                Att=Att*(z_w(i,j,k)-z_w(i,j,k-1))              
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
 
!                Att=(AttSW(ng)+                                         &
!     &               AttChl(ng)*Bio(i,k,iChlo)+                         &
!     &               AttFac)*                                           &
!     &               (z_w(i,j,k)-z_w(i,j,k-1))


                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                cff=PhyCN(ng)*12.0_r8
                cffi=Bio(i,k,iPhyt)+Bio(i,k,iDiaz)
                Chl2C=MIN(Bio(i,k,iChlo)/(cffi*cff+eps),                &
     &                    Chl2C_m(ng))

!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.

            !    Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
#ifdef FIXED_Vp0

                Vp=Vp0(ng)

#else
!
!  Temperature-dependent growth rate
!  This replaces the Eppley-1972 parameterization of Fennel et al. 2006
!  Rate is inspired from Cerco&Noel 2004 and Lomas et al. 2002 (Q10).
!  psl20190904: Use linear Q10 of Lomas et al.2002; lower NPP in winter.

                Vp = max( 1.50_r8,                                      &
     &                    0.55_r8 * exp( 0.08065_r8 * Bio(i,k,itemp) ) )
!psl            Vp = max( 2.15_r8,                                      &
!psl &                    0.60_r8 * exp( 0.07800_r8 * Bio(i,k,itemp) ) )

#endif   
                fac1=PAR*PhyIS(ng)
                Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                t_PPmax=Epp*fac1

!
!  Nitrogen fixation (Fennel, 2002).
!
#ifdef BULK_FLUXES
                ws=(Uwind(i,j)**2+Vwind(i,j)**2)*1.22_r8*0.0013_r8
#else
                ws=SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+          &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
#endif
                wscrit=0.062_r8
                Tcrit=24.75_r8

             IF (ws.le.wscrit) then
              amp=(TANH(2.0_r8*(Bio(i,k,itemp)-Tcrit))+1.0_r8)/3.0_r8+  &
     &            (1.0_r8/3.0_r8)
             ELSE
              amp=(TANH(2.0_r8*(Bio(i,k,itemp)-Tcrit))+1.0_r8)/6.0_r8+  &
     &            (1.0_r8/6.0_r8)
             END IF


              !  Vp2=2.2_r8*0.085_r8*(1.066_r8**Bio(i,k,itemp))       ! FeFac == 1.0nm Fe open ocean
                Vp2=2.2_r8*0.0_r8*(1.066_r8**Bio(i,k,itemp))       ! change Vp2 to 0
                Epp2=Vp2/SQRT(Vp2*Vp2+fac1*fac1)
                t_PPmax2=Epp2*fac1                                   ! xamp

                fac1=dtdays*t_PPmax
                fac12=dtdays*t_PPmax2
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
                K_PO41=K_NO3(ng)*8.0_r8            !  16
                K_PO42=K_NO3(ng)*8.0_r8            !  160
            
                cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                cff22=Bio(i,k,iPO4_)*K_PO41        ! uptake P by plankt

                inhNH4=1.0_r8/(1.0_r8+cff1)
                L_NH4=cff1/(1.0_r8+cff1)
                L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                L_PO4=cff22/(1.0_r8+cff22) 
                LTOT=MIN(L_NO3+L_NH4,L_PO4)
!
!  Nitrate and ammonium uptake by Phytoplankton.
!  
!  Updated from Azhar et al. (2014) as follows: the sink terms for
!  NO3, NH4 and PO4 are each assigned initial scratch values; the
!  scratch value for phosphorus uptake (N_Flux_NewProd_P) assumes
!  that nitrogen is not limiting, and vice versa for the scratch
!  nitrogen uptake (N_Flux_N). Then, the smaller of these two values
!  is added to the phytoplankton biomass pool. Once this is done,
!  the actual sinks of PO4, NO3, and NH4 are computed by multiplying
!  the scratch values by the value PNfrac1 or PNfrac2. For example,
!  if nitrogen is limiting, then the scratch nitrogen sink flux is
!  made smaller by multiplying it by PNfrac1.  -KGH, Mar. 2018.


!  WC_NITRIFICATION && defined WC_DENITRIFICATION
                f_NTR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg) + kinhO2dn)
                f_DNF = kinhO2dn/(Bio(i,k,iOxyg) + kinhO2dn)

                cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*Bio(i,k,iPhyt)
                cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio(i,k,iPhyt)
                cff6=fac1*K_PO41/(1.0_r8+cff22)*(Bio(i,k,iPhyt))

                NO3_scratch=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                NH4_scratch=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                PO4_scratch=Bio(i,k,iPO4_)/(1.0_r8+cff6)

                N_Flux_NewProd_NO3=NO3_scratch*cff4
                N_Flux_RegProd=NH4_scratch*cff5
                N_Flux_NewProd_P=PO4_scratch*cff6*rNxP_P
                N_Flux_N=N_Flux_NewProd_NO3+N_Flux_RegProd
                N_Flux_Lim=MIN(N_Flux_N,N_Flux_NewProd_P)  
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                           &
     &                       (1.0_r8-EsDON(ng)-                          &          
     &              (f_NTR+f_DNF)*ElDON(ng))*N_Flux_Lim  

                PNfrac1=N_Flux_NewProd_P/N_Flux_N
                PNfrac2=1.0_r8/PNfrac1
                N_Flux_NewProd_P_lim=N_Flux_NewProd_P*MIN(PNfrac2,1.0_r8)
                N_Flux_NewProd_NO3_lim=N_Flux_NewProd_NO3*MIN(PNfrac1,1.0_r8)
                N_Flux_RegProd_lim=N_Flux_RegProd*MIN(PNfrac1,1.0_r8)
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-N_Flux_NewProd_P_lim*NxP_P
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd_NO3_lim
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd_lim
!
!  Nitrogen fixation by diazotrophs.
!
                cff222=Bio(i,k,iPO4_)*K_PO42
                cff55=fac12*K_PO42/(1.0_r8+cff222)*(Bio(i,k,iDiaz)/rNxP_D)
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff55)
                N_Flux_NewProd_NFix=Bio(i,k,iPO4_)*cff55*rNxP_D
                Bio(i,k,iDiaz)=Bio(i,k,iDiaz)+N_Flux_NewProd_NFix
                Bio(i,k,iN2__)=Bio(i,k,iN2__)-N_Flux_NewProd_NFix
!
!  End of Nitrogen fixation process
!
                Bio(i,k,iChlo)=Bio(i,k,iChlo)+                          &
     &                         (1.0_r8-EsDON(ng)-                       &
     &                             (f_NTR+f_DNF)*ElDON(ng))*            &
     &                         (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*       &
     &                          Chl2C_m(ng)*Bio(i,k,iChlo))/            &
     &                         (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)

#ifdef DIAGNOSTICS_BIO
                DiaBio3d(i,j,k,iPPro)=DiaBio3d(i,j,k,iPPro)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  &
# endif
     &                                (N_Flux_Lim+                      & 
     &                                N_Flux_NewProd_NFix)*fiter     
                DiaBio3d(i,j,k,iNO3u)=DiaBio3d(i,j,k,iNO3u)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  & 
# endif
     &                                N_Flux_NewProd_NO3_lim*fiter
                DiaBio3d(i,j,k,iPO4u)=DiaBio3d(i,j,k,iPO4u)+            &
# ifdef WET_DRY
     &                                rmask_full(i,j)*                  & 
# endif
     &                                N_Flux_NewProd_P_lim*NxP_P*fiter    
   	        DiaBio3d(i,j,k,iNFix)=DiaBio3d(i,j,k,iNFix)+            &
# ifdef WET_DRY
     &		                      rmask_full(i,j)*                  &
# endif
     &                                N_Flux_NewProd_NFix*              &
     &                                fiter
#endif

!  WC_NITRIFICATION && defined WC_DENITRIFICATION
                f_NTR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg) + kinhO2dn)
                f_DNF = kinhO2dn/(Bio(i,k,iOxyg) + kinhO2dn)

!  Phytoplankton exudation of semilabile DON to DON
                Bio(i,k,iDON_)=Bio(i,k,iDON_)+EsDON(ng)*                &
     &                        (N_Flux_NewProd_NO3_lim+N_Flux_RegProd_lim)
!  Phytoplankton exudation of labile DON to NH4
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+(f_NTR+f_DNF)*ElDON(ng)*  &
     &                        (N_Flux_NewProd_NO3_lim+N_Flux_RegProd_lim)

!  Phytoplankton exudation of semilabile DOP to DOP
                Bio(i,k,iDOP_)=Bio(i,k,iDOP_)+EsDON(ng)*                &
     &                         N_Flux_NewProd_P_lim*NxP_P  
!  Phytoplankton exudation of labile DOP to PO4
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ElDON(ng)*                &
     &                         N_Flux_NewProd_P_lim*NxP_P
      

#if defined OXYGEN || defined H_SULF

! WC_DENITRIFICATION
                f_NTR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg) + kinhO2dn)
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
     &                                N_Flux_NewProd_NO3_lim*rOxNO3+    &
     &                                N_Flux_NewProd_NFix*rOxNH4+       & 
     &                                N_Flux_RegProd_lim*rOxNH4         &
     &                         -rOxNH4*f_NTR*ElDON(ng)*                 &
     &                         (N_Flux_NewProd_NO3_lim+N_Flux_RegProd_lim)
                Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
#endif
#ifdef CARBON
!
!  Total inorganic carbon (CO2) uptake during phytoplankton growth.
!
                cff1=PhyCN(ng)*(N_Flux_Lim+N_Flux_NewProd_NFix)
                Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-cff1

# ifdef TALK_NONCONSERV
!
!  Account for the uptake of NO3 on total alkalinity.
!
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+N_Flux_NewProd_NO3_lim+   &
     &                         N_Flux_NewProd_NFix
# endif
#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2 ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2- + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
#if defined OXYGEN || defined H_SULF
                fac2=MAX(Bio(i,k,iOxyg),0.0_r8)     ! O2 max
                fac3=MAX(fac2/(3.0_r8+fac2),0.0_r8) ! MM for O2 dependence
                fac1=dtdays*NitriR(ng)*fac3
                fac12=dtdays*NitriR2*fac3
#else
                fac1=dtdays*NitriR(ng)
#endif
!
!...........Through Nitrosomonas
!          
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                cff2=1.0_r8-MAX(0.0_r8,cff1)
                cff3=fac1*cff2
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO2_)=Bio(i,k,iNO2_)+N_Flux_Nitrifi
!
!...........Through Nitrobacter
!
                cff4=(PAR-I_thNO2)/(D_p5NO2+PAR-2.0_r8*I_thNO2)       
                cff5=1.0_r8-MAX(0.0_r8,cff4)
                cff6=fac12*cff5
                Bio(i,k,iNO2_)=Bio(i,k,iNO2_)/(1.0_r8+cff6)
                N_Flux_Nitrifi2=Bio(i,k,iNO2_)*cff6
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi2

#if defined OXYGEN || defined H_SULF
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-1.5_r8*N_Flux_Nitrifi-    &
     &                         0.5_r8*N_Flux_Nitrifi2
#endif
#ifdef DIAGNOSTICS_BIO
                DiaBio3d(i,j,k,iNit1)=DiaBio3d(i,j,k,iNit1)+            &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &   
# endif
     &                              N_Flux_Nitrifi*fiter              
                DiaBio3d(i,j,k,iNit2)=DiaBio3d(i,j,k,iNit2)+            &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &     
# endif     
     &                              N_Flux_Nitrifi2*fiter
#endif
#if defined CARBON && defined TALK_NONCONSERV
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)-N_Flux_Nitrifi
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
            ELSE
              DO k=N(ng),1,-1
#if defined OXYGEN || defined H_SULF
                fac2=MAX(Bio(i,k,iOxyg)-0.3_r8,0.0_r8)            ! O2 max
                fac3=MAX(fac2/(1.0_r8+fac2),0.0_r8)        ! MM for O2 dependence
                cff3=dtdays*NitriR(ng)*fac3
                cff32=dtdays*NitriR2*fac3
#else
                cff3=dtdays*NitriR(ng)
                cff32=dtdays*NitriR2
#endif
!
!.......Through Nitrosomonas
!
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO2_)=Bio(i,k,iNO2_)+N_Flux_Nitrifi
!
!.......Through Nitrobacter
!
                Bio(i,k,iNO2_)=Bio(i,k,iNO2_)/(1.0_r8+cff32)
                N_Flux_Nitrifi2=Bio(i,k,iNO2_)*cff32
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi2
#if defined OXYGEN || defined H_SULF
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-1.5_r8*N_Flux_Nitrifi-    &
     &                                0.5_r8*N_Flux_Nitrifi2
#endif
#ifdef DIAGNOSTICS_BIO
              DiaBio3d(i,j,k,iNit1)=DiaBio3d(i,j,k,iNit1)+              &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &     
# endif 
     &                              N_Flux_Nitrifi*fiter
              DiaBio3d(i,j,k,iNit2)=DiaBio3d(i,j,k,iNit2)+              &
# ifdef WET_DRY
     &                              rmask_full(i,j)*                    &      
# endif
     &                              N_Flux_Nitrifi2*fiter
#endif
#if defined CARBON && defined TALK_NONCONSERV
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)-N_Flux_Nitrifi
#endif
              END DO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!
!  Update from Azhar et al. (2014): diazotroph mortality coefficient
!  was previously set to a value of 0.05, but it is now set to equal
!  the phytoplankton mortality rate coefficient. -KGH, Mar. 2018.
!-----------------------------------------------------------------------
!
          fac1=dtdays*ZooGR(ng)
          cff2=dtdays*PhyMR(ng)                    ! phytopl. mortality
          cff22=dtdays*PhyMR(ng)                   ! diazotroph mortality
          DO k=1,N(ng)
            DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
!
!             Re-define fac1 for temperature-dependent grazing (Fei Da, 201705).
!             0.0875 is to achieve a `Q_10' of 2.4; see Lomas et al. 2002.
              fac1=dtdays * ZooGR(ng) * exp( 0.0875_r8 * Bio(i,k,itemp))
!psl Q10=2.1  fac1=dtdays * ZooGR(ng) * exp( 0.0742_r8 * Bio(i,k,itemp))
              cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff11=0.84_r8*fac1*Bio(i,k,iZoop)*Bio(i,k,iDiaz)/         & 
     &             (K_Phy(ng)+Bio(i,k,iDiaz)*Bio(i,k,iDiaz))
              cff3=1.0_r8/(1.0_r8+cff1)
              cff33=1.0_r8/(1.0_r8+cff11)
              Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
              Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
              Bio(i,k,iDiaz)=cff33*Bio(i,k,iDiaz)
!
!  Phytoplankton assimilated to zooplankton and egested to small
!  detritus.
!
!  Update from Azhar et al. (2014):
!  Diazotrophs have a higher N:P ratio than phytoplankton. In order to
!  simplify calculations, the uptake of nitrogen from diazotrophs into
!  zooplankton is limited according to the N:P ratio of phytoplankton.
!  Thus, zooplankton biomass will maintain a constant N:P ratio equal to
!  that of phytoplankton while excess nitrogen (N_Assim_Excess) is 
!  added to the small detritus nitrogen pool (SDeN). The corresponding
!  amount of carbon is also added to the small detritus carbon pool (SDeC).
!  
!  Additionally, the value N_Flux_Egest (the nitrogen egested from
!  phytoplankton) is now split into four distinct values: N_Flux_EgestP
!  (nitrogen egested from phytoplankton), N_Flux_EgestD (nitrogen egested
!  from diazotrophs), P_Flux_EgestP (phosphorus egested from phytoplankton),
!  and P_Flux_EgestD (phosphorus egested from diazotrophs). As with the
!  N_Assim_Excess calculation, this splitting is needed to account for the
!  different N:P ratios of phytoplankton and diazotrophs. -KGH, Mar. 2018.
!
              N_Flux_Assim=(cff1*Bio(i,k,iPhyt)+cff11*Bio(i,k,iDiaz))*  &
     &                     ZooAE_N(ng)
              P_Flux_Assim=(cff1*Bio(i,k,iPhyt)*NxP_P+                  &
     &                     cff11*(Bio(i,k,iDiaz)/rNxP_D))*ZooAE_N(ng)
              N_Assim_Excess=N_Flux_Assim-rNxP_P*P_Flux_Assim
              N_Flux_EgestP=cff1*Bio(i,k,iPhyt)*(1.0_r8-ZooAE_N(ng))
              N_Flux_EgestD=cff11*Bio(i,k,iDiaz)*(1.0_r8-ZooAE_N(ng))
              P_Flux_EgestP=cff1*Bio(i,k,iPhyt)*(1.0_r8-ZooAE_N(ng))*NxP_P
              P_Flux_EgestD=cff11*(Bio(i,k,iDiaz)/rNxP_D)*              &
     &                      (1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+P_Flux_Assim*rNxP_P
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                            &
     &                       N_Flux_EgestP+N_Flux_EgestD+N_Assim_Excess
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
     &                       P_Flux_EgestP+P_Flux_EgestD
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              N_Flux_PmortalD=cff22*MAX(Bio(i,k,iDiaz)-                 & 
     &                        PhyMin(ng),0.0_r8)   

              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
              Bio(i,k,iDiaz)=Bio(i,k,iDiaz)-N_Flux_PmortalD

              Bio(i,k,iChlo)=Bio(i,k,iChlo)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                            &
     &                       N_Flux_Pmortal+N_Flux_PmortalD
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
     &                      NxP_P*N_Flux_Pmortal+(N_Flux_PmortalD/rNxP_D)
#ifdef CARBON
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)+                            &
     &                       PhyCN(ng)*(N_Flux_EgestP+N_Flux_EgestD+    &
     &                       N_Flux_Pmortal+N_Flux_PmortalD+            &
     &                       N_Assim_Excess)+                           &
     &                       (PhyCN(ng)-ZooCN(ng))*P_Flux_Assim*rNxP_P
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!
!  Update from Azhar et al. (2014): nitrogen and phosphorus fluxes are
!  now separated, and fluxes associated with phytoplankton are separated
!  from fluxes associated with diazotrophs due to their different N:P
!  ratios. For example, the value that was previously N_Flux_Excret is
!  now four values: N_Flux_Zexcret_Phyt, N_Flux_Zexcret_Diaz,
!  P_Flux_Zexcret_Phyt, and P_Flux_Zexcret_Diaz. Additionally, a flux
!  of phosphorus from zooplankton metabolism is now added to the small
!  detritus phosphorus pool. This missing source had caused a mass
!  imbalance in phosphorus in the 2014 version of the code. Finally,
!  the carbon flux associated with zooplankton excretion is now added
!  to the small detritus carbon pool rather than the dissolved inorganic
!  carbon pool. -KGH, Mar. 2018.
!-----------------------------------------------------------------------
!
          cff1=dtdays*ZooBM(ng)
          fac2=dtdays*ZooMR(ng)
          fac3=dtdays*ZooER(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  Zooplankton excretion
!
              fac33=Bio(i,k,iPhyt)
              fac34=Bio(i,k,iDiaz)
              fac1=ZooAE_N(ng)*fac3*fac33*fac33/(K_Phy(ng)+fac33*fac33)
              fac12=ZooAE_N(ng)*fac3*fac34*fac34/(K_Phy(ng)+fac34*fac34)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+fac1)
              N_Flux_Zexcret_Phyt=Bio(i,k,iZoop)*fac1
              P_Flux_Zexcret_Phyt=N_Flux_Zexcret_Phyt*NxP_P
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+fac12)
              N_Flux_Zexcret_Diaz=Bio(i,k,iZoop)*fac12
              P_Flux_Zexcret_Diaz=N_Flux_Zexcret_Diaz/rNxP_P
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                            &
     &                       N_Flux_Zexcret_Phyt+N_Flux_Zexcret_Diaz
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
     &                       P_Flux_Zexcret_Phyt+P_Flux_Zexcret_Diaz
!
!  Zooplankton mortality
!
              cff2=fac2*Bio(i,k,iZoop)
              fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff3=fac1*ZooAE_N(ng)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+cff2+cff3)
              N_Flux_Zmortal=Bio(i,k,iZoop)*cff2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+N_Flux_Zmortal*NxP_P
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
              N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+N_Flux_Zmetabo*NxP_P
#if defined OXYGEN || defined H_SULF
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret_Phyt+&
     &                       N_Flux_Zexcret_Diaz)
#endif
#ifdef CARBON
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)+                            &
     &                       ZooCN(ng)*(N_Flux_Zmortal+                 &
     &                       N_Flux_Zexcret_Phyt+N_Flux_Zexcret_Diaz)
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       ZooCN(ng)*(N_Flux_Zmetabo)
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!
!  Update from Azhar et al. (2014): the values cff1P, cff1C, cff2P, and
!  cff2C are now included to maintain separate coagulation fluxes of
!  phosphorus and carbon.
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*(Bio(i,k,iSDeN)+Bio(i,k,iPhyt))
              cff1P=fac1*(Bio(i,k,iSDeP)+NxP_P*Bio(i,k,iPhyt))
              cff2=1.0_r8/(1.0_r8+cff1)
              cff2P=1.0_r8/(1.0_r8+cff1P)
#ifdef CARBON
              cff1C=fac1*(Bio(i,k,iSDeC)+PhyCN(ng)*Bio(i,k,iPhyt))
              cff2C=1.0_r8/(1.0_r8+cff1C)
#endif
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
              Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2P              
              N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
              N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
              P_Flux_CoagP=N_Flux_CoagP*NxP_P
              P_Flux_CoagD=Bio(i,k,iSDeP)*cff1P

              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
     &                       N_Flux_CoagP+N_Flux_CoagD
              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                            &
     &                       P_Flux_CoagP+P_Flux_CoagD 

#ifdef CARBON
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)*cff2C
              C_Flux_CoagP=N_Flux_CoagP*PhyCN(ng)
              C_Flux_CoagD=Bio(i,k,iSDeC)*cff1C
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)+C_Flux_CoagP+C_Flux_CoagD
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
#if defined OXYGEN || defined H_SULF
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=MAX(Bio(i,k,iOxyg)-1.0_r8,0.0_r8)     ! O2 off max
              fac2=MAX(Bio(i,k,iNO3_)-1.0_r8,0.0_r8)     ! NO3_off max
!              fac3=MAX(Bio(i,k,iNO2_)-2.0_r8,0.0_r8)     ! NO2_off max
              fac3=MAX(Bio(i,k,iNO2_)-0.0_r8,0.0_r8)     ! NO2_off max
              oxlim=MAX(fac1/(kO2+fac1),0.0_r8)          ! MM for O2 dependence
              dlim1=MAX((fac2/(kNO3+fac2))*LNO3*                        &
     &              (kinhO2dn/(fac1+kinhO2dn)),0.0_r8)   ! NO3 reduction
!
!  denit:loss NO2 15% (Canfield, 2010).
!
              dlim2=MAX((fac3/(kNO2+fac3))*LNO2*eNO2*                   &
     &              (kinhO2dn/(fac1+kinhO2dn)),0.0_r8)
 
              srlim=MAX((kinhO2an/(kinhO2an+fac1))*                     &
     &              (kinhNO3an/(kinhNO3an+fac2)),0.0_r8) ! MM SR dependence 
                
              tolim=oxlim+dlim1+dlim2+srlim

              fox=oxlim/tolim
              fd1=dlim1/tolim
              fd2=dlim2/tolim
              fsr=srlim/tolim

              cff1=dtdays*SDeRRN(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)
              cff4=1.0_r8/(1.0_r8+cff3)
              cff5=106.0_r8/8.0_r8                     ! rNO3:NH3
              cff6=53.0_r8/16.0_r8                     ! rH2S:N
              cff7=106.0_r8/12.0_r8                    ! rNO2:NH3
              cff8=53.0_r8/12.0_r8                     ! rdenit N2:NH3

!  WC_NITRIFICATION && defined WC_DENITRIFICATION
              f_NTR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg) + kinhO2dn)
              f_DNF = kinhO2dn/(Bio(i,k,iOxyg) + kinhO2dn)

              cff1  = dtdays*SDeNSR(ng)
!             Temperature-dependent solubilization (Fei Da 201705).
!             0.0875 corresponds to a `Q_10' of 2.4; Lomas et al. 2002.
              cff1  = cff1 * exp( 0.0875_r8 * Bio(i,k,itemp) )
!psl Q10=2.1  cff1  = cff1 * exp( 0.0742_r8 * Bio(i,k,itemp) )
              cff2  = deltN(ng)*cff1
              cff3a = (1-deltN(ng))*cff1*f_NTR
              cff3b = (1-deltN(ng))*cff1*f_DNF
              cff3  = cff3a + cff3b
              cff14  = 1.0_r8/(1.0_r8+cff2+cff3)

              cff4  = dtdays*LDeNSR(ng)
!             Temperature-dependent solubilization (Fei Da 201705).
!             0.0875 corresponds to a `Q_10' of 2.4; Lomas et al. 2002.
              cff4  = cff4 * exp( 0.0875_r8 * Bio(i,k,itemp) )
!psl Q10=2.1  cff4  = cff4 * exp( 0.0742_r8 * Bio(i,k,itemp) )
              cff12  = deltN(ng)*cff4
              cff13a = (1-deltN(ng))*cff4*f_NTR
              cff13b = (1-deltN(ng))*cff4*f_DNF
              cff13  = cff13a + cff13b
              cff15  = 1.0_r8/(1.0_r8+cff12+cff13)

              rDON  = dtdays * a0N(ng) * exp(0.0875_r8 * Bio(i,k,itemp))
!psl Q10=2.1  rDON  = dtdays * a0N(ng) * exp(0.0742_r8 * Bio(i,k,itemp))
              cff16a = rDON*f_NTR
              cff16b = rDON*f_DNF
              cff16  = cff16a + cff16b
              cff17 = 1.0_r8/(1.0_r8+cff16)

              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff14
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff15
              Bio(i,k,iDON_)=Bio(i,k,iDON_)*cff17

              N_Flux_Remine =Bio(i,k,iSDeN)*cff3+Bio(i,k,iLDeN)*cff13
              N_Flux_RemineD=Bio(i,k,iDON_)*cff16

              N_Flux_SolubS =Bio(i,k,iSDeN)*cff2
              N_Flux_SolubL =Bio(i,k,iLDeN)*cff12

              Bio(i,k,iNH4_)=Bio(i,k,iNH4_) + N_Flux_RemineD            &
     &                      +N_Flux_Remine
              Bio(i,k,iDON_)=Bio(i,k,iDON_)                             &
     &                      +N_Flux_SolubS + N_Flux_SolubL

              N_Flux_RemineSa=Bio(i,k,iSDeN)*cff3a 
              N_Flux_RemineLa=Bio(i,k,iLDeN)*cff13a
              N_Flux_RemineDa=Bio(i,k,iDON_)*cff16a

              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &      (N_Flux_RemineDa + N_Flux_RemineSa + N_Flux_RemineLa)*rOxNH4*fox
              Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
              
!              cff1=dtdays*SDeRRN(ng)
!              cff2=1.0_r8/(1.0_r8+cff1)
!
! Remineralization rate of N detrital
!
!              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
!              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
!              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3

! Temperature dependent remineralization of semi-labile DON
!              rDON=a0N(ng)*EXP(0.07_r8*Bio(i,k,itemp))
!              !cff5=rDON*dtdays
!              !cff6=1.0_r8/(1.0_r8+cff5)
!
!
!              Bio(i,k,iDON_)=Bio(i,k,iDON_) &
!            &  *(1.0_r8/(1.0_r8+(rDON*dtdays)))
!
!              N_Flux_RemineD=Bio(i,k,iDON_)*(rDON*dtdays)
!              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+                            &
!     &                       N_Flux_RemineD                             &
!     &                       +(1.0_r8-deltN(ng))*N_Flux_Remine
!              Bio(i,k,iDON_)=Bio(i,k,iDON_)+deltN(ng)*N_Flux_Remine
!              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
!     &                       N_Flux_RemineD*rOxNH4                      &
!     &                       -(1.0_r8-deltN(ng))*                       &
!     &                       N_Flux_Remine*rOxNH4
!              Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)     

                          
!
! Remineralization rate of P detrital
!
!              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2             
!              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
!              P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3 
!              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
!              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
!              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4*fox
!              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_Remine*fd1*cff5

! Temperature dependent remineralization of semi-labile DOP
!              rDON=0.0106_r8*EXP(0.07_r8*Bio(i,k,itemp))
!              !cff5=rDON*dtdays
!              !cff6=1.0_r8/(1.0_r8+cff5)
!
!
!              Bio(i,k,iDOP_)=Bio(i,k,iDOP_) &
!            &  *(1.0_r8/(1.0_r8+(rDON*dtdays)))
!
!              P_Flux_RemineD=Bio(i,k,iDOP_)*(rDON*dtdays)
!              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                            &
!     &                       P_Flux_RemineD                             &
!     &                       +(1.0_r8-deltN(ng))*P_Flux_Remine
!              Bio(i,k,iDOP_)=Bio(i,k,iDOP_)+deltN(ng)*P_Flux_Remine

! Temperature dependent solubilization of semi-labile DOP
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff14
              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff15
              Bio(i,k,iDOP_)=Bio(i,k,iDOP_)*cff17

              P_Flux_Remine =Bio(i,k,iSDeP)*cff3+Bio(i,k,iLDeP)*cff13
              P_Flux_RemineD=Bio(i,k,iDOP_)*cff16

              P_Flux_SolubS =Bio(i,k,iSDeP)*cff2
              P_Flux_SolubL =Bio(i,k,iLDeP)*cff12

              Bio(i,k,iPO4_)=Bio(i,k,iPO4_) + P_Flux_RemineD            &
     &                      +P_Flux_Remine
              Bio(i,k,iDOP_)=Bio(i,k,iDOP_)                             &
     &                      +P_Flux_SolubS + P_Flux_SolubL

              P_Flux_RemineSa=Bio(i,k,iSDeP)*cff3a 
              P_Flux_RemineLa=Bio(i,k,iLDeP)*cff13a
              P_Flux_RemineDa=Bio(i,k,iDOP_)*cff16a

!              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
!     &      (P_Flux_RemineDa + P_Flux_RemineSa + P_Flux_RemineLa)*rOxPO4*fox
!              Bio(i,k,iOxyg)=MAX(Bio(i,k,iOxyg),0.0_r8)
!              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_Remine*fd1*cff5

!  WC_NITRIFICATION && defined WC_DENITRIFICATION
! NO3 lost
              f_DNF = kinhO2dn/(Bio(i,k,iOxyg)+kinhO2dn)
!psl20180225  The role of f_WC is to bound the NO3 loss and to prevent negative
!               NO3 values (Feng et al. 2015, page-6). However, f_WC is also
!               introducing an artificial non-linear response to increases/
!               decreases of NO3 whenever F_WC<f_DNF is satisfied. So I'm
!               getting rid of f_WC and use a different approach to bound NO3.
!psl20180225  f_WC  = Bio(i,k,iNO3_)/(Bio(i,k,iNO3_) + K_WNO3)
              cff18  = 5.30_r8*    f_DNF       ! eta_DNF=5.3, Table-A4.
!psl20180225  cff18  = 5.30_r8*min(f_DNF,f_WC)

!             Use rates `cff1,cff4' defined above to obtain a
!               temperature-dependent remineralization (Fei Da 201705).
!               0.0875 corresponds to a `Q_10' of 2.4; Lomas et al. 2002.
              cff3c = (1._r8-deltN(ng))*cff1*cff18
              cff6c = (1._r8-deltN(ng))*cff4*cff18

              rDON  = a0N(ng)*exp(0.0875_r8*Bio(i,k,itemp))
!psl Q10=2.1  rDON  = a0N(ng)*exp(0.0742_r8*Bio(i,k,itemp))
              cff9c = rDON*dtdays*cff18

              N_Flux_RemineSc=Bio(i,k,iSDeN)*cff3c 
              N_Flux_RemineLc=Bio(i,k,iLDeN)*cff6c
              N_Flux_RemineDc=Bio(i,k,iDON_)*cff9c

!psl20180225  Bound the NO3 loss to prevent negative NO3 values:
              N_Flux_RemineDc = N_Flux_RemineSc + N_Flux_RemineLc       &
     &                        + N_Flux_RemineDc
              N_Flux_RemineSc = 0._r8
              N_Flux_RemineLc = 0._r8
              N_Flux_RemineDc = min(N_Flux_RemineDc, Bio(i,k,iNO3_)-eps)

              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-                            &
     &             N_Flux_RemineDc-N_Flux_RemineSc-N_Flux_RemineLc

!
! Anammox
!
              fac4=MAX(Bio(i,k,iNH4_),0.0_r8)          ! NH4
              fac5=MAX(Bio(i,k,iNO2_),0.0_r8)          ! NO2
              fanmx=dtdays*KANMX

              ranmx=1.00_r8*fanmx*fac4*fac5*(kinhO2dn/(fac1+kinhO2dn))  ! Anammox 85% NO2 loss

              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-ranmx
              Bio(i,k,iNO2_)=Bio(i,k,iNO2_)+N_Flux_Remine*fd1*cff5-     &             
     &                       N_Flux_Remine*fd2*cff7-ranmx       
              Bio(i,k,iN2__)=Bio(i,k,iN2__)+2.0_r8*ranmx+               &
     &                       N_Flux_Remine*fd2*cff7
# ifdef H_SULF
              Bio(i,k,iH2S_)=Bio(i,k,iH2S_)+N_Flux_Remine*fsr*cff6
              Bio(i,k,iSO4_)=Bio(i,k,iSO4_)-N_Flux_Remine*fsr*cff6
#  ifdef DIAGNOSTICS_BIO
              DiaBio3d(i,j,k,iOxic)=DiaBio3d(i,j,k,iOxic)+              &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                              N_Flux_Remine*fox*fiter            ! in N (NH4 production)
              DiaBio3d(i,j,k,iSRRa)=DiaBio3d(i,j,k,iSRRa)+              &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif     
     &                              N_Flux_Remine*fsr*cff6*fiter       ! in S (H2S production)
              DiaBio3d(i,j,k,iDno2)=DiaBio3d(i,j,k,iDno2)+              &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif     
     &                              N_Flux_Remine*fd2*cff8*fiter       ! in N (N2 production)
              DiaBio3d(i,j,k,iDno3)=DiaBio3d(i,j,k,iDno3)+              &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    & 
#   endif    
     &                              N_Flux_Remine*fd1*cff5*fiter       ! in N (NO2 production)
              DiaBio3d(i,j,k,iAnx_)=DiaBio3d(i,j,k,iAnx_)+              &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif     
     &                              ranmx*fiter                        ! in N (N2 production)
#  endif
# endif
            END DO
          END DO 
             
#else

          cff1=dtdays*SDeRRN(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRN(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
            END DO
          END DO            
            
#endif 
#if defined OXYGEN || defined H_SULF
!
!-----------------------------------------------------------------------
!  Sulfide oxidation by O2 and NO3 (Millero, 1991), (Konovalov,2006).
!  Al Azhar - Halifax, NS. (March 2011).
!  Modifications by F. Daryabor, 2017. IGN. UCPH.
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=dtdays*KH2SO
              cff2=dtdays*KH2SN1
              cff3=2.0_r8/8.0_r8              ! NO3/H2S
              cff4=dtdays*KH2SN2
              cff5=3.0_r8/8.0_r8              ! NO2/H2S
              fac1=MAX(Bio(i,k,iH2S_),0.0_r8)
!              fac12=MAX(Bio(i,k,iSO4_),0.0_r8)
              fac2=MAX(Bio(i,k,iOxyg),0.0_r8)
              fac22=MAX(Bio(i,k,iOxyg)/(Bio(i,k,iOxyg)+kSO2),0.0_r8)
              fac3=MAX(Bio(i,k,iNO3_)/(Bio(i,k,iNO3_)+kSNO3),0.0_r8)
              fac4=MAX(Bio(i,k,iNO2_)/(Bio(i,k,iNO2_)+kSNO2),0.0_r8)
!
              Bio(i,k,iH2S_)=Bio(i,k,iH2S_)-(cff1*fac1*fac22)-          &
     &                       ((cff2*fac1*fac3*cff3)*                    &
     &                       (kinhO2dn/(fac2+kinhO2dn)))-               &
     &                       ((cff4*fac1*fac4*cff5)*                    &
     &                       (kinhO2dn/(fac2+kinhO2dn)))
              Bio(i,k,iSO4_)=Bio(i,k,iSO4_)+(cff1*fac1*fac22)+         &
     &                       ((cff2*fac1*fac3*cff3)*                   &
     &                       (kinhO2dn/(fac2+kinhO2dn)))+               &
     &                       ((cff4*fac1*fac4*cff5)*                    &
     &                       (kinhO2dn/(fac2+kinhO2dn)))
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &                       (2.0_r8*cff1*fac1*fac22)             
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-((cff2*fac1*fac3)*          &
     &                       (kinhO2dn/(fac2+kinhO2dn)))
              Bio(i,k,iNO2_)=Bio(i,k,iNO2_)-((cff4*fac1*fac4)*          &
     &                       (kinhO2dn/(fac2+kinhO2dn)))+               &
     &                       ((cff2*fac1*fac3)*                         &
     &                       (kinhO2dn/(fac2+kinhO2dn)))
              Bio(i,k,iN2__)=Bio(i,k,iN2__)+((cff4*fac1*fac4)*          &
     &                       (kinhO2dn/(fac2+kinhO2dn)))

# ifdef DIAGNOSTICS_BIO
             DiaBio3d(i,j,k,iSOx_)=DiaBio3d(i,j,k,iSOx_)+               &
#  ifdef WET_DRY  
     &                              rmask_full(i,j)*                    &
#  endif   
     &                              cff1*fac1*fac22*fiter                     ! in S
              DiaBio3d(i,j,k,iSNx_)=DiaBio3d(i,j,k,iSNx_)+              &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif     
     &                              ((cff2*fac1*fac3)*                  &
     &                              (kinhO2dn/(fac2+kinhO2dn)))*fiter         ! in N
              DiaBio3d(i,j,k,iSNy_)=DiaBio3d(i,j,k,iSNy_)+              &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &     
#  endif     
     &                              ((cff4*fac1*fac4)*                  &
     &                              (kinhO2dn/(fac2+kinhO2dn)))*fiter           ! in N
# endif
            END DO
          END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          k=N(ng)
          DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
            SchmidtN_Ox=1638.0_r8-                                      &
     &                  Bio(i,k,itemp)*(81.83_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (1.483_r8-                      &
     &                                   Bio(i,k,itemp)*0.008004_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
            SchmidtN_Ox=1953.4_r8-                                      &
     &                  Bio(i,k,itemp)*(128.0_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (3.9918_r8-                     &
     &                                   Bio(i,k,itemp)*0.050091_r8))
# endif

            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
            TS=LOG((298.15_r8-Bio(i,k,itemp))/                          &
     &             (273.15_r8+Bio(i,k,itemp)))
            AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+          &
     &             Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
            O2satu=l2mol*EXP(AA)*O2palFrac(ng)
!
!  Add in O2 gas exchange.
!
            O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
     &                     O2_Flux*Hz_inv(i,k)
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iO2fx)=DiaBio2d(i,j,iO2fx)+                    &
#  ifdef WET_DRY
     &                          rmask_full(i,j)*                        &
#  endif
     &                          O2_Flux*fiter
# endif

          END DO
#endif

#ifdef CARBON
!
!-----------------------------------------------------------------------
!  Allow different remineralization rates for detrital C and detrital N.
!-----------------------------------------------------------------------
!
          cff1=dtdays*SDeRRC(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRC(ng)
          cff4=1.0_r8/(1.0_r8+cff3)

          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDeC)=Bio(i,k,iSDeC)*cff2
              Bio(i,k,iLDeC)=Bio(i,k,iLDeC)*cff4
              C_Flux_Remine=Bio(i,k,iSDeC)*cff1+Bio(i,k,iLDeC)*cff3
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+C_Flux_Remine
            END DO
          END DO

!
!-----------------------------------------------------------------------
!  Surface CO2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface.   
!
          k=N(ng)
# ifdef pCO2_RZ
          CALL pCO2_water_RZ (Istr, Iend, LBi, UBi, LBj, UBj,           &
     &                        IminS, ImaxS, j, DoNewton,                &
#  ifdef MASKING
     &                        rmask,                                    &
#  endif
     &                        Bio(IminS:,k,itemp), Bio(IminS:,k,isalt), &
     &                        Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk), &
     &                        pH, pCO2)
# else
          CALL pCO2_water (Istr, Iend, LBi, UBi, LBj, UBj,              &
     &                     IminS, ImaxS, j, DoNewton,                   &
#  ifdef MASKING
     &                     rmask,                                       &
#  endif
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk),    &
     &                     Bio(IminS:,k,iPO4_), 0.0_r8, pH, pCO2)  
# endif
!
!  Compute surface CO2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          DO i=Istr,Iend
!
!  Compute CO2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)**2+Vwind(i,j)**2
# else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
            SchmidtN=Acoef-                                             &
     &               Bio(i,k,itemp)*(Bcoef-                             &
     &                               Bio(i,k,itemp)*(Ccoef-             &
     &                               Bio(i,k,itemp)*Dcoef))
            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN)
!
!  Calculate CO2 solubility [mol/(kg.atm)] using Weiss (1974) formula.
!
            TempK=0.01_r8*(Bio(i,k,itemp)+273.15_r8)
            CO2_sol=EXP(A1+                                             &
     &                  A2/TempK+                                       &
     &                  A3*LOG(TempK)+                                  &
     &                  Bio(i,k,isalt)*(B1+TempK*(B2+B3*TempK)))
!
!  Add in CO2 gas exchange.
!
!            CALL caldate (r_date, tdays(ng), year, yday, month, iday,   &
!     &                    hour)
             CALL caldate (tdays(ng), yy_i=year, yd_r8=yday, mm_i=month, dd_i=iday, h_r8=hour)
            pmonth=2003.0_r8-1951.0_r8+yday/365.0_r8
            CO2_Flux=cff3*CO2_sol*(pCO2air(ng)-pCO2(i))
            Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                              &
     &                     CO2_Flux*Hz_inv(i,k)
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iCOfx)=DiaBio2d(i,j,iCOfx)+                    &
#  ifdef WET_DRY
     &                          rmask_full(i,j)*                        &
#  endif
     &                          CO2_Flux*fiter
            DiaBio2d(i,j,ipCO2)=pCO2(i)
#  ifdef WET_DRY
            DiaBio2d(i,j,ipCO2)=DiaBio2d(i,j,ipCO2)*rmask_full(i,j)
#  endif
# endif
          END DO
#endif
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink

            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
!  Update from Azhar et al. (2014): a benthic remineralization flux from
!  diazotrophs to the dissolved inorganic carbon pool is now included;
!  this had been omitted in the 2014 version of the code. -KGH, Mar. 2018.
!
            cff2=4.0_r8/16.0_r8
# if defined OXYGEN || defined H_SULF
            cff3=115.0_r8/16.0_r8          ! R_oxic O2
            cff4=106.0_r8/16.0_r8
            cff5=16.0_r8/16.0_r8
            cff6=106.0_r8/8.0_r8           ! R_denit NO3
            cff66=212.0_r8/16.0_r8         ! R_anamx Bo
            cff7=53.0_r8/16.0_r8           ! R_S red 
            cff8=196.0_r8/12.0_r8          ! R_anamx Bo
            cff9=106.0_r8/12.0_r8          ! R_denit NO2
            cff10=53.0_r8/12.0_r8          ! R_denit N2          
# endif

          DO i=Istr,Iend
            IF ((ibio.eq.iPhyt).or.                                     &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN).or.                                     &
     &          (ibio.eq.iSDeC).or.                                     &
     &          (ibio.eq.iLDeC).or.                                     &
     &          (ibio.eq.iSDeP).or.                                     &
     &          (ibio.eq.iLDeP).or.                                     &
     &          (ibio.eq.iDiaz)) THEN

#ifdef USECOS_BURIAL
!
!  The POM flux reaching the seafloor is partially resuspended in the 
!  lower water column level within the small detritus compartment,
!  assuming fragmentation of large particles. 
!  Druon 1999 (pg 32) assumes the resuspension rate (ReSuspR: 
!  nondimensional) is a function of the bottom friction velocity and 
!  the critical friction velocity for resuspension of Particulate 
!  Organic Matter: ReSuspR = (Ustarb/Ustar_crit)^2 where Ustarb is 
!  computed from the model bottom stress, bustr. The critical stress 
!  proposed by Peterson 1999 (Aquacultural Engineering, 21, 85-111) 
!  is 0.01 Pa equals Ustar_crit=3.1E-3 m s-1 for seawater rho = 1027. 
!  The factor 9.61E-6 is 3.1E-3^2.
!
!  The sinking flux FC(i,0) is in units millimole m-2 because it is 
!  concentration times layer thickness and has been integrated over 
!  model dt. When converting this to a quantity to be added to a 
!  tracer variable (concentration in millimole m-3) must divide
!  by layer thickness
!
                Ustarb=SQRT(SQRT((0.5_r8*(bustr(i,j)+bustr(i+1,j)))**2+ &
     &                 (0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))**2))
                ReSuspR=MIN(Ustarb**2/9.61E-6_r8,1.0_r8)
!
!  The particulate flux reaching the seafloor that will 
!  be remineralized/denitrified
                cff1=FC(i,0)*Hz_inv(i,1)*(1.0_r8-ReSuspR)      !mmol m-3
!
!  The complimentary resuspended fraction goes back to small detritus
                cff5=FC(i,0)*Hz_inv(i,1)*ReSuspR
!
                IF ((ibio.eq.iPhyt).or.                                 &
     &            (ibio.eq.iSDeN).or.                                   &
     &            (ibio.eq.iLDeN)) THEN
                  Bio(i,1,iSDeN)=Bio(i,1,iSDeN)+cff5
                ENDIF
                IF ((ibio.eq.iSDeP).or.(ibio.eq.iLDeP)) THEN 
                  Bio(i,1,iSDeP)=Bio(i,1,iSDeP)+cff5
                ENDIF
                IF ((ibio.eq.iSDeC).or.(ibio.eq.iLDeC)) THEN
                  Bio(i,1,iSDeC)=Bio(i,1,iSDeC)+cff5
                END IF
                IF (ibio.eq.iPhyt) THEN
                  Bio(i,1,iSDeC)=Bio(i,1,iSDeC)+cff5*PhyCN(ng)
                  Bio(i,1,iSDeP)=Bio(i,1,iSDeP)+cff5*NxP_P
                END IF    
                     
#else
  Original Fennel code:

  The particulate flux reaching the seafloor that will 
  be remineralized/denitrified
                cff1=FC(i,0)*Hz_inv(i,1)                      ! mmol m-3
#endif         
                END IF

              IF ((ibio.eq.iPhyt).or.                                   &
     &            (ibio.eq.iSDeN).or.                                   &
     &            (ibio.eq.iLDeN)) THEN
           
#ifdef USECOS_BURIAL
                IF (ibio.eq.iPhyt) THEN
                  CNbur = PhyCN(ng)
                ELSE
                  CNbur = 9.3_r8
                ENDIF
          !       cff5 = 0.01_r8          ! 1 percent of flux goes to DON
                fac1=CNbur*                                             &
     &               12.0_r8/1000.0_r8*                                 &
     &               cff1*Hz(i,j,1)*(365.0_r8/dtdays) ! gC m-2 year-1
                Cbe=MIN(0.75_r8,0.023_r8*fac1**0.5797_r8 )
# ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iNbur)=DiaBio2d(i,j,iNbur)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              Cbe*cff1*Hz(i,j,1)*fiter
!  The flux of nitrogen that reaches the seafloor
                DiaBio2d(i,j,iNbot)=DiaBio2d(i,j,iNbot)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              cff1*Hz(i,j,1)*fiter     
# endif   
!
!  cff1 will now be the unburied flux that is to be remineralized
!  Must multiply cff1 by Hz to convert to flux in meter-2 when saving
!  the denitrification diagnostic term. Leave in meter-3 when returning
!  C or N to dissolved tracer concentrations
                cff1=(1.0_r8-Cbe)*cff1                        ! mmol m-3
!
#endif /* BURIAL */             

# ifdef DENITRIFICATION
#  if defined OXYGEN || defined H_SULF
                fac1=MAX(Bio(i,1,iOxyg)-1.0_r8,0.0_r8)    ! O2 off max
                fac2=MAX(Bio(i,1,iNO3_)-1.0_r8,0.0_r8)    ! NO3_off max
!                fac3=MAX(Bio(i,1,iNO2_)-2.0_r8,0.0_r8)    ! NO2_off max
                fac3=MAX(Bio(i,1,iNO2_)-0.0_r8,0.0_r8)    ! NO2_off max
                oxlim=MAX(fac1/(kO2+fac1),0.0_r8)         ! MM O2 dependence
                dlim1=MAX((fac2/(kNO3+fac2))*LNO3*                      &
     &                (kinhO2dn/(fac1+kinhO2dn)),0.0_r8)  ! NO3 reduction
                dlim2=MAX((fac3/(kNO2+fac3))*LNO2*eNO2*                 &
     &                (kinhO2dn/(fac1+kinhO2dn)),0.0_r8)  ! denit:loss NO2
!
! MM SR dependence
!
                srlim=MAX((kinhO2an/(kinhO2an+fac1))*                   &
     &                (kinhNO3an/(kinhNO3an+fac2)),0.0_r8)                 
              
                tolim=oxlim+dlim1+dlim2+srlim

                fox=oxlim/tolim
                fd1=dlim1/tolim
                fd2=dlim2/tolim
                fsr=srlim/tolim

                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*(fox+fd1+fd2+fsr)         ! from ox,sr,df1,df2
                Bio(i,1,iNO3_)=Bio(i,1,iNO3_)-cff1*fd1*cff6                  ! consumed by denit+anamx
                Bio(i,1,iNO2_)=Bio(i,1,iNO2_)+cff1*fd1*cff6-cff1*fd2*cff9    ! produced by anamx
                Bio(i,1,iH2S_)=Bio(i,1,iH2S_)+cff1*fsr*cff7
                Bio(i,1,iSO4_)=Bio(i,1,iSO4_)-cff1*fsr*cff7
                Bio(i,1,iN2__)=Bio(i,1,iN2__)+cff1*fd2*cff9
                Bio(i,1,iDON_)=Bio(i,1,iDON_)+cff1*0.01_r8 ! 1/100
                Bio(i,1,iDOP_)=Bio(i,1,iDOP_)+cff1*0.01_r8 ! 1/100
#  endif
#  ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iDNIT)=DiaBio2d(i,j,iDNIT)+                &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                              fd2*cff9*cff1*Hz(i,j,1)*fiter
                DiaBio2d(i,j,iSRR_)=DiaBio2d(i,j,iSRR_)+                &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif   
     &                              fsr*cff7*cff1*Hz(i,j,1)*fiter
                DiaBio2d(i,j,iBFlx)=DiaBio2d(i,j,iBFlx)+                &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                              cff1*Hz(i,j,1)*fiter
                DiaBio2d(i,j,iSRRb)=DiaBio2d(i,j,iSRRb)+                &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                              cff1*Hz(i,j,1)*fsr*cff7*fiter
                DiaBio3d(i,j,1,iOxic)=DiaBio3d(i,j,1,iOxic)+            &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif 
     &                              cff1*fox*fiter
                DiaBio3d(i,j,1,iSRRa)=DiaBio3d(i,j,1,iSRRa)+            &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                              cff1*fsr*cff7*fiter
                DiaBio3d(i,j,1,iDno2)=DiaBio3d(i,j,1,iDno2)+            &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                               cff1*fd2*cff10*fiter
                DiaBio3d(i,j,1,iDno3)=DiaBio3d(i,j,1,iDno3)+            &
#   ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#   endif
     &                              cff1*fd1*cff6*fiter
#  endif
#  if defined OXYGEN || defined H_SULF
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3*fox 
#  endif
# else
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*(fox+fd1+fd2+fsr)
#  if defined OXYGEN || defined H_SULF
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
#  endif
# endif

            END IF


              IF ((ibio.eq.iPhyt).or.                                   &
     &            (ibio.eq.iSDeC).or.                                   &
     &            (ibio.eq.iLDeC)) THEN


#ifdef USECOS_BURIAL
                if ( ibio .eq. iPhyt ) then   !psl20170429
!                 cff1 currently represents the non-resuspended, non-buried
!                   phytoplankton flux. However, the Cbe calculation below
!                   requires the non-resuspended phytoplankton flux. Therefore
!                   you must rescale cff1:
                  cff1 = cff1 / (1._r8 - Cbe) !psl20170429
                end if
!  In the iPhyt case we are repeating the Cbe calculation above.
!  CNbur coefficient not required here because we are already in carbon
!  units except for phytoplankton, in which case use PhyCN ratio
                fac1=12.0_r8/1000.0_r8*                                 &
     &               cff1*Hz(i,j,1)*(365.0_r8/dtdays) ! gC m-2 year-1
                IF (ibio.eq.iPhyt) THEN
                  fac1=fac1*PhyCN(ng)
                ENDIF
                Cbe=MIN(0.75_r8,0.023_r8*fac1**0.5797_r8 )
!  Factor for converting phytoplankton nitrogen flux to carbon
                IF (ibio.eq.iPhyt) THEN
                  fac2=PhyCN(ng)
                ELSE
                  fac2=1.0_r8
                ENDIF
# ifdef DIAGNOSTICS_BIO
!  The flux of carbon that reaches the seafloor
                DiaBio2d(i,j,iCbot)=DiaBio2d(i,j,iCbot)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              fac2*cff1*Hz(i,j,1)*fiter
!  The flux of carbon that is buried
                DiaBio2d(i,j,iCbur)=DiaBio2d(i,j,iCbur)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              Cbe*fac2*cff1*Hz(i,j,1)*fiter
# endif
!
!  The fraction of the flux that is not buried
!  will be remineralized
                cff1=(1.0_r8-Cbe)*cff1                        ! mmol m-3
              END IF  
#endif /* BURIAL */
!             
              IF ((ibio.eq.iPhyt).or.                                   &
     &            (ibio.eq.iSDeP).or.                                   &
     &            (ibio.eq.iLDeP)) THEN
      if ( ibio .eq. iPhyt ) then   !psl20170429
!                We must rescale cff1 again because the Cbe calculation
!                below requires the non-resuspended phytoplankton flux. 
                  cff1 = cff1 / (1._r8 - Cbe) 
                end if
!  In the iPhyt case we are repeating the Cbe calculation above.
!  To estimate particulate phosphate burial efficiency, we convert the
!  particulate phosphate fluxes to equivalent carbon fluxes by assuming
!  a fixed CP ratio, For now I use CNbur to do a test.
                fac1=12.0_r8/1000.0_r8*                                 &
     &               cff1*Hz(i,j,1)*(365.0_r8/dtdays) ! gC m-2 year-1
                IF (ibio.eq.iPhyt) THEN
                  fac1=fac1*PhyCN(ng)
                ELSE
                  fac1=fac1*CNbur*rNxP_P  
                ENDIF
                Cbe=MIN(0.75_r8,0.023_r8*fac1**0.5797_r8 )
!  Factor for converting phytoplankton nitrogen flux to phosphate
                IF (ibio.eq.iPhyt) THEN
                  fac2=NxP_P   
                ELSE
                  fac2=1.0_r8      
                ENDIF                
# ifdef DIAGNOSTICS_BIO
!  The flux of phosphate that reaches the seafloor
                DiaBio2d(i,j,iPbot)=DiaBio2d(i,j,iPbot)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              fac2*cff1*Hz(i,j,1)*fiter
!  The flux of phosphate that is buried
                DiaBio2d(i,j,iPbur)=DiaBio2d(i,j,iPbur)+                &
#  ifdef WET_DRY
     &                              rmask_full(i,j)*                    &
#  endif
     &                              Cbe*fac2*cff1*Hz(i,j,1)*fiter
# endif                
!  cff1 will now be the unburied flux that is to be remineralized
                cff1=(1.0_r8-Cbe)*cff1                        ! mmol m-3 
              END IF                

# ifdef CARBON
#  ifdef DENITRIFICATION

            cff3=12.0_r8
            cff4=0.74_r8
            cff5=16.0_r8/16.0_r8
#  endif

            IF ((ibio.eq.iSDeC).or.                                     &
     &          (ibio.eq.iLDeC))THEN
          !    DO i=Istr,Iend
          !      cff1=cff2*FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+cff1*(1.0_r8-0.01_r8)
          !    END DO
            END IF
            IF (ibio.eq.iPhyt)THEN
          !    DO i=Istr,Iend
          !      cff1=cff2*FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+                          &
     &                          cff1*PhyCN(ng)*(1.0_r8-0.01_r8)
          !    END DO
            END IF
            IF (ibio.eq.iDiaz)THEN
          !     DO i=Istr,Iend
          !      cff1=cff2*FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+                          &
     &                         cff1*PhyCN(ng)*(1.0_r8-0.01_r8)
          !     END DO
            END IF
# endif
# ifdef H_SULF

            IF ((ibio.eq.iSDeP).or.                                     &
     &          (ibio.eq.iLDeP))THEN
          !    DO i=Istr,Iend
          !      cff1=cff2*FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+                          &
     &                         cff1*(fox+fd1+fd2+fsr)
          !    END DO
            END IF
            IF (ibio.eq.iPhyt)THEN
          !   DO i=Istr,Iend
          !      cff1=cff2*FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+                          &
     &                         NxP_P*cff1*(fox+fd1+fd2+fsr)
          !   END DO
            END IF
            IF (ibio.eq.iDiaz)THEN
          !   DO i=Istr,Iend
          !      cff1=cff2*FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+                          &
     &                         cff1*(fox+fd1+fd2+fsr)/rNxP_D
          !   END DO
            END IF
# endif
          END DO
# endif

          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012)
!-----------------------------------------------------------------------
!
      
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      
      
      END DO J_LOOP
      RETURN
      END SUBROUTINE biology_tile

#ifdef CARBON
# ifdef pCO2_RZ
      SUBROUTINE pCO2_water_RZ (Istr, Iend,                             &
     &                          LBi, UBi, LBj, UBj, IminS, ImaxS,       &
     &                          j, DoNewton,                            &
#  ifdef MASKING
     &                          rmask,                                  &
#  endif
     &                          T, S, TIC, TAlk, pH, pCO2)       
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4=0,            !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2= ppmv  (DoNewton=0)                              !
!                pCO2= ppmv  (DoNewton=1)                              !
!                                                                      !
!  This subroutine was adapted by Katja Fennel (Nov 2005) from         !
!  Zeebe and Wolf-Gladrow (2001).                                      !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Zeebe, R.E. and D. Wolf-Gladrow,  2005:  CO2 in Seawater:         !
!      Equilibrium, kinetics, isotopes, Elsevier Oceanographic         !
!      Series, 65, pp 346.                                             !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30

      integer :: Hstep, Ibrack, Inewton, i

      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: scl, sqrtS
      real(r8) :: borate, alk, dic
      real(r8) :: ff, K1, K2, K12, Kb, Kw
      real(r8) :: p5, p4, p3, p2, p1, p0
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        scl=S(i)/1.80655_r8

        alk= TAlk(i)*0.000001_r8
        dic = TIC(i)*0.000001_r8
        
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carboinic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!-----------------------------------------------------------------------
!
        K1=10.0_r8**(62.008_r8-                                         &
     &               invTk*3670.7_r8-                                   &
     &               logTk*9.7944_r8+                                   &
     &               S(i)*(0.0118_r8-                                   &
     &                     S(i)*0.000116_r8))
        K2=10.0_r8**(-4.777_r8-                                         &
     &               invTk*1394.7_r8+                                   &
     &               S(i)*(0.0184_r8-                                   &
     &                     S(i)*0.000118_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Solve for h in fifth-order polynomial. First calculate
!  polynomial coefficients.
!
        K12 = K1*K2

        p5 = -1.0_r8;
        p4 = -alk-Kb-K1;
        p3 = dic*K1-alk*(Kb+K1)+Kb*borate+Kw-Kb*K1-K12
        p2 = dic*(Kb*K1+2*K12)-alk*(Kb*K1+K12)+Kb*borate*K1             &
     &       +(Kw*Kb+Kw*K1-Kb*K12)
        p1 = 2.0_r8*dic*Kb*K12-alk*Kb*K12+Kb*borate*K12                 &
     &       +Kw*Kb*K1+Kw*K12
        p0 = Kw*Kb*K12;
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN
          X=X_guess
!
          DO Inewton=1,InewtonMax
!
!  Evaluate f([H+]) = p5*x^5+...+p1*x+p0
!
            fn=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
!
!  Evaluate derivative, df([H+])/dx:
!
!     df= d(fn)/d(X)
!
            df=(((5*p5*X+4*p4)*X+3*p3)*X+2*p2)*X+p1
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0) THEN
               EXIT BRACK_IT
            ELSE
               ftest=fni(1)/fni(3)
               IF (ftest.gt.0) THEN
                 X_hi=X_mid
               ELSE
                 X_lo=X_mid
               END IF
               X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff

#  ifdef MASKING
      ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF
#  endif

      END DO I_LOOP

      RETURN
      END SUBROUTINE pCO2_water_RZ
# else
      SUBROUTINE pCO2_water (Istr, Iend,                                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, j, DoNewton,                 &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
     &                       T, S, TIC, TAlk, PO4, SiO3, pH, pCO2)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     PO4        Inorganic phosphate (millimol/m3).                    !
!     SiO3       Inorganic silicate (millimol/m3).                     !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4=0,            !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2=0.35074945E+03 ppmv  (DoNewton=0)                !
!                pCO2=0.35073560E+03 ppmv  (DoNewton=1)                !
!                                                                      !
!  This subroutine was adapted by Mick Follows (Oct 1999) from OCMIP2  !
!  code CO2CALC. Modified for ROMS by Hernan Arango (Nov 2003).        !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(in) :: PO4(IminS:)     
      real(r8), intent(inout) :: pH(LBi:,LBj:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(in) :: PO4(IminS:ImaxS)   
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(in) :: SiO3

      real(r8), intent(out) :: pCO2(IminS:ImaxS)

!
!  Local variable declarations.
!

      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30

      integer :: Hstep, Ibrack, Inewton, i

      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: SO4, scl, sqrtS, sqrtSO4
      real(r8) :: alk, dic, phos, sili
      real(r8) :: borate, sulfate, fluoride
      real(r8) :: ff, K1, K2, K1p, K2p, K3p, Kb, Kf, Ks, Ksi, Kw
      real(r8) :: K12, K12p, K123p, invKb, invKs, invKsi
      real(r8) :: A, A2, B, B2, C, dA, dB
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        SO4=19.924_r8*S(i)/(1000.0_r8-1.005_r8*S(i))
        sqrtSO4=SQRT(SO4)
        scl=S(i)/1.80655_r8

        alk=TAlk(i)*0.000001_r8
        dic=TIC(i)*0.000001_r8
        phos=PO4(i)*0.000001_r8
        sili=SiO3*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carboinic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!-----------------------------------------------------------------------
!
        K1=10.0_r8**(62.008_r8-                                         &
     &               invTk*3670.7_r8-                                   &
     &               logTk*9.7944_r8+                                   &
     &               S(i)*(0.0118_r8-                                   &
     &                     S(i)*0.000116_r8))
        K2=10.0_r8**(-4.777_r8-                                         &
     &               invTk*1394.7_r8+                                   &
     &               S(i)*(0.0184_r8-                                   &
     &                     S(i)*0.000118_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute first (K1p), second (K2p), and third (K3p) dissociation
!  constant of phosphoric acid:
!
!           K1p = [H][H2PO4]/[H3PO4]
!           K2p = [H][HPO4]/[H2PO4]
!           K3p = [H][PO4]/[HPO4]
!
!  From DOE (1994) equations 7.2.20, 7.2.23, and 7.2.26, respectively.
!  With footnote using data from Millero (1974).
!-----------------------------------------------------------------------
!
        K1p=EXP(115.525_r8-                                             &
     &          invTk*4576.752_r8-                                      &
     &          logTk*18.453_r8+                                        &
     &          sqrtS*(0.69171_r8-invTk*106.736_r8)-                    &
     &          S(i)*(0.01844_r8+invTk*0.65643_r8))
        K2p=EXP(172.0883_r8-                                            &
     &          invTk*8814.715_r8-                                      &
     &          logTk*27.927_r8+                                        &
     &          sqrtS*(1.3566_r8-invTk*160.340_r8)-                     &
     &          S(i)*(0.05778_r8-invTk*0.37335_r8))
        K3p=EXP(-18.141_r8-                                             &
     &          invTk*3070.75_r8+                                       &
     &          sqrtS*(2.81197_r8+invTk*17.27039_r8)-                   &
     &          S(i)*(0.09984_r8+invTk*44.99486_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of silica, Ksi=[H][SiO(OH)3]/[Si(OH)4].
!  From Millero (1995; page 671) using data from Yao and Millero (1995).
!-----------------------------------------------------------------------
!
        Ksi=EXP(117.385_r8-                                             &
     &          invTk*8904.2_r8-                                        &
     &          logTk*19.334_r8+                                        &
     &          sqrtSO4*(3.5913_r8-invTk*458.79_r8)-                    &
     &          SO4*(1.5998_r8-invTk*188.74_r8-                         &
     &               SO4*(0.07871_r8-invTk*12.1652_r8))+                &
     &          LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!------------------------------------------------------------------------
!  Compute salinity constant of hydrogen sulfate, Ks = [H][SO4]/[HSO4].
!  From Dickson (1990, J. chem. Thermodynamics 22, 113)
!------------------------------------------------------------------------
!
        Ks=EXP(141.328_r8-                                              &
     &         invTk*4276.1_r8-                                         &
     &         logTk*23.093_r8+                                         &
     &         sqrtSO4*(324.57_r8-invTk*13856.0_r8-logTk*47.986_r8-     &
     &                  SO4*invTk*2698.0_r8)-                           &
     &         SO4*(771.54_r8-invTk*35474.0_r8-logTk*114.723_r8-        &
     &              SO4*invTk*1776.0_r8)+                               &
     &         LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute stability constant of hydrogen fluorid, Kf = [H][F]/[HF].
!  From Dickson and Riley (1979) -- change pH scale to total.
!-----------------------------------------------------------------------
!
        Kf=EXP(-12.641_r8+                                              &
     &         invTk*1590.2_r8+                                         &
     &         sqrtSO4*1.525_r8+                                        &
     &         LOG(1.0_r8-0.001005_r8*S(i))+                            &
     &         LOG(1.0_r8+0.1400_r8*scl/(96.062_r8*Ks)))
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974), sulfate (Morris
! and Riley, 1966), and fluoride (Riley, 1965).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
        sulfate=0.14_r8*scl/96.062_r8
        fluoride=0.000067_r8*scl/18.9984_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN
          X=X_guess
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          DO Inewton=1,InewtonMax
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
            X2=X*X
            X3=X2*X
            invX=1.0_r8/X
            invX2=1.0_r8/X2

            A=X*(K12p+X*(K1p+X))
            B=X*(K1+X)+K12
            C=1.0_r8/(1.0_r8+sulfate*invKs)

            A2=A*A
            B2=B*B
            dA=X*(2.0_r8*K1p+3.0_r8*X)+K12p
            dB=2.0_r8*X+K1
!
!  Evaluate f([H+]):
!
!     fn=HCO3+CO3+borate+OH+HPO4+2*PO4+H3PO4+silicate+Hfree+HSO4+HF-TALK
!
            fn=dic*K1*(X+2.0_r8*K2)/B+                                  &
     &         borate/(1.0_r8+X*invKb)+                                 &
     &         Kw*invX+                                                 &
     &         phos*(K12p*X+2.0_r8*K123p-X3)/A+                         &
     &         sili/(1.0_r8+X*invKsi)-                                  &
     &         X*C-                                                     &
     &         sulfate/(1.0_r8+Ks*invX*C)-                              &
     &         fluoride/(1.0_r8+Kf*invX)-                               &
     &         alk
!
!  Evaluate derivative, f(prime)([H+]):
!
!     df= d(fn)/d(X)
!
            df=dic*K1*(B-dB*(X+2.0_r8*K2))/B2-                          &
     &         borate/(invKb*(1.0+X*invKb)**2)-                         &
     &         Kw*invX2+                                                &
     &         phos*(A*(K12p-3.0_r8*X2)-dA*(K12p*X+2.0_r8*K123p-X3))/A2-&
     &         sili/(invKsi*(1.0_r8+X*invKsi)**2)+                      &
     &         C+                                                       &
     &         sulfate*Ks*C*invX2/((1.0_r8+Ks*invX*C)**2)+              &
     &         fluoride*Kf*invX2/((1.0_r8+Kf*invX)**2)
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
              X2=X*X
              X3=X2*X
              invX=1.0_r8/X

              A=X*(K12p+X*(K1p+X))+K123p
              B=X*(K1+X)+K12
              C=1.0_r8/(1.0_r8+sulfate*invKs)

              A2=A*A
              B2=B*B
              dA=X*(K1p*2.0_r8+3.0_r8*X2)+K12p
              dB=2.0_r8*X+K1
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=dic*(K1*X+2.0_r8*K12)/B+                       &
     &                   borate/(1.0_r8+X*invKb)+                       &
     &                   Kw*invX+                                       &
     &                   phos*(K12p*X+2.0_r8*K123p-X3)/A+               &
     &                   sili/(1.0_r8+X*invKsi)-                        &
     &                   X*C-                                           &
     &                   sulfate/(1.0_r8+Ks*invX*C)-                    &
     &                   fluoride/(1.0_r8+Kf*invX)-                     &
     &                   alk
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0.0_r8) THEN
              EXIT BRACK_IT
            ELSE
              ftest=fni(1)/fni(3)
              IF (ftest.gt.0.0) THEN
                X_hi=X_mid
              ELSE
                X_lo=X_mid
              END IF
              X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff

#  ifdef MASKING
      ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF
#  endif

      END DO I_LOOP

      RETURN
      END SUBROUTINE pCO2_water
# endif
#endif
