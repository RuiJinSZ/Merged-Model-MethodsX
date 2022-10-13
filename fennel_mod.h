!
!svn $Id: redoxh_mod.h 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Fennel et al. (2006) model:                          !
!                                                                      !
!   AttSW    Light attenuation due to sea water [1/m].                 !
!   AttChl   Light attenuation by Chlorophyll [1/(mg_Chl m2)].         !
!   BioIter  Maximum number of iterations to achieve convergence       !
!              of the nonlinear solution.                              !
!   Chl2C_m  Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].        !
!   ChlMin   Chlorophill minimum threshold value [mg_Chl/m3].          !
!   CoagR    Coagulation rate: agregation rate of SDeN + Phyt ==> LDeN !
!              [1/day].                                                !
!   D_p5NH4  Half-saturation radiation for nitrification inhibition    !
!              [Watts/m2].                                             !
!   I_thNH4  Radiation threshold for nitrification inhibition          !
!              [Watts/m2].                                             !
!   K_NH4    Inverse half-saturation for Phytoplankton NH4 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_NO3    Inverse half-saturation for Phytoplankton NO3 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_Phy    Zooplankton half-saturation, squared constant for         !
!              ingestion [mmol_N/m3]^2.                                !
!   LDeRR    Large Detrital re-mineralization rate [1/day].            !
!   NitriR   Nitrification rate: oxidation of NH4 to NO3 [1/day].      !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   PhyCN    Phytoplankton Carbon:Nitrogen ratio [mol_C/mol_N].        !
!   PhyIP    Phytoplankton NH4 inhibition parameter [1/(mmol_N)].      !
!   PhyIS    Phytoplankton, initial slope of the P-I curve             !
!              [1/(W m-2 day)].                                        !
!   ZooMin   Phytoplankton minimum threshold value [mmol_N/m3].        !
!   PhyMR    Phytoplankton mortality rate [1/day] to small detritus.   !
!   SDeAR    Small detritus aggregation rate into Large detritus       !
!              [1/day].                                                !
!   SDeBR    Small Detrital breakdown to NH4 rate [1/day].             !
!   SDeRR    Large Detrital re-mineralization rate [1/day].            !
!   Vp0      Eppley temperature-limited and light-limited growth       !
!              tuning parameter [nondimensional].                      !
!   wLDet    Vertical sinking velocities for Large Detritus            !
!              fraction [m/day].                                       !
!   wPhy     Vertical sinking velocity for Phytoplankton               !
!              fraction [m/day].                                       !
!   wSDet    Vertical sinking velocities for Small Detritus            !
!              fraction [m/day].                                       !
!   ZooAE_N  Zooplankton nitrogen assimilation efficiency fraction     !
!              [nondimensional].                                       !
!   ZooBM    Zooplankton basal metabolism [1/day].                     !
!   ZooCN    Zooplankton Carbon:Nitrogen ratio [mol_C/mol_N].          !
!   ZooER    Zooplankton specific excretion rate [1/day].              !
!   ZooGR    Zooplankton maximum growth rate [1/day].                  !
!   ZooMin   Zooplankton minimum threshold value [mmol_N/m3].          !
!   ZooMR    Zooplankton mortality to Detritus [1/day].                !
!   pCO2air  CO2 partial pressure in the air [ppmv].                   !
!   O2palFrac  Fraction of modern atmospheric pO2 [nondim.]
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iNO3_                  ! Nitrate concentration
      integer :: iNH4_                  ! Ammonium concentration
      integer :: iPhyt                  ! Phytoplankton concentration
      integer :: iZoop                  ! Zooplankton concentration
      integer :: iLDeN                  ! Large detritus N-concentration
      integer :: iSDeN                  ! Small detritus N-concentration
      integer :: iChlo                  ! Chlorophyll concentration
      integer :: iN2__                  ! Dinitrogen concentration
      integer :: iDON_                  ! Semilabile dissolved organic N
      integer :: irfDON_                ! Refractory dissolved organic N
      integer :: iDOP_                  ! Semilabile dissolved organic P
#ifdef CARBON
      integer :: iTIC_                  ! Total inorganic carbon
      integer :: iTAlk                  ! Total alkalinity
      integer :: iLDeC                  ! Large detritus C-concentration
      integer :: iSDeC                  ! Small detritus C-concentration
#endif
#ifdef OXYGEN
      integer :: iOxyg                  ! Dissolved oxygen concentration
#endif
#ifdef H_SULF
      integer :: iH2S_                  ! Hydrogen sulfide concentration
      integer :: iSO4_                  ! Sulfate concentration
      integer :: iNO2_                  ! Nitrite concentration
      integer :: iPO4_                  ! Phosphate concentration
      integer :: iDiaz                  ! diazotroph concentration
      integer :: iLDeP                  ! Large detritus P-concentration
      integer :: iSDeP                  ! Small detritus P-concentration
#endif

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Biological 2D diagnostic variable IDs.
!
      integer, allocatable :: iDbio2(:)       ! 2D biological terms

      integer  :: iCOfx                       ! air-sea CO2 flux
      integer  :: iDNIT                       ! denitrification flux
      integer  :: ipCO2                       ! partial pressure of CO2
      integer  :: iO2fx                       ! air-sea O2 flux
      integer  :: iSRR_                       ! sulfate reduction rate
      integer  :: iBFlx                       ! bottom flux rate
      integer  :: iSRRb                       ! bottom H2S flux rate
      integer  :: iNbur                       ! particulate N flux that's buried
      integer  :: iNbot                       ! total N flux reaching seafloor
      integer  :: iCbur                       ! particulate C flux that's buried
      integer  :: iCbot                       ! total C flux reaching seafloor
      integer  :: iPbur                       ! particulate P flux that's buried
      integer  :: iPbot                       ! total P flux reaching seafloor

!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms

      integer  :: iPPro = 1                   ! primary productivity
      integer  :: iNO3u = 2                   ! NO3 uptake
      integer  :: iNFix = 3                   ! N Fixation
      integer  :: iSOx_ = 4                   ! S oxidation O2
      integer  :: iOxic = 5                   ! oxic remineralization
      integer  :: iSRRa = 6                   ! S reduction
      integer  :: iDno2 = 7                   ! Denitrification (NO2 reduction)
      integer  :: iAnx_ = 8                   ! Anammox
      integer  :: iNit1 = 9                   ! nitrification (NH4 -> NO2)
      integer  :: iNit2 = 10                  ! nitrification (NO2 -> NO3)
      integer  :: iDno3 = 11                  ! NO3 reduction to NO2
      integer  :: iSNx_ = 12                  ! S oxidation NO3
      integer  :: iSNy_ = 13                  ! S oxidation NO2
      integer  :: iPO4u = 14                  ! PO4 uptake
#endif
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: AttSW(:)              ! 1/m
      real(r8), allocatable :: AttChl(:)             ! 1/(mg_Chl m2)
      real(r8), allocatable :: Chl2C_m(:)            ! mg_Chl/mg_C
      real(r8), allocatable :: ChlMin(:)             ! mg_Chl/m3
      real(r8), allocatable :: CoagR(:)              ! 1/day
      real(r8), allocatable :: D_p5NH4(:)            ! Watts/m2
      real(r8), allocatable :: I_thNH4(:)            ! Watts/m2
      real(r8), allocatable :: K_NH4(:)              ! m3/mmol_N
      real(r8), allocatable :: K_NO3(:)              ! m3/mmol_N
      real(r8), allocatable :: K_Phy(:)              ! (mmol_N/m3)^2
      real(r8), allocatable :: LDeRRN(:)             ! 1/day
      real(r8), allocatable :: LDeNSR(:)             ! 1/day
      real(r8), allocatable :: LDeRRC(:)             ! 1/day
      real(r8), allocatable :: NitriR(:)             ! 1/day
      real(r8), allocatable :: PARfrac(:)            ! nondimensional
      real(r8), allocatable :: PhyCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: PhyIP(:)              ! 1/mmol_N
      real(r8), allocatable :: PhyIS(:)              ! 1/(Watts m-2 day)
      real(r8), allocatable :: PhyMin(:)             ! mmol_N/m3
      real(r8), allocatable :: PhyMR(:)              ! 1/day
      real(r8), allocatable :: SDeAR(:)              ! 1/day
      real(r8), allocatable :: SDeBR(:)              ! 1/day
      real(r8), allocatable :: SDeRRN(:)             ! 1/day
      real(r8), allocatable :: SDeNSR(:)             ! 1/day
      real(r8), allocatable :: SDeRRC(:)             ! 1/day
      real(r8), allocatable :: Vp0(:)                ! nondimensional
      real(r8), allocatable :: wLDet(:)              ! m/day
      real(r8), allocatable :: wPhy(:)               ! m/day
      real(r8), allocatable :: wSDet(:)              ! m/day
      real(r8), allocatable :: ZooAE_N(:)            ! nondimensional
      real(r8), allocatable :: ZooBM(:)              ! 1/day
      real(r8), allocatable :: ZooCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: ZooER(:)              ! 1/day
      real(r8), allocatable :: ZooGR(:)              ! 1/day
      real(r8), allocatable :: ZooMin(:)             ! mmol_N/m3
      real(r8), allocatable :: ZooMR(:)              ! 1/day
      real(r8), allocatable :: pCO2air(:)            ! ppmv
      real(r8), allocatable :: ElDON(:)              ! nondimensional
      real(r8), allocatable :: EsDON(:)              ! nondimensional
      real(r8), allocatable :: deltN(:)              ! nondimensional
      real(r8), allocatable :: a0N(:)                ! 1/day
      real(r8), allocatable :: O2palFrac(:)          ! nondimensional
      real(r8), allocatable :: rkd1(:)               ! m-1
      real(r8), allocatable :: rkdChl1(:)            ! m2 mg_Chl-1
      real(r8), allocatable :: rkdTSS1(:)            ! m2 g-1
      real(r8), allocatable :: rkdS1(:)              ! m-1 PSU-1

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
#ifdef CARBON
# ifdef OXYGEN
#  ifdef H_SULF
      NBT=23       ! carbon + oxygen + sulf
#  else
      NBT=13       ! carbon + oxygen
#  endif
# else
#  ifdef H_SULF
      NBT=19       ! carbon + sulf
#  else
      NBT=12       ! carbon
#  endif
# endif
#else
# ifdef OXYGEN
#  ifdef H_SULF
      NBT=16       ! oxygen + sulf
#  else
      NBT=9        ! oxygen
#  endif
# else
#  ifdef H_SULF
      NBT=15       ! sulf
#  else
      NBT=8        ! basic bio_redoxh
#  endif
# endif
#endif


#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio3d=14
      NDbio2d=0
# ifdef DENITRIFICATION
      NDbio2d=NDbio2d+1
# endif
# ifdef CARBON
      NDbio2d=NDbio2d+2
# endif
# ifdef OXYGEN
      NDbio2d=NDbio2d+1
# endif
# ifdef H_SULF
      NDbio2d=NDbio2d+3
# endif
# ifdef USECOS_BURIAL
      NDbio2d=NDbio2d+6
#endif      
!
!  Initialize biology diagnostic indices (2D-variables).
!
      ic=0
# ifdef DENITRIFICATION
      iDNIT=ic+1
      ic=ic+1
# endif
# ifdef CARBON
      iCOfx=ic+1
      ipCO2=ic+2
      ic=ic+2
# endif
# ifdef OXYGEN
      iO2fx=ic+1
      ic=ic+1
# endif
# ifdef H_SULF
      iSRR_=ic+1
      iBFlx=ic+2
      iSRRb=ic+3
      ic=ic+3
# endif
# ifdef USECOS_BURIAL
      iNbur = ic+1
      iCbot = ic+2
      iCbur = ic+3
      iNbot = ic+4
      iPbot = ic+5
      iPbur = ic+6
      ic = ic+6
# endif
# endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(AttChl)) THEN
        allocate ( AttChl(Ngrids) )
      END IF
      IF (.not.allocated(Chl2C_m)) THEN
        allocate ( Chl2C_m(Ngrids) )
      END IF
      IF (.not.allocated(ChlMin)) THEN
        allocate ( ChlMin(Ngrids) )
      END IF
      IF (.not.allocated(CoagR)) THEN
        allocate ( CoagR(Ngrids) )
      END IF
      IF (.not.allocated(D_p5NH4)) THEN
        allocate ( D_p5NH4(Ngrids) )
      END IF
      IF (.not.allocated(I_thNH4)) THEN
        allocate ( I_thNH4(Ngrids) )
      END IF
      IF (.not.allocated(K_NH4)) THEN
        allocate ( K_NH4(Ngrids) )
      END IF
      IF (.not.allocated(K_NO3)) THEN
        allocate ( K_NO3(Ngrids) )
      END IF
      IF (.not.allocated(K_Phy)) THEN
        allocate ( K_Phy(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRN)) THEN
        allocate ( LDeRRN(Ngrids) )
      END IF
      IF (.not.allocated(LDeNSR)) THEN
        allocate ( LDeNSR(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRC)) THEN
        allocate ( LDeRRC(Ngrids) )
      END IF
      IF (.not.allocated(NitriR)) THEN
        allocate ( NitriR(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(PhyCN)) THEN
        allocate ( PhyCN(Ngrids) )
      END IF
      IF (.not.allocated(PhyIP)) THEN
        allocate ( PhyIP(Ngrids) )
      END IF
      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
      END IF
      IF (.not.allocated(PhyMin)) THEN
        allocate ( PhyMin(Ngrids) )
      END IF
      IF (.not.allocated(PhyMR)) THEN
        allocate ( PhyMR(Ngrids) )
      END IF
      IF (.not.allocated(SDeAR)) THEN
        allocate ( SDeAR(Ngrids) )
      END IF
      IF (.not.allocated(SDeBR)) THEN
        allocate ( SDeBR(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRN)) THEN
        allocate ( SDeRRN(Ngrids) )
      END IF
      IF (.not.allocated(SDeNSR)) THEN
        allocate ( SDeNSR(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRC)) THEN
        allocate ( SDeRRC(Ngrids) )
      END IF
      IF (.not.allocated(Vp0)) THEN
        allocate ( Vp0(Ngrids) )
      END IF
      IF (.not.allocated(wLDet)) THEN
        allocate ( wLDet(Ngrids) )
      END IF
      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
      END IF
      IF (.not.allocated(wSDet)) THEN
        allocate ( wSDet(Ngrids) )
      END IF
      IF (.not.allocated(ZooAE_N)) THEN
        allocate ( ZooAE_N(Ngrids) )
      END IF
      IF (.not.allocated(ZooBM)) THEN
        allocate ( ZooBM(Ngrids) )
      END IF
      IF (.not.allocated(ZooCN)) THEN
        allocate ( ZooCN(Ngrids) )
      END IF
      IF (.not.allocated(ZooER)) THEN
        allocate ( ZooER(Ngrids) )
      END IF
      IF (.not.allocated(ZooGR)) THEN
        allocate ( ZooGR(Ngrids) )
      END IF
      IF (.not.allocated(ZooMin)) THEN
        allocate ( ZooMin(Ngrids) )
      END IF
      IF (.not.allocated(ZooMR)) THEN
        allocate ( ZooMR(Ngrids) )
      END IF
      IF (.not.allocated(pCO2air)) THEN
        allocate ( pCO2air(Ngrids) )
      END IF
      IF (.not.allocated(ElDON)) THEN
        allocate ( ElDON(Ngrids) )
      END IF
      IF (.not.allocated(EsDON)) THEN
        allocate ( EsDON(Ngrids) )
      END IF
      IF (.not.allocated(deltN)) THEN
        allocate ( deltN(Ngrids) )
      END IF
      IF (.not.allocated(a0N)) THEN
        allocate ( a0N(Ngrids) )
      END IF  
      IF (.not.allocated(O2palFrac)) THEN
        allocate ( O2palFrac(Ngrids) )
      END IF
      IF (.not.allocated(rkd1)) THEN
        allocate ( rkd1(Ngrids) )
      END IF
      IF (.not.allocated(rkdChl1)) THEN
        allocate ( rkdChl1(Ngrids) )
      END IF  
      IF (.not.allocated(rkdTSS1)) THEN
        allocate ( rkdTSS1(Ngrids) )
      END IF
      IF (.not.allocated(rkdS1)) THEN
        allocate ( rkdS1(Ngrids) )
      END IF
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iNO3_=ic+1
      iNH4_=ic+2
      iPhyt=ic+3
      iZoop=ic+4
      iLDeN=ic+5
      iSDeN=ic+6
      iChlo=ic+7
      iN2__=ic+8
      iDON_=ic+9
      irfDON_=ic+10
      iDOP_=ic+11
      ic=ic+11
# ifdef CARBON
      iTIC_=ic+1
      iTAlk=ic+2
      iLDeC=ic+3
      iSDeC=ic+4
      ic=ic+4
# endif
# ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
# endif
# ifdef H_SULF
      iH2S_=ic+1
      iSO4_=ic+2
      iNO2_=ic+3
      iPO4_=ic+4
      iDiaz=ic+5
      iLDeP=ic+6
      iSDeP=ic+7
      ic=ic+7
# endif
      RETURN
      END SUBROUTINE initialize_biology
