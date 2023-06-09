module CNGrazingMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Routines different grazing methods
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_infnan_mod          , only : nan => shr_infnan_nan, assignment(=)
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use perf_mod                , only : t_startf, t_stopf
  use clm_varctl              , only : iulog
  use clm_time_manager        , only : get_curr_date
  use clm_varcon              , only : c_to_b, secspday
  use GridcellType            , only : grc
  use LandunitType            , only : lun
  use ColumnType              , only : col
  use PatchType               , only : patch
  use pftconMod               , only : pftcon
  use CNVegCarbonStateType    , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType     , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType  , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType   , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType , only : soilbiogeochem_state_type
  use TemperatureType         , only : temperature_type
  use CNGrazerType            , only : grazer_type
  use CNSharedParamsMod       , only : use_matrixcn
  use SurfaceAlbedoType       , only : surfalb_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNGrazing 
  public :: CNGrazingSeligman92
  public :: CNGrazingPftToColumn
  !
  !
  ! The following is only public for the sake of unit testing; it should not be called
  ! directly by CLM code outside this module

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  


  real(r8), private, parameter :: test_param = 10.0_r8   ! testing parameter
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  
  subroutine CNGrazing(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, grazer_inst,                 &
       temperature_inst, surfalb_inst)

  ! USES:
  use pftconMod       , only : noveg, nbrdlf_evr_shrub
  ! ARGUMENTS:

    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(grazer_type)               , intent(inout) :: grazer_inst
    type(temperature_type)          , intent(in)    :: temperature_inst
    type(surfalb_type)              , intent(in)    :: surfalb_inst




  ! LOCAL VARIABLES:

    call t_startf( 'CNGrazingSeligman92' )

      call CNGrazingSeligman92(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
         cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, grazer_inst,                 &
         temperature_inst, surfalb_inst)

    call t_stopf( 'CNGrazingSeligman92' )

  end subroutine CNGrazing
  
  
  

  !-----------------------------------------------------------------------
  subroutine CNGrazingSeligman92 (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, grazer_inst,                 &
       temperature_inst, surfalb_inst)
    !
    ! !DESCRIPTION:
    ! Grazing mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use pftconMod
    use clm_varcon      , only : secspday, c_to_b
    use clm_time_manager, only : get_step_size_real, is_near_local_noon
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(grazer_type)               , intent(inout) :: grazer_inst
    type(temperature_type)          , intent(in)    :: temperature_inst
    type(surfalb_type)              , intent(in)    :: surfalb_inst



    !
    ! !LOCAL VARIABLES:
    integer :: yr,mon,day,mcsec        ! outputs from get_curr_date
    integer :: p                         ! patch index
    integer :: g                         ! gridcell index
    integer :: c                         ! column index
    integer :: fp                        ! patch filter index
    real(r8) :: grassbioc                 ! carbon in grass for calculating grazed fraction (gC/m2)
    real(r8) :: grassbion                 ! nitrogen in grass for calculating grazed fraction (gN/m2)
    real(r8) :: mc                        ! rate for fractional grazing mortality carbon (1/s)
    real(r8) :: mn                        ! rate for fractional grazing mortality nitrogen (1/s)
    real(r8) :: vr                        ! residual carbon unavailible for grazing (gC/m2)
    real(r8) :: grasscn                   ! c:n ratio for grass for calculating grazed fraction
    real(r8) :: cx                        ! satiation consumption rate (gC/m2)
    real(r8) :: herbdens                  ! herbivore density   (animal/m2)
    real(r8) :: s                         ! grazing efficiency  (m2/animal)
    real(r8) :: ch4_frac                  ! ch4 output by grazers constant in CO2 equivalent
    real(r8) :: crherbfrac                 ! fraction of herbivore respiration
    real(r8) :: cfaecesfrac               ! fraction of feaces to litter
    real(r8) :: cfrac2litr(bounds%begp:bounds%endp)                ! fraction of grazed biomass going to litter C aboveground
    real(r8) :: nfrac2litr(bounds%begp:bounds%endp)                ! fraction of grazed biomass going to litter N aboveground
    real(r8) :: curinefrac                ! carbon fraction from urine to soil
    real(r8) :: dtime                     ! model time step (s)
    real(r8) :: coszen                    ! coszen to determine grazing time
    !-----------------------------------------------------------------------

    associate(& 
         ivt                                 =>    patch%itype                                                    , & ! Input:  [integer (:)]  pft vegetation type                                
         is_active                           =>    patch%active                                                   , &

         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                             , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                            , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
         livestemc                           =>    cnveg_carbonstate_inst%livestemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
         deadstemc                           =>    cnveg_carbonstate_inst%deadstemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
         livecrootc                          =>    cnveg_carbonstate_inst%livecrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
         deadcrootc                          =>    cnveg_carbonstate_inst%deadcrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
         gresp_xfer                          =>    cnveg_carbonstate_inst%gresp_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
         
         leafn                               =>    cnveg_nitrogenstate_inst%leafn_patch                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
         frootn                              =>    cnveg_nitrogenstate_inst%frootn_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
         livestemn                           =>    cnveg_nitrogenstate_inst%livestemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
         deadstemn                           =>    cnveg_nitrogenstate_inst%deadstemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
         livecrootn                          =>    cnveg_nitrogenstate_inst%livecrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
         deadcrootn                          =>    cnveg_nitrogenstate_inst%deadcrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
         retransn                            =>    cnveg_nitrogenstate_inst%retransn_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               

        coszen                               =>    surfalb_inst%coszen_col                                        , & ! Input:  [real(r8) (:)] column cosine of solar zenith angle
        t_veg10_day                         =>    temperature_inst%t_veg10_day_patch                             , & ! Input:  [real(r8) (:)] 10 day running mean of patch daytime time vegetation temperature (Kelvin), LUNA specific, but can be reused
        leafcn                               =>    pftcon%leafcn                                                  , & ! Input:  leaf C:N (gC/gN)


         grz_leafc_to_litter                 =>    cnveg_carbonflux_inst%grz_leafc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
         grz_frootc_to_litter                =>    cnveg_carbonflux_inst%grz_frootc_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
         grz_livestemc_to_litter             =>    cnveg_carbonflux_inst%grz_livestemc_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
         grz_deadstemc_to_litter             =>    cnveg_carbonflux_inst%grz_deadstemc_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
         grz_livecrootc_to_litter            =>    cnveg_carbonflux_inst%grz_livecrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         grz_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%grz_deadcrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         grz_to_atm                          =>    cnveg_carbonflux_inst%grz_to_atm_patch                         , & ! Output: [real(r8) (:)]                                                    
         grz_leafc_storage_to_litter         =>    cnveg_carbonflux_inst%grz_leafc_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         grz_frootc_storage_to_litter        =>    cnveg_carbonflux_inst%grz_frootc_storage_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         grz_livestemc_storage_to_litter     =>    cnveg_carbonflux_inst%grz_livestemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         grz_deadstemc_storage_to_litter     =>    cnveg_carbonflux_inst%grz_deadstemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         grz_livecrootc_storage_to_litter    =>    cnveg_carbonflux_inst%grz_livecrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         grz_deadcrootc_storage_to_litter    =>    cnveg_carbonflux_inst%grz_deadcrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         grz_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%grz_gresp_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         grz_leafc_xfer_to_litter            =>    cnveg_carbonflux_inst%grz_leafc_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         grz_frootc_xfer_to_litter           =>    cnveg_carbonflux_inst%grz_frootc_xfer_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         grz_livestemc_xfer_to_litter        =>    cnveg_carbonflux_inst%grz_livestemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         grz_deadstemc_xfer_to_litter        =>    cnveg_carbonflux_inst%grz_deadstemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         grz_livecrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%grz_livecrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         grz_deadcrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%grz_deadcrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         grz_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%grz_gresp_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         
         grz_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%grz_leafn_to_litter_patch              , & ! Output: [real(r8) (:)]                                                    
         grz_frootn_to_litter                =>    cnveg_nitrogenflux_inst%grz_frootn_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
         grz_livestemn_to_litter             =>    cnveg_nitrogenflux_inst%grz_livestemn_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         grz_deadstemn_to_litter             =>    cnveg_nitrogenflux_inst%grz_deadstemn_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         grz_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%grz_livecrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         grz_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%grz_deadcrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         grz_retransn_to_litter              =>    cnveg_nitrogenflux_inst%grz_retransn_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         grz_leafn_storage_to_litter         =>    cnveg_nitrogenflux_inst%grz_leafn_storage_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         grz_frootn_storage_to_litter        =>    cnveg_nitrogenflux_inst%grz_frootn_storage_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         grz_livestemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%grz_livestemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         grz_deadstemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%grz_deadstemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         grz_livecrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%grz_livecrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         grz_deadcrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%grz_deadcrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         grz_leafn_xfer_to_litter            =>    cnveg_nitrogenflux_inst%grz_leafn_xfer_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         grz_frootn_xfer_to_litter           =>    cnveg_nitrogenflux_inst%grz_frootn_xfer_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         grz_livestemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%grz_livestemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         grz_deadstemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%grz_deadstemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         grz_livecrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%grz_livecrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         grz_deadcrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%grz_deadcrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    

         grazed_totc                         =>    grazer_inst%grazed_totc_patch                                  , & ! Output: [real(r8) (:)]
         grazed_closs                        =>    grazer_inst%grazed_closs_patch                                 , & ! Output: [real(r8) (:)]
         grazed_totn                         =>    grazer_inst%grazed_totn_patch                                  , & ! Output: [real(r8) (:)]
         grazed_nloss                        =>    grazer_inst%grazed_nloss_patch                                 , & ! Output: [real(r8) (:)]
         grazed_pfts                         =>    grazer_inst%grazed_pfts                                          &  !Input:  [integer  (:)]
         
                                                    
         )



      dtime = get_step_size_real()
      ! valuse for Luo(2012)
      vr = 150.0_r8 * 1000 / 10000.0_r8 / c_to_b  ! convert from kg/ha biomass to gC/m2/sheep
      cx = 10.4_r8 / c_to_b * 1000.0_r8 / secspday ! convert from kg/day biomass to gC/sec per animal
      s  = 0.067 * 10000.0_r8 / secspday  ! convert from ha/day to m2/sec 
      herbdens  = 10.0_r8 / 10000.0_r8           ! convert from animal/ha to animal/m2
      ch4_frac = 0.03_r8 ! ch4 output by grazers constant in C
      crherbfrac = 0.5_r8       ! fraction of herbivore respiration
      cfaecesfrac = 0.3_r8         ! fraction of feaces to litter
      curinefrac = 0.8_r8 * 0.6_r8 * (12._r8 / 28._r8)           ! carbon fraction from urine to soil
      cfrac2litr(bounds%begp:bounds%endp) = 0.0_r8
      nfrac2litr(bounds%begp:bounds%endp) = 0.0_r8
      do fp = 1,num_soilp
        p = filter_soilp(fp)
        g = patch%gridcell(p)
        c = patch%column(p)
        grassbioc = 0.0_r8
        grassbion = 0.0_r8
        grasscn   = 1.0_r8
        mc = 0.0_r8
        mn = 0.0_r8
        grazed_totc(p) = 0.0_r8
        grazed_totn(p) = 0.0_r8
        grazed_closs(p) = 0.0_r8
        grazed_nloss(p) = 0.0_r8
        if (any(grazed_pfts  == ivt(p))) then
          write( iulog , * ) 'grz_check1'
          if ( coszen(c) > 0.0_r8 .and. t_veg10_day(p) >= 277.15_r8 ) then
            !get aboveground grazable biomass
            write( iulog , * ) 'grz_check2'
            call get_curr_date(yr,mon,day,mcsec)
            write( iulog , * ) 'grzdate ', yr,' ', mon,' ', day,' ', mcsec/3600 
              grassbioc =  leafc(p)              + &
                           leafc_storage(p)      + &
                           leafc_xfer(p)         + &
                           livestemc(p)          + &
                           livestemc_storage(p)  + &
                           livestemc_xfer(p)     + &
                           deadstemc(p)          + &
                           deadstemc_storage(p)  + &
                           deadstemc_xfer(p)     + &
                           gresp_storage(p)      + &
                           gresp_xfer(p)         
              grassbion =  leafn(p)              + &
                           leafn_storage(p)      + &
                           leafn_xfer(p)         + &
                           livestemn(p)          + &
                           livestemn_storage(p)  + &
                           livestemn_xfer(p)     + &
                           deadstemn(p)          + &
                           deadstemn_storage(p)  + &
                           deadstemn_xfer(p)     + &
                           retransn(p)
              write( iulog , * ) 'grass_c:', grassbioc , 'gC/m2 coszen',coszen(c), 'temp', t_veg10_day(p) - 273.15
              write( iulog , * ) 'dgrz_c:', s * herbdens * (grassbioc - vr) , 'gC/m2/s ' , herbdens * cx
              if ( s * herbdens * (grassbioc - vr) > 0.0_r8 ) then
                grazed_totc(p) = min(s * herbdens * (grassbioc - vr), herbdens * cx)
                mc = grazed_totc(p) / grassbioc
                if (grassbion <= 0.0_r8) then
                  grasscn = leafcn(ivt(p))
                  grazed_totn(p) = 0.0_r8
                  mn = 0.0_r8
                else
                  grasscn = grassbioc / grassbion
                  grazed_totn(p) = grazed_totc(p) / grasscn
                  mn = grazed_totn(p) / grassbion 
                endif
                write( iulog , * ) 'grass_n:', grassbion , 'gC/m2 grasscn', grasscn 
                write( iulog , * ) 'dgrz_n:', grazed_totn(p)
                ! just to make sure there are no crap values
                mc = max(0.0_r8 , mc )
                mn = max(0.0_r8 , mc )
                grazed_totc(p) = max(0.0_r8 , grazed_totc(p))
                grazed_totn(p) = max(0.0_r8 , grazed_totn(p))
              end if
            end if
        end if
        write( iulog , * ) 'grazed_c:', grazed_totc(p) , 'gC/m2'
        write( iulog , * ) 'grazed_n:', grazed_totn(p) , 'gC/m2'
        ! set litter fractions
        cfrac2litr(p) = min(cfaecesfrac+curinefrac , 1.0_r8)
        nfrac2litr(p) = cfrac2litr(p) / grasscn
        ! patch-level grazing carbon fluxes
        ! displayed pools
        if(.not. use_matrixcn)then
          grz_leafc_to_litter(p)               = leafc(p)               * mc 
          grz_frootc_to_litter(p)              = frootc(p)              * mc
          grz_livestemc_to_litter(p)           = livestemc(p)           * mc 
          grz_deadstemc_to_litter(p)           = deadstemc(p)           * mc 
          grz_livecrootc_to_litter(p)          = livecrootc(p)          * mc
          grz_deadcrootc_to_litter(p)          = deadcrootc(p)          * mc
          grz_to_atm(p)                        = grazed_totc(p) * ch4_frac

          ! storage pools
          grz_leafc_storage_to_litter(p)       = leafc_storage(p)       * mc 
          grz_frootc_storage_to_litter(p)      = frootc_storage(p)      * mc
          grz_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * mc 
          grz_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * mc 
          grz_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * mc
          grz_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * mc
          grz_gresp_storage_to_litter(p)       = gresp_storage(p)       * mc 

          ! transfer pools
          grz_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * mc 
          grz_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * mc
          grz_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * mc
          grz_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * mc 
          grz_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * mc
          grz_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * mc
          grz_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * mc 

          ! calculate losses carbon (what is not going to the litter)
          grazed_closs(p) = mc * ((leafc(p) + livestemc(p) + deadstemc(p) + &
          leafc_storage(p) + livestemc_storage(p) + deadstemc_storage(p) + gresp_storage(p) + &
          leafc_xfer(p) + livestemc_xfer(p) + deadstemc_xfer(p) + gresp_xfer(p))*(1.0_r8 -cfrac2litr(p)))

          ! patch-level grazing mortality nitrogen fluxes
          ! displayed pools
          grz_leafn_to_litter(p)               = leafn(p)               * mn
          grz_frootn_to_litter(p)              = frootn(p)              * mn
          grz_livestemn_to_litter(p)           = livestemn(p)           * mn
          grz_deadstemn_to_litter(p)           = livestemn(p)           * mn
          grz_livecrootn_to_litter(p)          = livecrootn(p)          * mn
          grz_deadcrootn_to_litter(p)          = deadcrootn(p)          * mn
          grz_retransn_to_litter(p)            = retransn(p)            * mn

          ! storage pools
          grz_leafn_storage_to_litter(p)       = leafn_storage(p)       * mn
          grz_frootn_storage_to_litter(p)      = frootn_storage(p)      * mn
          grz_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * mn
          grz_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * mn
          grz_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * mn
          grz_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * mn

          ! transfer pools
          grz_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * mn
          grz_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * mn
          grz_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * mn
          grz_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * mn
          grz_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * mn
          grz_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * mn

          ! calculate losses nitrogen  (what is not going to the litter)
          grazed_nloss(p) = mc * ((leafn(p) + livestemn(p) + deadstemn(p) + retransn(p) + &
          leafn_storage(p) + livestemn_storage(p) + deadstemc_storage(p) + &
          leafn_xfer(p) + livestemn_xfer(p) + deadstemn_xfer(p)) * (1.0_r8 -nfrac2litr(p)))

        ! NOTE: The non-matrix part of this update is in CNCStatUpdate2 CStateUpdate2grz (EBK 11/25/2019)
        !   and for Nitrogen The non-matrix part of this update is in CNNStatUpdate2 NStateUpdate2grz (EBK 11/25/2019)
        else ! use_matrixcn
        end if
      end do ! end of pft loop

      ! gather all patch-level litterfall fluxes from grazing to the column
      ! for litter C and N inputs

      call CNGrazingPftToColumn(num_soilc, filter_soilc, &
           soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
           grazer_inst, cfrac2litr, nfrac2litr)

    end associate 

  end subroutine CNGrazingSeligman92

 !-----------------------------------------------------------------------
 subroutine CNGrazingPftToColumn (num_soilc, filter_soilc, &
      soilbiogeochem_state_inst, CNVeg_carbonflux_inst, cnveg_nitrogenflux_inst, &
      grazer_inst, cfrac2litr, nfrac2litr)
   !
   ! !DESCRIPTION:
   ! called at the end of CNGrazing to gather all patch-level grazign litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clm_varpar , only : nlevdecomp, maxsoil_patches, i_litr_min, i_litr_max, i_met_lit
   !
   ! !ARGUMENTS:
   integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                         , intent(in)    :: filter_soilc(:) ! soil column filter
   type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
   type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
   type(grazer_type)               , intent(inout) :: grazer_inst
   real(r8)                        , intent(in)    :: cfrac2litr(:)
   real(r8)                        , intent(in)    :: nfrac2litr(:)
   !
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p,j,i  ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                                   & 
        ivt                              =>    patch%itype                                                      , & ! Input:  [integer  (:)   ]  pft vegetation type                                
        wtcol                            =>    patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
        
        lf_f                             =>    pftcon%lf_f                                                    , & ! Input:  leaf litter fraction
        fr_f                             =>    pftcon%fr_f                                                    , & ! Input:  fine root litter fraction
        
        leaf_prof                        =>    soilbiogeochem_state_inst%leaf_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    soilbiogeochem_state_inst%froot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    soilbiogeochem_state_inst%croot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    soilbiogeochem_state_inst%stem_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        grz_leafc_to_litter                 =>    cnveg_carbonflux_inst%grz_leafc_to_litter_patch                , & ! Input: [real(r8) (:)]                                                    
        grz_frootc_to_litter                =>    cnveg_carbonflux_inst%grz_frootc_to_litter_patch               , & ! Input: [real(r8) (:)]                                                    
        grz_livestemc_to_litter             =>    cnveg_carbonflux_inst%grz_livestemc_to_litter_patch            , & ! Input: [real(r8) (:)]                                                    
        grz_deadstemc_to_litter             =>    cnveg_carbonflux_inst%grz_deadstemc_to_litter_patch            , & ! Input: [real(r8) (:)]                                                    
        grz_livecrootc_to_litter            =>    cnveg_carbonflux_inst%grz_livecrootc_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        grz_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%grz_deadcrootc_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        grz_to_atm                          =>    cnveg_carbonflux_inst%grz_to_atm_patch                         , & ! Input: [real(r8) (:)]                                                    
        grz_leafc_storage_to_litter         =>    cnveg_carbonflux_inst%grz_leafc_storage_to_litter_patch        , & ! Input: [real(r8) (:)]                                                    
        grz_frootc_storage_to_litter        =>    cnveg_carbonflux_inst%grz_frootc_storage_to_litter_patch       , & ! Input: [real(r8) (:)]                                                    
        grz_livestemc_storage_to_litter     =>    cnveg_carbonflux_inst%grz_livestemc_storage_to_litter_patch    , & ! Input: [real(r8) (:)]                                                    
        grz_deadstemc_storage_to_litter     =>    cnveg_carbonflux_inst%grz_deadstemc_storage_to_litter_patch    , & ! Input: [real(r8) (:)]                                                    
        grz_livecrootc_storage_to_litter    =>    cnveg_carbonflux_inst%grz_livecrootc_storage_to_litter_patch   , & ! Input: [real(r8) (:)]                                                    
        grz_deadcrootc_storage_to_litter    =>    cnveg_carbonflux_inst%grz_deadcrootc_storage_to_litter_patch   , & ! Input: [real(r8) (:)]                                                    
        grz_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%grz_gresp_storage_to_litter_patch        , & ! Input: [real(r8) (:)]                                                    
        grz_leafc_xfer_to_litter            =>    cnveg_carbonflux_inst%grz_leafc_xfer_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        grz_frootc_xfer_to_litter           =>    cnveg_carbonflux_inst%grz_frootc_xfer_to_litter_patch          , & ! Input: [real(r8) (:)]                                                    
        grz_livestemc_xfer_to_litter        =>    cnveg_carbonflux_inst%grz_livestemc_xfer_to_litter_patch       , & ! Input: [real(r8) (:)]                                                    
        grz_deadstemc_xfer_to_litter        =>    cnveg_carbonflux_inst%grz_deadstemc_xfer_to_litter_patch       , & ! Input: [real(r8) (:)]                                                    
        grz_livecrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%grz_livecrootc_xfer_to_litter_patch      , & ! Input: [real(r8) (:)]                                                    
        grz_deadcrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%grz_deadcrootc_xfer_to_litter_patch      , & ! Input: [real(r8) (:)]                                                    
        grz_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%grz_gresp_xfer_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    

        grazed_totcp                        =>    grazer_inst%grazed_totc_patch                                  , & ! Input: [real(r8) (:)]
        grazed_clossp                       =>    grazer_inst%grazed_closs_patch                                 , & ! Input: [real(r8) (:)]

        grazed_totcc                        =>    grazer_inst%grazed_totc_col                                    , & ! Inout: [real(r8) (:)]
        grazed_clossc                       =>    grazer_inst%grazed_closs_col                                   , & ! Inout: [real(r8) (:)]
        grazing_c_to_litr_c                 =>    cnveg_carbonflux_inst%grazing_c_to_litr_c_col                  , & ! InOut: [real(r8) (:,:,:) ]  C fluxes associated with grazing to litter pools (gC/m3/s)
        grazing_c_to_cwdc                   =>    cnveg_carbonflux_inst%grazing_c_to_cwdc_col                    , & ! InOut: [real(r8) (:,:,:) ]  C fluxes associated with grazing to cwd pools (gC/m3/s)

        grz_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%grz_leafn_to_litter_patch              , & ! Input: [real(r8) (:)]                                                    
        grz_frootn_to_litter                =>    cnveg_nitrogenflux_inst%grz_frootn_to_litter_patch             , & ! Input: [real(r8) (:)]                                                    
        grz_livestemn_to_litter             =>    cnveg_nitrogenflux_inst%grz_livestemn_to_litter_patch          , & ! Input: [real(r8) (:)]                                                    
        grz_deadstemn_to_litter             =>    cnveg_nitrogenflux_inst%grz_livestemn_to_litter_patch          , & ! Input: [real(r8) (:)]                                                    
        grz_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%grz_livecrootn_to_litter_patch         , & ! Input: [real(r8) (:)]                                                    
        grz_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%grz_deadcrootn_to_litter_patch         , & ! Input: [real(r8) (:)]                                                    
        grz_retransn_to_litter              =>    cnveg_nitrogenflux_inst%grz_retransn_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        grz_leafn_storage_to_litter         =>    cnveg_nitrogenflux_inst%grz_leafn_storage_to_litter_patch      , & ! Input: [real(r8) (:)]                                                    
        grz_frootn_storage_to_litter        =>    cnveg_nitrogenflux_inst%grz_frootn_storage_to_litter_patch     , & ! Input: [real(r8) (:)]                                                    
        grz_livestemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%grz_livestemn_storage_to_litter_patch  , & ! Input: [real(r8) (:)]                                                    
        grz_deadstemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%grz_deadstemn_storage_to_litter_patch  , & ! Input: [real(r8) (:)]                                                    
        grz_livecrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%grz_livecrootn_storage_to_litter_patch , & ! Input: [real(r8) (:)]                                                    
        grz_deadcrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%grz_deadcrootn_storage_to_litter_patch , & ! Input: [real(r8) (:)]                                                    
        grz_leafn_xfer_to_litter            =>    cnveg_nitrogenflux_inst%grz_leafn_xfer_to_litter_patch         , & ! Input: [real(r8) (:)]                                                    
        grz_frootn_xfer_to_litter           =>    cnveg_nitrogenflux_inst%grz_frootn_xfer_to_litter_patch        , & ! Input: [real(r8) (:)]                                                    
        grz_livestemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%grz_livestemn_xfer_to_litter_patch     , & ! Input: [real(r8) (:)]                                                    
        grz_deadstemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%grz_deadstemn_xfer_to_litter_patch     , & ! Input: [real(r8) (:)]                                                    
        grz_livecrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%grz_livecrootn_xfer_to_litter_patch    , & ! Input: [real(r8) (:)]                                                    
        grz_deadcrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%grz_deadcrootn_xfer_to_litter_patch    , & ! Input: [real(r8) (:)]                                                    

        grazed_totnp                        =>    grazer_inst%grazed_totn_patch                                  , & ! Input: [real(r8) (:)]
        grazed_nlossp                       =>    grazer_inst%grazed_nloss_patch                                 , & ! Input: [real(r8) (:)]

        grazed_totnc                        =>    grazer_inst%grazed_totn_col                                    , & ! Input: [real(r8) (:)]
        grazed_nlossc                       =>    grazer_inst%grazed_nloss_col                                   , & ! Input: [real(r8) (:)]
        grazing_n_to_litr_n                 =>    cnveg_nitrogenflux_inst%grazing_n_to_litr_n_col                , & ! InOut: [real(r8) (:,:,:) ]  N fluxes associated with grazing to litter pools (gC/m3/s)
        grazing_n_to_cwdn                   =>    cnveg_nitrogenflux_inst%grazing_n_to_cwdn_col                    & ! InOut: [real(r8) (:,:,:) ]  N fluxes associated with grazing to cwd pools (gC/m3/s)
        )

     do j = 1, nlevdecomp
        do pi = 1,maxsoil_patches
           do fc = 1,num_soilc
              c = filter_soilc(fc)

              if (pi <=  col%npatches(c)) then
                 p = col%patchi(c) + pi - 1

                 if (patch%active(p)) then

                    do i = i_litr_min, i_litr_max
                       ! leaf grazing mortality carbon fluxes
                       grazing_c_to_litr_c(c,j,i) = &
                          grazing_c_to_litr_c(c,j,i) + &
                          grz_leafc_to_litter(p) * cfrac2litr(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)

                       ! fine root grazing mortality carbon fluxes
                       grazing_c_to_litr_c(c,j,i) = &
                          grazing_c_to_litr_c(c,j,i) + &
                          grz_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
                    end do

                    ! wood grazing mortality carbon fluxes
                    grazing_c_to_cwdc(c,j)  = grazing_c_to_cwdc(c,j)  + &
                         (grz_livestemc_to_litter(p) + grz_deadstemc_to_litter(p)) * &
                          cfrac2litr(p) * wtcol(p) * stem_prof(p,j) 
                    grazing_c_to_cwdc(c,j) = grazing_c_to_cwdc(c,j) + &
                         grz_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    grazing_c_to_cwdc(c,j) = grazing_c_to_cwdc(c,j) + &
                         grz_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 

                    ! storage grazing mortality carbon fluxes
                    ! Metabolic litter is treated differently than other types
                    ! of litter, so it gets this additional line after the
                    ! most recent loop over all litter types
                    grazing_c_to_litr_c(c,j,i_met_lit) = &
                       grazing_c_to_litr_c(c,j,i_met_lit) + &
                       grz_leafc_storage_to_litter(p) * cfrac2litr(p) * wtcol(p) * leaf_prof(p,j) + &
                       grz_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                       grz_livestemc_storage_to_litter(p) * cfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_deadstemc_storage_to_litter(p) * cfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       grz_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       grz_gresp_storage_to_litter(p) * cfrac2litr(p) * wtcol(p) * leaf_prof(p,j) + &

                    ! transfer grazing mortality carbon fluxes
                       grz_leafc_xfer_to_litter(p) * cfrac2litr(p) * wtcol(p) * leaf_prof(p,j) + &
                       grz_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                       grz_livestemc_xfer_to_litter(p) * cfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_deadstemc_xfer_to_litter(p) * cfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       grz_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       grz_gresp_xfer_to_litter(p) * cfrac2litr(p) * wtcol(p) * leaf_prof(p,j)

                    do i = i_litr_min, i_litr_max
                       grazing_n_to_litr_n(c,j,i) = &
                          grazing_n_to_litr_n(c,j,i) + &
                          ! leaf grazing mortality nitrogen fluxes
                          grz_leafn_to_litter(p) * nfrac2litr(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                          ! fine root litter nitrogen fluxes
                          grz_frootn_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
                    end do

                    ! wood grazing mortality nitrogen fluxes
                    grazing_n_to_cwdn(c,j)  = grazing_n_to_cwdn(c,j)  + &
                         (grz_livestemn_to_litter(p) + grz_deadstemn_to_litter(p)) * nfrac2litr(p) * &
                          wtcol(p) * stem_prof(p,j)
                    grazing_n_to_cwdn(c,j) = grazing_n_to_cwdn(c,j) + &
                         grz_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
                    grazing_n_to_cwdn(c,j) = grazing_n_to_cwdn(c,j) + &
                         grz_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)

                    ! Metabolic litter is treated differently than other types
                    ! of litter, so it gets this additional line after the
                    ! most recent loop over all litter types
                    grazing_n_to_litr_n(c,j,i_met_lit) = &
                       grazing_n_to_litr_n(c,j,i_met_lit) + &
                       ! retranslocated N pool grazing mortality fluxes
                       grz_retransn_to_litter(p) * nfrac2litr(p) * wtcol(p) * leaf_prof(p,j) + &
                       ! storage grazing mortality nitrogen fluxes
                       grz_leafn_storage_to_litter(p) * nfrac2litr(p) * wtcol(p) * leaf_prof(p,j) + &
                       grz_frootn_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                       grz_livestemn_storage_to_litter(p) * nfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_deadstemn_storage_to_litter(p) * nfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_livecrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       grz_deadcrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       ! transfer grazing mortality nitrogen fluxes
                       grz_leafn_xfer_to_litter(p) * nfrac2litr(p) * wtcol(p) * leaf_prof(p,j) + &
                       grz_frootn_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                       grz_livestemn_xfer_to_litter(p) * nfrac2litr(p) * wtcol(p) * stem_prof(p,j) + &
                       grz_deadstemn_xfer_to_litter(p)  * nfrac2litr(p)* wtcol(p) * stem_prof(p,j) + &
                       grz_livecrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                       grz_deadcrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)

                 end if
              end if

           end do

        end do
     end do
   
        do fc = 1,num_soilc
           c = filter_soilc(fc)
              if (col%active(c)) then
                 ! total grazed biomass and losses
                  grazed_totcc(c) = 0.0_r8
                  grazed_clossc(c) = 0.0_r8
                  grazed_totnc(c) = 0.0_r8
                  grazed_nlossc(c) = 0.0_r8
              end if
        end do


     do pi = 1,maxsoil_patches
        do fc = 1,num_soilc
           c = filter_soilc(fc)

           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1

              if (patch%active(p)) then
                 ! total grazed biomass and losses
                  grazed_totcc(c) = grazed_totcc(c)  + &
                      grazed_totcp(p) * wtcol(p)
                  grazed_clossc(c) = grazed_clossc(c)  + &
                      grazed_clossp(p) * wtcol(p)
                  grazed_totnc(c) = grazed_totnc(c)  + &
                      grazed_totnp(p) * wtcol(p)
                  grazed_nlossc(c) = grazed_nlossc(c)  + &
                      grazed_nlossp(p) * wtcol(p)
              end if
           end if

        end do

     end do

   end associate 

 end subroutine CNGrazingPftToColumn

end module CNGrazingMod

