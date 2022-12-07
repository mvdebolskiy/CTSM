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
  use filterColMod            , only : col_filter_from_ltypes
  use GridcellType            , only : grc
  use LandunitType            , only : lun
  use ColumnType              , only : col
  use PatchType               , only : patch
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNGrazing 
  !
  !
  ! The following is only public for the sake of unit testing; it should not be called
  ! directly by CLM code outside this module
  !public :: ComputeGroundHeatFluxAndDeriv       ! Computes G and dG/dT on surface of standing water, snow and soil

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: TestGrazRoutines       ! soubroutine to test this module


  real(r8), private, parameter :: test_param = 10.0_r8   ! 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  contains
  
  
  
  
  subroutine CNGrazing(bounds)

    filter_grazc = col_filter_from_ltypes( &
                 bounds = bounds, &
                 ltypes = [istsoil,istcrop], &
                 include_inactive = .true.)



  call t_startf( 'CNGrqazing' )

  call TestGrazRoutines()


  call t_startf( 'CNGrqazing' )
  end subroutine CNGrazing
  
  
  
  subroutine TestGrazRoutines()

  call t_startf( 'TestGrazRoutines' )


  call t_startf( 'CNGrqazing' )
  end subroutine TestGrazRoutines

end module CNGrazingMod

