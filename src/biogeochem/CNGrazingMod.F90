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
  public :: TestGrazRoutines       ! soubroutine to test this module
  !
  !
  ! The following is only public for the sake of unit testing; it should not be called
  ! directly by CLM code outside this module

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  


  real(r8), private, parameter :: test_param = 10.0_r8   ! 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  
  subroutine CNGrazing(bounds)

  ! ARGUMENTS

  type(bounds_type)              ,  intent(in)    :: bounds

  call t_startf( 'CNGrqazing' )

  call TestGrazRoutines()


  call t_stopf( 'CNGrqazing' )
  end subroutine CNGrazing
  
  
  
  subroutine TestGrazRoutines()

  call t_startf( 'TestGrazRoutines' )


  call t_stopf( 'TestGrazRoutines' )
  end subroutine TestGrazRoutines

end module CNGrazingMod

