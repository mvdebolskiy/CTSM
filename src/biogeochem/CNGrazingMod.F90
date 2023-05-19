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
  use CNGrazerType            , only : grazer_type
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
  


  real(r8), private, parameter :: test_param = 10.0_r8   ! testing parameter
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  
  subroutine CNGrazing(bounds, grazer_inst)

  ! USES:

  ! ARGUMENTS:

  type(bounds_type)              ,  intent(in)    :: bounds
  type(grazer_type)              ,  intent(inout) :: grazer_inst


  ! LOCAL VARIABLES:
  integer p

  call t_startf( 'CNGrqazing' )

  associate( &
            is_grazed               =>    grazer_inst%is_grazed                , & ! Input   [logical  (:)   ] patches that are grazed
            totc_grazedp            =>    grazer_inst%grazed_totc_patch        , & ! In/out  [real(r8) (:)   ] patch grazed biomass gC/m2s
            begc                    =>    bounds%begc                          , & ! Input   [integer        ] beginning column index
            begp                    =>    bounds%begp                          , & ! Input   [integer        ] beginning patch index
            endc                    =>    bounds%endc                          , & ! Input   [integer        ] ending column index
            endp                    =>    bounds%endp                            & ! Input   [integer        ] ending patch index
           )

    do p = begp, endp
      if (is_grazed(p)) then
        totc_grazedp = 1.0_r8
      write(iulog,*)
      endif
    enddo

  call TestGrazRoutines()

  end associate
  call t_stopf( 'CNGrqazing' )
  end subroutine CNGrazing
  
  
  
  subroutine TestGrazRoutines()

  call t_startf( 'TestGrazRoutines' )


  call t_stopf( 'TestGrazRoutines' )
  end subroutine TestGrazRoutines

end module CNGrazingMod

