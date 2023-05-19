module CNGrazerType

#include "shr_assert.h"

  !------------------------------------------------------------------------------

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog, use_grazing
  use clm_varpar      , only : mxpft
  use clm_varcon      , only : spval, ispval
  use GridcellType    , only : grc
  use LandunitType    , only : lun
  use ColumnType      , only : col
  use PatchType       , only : patch

 implicit none
  save
  private
  !
  type, public :: grazer_type

    integer           :: num_grazed_pfts                 = 3 ! number of pfts to graze
    integer           :: grazed_pfts(3)                  = ispval ! pfts to be grazed
    logical,  pointer :: is_grazed                       (:) !
    real(r8), pointer :: grazed_totc_patch             (:) ! patch grazed biomass gC/m2s
    real(r8), pointer :: grazed_totc_col               (:) ! column grazed biomass gC/m2s

  contains

    procedure, public  :: Init
    procedure, public  :: Restart
    procedure, private :: InitAllocate
    procedure, private :: InitHistory
    procedure, private :: InitCold

  end type grazer_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  subroutine Init( this, bounds)

    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start.
    !
    class(grazer_type)             :: this
    type(bounds_type) , intent(in) :: bounds

    call this%InitAllocate (bounds)
    if (use_grazing) then
      call this%InitHistory (bounds)
      call this%InitCold (bounds)
    endif

  end subroutine Init

  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start.
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    class(grazer_type)        :: this
    type(bounds_type) , intent(in) :: bounds

    ! !LOCAL VARIABLES:

    !------------------------------------------------------------------------

    !allocate(this%grazed_pfts  (1:this%num_grazed_pfts))                      ; this%grazed_pfts                     (:)   = ispval
    allocate(this%grazed_totc_patch (bounds%begp:bounds%endp))                ; this%grazed_totc_patch      (:) = nan
    allocate(this%grazed_totc_col   (bounds%begc:bounds%endc))                ; this%grazed_totc_col        (:) = nan
    allocate(this%is_grazed           (bounds%begp:bounds%endp))                ; this%is_grazed                (:) = .false.

  end subroutine InitAllocate

  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup the fields that can be output on history files.
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_cn, use_cndv
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(grazer_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg


    this%grazed_totc_col(begc:endc) = spval
    call hist_addfld1d (fname='GRAZED_TOTC_COL',  units='gC/m2s',  &
         avgflag='A', long_name='grazed aboveground biomass', &
         ptr_col=this%grazed_totc_col)

    this%grazed_totc_patch(begp:endp) = spval
    call hist_addfld1d (fname='GRAZED_TOTC_PATCH', units='gC/m2s',  &
         avgflag='A', long_name='grazed aboveground biomass per patch', &
         ptr_patch=this%grazed_totc_patch)


  end subroutine InitHistory

  subroutine InitCold(this, bounds)

    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start.
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_cn, use_cndv
    use pftconMod      , only : pftcon

    ! !ARGUMENTS:
    class(grazer_type)        :: this
    type(bounds_type) , intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer :: l,c,p,j
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    j = 0
              write(iulog, *) 'num grazed pfts:   ', this%num_grazed_pfts
    do p = 0, mxpft
      if ( pftcon%is_grass(p) == .true. ) then
          j = j +1
          write(iulog, *) 'grazed pft:   ', p
          this%grazed_pfts(j) = p
          
      endif
    end do

      do p = bounds%begp, bounds%endp
        c = patch%column(p)
        l = patch%landunit(p)
        if ( ANY(this%grazed_pfts == patch%itype(p)) ) then
          this%grazed_totc_patch(p) = 1.0_r8
          this%grazed_totc_col(c) = 1.0_r8
          this%is_grazed(p) = .true.
        else
          this%grazed_totc_patch(p) = 0.0_r8
          this%grazed_totc_col(c) = 0.0_r8
        endif
      end do

  end subroutine InitCold

  subroutine Restart(this, bounds)

    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start.
    !

    ! !USES:
    use clm_time_manager, only : is_restart, is_first_step

    ! !ARGUMENTS:
    class(grazer_type)        :: this
    type(bounds_type) , intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg





  end subroutine Restart


end module CNGrazerType