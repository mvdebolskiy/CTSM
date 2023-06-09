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
  use landunit_varcon , only : istsoil, istcrop, istdlak 
  use GridcellType    , only : grc
  use LandunitType    , only : lun
  use ColumnType      , only : col
  use PatchType       , only : patch
  use pftconMod       , only : pftcon

 implicit none
  save
  private
  !
  type, public :: grazer_type

    integer              num_grazed_pfts                 ! number of pfts to graze
    integer,  pointer :: grazed_pfts(:)                  ! pfts to be grazed
    real(r8), pointer :: grazed_totc_patch               (:) ! patch grazed carbon flux gC/m2/s
    real(r8), pointer :: grazed_closs_patch              (:) ! patch carbon loss from grazing gC/m2/s
    real(r8), pointer :: grazed_totc_col                 (:) ! column grazed carbon flux gC/m2/s
    real(r8), pointer :: grazed_closs_col                (:) ! column carbon loss flux from grazing gC/m2/s
    real(r8), pointer :: grazed_totn_patch               (:) ! patch grazed nitrogen flux gN/m2/s
    real(r8), pointer :: grazed_nloss_patch              (:) ! patch nitrogen loss flux from grazing gN/m2/s
    real(r8), pointer :: grazed_totn_col                 (:) ! column grazed nitrogen flux gN/m2/s
    real(r8), pointer :: grazed_nloss_col                (:) ! col nitrogen loss flux from grazing gN/m2/s
    

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
    call this%InitHistory (bounds)
    call this%InitCold (bounds)


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
    integer :: p,j
    !------------------------------------------------------------------------

    !allocate(this%grazed_pfts  (1:this%num_grazed_pfts))                      ; this%grazed_pfts                     (:)   = ispval
    allocate(this%grazed_totc_patch (bounds%begp:bounds%endp))                ; this%grazed_totc_patch      (:) = nan
    allocate(this%grazed_closs_patch (bounds%begp:bounds%endp))                ; this%grazed_closs_patch      (:) = nan
    allocate(this%grazed_totn_patch (bounds%begp:bounds%endp))                ; this%grazed_totn_patch      (:) = nan
    allocate(this%grazed_nloss_patch (bounds%begp:bounds%endp))                ; this%grazed_nloss_patch      (:) = nan
    allocate(this%grazed_totc_col   (bounds%begc:bounds%endc))                ; this%grazed_totc_col        (:) = nan
    allocate(this%grazed_closs_col   (bounds%begc:bounds%endc))                ; this%grazed_closs_col        (:) = nan
    allocate(this%grazed_totn_col   (bounds%begc:bounds%endc))                ; this%grazed_totn_col        (:) = nan
    allocate(this%grazed_nloss_col   (bounds%begc:bounds%endc))                ; this%grazed_nloss_col        (:) = nan

    j = 0
    do p = 0, mxpft
      if ( pftcon%is_grass(p) == .true. ) then
          j = j +1
      endif
    end do
    this%num_grazed_pfts = j
    write(iulog, *) 'num grazed pfts:   ', this%num_grazed_pfts
    allocate(this%grazed_pfts(1:this%num_grazed_pfts)) ; this%grazed_pfts(:) = ispval

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
         avgflag='A', long_name='grazed aboveground carbon', &
         ptr_col=this%grazed_totc_col)

    this%grazed_totc_patch(begp:endp) = spval
    call hist_addfld1d (fname='GRAZED_TOTC_PATCH', units='gC/m2s',  &
         avgflag='A', long_name='grazed aboveground carbon per patch', &
         ptr_patch=this%grazed_totc_patch)
    this%grazed_totn_col(begc:endc) = spval
    call hist_addfld1d (fname='GRAZED_TOTN_COL',  units='gN/m2s',  &
         avgflag='A', long_name='grazed aboveground nitrogen', &
         ptr_col=this%grazed_totc_col)
    this%grazed_totn_patch(begp:endp) = spval
    call hist_addfld1d (fname='GRAZED_TOTN_PATCH', units='gN/m2s',  &
         avgflag='A', long_name='grazed aboveground nitrogen per patch', &
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

    ! !ARGUMENTS:
    class(grazer_type)        :: this
    type(bounds_type) , intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer :: p, c, l, j, i
    integer :: fi                                        ! filter index
    !-----------------------------------------------------------------------

    ! Set column filters

    j = 0
    do p = 0, mxpft
      if ( pftcon%is_grass(p) == .true. ) then
          j = j +1
          this%grazed_pfts(j) = p
          write(iulog, *) 'grazed pft:   ', this%grazed_pfts(j)
      endif
    end do

    write(iulog, *) 'CNgrazed pfts:   ', this%grazed_pfts

    do p = bounds%begp, bounds%endp
       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%grazed_totc_patch(p) = 0.0_r8
          this%grazed_closs_patch(p) = 0.0_r8
          this%grazed_totn_patch(p) = 0.0_r8
          this%grazed_nloss_patch(p) = 0.0_r8
       endif
    enddo

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop ) then
          this%grazed_totc_col(c) = 0.0_r8
          this%grazed_closs_col(c) = 0.0_r8
          this%grazed_totn_col(c) = 0.0_r8
          this%grazed_nloss_col(c) = 0.0_r8
        endif
    enddo
    

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          write(iulog, * ) 'graz_totc_col', this%grazed_totc_col(c)
       end if
    end do
    do p = bounds%begp, bounds%endp
       l = patch%landunit(l)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          write(iulog, * ) 'graz_totc_patch', this%grazed_totc_patch(p)
       end if
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



  end subroutine Restart


end module CNGrazerType