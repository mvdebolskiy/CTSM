module ExcessIceStreamMod

    ! **********************************************************************
    ! --------------------------- IMPORTANT NOTE ---------------------------
    !
    ! This file is here temporarily in order to exercise the CDEPS stream code for this 3-d
    ! stream. In cases using the NUOPC driver/mediator, this version is used instead of the
    ! version in src/biogeophys. Changes to the science here should also be made in the
    ! similar file in src/biogeophys. Once we start using CDEPS by default, we can remove
    ! the version in src/biogeophys and move this version into there.
    ! **********************************************************************
  
#include "shr_assert.h"

!-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in soil moisture from data stream
  !
  ! !USES:
  use ESMF
  use dshr_strdata_mod   , only : shr_strdata_type, shr_strdata_print
  use dshr_strdata_mod   , only : shr_strdata_init_from_inline, shr_strdata_advance
  use dshr_methods_mod   , only : dshr_fldbun_getfldptr
  use dshr_stream_mod    , only : shr_stream_file_null
  use shr_kind_mod       , only : r8 => shr_kind_r8, cl => shr_kind_CL, cxx => shr_kind_CXX
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use shr_mpi_mod        , only : shr_mpi_bcast
  use decompMod          , only : bounds_type
  use abortutils         , only : endrun
  use clm_varctl         , only : iulog, use_soil_moisture_streams
  use controlMod         , only : NLFilename
  use LandunitType       , only : lun
  use ColumnType         , only : col                
  use SoilStateType      , only : soilstate_type
  use WaterStateBulkType , only : waterstatebulk_type
  use perf_mod           , only : t_startf, t_stopf
  use spm  dMod            , only : masterproc, mpicom, iam
  use clm_time_manager   , only : get_calendar, get_curr_date
  use clm_nlUtilsMod     , only : find_nlgroup_name
  use clm_varpar         , only : nlevsoi
  use clm_varcon         , only : denh2o, denice, watmin, spval
  use landunit_varcon    , only : istsoil, istcrop
  use lnd_comp_shr       , only : mesh, model_meshfile, model_clock
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ExcessIceInitStream ! excess ice initialization and file
  public :: ExcessIceInit

  ! !PRIVATE MEMBER DATA:
  type(shr_strdata_type)  :: sdat_exice                    ! excess ice input data stream
  character(*), parameter :: stream_var_name = "EXICE"    ! base string for field string
  integer, allocatable    :: g_to_ig(:)                    ! Array matching gridcell index to data index
  !PRIVATE TYPES:
  character(len=*), parameter, private :: sourcefile = &
         __FILE__
  character(*), parameter :: u_FILE_u = &
         __FILE__

contains
  
  !-----------------------------------------------------------------------
  !
  ! excess_ice_init_stream
  ! initializes the data stream
  !-----------------------------------------------------------------------
  subroutine ExcessIceInitStream(bounds)
    integer            :: rc                         ! error coe
    integer            :: i                          ! index
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    character(len=CL)  :: stream_fldfilename_exice   ! ustar stream filename to read
    character(len=CL)  :: stream_mapalgo = 'bilinear' !interpolation algorithm
    character(*), parameter :: subName = "('ExcessIceInitStream')"

    if (masterproc) then
        write(iulog,*) ' '
        write(iulog,*) 'Excess ice stream initialization started'
    endif
    
  end subroutine ExcessIceInitStream


  subroutine ExcessIceInit(bounds,soilstate_inst, waterstatebulk_inst)
    !
    ! Assign data stream information for prescribed soil moisture.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(soilstate_type)      , intent(in)    :: soilstate_inst
    type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: rc
    integer  :: c, g, j, ig, n


    character(*), parameter :: subName = "('ExcessIceInit')"

    if (masterproc) then
        write(iulog,*) ' '
        write(iulog,*) 'Excess ice initialization started'
    endif


  end subroutine ExcessIceInit
end module ExcessIceStreamMod