from subprocess import run
import os
import glob
import argparse
import sys
import xarray as xr
import numpy as np
import logging

# -- add python/ctsm  to path (needed if we want to run regrid_ggcmi_shdates stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm.utils import abort, import_coord_1d, import_coord_2d
from ctsm import ctsm_logging

logger = logging.getLogger(__name__)


def main():
    """
    Description
    -----------
    Calls function that regrids GGCMI sowing and harvest dates.
    """

    args = regrid_ggcmi_shdates_arg_process()
    regrid_ggcmi_shdates(
        args.regrid_resolution,
        args.regrid_template_file,
        args.regrid_input_directory,
        args.regrid_output_directory,
        args.regrid_extension,
        args.crop_list,
    )


def run_and_check(cmd):
    result = run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        abort(f"Trouble running `{result.args}` in shell:\n{result.stdout}\n{result.stderr}")


# Functionized because these are shared by process_ggcmi_shdates
def define_arguments(parser):
    # Required
    parser.add_argument(
        "-rr",
        "--regrid-resolution",
        help="Target CLM resolution, to be saved in output filenames.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-rt",
        "--regrid-template-file",
        help="Template netCDF file to be used in regridding of inputs. This can be a CLM output file (i.e., something with 1-d lat and lon variables) or a CLM surface dataset (i.e., something with 2-d LATIXY and LONGXY variables).",
        type=str,
        required=True,
    )
    
    default = ".nc"
    parser.add_argument(
        "-x",
        "--regrid-extension",
        help=f"File regrid_extension of raw GGCMI sowing/harvest date files (default {default}).",
        default=default,
    )
    parser.add_argument(
        "-c",
        "--crop-list",
        help="List of GGCMI crops to process; e.g., '--crop-list mai_rf,mai_ir'. If not provided, will process all GGCMI crops.",
        default=None,
    )
    return parser


def regrid_ggcmi_shdates(
    regrid_resolution,
    regrid_template_file_in,
    regrid_input_directory,
    regrid_output_directory,
    regrid_extension,
    crop_list,
):
    logger.info(f"Regridding GGCMI crop calendars to {regrid_resolution}:")

    # Ensure we can call necessary shell script(s)
    for cmd in ["module load cdo; cdo"]:
        run_and_check(f"{cmd} --help")

    previous_dir = os.getcwd()
    os.chdir(regrid_input_directory)
    if not os.path.exists(regrid_output_directory):
        os.makedirs(regrid_output_directory)

    templatefile = os.path.join(regrid_output_directory, "template.nc")
    if os.path.exists(templatefile):
        os.remove(templatefile)

    template_ds_in = xr.open_dataset(regrid_template_file_in)

    # Process inputs
    if crop_list is not None:
        crop_list = crop_list.split(",")
    if regrid_extension[0] != ".":
        regrid_extension = "." + regrid_extension

    # Import and format latitude
    if "lat" in template_ds_in:
        lat, Nlat = import_coord_1d(template_ds_in, "lat")
    elif "LATIXY" in template_ds_in:
        lat, Nlat = import_coord_2d(template_ds_in, "lat", "LATIXY")
        lat.attrs["axis"] = "Y"
    else:
        abort("No latitude variable found in regrid template file")

    # Flip latitude, if needed
    if lat.values[0] < lat.values[1]:
        lat = lat.reindex(lat=list(reversed(lat["lat"])))

    # Import and format longitude
    if "lon" in template_ds_in:
        lon, Nlon = import_coord_1d(template_ds_in, "lon")
    elif "LONGXY" in template_ds_in:
        lon, Nlon = import_coord_2d(template_ds_in, "lon", "LONGXY")
        lon.attrs["axis"] = "Y"
    else:
        abort("No longitude variable found in regrid template file")
    template_da_out = xr.DataArray(
        data=np.full((Nlat, Nlon), 0.0),
        dims={"lat": lat, "lon": lon},
        name="area",
    )

    # Save template Dataset for use by cdo
    template_ds_out = xr.Dataset(
        data_vars={
            "planting_day": template_da_out,
            "maturity_day": template_da_out,
            "growing_season_length": template_da_out,
        },
        coords={"lat": lat, "lon": lon},
    )
    template_ds_out.to_netcdf(templatefile, mode="w")

    # Loop through original crop calendar files, interpolating using cdo with nearest-neighbor
    pattern = "*" + regrid_extension
    input_files = glob.glob(pattern)
    if len(input_files) == 0:
        abort(f"No files found matching {os.path.join(os.getcwd(), pattern)}")
    input_files.sort()
    for f in input_files:
        this_crop = f[0:6]
        if crop_list is not None and this_crop not in crop_list:
            continue

        logger.info("    " + this_crop)
        f2 = os.path.join(regrid_output_directory, f)
        f3 = f2.replace(regrid_extension, f"_nninterp-{regrid_resolution}{regrid_extension}")

        if os.path.exists(f3):
            os.remove(f3)

        # Sometimes cdo fails for no apparent reason. In testing this never happened more than 3x in a row.
        try:
            run_and_check(f"module load cdo; cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")
        except:
            try:
                run_and_check(f"module load cdo; cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")
            except:
                try:
                    run_and_check(f"module load cdo; cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")
                except:
                    run_and_check(f"module load cdo; cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")

    # Delete template file, which is no longer needed
    os.remove(templatefile)
    os.chdir(previous_dir)


def regrid_ggcmi_shdates_arg_process():
    """Process input arguments

    Returns:
        argparse.ArgumentParser: Arguments/options
    """

    # set up logging allowing user control
    ctsm_logging.setup_logging_pre_config()

    parser = argparse.ArgumentParser(
        description="Regrids raw sowing and harvest date files provided by GGCMI to a target CLM resolution."
    )

    # Define arguments
    parser = define_arguments(parser)
    parser.add_argument(
        "-i",
        "--regrid-input-directory",
        help="Directory containing the raw GGCMI sowing/harvest date files.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--regrid-output-directory",
        help="Directory where regridded output files should be saved.",
        type=str,
        required=True,
    )
    ctsm_logging.add_logging_args(parser)

    # Get arguments
    args = parser.parse_args(sys.argv[1:])
    ctsm_logging.process_logging_args(args)

    # Process arguments
    args.regrid_template_file = os.path.realpath(args.regrid_template_file)
    args.regrid_input_directory = os.path.realpath(args.regrid_input_directory)
    args.regrid_output_directory = os.path.realpath(args.regrid_output_directory)

    return args

