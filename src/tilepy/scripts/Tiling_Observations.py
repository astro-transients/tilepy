########################################################################
#    Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler    #
#  Script to obtained the pointing observations of a GW/GRB follow-up  #
#  ------------------------- Version 1.3.0 -------------------------   #
########################################################################

import argparse
import logging
import sys
import time
import traceback

from astropy.time import Time

from pathlib import Path

from tilepy.include.CampaignDefinition import ObservationParameters
from tilepy.include.ObservationScheduler import GetSchedule
from tilepy.include.PointingTools import getdate

__all__ = ["Tiling_Observations"]

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


def Tiling_Observations(obspar):
    return_code = 0
    try:
        GetSchedule(obspar)
    except Exception:
        logger.error(
            f"An error occurred during the execution:\n{traceback.format_exc()}"
        )
        return_code = 1

    return return_code


def main():
    start = time.time()

    time_now = Time.now()
    time_now_str = time_now.datetime.strftime("%Y-%m-%d %H:%M:%S")

    parser = argparse.ArgumentParser(
        description="Start the pointing observation of a GW event"
    )
    parser.add_argument(
        "-skymap",
        metavar="skymap",
        default=None,
        help="FITS file with the sky localization, e.g.for GW https://urlpath/Bayestar.fits.gz",
    )
    parser.add_argument(
        "-mode",
        metavar="mode",
        default="healpix",
        choices=["gaussian", "locprob", "healpix"],
        help="Mode for reading maps.",
    )
    parser.add_argument(
        "-ra",
        metavar="ra",
        default=None,
        help="Right ascension of the target, in degrees.",
    )
    parser.add_argument(
        "-dec",
        metavar="dec",
        default=None,
        help="Declination of the target, in degrees",
    )
    parser.add_argument(
        "-sigma",
        metavar="sigma",
        default=None,
        help="Sigma for the creation of the gaussian map, in degrees",
    )
    parser.add_argument(
        "-nside",
        metavar="nside",
        default=None,
        help="Nside for the creation of the gaussian map.",
    )
    parser.add_argument(
        "-time",
        metavar='"YYYY-MM-DD HH:MM:SS"',
        default=time_now_str,
        help="optional: date and time of the event (default: NOW, i.e. %(default)s)",
    )
    parser.add_argument(
        "-i",
        metavar="input path",
        help="Path to the input datasets (where galaxy cat should be for GW case)",
        default="../../dataset/",
    )
    parser.add_argument(
        "-o",
        metavar="output path",
        help="Path to the output folder",
        default="./output",
    )
    parser.add_argument(
        "-cfg",
        metavar="config file",
        help="Config file for the tiling scheduling",
        default="../config/FollowupParameters_LST.ini",
    )
    parser.add_argument(
        "-galcatName", metavar="galaxy catalog name", default="Gladeplus.h5"
    )
    parser.add_argument("-tiles", metavar="tiles already observed", default=None)
    parser.add_argument(
        "-eventName", metavar="Name of the observed event", default=None
    )
    parser.add_argument(
        "-logname",
        metavar="Name of the output log file.",
        default="tiling_observations.log",
    )

    args = parser.parse_args()
    skymap = args.skymap
    mode = args.mode
    ra = args.ra
    dec = args.dec
    sigma = args.sigma
    nside = args.nside
    obsTime = getdate(args.time)
    datasetDir = args.i
    outDir = args.o
    cfgFile = args.cfg
    galcatName = args.galcatName
    pointingsFile = args.tiles
    eventName = args.eventName
    logname = args.logname

    if not Path(datasetDir).exists():
        raise FileNotFoundError(f"Dataset directory {datasetDir} not found.")

    galaxy_catalog = Path(f"{datasetDir}/{galcatName}")

    if not galaxy_catalog.exists():
        raise FileNotFoundError(f"Galaxy catalog file {galaxy_catalog} not found.")

    if not Path(cfgFile).exists():
        raise FileNotFoundError(f"Configuration file {cfgFile} not found.")

    if not Path(outDir).exists():
        Path(outDir).mkdir(parents=True)

    logging.basicConfig(filename=logname)

    if skymap is None and mode not in ["gaussian"]:
        raise ValueError(
            f"The skymap argument is mandatory when setting mode to {mode}."
        )

    attr_list_gaussian = ["ra", "dec", "sigma", "nside"]
    if mode == "gaussian":
        for i, attr in enumerate([ra, dec, sigma, nside]):
            if attr is None:
                raise ValueError(
                    f'"{attr_list_gaussian[i]}" must be defined in "gaussian" mode.'
                )

        if skymap is not None:
            raise ValueError('Cannot specify a skymap URL when mode is "gaussian".')

    ################################################

    obspar = ObservationParameters()
    obspar.add_parsed_args(
        skymap,
        obsTime,
        datasetDir,
        galcatName,
        outDir,
        pointingsFile,
        None,
        eventName,
        mode,
        ra,
        dec,
        sigma,
        nside,
    )
    obspar.from_configfile(cfgFile)

    return_code = Tiling_Observations(obspar)

    end = time.time()
    logger.info(f"Execution time: {end - start:.0f} [sec]")
    logger.info(f"Return code: {return_code}")
    sys.exit(return_code)


if __name__ == "__main__":
    main()
