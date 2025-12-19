import argparse
import logging
import sys
import time
import traceback
from pathlib import Path

from tilepy.include.CampaignDefinition import ObservationParameters
from tilepy.tools.VisualizationTools import CompareTwoTilings

__all__ = ["PlottingTwoCampaigns"]

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


def PlottingTwoCampaigns(obspar, PointingsFile1, PointingsFile2):
    return_code = 0
    try:
        plotType = "gnomonic"
        CompareTwoTilings(
            obspar.skymap, PointingsFile1, PointingsFile2, obspar.FOV, plotType
        )
        plotType = "mollweide"
        CompareTwoTilings(
            obspar.skymap, PointingsFile1, PointingsFile2, obspar.FOV, plotType
        )
    except Exception:
        logger.error(
            f"An error occurred during the execution:\n{traceback.format_exc()}"
        )
        return_code = 1

    return return_code


def main():
    start = time.time()

    parser = argparse.ArgumentParser(
        description="Start the LST pointing observation of a GW event"
    )
    parser.add_argument(
        "-skymap",
        metavar="skymap",
        default="https://gracedb.ligo.org/api/superevents/MS230522j/files/bayestar.fits.gz",
        help="FITS file with the sky localization, e.g.for GW https://urlpath/Bayestar.fits.gz",
    )
    parser.add_argument(
        "-time",
        metavar='"YYYY-MM-DD HH:MM:SS"',
        default="2023-07-27 08:30:10",
        help="optional: date and time of the event (default: NOW, i.e. %(default)s)",
    )
    parser.add_argument(
        "-file1",
        metavar="first scheduling",
        help="file with the first schedule to plot",
    )
    parser.add_argument(
        "-file2",
        metavar="second scheduling",
        help="file with the second schedule to plot",
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
        default="plotting_two_campaigns.log",
    )
    parser.add_argument(
        "-vetoWindowsFile",
        help="File containing time windows to exclude from the computation..",
        default=None,
    )

    args = parser.parse_args()
    skymap = args.skymap
    obsTime = args.time
    PointingsFile1 = args.file1
    PointingsFile2 = args.file2
    datasetDir = args.i
    outDir = args.o
    cfgFile = args.cfg
    galcatName = args.galcatName
    pointingsFile = args.tiles
    eventName = args.eventName
    logname = args.logname
    vetoWindowsFile = args.vetoWindowsFile

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
        "healpix",
        None,
        None,
        None,
        None,
        vetoWindowsFile,
    )
    obspar.from_configfile(cfgFile)

    return_code = PlottingTwoCampaigns(obspar, PointingsFile1, PointingsFile2)

    end = time.time()
    logger.info(f"Execution time: {end - start:.0f} [sec]")
    logger.info(f"Return code: {return_code}")
    sys.exit(return_code)


if __name__ == "__main__":
    main()
