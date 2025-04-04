import argparse
import os
import time

from tilepy.include.CampaignDefinition import ObservationParameters

from tilepy.tools.VisualizationTools import CompareTwoTilings

__all__ = ["PlottingTwoCampaigns"]


def PlottingTwoCampaigns(obspar, PointingsFile1, PointingsFile2):
    plotType = "gnomonic"
    CompareTwoTilings(
        obspar.skymap, PointingsFile1, PointingsFile2, obspar.FOV, plotType
    )
    plotType = "mollweide"
    CompareTwoTilings(
        obspar.skymap, PointingsFile1, PointingsFile2, obspar.FOV, plotType
    )


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
        "-locCut",
        metavar="limit on skyloc to perform a followup",
        help="Options are: loose or std",
        default=None,
    )
    parser.add_argument(
        "-eventName", metavar="Name of the observed event", default=None
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

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    obspar = ObservationParameters()
    obspar.add_parsed_args(
        skymap, obsTime, datasetDir, galcatName, outDir, pointingsFile, eventName
    )
    obspar.from_configfile(cfgFile)

    PlottingTwoCampaigns(obspar, PointingsFile1, PointingsFile2)

    end = time.time()
    print("Execution time: ", end - start)


if __name__ == "__main__":
    main()
