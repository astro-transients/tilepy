########################################################################
#    Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler    #
#  Script to obtained the pointing observations of a GW/GRB follow-up  #
#  ------------------------- Version 1.3.0 -------------------------   #
########################################################################

import argparse
import datetime
import os
import time

from tilepy.include.ObservationScheduler import GetSchedule
from tilepy.include.CampaignDefinition import ObservationParameters

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
    "-i",
    metavar="input path",
    help="Path to the input datasets (where galaxy cat should be for GW case)",
    default="../../dataset/",
)
parser.add_argument(
    "-o", metavar="output path", help="Path to the output folder", default="./output"
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
parser.add_argument("-eventName", metavar="Name of the observed event", default=None)
args = parser.parse_args()
url = args.skymap
obsTime = datetime.datetime.fromisoformat(args.time)
datasetDir = args.i
outDir = args.o
cfgFile = args.cfg
galcatName = args.galcatName
pointingsFile = args.tiles
eventName = args.eventName

if not os.path.exists(outDir):
    os.makedirs(outDir)

skymap = "https://gracedb.ligo.org/api/superevents/S190728q/files/GW190728_064510_PublicationSamples_flattened.fits.gz,0"

obspar = ObservationParameters()
obspar.add_parsed_args(
    skymap, obsTime, datasetDir, galcatName, outDir, pointingsFile, eventName
)
obspar.from_configfile(cfgFile)

GetSchedule(obspar)

end = time.time()
print("Execution time: ", end - start)
