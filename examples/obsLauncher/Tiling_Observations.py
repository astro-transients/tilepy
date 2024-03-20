########################################################################
#    Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler    #
#  Script to obtained the pointing observations of a GW/GRB follow-up  #
#  ------------------------- Version 1.3.0 -------------------------   #
########################################################################

from tilepy.include.ObservationScheduler import getSchedule
from tilepy.include.PointingTools import ObservationParameters, NextWindowTools, getdate
import time
import argparse
import os
from pathlib import Path



start = time.time()

###########################
#####    Parsing  ######
###########################

parser = argparse.ArgumentParser(description='Start the LST pointing observation of a GW event')
parser.add_argument('-alertType',metavar = 'type of followup', help = 'options: gbm, gbmpng or gw')
parser.add_argument('-url', metavar='skymap', default = 'https://gracedb.ligo.org/api/superevents/MS230522j/files/bayestar.fits.gz',
                    help='the url to the FITS file with the sky localization, e.g.for GW https://urlpath/Bayestar.fits.gz')
parser.add_argument('-time', metavar='\"YYYY-MM-DD HH:MM:SS\"', default= "2022-12-14 14:30:10",
                    help='optional: date and time of the event (default: NOW, i.e. %(default)s)')
parser.add_argument('-i',metavar = 'input path', help='Path to the input datasets (where galaxy cat should be for GW case)', default = "../dataset/")
parser.add_argument('-o',metavar = 'output path', help='Path to the output folder',default='./output')
parser.add_argument('-cfg',metavar = 'config file', help='Config file for the tiling scheduling',default='./config/FollowupParameters.ini')
parser.add_argument('-galcatName', metavar='galaxy catalog name', default=None)
parser.add_argument('-tiles', metavar='tiles already observed', default= "False")
parser.add_argument('-locCut', metavar='limit on skyloc to perform a followup', help='Options are: loose or std', default=None)

args = parser.parse_args()
alertType = args.alertType
url = args.url
obsTime = getdate(args.time)
datasetDir = args.i
outDir = args.o
cfgFile = args.cfg
galcatName = args.galcatName
pointingsFile = args.tiles
locCut = args.locCut

if not os.path.exists(outDir):
    os.makedirs(outDir)

################################################
#url = 'https://gracedb.ligo.org/api/superevents/S190728q/files/GW190728_064510_PublicationSamples_flattened.fits.gz,0'

obspar = ObservationParameters()
obspar.add_parsed_args(url,obsTime,datasetDir,galcatName,outDir,pointingsFile,alertType,locCut)
obspar.from_configfile(cfgFile)

#GetSchedule_confile(url,ObsTime,datasetDir,galcatname,outDir,cfgFile,PointingsFile,type)
getSchedule(obspar)

end = time.time()
print('Execution time: ', end - start)