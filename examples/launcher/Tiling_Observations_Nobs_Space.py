#####################################################################
#  Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler  #
#  Script to obtained the Pointing observations of a GW follow-up  #
#  with several sites  #
#####################################################################

from tilepy.include.ObservationScheduler import getdate, GetUniversalSchedule
from tilepy.include.PointingTools import ObservationParameters
import time
import argparse
import os


start = time.time()

###########################
#####    Parsing  ######
###########################

parser = argparse.ArgumentParser(description='Start the LST pointing observation of a GW event')
parser.add_argument('-alertType',metavar = 'type of followup', help = 'options: gbm, gbmpng or gw', default= 'gw')
parser.add_argument('-skymap', metavar='skymap', default = 'https://gracedb.ligo.org/api/superevents/MS230522j/files/bayestar.fits.gz',
                    help='FITS file with the sky localization, e.g.for GW https://urlpath/Bayestar.fits.gz')
parser.add_argument('-time', metavar='\"YYYY-MM-DD HH:MM:SS\"', default= "2023-07-27 01:30:10",
                    help='optional: date and time of the event (default: NOW, i.e. %(default)s)')
parser.add_argument('-i',metavar = 'input path', help='Path to the input datasets (where galaxy cat should be for GW case)', default = "../../dataset/")
parser.add_argument('-o',metavar = 'output path', help='Path to the output folder',default='./output')
parser.add_argument('-cfg',metavar = 'config file', help='Config file for the tiling scheduling',default='../config/FollowupParameters.ini')
parser.add_argument('-galcatName', metavar='galaxy catalog name', default="Gladeplus.h5")
parser.add_argument('-tiles', metavar='tiles already observed', default= None)
#parser.add_argument('-base', metavar='ground or space', default= "ground")

args = parser.parse_args()
alertType = args.alertType
skymap = args.skymap
obsTime = getdate(args.time)
datasetDir = args.i
outDir = args.o
cfgFile = args.cfg
galcatName = args.galcatName
pointingsFile = args.tiles

if not os.path.exists(outDir):
    os.makedirs(outDir)

################################################

skymap = 'https://gracedb.ligo.org/api/superevents/S190728q/files/GW190728_064510_PublicationSamples_flattened.fits.gz,0'

ObsArray = ['HESS', "SWIFT"]
parameters = []

for i in ObsArray:
    parameters.append("../config/FollowupParameters_%s.ini" % i)
print("===========================================================================================")
print('parameters', parameters)
obsparameters = []

for j in range(len(parameters)):
    obspar = ObservationParameters()
    obspar.add_parsed_args(skymap, obsTime, datasetDir, galcatName, outDir, pointingsFile, alertType)
    obspar.from_configfile(parameters[j])
    obsparameters.append(obspar)
    print("obsparameters.base",obsparameters[0].base)

GetUniversalSchedule(obsparameters)

end = time.time()
print('Execution time: ', end - start)





