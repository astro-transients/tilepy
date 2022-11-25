############################################################################
#  Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler          #
#  Script to obtained the Pointing observations of a GBM tiling follow-up  #
############################################################################
# The script should be called as python Launch_SingleTelescope_Tiling_GBM.py -url name.gz -time "YYYY-MM-DD HH:MM:SS" -i input_path -o output_path

from tilepy.include.ObservationScheduler import GetSchedule_GBMfromPNG, getdate
import time
import argparse
import os
#import datetime
#from datetime import timezone


start = time.time()

###########################
#####    Parsing  ######
###########################

parser = argparse.ArgumentParser(description='Start the LST pointing observation of a GBM event')
parser.add_argument('-url', metavar='GBMmap', default = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2021/bn211130636/quicklook/glg_skymap_all_bn211130636.png',
                    help='the url to the png file with the GBM localization, e.g. https://urlpath/skymap.png')
parser.add_argument('-time', metavar='\"YYYY-MM-DD HH:MM:SS\"', default= "2022-12-13 10:30:10",
                    help='optional: date and time of the event (default: NOW, i.e. %(default)s)')
parser.add_argument('-i',metavar = 'input path', help='Path to the input datasets', default = "./datasets")
parser.add_argument('-o',metavar = 'output path', help='Path to the output folder', default="./output")
parser.add_argument('-cfg',metavar = 'config file', help='Config file for the tiling scheduling',default='./config/testConfig.ini')


#default= "2021-12-05 10:10:10"
args = parser.parse_args()
url = args.url
ObsTime = getdate(args.time)
datasetDir = args.i
outDir = args.o
cfgFile = args.cfg

if not os.path.exists(outDir):
    os.makedirs(outDir)


################################################

GetSchedule_GBMfromPNG(url, ObsTime, datasetDir,outDir,cfgFile)

end = time.time()
print('Execution time: ', end - start)










