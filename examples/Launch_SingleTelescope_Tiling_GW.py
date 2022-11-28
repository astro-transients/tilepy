####################################################################
#  Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler  #
#  Script to obtained the Pointing observations of a GW follow-up  #
####################################################################
# The script should be called as python Launch_SingleTelescope_Tiling_GW.py -url name.gz -time "YYYY-MM-DD HH:MM:SS" -i input_path -o output_path

from tilepy.include.ObservationScheduler import GetSchedule_GW, getdate
import time
import argparse
import os


start = time.time()

###########################
#####    Parsing  ######
###########################


parser = argparse.ArgumentParser(description='Start the pointing observation of a GW event')
parser.add_argument('-url', metavar='GWmap', default = 'https://gracedb.ligo.org/api/superevents/MS221112d/files/bayestar.fits.gz,1',
                    help='the url to the FITS file with the GW localization, e.g. https://urlpath/Bayestar.fits.gz')
parser.add_argument('-time', metavar='\"YYYY-MM-DD HH:MM:SS\"', default= "2022-12-14 00:30:10",
                    help='optional: date and time of the event (default: NOW, i.e. %(default)s)')
parser.add_argument('-i',metavar = 'input path', help='Path to the input datasets', default = "./datasets")
parser.add_argument('-o',metavar = 'output path', help='Path to the output folder',default='./output')
parser.add_argument('-cfg',metavar = 'config file', help='Config file for the tiling scheduling',default='./config/ExampleConfig.ini')

args = parser.parse_args()
url = args.url
ObsTime = getdate(args.time)
datasetDir = args.i
outDir = args.o
cfgFile = args.cfg

if not os.path.exists(outDir):
    os.makedirs(outDir)

################################################
#url = 'https://gracedb.ligo.org/api/superevents/S190728q/files/GW190728_064510_PublicationSamples_flattened.fits.gz,0'

GetSchedule_GW(url, ObsTime,datasetDir,outDir, cfgFile)

end = time.time()
print('Execution time: ', end - start)





