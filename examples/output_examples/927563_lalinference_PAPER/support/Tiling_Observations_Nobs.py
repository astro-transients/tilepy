#####################################################################
#  Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler  #
#  Script to obtained the Pointing observations of a GW follow-up  #
#  with several sites  #
#####################################################################


from tilepy.include.ObservationScheduler import getdate, GetUniversalSchedule
from tilepy.include.PointingTools import ObservationParameters, NextWindowTools
import time
import argparse
import os


start = time.time()

###########################
#####    Parsing  ######
###########################

parser = argparse.ArgumentParser(description='Start the LST pointing observation of a GW event')
parser.add_argument('-alertType',metavar = 'type of followup', help = 'options: gbm, gbmpng or gw', default= 'gw')
parser.add_argument('-url', metavar='skymap', default = 'https://gracedb.ligo.org/api/superevents/MS230522j/files/bayestar.fits.gz',
                    help='the url to the FITS file with the sky localization, e.g.for GW https://urlpath/Bayestar.fits.gz')
parser.add_argument('-time', metavar='\"YYYY-MM-DD HH:MM:SS\"', default= "2023-03-15 10:30:10",
                    help='optional: date and time of the event (default: NOW, i.e. %(default)s)')
parser.add_argument('-i',metavar = 'input path', help='Path to the input datasets (where galaxy cat should be for GW case)', default = "./dataset/")
parser.add_argument('-o',metavar = 'output path', help='Path to the output folder',default='./output')
parser.add_argument('-cfg',metavar = 'config file', help='Config file for the tiling scheduling',default='./config/FollowupParameters.ini')
parser.add_argument('-galcatName', metavar='galaxy catalog name', default="converted_GLADE.h5")
parser.add_argument('-tiles', metavar='tiles already observed', default= None)
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

#url = 'https://gracedb.ligo.org/api/superevents/S190814bv/files/LALInference.v1.fits.gz'
#url = 'https://dcc.ligo.org/public/0119/P1500071/007/939340_lalinference.fits.gz'
#url = 'https://dcc.ligo.org/public/0119/P1500071/007/644826_lalinference.fits.gz'
#url = 'https://dcc.ligo.org/public/0119/P1500071/007/644826_lalinference.fits.gz'
#url = 'https://dcc.ligo.org/public/0119/P1500071/007/790258_lalinference.fits.gz'
#url = 'https://dcc.ligo.org/public/0119/P1500071/007/939340_lalinference.fits.gz'
#url = 'https://gracedb.ligo.org/api/superevents/MS231107m/files/bayestar.multiorder.fits,1'
url = 'https://dcc.ligo.org/public/0119/P1500071/007/927563_lalinference.fits.gz'
#url = 'https://gracedb.ligo.org/api/superevents/S190728q/files/GW190728_064510_PublicationSamples_flattened.fits.gz,0'
#url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2023/bn231012231/quicklook/glg_healpix_all_bn231012231.fit'
#url = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn221126547/quicklook/glg_locplot_all_bn221126547.png'
#url = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn221119627/quicklook/glg_locplot_all_bn221119627.png'
#url = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn220424481/quicklook/glg_locplot_all_bn220424481.png'
#url = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn220627901/quicklook/glg_locplot_all_bn220627901.png'
#url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn220627890/quicklook/glg_locplot_all_bn220627890.png'
#url = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn220424481/quicklook/glg_locplot_all_bn220424481.png'
#url = 'https://gracedb.ligo.org/api/superevents/S230924an/files/bayestar.multiorder.fits,1'
#url = 'https://gracedb.ligo.org/api/superevents/MS230826n/files/bayestar.multiorder.fits,1'
#url = 'https://dcc.ligo.org/public/0146/G1701985/001/LALInference_v2.fits.gz'
ObsArray = ['CTAN', 'CTAS']
parameters = []

for i in ObsArray:
    parameters.append("./config/FollowupParameters_%s.ini" % i)
print("===========================================================================================")
print('parameters', parameters)
obsparameters = []

for j in range(len(parameters)):
    obspar = ObservationParameters()
    obspar.add_parsed_args(url, obsTime, datasetDir, galcatName, outDir, pointingsFile, alertType, locCut)
    obspar.from_configfile(parameters[j])
    obsparameters.append(obspar)

GetUniversalSchedule(obsparameters)



end = time.time()
print('Execution time: ', end - start)





