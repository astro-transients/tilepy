import argparse
from ..include/TilingDetermination import *
import os

parser = argparse.ArgumentParser(description='Start the CTA GW pointing simulation and analysis (PgwinFov)')
parser.add_argument('index', metavar='GWmap',
                    help='the FITS file with the GW localization, e.g. Bayestar.fits.gz')
parser.add_argument('-params', metavar='parameters.ini', default="./configs/PGWinFoV.ini",
                    help='file holding the main parameters (default: %(default)s)')
args = parser.parse_args()

j=int(args.index)
parameters=args.params

# GW file
InputFileName = '../../dataset/BNS-GW_onAxis5deg.txt'
InputList = TableImportCTA(InputFileName)

GWFile = "../../dataset/skymaps/" + InputList['run'][j] + '_' + InputList['MergerID'][j] + "_skymap.fits.gz"
name =  InputList['run'][j] + '_' + InputList['MergerID'][j]

dirName='../output/PGWonFoV'

InjectTimeFile = '../../dataset/BNS-GW-Time_onAxis5deg_postRome.txt'
InputTimeList = TableImportCTA_Time(InjectTimeFile)

pointingsFile= '../../gw-follow-up-simulations-side-results/TestPointings/run0017_MergerID000261_GWOptimisation.txt'
FOV= 2.0

PointingPlottingGWCTA(GWFile,name,dirName,pointingsFile,FOV,InputTimeList['Observatory'][j])
