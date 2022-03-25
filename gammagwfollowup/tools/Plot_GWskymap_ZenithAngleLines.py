import sys
import argparse
import os
from gammagwfollowup.include.PointingPlotting import PointingPlottingGW_ZenithSteps
from gammagwfollowup.include.PointingTools import TableImportCTA, TableImportCTA_Time


parser = argparse.ArgumentParser(description='Start the CTA GW pointing simulation and analysis (PgwinFov)')
parser.add_argument('index', metavar='GWmap',
                    help='the FITS file with the GW localization, e.g. Bayestar.fits.gz')
parser.add_argument('-i',metavar = 'input path', help='Path to the input datasets')
args = parser.parse_args()
j=int(args.index)
datasetDir = args.i

# GW file
InputFileName = datasetDir +'/BNS-GW_onAxis5deg.txt'
InputList = TableImportCTA(InputFileName)

GWFile = datasetDir +"/skymaps/" + InputList['run'][j] + '_' + InputList['MergerID'][j] + "_skymap.fits.gz"
name =  InputList['run'][j] + '_' + InputList['MergerID'][j]

dirName='../output/PGWonFoV'

InjectTimeFile = datasetDir + '/BNS-GW-Time_onAxis5deg_postRome.txt'
InputTimeList = TableImportCTA_Time(InjectTimeFile)

FOV= 2.0

PointingPlottingGW_ZenithSteps(GWFile,name,dirName,FOV,InputTimeList[j])
