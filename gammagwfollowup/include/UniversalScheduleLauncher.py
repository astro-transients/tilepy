############################################################################
#       Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler     #
#           LST observation scheduler of GBM alerts and GW events          #
############################################################################


from .UniversalObservationScheduler import PGWinFoV_NObs, PGalinFoV_NObs
from .RankingObservationTimes import RankingTimes, RankingTimes_SkyMapInput_2D
from .PointingPlotting import PointingPlotting
from astropy.coordinates import SkyCoord
from .PointingTools import Tools, LoadGalaxies, getdate, GetGBMMap, GetGWMap, Check2Dor3D
from astropy.io import fits, ascii
import time
import healpy as hp
import numpy as np
from astropy import units as u
import datetime
import os


def GetUniversalSchedule(URL, date, datasetDir, outDir, Type, ObsArray):
    targetType = 'Tiling'
    
    if Type == 'gbm':
        fitsMap, filename = GetGBMMap(URL)

    else: 
        fitsMap, filename = GetGWMap(URL)

    prob, has3D = Check2Dor3D(fitsMap,filename)

    
    print("===========================================================================================")
    name = URL.split('/')[-3]
    PointingsFile = "False"
    parameters = []
    ObservationTime = date
    outputDir = "%s/%s" % (outDir, name)
    galaxies = datasetDir + "/GLADE.txt"

    for i in ObsArray:
        parameters.append("./configs/FollowupParameters_%s.ini" %i)
    print("===========================================================================================")
    print("Starting the IACTs GW - 2D pointing calculation with the following parameters\n")
    print("Filename: ", name)
    print("Date: ", ObservationTime)
    print("Previous pointings: ", PointingsFile)
    print("Parameters: ", parameters)
    print("Dataset: ", datasetDir)
    print("Output: ", outputDir)


    if has3D:
        dirName = '%s/PGallinFoV' % outputDir
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, cat, ObsParameters = PGalinFoV_NObs(filename, ObservationTime, PointingsFile, galaxies, parameters, dirName, ObsArray)
    else: 
        dirName = '%s/PGWinFoV_NObs' % outputDir
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, ObsParameters = PGWinFoV_NObs(filename, ObservationTime, PointingsFile, parameters, dirName, ObsArray)
    if (len(SuggestedPointings) != 0):
        print(SuggestedPointings)
        FOLLOWUP = True
        outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
        ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
        print()
        
        #for obspar in parameters:
        for j in range(len(parameters)):
            obspar1 = ObsParameters[j]
            SuggestedPointings_1 = SuggestedPointings[SuggestedPointings['ObsName'] == obspar1.name]
            print(SuggestedPointings_1)
            ascii.write(SuggestedPointings_1, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, ObsArray[j]) , overwrite=True, fast_writer=False)
            RankingTimes_SkyMapInput_2D(ObservationTime, prob, parameters[j], targetType, dirName,'%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, ObsArray[j]))
            PointingPlotting(prob, parameters[j], name, dirName, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, ObsArray[j]))
    else:
        FOLLOWUP = False
        print('No observations are scheduled')

