############################################################################
#       Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler     #
#           LST observation scheduler of GBM alerts and GW events          #
############################################################################


from .UniversalObservationScheduler import PGWinFoV_NObs
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
    dirName = '%s/PGWinFoV_NObs' % outputDir

    if not os.path.exists(dirName):
        os.makedirs(dirName)

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
        print("Will do that later")
    else: 
        print("Will do that now")
        SuggestedPointings = PGWinFoV_NObs(filename, ObservationTime, PointingsFile, parameters, dirName, ObsArray)

    if (len(SuggestedPointings) != 0):
        print(SuggestedPointings)
        FOLLOWUP = True
        outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
        ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
        print()
        #for obspar in parameters:
            #RankingTimes_SkyMapInput_2D(ObservationTime, prob, obspar, targetType, dirName,'%s/SuggestedPointings_GWOptimisation.txt' % dirName)
            #PointingPlotting(prob, obspar, name, dirName, '%s/SuggestedPointings_GWOptimisation.txt' % dirName)
    else:
        FOLLOWUP = False
        print('No observations are scheduled')
