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


    name = URL.split('/')[-3]
    print("===========================================================================================")
    PointingsFile = "False"
    galaxies = datasetDir + "/GLADE.txt"
    parameters = []

    j = 0
    for i in ObsArray:
        parameters.append("./configs/FollowupParameters_%s.ini" %i)
    print(parameters)

    if has3D:
        print("Will do that later")
    else: 
        print("Will do that now")
        #PGWinFoV_NObs(parameters,parameters,parameters,parameters,ObsArray)



    #SuggestedPointings, t0 = PGWinFoV(filename, ObservationTime, PointingsFile, parameters, dirName)