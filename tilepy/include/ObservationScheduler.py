############################################################################
#       Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler     #
#           LST observation scheduler of GBM alerts and GW events          #
############################################################################


from .TilingDetermination import PGWinFoV, PGalinFoV
from .RankingObservationTimes import RankingTimes, RankingTimes_SkyMapInput_2D
from .PointingPlotting import PointingPlotting
from astropy.coordinates import SkyCoord
from .PointingTools import Tools, LoadGalaxies, getdate, GetGBMMap, GetGWMap, Check2Dor3D, ObservationParameters
from astropy.io import fits, ascii
import time
import healpy as hp
import numpy as np
from astropy import units as u
import datetime
import os
import json



def GetSchedule_confile(URL, date,datasetDir,outDir,cfgFile, Type):

    if Type == 'gbmpng':
        targetType = 'GBM_Pointing'
        fitsMap, filename = GetGBMMap(URL)
        name = URL.split('/')[-3]
    elif Type == 'gbm':
        targetType = 'GBM_Pointing'
        fitsMap = fits.open(URL)
        filename = URL
        name = URL.split('all_')[1].split('_v00')[0]
    else: 
        targetType = 'GW_Pointing'
        fitsMap, filename = GetGWMap(URL)
        name = URL.split('/')[-3]

    prob, has3D = Check2Dor3D(fitsMap,filename)

    
    print("===========================================================================================")
    PointingsFile = "False"
    galaxies = datasetDir + "/GLADE.txt"
    #cfgFile = "./configs/FollowupParameters.ini"

    if has3D:

        ObservationTime = date
        outputDir = "%s/%s" % (outDir, name)
        dirName = '%s/PGallinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Starting the LST GW - 3D pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galaxies)
        print("Config parameters: ", cfgFile)
        print("Dataset: ", datasetDir)
        print("Output: ", outputDir)
        
        obspar = ObservationParameters()
        obspar.from_configfile(cfgFile)

        SuggestedPointings, cat = PGalinFoV(filename, ObservationTime, PointingsFile, galaxies, obspar, dirName)

        print(SuggestedPointings)
        print("===========================================================================================")
        print()

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
            ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
            print()
            RankingTimes(ObservationTime, filename, cat, obspar, targetType, dirName,
                         '%s/SuggestedPointings_GWOptimisation.txt' % dirName, obspar.name)
            PointingPlotting(prob, obspar, name, dirName,
                             '%s/SuggestedPointings_GWOptimisation.txt' % dirName, obspar.name, filename)
        else:
            FOLLOWUP = False
            print('No observations are scheduled')

    else:

        ObservationTime = date
        outputDir = "%s/%s" % (outDir, name)
        dirName = '%s/PGWinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Starting the LST GW - 2D pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galaxies)
        print("Parameters: ", cfgFile)
        print("Dataset: ", datasetDir)
        print("Output: ", outputDir)
        
        obspar = ObservationParameters()
        obspar.from_configfile(cfgFile)

        SuggestedPointings, t0 = PGWinFoV(filename, ObservationTime, PointingsFile, obspar, dirName)

        print(SuggestedPointings)
        print("===========================================================================================")
        print()

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
            ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
            print()
            cat = LoadGalaxies(galaxies)
            RankingTimes(ObservationTime, filename, cat, obspar, targetType, dirName,
                         '%s/SuggestedPointings_GWOptimisation.txt' % dirName, obspar.name)
            PointingPlotting(prob, obspar, name, dirName, '%s/SuggestedPointings_GWOptimisation.txt' % dirName, obspar.name, filename)
        else:
            FOLLOWUP = False
            print('No observations are scheduled')


def GetSchedule_funcarg(Type, URL, date,datasetDir,outDir, name, Lat, Lon, Height, gSunDown, HorizonSun, gMoonDown,
                 HorizonMoon, gMoonGrey, gMoonPhase, MoonSourceSeparation,
                 MaxMoonSourceSeparation, max_zenith, FOV, MaxRuns, MaxNights,
                 Duration, MinDuration, UseGreytime, MinSlewing, online,
                 MinimumProbCutForCatalogue, MinProbCut, doplot, SecondRound ,
                 FulFillReq_Percentage, PercentCoverage, ReducedNside, HRnside,
                 Mangrove):

    if Type == 'gbmpng':
        targetType = 'GBM_Pointing'
        fitsMap, filename = GetGBMMap(URL)
        name = URL.split('/')[-3]
    elif Type == 'gbm':
        targetType = 'GBM_Pointing'
        fitsMap = fits.open(URL)
        filename = URL
        name = URL.split('all_')[1].split('_v00')[0]
    else: 
        targetType = 'GW_Pointing'
        fitsMap, filename = GetGWMap(URL)
        name = URL.split('/')[-3]

    prob, has3D = Check2Dor3D(fitsMap,filename) 

    print("===========================================================================================")
    PointingsFile = "False"
    galaxies = datasetDir + "/GLADE.txt"
    #cfgFile = "./configs/FollowupParameters.ini"

    obspar = ObservationParameters()
    obspar.from_args(name, Lat, Lon, Height, gSunDown, HorizonSun, gMoonDown,
                 HorizonMoon, gMoonGrey, gMoonPhase, MoonSourceSeparation,
                 MaxMoonSourceSeparation, max_zenith, FOV, MaxRuns, MaxNights,
                 Duration, MinDuration, UseGreytime, MinSlewing, online,
                 MinimumProbCutForCatalogue, MinProbCut, doplot, SecondRound ,
                 FulFillReq_Percentage, PercentCoverage, ReducedNside, HRnside,
                 Mangrove)

    if has3D:

        ObservationTime = date
        outputDir = "%s/%s" % (outDir, name)
        dirName = '%s/PGallinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galaxies)
        print("Dataset: ", datasetDir)
        print("Output: ", outputDir)
        

        SuggestedPointings, cat = PGalinFoV(filename, ObservationTime, PointingsFile, galaxies, obspar, dirName)

        print(SuggestedPointings)
        print("===========================================================================================")
        print()

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            df = SuggestedPointings.to_pandas()
            table_dict = df.to_dict()
            SuggestedPointings_AstroCOLIBRI = json.dumps(table_dict)
            print()
            return SuggestedPointings_AstroCOLIBRI
        else:
            FOLLOWUP = False
            print('No observations are scheduled')
            return None

    else:

        ObservationTime = date
        outputDir = "%s/%s" % (outDir, name)
        dirName = '%s/PGWinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galaxies)
        print("Dataset: ", datasetDir)
        print("Output: ", outputDir)
        

        SuggestedPointings, t0 = PGWinFoV(filename, ObservationTime, PointingsFile, obspar, dirName)

        print(SuggestedPointings)
        print("===========================================================================================")
        print()

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            df = SuggestedPointings.to_pandas()
            table_dict = df.to_dict()
            SuggestedPointings_AstroCOLIBRI = json.dumps(table_dict)
            print()
            return SuggestedPointings_AstroCOLIBRI
        else:
            FOLLOWUP = False
            print('No observations are scheduled')
            return None



