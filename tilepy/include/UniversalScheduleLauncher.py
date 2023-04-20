############################################################################
#       Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler     #
#           LST observation scheduler of GBM alerts and GW events          #
############################################################################


from .UniversalObservationScheduler import PGWinFoV_NObs, PGalinFoV_NObs
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


def GetUniversalSchedule(URL, date, datasetDir, galcatname, outDir, Type, ObsArray):

    '''
    Description: Top level function that is called by the user with specific arguments and creates a folder with the tiling schedules for several telescopes working together and visibility plots.  
    Args:
        URL: the url of the probability fits or  png map
        date: the desired time for scheduling to start 
        datasetDir: Path to the directory containting the datset like the galaxy catalog
        outDir: Path to the output directory where the schedules and plots will eb saved 
        Type: The type of the url given. gw if fits GW map, gbm if fits GBM map and gbmpng if PNG GBM map
        ObsArray: array of strings containing the name of the configuration files of telescopes. 
    '''

    targetType = 'Tiling'
    

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
    parameters = []
    ObservationTime = date
    outputDir = "%s/%s" % (outDir, name)
    galaxies = datasetDir + galcatname

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

    print('parameters', parameters)
    ObsParameters = []

    for j in range(len(parameters)):
        obspar = ObservationParameters()
        obspar.from_configfile(parameters[j])
        ObsParameters.append(obspar)


    if has3D:
        dirName = '%s/PGalinFoV_NObs' % outputDir
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, cat, ObsParameters = PGalinFoV_NObs(filename, ObservationTime, PointingsFile, galaxies, parameters, dirName, ObsArray, ObsParameters)
    else: 
        dirName = '%s/PGWinFoV_NObs' % outputDir
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, ObsParameters = PGWinFoV_NObs(filename, ObservationTime, PointingsFile, parameters, dirName, ObsArray, ObsParameters)
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
            if (len(SuggestedPointings_1) != 0):
                ascii.write(SuggestedPointings_1, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, ObsArray[j]) , overwrite=True, fast_writer=False)
                RankingTimes_SkyMapInput_2D(ObservationTime, prob, ObsParameters[j], targetType, dirName,'%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, ObsArray[j]), ObsArray[j])
                PointingPlotting(prob, ObsParameters[j], name, dirName, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, ObsArray[j]), ObsArray[j], filename)
    else:
        FOLLOWUP = False
        print('No observations are scheduled')

