##################################################################################################
#        Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler, Mathieu de Bony         #
#                        Observation scheduler of GBM alerts and GW events                       #
##################################################################################################


from .TilingDetermination import PGWinFoV, PGalinFoV
from .TilingDetermination import PGWinFoV_NObs, PGalinFoV_NObs
from .RankingObservationTimes import RankingTimes, RankingTimes_2D
from .PointingPlotting import PointingPlotting
from .PointingTools import GetGBMMap, GetGWMap,getdate, Check2Dor3D, ObservationParameters, GetAreaSkymap5090, GetAreaSkymap5090_Flat
from astropy.io import fits, ascii
from astropy.table import QTable
from astropy import units as u
import os
import json
import numpy as np
import healpy as hp
import ligo.skymap.postprocess as lsp
from astropy.coordinates import SkyCoord

import time
import datetime


def GetSchedule(obspar):
    """
    Top level function that is called by the user with specific arguments and creates a folder 
    with the tiling schedules for a single telescope and visibility plots.  

    :param obspar: the set of parameters needed to launch the tiling scheduler
    :type obspar: class ObservationParameters
    """

    URL = obspar.url

    if obspar.alertType == 'gbmpng':
        fitsMap, filename = GetGBMMap(URL)
        if fitsMap is None and filename is None:
            print('The localization map is not available, returning.')
            return
        name = URL.split('/')[-3]
    elif obspar.alertType == 'gbm':
        fitsMap = fits.open(URL)
        if fitsMap is None:
            print('The localization map is not available, returning.')
            return
        filename = URL
        name = URL.split('all_')[1].split('_v00')[0]
    else:
        fitsMap, filename = GetGWMap(URL)
        name = URL.split('/')[-3]

    
    prob, has3D, origNSIDE = Check2Dor3D(fitsMap, filename, obspar)

    # adapting the resolutions to the one provided in the original map
    if (obspar.HRnside > origNSIDE) :
        print("Reducing HRnside to the value from the original map: NSIDE=",origNSIDE)
        obspar.HRnside = origNSIDE
    if (obspar.reducedNside > obspar.HRnside):
        obspar.reducedNside = obspar.HRnside

    if obspar.locCut != None:
        if(obspar.MO==True):
            area_50, area_90 = GetAreaSkymap5090(filename)
        if(obspar.MO==False):
            area_50, area_90 = GetAreaSkymap5090_Flat(filename)
        if (obspar.locCut== 'loose' and area_90 > 10000) or (obspar.locCut== 'std' and area_50 > 1000) or (obspar.locCut== 'tight' and area_90 > 650) :
            return

    print("===========================================================================================")

    ObservationTime = obspar.obsTime
    outputDir = "%s/%s" % (obspar.outDir, name)

    if has3D:
        dirName = f"{outputDir}/PGallinFoV{obspar.strategy}"
        galaxies = obspar.datasetDir + obspar.galcatName

    else:
        dirName = f"{outputDir}/PGinFoV"

    if not os.path.exists(dirName):
        os.makedirs(dirName)


    if has3D:
        print("===========================================================================================")
        print("Starting the 3D pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", obspar.obsTime)
        print("Previous pointings: ", obspar.pointingsFile)
        print("Catalog: ", galaxies)
        print("Dataset: ", obspar.datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()

        SuggestedPointings, cat = PGalinFoV(filename, galaxies, obspar, dirName)

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            outfilename = '%s/SuggestedPointings_GalProbOptimisation.txt' % dirName
            ascii.write(SuggestedPointings, outfilename,
                        overwrite=True, fast_writer=False)
            print()
            if(obspar.doRank):
                RankingTimes(obspar.obsTime, filename, cat, obspar, obspar.alertType, dirName,
                            '%s/SuggestedPointings_GalProbOptimisation.txt' % dirName, obspar.name)
            if(obspar.doPlot):
                PointingPlotting(prob, obspar, name, dirName,
                                '%s/SuggestedPointings_GalProbOptimisation.txt' % dirName, obspar.name, filename)
        else:
            FOLLOWUP = False
            print('No observations are scheduled')

    else:

        print("===========================================================================================")
        print("Starting the 2D pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", obspar.obsTime)
        print("Previous pointings: ", obspar.pointingsFile)
        print("Dataset: ", obspar.datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()

        SuggestedPointings, t0 = PGWinFoV(filename, obspar, dirName)

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            outfilename = '%s/SuggestedPointings_2DProbOptimisation.txt' % dirName
            ascii.write(SuggestedPointings, outfilename,
                        overwrite=True, fast_writer=False)
            print()
            if(obspar.doRank):
                RankingTimes_2D(obspar.obsTime, prob, obspar, obspar.alertType, dirName,
                                '%s/SuggestedPointings_2DProbOptimisation.txt' % dirName, obspar.name)
            if(obspar.doPlot):
                PointingPlotting(prob, obspar, name, dirName,
                                '%s/SuggestedPointings_2DProbOptimisation.txt' % dirName, obspar.name, filename)
        else:
            FOLLOWUP = False
            print('No observations are scheduled')



def GetUniversalSchedule(obspar):
    '''
    Top level function that is called by the user with specific arguments and creates a folder 
    with the tiling schedules for multiple telescopes/observartories and visibility plots.  

    :param obspar: a list of sets of parameters for each observatory needed to launch the tiling scheduler
    :type obsparameters: list of class ObservationParameters
    '''

    URL = obspar[0].url
    print(URL)

    if obspar[0].alertType == 'gbmpng':
        fitsMap, filename = GetGBMMap(URL)
        if fitsMap is None and filename is None:
            print('The localization map is not available, returning.')
            return
        name = URL.split('/')[-3]
    elif obspar[0].alertType == 'gbm':
        fitsMap = fits.open(URL)
        if fitsMap is None:
            print('The localization map is not available, returning.')
            return
        filename = URL
        name = URL.split('all_')[1].split('_v00')[0]
    else:
        fitsMap, filename = GetGWMap(URL)
        name = URL.split('/')[-3]

    
    prob, has3D, origNSIDE = Check2Dor3D(fitsMap, filename, obspar[0])

    print("===========================================================================================")
    ObservationTime = obspar[0].obsTime
    outputDir =  "%s/%s" % (obspar[0].outDir, name)
    print(obspar[0].datasetDir)
    print(obspar[0].galcatName)
    galaxies = obspar[0].datasetDir + obspar[0].galcatName

    if has3D:
        dirName = '%s/PGalinFoV_NObs' % outputDir
        galaxies = obspar[0].datasetDir + obspar[0].galcatName
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, cat, obspar = PGalinFoV_NObs(
            filename, ObservationTime, obspar[0].pointingsFile, galaxies, obspar, dirName)
    else:
        dirName = '%s/PGWinFoV_NObs' % outputDir
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, obspar = PGWinFoV_NObs(
            filename, ObservationTime, obspar[0].pointingsFile, obspar, dirName)
    if (len(SuggestedPointings) != 0):
        print(SuggestedPointings)
        FOLLOWUP = True
        outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
        ascii.write(SuggestedPointings, outfilename,
                    overwrite=True, fast_writer=False)
        print()

        # for obspar in parameters:
        for j in range(len(obspar)):
            obspar1 = obspar[j]
            SuggestedPointings_1 = SuggestedPointings[SuggestedPointings['ObsName'] == obspar1.name]
            print(SuggestedPointings_1)
            if (len(SuggestedPointings_1) != 0):
                ascii.write(SuggestedPointings_1, '%s/SuggestedPointings_GWOptimisation_%s.txt' %
                            (dirName, obspar[j].name), overwrite=True, fast_writer=False)
                RankingTimes_2D(ObservationTime, prob, obspar[j], obspar[j].alertType, dirName,
                                '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, obspar[j].name), obspar[j].name)
                PointingPlotting(prob, obspar[j], obspar[j].name, dirName, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (
                    dirName, obspar[j].name), obspar[j].name, filename)
        PointingPlotting(prob, obspar[0], "all", dirName,
                             '%s/SuggestedPointings_GWOptimisation.txt' % dirName, "all", filename)
    else:
        FOLLOWUP = False
        print('No observations are scheduled')
