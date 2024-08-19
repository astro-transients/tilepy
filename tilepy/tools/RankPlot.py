
from tilepy.include.PointingTools import ObservationParameters
import time
import argparse
import os
from tilepy.include.RankingObservationTimes import RankingTimes_2D
from tilepy.include.PointingPlotting import PointingPlotting
from tilepy.include.PointingTools import Check2Dor3D, ObservationParameters, LoadHealpixMap, LoadGalaxies, LoadGalaxies_SteMgal, CorrelateGalaxies_LVC, CorrelateGalaxies_LVC_SteMass
from tilepy.include.MapReader import GetSkymap
import os
import time
import healpy as hp
import numpy as np
from astropy.io import fits


def RankPlot(skymap, alertType, obsTime, configDir, datasetDir, outDir, galcatName, pointingsFile, locCut, ObsArray):

    parameters = []

    for i in ObsArray:
        parameters.append("%s/FollowupParameters_%s.ini" %(configDir ,i))
    print("===========================================================================================")
    print('parameters', parameters)
    obsparameters = []



    for j in range(len(parameters)):
        obspar = ObservationParameters()
        obspar.add_parsed_args(skymap, obsTime, datasetDir, galcatName, outDir, pointingsFile, alertType, locCut)
        obspar.from_configfile(parameters[j])
        print(obspar)
        obsparameters.append(obspar)
    
    
    galaxies = datasetDir + galcatName

    skymap = obsparameters[0].skymap
    print(skymap)

    fitsMap, filename, name = GetSkymap(skymap) 
    #if obsparameters[0].alertType == 'gbmpng':
    #    fitsMap, filename = GetGBMMap(skymap)
    #    if fitsMap is None and filename is None:
    #        print('The localization map is not available, returning.')
    #    name = URL.split('/')[-3]
    #elif obsparameters[0].alertType == 'gbm':
    #    fitsMap = fits.open(URL)
    #    if fitsMap is None:
    #        print('The localization map is not available, returning.')
    #    filename = URL
    #    name = URL.split('all_')[1].split('_v00')[0]
    #else:
    #    fitsMap, filename = GetGWMap(URL)
    #    name = URL.split('/')[-3]

    prob, has3D, origNSIDE = Check2Dor3D(fitsMap, filename, obspar)
    print("Original NSIDE", origNSIDE, has3D)
    #has3D = True

    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    # load galaxy catalogue
    if not obspar.mangrove:
        cat = LoadGalaxies(galaxies)
    else:
        cat = LoadGalaxies_SteMgal(galaxies)


    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(
            prob, distmu, distsigma, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(
            prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)


    print("===========================================================================================")
    ObservationTime = obsparameters[0].obsTime
    outputDir =  "%s/%s" % (obsparameters[0].outDir, name)
    print(obsparameters[0].datasetDir)
    print(obsparameters[0].galcatName)
    print(obsparameters[0].name)
    
    
    if has3D:
        dirName = '%s/PGalinFoV_NObs' % outputDir
    else:
        dirName = '%s/PGWinFoV_NObs' % outputDir


    for j in range(len(obsparameters)):
        #RankingTimes_2D(ObservationTime, prob, obsparameters[j], obsparameters[j].alertType, dirName,
        #                '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, obsparameters[j].name), obsparameters[j].name)
        PointingPlotting(prob, obsparameters[j], obsparameters[j].name, dirName, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (
            dirName, obsparameters[j].name), obsparameters[j].name, tGals0)
        
    PointingPlotting(prob, obsparameters[0], "all", dirName,
                            '%s/SuggestedPointings_GWOptimisation.txt' % dirName, "all", tGals0)

