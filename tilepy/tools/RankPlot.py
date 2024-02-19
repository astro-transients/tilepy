
from tilepy.include.PointingTools import ObservationParameters
import time
import argparse
import os
from tilepy.include.RankingObservationTimes import RankingTimes_2D
from tilepy.include.PointingPlotting import PointingPlotting
from tilepy.include.PointingTools import GetGBMMap, GetGWMap, Check2Dor3D, ObservationParameters
import os
import time


def RankPlot(url, alertType, obsTime, configDir, datasetDir, outDir, galcatName, pointingsFile, locCut, ObsArray):

    parameters = []

    for i in ObsArray:
        parameters.append("%s/FollowupParameters_%s.ini" %(configDir ,i))
    print("===========================================================================================")
    print('parameters', parameters)
    obsparameters = []

    for j in range(len(parameters)):
        obspar = ObservationParameters()
        obspar.add_parsed_args(url, obsTime, datasetDir, galcatName, outDir, pointingsFile, alertType, locCut)
        obspar.from_configfile(parameters[j])
        print(obspar)
        obsparameters.append(obspar)


    URL = obsparameters[0].url
    print(URL)

    if obsparameters[0].alertType == 'gbmpng':
        fitsMap, filename = GetGBMMap(URL)
        if fitsMap is None and filename is None:
            print('The localization map is not available, returning.')


        name = URL.split('/')[-3]
    elif obsparameters[0].alertType == 'gbm':
        fitsMap = fits.open(URL)
        if fitsMap is None:
            print('The localization map is not available, returning.')


        filename = URL
        name = URL.split('all_')[1].split('_v00')[0]
    else:
        fitsMap, filename = GetGWMap(URL)
        name = URL.split('/')[-3]


    prob, has3D, origNSIDE = Check2Dor3D(fitsMap, filename, obsparameters[0])


    print("===========================================================================================")
    ObservationTime = obsparameters[0].obsTime
    outputDir =  "%s/%s" % (obsparameters[0].outDir, name)
    print(obsparameters[0].datasetDir)
    print(obsparameters[0].galcatName)
    print(obsparameters[0].name)
    galaxies = obsparameters[0].datasetDir + obsparameters[0].galcatName



    if has3D:
        dirName = '%s/PGalinFoV_NObs' % outputDir
    else:
        dirName = '%s/PGWinFoV_NObs' % outputDir


    for j in range(len(obsparameters)):
        RankingTimes_2D(ObservationTime, prob, obsparameters[j], obsparameters[j].alertType, dirName,
                        '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, obsparameters[j].name), obsparameters[j].name)
        PointingPlotting(prob, obsparameters[j], obsparameters[j].name, dirName, '%s/SuggestedPointings_GWOptimisation_%s.txt' % (
            dirName, obsparameters[j].name), obsparameters[j].name, filename, galaxies)
        
    PointingPlotting(prob, obsparameters[0], "all", dirName,
                            '%s/SuggestedPointings_GWOptimisation.txt' % dirName, "all", filename, galaxies)

