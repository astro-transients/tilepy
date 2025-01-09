#
# Copyright (C) 2016-2024  tilepy developers (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

##################################################################################################
#                        Observation scheduler of GBM alerts and GW events                       #
##################################################################################################


import os

import astropy.units as u
from astropy.io import ascii

from .MapManagement import MapReader
from .MapManagement import SkyMap
from .PointingPlotting import PointingPlotting
from .RankingObservationTimes import RankingTimes, RankingTimes_2D
from .TilingDetermination import PGWinFoV, PGalinFoV
from .TilingDetermination import PGWinFoV_NObs, PGalinFoV_NObs

__all__ = [
    "GetSchedule",
    "GetUniversalSchedule",
]


def GetSchedule(obspar):
    """
    Top level function that is called by the user with specific arguments and creates a folder 
    with the tiling schedules for a single telescope and visibility plots.  

    :param obspar: the set of parameters needed to launch the tiling scheduler
    :type obspar: class ObservationParameters
    """

    raw_map = MapReader(obspar)
    skymap = SkyMap(obspar, raw_map)

    if obspar.locCut90 != None:
        area_90 = skymap.getArea(0.9).to_value(u.deg*u.deg)
        if obspar.locCut90 < area_90:
            print('The 90% area (' + str(area_90) + ' deg^2) is larger than the maximum allowed in the configuration (' + str(obspar.locCut90) + ' deg^2)')
            return

    print("===========================================================================================")

    outputDir = "%s/%s" % (obspar.outDir, raw_map.name_event)

    if skymap.is3D:
        dirName = f"{outputDir}/PGallinFoV{obspar.strategy}"
        galaxies = obspar.datasetDir + obspar.galcatName
    else:
        dirName = f"{outputDir}/PGinFoV"

    if not os.path.exists(dirName):
        os.makedirs(dirName)

    if skymap.is3D:
        print("===========================================================================================")
        print("Starting the 3D pointing calculation with the following parameters\n")
        print("Filename: ", raw_map.name_event)
        print("Date: ", obspar.obsTime)
        print("Previous pointings: ", obspar.pointingsFile)
        print("Catalog: ", galaxies)
        print("Dataset: ", obspar.datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()

        SuggestedPointings, cat = PGalinFoV(skymap, raw_map.name_event, galaxies, obspar, dirName)

        if (len(SuggestedPointings) != 0):
            FOLLOWUP = True
            outfilename = '%s/SuggestedPointings_GalProbOptimisation.txt' % dirName
            ascii.write(SuggestedPointings, outfilename,
                        overwrite=True, fast_writer=False)
            print()
            print(f"Resulting pointings file is {outfilename}")
            if (obspar.doRank):
                RankingTimes(obspar.obsTime, skymap, cat, obspar, dirName,
                             '%s/SuggestedPointings_GalProbOptimisation.txt' % dirName, obspar.obs_name)
            if (obspar.doPlot):
                PointingPlotting(skymap.getMap('prob', obspar.HRnside), obspar, raw_map.name_event, dirName,
                                 '%s/SuggestedPointings_GalProbOptimisation.txt' % dirName, obspar.obs_name, cat)
        else:
            FOLLOWUP = False
            print('No observations are scheduled')

    else:

        print("===========================================================================================")
        print("Starting the 2D pointing calculation with the following parameters\n")
        print("Filename: ", raw_map.name_event)
        print("Date: ", obspar.obsTime)
        print("Previous pointings: ", obspar.pointingsFile)
        print("Dataset: ", obspar.datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()

        SuggestedPointings, t0 = PGWinFoV(skymap, raw_map.name_event, obspar, dirName)

        if (len(SuggestedPointings) != 0):
            gal = []
            FOLLOWUP = True
            outfilename = '%s/SuggestedPointings_2DProbOptimisation.txt' % dirName
            ascii.write(SuggestedPointings, outfilename,
                        overwrite=True, fast_writer=False)
            print()
            print(f"Resulting pointings file is {outfilename}")
            if (obspar.doRank):
                RankingTimes_2D(obspar.obsTime, skymap.getMap('prob', obspar.HRnside), obspar, dirName,
                                '%s/SuggestedPointings_2DProbOptimisation.txt' % dirName, obspar.obs_name)
            if (obspar.doPlot):
                PointingPlotting(skymap.getMap('prob', obspar.HRnside), obspar, raw_map.name_event, dirName,
                                 '%s/SuggestedPointings_2DProbOptimisation.txt' % dirName, obspar.obs_name, gal)
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

    raw_map = MapReader(obspar[0])
    skymap = SkyMap(obspar[0], raw_map)

    print("===========================================================================================")
    ObservationTime = obspar[0].obsTime
    outputDir = "%s/%s" % (obspar[0].outDir, raw_map.name_event)
    galaxies = obspar[0].datasetDir + obspar[0].galcatName
    cat = None

    if skymap.is3D:
        print("===========================================================================================")
        print("Starting the 3D pointing calculation with the following parameters\n")
        print("Filename: ", raw_map.name_event)
        print("Date: ", obspar[0].obsTime)
        print("Catalog: ", obspar[0].galcatName)
        print("Dataset: ", obspar[0].datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()
        dirName = '%s/PGalinFoV_NObs' % outputDir
        galaxies = obspar[0].datasetDir + obspar[0].galcatName
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, cat, obspar = PGalinFoV_NObs(
            skymap, raw_map.name_event, ObservationTime, obspar[0].pointingsFile, galaxies, obspar, dirName)
    else:
        print("===========================================================================================")
        print("Starting the 2D pointing calculation with the following parameters\n")
        print("Filename: ", raw_map.name_event)
        print("Date: ", obspar[0].obsTime)
        print("Dataset: ", obspar[0].datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()
        dirName = '%s/PGWinFoV_NObs' % outputDir
        if not os.path.exists(dirName):
            os.makedirs(dirName)
        SuggestedPointings, obspar = PGWinFoV_NObs(
            skymap, raw_map.name_event, ObservationTime, obspar[0].pointingsFile, obspar, dirName)
    if (len(SuggestedPointings) != 0):
        print(SuggestedPointings)
        FOLLOWUP = True
        outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
        ascii.write(SuggestedPointings, outfilename,
                    overwrite=True, fast_writer=False)
        print()
        print(f"Resulting pointings file is {outfilename}")
        # for obspar in parameters:
        for j in range(len(obspar)):
            obspar1 = obspar[j]
            SuggestedPointings_1 = SuggestedPointings[SuggestedPointings['ObsName'] == obspar1.obs_name]
            print(SuggestedPointings_1)
            if (len(SuggestedPointings_1) != 0):
                ascii.write(SuggestedPointings_1, '%s/SuggestedPointings_GWOptimisation_%s.txt' %
                            (dirName, obspar[j].obs_name), overwrite=True, fast_writer=False)
                RankingTimes_2D(ObservationTime, skymap.getMap('prob', obspar[j].HRnside), obspar[j], dirName,
                                '%s/SuggestedPointings_GWOptimisation_%s.txt' % (dirName, obspar[j].obs_name),
                                obspar[j].obs_name)
                PointingPlotting(skymap.getMap('prob', obspar[j].HRnside), obspar[j], obspar[j].obs_name, dirName,
                                 '%s/SuggestedPointings_GWOptimisation_%s.txt' % (
                                     dirName, obspar[j].obs_name), obspar[j].obs_name, cat)
        PointingPlotting(skymap.getMap('prob', obspar[j].HRnside), obspar[0], "all", dirName, '%s/SuggestedPointings_GWOptimisation.txt' % dirName, "all",
                         cat)
    else:
        FOLLOWUP = False
        print('No observations are scheduled')


def GetObservationCoordinates(obspar):
    """
    Top level function that is called by the user with specific arguments and creates a folder 
    with potentially observable coordinates (generally speaking, not linked to any telescope)  

    :param obspar: the set of parameters needed to launch the tiling scheduler
    :type obspar: class ObservationParameters
    """

    raw_map = MapReader(obspar)
    skymap = SkyMap(obspar, raw_map)

    if obspar.locCut90 != None:
        area_90 = skymap.getArea(0.9).to_value(u.deg*u.deg)
        if obspar.locCut90 < area_90:
            print('The 90% area (' + str(area_90) + ' deg^2) is larger than the maximum allowed in the configuration (' + str(obspar.locCut90) + ' deg^2)')
            return

    print("===========================================================================================")

    outputDir = "%s/%s" % (obspar.outDir, raw_map.name_event)

    if skymap.is3D:
        dirName = f"{outputDir}/PGalRank{obspar.strategy}"
        galaxies = obspar.datasetDir + obspar.galcatName
    else:
        dirName = f"{outputDir}/PGRank"

    if not os.path.exists(dirName):
        os.makedirs(dirName)

    if skymap.is3D:
        print("===========================================================================================")
        print("Starting the identification of the most probable coordinates the following parameters\n")
        print("Filename: ", raw_map.name_event)
        print("Catalog: ", galaxies)
        print("Dataset: ", obspar.datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()

        SuggestedCoordinate, cat = ObtainMostProbableCoordinatesGal(skymap, raw_map.name_event, galaxies, obspar, dirName)

        if (len(SuggestedCoordinates) != 0):
            outfilename = '%s/SuggestedCoordinates_GalProbOptimisation.txt' % dirName
            ascii.write(SuggestedCoordinates, outfilename,
                        overwrite=True, fast_writer=False)
            print()
            print(f"Resulting pointings file is {outfilename}")

            if (obspar.doPlot):
                PointingPlotting(skymap.getMap('prob', obspar.HRnside), obspar, raw_map.name_event, dirName,
                                 '%s/SuggestedCoordinates_GalProbOptimisation.txt' % dirName, obspar.obs_name, cat)
        else:
            print('No observations are scheduled')

    else:

        print("===========================================================================================")
        print("Starting the 2D pointing calculation with the following parameters\n")
        print("Filename: ", raw_map.name_event)
        print("Dataset: ", obspar.datasetDir)
        print("Output: ", outputDir)
        print("===========================================================================================")
        print()

        SuggestedCoordinates = ObtainMostProbableCoordinates(skymap, raw_map.name_event, obspar, dirName)

        if (len(SuggestedCoordinates) != 0):
            gal = []
            outfilename = '%s/SuggestedCoordinates_GWProbOptimisation.txt' % dirName
            ascii.write(SuggestedCoordinates, outfilename,
                        overwrite=True, fast_writer=False)
            print()
            print(f"Resulting pointings file is {outfilename}")

            if (obspar.doPlot):
                PointingPlotting(skymap.getMap('prob', obspar.HRnside), obspar, raw_map.name_event, dirName,
                                 '%s/SuggestedCoordinates_GWProbOptimisation.txt' % dirName, obspar.obs_name, gal)
        else:
            print('No observations are scheduled')
