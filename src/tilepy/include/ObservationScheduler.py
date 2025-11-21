# Copyright (C) 2016-2025  tilepy developers
# (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

##################################################################################################
#                        Observation scheduler of GBM alerts and GW events                       #
##################################################################################################


import logging

import astropy.units as u
from astropy.io import ascii
from astropy.table import Table

from pathlib import Path

from .MapManagement import SkyMap, create_map_reader
from .PointingPlotting import PlotAccRegion, PointingPlotting
from .RankingObservationTimes import (
    PlotAccRegionTimePix,
    PlotAccRegionTimeRadec,
    Ranking_Space,
    Ranking_Space_AI,
    RankingTimes,
    RankingTimes_2D,
)
from .TilingDetermination import (
    GetBestTiles2D,
    GetBestTiles3D,
    PGalinFoV,
    PGalinFoV_NObs,
    PGalinFoV_Space_NObs,
    PGWinFoV,
    PGWinFoV_NObs,
    PGWinFoV_Space_NObs,
)

__all__ = [
    "GetSchedule",
    "GetUniversalSchedule",
]

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


def GetSchedule(obspar):
    """
    Generates a tiling schedule and visibility plots for a single telescope.

    This top level function takes the observation parameters, computes the optimal tiling
    and creates a dedicated output folder containing the scheduled pointings and optional plots.

    Parameters
    ----------
    obspar : ObservationParameters
        The set of parameters needed to launch the tiling scheduler.

    Returns
    -------
    None

    Notes
    -----
    If the 90% localization area is larger than the configured maximum (`locCut`), no schedule is generated.

    Examples
    --------
    >>> GetSchedule(obspar)

    """

    raw_map = create_map_reader(obspar)
    skymap = SkyMap(obspar, raw_map)

    area_90 = skymap.getArea(0.9).to_value(u.deg * u.deg)
    area_50 = skymap.getArea(0.5).to_value(u.deg * u.deg)
    area_percentage = skymap.getArea(obspar.percentageMOC).to_value(u.deg * u.deg)

    if obspar.locCut is not None:
        # FIXME : a repetitive calculation
        # area_90 = skymap.getArea(0.9).to_value(u.deg * u.deg)
        if obspar.locCut < area_percentage:
            logger.info(
                f"The {obspar.percentageMOC * 100:.1f}% area ({area_percentage:.2f} deg^2) is larger than the maximum allowed in the configuration ({str(obspar.locCut)} deg^2)"
            )
            return

    # FIXME : Not necessary but could ovoid to have along "==" in the code
    # print("=" * 91)
    logger.info(
        "==========================================================================================="
    )

    outputDir = Path(f"{obspar.outDir}/{raw_map.name_event}")

    if skymap.is3D:
        dirName = outputDir / f"PGallinFoV{obspar.strategy}"
        galaxies = Path(f"{obspar.datasetDir}/{obspar.galcatName}")
    else:
        dirName = outputDir / "PGinFoV"

    if not dirName.exists():
        dirName.mkdir(parents=True)

    if skymap.is3D:
        logger.info(
            "==========================================================================================="
        )
        logger.info(
            "Starting the 3D pointing calculation with the following parameters\n"
        )
        logger.info(f"Filename:  {raw_map.name_event}")
        logger.info(f"Date: {obspar.obsTime}")
        logger.info(f"Previous pointings: {obspar.pointingsFile}")
        logger.info(f"Catalog: {str(galaxies)}")
        logger.info(f"Dataset: {obspar.datasetDir}")
        logger.info(f"Output: {str(outputDir)}")
        logger.info(f"90% area = {area_90}. 50% area = {area_50}")
        logger.info(
            "===========================================================================================\n"
        )

        SuggestedPointings, cat = PGalinFoV(
            skymap, raw_map.name_event, str(galaxies), obspar, str(dirName)
        )

        if len(SuggestedPointings) != 0:
            outfilename = f"{str(dirName)}/SuggestedPointings_GalProbOptimisation.txt"
            ascii.write(
                SuggestedPointings, outfilename, overwrite=True, fast_writer=False
            )
            logger.info(f"\nResulting pointings file is {outfilename}")
            if obspar.doRank:
                RankingTimes(
                    obspar,
                    skymap,
                    cat,
                    str(dirName),
                    f"{str(dirName)}/SuggestedPointings_GalProbOptimisation.txt",
                )
            if obspar.doPlot:
                PointingPlotting(
                    skymap.getMap("prob", obspar.HRnside),
                    obspar,
                    raw_map.name_event,
                    str(dirName),
                    f"{str(dirName)}/SuggestedPointings_GalProbOptimisation.txt",
                    obspar.obs_name,
                    cat,
                )
        else:
            logger.info("No observations are scheduled")

    else:
        logger.info(
            "==========================================================================================="
        )
        logger.info(
            "Starting the 2D pointing calculation with the following parameters\n"
        )
        logger.info(f"Filename:  {raw_map.name_event}")
        logger.info(f"Date: {obspar.obsTime}")
        logger.info(f"Previous pointings: {obspar.pointingsFile}")
        logger.info(f"Dataset: {obspar.datasetDir}")
        logger.info(f"Output: {str(outputDir)}")
        logger.info(f"90% area = {area_90}. 50% area = {area_50}")
        logger.info(
            "===========================================================================================\n"
        )

        SuggestedPointings, t0 = PGWinFoV(skymap, raw_map.name_event, obspar, str(dirName))

        if len(SuggestedPointings) != 0:
            gal = []
            outfilename = f"{str(dirName)}/SuggestedPointings_2DProbOptimisation.txt"
            ascii.write(
                SuggestedPointings, outfilename, overwrite=True, fast_writer=False
            )
            logger.info(f"Resulting pointings file is {outfilename}")
            if obspar.doRank:
                RankingTimes_2D(
                    obspar,
                    skymap.getMap("prob", obspar.HRnside),
                    str(dirName),
                    f"{str(dirName)}/SuggestedPointings_2DProbOptimisation.txt",
                )
            if obspar.doPlot:
                PointingPlotting(
                    skymap.getMap("prob", obspar.HRnside),
                    obspar,
                    raw_map.name_event,
                    str(dirName),
                    f"{str(dirName)}/SuggestedPointings_2DProbOptimisation.txt",
                    obspar.obs_name,
                    gal,
                )
        else:
            logger.info("No observations are scheduled")


def GetUniversalSchedule(obspar):
    """
    Generates tiling schedules and visibility plots for multiple telescopes/observatories.

    This top-level function takes as input a list of observation parameter sets (one per observatory),
    computes optimal tilings according to the skymap, and saves the results in output folders.
    Optionally, it produces ranked pointings and plots if requested in the input parameters.

    Parameters
    ----------
    obspar : list of ObservationParameters
        List of parameter objects, one for each observatory, defining the scheduling configuration.

    Returns
    -------
    None

    Notes
    -----
    The function creates necessary output folders and writes results to disk. If no suitable
    schedule can be generated for any observatory, no output file is created for that observatory.

    Examples
    --------
    >>> GetUniversalSchedule([obspar1, obspar2])

    """

    raw_map = create_map_reader(obspar[0])
    skymap = SkyMap(obspar[0], raw_map)

    area_90 = skymap.getArea(0.9).to_value(u.deg * u.deg)
    area_50 = skymap.getArea(0.5).to_value(u.deg * u.deg)

    base = obspar[0].base

    ObservationTime = obspar[0].obsTime

    outputDir = Path(f"{obspar[0].outDir}/{raw_map.name_event}")

    if skymap.is3D:
        galaxies = Path(f"{obspar[0].datasetDir}/{obspar[0].galcatName}")

    cat = None

    if base == "grid":
        if skymap.is3D:
            dirName = outputDir / "GetBestTiles3D"
            if not dirName.exists():
                dirName.mkdir(parents=True)
            SuggestedPointings = GetBestTiles3D(
                skymap,
                raw_map.name_event,
                obspar[0].pointingsFile,
                str(galaxies),
                obspar,
                str(dirName),
            )
            # print(SuggestedPointings)
        else:
            dirName = outputDir / "GetBestTiles2D"
            if not dirName.exists():
                dirName.mkdir(parents=True)
            SuggestedPointings = GetBestTiles2D(
                skymap, raw_map.name_event, obspar[0].pointingsFile, obspar, str(dirName)
            )

    # SPACE
    elif base == "space":
        if skymap.is3D:
            logger.info(
                "==========================================================================================="
            )
            logger.info(
                "Starting the 3D pointing calculation with the following parameters\n"
            )
            logger.info(f"Filename:  {raw_map.name_event}")
            logger.info(f"Date: {obspar[0].obsTime}")
            logger.info(f"Previous pointings: {obspar[0].pointingsFile}")
            logger.info(f"Catalog: {str(galaxies)}")
            logger.info(f"Dataset: {obspar[0].datasetDir}")
            logger.info(f"Output: {str(outputDir)}")
            logger.info(f"90% area = {area_90}. 50% area = {area_50}")
            logger.info(
                "===========================================================================================\n"
            )
            dirName = outputDir / "PGalinFoV_NObs_Space"
            if not dirName.exists():
                dirName.mkdir(parents=True)
            SuggestedPointings, result = PGalinFoV_Space_NObs(
                skymap,
                raw_map.name_event,
                ObservationTime,
                obspar[0].pointingsFile,
                str(galaxies),
                obspar,
                str(dirName),
            )
            if obspar[0].doPlot and len(result["first_values1"]) > 0:
                PlotAccRegion(
                    skymap,
                    str(dirName),
                    obspar[0].reducedNside,
                    result["Occultedpixels"],
                    result["first_values"],
                )
            if obspar[0].doRank:
                PlotAccRegionTimePix(
                    str(dirName),
                    result["AvailablePixPerTime"],
                    result["ProbaTime"],
                    result["TestTime"],
                )
                PlotAccRegionTimeRadec(
                    str(dirName),
                    result["AvailablePixPerTime"],
                    result["ProbaTime"],
                    result["TestTime"],
                    obspar[0].reducedNside,
                )

        else:
            logger.info(
                "==========================================================================================="
            )
            logger.info(
                "Starting the 2D pointing calculation with the following parameters\n"
            )
            logger.info(f"Filename:  {raw_map.name_event}")
            logger.info(f"Date: {obspar[0].obsTime}")
            logger.info(f"Previous pointings: {obspar[0].pointingsFile}")
            logger.info(f"Dataset: {obspar[0].datasetDir}")
            logger.info(f"Output: {str(outputDir)}")
            logger.info(f"90% area = {area_90}. 50% area = {area_50}")
            logger.info(
                "===========================================================================================\n"
            )
            dirName = outputDir / "PGWinFoV_Space_NObs"
            if not dirName.exists():
                dirName.mkdir(parents=True)
            SuggestedPointings, result = PGWinFoV_Space_NObs(
                skymap,
                raw_map.name_event,
                ObservationTime,
                obspar[0].pointingsFile,
                obspar,
                str(dirName),
            )
            if obspar[0].doPlot and len(result["first_values1"]) > 0:
                PlotAccRegion(
                    skymap,
                    str(dirName),
                    obspar[0].reducedNside,
                    result["Occultedpixels"],
                    result["first_values"],
                )
            if obspar[0].doRank:
                PlotAccRegionTimePix(
                    str(dirName),
                    result["AvailablePixPerTime"],
                    result["ProbaTime"],
                    result["TestTime"],
                )
                PlotAccRegionTimeRadec(
                    str(dirName),
                    result["AvailablePixPerTime"],
                    result["ProbaTime"],
                    result["TestTime"],
                    obspar[0].reducedNside,
                )
                # PlotAccRegionTimePix(dirName, AvailablePixPerTime, ProbaTime, TestTime)
                # PlotAccRegionTimeRadec(
                #    dirName, AvailablePixPerTime, ProbaTime, result["TestTime"], reducedNside
                # )

    else:
        # GROUND
        if skymap.is3D:
            logger.info(
                "==========================================================================================="
            )
            logger.info(
                "Starting the 3D pointing calculation with the following parameters\n"
            )
            logger.info(f"Filename:  {raw_map.name_event}")
            logger.info(f"Date: {obspar[0].obsTime}")
            logger.info(f"Previous pointings: {obspar[0].pointingsFile}")
            logger.info(f"Catalog: {str(galaxies)}")
            logger.info(f"Dataset: {obspar[0].datasetDir}")
            logger.info(f"Output: {str(outputDir)}")
            logger.info(f"90% area = {area_90}. 50% area = {area_50}")
            logger.info(
                "===========================================================================================\n"
            )
            dirName = outputDir / "PGalinFoV_NObs"
            if not dirName.exists():
                dirName.mkdir(parents=True)
            SuggestedPointings, cat, obspar = PGalinFoV_NObs(
                skymap,
                raw_map.name_event,
                ObservationTime,
                obspar[0].pointingsFile,
                str(galaxies),
                obspar,
                str(dirName),
            )

        else:
            logger.info(
                "==========================================================================================="
            )
            logger.info(
                "Starting the 2D pointing calculation with the following parameters\n"
            )
            logger.info(f"Filename:  {raw_map.name_event}")
            logger.info(f"Date: {obspar[0].obsTime}")
            logger.info(f"Previous pointings: {obspar[0].pointingsFile}")
            logger.info(f"Dataset: {obspar[0].datasetDir}")
            logger.info(f"Output: {str(outputDir)}")
            logger.info(f"90% area = {area_90}. 50% area = {area_50}")
            logger.info(
                "===========================================================================================\n"
            )
            dirName = outputDir / "PGWinFoV_NObs"
            if not dirName.exists():
                dirName.mkdir(parents=True)
            SuggestedPointings, obspar = PGWinFoV_NObs(
                skymap,
                raw_map.name_event,
                ObservationTime,
                obspar[0].pointingsFile,
                obspar,
                str(dirName),
            )

    if len(SuggestedPointings) != 0:
        logger.info(f"Suggested pointings: {SuggestedPointings}")
        outfilename = str( dirName / "SuggestedPointings_GWOptimisation.txt")
        ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
        logger.info(f"Resulting pointings file is {outfilename}")

        if base in ["space"]:
            for j, obs in enumerate(obspar):
                obspar1 = obspar[j]
                SuggestedPointings_1 = SuggestedPointings[
                    SuggestedPointings["ObsName"] == obspar[j].obs_name
                ]
                logger.info(
                    f"Suggested pointings for {obspar[j].obs_name}: {SuggestedPointings_1}"
                )
                if base == "space":
                    time_table = Table(
                        [result["SatTimes"], result["saa"]], names=("SatTimes", "SAA")
                    )
                    ascii.write(
                        time_table,
                        f"{str(dirName)}/SAA_Times_{obspar[j].obs_name}.txt",
                        overwrite=True,
                        fast_writer=False,
                    )
                if len(SuggestedPointings_1) != 0:
                    ascii.write(
                        SuggestedPointings_1,
                        f"{str(dirName)}/SuggestedPointings_GWOptimisation_{obspar[j].obs_name}.txt",
                        overwrite=True,
                        fast_writer=False,
                    )
                    if obspar[j].doRank:
                        Ranking_Space(
                            dirName,
                            f"{str(dirName)}/SuggestedPointings_GWOptimisation_{obspar[j].obs_name}.txt",
                            obspar[j],
                            obspar[j].alphaR,
                            obspar[j].betaR,
                            skymap,
                        )
                        Ranking_Space_AI(
                            dirName,
                            f"{str(dirName)}/SuggestedPointings_GWOptimisation_{obspar[j].obs_name}.txt",
                            obspar[j],
                            skymap,
                        )

        elif base in ["grid"]:
            logger.info("This is a grid and not observatory dependent")

        else:
            # for obspar in parameters:
            for j, obs in enumerate(obspar):
                obspar1 = obspar[j]
                SuggestedPointings_1 = SuggestedPointings[
                    SuggestedPointings["ObsName"] == obspar1.obs_name
                ]
                logger.info(
                    f"Suggested pointings for {obspar1.obs_name}: {SuggestedPointings_1}"
                )
                if len(SuggestedPointings_1) != 0:
                    ascii.write(
                        SuggestedPointings_1,
                        f"{str(dirName)}/SuggestedPointings_GWOptimisation_{obspar[j].obs_name}.txt",
                        overwrite=True,
                        fast_writer=False,
                    )
                    if obspar[j].doRank:
                        RankingTimes_2D(
                            obspar[j],
                            skymap.getMap("prob", obspar[j].HRnside),
                            dirName,
                            f"{str(dirName)}/SuggestedPointings_GWOptimisation_{obspar[j].obs_name}.txt",
                        )
                    if obspar[j].doPlot:
                        PointingPlotting(
                            skymap.getMap("prob", obspar[j].HRnside),
                            obspar[j],
                            obspar[j].obs_name,
                            dirName,
                            f"{str(dirName)}/SuggestedPointings_GWOptimisation_{obspar[j].obs_name}.txt",
                            obspar[j].obs_name,
                            cat,
                        )
            if obspar[j].doPlot:
                PointingPlotting(
                    skymap.getMap("prob", obspar[j].HRnside),
                    obspar[0],
                    "all",
                    dirName,
                    f"{str(dirName)}/SuggestedPointings_GWOptimisation.txt",
                    "all",
                    cat,
                )

    else:
        logger.info("No observations are scheduled")
