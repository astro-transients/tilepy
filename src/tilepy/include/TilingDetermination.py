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
#       Mid-level functions use to obtain the tiling observation scheduling
##################################################################################################


import datetime
import random

import astropy.coordinates as co
import healpy as hp
import numpy as np
import pytz
import six
from astropy import units as u
from astropy.table import Table
from six.moves import configparser
from .MaskingTools import (
    ZenithAngleCut,
    VisibleAtTime,
    FulfillsRequirement,
    FulfillsRequirementGreyObservations,
    GetBestGridPos2D,
    GetBestGridPos3D,
    OccultationCut,
    SAA_Times,
)
from .PointingTools import (
    ComputeProbability2D,
    ComputeProbGalTargeted,
    ComputeProbPGALIntegrateFoV,
    FilterGalaxies,
    Get90RegionPixReduced,
    GetBestNSIDE,
    GetSatelliteName,
    GetSatellitePositions,
    GetSatelliteTime,
    LoadGalaxies,
    LoadGalaxies_SteMgal,
    MangroveGalaxiesProbabilities,
    ModifyCatalogue,
    NextWindowTools,
    NightDarkObservation,
    NightDarkObservationwithGreyTime,
    SubstractPointings,
    SubstractPointings2D,
    Tools,
    TransformRADecToPix,
    TransformPixToRaDec,
    FindMatchingPixList,
    FindMatchingCoords,
)


if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


utc = pytz.UTC
############################################

#              General definitions              #

############################################

__all__ = [
    "PGWinFoV",
    "PGalinFoV",
    "ObservationStartperObs",
    "PGWinFoV_NObs",
    "PGalinFoV_NObs",
    "PGWinFoV_Space_NObs",
]


def PGWinFoV(skymap, nameEvent, obspar, dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule based on a 2D method.

    Parameters
    ----------
    skymap : SkyMap
        The object containing the sky maps.
    nameEvent : str
        The name of the event.
    obspar : ObservationParameters
        Observation parameters, including:
        - obsTime (datetime): Desired time for scheduling to start.
        - pointingsFile (str): Path to file with previous pointings.
        - ... (other parameters as needed)
        Class containing the telescope configuration parameters.
    dirName : str
        Output directory for schedules and plots.

    Returns
    -------
    SuggestedPointings : astropy.table.Table
        Table of scheduled pointings.
    ObservationTime0 : str
        the desired time for scheduling to start.

    """

    ObservationTime0 = obspar.obsTime
    PointingFile = obspar.pointingsFile
    # Main parameters

    print(obspar)

    random.seed()

    RAarray = []
    DECarray = []
    pixlist = []
    pixlistHR = []
    pixlist1 = []
    pixlistHR1 = []
    P_GWarray = []
    ObservationTimearray = []
    Round = []
    Fov_obs = []
    ObsName = []
    Duration = []

    print()
    print("-------------------   NEW EVENT   --------------------")
    print()

    # Retrieve maps
    prob = skymap.getMap("prob", obspar.reducedNside)
    highres = skymap.getMap("prob", obspar.HRnside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside
    )
    radecs = co.SkyCoord(rapix, decpix, frame="icrs", unit=(u.deg, u.deg))

    # Add observed pixels to pixlist
    maxRuns = obspar.maxRuns
    if PointingFile is not None:
        print(
            "==========================================================================================="
        )
        pixlist, pixlistHR, sumPGW, doneObs = SubstractPointings2D(
            PointingFile, prob, obspar, pixlist, pixlistHR
        )
        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs

        print(f"Total GW probability already covered: {sumPGW}")
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )

    #######################################################

    print("----------   NEW FOLLOW-UP ATTEMPT   ----------")

    if obspar.useGreytime:
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0, obspar)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0, obspar)

    counter = 0
    for j, NightDarkRun in enumerate(NightDarkRuns):
        if len(ObservationTimearray) < maxRuns:
            ObservationTime = NightDarkRun
            ObsBool, yprob = ZenithAngleCut(prob, ObservationTime, obspar)
            if ObsBool:
                # Round 1
                P_GW, TC, pixlist, pixlistHR = ComputeProbability2D(
                    obspar,
                    prob,
                    highres,
                    radecs,
                    ObservationTime,
                    pixlist,
                    pixlistHR,
                    counter,
                    dirName,
                )
                if (P_GW <= obspar.minProbcut) and obspar.secondRound:
                    # Try Round 2
                    # print('The minimum probability cut being', minProbcut * 100, '% is, unfortunately, not reached.')
                    yprob1 = highres
                    P_GW, TC, pixlist1, pixlistHR1 = ComputeProbability2D(
                        obspar,
                        prob,
                        yprob1,
                        radecs,
                        ObservationTime,
                        pixlist1,
                        pixlistHR1,
                        counter,
                        dirName,
                    )
                    if P_GW <= obspar.minProbcut:
                        print(
                            "Tile Pgw= ",
                            P_GW,
                            " is smaller than the minProbCut (",
                            obspar.minProbcut,
                            ") => skip this tile",
                        )
                    else:
                        Round.append(2)
                        P_GWarray.append(float("{:1.4f}".format(float(P_GW))))
                        RAarray.append(float("{:3.4f}".format(float(TC.ra.deg))))
                        DECarray.append(float("{:3.4f}".format(float(TC.dec.deg))))
                        ObservationTimearray.append(
                            ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                        )
                        ObsName.append(obspar.obs_name)
                        Duration.append(obspar.duration)
                        Fov_obs.append(obspar.FOV)
                        counter = counter + 1
                elif P_GW >= obspar.minProbcut:
                    Round.append(1)
                    P_GWarray.append(float("{:1.4f}".format(float(P_GW))))
                    RAarray.append(float("{:3.4f}".format(float(TC.ra.deg))))
                    DECarray.append(float("{:3.4f}".format(float(TC.dec.deg))))
                    ObservationTimearray.append(
                        ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                    )
                    ObsName.append(obspar.obs_name)
                    Duration.append(obspar.duration)
                    Fov_obs.append(obspar.FOV)
                    counter = counter + 1
        else:
            break

    print(
        "\nTotal GW probability covered: ",
        float("{:1.4f}".format(float(sum(P_GWarray)))),
        "Number of runs that fulfill darkness condition  :",
        len(NightDarkRuns),
        "Number of effective pointings: ",
        len(ObservationTimearray),
    )

    # List of suggested pointings
    SuggestedPointings = Table(
        [
            ObservationTimearray,
            RAarray,
            DECarray,
            P_GWarray,
            Round,
            ObsName,
            Duration,
            Fov_obs,
        ],
        names=[
            "Time[UTC]",
            "RA[deg]",
            "DEC[deg]",
            "PGW",
            "Round",
            "ObsName",
            "Duration",
            "FoV",
        ],
    )

    if len(SuggestedPointings) != 0:
        print(
            "\n================================= Tiling found ============================================="
        )
        print(SuggestedPointings)
        print(
            "============================================================================================\n"
        )
        print(f"The total probability PGW: {np.sum(P_GWarray):.4f}")
    return (SuggestedPointings, ObservationTime0)


def PGalinFoV(skymap, nameEvent, galFile, obspar, dirName):
    """
    Compute an observation schedule based on a 3D (galaxy-weighted) method.

    This mid-level function is called by :func:`tilepy.include.observationschedule.GetSchedule`
    and produces a suggested schedule of pointings for the given observatory,
    using either the target galaxy strategy or the integrated galaxy probability strategy,
    depending on user input and the telescope field of view.

    Parameters
    ----------
    skymap : SkyMap
        The object storing sky maps.
    nameEvent : str
        The name of the event.
    galFile : str
        Path to the galaxy catalog.
    obspar : ObservationParameters
        Telescope configuration parameters used in the scheduling.
    dirName : str
        Path to the output directory where the schedules and plots will be saved.

    Returns
    -------
    SuggestedPointings : astropy.table.Table
        Table of suggested pointings (with time, coordinates, probability, etc.).
    tGals0 : astropy.table.Table
        Filtered and ranked list of galaxies for scheduling.

    """

    # The desired time for scheduling to start
    ObservationTime0 = obspar.obsTime
    # The path to the text file containing the pointings that have already been performed before the scheduling
    PointingFile = obspar.pointingsFile

    # Main Parameters
    print(obspar)

    # load galaxy catalog from local file
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)

    nside = obspar.HRnside
    prob = skymap.getMap("prob", obspar.HRnside)

    if skymap.is3D:
        print("Skymap is 3D")
    else:
        print("Skymap is 2D")

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        sum_dP_dV = cat["dp_dV"].sum()
    else:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        tGals0 = MangroveGalaxiesProbabilities(tGals0)
        sum_dP_dV = cat["dp_dV"].sum()

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []

    #########################
    maxRuns = obspar.maxRuns
    if PointingFile is None:
        tGals = tGals0
        print("No pointings were given to be substracted")
    else:

        (
            ra,
            dec,
            tGals,
            AlreadyObservedPgw,
            AlreadyObservedPgal,
            alreadysumipixarray1,
            doneObs,
        ) = SubstractPointings(
            PointingFile, tGals0, alreadysumipixarray1, sum_dP_dV, prob, obspar, nside
        )
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)
        # for second round
        # ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray2, doneObs = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,prob, obspar, nside)
        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(
            f"Total GW probability already covered: {sumPGW}, "
            f"Total Gal probability already covered: {sumPGAL}"
        )
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )
    ##########################

    tGals_aux = tGals
    tGals_aux2 = tGals

    P_GALarray = []
    P_GWarray = []
    ObservationTimearray = []
    RAarray = []
    DECarray = []
    ObsName = []
    Duration = []
    Fov_obs = []

    Round = []
    print("----------   NEW FOLLOW-UP ATTEMPT   ----------")
    print(
        "maxRuns: ",
        maxRuns,
        "MinimumProbCutForCatalogue: ",
        obspar.minimumProbCutForCatalogue,
    )

    if obspar.useGreytime:
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0, obspar)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0, obspar)

    counter = 0
    if obspar.strategy == "integrated":
        for j, NightDarkRun in enumerate(NightDarkRuns):
            if len(ObservationTimearray) < maxRuns:
                ObservationTime = NightDarkRun
                visible, altaz, tGals_aux = VisibleAtTime(
                    ObservationTime, tGals_aux, obspar
                )
                if visible:
                    # select galaxies within the slightly enlarged visiblity window
                    visiMask = altaz.alt.value > 90 - (obspar.maxZenith + obspar.FOV)
                    visiGals = tGals_aux[visiMask]
                    visiGals = ModifyCatalogue(
                        prob, visiGals, obspar.FOV, sum_dP_dV, nside
                    )
                    mask, minz = FulfillsRequirement(visiGals, obspar, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar
                        )
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]
                    if finalGals["dp_dV_FOV"][0] > obspar.minProbcut:
                        # final galaxies within the FoV
                        # print(f"Condition met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} is greater than {obspar.minProbcut}")
                        if (
                            (finalGals["dp_dV_FOV"][0] < (2 * obspar.minProbcut))
                            and (sum(P_GWarray) > 0.40)
                            and obspar.secondRound
                        ):
                            visible, altaz, tGals_aux2 = VisibleAtTime(
                                ObservationTime, tGals_aux2, obspar
                            )
                            if visible:
                                visiMask = altaz.alt.value > 90 - (
                                    obspar.maxZenith + obspar.FOV
                                )
                                visiGals2 = tGals_aux2[visiMask]
                                visiGals2 = ModifyCatalogue(
                                    prob, visiGals2, obspar.FOV, sum_dP_dV, nside
                                )

                                mask, minz = FulfillsRequirement(
                                    visiGals2, obspar, UsePix=False
                                )

                                if obspar.useGreytime:
                                    maskgrey = FulfillsRequirementGreyObservations(
                                        ObservationTime, visiGals2, obspar
                                    )
                                    finalGals2 = visiGals2[mask & maskgrey]
                                if not obspar.useGreytime:
                                    finalGals2 = visiGals2[mask]

                                p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = (
                                    ComputeProbPGALIntegrateFoV(
                                        prob,
                                        ObservationTime,
                                        obspar.location,
                                        finalGals2,
                                        False,
                                        visiGals2,
                                        tGals_aux2,
                                        sum_dP_dV,
                                        alreadysumipixarray2,
                                        nside,
                                        minz,
                                        obspar,
                                        counter,
                                        nameEvent,
                                        dirName,
                                        obspar.doPlot,
                                    )
                                )

                                RAarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals2["RAJ2000"][:1])
                                        )
                                    )
                                )
                                DECarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals2["DEJ2000"][:1])
                                        )
                                    )
                                )
                                Round.append(2)
                                P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                ObservationTimearray.append(
                                    ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                )
                                ObsName.append(obspar.obs_name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1

                            else:
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                    ComputeProbPGALIntegrateFoV(
                                        prob,
                                        ObservationTime,
                                        obspar.location,
                                        finalGals,
                                        False,
                                        visiGals,
                                        tGals_aux,
                                        sum_dP_dV,
                                        alreadysumipixarray1,
                                        nside,
                                        minz,
                                        obspar,
                                        counter,
                                        nameEvent,
                                        dirName,
                                        obspar.doPlot,
                                    )
                                )
                                RAarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["RAJ2000"][:1])
                                        )
                                    )
                                )
                                DECarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["DEJ2000"][:1])
                                        )
                                    )
                                )
                                Round.append(1)
                                P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                ObservationTimearray.append(
                                    ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                )
                                ObsName.append(obspar.obs_name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1
                        else:
                            # print("We are in round 1")
                            # print("\n=================================")
                            # print("TARGET COORDINATES AND DETAILS...")
                            # print("=================================")
                            p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                ComputeProbPGALIntegrateFoV(
                                    prob,
                                    ObservationTime,
                                    obspar.location,
                                    finalGals,
                                    False,
                                    visiGals,
                                    tGals_aux,
                                    sum_dP_dV,
                                    alreadysumipixarray1,
                                    nside,
                                    minz,
                                    obspar,
                                    counter,
                                    nameEvent,
                                    dirName,
                                    obspar.doPlot,
                                )
                            )

                            RAarray.append(
                                float("{:3.4f}".format(float(finalGals["RAJ2000"][:1])))
                            )
                            DECarray.append(
                                float("{:3.4f}".format(float(finalGals["DEJ2000"][:1])))
                            )
                            Round.append(1)
                            P_GALarray.append(float("{:1.4f}".format(p_gal)))
                            P_GWarray.append(float("{:1.4f}".format(p_gw)))
                            ObservationTimearray.append(
                                ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                            )
                            ObsName.append(obspar.obs_name)
                            Duration.append(obspar.duration)
                            Fov_obs.append(obspar.FOV)
                            counter = counter + 1

                    else:
                        print(
                            f"Condition not met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} must be greater than {obspar.minProbcut}"
                        )
            else:
                break
    if obspar.strategy == "targeted":
        for j, NightDarkRun in enumerate(NightDarkRuns):
            if len(ObservationTimearray) < maxRuns:
                ObservationTime = NightDarkRun
                visible, altaz, tGals_aux = VisibleAtTime(
                    ObservationTime, tGals_aux, obspar
                )
                if visible:
                    # select galaxies within the slightly enlarged visiblity window
                    visiMask = altaz.alt.value > 90 - (obspar.maxZenith + obspar.FOV)
                    visiGals = tGals_aux[visiMask]
                    mask, minz = FulfillsRequirement(visiGals, obspar, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar
                        )
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]
                    # print('finalGals', finalGals,tGals['dp_dV'][:1]*obspar.minProbcut)
                    if finalGals["dp_dV"][0] > tGals["dp_dV"][0] * obspar.minProbcut:
                        # print(f"Condition met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV'][0]} is greater than {obspar.minProbcut}")
                        if (
                            (finalGals["dp_dV"][0] < (2 * obspar.minProbcut))
                            and (sum(P_GWarray) > 0.40)
                            and obspar.secondRound
                        ):
                            visible, altaz, tGals_aux2 = VisibleAtTime(
                                ObservationTime, tGals_aux2, obspar
                            )
                            if visible:
                                visiMask = altaz.alt.value > 90 - (
                                    obspar.maxZenith + obspar.FOV
                                )
                                visiGals2 = tGals_aux2[visiMask]
                                mask, minz = FulfillsRequirement(
                                    visiGals2, obspar, UsePix=False
                                )

                                if obspar.useGreytime:
                                    maskgrey = FulfillsRequirementGreyObservations(
                                        ObservationTime, visiGals2, obspar
                                    )
                                    finalGals2 = visiGals2[mask & maskgrey]
                                if not obspar.useGreytime:
                                    finalGals2 = visiGals2[mask]
                                p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = (
                                    ComputeProbGalTargeted(
                                        prob,
                                        ObservationTime,
                                        finalGals2,
                                        visiGals2,
                                        tGals_aux2,
                                        sum_dP_dV,
                                        alreadysumipixarray2,
                                        nside,
                                        minz,
                                        obspar,
                                        counter,
                                        dirName,
                                    )
                                )

                                RAarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals2["RAJ2000"][:1])
                                        )
                                    )
                                )
                                DECarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals2["DEJ2000"][:1])
                                        )
                                    )
                                )
                                Round.append(2)
                                P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                ObservationTimearray.append(
                                    ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                )
                                ObsName.append(obspar.obs_name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1

                            else:
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                    ComputeProbGalTargeted(
                                        prob,
                                        ObservationTime,
                                        finalGals,
                                        visiGals,
                                        tGals_aux,
                                        sum_dP_dV,
                                        alreadysumipixarray1,
                                        nside,
                                        minz,
                                        obspar,
                                        counter,
                                        dirName,
                                    )
                                )
                                RAarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["RAJ2000"][:1])
                                        )
                                    )
                                )
                                DECarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["DEJ2000"][:1])
                                        )
                                    )
                                )
                                Round.append(1)
                                P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                ObservationTimearray.append(
                                    ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                )
                                ObsName.append(obspar.obs_name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1
                        else:
                            # print("We are in round 1")
                            p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                ComputeProbGalTargeted(
                                    prob,
                                    ObservationTime,
                                    finalGals,
                                    visiGals,
                                    tGals_aux,
                                    sum_dP_dV,
                                    alreadysumipixarray1,
                                    nside,
                                    minz,
                                    obspar,
                                    counter,
                                    dirName,
                                )
                            )
                            RAarray.append(
                                float("{:3.4f}".format(float(finalGals["RAJ2000"][:1])))
                            )
                            DECarray.append(
                                float("{:3.4f}".format(float(finalGals["DEJ2000"][:1])))
                            )
                            Round.append(1)
                            P_GALarray.append(float("{:1.4f}".format(p_gal)))
                            P_GWarray.append(float("{:1.4f}".format(p_gw)))
                            ObservationTimearray.append(
                                ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                            )
                            ObsName.append(obspar.obs_name)
                            Duration.append(obspar.duration)
                            Fov_obs.append(obspar.FOV)
                            counter = counter + 1

                    else:
                        print(
                            f"Condition not met at {ObservationTime}: dp/dV = {finalGals['dp_dV'][0]} must be greater than {obspar.minProbcut}"
                        )

            else:
                break

    # List of suggested pointings
    SuggestedPointings = Table(
        [
            ObservationTimearray,
            RAarray,
            DECarray,
            P_GWarray,
            P_GALarray,
            Round,
            ObsName,
            Duration,
            Fov_obs,
        ],
        names=[
            "Time[UTC]",
            "RA[deg]",
            "DEC[deg]",
            "PGW",
            "Pgal",
            "Round",
            "ObsName",
            "Duration",
            "FoV",
        ],
    )

    if len(SuggestedPointings) != 0:
        print(
            "\n================================= Tiling found ============================================="
        )
        print(SuggestedPointings)
        print(
            "============================================================================================\n"
        )
        print(f"The total probability PGal: {np.sum(P_GALarray):.4f}")
        print(f"The total probability PGW: {np.sum(P_GWarray):.4f}")
    return SuggestedPointings, tGals0


def ObservationStartperObs(obsparameters, ObservationTime0):
    """
    Compute the first observation time for each observatory involved in the scheduling.

    This mid-level function is called by Nobs Tiling functions to determine the first available observation
    time for each observatory and whether each observatory can observe on the same night.

    Parameters
    ----------
    obsparameters : list of ObservationParameters
        A list of sets of parameters for each observatory needed to launch the tiling scheduler.
    ObservationTime0 : sr
        The desired start time for scheduling to begin.

    Returns
    -------
    obs_time : datetime
        The current observation time, possibly adjusted.
    SameNight : numpy.ndarray of bool
        An array indicating whether each observatory is available for observation on the same night.
    NewActiveObs : list of ObservationParameters
        A list of observatories that are available to observe.
    NewActiveObsStart : numpy.ndarray of datetime
        A sorted list of the first available observation times for each observatory.

    """

    print("ObservationTime0", ObservationTime0)

    print("obsparameters", len(obsparameters))
    # Finding the start time for each observatory and checking if it's now
    FirstDark = np.full(len(obsparameters), False, dtype=bool)
    FirstDark_Flag = np.full(len(obsparameters), False, dtype=bool)
    # print(len(ObsFirstTime))
    obs_time = ObservationTime0
    if obs_time.tzinfo is None:
        obs_time = utc.localize(obs_time)
    ObsFirstTime = []

    j = 0
    for obspar1 in obsparameters:
        if obsparameters[j].base == "space":
            dark_at_start = True
            FirstDark[j] = dark_at_start

        else:
            dark_at_start = False

            if obsparameters[j].useGreytime:
                dark_at_start = Tools.CheckWindowGrey(obs_time, obsparameters[j])
            if not obsparameters[j].useGreytime:
                dark_at_start = Tools.CheckWindow(obs_time, obsparameters[j])
            FirstDark[j] = dark_at_start

        # THIS WILL CREATE A DATETIME OBJECT WITH IN THE FORM XX+00:00 WITH NO DOTS
        if FirstDark[j]:
            FirstDark_Flag[j] = True
            if obs_time.tzinfo is None:
                obs_time = utc.localize(obs_time)
            ObsFirstTime.append(obs_time)
        else:  # THIS WILL CREATE A DATETIME OBJECT WITH IN THE FORM .XX+00:00
            if obsparameters[j].useGreytime:
                ObsFirstTime1 = NextWindowTools.NextObservationWindowGrey(
                    time=obs_time, obspar=obsparameters[j]
                )
                ObsFirstTime.append(ObsFirstTime1)
            if not obsparameters[j].useGreytime:
                ObsFirstTime1 = NextWindowTools.NextObservationWindow(
                    time=obs_time, obspar=obsparameters[j]
                )
                ObsFirstTime.append(ObsFirstTime1)
            if ObsFirstTime1:
                if ObsFirstTime1.tzinfo is None:
                    ObsFirstTime1 = utc.localize(ObsFirstTime1)
                if obs_time.tzinfo is None:
                    obs_time = utc.localize(obs_time)
                if ObsFirstTime1 < obs_time + datetime.timedelta(hours=24):
                    FirstDark_Flag[j] = True
        j += 1

    # Checking which observatories are availabe for observations and saving their start time
    ActiveObsStart = []
    ActiveObs = []
    SameNight = np.full(len(obsparameters), False, dtype=bool)

    j = 0
    for obspar in obsparameters:
        if FirstDark_Flag[j]:
            if ObsFirstTime[j].tzinfo is None:
                ObsFirstTime = utc.localize(ObsFirstTime[j])
            ActiveObsStart.append(ObsFirstTime[j])
            ActiveObs.append(obsparameters[j])
            SameNight[j] = True
        j += 1

    # Sorting observatories according to their first obsevation time available
    NewActiveObsStart = np.sort(ActiveObsStart)
    NewActiveObs = ActiveObs
    ind = np.argsort(ActiveObsStart)
    ind = np.array(ind)
    NewActiveObs = np.take(ActiveObs, ind)

    return obs_time, SameNight, NewActiveObs, NewActiveObsStart


def PGWinFoV_NObs(
    skymap, nameEvent, ObservationTime0, PointingFile, obsparameters, dirName
):
    """
    Compute an observation schedule for multiple telescopes/observatories based on a 2D method.

    This mid-level function is called by `GetSchedule` (see :func:`tilepy.include.ObservationScheduler.GetSchedule`)
    and produces a suggested schedule of pointings for each observatory, given the input sky map,
    pointings already performed, and observation parameters.

    Parameters
    ----------
    skymap : SkyMap
        The object storing sky maps.
    nameEvent : str
        The name of the event.
    ObservationTime0 : str
        The desired start time for scheduling.
    PointingFile : str
        Path to the text file containing pointings already performed before scheduling.
    obsparameters : list of ObservationParameters
        Parameters for each observatory needed to launch the tiling scheduler.
    dirName : str
        Path to the output directory where the schedules and plots will be saved.

    Returns
    -------
    SuggestedPointings : astropy.table.Table
        Table of suggested pointings for all observatories.
    obsparameters : list of ObservationParameters
        (Possibly updated) list of parameters for each observatory.

    """

    obs_time, SameNight, NewActiveObs, NewActiveObsStart = ObservationStartperObs(
        obsparameters, ObservationTime0
    )

    # START
    #################################################################################################################################################
    random.seed()
    RAarray = []
    DECarray = []
    pixlist = []
    pixlistHR = []
    pixlist1 = []
    pixlistHR1 = []
    P_GWarray = []
    ObservationTimearray = []
    Round = []
    ObsName = []
    Duration = []
    Fov_obs = []
    #################################################################################################################################################
    obspar = obsparameters[0]

    # Retrieve maps
    prob = skymap.getMap("prob", obspar.reducedNside)
    highres = skymap.getMap("prob", obspar.HRnside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside
    )
    radecs = co.SkyCoord(rapix, decpix, frame="icrs", unit=(u.deg, u.deg))
    maxRuns = obspar.maxRuns
    # Add observed pixels to pixlist
    if PointingFile is not None:
        print(PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist)
        pixlist, pixlistHR, sumPGW, doneObs = SubstractPointings2D(
            PointingFile, prob, obspar, pixlist, pixlistHR
        )

        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(f"Total GW probability already covered: {sumPGW}")
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )
    #################################################################################################################################################
    ITERATION_OBS = 0
    TIME_MIN_ALL = []
    TIME_MIN = obs_time + datetime.timedelta(hours=12)
    if TIME_MIN.tzinfo is None:
        TIME_MIN = utc.localize(TIME_MIN)
    NewActiveObsTime = NewActiveObsStart
    NUMBER_OBS = np.zeros(len(NewActiveObs))
    #################################################################################################################################################
    counter = 0
    i = 0
    couter_per_obs = np.zeros(len(NewActiveObs))
    print("------NewActiveObsTime--------", NewActiveObs[0].obs_name)
    while (i < 500) & any(SameNight):
        for j, obs in enumerate(NewActiveObs):
            obspar = NewActiveObs[j]
            print("Observatory: ", obspar.obs_name)
            ObservationTime = NewActiveObsTime[j]
            if ITERATION_OBS == len(obsparameters):
                TIME_MIN_ALL = []
                ITERATION_OBS = 0
            ITERATION_OBS += 1

            if obspar.base == "space":
                SatelliteName = GetSatelliteName(obspar.obs_name, obspar.stationsurl)
                print(obspar.obs_name)

            if couter_per_obs[j] >= maxRuns:
                SameNight[j] = False
            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
                pixlistHROcc = None
                if obspar.base == "space":
                    SatelliteTime = GetSatelliteTime(SatelliteName, ObservationTime)
                    satellitePosition, satelliteLocation = GetSatellitePositions(
                        SatelliteName, SatelliteTime
                    )
                    ObsBool, yprob, pixlistHROcc = OccultationCut(
                        prob,
                        obspar.reducedNside,
                        ObservationTime,
                        obspar,
                        satellitePosition,
                        satelliteLocation,
                    )
                else:
                    ObsBool, yprob = ZenithAngleCut(prob, ObservationTime, obspar)

                if ObsBool:
                    # Round 1
                    P_GW, TC, pixlist, pixlistHR = ComputeProbability2D(
                        obspar,
                        prob,
                        highres,
                        radecs,
                        ObservationTime,
                        pixlist,
                        pixlistHR,
                        counter,
                        dirName,
                        pixlistHROcc,
                    )
                    # print(P_GW, obspar.name)
                    if (P_GW <= obspar.minProbcut) and obspar.secondRound:
                        # Try Round 2
                        # print('The minimum probability cut being', minProbcut * 100, '% is, unfortunately, not reached.')
                        yprob1 = highres
                        P_GW, TC, pixlist1, pixlistHR1 = ComputeProbability2D(
                            prob,
                            yprob1,
                            radecs,
                            ObservationTime,
                            pixlist1,
                            pixlistHR1,
                            counter,
                            dirName,
                            pixlistHROcc,
                        )
                        if P_GW <= obspar.minProbcut:
                            print(
                                "Tile Pgw= ",
                                P_GW,
                                " is smaller than the minProbCut (",
                                obspar.minProbCut,
                                ") => skip this tile",
                            )
                        else:
                            Round.append(2)
                            P_GWarray.append(float("{:1.4f}".format(float(P_GW))))
                            RAarray.append(float("{:3.4f}".format(float(TC.ra.deg))))
                            DECarray.append(float("{:3.4f}".format(float(TC.dec.deg))))
                            ObservationTime = str(ObservationTime).split("+")[0]
                            try:
                                ObservationTime = datetime.datetime.strptime(
                                    ObservationTime, "%Y-%m-%d %H:%M:%S"
                                )
                            except ValueError:
                                ObservationTime = str(ObservationTime).split(".")[0]
                                ObservationTime = datetime.datetime.strptime(
                                    ObservationTime, "%Y-%m-%d %H:%M:%S"
                                )
                            ObservationTimearray.append(ObservationTime)
                            ObsName.append(obspar.obs_name)
                            counter = counter + 1
                            couter_per_obs[j] += 1
                            Duration.append(obspar.duration)
                            Fov_obs.append(obspar.FOV)

                    elif P_GW >= obspar.minProbcut:
                        Round.append(1)
                        P_GWarray.append(float("{:1.4f}".format(float(P_GW))))
                        RAarray.append(float("{:3.4f}".format(float(TC.ra.deg))))
                        DECarray.append(float("{:3.4f}".format(float(TC.dec.deg))))
                        ObservationTime = str(ObservationTime).split("+")[0]
                        try:
                            ObservationTime = datetime.datetime.strptime(
                                ObservationTime, "%Y-%m-%d %H:%M:%S"
                            )
                        except ValueError:
                            ObservationTime = str(ObservationTime).split(".")[0]
                            ObservationTime = datetime.datetime.strptime(
                                ObservationTime, "%Y-%m-%d %H:%M:%S"
                            )
                        ObservationTimearray.append(ObservationTime)
                        ObsName.append(obspar.obs_name)
                        counter = counter + 1
                        couter_per_obs[j] += 1
                        Duration.append(obspar.duration)
                        Fov_obs.append(obspar.FOV)

                # HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
                NewActiveObsTime[j] = NewActiveObsTime[j] + datetime.timedelta(
                    minutes=obspar.duration
                )

                # HERE WE DETERMINE IF WE ARE STILL IN THE SAME NIGHT FOR THIS OBSERVATORY
                # if (NewActiveObsTime[j] > Tools.NextSunrise(obsstart, obspar)) | (obsstart > Tools.NextMoonrise(obsstart, obspar)):

                if obspar.base == "space":
                    if NewActiveObsTime[j] > obs_time + datetime.timedelta(hours=24):
                        SameNight[j] = False
                else:
                    if not obsparameters[j].useGreytime:
                        if not Tools.CheckWindow(NewActiveObsTime[j], obspar):
                            SameNight[j] = False
                    if obsparameters[j].useGreytime:
                        if not Tools.CheckWindowGrey(NewActiveObsTime[j], obspar):
                            SameNight[j] = False
                    if obspar.sunDown > 10:
                        if NewActiveObsTime[j] > obs_time + datetime.timedelta(
                            hours=24
                        ):
                            SameNight[j] = False

                NUMBER_OBS[j] += 1

            if SameNight[j]:
                TIME_MIN = NewActiveObsTime[j]
                TIME_MIN_ALL.append(TIME_MIN)
                TIME_MIN = np.min(TIME_MIN_ALL)
            else:
                TIME_MIN = TIME_MIN + datetime.timedelta(hours=12)

        i += 1

    SuggestedPointings = Table(
        [
            ObservationTimearray,
            RAarray,
            DECarray,
            P_GWarray,
            Round,
            ObsName,
            Duration,
            Fov_obs,
        ],
        names=[
            "Time[UTC]",
            "RA[deg]",
            "DEC[deg]",
            "PGW",
            "Round",
            "ObsName",
            "Duration",
            "FoV",
        ],
    )
    print("The total probability PGW: ", np.sum(P_GWarray))

    return SuggestedPointings, obsparameters


def PGalinFoV_NObs(
    skymap, nameEvent, ObservationTime0, PointingFile, galFile, obsparameters, dirName
):
    """
    Compute the optimal observation schedule based on galaxy probability and gravitational wave data for multiple observatories.
    """

    obs_time, SameNight, NewActiveObs, NewActiveObsStart = ObservationStartperObs(
        obsparameters, ObservationTime0
    )
    # START
    #################################################################################################################################################
    random.seed()
    P_GALarray = []
    P_GWarray = []
    ObservationTimearray = []
    RAarray = []
    DECarray = []
    alreadysumipixarray1 = []
    alreadysumipixarray2 = []
    Round = []
    ObsName = []
    Duration = []
    Fov_obs = []
    #################################################################################################################################################
    obspar = obsparameters[0]

    nside = obspar.HRnside
    prob = skymap.getMap("prob", obspar.HRnside)

    # load galaxy catalogue
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        sum_dP_dV = cat["dp_dV"].sum()
    else:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        tGals0 = MangroveGalaxiesProbabilities(tGals0)
        sum_dP_dV = cat["dp_dV"].sum()

    # Add observed pixels to pixlist
    maxRuns = obspar.maxRuns
    if PointingFile is None:
        tGals = tGals0
        print("No pointings were given to be substracted")
    else:
        # tGals_aux = tGals
        (
            ra,
            dec,
            tGals,
            AlreadyObservedPgw,
            AlreadyObservedPgal,
            alreadysumipixarray1,
            doneObs,
        ) = SubstractPointings(
            PointingFile, tGals0, alreadysumipixarray1, sum_dP_dV, prob, obspar, nside
        )
        maxRuns = obspar.maxRuns - len(ra)
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)
        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(
            f"Total GW probability already covered: {sumPGW}, "
            f"Total Gal probability already covered: {sumPGAL}"
        )
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )

    tGals_aux = tGals
    tGals_aux2 = tGals
    #################################################################################################################################################
    ITERATION_OBS = 0
    TIME_MIN_ALL = []
    TIME_MIN = obs_time + datetime.timedelta(hours=12)
    if TIME_MIN.tzinfo is None:
        TIME_MIN = utc.localize(TIME_MIN)
    NewActiveObsTime = NewActiveObsStart
    NUMBER_OBS = np.zeros(len(NewActiveObs))
    #################################################################################################################################################

    counter = 0
    i = 0
    couter_per_obs = np.zeros(len(NewActiveObs))
    while (i < 5000) & any(SameNight):
        for j, obs in enumerate(NewActiveObs):
            obspar = NewActiveObs[j]
            ObservationTime = NewActiveObsTime[j]
            if ITERATION_OBS == len(obsparameters):
                TIME_MIN_ALL = []
                ITERATION_OBS = 0

            ITERATION_OBS += 1
            if couter_per_obs[j] >= maxRuns:
                SameNight[j] = False
            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:

                if obspar.strategy == "integrated":
                    visible, altaz, tGals_aux = VisibleAtTime(
                        ObservationTime, tGals_aux, obspar
                    )

                    if visible:

                        # select galaxies within the slightly enlarged visiblity window
                        visiMask = altaz.alt.value > 90 - (
                            obspar.maxZenith + obspar.FOV
                        )
                        visiGals = tGals_aux[visiMask]
                        visiGals = ModifyCatalogue(
                            prob, visiGals, obspar.FOV, sum_dP_dV, nside
                        )

                        mask, minz = FulfillsRequirement(visiGals, obspar, UsePix=False)
                        if obspar.useGreytime:
                            maskgrey = FulfillsRequirementGreyObservations(
                                ObservationTime, visiGals, obspar
                            )
                            finalGals = visiGals[mask & maskgrey]
                        if not obspar.useGreytime:
                            finalGals = visiGals[mask]

                        if finalGals["dp_dV_FOV"][0] > obspar.minProbcut:
                            # final galaxies within the FoV
                            # print(f"Condition met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} is greater than {obspar.minProbcut}")
                            if (
                                (finalGals["dp_dV_FOV"][0] < (2 * obspar.minProbcut))
                                and (sum(P_GWarray) > 0.40)
                                and obspar.secondRound
                            ):
                                print("probability", finalGals["dp_dV_FOV"][:1])
                                visible, altaz, tGals_aux2 = VisibleAtTime(
                                    ObservationTime, tGals_aux2, obspar
                                )
                                if visible:
                                    visiMask = altaz.alt.value > 90 - (
                                        obspar.maxZenith + obspar.FOV
                                    )
                                    visiGals2 = tGals_aux2[visiMask]
                                    visiGals2 = ModifyCatalogue(
                                        prob, visiGals2, obspar.FOV, sum_dP_dV, nside
                                    )

                                    mask, minz = FulfillsRequirement(
                                        visiGals2, obspar, UsePix=False
                                    )

                                    if obspar.useGreytime:
                                        maskgrey = FulfillsRequirementGreyObservations(
                                            ObservationTime, visiGals2, obspar
                                        )
                                        finalGals2 = visiGals2[mask & maskgrey]
                                    if not obspar.useGreytime:
                                        finalGals2 = visiGals2[mask]

                                    p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = (
                                        ComputeProbPGALIntegrateFoV(
                                            prob,
                                            ObservationTime,
                                            obspar.location,
                                            finalGals2,
                                            False,
                                            visiGals2,
                                            tGals_aux2,
                                            sum_dP_dV,
                                            alreadysumipixarray2,
                                            nside,
                                            minz,
                                            obspar,
                                            counter,
                                            nameEvent,
                                            dirName,
                                            obspar.doPlot,
                                        )
                                    )
                                    RAarray.append(
                                        float(
                                            "{:3.4f}".format(
                                                float(finalGals2["RAJ2000"][:1])
                                            )
                                        )
                                    )
                                    DECarray.append(
                                        float(
                                            "{:3.4f}".format(
                                                float(finalGals2["DEJ2000"][:1])
                                            )
                                        )
                                    )
                                    Round.append(2)
                                    P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                    P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                    ObservationTime = str(ObservationTime).split(".")[0]
                                    try:
                                        ObservationTime = datetime.datetime.strptime(
                                            ObservationTime, "%Y-%m-%d %H:%M:%S"
                                        )
                                    except Exception:
                                        ObservationTime = str(ObservationTime).split(
                                            "+"
                                        )[0]
                                        ObservationTime = datetime.datetime.strptime(
                                            ObservationTime, "%Y-%m-%d %H:%M:%S"
                                        )
                                    ObservationTimearray.append(ObservationTime)
                                    ObsName.append(obspar.obs_name)
                                    counter = counter + 1
                                    couter_per_obs[j] += 1
                                    Duration.append(obspar.duration)
                                    Fov_obs.append(obspar.FOV)

                            else:
                                # print("\n=================================")
                                # print("TARGET COORDINATES AND DETAILS...")
                                # print("=================================")
                                # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                    ComputeProbPGALIntegrateFoV(
                                        prob,
                                        ObservationTime,
                                        obspar.location,
                                        finalGals,
                                        False,
                                        visiGals,
                                        tGals_aux,
                                        sum_dP_dV,
                                        alreadysumipixarray1,
                                        nside,
                                        minz,
                                        obspar,
                                        counter,
                                        nameEvent,
                                        dirName,
                                        obspar.doPlot,
                                    )
                                )
                                RAarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["RAJ2000"][:1])
                                        )
                                    )
                                )
                                DECarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["DEJ2000"][:1])
                                        )
                                    )
                                )
                                Round.append(1)
                                P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                ObservationTime = str(ObservationTime).split(".")[0]
                                try:
                                    ObservationTime = datetime.datetime.strptime(
                                        ObservationTime, "%Y-%m-%d %H:%M:%S"
                                    )
                                except Exception:
                                    ObservationTime = str(ObservationTime).split("+")[0]
                                    ObservationTime = datetime.datetime.strptime(
                                        ObservationTime, "%Y-%m-%d %H:%M:%S"
                                    )
                                ObservationTimearray.append(ObservationTime)
                                ObsName.append(obspar.obs_name)
                                counter = counter + 1
                                couter_per_obs[j] += 1
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                            # ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))

                        else:
                            print(
                                f"Condition NOT met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} is greater than {obspar.minProbcut}"
                            )

                if obspar.strategy == "targeted":
                    visible, altaz, tGals_aux = VisibleAtTime(
                        ObservationTime, tGals_aux, obspar
                    )
                    if visible:
                        # select galaxies within the slightly enlarged visiblity window
                        visiMask = altaz.alt.value > 90 - (
                            obspar.maxZenith + obspar.FOV
                        )
                        visiGals = tGals_aux[visiMask]
                        mask, minz = FulfillsRequirement(visiGals, obspar, UsePix=False)
                        if obspar.useGreytime:
                            maskgrey = FulfillsRequirementGreyObservations(
                                ObservationTime,
                                visiGals,
                                obspar.location,
                                obspar.minMoonSourceSeparation,
                            )
                            finalGals = visiGals[mask & maskgrey]
                        if not obspar.useGreytime:
                            finalGals = visiGals[mask]
                        # print('finalGals', finalGals,tGals['dp_dV'][:1]*obspar.minProbcut)
                        if (
                            finalGals["dp_dV"][0]
                            > tGals["dp_dV"][:1] * obspar.minProbcut
                        ):
                            print(
                                f"Condition met at {ObservationTime}: dp/dV = {finalGals['dp_dV'][0]} is greater than {obspar.minProbcut}"
                            )
                            if (
                                (finalGals["dp_dV"][0] < (2 * obspar.minProbcut))
                                and (sum(P_GWarray) > 0.40)
                                and obspar.secondRound
                            ):
                                visible, altaz, tGals_aux2 = VisibleAtTime(
                                    ObservationTime, tGals_aux2, obspar
                                )
                                if visible:
                                    visiMask = altaz.alt.value > 90 - (
                                        obspar.maxZenith + obspar.FOV
                                    )
                                    visiGals2 = tGals_aux2[visiMask]
                                    mask, minz = FulfillsRequirement(
                                        visiGals2, obspar, UsePix=False
                                    )

                                    if obspar.useGreytime:
                                        maskgrey = FulfillsRequirementGreyObservations(
                                            ObservationTime, visiGals2, obspar
                                        )
                                        finalGals2 = visiGals2[mask & maskgrey]
                                    if not obspar.useGreytime:
                                        finalGals2 = visiGals2[mask]
                                    p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = (
                                        ComputeProbGalTargeted(
                                            prob,
                                            ObservationTime,
                                            finalGals2,
                                            visiGals2,
                                            tGals_aux2,
                                            sum_dP_dV,
                                            alreadysumipixarray2,
                                            nside,
                                            minz,
                                            obspar,
                                            counter,
                                            dirName,
                                        )
                                    )

                                    RAarray.append(
                                        float(
                                            "{:3.4f}".format(
                                                float(finalGals2["RAJ2000"][:1])
                                            )
                                        )
                                    )
                                    DECarray.append(
                                        float(
                                            "{:3.4f}".format(
                                                float(finalGals2["DEJ2000"][:1])
                                            )
                                        )
                                    )
                                    Round.append(2)
                                    P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                    P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                    ObservationTimearray.append(
                                        ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                    )
                                    ObsName.append(obspar.obs_name)
                                    counter = counter + 1
                                    couter_per_obs[j] += 1
                                    Duration.append(obspar.duration)
                                    Fov_obs.append(obspar.FOV)

                                else:
                                    p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                        ComputeProbGalTargeted(
                                            prob,
                                            ObservationTime,
                                            finalGals,
                                            visiGals,
                                            tGals_aux,
                                            sum_dP_dV,
                                            alreadysumipixarray1,
                                            nside,
                                            minz,
                                            obspar,
                                            counter,
                                            dirName,
                                        )
                                    )
                                    RAarray.append(
                                        float(
                                            "{:3.4f}".format(
                                                float(finalGals["RAJ2000"][:1])
                                            )
                                        )
                                    )
                                    DECarray.append(
                                        float(
                                            "{:3.4f}".format(
                                                float(finalGals["DEJ2000"][:1])
                                            )
                                        )
                                    )
                                    Round.append(1)
                                    P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                    P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                    ObservationTimearray.append(
                                        ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                    )
                                    ObsName.append(obspar.obs_name)
                                    counter = counter + 1
                                    couter_per_obs[j] += 1
                                    Duration.append(obspar.duration)
                                    Fov_obs.append(obspar.FOV)
                            else:
                                # print("We are in round 1")
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = (
                                    ComputeProbGalTargeted(
                                        prob,
                                        ObservationTime,
                                        finalGals,
                                        visiGals,
                                        tGals_aux,
                                        sum_dP_dV,
                                        alreadysumipixarray1,
                                        nside,
                                        minz,
                                        obspar,
                                        counter,
                                        dirName,
                                    )
                                )
                                RAarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["RAJ2000"][:1])
                                        )
                                    )
                                )
                                DECarray.append(
                                    float(
                                        "{:3.4f}".format(
                                            float(finalGals["DEJ2000"][:1])
                                        )
                                    )
                                )
                                Round.append(1)
                                P_GALarray.append(float("{:1.4f}".format(p_gal)))
                                P_GWarray.append(float("{:1.4f}".format(p_gw)))
                                ObservationTimearray.append(
                                    ObservationTime.strftime("%Y-%m-%d %H:%M:%S")
                                )
                                ObsName.append(obspar.obs_name)
                                counter = counter + 1
                                couter_per_obs[j] += 1
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)

                        else:
                            print(
                                f"Condition NOT met at {ObservationTime}: dp/dV = {finalGals['dp_dV'][0]} is greater than {obspar.minProbcut}"
                            )

                # HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
                NewActiveObsTime[j] = NewActiveObsTime[j] + datetime.timedelta(
                    minutes=obspar.duration
                )

                # HERE WE DETERMINE IF WE ARE STILL IN THE SAME NIGHT FOR THIS OBSERVATORY
                # if (NewActiveObsTime[j] > Tools.NextSunrise(NewActiveObsStart[j], NewActiveObs[j])) | (NewActiveObsStart[j] > Tools.NextMoonrise(obs_time, NewActiveObs[j])):
                if not obsparameters[j].useGreytime:
                    if not Tools.CheckWindow(NewActiveObsTime[j], obspar):
                        SameNight[j] = False
                if obsparameters[j].useGreytime:
                    if not Tools.CheckWindowGrey(NewActiveObsTime[j], obspar):
                        SameNight[j] = False
                if obspar.sunDown > 10:
                    if NewActiveObsTime[j] > obs_time + datetime.timedelta(hours=24):
                        SameNight[j] = False

                NUMBER_OBS[j] += 1

            if SameNight[j]:
                TIME_MIN = NewActiveObsTime[j]
                TIME_MIN_ALL.append(TIME_MIN)
                TIME_MIN = np.min(TIME_MIN_ALL)
            else:
                TIME_MIN = TIME_MIN + datetime.timedelta(hours=12)

        i += 1

    SuggestedPointings = Table(
        [
            ObservationTimearray,
            RAarray,
            DECarray,
            P_GWarray,
            P_GALarray,
            Round,
            ObsName,
            Duration,
            Fov_obs,
        ],
        names=[
            "Time[UTC]",
            "RA[deg]",
            "DEC[deg]",
            "PGW",
            "Pgal",
            "Round",
            "ObsName",
            "Duration",
            "FoV",
        ],
    )
    print("The total probability PGal: ", np.sum(P_GALarray))
    print("The total probability PGW: ", np.sum(P_GWarray))
    return SuggestedPointings, tGals0, obsparameters


def GetBestTiles2D(skymap, nameEvent, PointingFile, obsparameters, dirName):
    """
    Compute the best observation tiles based on a 2D method, considering galaxy probabilities and observatory constraints.
    """

    random.seed()
    RAarray = []
    DECarray = []
    P_GWarray = []
    Occultedpixels = []
    pixlistHR = []
    obspar = obsparameters[0]

    radius = obspar.FOV
    HRnside, reducedNside = GetBestNSIDE(obspar.reducedNside, obspar.HRnside, radius)

    # Retrieve maps
    prob = skymap.getMap("prob", reducedNside)
    highres = skymap.getMap("prob", HRnside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, _ = Get90RegionPixReduced(prob, obspar.percentageMOC, reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame="icrs", unit=(u.deg, u.deg))
    maxRuns = obspar.maxRuns

    doPlot = obspar.doPlot

    # Add observed pixels to pixlist
    if PointingFile is not None:

        # FIXME: pixlist is undefined in this scope
        # The program will crash is the if branch is executed
        print(
            PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist  # noqa: F821
        )
        pixlist, pixlistHR, sumPGW, doneObs = SubstractPointings2D(
            PointingFile, prob, obspar, pixlist, pixlistHR  # noqa: F821
        )

        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(f"Total GW probability already covered: {sumPGW}")
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )

    ipix = TransformRADecToPix(radecs, reducedNside)
    newpix = ipix

    first_values = GetBestGridPos2D(
        prob,
        highres,
        HRnside,
        reducedNside,
        newpix,
        radius,
        maxRuns,
        Occultedpixels,
        doPlot,
        dirName,
        obspar.numberSides,
        pixlistHR,
        obspar.minProbcut,
    )

    # ObsName = [obspar.name for j in range(len(first_values))]
    RAarray = [row["PIXRA"] for row in first_values]
    DECarray = [row["PIXDEC"] for row in first_values]
    P_GWarray = [row["PIXFOVPROB"] for row in first_values]

    SuggestedPointings = Table(
        [RAarray, DECarray, P_GWarray], names=["RA[deg]", "DEC[deg]", "PGW"]
    )

    return SuggestedPointings


def GetBestTiles3D(skymap, nameEvent, PointingFile, galFile, obsparameters, dirName):
    """
    Compute the best observation tiles based on a 3D method for space-based observatories, considering galaxy probabilities and satellite constraints.
    """

    random.seed()
    RAarray = []
    DECarray = []
    pixlist = []
    pixlistHR = []
    ObsName = []
    Occultedpixels = []

    obspar = obsparameters[0]
    radius = obspar.FOV
    HRnside, reducedNside = GetBestNSIDE(obspar.reducedNside, obspar.HRnside, radius)

    maxRuns = obspar.maxRuns

    # Retrieve maps
    prob = skymap.getMap("prob", reducedNside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, reducedNside
    )
    radecs = co.SkyCoord(rapix, decpix, frame="icrs", unit=(u.deg, u.deg))
    maxRuns = obspar.maxRuns

    doPlot = obspar.doPlot

    # load galaxy catalogue
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        sum_dP_dV = cat["dp_dV"].sum()
    else:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        tGals0 = MangroveGalaxiesProbabilities(tGals0)
        sum_dP_dV = cat["dp_dV"].sum()

    # Add observed pixels to pixlist
    if PointingFile is not None:
        print(PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist)
        pixlist, pixlistHR, sumPGW, doneObs = SubstractPointings2D(
            PointingFile, prob, obspar, pixlist, pixlistHR
        )

        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(f"Total GW probability already covered: {sumPGW}")
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print("========")

    ipix = TransformRADecToPix(radecs, reducedNside)
    newpix = ipix

    # CONVERTING newpix to angles on the coordinate grid
    pixradec = TransformPixToRaDec(newpix, reducedNside)

    first_values = GetBestGridPos3D(
        prob,
        tGals0,
        pixradec,
        newpix,
        radius,
        sum_dP_dV,
        HRnside,
        obspar.numberSides,
        maxRuns,
        doPlot,
        dirName,
        reducedNside,
        Occultedpixels,
        obspar.minProbcut,
    )

    ObsName = [obspar.obs_name for j in range(len(first_values))]
    RAarray = [row["PIXRA"] for row in first_values]
    DECarray = [row["PIXDEC"] for row in first_values]
    P_Galarray = [row["PIXFOVPROB"] for row in first_values]

    SuggestedPointings = Table(
        [ObsName, RAarray, DECarray, P_Galarray],
        names=["ObsName", "RA[deg]", "DEC[deg]", "PGal"],
    )
    return SuggestedPointings


def PGWinFoV_Space_NObs(
    skymap, nameEvent, ObservationTime0, PointingFile, obsparameters, dirName
):
    """
    Compute an observation schedule for space-based observatories using a 3D method.

    It calculates optimal observation times considering the galaxy catalog, GW probability map, and satellite constraints.
    The function returns suggested pointings along with satellite visibility times and SAA status.
    """

    random.seed()
    RAarray = []
    DECarray = []
    pixlist = []
    pixlistHR = []
    P_GWarray = []
    ObsName = []
    Occultedpixels = []

    obspar = obsparameters[0]
    radius = obspar.FOV
    HRnside, reducedNside = GetBestNSIDE(obspar.reducedNside, obspar.HRnside, radius)
    #################################################################################################################################################

    # Retrieve maps
    prob = skymap.getMap("prob", reducedNside)
    highres = skymap.getMap("prob", HRnside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, reducedNside
    )
    radecs = co.SkyCoord(rapix, decpix, frame="icrs", unit=(u.deg, u.deg))
    maxRuns = obspar.maxRuns

    doPlot = obspar.doPlot

    # Add observed pixels to pixlist
    if PointingFile is not None:
        print(PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist)
        pixlist, pixlistHR, sumPGW, doneObs = SubstractPointings2D(
            PointingFile, prob, obspar, pixlist, pixlistHR
        )

        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(f"Total GW probability already covered: {sumPGW}")
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )

    ipix = TransformRADecToPix(radecs, reducedNside)
    newpix = ipix

    first_values1 = GetBestGridPos2D(
        prob,
        highres,
        HRnside,
        reducedNside,
        newpix,
        radius,
        maxRuns,
        Occultedpixels,
        doPlot,
        dirName,
        obspar.numberSides,
        pixlistHR,
        obspar.minProbcut,
    )

    # FOR SPACE ########################################################
    # Computing the satellite name
    SatelliteName = GetSatelliteName(obspar.obs_name, obspar.stationsurl)

    saa = np.empty(obspar.maxRuns + 1, dtype=bool)
    SatTimes = np.empty(obspar.maxRuns + 1, dtype=bool)

    # Iterate through time steps within the specified duration -> to be modified
    current_time = ObservationTime0
    start_time = ObservationTime0
    step = int(obspar.duration / obspar.maxRuns)
    step = datetime.timedelta(minutes=step)

    duration = obspar.duration

    SatTimes, saa = SAA_Times(
        duration,
        start_time,
        current_time,
        SatelliteName,
        saa,
        SatTimes,
        step,
        doPlot,
        dirName,
        obspar.datasetDir,
        obspar.SAAThreshold,
    )

    i = 0
    current_time = ObservationTime0
    start_time = ObservationTime0
    AvailablePixPerTime = []
    TestTime = []
    RadecsVsTimes = []
    matching_tables = []
    ProbaTime = []
    while current_time <= start_time + datetime.timedelta(minutes=duration):
        # Need to get a list of highest pixels
        SatelliteTime = GetSatelliteTime(SatelliteName, current_time)
        satellitePosition, satelliteLocation = GetSatellitePositions(
            SatelliteName, SatelliteTime
        )
        ObsBool, yprob, pixlistRROcc = OccultationCut(
            prob,
            reducedNside,
            current_time,
            obspar,
            satellitePosition,
            satelliteLocation,
        )

        # Let's get the list of pixels available at each iteration
        firstvalue1 = first_values1

        matching_rows1 = FindMatchingCoords(1, firstvalue1, pixlistRROcc, reducedNside)
        matching_tables.append(matching_rows1)

        radectime = co.SkyCoord(
            ra=matching_rows1["PIXRA"] * u.deg, dec=matching_rows1["PIXDEC"] * u.deg
        )
        pix_idx = TransformRADecToPix(radectime, reducedNside)
        pix_proba = matching_rows1["PIXFOVPROB"]

        RadecsVsTimes.append(radectime)
        AvailablePixPerTime.append(pix_idx)
        TestTime.append(current_time)
        ProbaTime.append(pix_proba)

        # List of all cculted pixels
        Occultedpixels.append(pixlistRROcc)
        current_time += step
        i += 1

    # WE CAN GET THE LIST OF PIXELS AVAILABLE AT ALL TIMES
    Occultedpixels = [item for sublist in Occultedpixels for item in sublist]
    OldPix = ipix
    searchpix = np.isin(OldPix, Occultedpixels, invert=True)
    newpix = OldPix[searchpix]

    # Find common pixels
    first_values = FindMatchingPixList(newpix, first_values1)

    ObsName = [obspar.obs_name for j in range(len(first_values))]
    RAarray = [row["PIXRA"] for row in first_values]
    DECarray = [row["PIXDEC"] for row in first_values]
    P_GWarray = [row["PIXFOVPROB"] for row in first_values]

    SuggestedPointings = Table(
        [ObsName, RAarray, DECarray, P_GWarray],
        names=["ObsName", "RA(deg)", "DEC(deg)", "PGW"],
    )

    result = {
        "SatTimes": SatTimes,
        "saa": saa,
        "first_values": first_values,
        "first_values1": first_values1,
        "TestTime": TestTime,
        "matching_tables": matching_tables,
        "Occultedpixels": Occultedpixels,
    }

    return SuggestedPointings, result


def PGalinFoV_Space_NObs(
    skymap, nameEvent, ObservationTime0, PointingFile, galFile, obsparameters, dirName
):
    """
    Compute an observation schedule for space-based observatories using a 3D method.

    Called by :func:`tilepy.include.observationschedule.GetSchedule`, this function generates a schedule of pointings for each observatory,
    considering the field of view (FoV), galaxy catalog, and satellite position with occultation constraints.

    Parameters
    ----------
    skymap : SkyMap
        The object storing sky maps.
    nameEvent : str
        The name of the event.
    ObservationTime0 : datetime
        The desired start time for scheduling to begin.
    PointingFile : str, optional
        Path to the text file containing the pointings that have already been performed.
    galFile : str
        Path to the galaxy catalog file.
    obsparameters : list of ObservationParameters
        A list of sets of parameters for each observatory needed to launch the tiling scheduler.
    dirName : str
        Path to the output directory where the schedules and plots will be saved.

    Returns
    -------
    SuggestedPointings : astropy.table.Table
        Table of suggested pointings with their RA, DEC, and galaxy probability.
    SatTimes : numpy.ndarray
        Array of satellite observation times.
    saa : numpy.ndarray
        Array indicating the satellite's South Atlantic Anomaly (SAA) status at each time step.

    """

    random.seed()
    RAarray = []
    DECarray = []
    pixlist = []
    pixlistHR = []
    ObsName = []
    Occultedpixels = []

    obspar = obsparameters[0]
    radius = obspar.FOV
    HRnside, reducedNside = GetBestNSIDE(obspar.reducedNside, obspar.HRnside, radius)
    #################################################################################################################################################

    # Retrieve maps
    prob = skymap.getMap("prob", reducedNside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, reducedNside
    )
    radecs = co.SkyCoord(rapix, decpix, frame="icrs", unit=(u.deg, u.deg))
    maxRuns = obspar.maxRuns

    doPlot = obspar.doPlot

    # load galaxy catalogue
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        sum_dP_dV = cat["dp_dV"].sum()
    else:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        tGals0 = MangroveGalaxiesProbabilities(tGals0)
        sum_dP_dV = cat["dp_dV"].sum()

    # Add observed pixels to pixlist
    if PointingFile is not None:
        print(PointingFile, prob, reducedNside, radius, pixlist)
        pixlist, pixlistHR, sumPGW, doneObs = SubstractPointings2D(
            PointingFile, prob, obspar, pixlist, pixlistHR
        )

        if obspar.countPrevious:
            maxRuns = obspar.maxRuns - doneObs
        print(
            "==========================================================================================="
        )
        print()
        print(f"Total GW probability already covered: {sumPGW}")
        print(
            f"Count Previous = {obspar.countPrevious}, Number of pointings already done: {doneObs}, "
            f"Max Runs was {obspar.maxRuns}, now is {maxRuns}"
        )
        print(
            "==========================================================================================="
        )

    ipix = TransformRADecToPix(radecs, reducedNside)
    newpix = ipix
    pixradec = radecs

    first_values1 = GetBestGridPos3D(
        prob,
        tGals0,
        pixradec,
        newpix,
        radius,
        sum_dP_dV,
        HRnside,
        obspar.numberSides,
        maxRuns,
        doPlot,
        dirName,
        reducedNside,
        Occultedpixels,
        obspar.minProbcut,
    )

    # FOR SPACE ######################################################
    # Computing the satellite name
    SatelliteName = GetSatelliteName(obspar.obs_name, obspar.stationsurl)

    saa = np.empty(obspar.maxRuns + 1, dtype=bool)
    SatTimes = np.empty(obspar.maxRuns + 1, dtype=bool)

    # Iterate through time steps within the specified duration
    current_time = ObservationTime0
    start_time = ObservationTime0
    step = int(obspar.duration / obspar.maxRuns)
    step = datetime.timedelta(minutes=step)

    duration = obspar.duration

    SatTimes, saa = SAA_Times(
        duration,
        start_time,
        current_time,
        SatelliteName,
        saa,
        SatTimes,
        step,
        doPlot,
        dirName,
        obspar.datasetDir,
        obspar.SAAThreshold,
    )

    i = 0
    current_time = ObservationTime0
    start_time = ObservationTime0
    AvailablePixPerTime = []
    TestTime = []
    RadecsVsTimes = []
    matching_tables = []
    ProbaTime = []
    while current_time <= start_time + datetime.timedelta(minutes=duration):
        # Need to get a list of highest pixels
        SatelliteTime = GetSatelliteTime(SatelliteName, current_time)
        satellitePosition, satelliteLocation = GetSatellitePositions(
            SatelliteName, SatelliteTime
        )
        ObsBool, yprob, pixlistRROcc = OccultationCut(
            prob,
            reducedNside,
            current_time,
            obspar,
            satellitePosition,
            satelliteLocation,
        )

        # Let's get the list of pixels available at each iteration
        firstvalue1 = first_values1

        matching_rows1 = FindMatchingCoords(1, firstvalue1, pixlistRROcc, reducedNside)
        matching_tables.append(matching_rows1)

        radectime = co.SkyCoord(
            ra=matching_rows1["PIXRA"] * u.deg, dec=matching_rows1["PIXDEC"] * u.deg
        )
        theta = np.radians(90.0 - matching_rows1["PIXDEC"])
        phi = np.radians(matching_rows1["PIXRA"])  # phi = longitude
        pix_idx = hp.ang2pix(reducedNside, theta, phi, nest=False)

        pix_proba = matching_rows1["PIXFOVPROB"]

        RadecsVsTimes.append(radectime)
        AvailablePixPerTime.append(pix_idx)
        TestTime.append(current_time)
        ProbaTime.append(pix_proba)

        # List of all cculted pixels
        Occultedpixels.append(pixlistRROcc)
        current_time += step
        i += 1

    # WE CAN GET THE LIST OF PIXELS AVAILABLE AT ALL TIMES --> here we are getting them for all the 90% region... we can only get then for furst value if we want
    Occultedpixels = [item for sublist in Occultedpixels for item in sublist]
    OldPix = ipix
    searchpix = np.isin(OldPix, Occultedpixels, invert=True)
    newpix = OldPix[searchpix]

    # CONVERTING newpix to angles on the coordinate grid
    pixradec = TransformPixToRaDec(newpix, reducedNside)

    # Finding the common radec betweem visible pixels and the grid
    first_values_coords = co.SkyCoord(
        ra=first_values1["PIXRA"], dec=first_values1["PIXDEC"], unit="deg"
    )

    matching_rows = []
    for coord in pixradec:
        sep = first_values_coords.separation(coord)
        matches = np.where(sep < 1e-2 * u.deg)[0]  # adjust tolerance as needed

        matching_rows.extend(first_values1[matches])

    if matching_rows:
        first_values = Table(rows=matching_rows, names=first_values1.colnames)
    else:
        print("No coordinates matched within the tolerance.")
        first_values = Table(names=first_values1.colnames)

    # FOR TARGETED HERE TRY TO FIND OUT WHICH GALAXIES ARE IN THE VISIBLE PART. Then choose the highest 10 betwee nthem

    ObsName = [obspar.obs_name for j in range(len(first_values))]
    RAarray = [row["PIXRA"] for row in first_values]
    DECarray = [row["PIXDEC"] for row in first_values]
    P_Galarray = [row["PIXFOVPROB"] for row in first_values]

    SuggestedPointings = Table(
        [ObsName, RAarray, DECarray, P_Galarray],
        names=["ObsName", "RA[deg]", "DEC[deg]", "PGal"],
    )

    result = {
        "SatTimes": SatTimes,
        "saa": saa,
        "first_values": first_values,
        "first_values1": first_values1,
        "TestTime": TestTime,
        "matching_tables": matching_tables,
        "Occultedpixels": Occultedpixels,
    }

    return SuggestedPointings, result
