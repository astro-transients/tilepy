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
import numpy as np
import pytz
import six
from astropy import units as u
from astropy.table import Table
from six.moves import configparser

from .PointingTools import (NightDarkObservation,
                            NightDarkObservationwithGreyTime, Get90RegionPixReduced, ZenithAngleCut,
                            ComputeProbability2D,
                            FulfillsRequirement, VisibleAtTime, LoadGalaxies, SubstractPointings2D, Tools,
                            LoadGalaxies_SteMgal, SubstractPointings,
                            ModifyCatalogue, FulfillsRequirementGreyObservations, ComputeProbPGALIntegrateFoV,
                            ComputeProbGalTargeted,
                            NextWindowTools, FilterGalaxies, MangroveGalaxiesProbabilities)

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


utc = pytz.UTC
############################################

#              General definitions              #

############################################

def PGWinFoV(skymap, eventName, obspar, dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule based on a 2D method.  
    
    :param skymap: The object containing the skympas
    :type skymap: SkyMap
    :param eventName: The name of the event
    :type eventName: str
    :param ObservationTime0: the desired time for scheduling to start 
    :type ObservationTime0: str
    :param PointingFile: The path to the text file containing the pointings that have already been performed before the scheduling
    :type PointingFile: str 
    :param obspar: Class containing the telescope configuration parameters to be used in the scheduling
    :type obspar: Observation parameters class
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :type dirName: str
    :return: SuggestedPointings, cat
    rtype: ascii table, astropy table
    """

    ObservationTime0 = obspar.obsTime
    PointingFile = obspar.pointingsFile 
    # Main parameters

    print(obspar)

    random.seed()

    RAarray = []
    DECarray = []
    pixlist = []
    ipixlistHR = []
    pixlist1 = []
    ipixlistHR1 = []
    P_GWarray = []
    ObservationTimearray = []
    Round = []
    Fov_obs = []
    ObsName = []
    Duration = []

    print()
    print('-------------------   NEW EVENT   --------------------')
    print()

    # Retrieve maps
    nside = obspar.reducedNside
    prob = skymap.getMap('prob', obspar.reducedNside)
    highres = skymap.getMap('prob', obspar.HRnside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='icrs', unit=(u.deg, u.deg))

    # Add observed pixels to pixlist
    if (PointingFile != None):
        pixlist, P_GW = SubstractPointings2D(
            PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist)
        print('Already observed probability =', P_GW)

    #######################################################

    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    if (obspar.useGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(
            ObservationTime0, obspar)

    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0, obspar)

    counter = 0
    for j, NightDarkRun in enumerate(NightDarkRuns):
        if (len(ObservationTimearray) < obspar.maxRuns):
            ObservationTime = NightDarkRun
            ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, obspar.minProbcut,
                                            obspar.maxZenith, obspar.location, obspar.minMoonSourceSeparation, obspar.useGreytime)
            if ObsBool:
                # Round 1
                P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                     obspar.location, obspar.maxZenith, obspar.FOV, eventName, pixlist, ipixlistHR, counter, dirName, obspar.useGreytime, obspar.doPlot)
                if ((P_GW <= obspar.minProbcut) and obspar.secondRound):
                    # Try Round 2
                    # print('The minimum probability cut being', minProbcut * 100, '% is, unfortunately, not reached.')
                    yprob1 = highres
                    P_GW, TC, pixlist1, ipixlistHR1 = ComputeProbability2D(prob, yprob1, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                           obspar.location, obspar.maxZenith, obspar.FOV, eventName, pixlist1, ipixlistHR1, counter, dirName, obspar.useGreytime, obspar.doPlot)
                    if ((P_GW <= obspar.minProbcut)):
                        print('Tile Pgw= ',P_GW,' is smaller than the minProbCut (',obspar.minProbcut,') => skip this tile')
                    else:
                        Round.append(2)
                        P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                        RAarray.append(
                            float('{:3.4f}'.format(float(TC.ra.deg))))
                        DECarray.append(
                            float('{:3.4f}'.format(float(TC.dec.deg))))
                        ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                        ObsName.append(obspar.name)
                        Duration.append(obspar.duration)
                        Fov_obs.append(obspar.FOV)
                        counter = counter+1
                elif (P_GW >= obspar.minProbcut):
                    Round.append(1)
                    P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                    RAarray.append(float('{:3.4f}'.format(float(TC.ra.deg))))
                    DECarray.append(float('{:3.4f}'.format(float(TC.dec.deg))))
                    ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                    ObsName.append(obspar.name)
                    Duration.append(obspar.duration)
                    Fov_obs.append(obspar.FOV)
                    counter = counter+1
        else:
            break

    print()
    # print("===========================================================================================")

    # print()
    print("Total GW probability covered: ", float('{:1.4f}'.format(float(sum(P_GWarray)))), "Number of runs that fulfill darkness condition  :",
          len(NightDarkRuns), "Number of effective pointings: ", len(ObservationTimearray))

    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, Round, ObsName, Duration, Fov_obs], names=[
                               'Time[UTC]', 'RA[deg]', 'DEC[deg]', 'PGW', 'Round', 'ObsName', 'Duration', 'FoV'])


    print()
    print("================================= Tiling found =============================================")
    print(SuggestedPointings)
    print("===========================================================================================")
    print()
    print('The total probability PGW: ', np.sum(P_GWarray))
    return (SuggestedPointings, ObservationTime0)


def PGalinFoV(skymap, nameEvent, galFile,obspar,dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule based on a 3D method. 
    Depending on the user input in the configuration file and the telescope FoV the pointings use the targtet galaxy strategy or integrated galaxy probability strategy.
    
    :param skymap: The object storing sky maps
    :type skymap: SkyMap
    :param nameEvent: The name of the event
    :type nameEvent: str
    :param ObservationTime0: the desired time for scheduling to start 
    :type ObservationTime0: str
    :param PointingFile: The path to the text file containing the pointings that have already been performed before the scheduling
    :type PointingFile: str
    :param galFile: The path to the galaxy catalog
    :type galFile: str   
    :param obspar: Class containing the telescope configuration parameters to be used in the scheduling
    :type obspar: Observation parameters class
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :type dirName: str
    :return: SuggestedPointings, cat
    rtype: ascii table, astropy table
    """

    ObservationTime0 = obspar.obsTime
    PointingFile = obspar.pointingsFile 
    
    # Main Parameters
    print(obspar)

    # load galaxy catalog from local file
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)

    nside = obspar.HRnside
    prob = skymap.getMap('prob', obspar.HRnside)

    if skymap.is3D:
        print("Found a generic map without 3D information")
    else:
        print("Found a 3D reconstruction")

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        sum_dP_dV = cat['dp_dV'].sum()
    else:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat)
        tGals0 = MangroveGalaxiesProbabilities(tGals0)
        sum_dP_dV = cat['dp_dV'].sum()

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []

    #########################
    if (PointingFile == None):
        tGals = tGals0
        print('No pointings were given to be substracted')
    else:

        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal, alreadysumipixarray1 = SubstractPointings(
            PointingFile, tGals0, alreadysumipixarray1, sum_dP_dV, obspar.FOV, prob, nside)
        # for second round
        # ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray2 = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,obspar.FOV,prob,nside)
        maxRuns = obspar.maxRuns - len(np.atleast_1d(ra))
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)

        print("===========================================================================================")
        print()
        print(nameEvent, "Total GW probability already covered: ", sumPGW,
            "Total Gal probability already covered: ",
            sumPGAL, "Number of effective pointings already done: ", len(np.atleast_1d(ra)))

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
    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    print('maxRuns: ', obspar.maxRuns, 'MinimumProbCutForCatalogue: ',
          obspar.minimumProbCutForCatalogue)

    if (obspar.useGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(
            ObservationTime0, obspar)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0, obspar)

    totalProb = 0.
    counter = 0
    if(obspar.strategy == 'integrated'):
        for j, NightDarkRun in enumerate(NightDarkRuns):
            if (len(ObservationTimearray) < obspar.maxRuns):
                ObservationTime = NightDarkRun
                visible, altaz, tGals_aux = VisibleAtTime(
                    ObservationTime, tGals_aux, obspar.maxZenith, obspar.location)
                if (visible):
                    # select galaxies within the slightly enlarged visiblity window
                    visiMask = altaz.alt.value > 90 - (obspar.maxZenith+obspar.FOV)
                    visiGals = tGals_aux[visiMask]
                    visiGals = ModifyCatalogue(prob, visiGals, obspar.FOV, sum_dP_dV, nside)
                    mask, minz = FulfillsRequirement(
                        visiGals, obspar, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]
                    if (finalGals['dp_dV_FOV'][0] > obspar.minProbcut):
                        # final galaxies within the FoV
                        # print(f"Condition met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} is greater than {obspar.minProbcut}")
                        if ((finalGals['dp_dV_FOV'][0] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
                            visible, altaz, tGals_aux2 = VisibleAtTime(
                                ObservationTime, tGals_aux2, obspar.maxZenith, obspar.location)
                            if (visible):
                                visiMask = altaz.alt.value > 90 - \
                                    (obspar.maxZenith + obspar.FOV)
                                visiGals2 = tGals_aux2[visiMask]
                                visiGals2 = ModifyCatalogue(prob, visiGals2, obspar.FOV, sum_dP_dV, nside)

                                mask, minz = FulfillsRequirement(
                                    visiGals2, obspar, UsePix=False)

                                if obspar.useGreytime:
                                    maskgrey = FulfillsRequirementGreyObservations(
                                        ObservationTime, visiGals2, obspar.location, obspar.minMoonSourceSeparation)
                                    finalGals2 = visiGals2[mask & maskgrey]
                                if not obspar.useGreytime:
                                    finalGals2 = visiGals2[mask]
                                
                                p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(
                                    prob, ObservationTime, obspar.location, finalGals2, False, visiGals2, tGals_aux2, sum_dP_dV, alreadysumipixarray2, nside, minz,obspar, counter, nameEvent, dirName, obspar.doPlot)

                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals2['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals2['DEJ2000'][:1]))))
                                Round.append(2)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                ObsName.append(obspar.name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1

                            else: 
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(
                                    prob, ObservationTime, obspar.location, finalGals, False, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, nameEvent, dirName, obspar.doPlot)
                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals['DEJ2000'][:1]))))
                                Round.append(1)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                ObsName.append(obspar.name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1
                        else:
                            # print("We are in round 1")
                            # print("\n=================================")
                            # print("TARGET COORDINATES AND DETAILS...")
                            # print("=================================")
                            p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(
                                prob, ObservationTime, obspar.location, finalGals, False, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, nameEvent, dirName, obspar.doPlot)

                            RAarray.append(float('{:3.4f}'.format(
                                float(finalGals['RAJ2000'][:1]))))
                            DECarray.append(float('{:3.4f}'.format(
                                float(finalGals['DEJ2000'][:1]))))
                            Round.append(1)
                            P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                            P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                            ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                            ObsName.append(obspar.name)
                            Duration.append(obspar.duration)
                            Fov_obs.append(obspar.FOV)
                            counter = counter + 1

                    else:
                        print(f"Condition not met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} must be greater than {obspar.minProbcut}")
            else:
                break
    if(obspar.strategy == 'targeted'):
        for j, NightDarkRun in enumerate(NightDarkRuns):
            if (len(ObservationTimearray) < obspar.maxRuns):
                ObservationTime = NightDarkRun
                visible, altaz, tGals_aux = VisibleAtTime(
                    ObservationTime, tGals_aux, obspar.maxZenith, obspar.location)
                if (visible):
                    # select galaxies within the slightly enlarged visiblity window
                    visiMask = altaz.alt.value > 90 - \
                        (obspar.maxZenith+obspar.FOV)
                    visiGals = tGals_aux[visiMask]
                    mask, minz = FulfillsRequirement(visiGals, obspar, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]
                    #print('finalGals', finalGals,tGals['dp_dV'][:1]*obspar.minProbcut)
                    if (finalGals['dp_dV'][0] > tGals['dp_dV'][0]*obspar.minProbcut):
                        # print(f"Condition met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV'][0]} is greater than {obspar.minProbcut}")
                        if ((finalGals['dp_dV'][0] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
                            visible, altaz, tGals_aux2 = VisibleAtTime(
                                ObservationTime, tGals_aux2, obspar.maxZenith, obspar.location)
                            if (visible):
                                visiMask = altaz.alt.value > 90 - \
                                    (obspar.maxZenith + obspar.FOV)
                                visiGals2 = tGals_aux2[visiMask]
                                mask, minz = FulfillsRequirement(
                                    visiGals2, obspar, UsePix=False)

                                if obspar.useGreytime:
                                    maskgrey = FulfillsRequirementGreyObservations(
                                        ObservationTime, visiGals2, obspar.location, obspar.minMoonSourceSeparation)
                                    finalGals2 = visiGals2[mask & maskgrey]
                                if not obspar.useGreytime:
                                    finalGals2 = visiGals2[mask]
                                p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbGalTargeted(
                                    prob, ObservationTime, finalGals2, visiGals2, tGals_aux2, sum_dP_dV, alreadysumipixarray2, nside, minz, obspar, counter, dirName)

                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals2['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals2['DEJ2000'][:1]))))
                                Round.append(2)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                ObsName.append(obspar.name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1

                            else:  
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbGalTargeted(
                                    prob, ObservationTime, finalGals, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz,obspar,counter, dirName)
                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals['DEJ2000'][:1]))))
                                Round.append(1)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                ObsName.append(obspar.name)
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                                counter = counter + 1
                        else:
                            #print("We are in round 1")
                            p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbGalTargeted(
                                    prob, ObservationTime, finalGals, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, dirName)
                            RAarray.append(float('{:3.4f}'.format(
                                float(finalGals['RAJ2000'][:1]))))
                            DECarray.append(float('{:3.4f}'.format(
                                float(finalGals['DEJ2000'][:1]))))
                            Round.append(1)
                            P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                            P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                            ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                            ObsName.append(obspar.name)
                            Duration.append(obspar.duration)
                            Fov_obs.append(obspar.FOV)
                            counter = counter + 1

                    else:
                        print(f"Condition not met at {ObservationTime}: dp/dV = {finalGals['dp_dV'][0]} must be greater than {obspar.minProbcut}")

            else:
                break
        
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, Round, ObsName, Duration, Fov_obs], names=[
                               'Time[UTC]', 'RA[deg]', 'DEC[deg]', 'PGW', 'Pgal', 'Round', 'ObsName', 'Duration', 'FoV'])
    print()
    print()
    print('The total probability PGal: ', np.sum(P_GALarray))
    print('The total probability PGW: ', np.sum(P_GWarray))
    return SuggestedPointings, tGals0

def ObservationStartperObs(obsparameters, ObservationTime0):
    '''
    Mid-level function that is called by Nobs Tiling functions for multiple telescopes to find the first observation time for each observatory involved.  

    :param obsparameters: a list of sets of parameters for each observatory needed to launch the tiling scheduler
    :type obsparameters: list of class ObservationParameters
    :param ObservationTime0: the desired time for scheduling to start 
    :type ObservationTime0: str
    :return: obs_time, SameNight, NewActiveObs, NewActiveObsStart
    rtype: datetime, boolean, list, list
    '''


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

        dark_at_start = False

        if obsparameters[j].useGreytime:
            dark_at_start = Tools.CheckWindowGrey(obs_time, obsparameters[j])
        if not obsparameters[j].useGreytime:
            dark_at_start = Tools.CheckWindow(obs_time, obsparameters[j])
        FirstDark[j] = dark_at_start

        # THIS WILL CREATE A DATETIME OBJECT WITH IN THE FORM XX+00:00 WITH NO DOTS
        if FirstDark[j] == True:
            FirstDark_Flag[j] = True
            if obs_time.tzinfo is None:
                obs_time = utc.localize(obs_time)
            ObsFirstTime.append(obs_time)
        else:  # THIS WILL CREATE A DATETIME OBJECT WITH IN THE FORM .XX+00:00
            if obsparameters[j].useGreytime:
                ObsFirstTime1 = NextWindowTools.NextObservationWindowGrey(
                    time=obs_time, obspar=obsparameters[j])
                ObsFirstTime.append(ObsFirstTime1)
            if not obsparameters[j].useGreytime:
                ObsFirstTime1 = NextWindowTools.NextObservationWindow(
                    time=obs_time, obspar=obsparameters[j])
                ObsFirstTime.append(ObsFirstTime1)
            if ObsFirstTime1 != False:
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


def PGWinFoV_NObs(skymap, nameEvent, ObservationTime0, PointingFile, obsparameters, dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule for multiple telescopes/observartories based on a 2D method.   

    :param skymap: The object storing sky maps
    :type skymap: SkyMap
    :param nameEvent: The name of the event
    :type nameEvent: str
    :param ObservationTime0: the desired time for scheduling to start 
    :type ObservationTime0: str
    :param PointingFile: The path to the text file containing the pointings that have already been performed before the scheduling
    :type PointingFile: str
    :param obsparameters: a list of sets of parameters for each observatory needed to launch the tiling scheduler
    :type obsparameters: list of class ObservationParameters
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :type dirName: str
    :return: SuggestedPointings, cat
    rtype: ascii table, astropy table
    """
     
    obs_time, SameNight, NewActiveObs, NewActiveObsStart = ObservationStartperObs(obsparameters, ObservationTime0)


    # START
#################################################################################################################################################
    random.seed()
    RAarray = []
    DECarray = []
    pixlist = []
    ipixlistHR = []
    pixlist1 = []
    ipixlistHR1 = []
    P_GWarray = []
    ObservationTimearray = []
    Round = []
    ObsName = []
    Duration = []
    Fov_obs = []
#################################################################################################################################################
    obspar = obsparameters[0]

    # Retrieve maps
    nside = obspar.reducedNside
    prob = skymap.getMap('prob', obspar.reducedNside)
    highres = skymap.getMap('prob', obspar.HRnside)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='icrs', unit=(u.deg, u.deg))
    # Add observed pixels to pixlist
    if (PointingFile != None):
        print(PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist)
        pixlist, P_GW = SubstractPointings2D(
            PointingFile, prob, obspar.reducedNside, obspar.FOV, pixlist)
        print('Already observed probability =', P_GW)
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
    while (i < 500) & any(SameNight):
        for j in range(len(NewActiveObs)):
            obspar = NewActiveObs[j]
            # print(j)
            # print(NewActiveObs[0].name)
            # print(obspar.name)
            ObservationTime = NewActiveObsTime[j]
            if ITERATION_OBS == len(obsparameters):
                TIME_MIN_ALL = []
                ITERATION_OBS = 0
            ITERATION_OBS += 1
            
            if(couter_per_obs[j] >=  obspar.maxRuns):
                SameNight[j] = False
            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
                ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, obspar.minProbcut,
                                            obspar.maxZenith, obspar.location, obspar.minMoonSourceSeparation, obspar.useGreytime)
                if ObsBool:
                    # Round 1
                    P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                         obspar.location, obspar.maxZenith, obspar.FOV, nameEvent, pixlist, ipixlistHR, counter, dirName, obspar.useGreytime, obspar.doPlot)
                    # print(P_GW, obspar.name)
                    if ((P_GW <= obspar.minProbcut) and obspar.secondRound):
                        # Try Round 2
                        # print('The minimum probability cut being', minProbcut * 100, '% is, unfortunately, not reached.')
                        yprob1 = highres
                        P_GW, TC, pixlist1, ipixlistHR1 = ComputeProbability2D(prob, yprob1, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                               obspar.location, obspar.maxZenith, obspar.FOV, nameEvent, pixlist1, ipixlistHR1, counter, dirName, obspar.useGreytime, obspar.doPlot)
                        if ((P_GW <= obspar.minProbcut)):
                            print('Tile Pgw= ', P_GW, ' is smaller than the minProbCut (',
                                  obspar.minProbCut, ') => skip this tile')
                        else:
                            Round.append(2)
                            P_GWarray.append(
                                float('{:1.4f}'.format(float(P_GW))))
                            RAarray.append(
                                float('{:3.4f}'.format(float(TC.ra.deg))))
                            DECarray.append(
                                float('{:3.4f}'.format(float(TC.dec.deg))))
                            ObservationTime = str(
                                ObservationTime).split('+')[0]
                            try:
                                ObservationTime = datetime.datetime.strptime(
                                    ObservationTime, '%Y-%m-%d %H:%M:%S')
                            except ValueError:
                                ObservationTime = str(
                                    ObservationTime).split('.')[0]
                                ObservationTime = datetime.datetime.strptime(
                                    ObservationTime, '%Y-%m-%d %H:%M:%S')
                            ObservationTimearray.append(ObservationTime)
                            ObsName.append(obspar.name)
                            counter = counter + 1
                            couter_per_obs[j] += 1
                            Duration.append(obspar.duration)
                            Fov_obs.append(obspar.FOV)

                    elif (P_GW >= obspar.minProbcut):
                        Round.append(1)
                        P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                        RAarray.append(
                            float('{:3.4f}'.format(float(TC.ra.deg))))
                        DECarray.append(
                            float('{:3.4f}'.format(float(TC.dec.deg))))
                        ObservationTime = str(ObservationTime).split('+')[0]
                        try:
                            ObservationTime = datetime.datetime.strptime(
                                ObservationTime, '%Y-%m-%d %H:%M:%S')
                        except ValueError:
                            ObservationTime = str(
                                ObservationTime).split('.')[0]
                            ObservationTime = datetime.datetime.strptime(
                                ObservationTime, '%Y-%m-%d %H:%M:%S')
                        ObservationTimearray.append(ObservationTime)
                        ObsName.append(obspar.name)
                        counter = counter+1
                        couter_per_obs[j] += 1
                        Duration.append(obspar.duration)
                        Fov_obs.append(obspar.FOV)

                # HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
                NewActiveObsTime[j] = NewActiveObsTime[j] + \
                    datetime.timedelta(minutes=obspar.duration)

                # HERE WE DETERMINE IF WE ARE STILL IN THE SAME NIGHT FOR THIS OBSERVATORY
                # if (NewActiveObsTime[j] > Tools.NextSunrise(obsstart, obspar)) | (obsstart > Tools.NextMoonrise(obsstart, obspar)):

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

    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, Round, ObsName, Duration, Fov_obs], names=[
                               'Time[UTC]', 'RA(deg)', 'DEC(deg)', 'PGW', 'Round', 'ObsName', 'Duration', 'FoV'])
    print('The total probability PGW: ', np.sum(P_GWarray))

    return SuggestedPointings, obsparameters

def PGalinFoV_NObs(skymap, nameEvent, ObservationTime0, PointingFile, galFile, obsparameters, dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule for multiple telescopes/observartories based on a 3D method.  
     
    :param skymap: The object storing skympas
    :type skymap: SkyMap
    :param ObservationTime0: the desired time for scheduling to start 
    :type ObservationTime0: str
    :param PointingFile: The path to the text file containing the pointings that have already been performed before the scheduling
    :type PointingFile: str
    :param galFile: The path to the galaxy catalog
    :type galFile: str   
    :param obsparameters: a list of sets of parameters for each observatory needed to launch the tiling scheduler
    :type obsparameters: list of class ObservationParameters
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :type dirName: str
    :return: SuggestedPointings, cat, obsparameters
    rtype: ascii table, astropy table, list
    """

    obs_time, SameNight, NewActiveObs, NewActiveObsStart = ObservationStartperObs(obsparameters, ObservationTime0)
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
    has3D = True
    Duration = []
    Fov_obs = []
#################################################################################################################################################
    obspar = obsparameters[0]

    nside = obspar.HRnside
    prob = skymap.getMap('prob', obspar.HRnside)

    # load galaxy catalogue
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
        sum_dP_dV = cat['dp_dV'].sum()
    else:
        cat = skymap.computeGalaxyProbability(cat)
        tGals0 = FilterGalaxies(cat)
        tGals0 = MangroveGalaxiesProbabilities(tGals0)
        sum_dP_dV = cat['dp_dV'].sum()

    # Add observed pixels to pixlist
    if (PointingFile == None):
        tGals = tGals0
        print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal, alreadysumipixarray1 = SubstractPointings(
            PointingFile, tGals0, alreadysumipixarray1, sum_dP_dV, obspar.FOV, prob, nside)
        maxRuns = obspar.maxRuns - len(ra)
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)

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

    totalProb = 0.
    counter = 0
    # print(SameNight)
    # print(NewActiveObs[0].name, NewActiveObs[1].name, NewActiveObs[2].name)
    # print(NewActiveObsTime)
    i = 0
    couter_per_obs = np.zeros(len(NewActiveObs))
    while (i < 5000) & any(SameNight):
        for j in range(len(NewActiveObs)):
            obspar = NewActiveObs[j]
            ObservationTime = NewActiveObsTime[j]
            if ITERATION_OBS == len(obsparameters):
                TIME_MIN_ALL = []
                ITERATION_OBS = 0

            ITERATION_OBS += 1
            if(couter_per_obs[j] >=  obspar.maxRuns):
                SameNight[j] = False
            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:

                if(obspar.strategy == 'integrated'):
                    visible, altaz, tGals_aux = VisibleAtTime(
                        ObservationTime, tGals_aux, obspar.maxZenith, obspar.location)

                    if (visible):

                        # select galaxies within the slightly enlarged visiblity window
                        visiMask = altaz.alt.value > 90 - \
                            (obspar.maxZenith+obspar.FOV)
                        visiGals = tGals_aux[visiMask]
                        visiGals = ModifyCatalogue(
                            prob, visiGals, obspar.FOV, sum_dP_dV, nside)

                        mask, minz = FulfillsRequirement(
                            visiGals, obspar, UsePix=False)
                        if obspar.useGreytime:
                            maskgrey = FulfillsRequirementGreyObservations(
                                ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                            finalGals = visiGals[mask & maskgrey]
                        if not obspar.useGreytime:
                            finalGals = visiGals[mask]

                        if (finalGals['dp_dV_FOV'][0] > obspar.minProbcut):
                            # final galaxies within the FoV
                            #print(f"Condition met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} is greater than {obspar.minProbcut}")
                            if ((finalGals['dp_dV_FOV'][0] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
                                print('probability', finalGals['dp_dV_FOV'][:1])
                                visible, altaz, tGals_aux2 = VisibleAtTime(
                                    ObservationTime, tGals_aux2, obspar.maxZenith, obspar.location)
                                if (visible):
                                    visiMask = altaz.alt.value > 90 - \
                                        (obspar.maxZenith + obspar.FOV)
                                    visiGals2 = tGals_aux2[visiMask]
                                    visiGals2 = ModifyCatalogue(
                                        prob, visiGals2, obspar.FOV, sum_dP_dV, nside)

                                    mask, minz = FulfillsRequirement(
                                        visiGals2, obspar, UsePix=False)

                                    if obspar.useGreytime:
                                        maskgrey = FulfillsRequirementGreyObservations(
                                            ObservationTime, visiGals2, obspar.location, obspar.minMoonSourceSeparation)
                                        finalGals2 = visiGals2[mask & maskgrey]
                                    if not obspar.useGreytime:
                                        finalGals2 = visiGals2[mask]

                                    p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(
                                        prob, ObservationTime, obspar.location, finalGals2, False, visiGals2, tGals_aux2, sum_dP_dV, alreadysumipixarray2, nside, minz, obspar, counter, nameEvent, dirName, obspar.doPlot)
                                    RAarray.append(float('{:3.4f}'.format(
                                        float(finalGals2['RAJ2000'][:1]))))
                                    DECarray.append(float('{:3.4f}'.format(
                                        float(finalGals2['DEJ2000'][:1]))))
                                    Round.append(2)
                                    P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                    P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                    ObservationTime = str(ObservationTime).split('.')[0]
                                    try:
                                        ObservationTime = datetime.datetime.strptime(
                                            ObservationTime, '%Y-%m-%d %H:%M:%S')
                                    except:
                                        ObservationTime = str(
                                            ObservationTime).split('+')[0]
                                        ObservationTime = datetime.datetime.strptime(
                                            ObservationTime, '%Y-%m-%d %H:%M:%S')
                                    ObservationTimearray.append(ObservationTime)
                                    ObsName.append(obspar.name)
                                    counter = counter + 1
                                    couter_per_obs[j] += 1
                                    Duration.append(obspar.duration)
                                    Fov_obs.append(obspar.FOV)

                            else:
                                # print("\n=================================")
                                # print("TARGET COORDINATES AND DETAILS...")
                                # print("=================================")
                                # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(
                                    prob, ObservationTime, obspar.location, finalGals, False, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, nameEvent, dirName, obspar.doPlot)
                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals['DEJ2000'][:1]))))
                                Round.append(1)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTime = str(ObservationTime).split('.')[0]
                                try:
                                    ObservationTime = datetime.datetime.strptime(
                                        ObservationTime, '%Y-%m-%d %H:%M:%S')
                                except:
                                    ObservationTime = str(
                                        ObservationTime).split('+')[0]
                                    ObservationTime = datetime.datetime.strptime(
                                        ObservationTime, '%Y-%m-%d %H:%M:%S')
                                ObservationTimearray.append(ObservationTime)
                                ObsName.append(obspar.name)
                                counter = counter + 1
                                couter_per_obs[j] += 1
                                Duration.append(obspar.duration)
                                Fov_obs.append(obspar.FOV)
                            # ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))

                        else:
                            print(f"Condition NOT met at {ObservationTime}: dp/dV_FOV = {finalGals['dp_dV_FOV'][0]} is greater than {obspar.minProbcut}")


                if(obspar.strategy == 'targeted'):
                        visible, altaz, tGals_aux = VisibleAtTime(
                            ObservationTime, tGals_aux, obspar.maxZenith, obspar.location)
                        if (visible):
                            # select galaxies within the slightly enlarged visiblity window
                            visiMask = altaz.alt.value > 90 - \
                                (obspar.maxZenith+obspar.FOV)
                            visiGals = tGals_aux[visiMask]
                            mask, minz = FulfillsRequirement(
                                visiGals, obspar, UsePix=False)
                            if obspar.useGreytime:
                                maskgrey = FulfillsRequirementGreyObservations(
                                    ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                                finalGals = visiGals[mask & maskgrey]
                            if not obspar.useGreytime:
                                finalGals = visiGals[mask]
                            #print('finalGals', finalGals,tGals['dp_dV'][:1]*obspar.minProbcut)
                            if (finalGals['dp_dV'][0] > tGals['dp_dV'][:1]*obspar.minProbcut):
                                print(f"Condition met at {ObservationTime}: dp/dV = {finalGals['dp_dV'][0]} is greater than {obspar.minProbcut}")
                                if ((finalGals['dp_dV'][0] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
                                    visible, altaz, tGals_aux2 = VisibleAtTime(
                                        ObservationTime, tGals_aux2, obspar.maxZenith, obspar.location)
                                    if (visible):
                                        visiMask = altaz.alt.value > 90 - \
                                            (obspar.maxZenith + obspar.FOV)
                                        visiGals2 = tGals_aux2[visiMask]
                                        mask, minz = FulfillsRequirement(
                                            visiGals2, obspar, UsePix=False)

                                        if obspar.useGreytime:
                                            maskgrey = FulfillsRequirementGreyObservations(
                                                ObservationTime, visiGals2, obspar.location, obspar.minMoonSourceSeparation)
                                            finalGals2 = visiGals2[mask & maskgrey]
                                        if not obspar.useGreytime:
                                            finalGals2 = visiGals2[mask]
                                        p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbGalTargeted(
                                            prob, ObservationTime, finalGals2, visiGals2, tGals_aux2, sum_dP_dV, alreadysumipixarray2, nside, minz, obspar, counter, dirName)

                                        RAarray.append(float('{:3.4f}'.format(
                                            float(finalGals2['RAJ2000'][:1]))))
                                        DECarray.append(float('{:3.4f}'.format(
                                            float(finalGals2['DEJ2000'][:1]))))
                                        Round.append(2)
                                        P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                        P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                        ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                        ObsName.append(obspar.name)
                                        counter = counter + 1
                                        couter_per_obs[j] += 1
                                        Duration.append(obspar.duration)
                                        Fov_obs.append(obspar.FOV)

                                    else:  
                                        p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbGalTargeted(
                                            prob, ObservationTime, finalGals, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz,obspar,counter, dirName)
                                        RAarray.append(float('{:3.4f}'.format(
                                            float(finalGals['RAJ2000'][:1]))))
                                        DECarray.append(float('{:3.4f}'.format(
                                            float(finalGals['DEJ2000'][:1]))))
                                        Round.append(1)
                                        P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                        P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                        ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                        ObsName.append(obspar.name)
                                        counter = counter + 1
                                        couter_per_obs[j] += 1
                                        Duration.append(obspar.duration)
                                        Fov_obs.append(obspar.FOV)
                                else:
                                    #print("We are in round 1")
                                    p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbGalTargeted(
                                            prob, ObservationTime, finalGals, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, dirName)
                                    RAarray.append(float('{:3.4f}'.format(
                                        float(finalGals['RAJ2000'][:1]))))
                                    DECarray.append(float('{:3.4f}'.format(
                                        float(finalGals['DEJ2000'][:1]))))
                                    Round.append(1)
                                    P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                    P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                    ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                    ObsName.append(obspar.name)
                                    counter = counter + 1
                                    couter_per_obs[j] += 1
                                    Duration.append(obspar.duration)
                                    Fov_obs.append(obspar.FOV)

                            else:
                                print(f"Condition NOT met at {ObservationTime}: dp/dV = {finalGals['dp_dV'][0]} is greater than {obspar.minProbcut}")

        


                # HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
                NewActiveObsTime[j] = NewActiveObsTime[j] + \
                    datetime.timedelta(minutes=obspar.duration)

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

    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, Round, ObsName, Duration, Fov_obs], names=[
                               'Time[UTC]', 'RA(deg)', 'DEC(deg)', 'PGW', 'Pgal', 'Round', 'ObsName', 'Duration', 'FoV'])
    print('The total probability PGal: ', np.sum(P_GALarray))
    print('The total probability PGW: ', np.sum(P_GWarray))
    return SuggestedPointings, tGals0, obsparameters
