from .PointingTools import (NightDarkObservation, SelectObservatory_fromHotspot,
                            NightDarkObservationwithGreyTime, LoadHealpixMap, LoadHealpixUNIQMap,
                            Get90RegionPixReduced, ZenithAngleCut, ComputeProbability2D,
                            FulfillsRequirement, VisibleAtTime, LoadGalaxies, CorrelateGalaxies_LVC, SubstractPointings2D, SimpleGWprob, ComputeProbGalTargeted,
                            Tools, LoadGalaxies_SteMgal, CorrelateGalaxies_LVC_SteMass, SubstractPointings,
                            ModifyCatalogue, FulfillsRequirementGreyObservations, ComputeProbPGALIntegrateFoV,ComputeProbGalTargeted,
                            ModifyCataloguePIX, ObservationParameters, NextWindowTools,
                            ComputeProbability2D_SelectClusters, GiveProbToGalaxy, LoadGalaxiesSimulation)
from .Observatories import CTASouthObservatory, CTANorthObservatory
import numpy as np
from astropy import units as u
from astropy.table import Table
import datetime
import healpy as hp
import astropy.coordinates as co
import random
import pytz
from .ObservingTimes import ObtainObservingTimes
from six.moves import configparser
import six
if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


utc = pytz.UTC
############################################

#              General definitions              #

############################################

def PGWinFoV(filename,ObservationTime0,PointingFile,obspar,dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule based on a 2D method.  
    
    :param filename: the name of the probability map file
    :type filename: str
    :param ObservationTime0: the desired time for scheduling to start 
    :type date: str
    :param PointingFile: The text file containing the pointings that have already been performed before the scheduling
    :type datasetDir: str
    :param obspar: Class containing the telescope configuration parameters to be used in the scheduling
    :type galcatname: str
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :return: SuggestedPointings, ObservationTime0
    rtype: ascii table, str
    """
    # Main parameters

    print(obspar)

    # link to the GW map
    name = filename.split('.')[0].split('/')[-1]

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

    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading map from ', filename)
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    prob = hp.pixelfunc.ud_grade(tprob, obspar.reducedNside, power=-2)
    nside = obspar.reducedNside

    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='fk5', unit=(u.deg, u.deg))

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
                                                                     obspar.location, obspar.maxZenith, obspar.FOV, name, pixlist, ipixlistHR, counter, dirName, obspar.useGreytime, obspar.doPlot)
                if ((P_GW <= obspar.minProbcut) and obspar.secondRound):
                    # Try Round 2
                    # print('The minimum probability cut being', minProbcut * 100, '% is, unfortunately, not reached.')
                    yprob1 = highres
                    P_GW, TC, pixlist1, ipixlistHR1 = ComputeProbability2D(prob, yprob1, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                           obspar.location, obspar.maxZenith, obspar.FOV, name, pixlist1, ipixlistHR1, counter, dirName, obspar.useGreytime, obspar.doPlot)
                    if ((P_GW <= obspar.minProbcut)):
                        print('Fail')
                    else:
                        Round.append(2)
                        P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                        RAarray.append(
                            float('{:3.4f}'.format(float(TC.ra.deg))))
                        DECarray.append(
                            float('{:3.4f}'.format(float(TC.dec.deg))))
                        ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                        counter = counter+1
                elif (P_GW >= obspar.minProbcut):
                    Round.append(1)
                    P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                    RAarray.append(float('{:3.4f}'.format(float(TC.ra.deg))))
                    DECarray.append(float('{:3.4f}'.format(float(TC.dec.deg))))
                    ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                    counter = counter+1
        else:
            break

    print()
    # print("===========================================================================================")

    # print()
    print("Total GW probability covered: ", float('{:1.4f}'.format(float(sum(P_GWarray)))), "Number of runs that fulfill darkness condition  :",
          len(NightDarkRuns), "Number of effective pointings: ", len(ObservationTimearray))

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, Round], names=[
                               'Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'PGW', 'Round'])
    return (SuggestedPointings, ObservationTime0)


def PGalinFoV(filename,ObservationTime0,PointingFile,galFile,obspar,dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule based on a 2D method.  
    
    :param filename: the name of the probability map file
    :type filename: str
    :param ObservationTime0: the desired time for scheduling to start 
    :type date: str
    :param PointingFile: The text file containing the pointings that have already been performed before the scheduling
    :type datasetDir: str
    :param galFile: The path to the galaxy catalog
    :type datasetDir: str   
    :param obspar: Class containing the telescope configuration parameters to be used in the scheduling
    :type galcatname: str
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :return: SuggestedPointings, cat
    rtype: ascii table, astropy table
    """

    # Main Parameters
    print(obspar)

    # load galaxy catalog from local file
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)
    print('done loading galaxies')

    name = filename.split('.')[0].split('/')[-1]

    print('Loading map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    has3D = True
    if (len(distnorm) == 0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D = False
    else:
        print("Found a 3D reconstruction")

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(
            prob, distmu, distsigma, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(
            prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)

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
        print(
            name, "Total GW probability already covered: ", sumPGW,
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
                    visiMask = altaz.alt.value > 90 - \
                        (obspar.maxZenith+obspar.FOV)
                    visiGals = tGals_aux[visiMask]
                    visiGals = ModifyCatalogue(
                        prob, visiGals, obspar.FOV, sum_dP_dV, nside)

                    mask, minz = FulfillsRequirement(
                        visiGals, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]

                    if (finalGals['dp_dV_FOV'][:1] > obspar.minProbcut):
                        # final galaxies within the FoV
                        if ((finalGals['dp_dV_FOV'][:1] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
                            visible, altaz, tGals_aux2 = VisibleAtTime(
                                ObservationTime, tGals_aux2, obspar.maxZenith, obspar.location)
                            if (visible):
                                visiMask = altaz.alt.value > 90 - \
                                    (obspar.maxZenith + obspar.FOV)
                                visiGals2 = tGals_aux2[visiMask]
                                visiGals2 = ModifyCatalogue(prob, visiGals2, obspar.FOV, sum_dP_dV, nside)

                                mask, minz = FulfillsRequirement(
                                    visiGals2, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=False)

                                if obspar.useGreytime:
                                    maskgrey = FulfillsRequirementGreyObservations(
                                        ObservationTime, visiGals2, obspar.location, obspar.minMoonSourceSeparation)
                                    finalGals2 = visiGals2[mask & maskgrey]
                                if not obspar.useGreytime:
                                    finalGals2 = visiGals2[mask]
                                
                                p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(
                                    prob, ObservationTime, obspar.location, finalGals2, False, visiGals2, tGals_aux2, sum_dP_dV, alreadysumipixarray2, nside, minz,obspar, counter, name, dirName, obspar.doPlot)

                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals2['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals2['DEJ2000'][:1]))))
                                Round.append(2)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                counter = counter + 1

                            else: 
                                p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(
                                    prob, ObservationTime, obspar.location, finalGals, False, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, name, dirName, obspar.doPlot)
                                RAarray.append(float('{:3.4f}'.format(
                                    float(finalGals['RAJ2000'][:1]))))
                                DECarray.append(float('{:3.4f}'.format(
                                    float(finalGals['DEJ2000'][:1]))))
                                Round.append(1)
                                P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                                P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                                ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                                counter = counter + 1
                        else:
                            # print("We are in round 1")
                            # print("\n=================================")
                            # print("TARGET COORDINATES AND DETAILS...")
                            # print("=================================")
                            p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(
                                prob, ObservationTime, obspar.location, finalGals, False, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar, counter, name, dirName, obspar.doPlot)

                            RAarray.append(float('{:3.4f}'.format(
                                float(finalGals['RAJ2000'][:1]))))
                            DECarray.append(float('{:3.4f}'.format(
                                float(finalGals['DEJ2000'][:1]))))
                            Round.append(1)
                            P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                            P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                            ObservationTimearray.append(ObservationTime.strftime("%Y-%m-%d %H:%M:%S"))
                            counter = counter + 1

                    else:
                        print("NOT passing the cut on dp/dV_FOV > ", obspar.minProbcut)
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
                    mask, minz = FulfillsRequirement(
                        visiGals, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]
                    #print('finalGals', finalGals,tGals['dp_dV'][:1]*obspar.minProbcut)
                    if (finalGals['dp_dV'][:1] > tGals['dp_dV'][:1]*obspar.minProbcut):
                        if ((finalGals['dp_dV'][:1] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
                            visible, altaz, tGals_aux2 = VisibleAtTime(
                                ObservationTime, tGals_aux2, obspar.maxZenith, obspar.location)
                            if (visible):
                                visiMask = altaz.alt.value > 90 - \
                                    (obspar.maxZenith + obspar.FOV)
                                visiGals2 = tGals_aux2[visiMask]
                                mask, minz = FulfillsRequirement(
                                    visiGals2, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=False)

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
                            counter = counter + 1

                    else:
                        print("NOT passing the cut on dp/dV_FOV > ", obspar.minProbcut)
            else:
                break
        
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, Round], names=[
                               'Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'PGW', 'Pgal', 'Round'])
    print()
    print("================================= Tiling found =============================================")
    print(SuggestedPointings)
    print("===========================================================================================")
    print()
    return SuggestedPointings, cat


def PGalinFoV_PixRegion(filename,ObservationTime0,PointingFile,galFile, obspar,dirName):
    """
    Mid-level function that is called by GetSchedule to compute a observation schedule based on a 2D method.  
    
    :param filename: the name of the probability map file
    :type filename: str
    :param ObservationTime0: the desired time for scheduling to start 
    :type date: str
    :param PointingFile: The text file containing the pointings that have already been performed before the scheduling
    :type datasetDir: str
    :param galFile: The path to the galaxy catalog
    :type datasetDir: str   
    :param obspar: Class containing the telescope configuration parameters to be used in the scheduling
    :type galcatname: str
    :param dirName: Path to the output directory where the schedules and plots will eb saved
    :return: SuggestedPointings, cat, area
    rtype: ascii table, astropy table, float
    """


#(filename, ObservationTime, PointingsFile, galaxies, obspar, dirName)
    # Main parameters from config

    ###############################

    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    # galFile='./GLADE_2clean.txt'

    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)
    print('done loading galaxies')

    name = filename.split('.')[0].split('/')[-1]

    # name = "Default"
    # if ('G' in filename):
    #    names = filename.split("_")
    #    name = names[0]

    print()
    print('-------------------   NEW EVENT   --------------------')
    print()

    print('Loading map from ', filename)
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob,obspar.reducedNside,power=-2)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    has3D = True
    if (len(distnorm) == 0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D = False
    else:
        print("Found a 3D reconstruction")

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(tprob, distmu, distsigma, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(tprob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []
    #########################

    if (PointingFile == None):
        tGals = tGals0
        print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal, alreadysumipixarray1 = SubstractPointings(
            PointingFile, tGals0, alreadysumipixarray1, sum_dP_dV, obspar.FOV, prob, nside)
        maxRuns = obspar.maxRuns - len(np.atleast_1d(ra))
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)

        # ObservedPointings = Table([time, ra, dec, AlreadyObservedPgw, AlreadyObservedPgal],names=['Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'Covered GW probability','Pgal covered'])
        # print(ObservedPointings)
        print("===========================================================================================")
        print()
        print(
            name, "Total GW probability already covered: ", sumPGW,
            "Total Gal probability already covered: ",
            sumPGAL, "Number of effective pointings already done: ", len(np.atleast_1d(ra)))
    #########################

    tGals_aux = tGals
    tGals_aux2 = tGals

    P_GALarray = []
    P_GWarray = []
    ObservationTimearray = []
    RAarray = []
    DECarray = []

    # In case one wants to see if a next round would give us better results..
    # So the time will be j-1
    nextround = False
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
    n = 0

    ###############################

    # Get the RA & DEC of pixles of the pixels in an enclosed probability region (% precised by percentageMOC).
    # Reduce these RA DEC to angles in maps with smaller resolution (reducedNside)

    pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob,obspar.PercentCoverage,obspar.reducedNside)


    ##############################
    counter = 0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < obspar.maxRuns):
            ObservationTime = NightDarkRuns[j]
            if (nextround):
                ObservationTime = NightDarkRuns[j - 1]
                nextround = False
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, obspar.maxZenith,
                                                      obspar.location)

            if (visible):

                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - \
                    (obspar.maxZenith + obspar.FOV)
                visiGals = tGals_aux[visiMask]

                mask, minz = FulfillsRequirement(
                    visiGals, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=True)

                finalGals = visiGals[mask]
                visiPix = ModifyCataloguePIX(pix_ra1, pix_dec1, ObservationTime, obspar.maxZenith, tprob, finalGals, obspar.FOV,
                                             sum_dP_dV, nside, obspar.reducedNside, minz,obspar.location)

                if (visiPix['PIXFOVPROB'][:1] > obspar.minProbcut):
                    n = n + 1
                    # final galaxies within the FoV

                    # print("\n=================================")
                    # print("TARGET COORDINATES AND DETAILS...")
                    # print("=================================")
                    # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob, ObservationTime, obspar.location,
                                                                                               visiPix, True, visiGals,
                                                                                               tGals_aux, sum_dP_dV,
                                                                                               alreadysumipixarray1,
                                                                                               nside, minz, obspar, counter, name, dirName, obspar.doPlot)
                    RAarray.append(
                        float('{:3.4f}'.format(float(visiPix['PIXRA'][:1]))))
                    DECarray.append(
                        float('{:3.4f}'.format(float(visiPix['PIXDEC'][:1]))))
                    Round.append(1)
                    P_GALarray.append(float('{:1.4f}'.format(p_gal)))
                    P_GWarray.append(float('{:1.4f}'.format(p_gw)))
                    ObservationTimearray.append(
                        str(ObservationTime).split('.')[0])
                    counter = counter+1

                else:
                    print("\nOptimal pointing position is: ")
                    print(visiPix['PIXRA', 'PIXDEC', 'PIXFOVPROB'][:1])
                    print("NOT passing the cut on dp_dV_FOV > ", obspar.minProbcut)

        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, Round],
                               names=['Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'PGW',
                                      'Pgal', 'Round'], )
    print(SuggestedPointings)
    print("Name", name, "Total GW probability covered: ", float('{:1.4f}'.format(float(sum(P_GWarray)))), "Total Gal probability covered: ", sum(P_GALarray),
          "Number of runs that fulfill darkness condition  :", len(
              NightDarkRuns), "Number of effective pointings: ",
          len(ObservationTimearray))
    return SuggestedPointings, cat, area

    # print("===========================================================================================")
    # print()
    # print()


def PGWonFoV_WindowsfromIRFs(filename, InputChar, TC, parameters, dirName):

    UseObs = InputChar['Observatory']
    run = InputChar['run']
    mergerID = InputChar['MergerID']
    zenith = InputChar['Zenith']
    ObservationTime0 = datetime.datetime.strptime(
        InputChar['Time'], '%Y-%m-%d %H:%M:%S.%f')

    # Main parameters from config

    obspar = ObservationParameters()
    obspar.from_configfile(parameters)
    print(obspar)

    # Observatory
    if UseObs == 'South':
        print('Observed form the', UseObs)
        observatory = CTASouthObservatory()
    else:
        print('Observed from the', UseObs)
        observatory = CTANorthObservatory()

    # link to the GW map
    name = filename.split('.')[0].split('/')[-1]
    # if('G' in filename):
    #    names = filename.split("_")
    #    name= names[0]

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
    PreDefWindow = []
    Duration = []

    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading map from ', filename)
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    prob = hp.pixelfunc.ud_grade(tprob, obspar.reducedNside, power=-2)

    nside = obspar.reducedNside

    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='fk5', unit=(u.deg, u.deg))
    has3D = True
    if (len(distnorm) == 0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D = False
    else:
        print("Found a 3D reconstruction")

    # print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    # Computation of observation times considering Moon and Sun in 24 hours
    # If observations are possible right now or in less than XX seconds, Prompt=True
    # Check wheather the probability is enough to start observations (This can be done in the future)
    # The effective first window is obtained

    # ToDO need to change that to no max integration time (limited by visibility)
    # ToDo need to change the maximum latency to 48 hours
    totalTime = 7200

    followupDelay = 30
    SlewingTime = 210
    total_followupDelay = followupDelay + SlewingTime

    # Else, continue
    # total_followupDelay=0
    # Obtain the times from the comparison to a given spectrum and sensitivities computed by cssens
    # Minimum set to 10 seconds!
    # Starting of observation and duration of observation
    tstar, tobs = ObtainObservingTimes(
        totalTime, total_followupDelay, run, mergerID, observatory, obspar.maxZenith)

    counter = 0
    # Loop where probabilities are computed comes here

    for j in range(0, len(tobs), 1):
        if (j > obspar.maxRuns):
            break
        else:
            ObservationTime = ObservationTime0 + \
                datetime.timedelta(seconds=tstar[j])
            # print(ObservationTime)

            ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, obspar.minProbcut,
                                            obspar.maxZenith, observatory.location, obspar.minMoonSourceSeparation, False)
            # print(ObsBool)
            if ObsBool:
                # Round 1
                # print('Round ',j,'counter',counter)
                P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, obspar.reducedNside, obspar.HRnside,
                                                                     obspar.minProbcut, ObservationTime, observatory.location,
                                                                     obspar.maxZenith, obspar.FOV, name, pixlist, ipixlistHR,
                                                                     counter,
                                                                     dirName, False, obspar.doPlot)
                print('P_GW', P_GW)
                if ((P_GW >= obspar.minProbcut)):
                    P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                    RAarray.append(float('{:3.4f}'.format(float(TC.ra.deg))))
                    DECarray.append(float('{:3.4f}'.format(float(TC.dec.deg))))
                    ObservationTimearray.append(ObservationTime)
                    Duration.append(tobs[j])
                    # PreDefWindow.append(predefWind[j])
                    # ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar)
                    counter = counter + 1
                else:
                    print('Probability too low')

            # else:
            #    break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, Duration],
                               names=['Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'PGW', 'Duration[s]'])

    return (SuggestedPointings, ObservationTime0, obspar.FOV, nside, len(tobs))


def PGWonFoV_WindowOptimisation(filename, timeStr, TC, parameters, conf, datasetDir, outDir):

    # UseObs = InputChar['Observatory']
    UseObs = SelectObservatory_fromHotspot(filename)

    # ID retrieved from the filename
    ID = filename.split('/')[-1].split('.')[0]

    # zenith=InputChar['Zenith']

    ObservationTime0 = datetime.datetime.strptime(
        timeStr.split('.')[0], '%Y-%m-%d %H:%M:%S')
    ObservationTime0 = pytz.utc.localize(ObservationTime0)

    # Main parameters from config

    obspar = ObservationParameters()
    obspar.from_configfile(parameters)

    ##################

    # path = outDir + '/PointingPlotting/' + run + '_' + mergerID + '/EvolutionPlot/'
    path = outDir + '/PointingPlotting/' + ID + '/EvolutionPlot/'
    # Observatory
    if UseObs == 'South':
        print('Observed form the', UseObs)
        obspar.name = UseObs
        obspar.lon = CTASouthObservatory().Lon
        obspar.lat = CTASouthObservatory().Lat
        obspar.height = CTASouthObservatory().Height
        obspar.location = CTASouthObservatory().location
    else:
        print('Observed from the', UseObs)
        obspar.name = UseObs
        obspar.lon = CTANorthObservatory().Lon
        obspar.lat = CTANorthObservatory().Lat
        obspar.height = CTANorthObservatory().Height
        obspar.location = CTANorthObservatory().location

    print(obspar)
    # link to the GW map
    # name = filename.split('.')[0].split('/')[-1]
    # if('G' in filename):
    #    names = filename.split("_")
    #    name= names[0]

    random.seed()

    RAarray = []
    DECarray = []
    Obsarray = []
    pixlist = []
    ipixlistHR = []
    pixlist1 = []
    ipixlistHR1 = []
    P_GWarray = []
    ObservationTimearray = []
    Exposure = []
    Delay = []
    ZenIniarray = []
    ZenEndarray = []
    ObsBoolarray = []

    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading map from ', filename)
    # tprob, distmu, distsigma, distnorm  = LoadHealpixUNIQMap(filename)

    # Get the skymap as an ordered dict
    import ligo.skymap.io.fits as lf
    skymap_OD = lf.read_sky_map(filename)

    # Get the main parameters
    tprob = skymap_OD[0]
    has3D = False
    print('.fits file only has 2D information')
    # distmu = skymap_OD[1]['distmean']
    # distsigma = skymap_OD[1]['diststd']
    # distnorm = skymap_OD[1]['diststd']

    prob = hp.pixelfunc.ud_grade(tprob, obspar.reducedNside, power=-2)
    nside = obspar.reducedNside
    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)

    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='fk5', unit=(u.deg, u.deg))
    has3D = True

    # if (len(distnorm) == 0):
    #    print("Found a generic map without 3D information")
    #    # flag the event for special treatment
    #    has3D = False
    # else:
    #    print("Found a 3D reconstruction")

    # print('----------   NEW FOLLOW-UP ATTEMPT   ----------')

    # Set of delays and slewing times
    # Units: seconds
    totalTime = 172800  # 48h
    followupDelay = 30  # Minimum delay to start observaiton
    SlewingTime = 210  # Slewing to first position
    interObsSlew = 20  # Slewing time between observations
    total_followupDelay = followupDelay + SlewingTime

    # From the injection time, look for the next window. Time is the time of the first observation
    ObservationTime = ObservationTime0 + \
        datetime.timedelta(seconds=total_followupDelay)

    print("Main info of this scheduling:", totalTime,
          total_followupDelay, ID, TC, obspar.location)

    TotalNights = 2
    counter = 0
    counterTotalPossible = 0
    maxRuns = 20
    # GrbsensMax = 16384
    # 15 minutes.  Allows to check if 15 minutes later there is a better opportunity
    AuxTimeNextTry = 900
    # minProbcut = 0.005

    # case 1 = TstartNight == False
    # case 2 = "The source is not on the temporal FoV of the instrument"
    for nights in range(0, TotalNights):
        TstartNight = NextWindowTools.NextObservationWindow(
            ObservationTime, obspar)
        print("TstartNight", TstartNight)
        if (TstartNight == False):
            TendNight = False
            print("Night number", nights, " window starts at",
                  TstartNight, "finishes at", TendNight)
            ObsCase = 'NoDarknessFound'
            # print("===== RESULTS ========")
            print('++++ No Darkness Time Found +++++')
            P_GWarray.append(0)
            RAarray.append(0)
            DECarray.append(0)
            Obsarray.append(0)
            ZenIniarray.append(0)
            ZenEndarray.append(0)
            Exposure.append(0)
            Delay.append(0)
            ObsBoolarray.append(ObsCase)
            break
        else:
            print("Observations can be allocated")
            TendNight = NextWindowTools.EndObservationWindow(
                TstartNight, obspar)
            print("Night number", nights, " window starts at",
                  TstartNight, "finishes at", TendNight)

            # Look for the first time that the C.R. observable is > 5%
            TemporalBin = datetime.timedelta(seconds=600)  # Every 10 minutes
            for i in range(0, 24):
                checkTime = TstartNight + i * TemporalBin
                print(TendNight, TstartNight)
                print(TendNight.tzinfo, TstartNight.tzinfo)
                if ((TendNight - checkTime).total_seconds() < 0):
                    ObsBool = False
                    print("The source is not on the temporal FoV of the instrument")
                    break
                ObsBool, yprob = ZenithAngleCut(prob, nside, checkTime, obspar.minProbcut,
                                                obspar.maxZenith, obspar.location, obspar.minMoonSourceSeparation, useGreytime=False)
                print(checkTime, ObsBool)
                if (ObsBool == True):
                    StartObsTime = checkTime
                    print("The first time the GW is observed is", StartObsTime, "minProbcut", obspar.minProbcut, "maxZenith",
                          obspar.maxZenith)
                    break
            if (ObsBool):
                # Get the first 100 clusters of probability for the masked map
                for tt in range(0, 1000):
                    print("  ")
                    print('++++++++++++++++++++++++++++++++++')
                    print("---- NEW OBSERVATION ATTEMPT ----")
                    print('++++++++++++++++++++++++++++++++++')
                    print("  ")
                    print("Observation number ", counter)
                    print("Starting time", StartObsTime)
                    if (tt > maxRuns):
                        break
                    if ((TendNight - datetime.timedelta(seconds=interObsSlew) - StartObsTime).total_seconds() <= 0):
                        print("---- NIGHT IS OVER ----")
                        break
                    # Total Exposure for the night
                    TotalExposure = int(
                        (TendNight - StartObsTime).total_seconds())
                    DelayObs = int(
                        (StartObsTime - ObservationTime0).total_seconds())

                    print("ObservationTime0", ObservationTime0)
                    print("TotalExposure", TotalExposure)
                    print("DelayObs", DelayObs)
                    print("----------------------------")
                    P_GW, TC, ObsExp, ZenIni, ZenEnd, ObsCase, pixlist, ipixlistHR = ComputeProbability2D_SelectClusters(
                        prob, highres, radecs,  conf, StartObsTime, DelayObs, interObsSlew, obspar, ID, pixlist, ipixlistHR, counter, datasetDir, outDir, False, False)
                    print("=============")
                    print("P_GW, ObsExp, ZenIni, ZenEnd, ObsCase")
                    print(P_GW, ObsExp, ZenIni, ZenEnd, ObsCase)
                    print("=============")
                    counterTotalPossible = counterTotalPossible + 1
                    # print("PGW", P_GW,'COORDINATES:', TC,'ZENITH CHANGE', ZenIni,'->', ZenEnd)
                    # print("TotalExposure: ",TotalExposure,"DelayObs:",DelayObs, "Observation Number:",counter, "Exposure:", ObsExp)

                    ObservationTimearray.append(str(StartObsTime).split('.')[0].split('+')[0])

                    if (ObsCase == 'TimeNotEnoughIte' or ObsCase == 'TimeNotEnough'):
                        StartObsTime = StartObsTime + \
                            datetime.timedelta(
                                seconds=ObsExp + interObsSlew + AuxTimeNextTry)
                    else:
                        StartObsTime = StartObsTime + \
                            datetime.timedelta(seconds=ObsExp + interObsSlew)

                    # PreDefWindow.append(predefWind[j])
                    # ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar)
                    if ((P_GW >= obspar.minProbcut)):
                        # print("===== RESULTS ========")
                        print('++++ SCHEDULING OBS +++++')
                        P_GWarray.append(float('{:1.4f}'.format(float(P_GW))))
                        RAarray.append(
                            float('{:3.4f}'.format(float(TC.ra.deg))))
                        DECarray.append(
                            float('{:3.4f}'.format(float(TC.dec.deg))))
                        Obsarray.append(obspar.name)
                        ZenIniarray.append(
                            float('{:1.4f}'.format(float(ZenIni))))
                        ZenEndarray.append(
                            float('{:1.4f}'.format(float(ZenEnd))))
                        ObsBoolarray.append('True')
                        Exposure.append(ObsExp)
                        Delay.append(DelayObs)
                        counter = counter + 1
                    # Although OBS AND NIGHT WAS TRUE
                    elif ((P_GW <= obspar.minProbcut) and (P_GW > 0)):
                        # print("===== RESULTS ========")
                        print('++++ Probability too low +++++')
                        P_GWarray.append(
                            0)  # Careful with including the PGW in this case, afterwards is used to compute total PGW and yields to bad results
                        RAarray.append(0)
                        DECarray.append(0)
                        Obsarray.append(0)
                        ZenIniarray.append(0)
                        ZenEndarray.append(0)
                        Exposure.append(0)
                        Delay.append(0)
                        ObsBoolarray.append('ProbTooLow')
                    else:
                        # print("===== RESULTS ========")
                        print('++++ Probability is Zero +++++')
                        P_GWarray.append(0)
                        RAarray.append(0)
                        DECarray.append(0)
                        Obsarray.append(0)
                        ZenIniarray.append(0)
                        ZenEndarray.append(0)
                        Exposure.append(0)
                        Delay.append(0)
                        ObsBoolarray.append('ProbZero')
            else:
                # The event hasnt been found to be on the temporal FoV of the instrument
                # print("===== RESULTS ========")
                print('++++ The event is not in temporal FoV of the instrument +++++')
                ObservationTimearray.append(str(StartObsTime).split('.')[0].split('+')[0])
                P_GWarray.append(0)
                RAarray.append(0)
                DECarray.append(0)
                Obsarray.append(0)
                ZenIniarray.append(0)
                ZenEndarray.append(0)
                Exposure.append(0)
                Delay.append(0)
                ObsBoolarray.append('NoTemporalFoV')

        # Auxiliary jump of 12 hours to the next day
        ObservationTime = TendNight + datetime.timedelta(seconds=43200)
        # nights = nights + 1
    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    print("TOTAL POSSIBLE: ", counterTotalPossible, "DONE: ",
          counter, "Probability: ", np.sum(P_GWarray))
    print(len(ObservationTimearray), len(RAarray), len(Obsarray), len(P_GWarray), len(ObsBoolarray), len(ZenIniarray), len(ZenEndarray), len(Exposure),
         len(Delay))
    SuggestedPointings = Table(
        [ObservationTimearray, RAarray, DECarray, Obsarray, P_GWarray, ObsBoolarray, ZenIniarray, ZenEndarray, Exposure,
         Delay],
        names=['Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'Observatory', 'PGW', 'ObsInfo', 'ZenIni[deg]', 'ZenEnd[deg]',
               'Duration[s]', 'Delay[s]'])
    return (SuggestedPointings, ObservationTime0, obspar, counter)


def PGalonFoV_WindowsFromList(filename, galFile, InputObservationList, UseObs, distance, Edistance_max, Edistance_min,
                              ObservationTime0, parameters, dirName):
    # Main Parameters

    #########################
    cfg = parameters
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWGalaxyProbabilityIntegrated-Parameters'

    try:
        maxZenith = int(parser.get(section, 'maxZenith'))
        maxNights = int(parser.get(section, 'maxNights'))
        FOV = float(parser.get(section, 'FOV'))
        maxRuns = int(parser.get(section, 'maxRuns'))
        probCut = float(parser.get(section, 'probCut'))
        MinimumProbCutForCatalogue = float(
            parser.get(section, 'MinimumProbCutForCatalogue'))
        doPlot = (parser.getboolean(section, 'doPlot'))
        Duration = int(parser.get(section, 'Duration'))
        MinDuration = int(parser.get(section, 'MinDuration'))
        secondRound = (parser.getboolean(section, 'secondRound'))
        zenithWeighting = float(
            parser.get(section, 'zenithWeighting'))
        useGreytime = (parser.getboolean(section, 'useGreytime'))

    except Exception as x:
        print(x)

    print('GWGalaxyProbabilityIntegrated - Parameters:', maxZenith, maxNights, FOV, maxRuns, probCut,
          MinimumProbCutForCatalogue, doPlot,
          Duration, MinDuration, secondRound, zenithWeighting, dirName, useGreytime)

    #########################

    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    cat = LoadGalaxiesSimulation(galFile)
    print('done loading galaxies')

    # Observatory
    if UseObs == 'South':
        print('Observed form the', UseObs)
        observatory = CTASouthObservatory()
    else:
        print('Observed from the', UseObs)
        observatory = CTANorthObservatory()
    # name = filename.split('.')[0].split('/')[-1]

    # name = "Default"
    # if ('G' in filename):
    #    names = filename.split("_")
    #    name = names[0]

    # print()
    # print('-------------------   NEW LVC EVENT   --------------------')
    # print()

    print('Loading map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    has3D = True
    if (len(distnorm) == 0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D = False
    else:
        print("Found a 3D reconstruction")

    name = filename.split('.')[0].split('/')[-1]
    # correlate GW map with galaxy catalog, retrieve ordered list
    tGals, sum_dP_dV = GiveProbToGalaxy(
        prob, cat, distance, Edistance_max, Edistance_min, MinimumProbCutForCatalogue)
    tGals_aux = tGals
    tGals_aux2 = tGals

    P_GALarray = []
    P_GWarray = []
    ObservationTimearray = []
    RAarray = []
    DECarray = []
    alreadysumipixarray1 = []
    alreadysumipixarray2 = []
    PreDefWindow = []
    Round = []
    print('+++++++++++++++++++++++++++++++++++++++++++++++')
    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    print('+++++++++++++++++++++++++++++++++++++++++++++++')
    print('maxRuns: ', maxRuns, 'MinimumProbCutForCatalogue: ',
          MinimumProbCutForCatalogue)
    # SlewingTime=datetime.timedelta(seconds=210)
    ObservationTime = ObservationTime0
    time = NextWindowTools.NextObservationWindow(ObservationTime, observatory)
    WindowDurations = InputObservationList['Interval']
    # WindowDurations = [15, 17, 20, 23, 27, 33, 40, 50, 64, 85,119,178,296,595,1905]
    NightDarkRuns = NextWindowTools.CheckWindowCreateArray(
        time, observatory, WindowDurations)

    # print('EffectiveRunsTime',len(NightDarkRuns),'being',NightDarkRuns)

    predefWind = InputObservationList['pointingNumber']
    totalProb = 0.
    counter = 0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < maxRuns):
            ObservationTime = NightDarkRuns[j]
            visible, altaz, tGals_aux = VisibleAtTime(
                ObservationTime, tGals_aux, maxZenith, observatory.location)

            if (visible):

                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (maxZenith + FOV)
                visiGals = tGals_aux[visiMask]
                visiGals = ModifyCatalogue(
                    prob, visiGals, FOV, sum_dP_dV, nside)

                mask, minz = FulfillsRequirement(
                    visiGals, maxZenith, FOV, zenithWeighting, UsePix=False)

                finalGals = visiGals[mask]

                if (finalGals['dp_dV_FOV'][:1] > probCut):
                    # final galaxies within the FoV
                    if ((finalGals['dp_dV_FOV'][:1] < (2 * probCut)) and (
                            sum(P_GWarray) > 0.40) and secondRound):  # This notes LIGOVirgo type of signal
                        print('probability', finalGals['dp_dV_FOV'][:1])
                        visible, altaz, tGals_aux2 = VisibleAtTime(ObservationTime, tGals_aux2, maxZenith,
                                                                   observatory.location)
                        if (visible):
                            visiMask = altaz.alt.value > 90 - \
                                (maxZenith + FOV)
                            visiGals2 = tGals_aux2[visiMask]
                            visiGals2 = ModifyCatalogue(
                                prob, visiGals2, FOV, sum_dP_dV, nside)

                            mask, minz = FulfillsRequirement(visiGals2, maxZenith, FOV, zenithWeighting,
                                                             UsePix=False)

                            finalGals2 = visiGals2[mask]
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(prob,
                                                                                                        ObservationTime,
                                                                                                        observatory.location,
                                                                                                        finalGals2,
                                                                                                        False,
                                                                                                        visiGals2,
                                                                                                        tGals_aux2,
                                                                                                        sum_dP_dV,
                                                                                                        alreadysumipixarray2,
                                                                                                        nside, minz,
                                                                                                        maxZenith, FOV,
                                                                                                        counter, name,
                                                                                                        dirName, doPlot)

                            RAarray.append(finalGals2['RAJ2000'][:1])
                            DECarray.append(finalGals2['DEJ2000'][:1])
                            PreDefWindow.append(predefWind[j])
                            Round.append(2)
                    else:
                        # print("\n=================================")
                        # print("TARGET COORDINATES AND DETAILS...")
                        # print("=================================")
                        # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                        p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob,
                                                                                                   ObservationTime,
                                                                                                   observatory.location,
                                                                                                   finalGals, False,
                                                                                                   visiGals, tGals_aux,
                                                                                                   sum_dP_dV,
                                                                                                   alreadysumipixarray1,
                                                                                                   nside, minz,
                                                                                                   maxZenith, FOV,
                                                                                                   counter, name,
                                                                                                   dirName, doPlot)
                        RAarray.append(finalGals['RAJ2000'][:1])
                        DECarray.append(finalGals['DEJ2000'][:1])
                        PreDefWindow.append(predefWind[j])
                        Round.append(1)
                    P_GALarray.append(p_gal)
                    P_GWarray.append(p_gw)
                    ObservationTimearray.append(ObservationTime)
                    counter = counter + 1
                    # ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))

                else:
                    # print("Optimal pointing position is: ")
                    # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    print("NOT passing the cut on dp_dV_FOV > ", probCut, 'as', finalGals['dp_dV_FOV'][:1],
                          visiGals['dp_dV_FOV'][:1])
        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, PreDefWindow, Round],
                               names=['Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'PGW', 'Pgal', 'preDefWind[s]',
                                      'Round'])
    return SuggestedPointings, cat, FOV, nside




def ObservationStartperObs(obsparameters, ObservationTime0):
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


def PGWinFoV_NObs(filename, ObservationTime0, PointingFile, dirName, obsparameters):
    obs_time, SameNight, NewActiveObs, NewActiveObsStart = ObservationStartperObs(obsparameters, ObservationTime0)

    # START
#################################################################################################################################################
    name = filename.split('.')[0].split('/')[-1]
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
#################################################################################################################################################
    obspar = obsparameters[0]
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    prob = hp.pixelfunc.ud_grade(tprob, obspar.reducedNside, power=-2)
    nside = obspar.reducedNside
    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(
        prob, obspar.percentageMOC, obspar.reducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='fk5', unit=(u.deg, u.deg))
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
            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
                ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, obspar.minProbcut,
                                            obspar.maxZenith, obspar.location, obspar.minMoonSourceSeparation, obspar.useGreytime)
                if ObsBool:
                    # Round 1
                    P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                         obspar.location, obspar.maxZenith, obspar.FOV, name, pixlist, ipixlistHR, counter, dirName, obspar.useGreytime, obspar.doPlot)
                    # print(P_GW, obspar.name)
                    if ((P_GW <= obspar.minProbcut) and obspar.secondRound):
                        # Try Round 2
                        # print('The minimum probability cut being', minProbcut * 100, '% is, unfortunately, not reached.')
                        yprob1 = highres
                        P_GW, TC, pixlist1, ipixlistHR1 = ComputeProbability2D(prob, yprob1, radecs, obspar.reducedNside, obspar.HRnside, obspar.minProbcut, ObservationTime,
                                                                               obspar.location, obspar.maxZenith, obspar.FOV, name, pixlist1, ipixlistHR1, counter, dirName, obspar.useGreytime, obspar.doPlot)
                        if ((P_GW <= obspar.minProbcut)):
                            print('Fail')
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

    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, Round, ObsName], names=[
                               'Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'PGW', 'Round', 'ObsName'])

    return SuggestedPointings, obsparameters


def PGalinFoV_NObs(filename, ObservationTime0, PointingFile, dirName, obsparameters):

    obs_time, SameNight, NewActiveObs, NewActiveObsStart = ObservationStartperObs(obsparameters, ObservationTime0)
    # START
#################################################################################################################################################
    name = filename.split('.')[0].split('/')[-1]
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
#################################################################################################################################################
    obspar = obsparameters[0]
    print('Loading map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    # load galaxy catalogue
    galFile = obspar.datasetDir + obspar.galcatName
    
    if not obspar.mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)
    print('done loading galaxies')

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(
            prob, distmu, distsigma, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(
            prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.minimumProbCutForCatalogue)

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
    while (i < 500) & any(SameNight):
        for j in range(len(NewActiveObs)):
            obspar = NewActiveObs[j]
            ObservationTime = NewActiveObsTime[j]
            if ITERATION_OBS == len(obsparameters):
                TIME_MIN_ALL = []
                ITERATION_OBS = 0

            ITERATION_OBS += 1

            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
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
                        visiGals, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=False)
                    if obspar.useGreytime:
                        maskgrey = FulfillsRequirementGreyObservations(
                            ObservationTime, visiGals, obspar.location, obspar.minMoonSourceSeparation)
                        finalGals = visiGals[mask & maskgrey]
                    if not obspar.useGreytime:
                        finalGals = visiGals[mask]

                    if (finalGals['dp_dV_FOV'][:1] > obspar.minProbcut):
                        # final galaxies within the FoV
                        # This notes LIGOVirgo type of signal
                        if ((finalGals['dp_dV_FOV'][:1] < (2 * obspar.minProbcut)) and (sum(P_GWarray) > 0.40) and obspar.secondRound):
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
                                    visiGals2, obspar.maxZenith, obspar.FOV, obspar.zenithWeighting, UsePix=False)

                                if obspar.useGreytime:
                                    maskgrey = FulfillsRequirementGreyObservations(
                                        ObservationTime, visiGals2, obspar.location, obspar.minMoonSourceSeparation)
                                    finalGals2 = visiGals2[mask & maskgrey]
                                if not obspar.useGreytime:
                                    finalGals2 = visiGals2[mask]

                                p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(
                                    prob, ObservationTime, obspar.location, finalGals2, False, visiGals2, tGals_aux2, sum_dP_dV, alreadysumipixarray2, nside, minz, obspar.maxZenith, obspar.FOV, counter, name, dirName, obspar.doPlot)

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

                        else:
                            # print("\n=================================")
                            # print("TARGET COORDINATES AND DETAILS...")
                            # print("=================================")
                            # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                            p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(
                                prob, ObservationTime, obspar.location, finalGals, False, visiGals, tGals_aux, sum_dP_dV, alreadysumipixarray1, nside, minz, obspar.maxZenith, obspar.FOV, counter, name, dirName, obspar.doPlot)
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
                        # ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))

                    else:
                        # print("Optimal pointing position is: ")
                        # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                        print(finalGals['dp_dV_FOV'][:1])
                        print("NOT passing the cut on dp_dV_FOV > ",
                              obspar.minProbcut)

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

    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, Round, ObsName], names=[
                               'Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'PGW', 'Pgal', 'Round', 'ObsName'])
    return SuggestedPointings, cat, obsparameters