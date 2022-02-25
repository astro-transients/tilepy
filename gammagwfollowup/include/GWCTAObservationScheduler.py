import numpy as np
from astropy import units as u
from astropy.utils import iers
import astropy.coordinates as co
import random
from .GWHESSPointingTools import (LoadHealpixMap, Get90RegionPixReduced,
                                  ZenithAngleCut,ComputeProbability2D,
                                  FulfillsRequirement,ModifyCataloguePIX, ComputeProbPGALIntegrateFoV,
                                  CTASouthObservatory,CTANorthObservatory,
                                  VisibleAtTime,LoadGalaxies,CorrelateGalaxies_LVC,ModifyCatalogue
                                  )
from .GWCTAPointingTools import (NextWindowTools,
                                 ComputeProbability2D_SelectClusters,GiveProbToGalaxy, LoadGalaxiesSimulation
                                 )
from .ObservingTimes import ObtainObservingTimes
import os
from astropy.table import Table
import datetime
import healpy as hp

from six.moves import configparser
import six
if six.PY2:
  ConfigParser = configparser.SafeConfigParser
else:
  ConfigParser = configparser.ConfigParser

############################################

#              General definitions         #

############################################
#iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))
#iers_file='./finals2000A.all'
# iers_file = './dataset/finals2000A.all'
iers_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)

def PGWonFoV(filename,InputChar,TC,parameters,dirName):
    UseObs =InputChar['Observatory']
    run=InputChar['run']
    mergerID=InputChar['MergerID']
    zenith=InputChar['Zenith']
    ObservationTime0= datetime.datetime.strptime(InputChar['Time'], '%Y-%m-%d %H:%M:%S.%f')
    # Main parameters

    ##################
    cfg = parameters
    parser = ConfigParser()
    #print("parameters",parameters)
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWProbDensityIntegration-Parameters'

    try:
        max_zenith = int(parser.get(section, 'max_zenith'))
        MaxNights = int(parser.get(section, 'MaxNights'))
        FOV = float(parser.get(section, 'FOV'))
        MaxRuns = int(parser.get(section, 'MaxRuns'))
        MinProbCut = float(parser.get(section, 'MinProbCut'))
        doplot = (parser.getboolean(section, 'doplot'))
        Duration = int(parser.get(section, 'Duration'))
        MinDuration = int(parser.get(section, 'MinDuration'))
        MinSlewing = int(parser.get(section, 'MinSlewing'))
        SecondRound = (parser.getboolean(section, 'SecondRound'))
        PercentCoverage = float(parser.get(section, 'PercentCoverage'))
        ReducedNside = int(parser.get(section, 'ReducedNside'))
        HRnside = int(parser.get(section, 'HRnside'))
        UseGreytime = (parser.getboolean(section, 'UseGreytime'))

    except Exception as x:
        print(x)

    print('GWProbDensityIntegration-Parameters:', max_zenith,MaxNights,FOV, MaxRuns, MinProbCut, doplot, Duration, MinDuration, SecondRound,PercentCoverage,ReducedNside,HRnside,UseGreytime)

    ##################

    #Observatory
    if UseObs=='South':
        print('Observed form the',UseObs)
        observatory=CTASouthObservatory()
    else:
        print('Observed from the',UseObs)
        observatory =CTANorthObservatory()

    # link to the GW map
    name = filename.split('.')[0].split('/')[-1]
    #if('G' in filename):
    #    names = filename.split("_")
    #    name= names[0]

    random.seed()

    RAarray = []
    DECarray = []
    pixlist = []
    ipixlistHR=[]
    pixlist1 = []
    ipixlistHR1=[]
    P_GWarray = []
    ObservationTimearray= []
    Round = []
    PreDefWindow =[]
    Duration = []


    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading map from ', filename)
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob,ReducedNside,power=-2)

    nside = ReducedNside

    highres=hp.pixelfunc.ud_grade(prob, HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix,areapix=Get90RegionPixReduced(prob,PercentCoverage,ReducedNside)
    radecs= co.SkyCoord(rapix,decpix, frame='fk5', unit=(u.deg, u.deg))
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


    #ToDO need to change that to no max integration time (limited by visibility)
    #ToDo need to change the maximum latency to 48 hours
    totalTime = 7200

    followupDelay=30
    SlewingTime=210
    total_followupDelay=followupDelay+SlewingTime


    # Else, continue
    #total_followupDelay=0
    print(totalTime, total_followupDelay, run, mergerID, TC,observatory,zenith)
    # Obtain the times from the comparison to a given spectrum and sensitivities computed by cssens
    # Minimum set to 10 seconds!
    # Starting of observation and duration of observation
    tstar, tobs = ObtainObservingTimes(totalTime, total_followupDelay, run, mergerID,observatory,zenith)

    counter = 0
    # Loop where probabilities are computed comes here

    for j in range(0,len(tobs),1):
        if(j>MaxRuns):
            break
        else:
            ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar[j])
            #print(ObservationTime)

            ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, MinProbCut, max_zenith, observatory.Location,
                                        False)
            #print(ObsBool)
            if ObsBool:
                # Round 1
                #print('Round ',j,'counter',counter)
                P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, ReducedNside, HRnside,
                                                                 MinProbCut, ObservationTime, observatory.Location,
                                                                 max_zenith, FOV, name, pixlist, ipixlistHR, counter,
                                                              dirName, False, doplot)
                print('P_GW',P_GW)
                if ((P_GW >= MinProbCut)):
                    P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                    RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                    DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                    ObservationTimearray.append(ObservationTime)
                    Duration.append(tobs[j])
                    #PreDefWindow.append(predefWind[j])
                    #ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar)
                    counter = counter + 1
                else:
                    print('Probability too low')

            #else:
            #    break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,Duration], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','Duration'])

    return (SuggestedPointings, ObservationTime0, FOV, nside,len(tobs))

def PGWonFoV_WindowOptimisation(filename,InputChar,TC,parameters,datasetDir,outDir):
    UseObs =InputChar['Observatory']
    run=InputChar['run']
    mergerID=InputChar['MergerID']
    #zenith=InputChar['Zenith']
    ObservationTime0= datetime.datetime.strptime(InputChar['Time'], '%Y-%m-%d %H:%M:%S')
    # Main parameters

    ##################

    cfg = parameters
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWProbDensityIntegration-Parameters'

    try:
        max_zenith = int(parser.get(section, 'max_zenith'))
        MaxNights = int(parser.get(section, 'MaxNights'))
        FOV = float(parser.get(section, 'FOV'))
        MaxRuns = int(parser.get(section, 'MaxRuns'))
        MinProbCut = float(parser.get(section, 'MinProbCut'))
        doplot = (parser.getboolean(section, 'doplot'))
        Duration = int(parser.get(section, 'Duration'))
        MinDuration = int(parser.get(section, 'MinDuration'))
        MinSlewing = int(parser.get(section, 'MinSlewing'))
        SecondRound = (parser.getboolean(section, 'SecondRound'))
        PercentCoverage = float(parser.get(section, 'PercentCoverage'))
        ReducedNside = int(parser.get(section, 'ReducedNside'))
        HRnside = int(parser.get(section, 'HRnside'))
        UseGreytime = (parser.getboolean(section, 'UseGreytime'))

    except Exception as x:
        print(x)

    print('GWProbDensityIntegration-Parameters:', max_zenith,MaxNights,FOV, MaxRuns, MinProbCut, doplot, Duration, MinDuration, SecondRound,PercentCoverage,ReducedNside,HRnside,UseGreytime)

    ##################
    path = outDir + '/PointingPlotting/' + run + '_' + mergerID + '/EvolutionPlot/'
    #Observatory
    if UseObs=='South':
        print('Observed form the',UseObs)
        observatory=CTASouthObservatory()
    else:
        print('Observed from the',UseObs)
        observatory =CTANorthObservatory()

    # link to the GW map
    name = filename.split('.')[0].split('/')[-1]
    #if('G' in filename):
    #    names = filename.split("_")
    #    name= names[0]

    random.seed()

    RAarray = []
    DECarray = []
    Obsarray = []
    pixlist = []
    ipixlistHR=[]
    pixlist1 = []
    ipixlistHR1=[]
    P_GWarray = []
    ObservationTimearray= []
    Exposure = []
    Delay = []
    ZenIniarray = []
    ZenEndarray = []
    ObsBoolarray = []


    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading map from ', filename)
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob,ReducedNside,power=-2)

    nside = ReducedNside

    highres=hp.pixelfunc.ud_grade(prob, HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix,areapix=Get90RegionPixReduced(prob,PercentCoverage,ReducedNside)
    radecs= co.SkyCoord(rapix,decpix, frame='fk5', unit=(u.deg, u.deg))
    has3D = True
    if (len(distnorm) == 0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D = False
    else:
        print("Found a 3D reconstruction")


    # print('----------   NEW FOLLOW-UP ATTEMPT   ----------')

    # Set of delays and slewing times
    totalTime = 172800    # 48h
    followupDelay = 30    # Minimum delay to start observaiton
    SlewingTime = 210     # Slewing to first position
    interObsSlew = 20     # Slewing time between observations
    total_followupDelay = followupDelay + SlewingTime

    # From the injection time, look for the next window. Time is the time of the first observation
    ObservationTime = ObservationTime0 + datetime.timedelta(seconds=total_followupDelay)

    print("Main info of this scheduling:",totalTime, total_followupDelay, run, mergerID, TC,observatory)

    TotalNights = 2
    counter = 0
    counterTotalPossible = 0
    MaxRuns = 20
    #GrbsensMax = 16384
    AuxTimeNextTry = 900   # 15 minutes.  Allows to check if 15 minutes later there is a better opportunity
    #MinProbCut = 0.005


    # case 1 = TstartNight == False
    # case 2 = "The source is not on the temporal FoV of the instrument"
    for nights in range(0,TotalNights):
        TstartNight = NextWindowTools.NextObservationWindow(time = ObservationTime, obsSite = observatory)
        print("TstartNight",TstartNight)
        if(TstartNight == False):
            TendNight = False
            print("Night number", nights, " window starts at", TstartNight, "finishes at", TendNight)
            ObsCase = 'NoDarknessFound'
            print("===== RESULTS ========")
            print('++++ No Darkness Time Found +++++')
            ObservationTimearray.append(ObservationTime)
            P_GWarray.append(0)
            RAarray.append(0)
            DECarray.append(0)
            ZenIniarray.append(0)
            ZenEndarray.append(0)
            Exposure.append(0)
            ObsBoolarray.append(ObsCase)
            break
        else:
            print("Observations can be allocated")
            TendNight = NextWindowTools.EndObservationWindow(TstartNight,observatory)
            print("Night number", nights, " window starts at", TstartNight, "finishes at", TendNight)

            # Look for the first time that the C.R. observable is > 5%
            TemporalBin = datetime.timedelta(seconds = 600) # Every 10 minutes
            for i in range(0, 24):
                checkTime = TstartNight + i*TemporalBin
                if((TendNight-checkTime).total_seconds()<0):
                    ObsBool = False
                    print("The source is not on the temporal FoV of the instrument")
                    break
                ObsBool, yprob = ZenithAngleCut(prob, nside, checkTime, MinProbCut, max_zenith, observatory.Location,usegreytime=False)
                print(checkTime, ObsBool)
                if(ObsBool==True):
                    StartObsTime=checkTime
                    print("The first time the GW is observed is", StartObsTime, "MinProbCut", MinProbCut, "max_zenith", max_zenith)
                    break
            if(ObsBool):
                # Get the first 100 clusters of probability for the masked map
                for tt in range(0,1000):
                    print("  ")
                    print('++++++++++++++++++++++++++++++++++')
                    print("---- NEW OBSERVATION ATTEMPT ----")
                    print('++++++++++++++++++++++++++++++++++')
                    print("  ")
                    print("Observation number ", counter)
                    print("Starting time",StartObsTime)
                    if(tt>MaxRuns):
                        break
                    if((TendNight-datetime.timedelta(seconds =interObsSlew)-StartObsTime).total_seconds()<=0):
                        print("---- NIGHT IS OVER ----")
                        break
                    TotalExposure = int((TendNight-StartObsTime).total_seconds())     # Total Exposure for the night
                    DelayObs = int((StartObsTime-ObservationTime0).total_seconds())

                    print("ObservationTime0", ObservationTime0)
                    print("TotalExposure", TotalExposure)
                    print("DelayObs",DelayObs)
                    print("----------------------------")
                    # TO-DO: Need to have a class to set all this parameters for an observation, so there are less arguments passed.
                    P_GW, TC, ObsExp, ZenIni, ZenEnd, ObsCase, pixlist, ipixlistHR = ComputeProbability2D_SelectClusters(prob, highres, radecs, ReducedNside, HRnside,
                                                         MinProbCut,TotalExposure, StartObsTime, DelayObs, interObsSlew, observatory,
                                                         max_zenith, FOV, run, mergerID, pixlist, ipixlistHR, counter,datasetDir,outDir, False, False)
                    print("=============")
                    print("P_GW, ObsExp, ZenIni, ZenEnd, ObsCase")
                    print(P_GW, ObsExp, ZenIni, ZenEnd,ObsCase)
                    print("=============")
                    counterTotalPossible = counterTotalPossible+1
                    #print("PGW", P_GW,'COORDINATES:', TC,'ZENITH CHANGE', ZenIni,'->', ZenEnd)
                    #print("TotalExposure: ",TotalExposure,"DelayObs:",DelayObs, "Observation Number:",counter, "Exposure:", ObsExp)

                    ObservationTimearray.append(StartObsTime)

                    if (ObsCase == 'TimeNotEnoughIte' or ObsCase == 'TimeNotEnough'):
                        StartObsTime = StartObsTime + datetime.timedelta(seconds = ObsExp + interObsSlew + AuxTimeNextTry)
                    else:
                        StartObsTime = StartObsTime + datetime.timedelta(seconds=ObsExp + interObsSlew)


                    # PreDefWindow.append(predefWind[j])
                    # ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar)
                    if ((P_GW >= MinProbCut)):
                        print("===== RESULTS ========")
                        print('++++ SCHEDULING OBS +++++')
                        P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                        RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                        DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                        Obsarray.append(observatory.Name)
                        ZenIniarray.append(np.float('{:1.4f}'.format(np.float(ZenIni))))
                        ZenEndarray.append(np.float('{:1.4f}'.format(np.float(ZenEnd))))
                        ObsBoolarray.append('True')
                        Exposure.append(ObsExp)
                        Delay.append(DelayObs)
                        counter = counter + 1
                    elif((P_GW <= MinProbCut) and (P_GW>0)): # Although OBS AND NIGHT WAS TRUE
                        print("===== RESULTS ========")
                        print('++++ Probability too low +++++')
                        P_GWarray.append(0)  # Careful with including the PGW in this case, afterwards is used to compute total PGW and yields to bad results
                        RAarray.append(0)
                        DECarray.append(0)
                        Obsarray.append(0)
                        ZenIniarray.append(0)
                        ZenEndarray.append(0)
                        Exposure.append(0)
                        Delay.append(0)
                        ObsBoolarray.append('ProbTooLow')
                    else:
                        print("===== RESULTS ========")
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
                print("===== RESULTS ========")
                print('++++ The event is not in temporal FoV of the instrument +++++')
                ObservationTimearray.append(ObservationTime)
                P_GWarray.append(0)
                RAarray.append(0)
                DECarray.append(0)
                Obsarray.append(0)
                ZenIniarray.append(0)
                ZenEndarray.append(0)
                Exposure.append(0)
                Delay.append(0)
                ObsBoolarray.append('NoTemporalFoV')

        ObservationTime = TendNight + datetime.timedelta(seconds=43200)  # Auxiliary jump of 12 hours to the next day
        #nights = nights + 1
    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    print("TOTAL POSSIBLE: ",counterTotalPossible,"DONE: ",counter)
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray, Obsarray, P_GWarray,ObsBoolarray,ZenIniarray, ZenEndarray, Exposure, Delay], names=['Observation Time UTC','RA(deg)','DEC(deg)','Observatory','PGW','ObsInfo','ZenIni','ZenEnd','Duration', 'Delay'])
    return (SuggestedPointings, ObservationTime0, FOV, nside, counter)

def PGalonFoV(filename,galFile,InputObservationList,UseObs,distance,Edistance_max,Edistance_min,ObservationTime0,parameters,dirName):


    # Main Parameters

    #########################
    cfg = parameters
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWGalaxyProbabilityIntegrated-Parameters'

    try:
        max_zenith = int(parser.get(section, 'max_zenith'))
        MaxNights = int(parser.get(section, 'MaxNights'))
        FOV = float(parser.get(section, 'FOV'))
        MaxRuns = int(parser.get(section, 'MaxRuns'))
        probCut = float(parser.get(section, 'probCut'))
        MinimumProbCutForCatalogue = float(parser.get(section, 'MinimumProbCutForCatalogue'))
        doplot = (parser.getboolean(section, 'doplot'))
        Duration = int(parser.get(section, 'Duration'))
        MinDuration = int(parser.get(section, 'MinDuration'))
        SecondRound = (parser.getboolean(section, 'SecondRound'))
        FulFillReq_Percentage = float(parser.get(section, 'FulFillReq_Percentage'))
        UseGreytime = (parser.getboolean(section, 'UseGreytime'))

    except Exception as x:
        print(x)

    print('GWGalaxyProbabilityIntegrated - Parameters:', max_zenith, MaxNights, FOV, MaxRuns, probCut, MinimumProbCutForCatalogue, doplot,
          Duration, MinDuration, SecondRound, FulFillReq_Percentage, dirName, UseGreytime)

    #########################
    
    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    cat = LoadGalaxiesSimulation(galFile)
    print('done loading galaxies')

    #Observatory
    if UseObs=='South':
        print('Observed form the',UseObs)
        observatory=CTASouthObservatory()
    else:
        print('Observed from the',UseObs)
        observatory =CTANorthObservatory()
    #name = filename.split('.')[0].split('/')[-1]
    
    #name = "Default"
    #if ('G' in filename):
    #    names = filename.split("_")
    #    name = names[0]
    
    
    #print()
    #print('-------------------   NEW LVC EVENT   --------------------')
    #print()
    
    print('Loading map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)
    
    has3D=True
    if(len(distnorm)==0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D=False
    else:
        print("Found a 3D reconstruction")

    name = filename.split('.')[0].split('/')[-1]
    # correlate GW map with galaxy catalog, retrieve ordered list
    tGals,sum_dP_dV=GiveProbToGalaxy(prob,cat,distance,Edistance_max,Edistance_min,MinimumProbCutForCatalogue)
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
    print('MaxRuns: ',MaxRuns,'MinimumProbCutForCatalogue: ',MinimumProbCutForCatalogue)
    #SlewingTime=datetime.timedelta(seconds=210)
    ObservationTime=ObservationTime0
    time = NextWindowTools.NextObservationWindow(ObservationTime, observatory)
    WindowDurations= InputObservationList['Interval']
    #WindowDurations = [15, 17, 20, 23, 27, 33, 40, 50, 64, 85,119,178,296,595,1905]
    NightDarkRuns = NextWindowTools.CheckWindowCreateArray(time, observatory, WindowDurations)

    # print('EffectiveRunsTime',len(NightDarkRuns),'being',NightDarkRuns)

    predefWind=InputObservationList['pointingNumber']
    totalProb=0.
    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, max_zenith,observatory.Location)
            
            if (visible):
                
                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (max_zenith+FOV)
                visiGals= tGals_aux[visiMask]
                visiGals = ModifyCatalogue(prob,visiGals,FOV,sum_dP_dV,nside)
                
                mask, minz = FulfillsRequirement(visiGals, max_zenith,FOV,FulFillReq_Percentage,UsePix=False)

                finalGals = visiGals[mask]

                if(finalGals['dp_dV_FOV'][:1] > probCut):
                    # final galaxies within the FoV
                    if ((finalGals['dp_dV_FOV'][:1] < (2 * probCut)) and (sum(P_GWarray) > 0.40) and SecondRound):  # This notes LIGOVirgo type of signal
                        print('probability', finalGals['dp_dV_FOV'][:1])
                        visible, altaz, tGals_aux2 = VisibleAtTime(ObservationTime, tGals_aux2, max_zenith,observatory.Location)
                        if (visible):
                            visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                            visiGals2 = tGals_aux2[visiMask]
                            visiGals2 = ModifyCatalogue(prob,visiGals2, FOV, sum_dP_dV,nside)
                            
                            mask, minz = FulfillsRequirement(visiGals2, max_zenith,FOV,FulFillReq_Percentage,UsePix=False)

                            finalGals2 = visiGals2[mask]
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,observatory.Location,finalGals2,False,visiGals2,tGals_aux2,sum_dP_dV, alreadysumipixarray2,nside, minz,max_zenith, FOV,counter,name,dirName,doplot)


                            RAarray.append(finalGals2['RAJ2000'][:1])
                            DECarray.append(finalGals2['DEJ2000'][:1])
                            PreDefWindow.append(predefWind[j])
                            Round.append(2)
                    else:
                        #print("\n=================================")
                        #print("TARGET COORDINATES AND DETAILS...")
                        #print("=================================")
                        #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                        p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,observatory.Location,finalGals,False, visiGals,tGals_aux, sum_dP_dV,alreadysumipixarray1,nside, minz,max_zenith, FOV, counter,name,dirName,doplot)
                        RAarray.append(finalGals['RAJ2000'][:1])
                        DECarray.append(finalGals['DEJ2000'][:1])
                        PreDefWindow.append(predefWind[j])
                        Round.append(1)
                    P_GALarray.append(p_gal)
                    P_GWarray.append(p_gw)
                    ObservationTimearray.append(ObservationTime)
                    counter = counter + 1
                    #ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))
            
                else:
                    #print("Optimal pointing position is: ")
                    #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    print("NOT passing the cut on dp_dV_FOV > ",probCut,'as',finalGals['dp_dV_FOV'][:1],visiGals['dp_dV_FOV'][:1] )
        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray ,P_GWarray,P_GALarray,PreDefWindow,Round], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','Pgal','preDefWind','Round'])
    return SuggestedPointings,cat,FOV,nside
