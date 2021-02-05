import sys
sys.path.append('./include')
from GWCTAPointingTools import *
from ObservingTimes import *
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

iers_file='./finals2000A.all'
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

    print('Loading GW map from ', filename)
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

def PGWonFoV_WindowOptimisation(filename,InputChar,TC,parameters,dirName):
    UseObs =InputChar['Observatory']
    run=InputChar['run']
    mergerID=InputChar['MergerID']
    zenith=InputChar['Zenith']
    ObservationTime0= datetime.datetime.strptime(InputChar['Time'], '%Y-%m-%d %H:%M:%S.%f')
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
    path = dirName + '/PointingPlotting/' + run + '_' + mergerID + '/EvolutionPlot/'
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
    Exposure = []
    ZenIniarray = []
    ZenEndarray = []
    ObsBoolarray = []


    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading GW map from ', filename)
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
    followupDelay = 30
    SlewingTime = 210
    interObsSlew = 10
    total_followupDelay = followupDelay + SlewingTime

    # From the injection time, look for the next window. Time is the time of the first observation
    ObservationTime = ObservationTime0 + datetime.timedelta(seconds=total_followupDelay)

    print("Main info of this scheduling:",totalTime, total_followupDelay, run, mergerID, TC,observatory,zenith)

    TotalNights=1
    counter = 0
    counterTotalPossible = 0
    MaxRuns = 100

    for nights in range(0,TotalNights):
        nights=nights+1
        TstartNight = NextWindowTools.NextObservationWindow(time = ObservationTime, obsSite = observatory)
        if(TstartNight != False):
            TendNight = NextWindowTools.EndObservationWindow(TstartNight,observatory)
        print("Night number",nights," window starts at", TstartNight,"finishes at", TendNight)

        # Look for the first time that the C.R. observable is > 5%
        TemporalBin = datetime.timedelta(seconds=1800) # Every 30 minutes
        MinProbCut = 0.00001

        for i in range(0, 12):
            checkTime=TstartNight+i*TemporalBin
            print((checkTime-TendNight).seconds)
            if((TendNight-checkTime).seconds<0):
                FistObsTime = False
                print("The source is not on the temporal FoV of the instrument")
                break
            ObsBool, yprob = ZenithAngleCut(prob, nside, checkTime, MinProbCut, max_zenith, observatory.Location,usegreytime=False)
            if(ObsBool==True):
                StartObsTime=checkTime
                break
        print("The first observing time is",StartObsTime, "MinProbCut",MinProbCut,"max_zenith",max_zenith)

        # Get the first 100 clusters of probability for the masked map
        for tt in range(0,1000):
            if(tt>MaxRuns):
                break
            if((TendNight-StartObsTime).seconds<0):
                print("---- NIGHT IS OVER ----")
                break
            print("Starting observation number ", counter)
            TotalExposure = (TendNight-StartObsTime).seconds
            DelayObs = (StartObsTime-ObservationTime0).seconds

            P_GW, TC, ObsExp, ZenIni, ZenEnd, pixlist, ipixlistHR = ComputeProbability2D_SelectClusters(prob, highres, radecs, ReducedNside, HRnside,
                                                         MinProbCut,TotalExposure, StartObsTime, DelayObs, observatory,
                                                         max_zenith, FOV, run, mergerID, pixlist, ipixlistHR, counter,
                                                         dirName, False, True)
            counterTotalPossible = counterTotalPossible+1
            print("PGW", P_GW,'COORDINATES:', TC,'ZENITH CHANGE', ZenIni,'->', ZenEnd)
            print("TotalExposure: ",TotalExposure,"DelayObs:",DelayObs, "Observation Number:",counter, "Exposure:", ObsExp)

            #print('P_GW', P_GW)

            ObservationTimearray.append(StartObsTime)
            Exposure.append(ObsExp)
            # PreDefWindow.append(predefWind[j])
            # ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar)
            if ((P_GW >= MinProbCut)):
                P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                ZenIniarray.append(np.float('{:1.4f}'.format(np.float(ZenIni))))
                ZenEndarray.append(np.float('{:1.4f}'.format(np.float(ZenEnd))))
                ObsBoolarray.append(True)
                counter = counter + 1
            else:
                P_GWarray.append(0)
                RAarray.append(0)
                DECarray.append(0)
                ZenIniarray.append(0)
                ZenEndarray.append(0)
                ObsBoolarray.append(False)
                print('Probability too low')

            StartObsTime = StartObsTime + datetime.timedelta(seconds=ObsExp+interObsSlew)

        ObservationTime = TendNight + datetime.timedelta(seconds=43200)  # Auxiliary jump of 12 hours to the next day

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    print("TOTAL POSSIBLE: ",counterTotalPossible,"DONE: ",counter)
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,ObsBoolarray,ZenIniarray, ZenEndarray, Exposure], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','ObsBoolarray','ZenIni','ZenEnd','Duration'])
    return (SuggestedPointings, ObservationTime0, FOV, nside, len(ObservationTimearray))

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
    
    print('Loading GW map from ', filename)
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

    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
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
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,finalGals2,False,visiGals2,tGals_aux2,sum_dP_dV, alreadysumipixarray2,nside, minz,max_zenith, FOV,counter,name,dirName,doplot)


                            RAarray.append(finalGals2['RAJ2000'][:1])
                            DECarray.append(finalGals2['DEJ2000'][:1])
                            PreDefWindow.append(predefWind[j])
                            Round.append(2)
                    else:
                        #print("\n=================================")
                        #print("TARGET COORDINATES AND DETAILS...")
                        #print("=================================")
                        #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                        p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,finalGals,False, visiGals,tGals_aux, sum_dP_dV,alreadysumipixarray1,nside, minz,max_zenith, FOV, counter,name,dirName,doplot)
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

def PGalonFoV_PixRegion(filename,ObservationTime0):

    # Main Parameters

    ###############################
    max_zenith = 60
    MaxNights = 1
    FOV = 2.5
    MaxRuns = 20  # Maximum number of pointings/runs
    probCut = 0.005
    MinimumProbCutForCatalogue = 0.01
    doplot = False
    FulFillReq_Percentage = 0.75
    NewNside = 64
    PercentCoverage = 0.99
    Duration=5
    MinDuration=1
    ###############################

    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    # galFile='./GLADE_2clean.txt'

    cat = LoadGalaxies(galFile)
    print('done loading galaxies')

    name = filename.split('.')[0].split('/')[-1]

    # name = "Default"
    # if ('G' in filename):
    #    names = filename.split("_")
    #    name = names[0]

    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()

    print('Loading GW map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
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
    tGals, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D, MinimumProbCutForCatalogue)
    tGals_aux = tGals
    tGals_aux2 = tGals

    P_GALarray = []
    P_GWarray = []
    ObservationTimearray = []
    RAarray = []
    DECarray = []
    alreadysumipixarray1 = []
    alreadysumipixarray2 = []


    # In case one wants to see if a next round would give us better results..
    # So the time will be j-1
    nextround = False
    Round = []
    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    print('MaxRuns: ', MaxRuns, 'MinimumProbCutForCatalogue: ', MinimumProbCutForCatalogue)

    SlewingTime = datetime.timedelta(minutes=210)
    ObservationTime = ObservationTime0 + SlewingTime
    time = NextWindowTools.NextObservationWindow(ObservationTime, CTANorthObservatory())
    WindowDurations = [15, 17, 20, 23, 27, 33, 40, 50, 64, 85,119,178,296,595,1905]
    NightDarkRuns = NextWindowTools.CheckWindowCreateArray(time, CTANorthObservatory(), WindowDurations)

    totalProb = 0.
    n = 0

    ###############################

    # Get the RA & DEC of pixles of the pixels in an enclosed probability region (% precised by PercentCoverage).
    # Reduce these RA DEC to angles in maps with smaller resolution (NewNside)


    pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob, PercentCoverage, NewNside)



    ##############################

    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            if (nextround):
                ObservationTime = NightDarkRuns[j - 1]
                nextround = False
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, max_zenith,
                                                      CTANorthObservatory().Location)

            if (visible):

                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                visiGals = tGals_aux[visiMask]

                mask, minz = FulfillsRequirement(visiGals, max_zenith, FOV, FulFillReq_Percentage, UsePix=True)

                finalGals = visiGals[mask]
                visiPix = ModifyCataloguePIX(pix_ra1, pix_dec1, ObservationTime, max_zenith, prob, finalGals, FOV,
                                             sum_dP_dV, nside, NewNside, minz)

                if (visiPix['PIXFOVPROB'][:1] > probCut):
                    n = n + 1
                    # final galaxies within the FoV

                    # print("\n=================================")
                    # print("TARGET COORDINATES AND DETAILS...")
                    # print("=================================")
                    # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob, ObservationTime,
                                                                                               visiPix, True, visiGals,
                                                                                               tGals_aux, sum_dP_dV,
                                                                                               alreadysumipixarray1,
                                                                                               nside, minz, max_zenith,
                                                                                               FOV, name, doplot)
                    RAarray.append(visiPix['PIXRA'][:1])
                    DECarray.append(visiPix['PIXDEC'][:1])
                    Round.append(1)
                    P_GALarray.append(p_gal)
                    P_GWarray.append(p_gw)
                    ObservationTimearray.append(ObservationTime)

                else:
                    print("Optimal pointing position is: ")
                    print(visiPix['PIXRA', 'PIXDEC', 'PIXFOVPROB'][:1])
                    print("NOT passing the cut on dp_dV_FOV > ", probCut)

        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, P_GALarray, Round],
                               names=['Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'PGW',
                                      'Pgal', 'Round'], )
    print(SuggestedPointings)
    print("Name", name, "Total GW probability covered: ", sum(P_GWarray), "Total Gal probability covered: ", sum(P_GALarray),
    "Number of runs that fulfill darkness condition  :", len(NightDarkRuns), "Number of effective pointings: ",
    len(ObservationTimearray))
    return SuggestedPointings,cat,area

    #print("===========================================================================================")
    #print()
    #print()
