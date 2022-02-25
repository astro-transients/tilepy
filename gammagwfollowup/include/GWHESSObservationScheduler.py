from .GWHESSPointingTools import (NightDarkObservation,
                                  NightDarkObservationwithGreyTime,
                                  LoadHealpixMap,LoadGalaxies,CorrelateGalaxies_LVC,
                                  Get90RegionPixReduced,
                                  SubstractPointings2D,ZenithAngleCut,ComputeProbability2D,
                                  VisibleAtTime, FulfillsRequirement, SimpleGWprob,ComputeProbBCFOV,
                                  Tools,LoadGalaxies_SteMgal, CorrelateGalaxies_LVC_SteMass, SubstractPointings,
                                  ModifyCatalogue,FulfillsRequirementGreyObservations,ComputeProbPGALIntegrateFoV,
                                  ModifyCataloguePIX)
import random
import numpy as np
from astropy.table import Table
import astropy.coordinates as co
from astropy import units as u
import healpy as hp
from six.moves import configparser
import six
if six.PY2:
  ConfigParser = configparser.SafeConfigParser
else:
  ConfigParser = configparser.ConfigParser
############################################

#              General definitions              #

############################################

def PGWinFoV(filename,ObservationTime0,PointingFile,galFile,parameters,dirName,Observatory):

    # Main parameters

    ##################
    cfg = parameters
    parser = ConfigParser()
    parser.read(cfg)
    parser.sections()
    section = 'GWBestGalaxyParameters'

    try:
        max_zenith = int(parser.get(section, 'max_zenith'))
        MaxNights = int(parser.get(section, 'MaxNights'))
        FOV = float(parser.get(section, 'FOV'))
        MaxRuns = int(parser.get(section, 'MaxRuns'))
        MinProbCut = float(parser.get(section, 'MinProbCut'))
        doplot = (parser.getboolean(section, 'doplot'))
        Duration = int(parser.get(section, 'Duration'))
        MinDuration = int(parser.get(section, 'MinDuration'))
        SecondRound = (parser.getboolean(section, 'SecondRound'))
        PercentCoverage = float(parser.get(section, 'PercentCoverage'))
        ReducedNside = int(parser.get(section, 'ReducedNside'))
        HRnside = int(parser.get(section, 'HRnside'))
        UseGreytime = (parser.getboolean(section, 'UseGreytime'))
        Mangrove = (parser.getboolean(section, 'Mangrove'))
    except Exception as x:
        print(x)

    print('Parameters:', max_zenith,MaxNights,FOV, MaxRuns, MinProbCut, doplot, Duration, MinDuration, SecondRound,ReducedNside,HRnside,UseGreytime)
    #########################

    # link to the GW map
    name = filename.split('.')[0].split('/')[-1]
    #if('G' in filename):
    #    names = filename.split("_")
    #    name= names[0]

    random.seed()

    RAarray = []
    DECarray = []
    pixlist = []
    ipixlistHR = []
    pixlist1 = []
    ipixlistHR1=[]
    P_GWarray = []
    ObservationTimearray= []
    Round = []

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
    ##########################
    #Add observed pixels to pixlist

    if (PointingFile=='False'):
       print('No pointings were given to be substracted')
    else:
        pixlist,P_GW=SubstractPointings2D(PointingFile,prob,ReducedNside,FOV,pixlist)
        print('Already observed probability =',P_GW)
    ##########################

    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    if(UseGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0,Observatory,MaxNights,Duration,MinDuration)

    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0,Observatory,MaxNights,Duration,MinDuration)

    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            ObsBool,yprob=ZenithAngleCut(prob,nside,ObservationTime,MinProbCut,max_zenith,Observatory.Location,UseGreytime)
            if ObsBool:
                # Round 1
                P_GW,TC,pixlist,ipixlistHR = ComputeProbability2D(prob,highres,radecs,ReducedNside,HRnside,MinProbCut,ObservationTime,Observatory.Location,max_zenith,FOV,name,pixlist,ipixlistHR,counter,dirName,UseGreytime,doplot)
                if ((P_GW <= MinProbCut)and SecondRound):
                    #Try Round 2
                    #print('The minimum probability cut being', MinProbCut * 100, '% is, unfortunately, not reached.')
                    yprob1=highres
                    P_GW, TC, pixlist1,ipixlistHR1 = ComputeProbability2D(prob,yprob1,radecs, nside,ReducedNside,HRnside,PercentCoverage, ObservationTime,Observatory.Location, max_zenith,FOV, name, pixlist1,ipixlistHR1, counter,dirName,UseGreytime,doplot)
                    if ((P_GW <= MinProbCut)):
                        print('Fail')
                    else:
                        Round.append(2)
                        P_GWarray.append(P_GW)
                        RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                        DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                        ObservationTimearray.append(ObservationTime)
                        counter = counter + 1
                elif(P_GW >= MinProbCut):
                    Round.append(1)
                    P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                    RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                    DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                    ObservationTimearray.append(ObservationTime)
                    counter=counter+1
        else:
            break


    print()
    print("===========================================================================================")
    print()
    print("===========================================================================================")
    print()
    print("Total GW probability covered: ", sum(P_GWarray), "Number of runs that fulfill darkness condition  :",
                  len(NightDarkRuns), "Number of effective pointings: ", len(ObservationTimearray))

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,Round], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','Round'])
    return(SuggestedPointings,ObservationTime0)

def BestCandidateonPGal(filename,ObservationTime0,galFile):
    
    # Main Parameters

    ####################
    # ToDo: put these parameters into the 'parameters.ini' file and extract them here
    '''
    cfg = './configs/BestCandidateonPGal.ini'
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWBestGalaxyParameters'

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


    except Exception, x:
        print x

    print('GWBestGalaxyParameters:', max_zenith,MaxNights,FOV, MaxRuns, probCut, MinimumProbCutForCatalogue, doplot, Duration, MinDuration, SecondRound, FulFillReq_Percentage)
    '''
    #########################
    max_zenith = 60
    FOV = 1.5
    MaxRuns = 20 # Maximum number of pointings/runs
    MaxNights = 3
    MinimumProbCutForCatalogue = 0.01
    probCut = 0.05
    doplot=False
    SecondRound=False
    Duration=28
    MinDuration=10
    FulFillReq_Percentage=0.75
    ####################
    
    # load galaxy catalog from local file
    cat = LoadGalaxies(galFile)
    print('done loading galaxies')

    name = filename.split('.')[0].split('/')[-1]
    #if('G' in filename):
    #    names = filename.split("_")
    #    name= names[0]
    
    
    print()
    print('-------------------   NEW LVC EVENT   --------------------')
    print()
    
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


    # correlate GW map with galaxy catalog, retrieve ordered list
    tGals,sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D,MinimumProbCutForCatalogue)
    tGals_aux = tGals
    tGals_aux2 = tGals


    P_GALarray = []
    P_GWarray = []
    ObservationTimearray = []
    ObservationTimearrayNamibia = []
    RAarray = []
    DECarray = []
    alreadysumipixarray1 = []
    alreadysumipixarray2 = []
    Round=[]
    savedcircle = co.SkyCoord([],[], frame='fk5', unit=(u.deg, u.deg))
    savedcircle2 = co.SkyCoord([], [], frame='fk5', unit=(u.deg, u.deg))
    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    
    #Time in UTC

    NightDarkRuns = NightDarkObservation(ObservationTime0,Observatory,MaxNights,Duration,MinDuration)

    totalProb=0.
    
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, max_zenith,Observatory.Location)
            
            if (visible):
                
                visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                visiGals = tGals_aux[visiMask]
                #visiGals = ModifyCatalogue(visiGals, FOV, sum_dP_dV)
                
                mask, minz = FulfillsRequirement(visiGals, max_zenith,FOV,FulFillReq_Percentage,UsePix=False)

                finalGals = visiGals[mask]
                probability = SimpleGWprob(prob,finalGals,alreadysumipixarray1,FOV,nside)
                
                if (probability > probCut):
                    # final galaxies within the FoV
                    if((probability< (2*probCut)) and (sum(P_GWarray)>0.40)and SecondRound): #This notes LIGOVirgo type of signal
                        #print('probability',probability)
                        visible, altaz, tGals_aux2 = VisibleAtTime(ObservationTime, tGals_aux2, max_zenith,Observatory.Location)
                        if (visible):
                            visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                            visiGals2 = tGals_aux2[visiMask]
                            # visiGals = ModifyCatalogue(visiGals, FOV, sum_dP_dV)
                            
                            mask, minz = FulfillsRequirement(visiGals2, max_zenith,FOV,FulFillReq_Percentage,UsePix=False)

                            finalGals2 = visiGals2[mask]
                            probability = SimpleGWprob(prob,finalGals2, alreadysumipixarray2,FOV,nside)
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2, savedcircle2 = ComputeProbBCFOV(prob,ObservationTime,
                                                                                                           finalGals2, visiGals2,
                                                                                                           tGals_aux2, sum_dP_dV,
                                                                                                           alreadysumipixarray2,
                                                                                                           nside, minz,max_zenith, FOV, name,
                                                                                                           savedcircle2, dirName,doplot)
                                
                            RAarray.append(finalGals2['RAJ2000'][:1])
                            DECarray.append(finalGals2['DEJ2000'][:1])
                            Round.append(2)
                    else:
                        #print("\n=================================")
                        #print("TARGET COORDINATES AND DETAILS...")
                        #print("=================================")
                        #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV'][:1])
                        
                        p_gal, p_gw, tGals_aux, alreadysumipixarray1,savedcircle = ComputeProbBCFOV(prob,ObservationTime, finalGals,visiGals, tGals_aux,sum_dP_dV,alreadysumipixarray1, nside, minz,max_zenith, FOV, name,savedcircle,dirName,doplot)
                        RAarray.append(finalGals['RAJ2000'][:1])
                        DECarray.append(finalGals['DEJ2000'][:1])
                        Round.append(1)
                    P_GALarray.append(p_gal)
                    #print('p_gal/sum_dP_dV=',p_gal/sum_dP_dV)
                    P_GWarray.append(p_gw)
                    ObservationTimearray.append(ObservationTime)
                    ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))
                
                
                else:
                    print("Optimal pointing position is: ")
                    print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV'][:1])
                    print("NOT passing the cut on dp_dV > ", probCut)
        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,P_GALarray,Round], names=['Observation Time UTC','RA(deg)','Dec(deg)','PGW','Pgal','Round'])
    return SuggestedPointings,cat


def PGalinFoV(filename,ObservationTime0,PointingFile,galFile,parameters,dirName,Observatory):
    
    # Main Parameters

    #########################
    cfg = parameters
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWBestGalaxyParameters'

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
        Mangrove = (parser.getboolean(section, 'Mangrove'))

    except Exception as x:
        print(x)

    print('GWBestGalaxyParameters:', max_zenith,MaxNights,FOV, MaxRuns, probCut, MinimumProbCutForCatalogue, doplot, Duration, MinDuration, SecondRound, FulFillReq_Percentage,dirName,UseGreytime, Mangrove)
    #########################
    
    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    if not Mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)
    print('done loading galaxies')

    name = filename.split('.')[0].split('/')[-1]
    
    #name = "Default"
    #if ('G' in filename):
    #    names = filename.split("_")
    #    name = names[0]
    
    
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


    # correlate GW map with galaxy catalog, retrieve ordered list
    if not Mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D,MinimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D,MinimumProbCutForCatalogue)

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []

    #########################
    if (PointingFile=='False'):
       tGals = tGals0
       print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray1 = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,FOV,prob,nside)
        MaxRuns = MaxRuns - len(ra)
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)

        #ObservedPointings = Table([time, ra, dec, AlreadyObservedPgw, AlreadyObservedPgal],names=['Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'Covered GW probability','Pgal covered'])
        #print(ObservedPointings)
        print("===========================================================================================")
        print()
        print(
            name, "Total GW probability already covered: ", sumPGW,
            "Total Gal probability already covered: ",
            sumPGAL, "Number of effective pointings already done: ",len(ra))

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
    print('MaxRuns: ',MaxRuns,'MinimumProbCutForCatalogue: ',MinimumProbCutForCatalogue)

    if(UseGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0,Observatory,MaxNights,Duration,MinDuration)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0,Observatory,MaxNights,Duration,MinDuration)
    
    # print('EffectiveRunsTime',len(NightDarkRuns),'being',NightDarkRuns)
    
    
    totalProb=0.
    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, max_zenith,Observatory.Location)
            
            if (visible):
                
                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (max_zenith+FOV)
                visiGals= tGals_aux[visiMask]
                visiGals = ModifyCatalogue(prob,visiGals,FOV,sum_dP_dV,nside)
                
                mask, minz = FulfillsRequirement(visiGals, max_zenith,FOV,FulFillReq_Percentage,UsePix=False)
                if UseGreytime:
                    maskgrey=FulfillsRequirementGreyObservations(ObservationTime,visiGals,Observatory.Location)
                    finalGals=visiGals[mask&maskgrey]
                if not UseGreytime:
                    finalGals = visiGals[mask]
                
                if(finalGals['dp_dV_FOV'][:1] > probCut):
                    # final galaxies within the FoV
                    if ((finalGals['dp_dV_FOV'][:1] < (2 * probCut)) and (sum(P_GWarray) > 0.40) and SecondRound):  # This notes LIGOVirgo type of signal
                        print('probability', finalGals['dp_dV_FOV'][:1])
                        visible, altaz, tGals_aux2 = VisibleAtTime(ObservationTime, tGals_aux2, max_zenith,Observatory.Location)
                        if (visible):
                            visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                            visiGals2 = tGals_aux2[visiMask]
                            visiGals2 = ModifyCatalogue(prob,visiGals2, FOV, sum_dP_dV,nside)
                            
                            mask, minz = FulfillsRequirement(visiGals2, max_zenith,FOV,FulFillReq_Percentage,UsePix=False)

                            finalGals2 = visiGals2[mask]
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,Observatory.Location,finalGals2,False,visiGals2,tGals_aux2,sum_dP_dV, alreadysumipixarray2,nside, minz,max_zenith, FOV,counter,name,dirName,doplot)

                            RAarray.append(np.float('{:3.4f}'.format(np.float(finalGals2['RAJ2000'][:1]))))
                            DECarray.append(np.float('{:3.4f}'.format(np.float(finalGals2['DEJ2000'][:1]))))
                            Round.append(2)

                    else:
                        #print("\n=================================")
                        #print("TARGET COORDINATES AND DETAILS...")
                        #print("=================================")
                        #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                        p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,Observatory.Location,finalGals,False, visiGals,tGals_aux, sum_dP_dV,alreadysumipixarray1,nside, minz,max_zenith, FOV, counter,name, dirName,doplot)
                        RAarray.append(np.float('{:3.4f}'.format(np.float(finalGals['RAJ2000'][:1]))))
                        DECarray.append(np.float('{:3.4f}'.format(np.float(finalGals['DEJ2000'][:1]))))
                        Round.append(1)
                    P_GALarray.append(np.float('{:1.4f}'.format(p_gal)))
                    P_GWarray.append(np.float('{:1.4f}'.format(p_gw)))
                    ObservationTimearray.append(ObservationTime)
                    counter = counter + 1
                    #ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))
            
                else:
                    #print("Optimal pointing position is: ")
                    #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    print("NOT passing the cut on dp_dV_FOV > ",probCut)
        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,P_GALarray,Round], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','Pgal','Round'])
    return SuggestedPointings,cat

def PGalinFoV_PixRegion(filename,ObservationTime0,PointingFile,galFile, parameters,dirName, Observatory):

    # Main Parameters

    ###############################
    # ToDo: put these parameters into the 'parameters.ini' file and extract them here
    cfg = parameters
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWBestGalaxyParameters'

    try:
        max_zenith = int(parser.get(section, 'max_zenith'))
        MaxNights = int(parser.get(section, 'MaxNights'))
        FOV = float(parser.get(section, 'FOV'))
        MaxRuns = int(parser.get(section, 'MaxRuns'))
        probCut = float(parser.get(section, 'probCut'))
        MinimumProbCutForCatalogue = float(parser.get(section, 'MinimumProbCutForCatalogue'))
        doplot = (parser.getboolean(section, 'doplot'))
        FulFillReq_Percentage = float(parser.get(section, 'FulFillReq_Percentage'))
        print('GWBestGalaxyParameters:', max_zenith,MaxNights,FOV, MaxRuns, probCut, MinimumProbCutForCatalogue, doplot, FulFillReq_Percentage)
        NewNside = int(parser.get(section, 'NewNside'))
        PercentCoverage = float(parser.get(section, 'PercentCoverage'))
        Duration = int(parser.get(section, 'Duration'))
        MinDuration = int(parser.get(section, 'MinDuration'))
        UseGreytime = (parser.getboolean(section, 'UseGreytime'))
        Mangrove = (parser.getboolean(section, 'Mangrove'))


    except Exception as x:
        print(x)

    print('GWBestGalaxyParameters:', max_zenith,MaxNights,FOV, MaxRuns, probCut, MinimumProbCutForCatalogue, doplot, Duration, MinDuration, FulFillReq_Percentage, NewNside,UseGreytime)

    ###############################

    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    # galFile='./GLADE_2clean.txt'

    if not Mangrove:
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
    if not Mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D,MinimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D,MinimumProbCutForCatalogue)

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []
    #########################

    if (PointingFile=='False'):
       tGals = tGals0
       print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray1 = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,FOV,prob,nside)
        MaxRuns = MaxRuns - len(ra)
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)

        #ObservedPointings = Table([time, ra, dec, AlreadyObservedPgw, AlreadyObservedPgal],names=['Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'Covered GW probability','Pgal covered'])
        #print(ObservedPointings)
        print("===========================================================================================")
        print()
        print(
            name, "Total GW probability already covered: ", sumPGW,
            "Total Gal probability already covered: ",
            sumPGAL, "Number of effective pointings already done: ",len(ra))
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
    print('MaxRuns: ', MaxRuns, 'MinimumProbCutForCatalogue: ', MinimumProbCutForCatalogue)

    if(UseGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0, Observatory,MaxNights,Duration,MinDuration)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0,Observatory, MaxNights,Duration,MinDuration)

    totalProb = 0.
    n = 0

    ###############################

    # Get the RA & DEC of pixles of the pixels in an enclosed probability region (% precised by PercentCoverage).
    # Reduce these RA DEC to angles in maps with smaller resolution (NewNside)


    pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob, PercentCoverage, NewNside)



    ##############################
    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            if (nextround):
                ObservationTime = NightDarkRuns[j - 1]
                nextround = False
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, max_zenith,
                                                      Observatory.Location)

            if (visible):

                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                visiGals = tGals_aux[visiMask]

                mask, minz = FulfillsRequirement(visiGals, max_zenith, FOV, FulFillReq_Percentage, UsePix=True)

                finalGals = visiGals[mask]
                visiPix = ModifyCataloguePIX(pix_ra1, pix_dec1, ObservationTime, max_zenith, prob, finalGals, FOV,
                                             sum_dP_dV, nside, NewNside, minz,Observatory.Location)

                if (visiPix['PIXFOVPROB'][:1] > probCut):
                    n = n + 1
                    # final galaxies within the FoV

                    # print("\n=================================")
                    # print("TARGET COORDINATES AND DETAILS...")
                    # print("=================================")
                    # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob, ObservationTime, Observatory.Location,
                                                                                               visiPix, True, visiGals,
                                                                                               tGals_aux, sum_dP_dV,
                                                                                               alreadysumipixarray1,
                                                                                               nside, minz, max_zenith,
                                                                                               FOV, counter, name,dirName, doplot)
                    RAarray.append(np.float('{:3.4f}'.format(np.float(visiPix['PIXRA'][:1]))))
                    DECarray.append(np.float('{:3.4f}'.format(np.float(visiPix['PIXDEC'][:1]))))
                    Round.append(1)
                    P_GALarray.append(np.float('{:1.4f}'.format(p_gal)))
                    P_GWarray.append(np.float('{:1.4f}'.format(p_gw)))
                    ObservationTimearray.append(ObservationTime)
                    counter=counter+1

                else:
                    print("\nOptimal pointing position is: ")
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

