from .PointingTools import (NightDarkObservation,
                                  NightDarkObservationwithGreyTime,LoadHealpixMap, Get90RegionPixReduced,
                                  ZenithAngleCut,ComputeProbability2D,
                                  FulfillsRequirement,
                                  VisibleAtTime,LoadGalaxies,CorrelateGalaxies_LVC, SubstractPointings2D,SimpleGWprob,ComputeProbBCFOV,
                                  Tools,LoadGalaxies_SteMgal, CorrelateGalaxies_LVC_SteMass, SubstractPointings,
                                  ModifyCatalogue,FulfillsRequirementGreyObservations,ComputeProbPGALIntegrateFoV,
                                  ModifyCataloguePIX, ObservationParameters, NextWindowTools,
                                  ComputeProbability2D_SelectClusters,GiveProbToGalaxy, LoadGalaxiesSimulation)
from .Observatories import CTASouthObservatory,CTANorthObservatory
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
############################################

#              General definitions              #

############################################

def PGWinFoV(filename,ObservationTime0,PointingFile,parameters,dirName):

    # Main parameters

    obspar = ObservationParameters()
    obspar.from_configfile(parameters)
    
    print(obspar)

    # link to the GW map
    name = filename.split('.')[0].split('/')[-1]

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
    prob = hp.pixelfunc.ud_grade(tprob,obspar.ReducedNside,power=-2)
    nside = obspar.ReducedNside

    highres=hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)

    # Create table for 2D probability at 90% containment
    rapix, decpix,areapix=Get90RegionPixReduced(prob,obspar.PercentCoverage,obspar.ReducedNside)
    radecs= co.SkyCoord(rapix,decpix, frame='fk5', unit=(u.deg, u.deg))

    #Add observed pixels to pixlist
    if (PointingFile=='False'):
       print('No pointings were given to be substracted')
    else:
        pixlist,P_GW=SubstractPointings2D(PointingFile,prob,obspar.ReducedNside,obspar.FOV,pixlist)
        print('Already observed probability =',P_GW)

    #######################################################

    print('----------   NEW FOLLOW-UP ATTEMPT   ----------')
    if(obspar.UseGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0,obspar)

    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0,obspar)

    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < obspar.MaxRuns):
            ObservationTime = NightDarkRuns[j]
            ObsBool,yprob=ZenithAngleCut(prob,nside,ObservationTime,obspar.MinProbCut,obspar.max_zenith,obspar.Location,obspar.UseGreytime)
            if ObsBool:
                # Round 1
                P_GW,TC,pixlist,ipixlistHR = ComputeProbability2D(prob,highres,radecs,obspar.ReducedNside,obspar.HRnside,obspar.MinProbCut,ObservationTime,obspar.Location,obspar.max_zenith,obspar.FOV,name,pixlist,ipixlistHR,counter,dirName,obspar.UseGreytime,obspar.doplot)
                if ((P_GW <= obspar.MinProbCut)and obspar.SecondRound):
                    #Try Round 2
                    #print('The minimum probability cut being', MinProbCut * 100, '% is, unfortunately, not reached.')
                    yprob1=highres
                    P_GW, TC, pixlist1,ipixlistHR1 = ComputeProbability2D(prob,yprob1,radecs, nside,obspar.ReducedNside,obspar.HRnside,obspar.PercentCoverage, ObservationTime,obspar.Location, obspar.max_zenith,FOV, name, pixlist1,ipixlistHR1, counter,dirName,obspar.UseGreytime,obspar.doplot)
                    if ((P_GW <= obspar.MinProbCut)):
                        print('Fail')
                    else:
                        Round.append(2)
                        P_GWarray.append(P_GW)
                        RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                        DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                        ObservationTimearray.append(str(ObservationTime).split('.')[0])
                        counter = counter + 1
                elif(P_GW >= obspar.MinProbCut):
                    Round.append(1)
                    P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                    RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                    DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                    ObservationTimearray.append(str(ObservationTime).split('.')[0])
                    counter=counter+1
        else:
            break


    print()
    #print("===========================================================================================")
    #print()
    print("Total GW probability covered: ", sum(P_GWarray), "Number of runs that fulfill darkness condition  :",
                  len(NightDarkRuns), "Number of effective pointings: ", len(ObservationTimearray))

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,Round], names=['Observation Time UTC','RA[deg]','DEC[deg]','PGW','Round'])
    return(SuggestedPointings,ObservationTime0)

def BestCandidateonPGal(filename,ObservationTime0,galFile):
    # THIS SCRIPT IS Outdated!!! AND NOT USED ANYMORE
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
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,P_GALarray,Round], names=['Observation Time UTC','RA[deg]','Dec[deg]','PGW','Pgal','Round'])
    return SuggestedPointings,cat

def PGalinFoV(filename,ObservationTime0,PointingFile,galFile,parameters,dirName):
    
    # Main Parameters

    obspar = ObservationParameters()
    obspar.from_configfile(parameters)
    print(obspar)

    
    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    if not obspar.Mangrove:
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
    if not obspar.Mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D,obspar.MinimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.MinimumProbCutForCatalogue)

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []

    #########################
    if (PointingFile=='False'):
       tGals = tGals0
       print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray1 = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,obspar.FOV,prob,nside)
        MaxRuns = obspar.MaxRuns - len(ra)
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
    print('MaxRuns: ',MaxRuns,'MinimumProbCutForCatalogue: ', obspar.MinimumProbCutForCatalogue)

    if(obspar.UseGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0,obspar)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0,obspar)
    
    # print('EffectiveRunsTime',len(NightDarkRuns),'being',NightDarkRuns)
    
    
    totalProb=0.
    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, obspar.max_zenith, obspar.Location)
            
            if (visible):
                
                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (obspar.max_zenith+obspar.FOV)
                visiGals= tGals_aux[visiMask]
                visiGals = ModifyCatalogue(prob,visiGals,obspar.FOV,sum_dP_dV,nside)
                
                mask, minz = FulfillsRequirement(visiGals, obspar.max_zenith,obspar.FOV,obspar.FulFillReq_Percentage,UsePix=False)
                if obspar.UseGreytime:
                    maskgrey=FulfillsRequirementGreyObservations(ObservationTime,visiGals,obspar.Location)
                    finalGals=visiGals[mask&maskgrey]
                if not obspar.UseGreytime:
                    finalGals = visiGals[mask]
                
                if(finalGals['dp_dV_FOV'][:1] > obspar.probCut):
                    # final galaxies within the FoV
                    if ((finalGals['dp_dV_FOV'][:1] < (2 * probCut)) and (sum(P_GWarray) > 0.40) and obspar.SecondRound):  # This notes LIGOVirgo type of signal
                        print('probability', finalGals['dp_dV_FOV'][:1])
                        visible, altaz, tGals_aux2 = VisibleAtTime(ObservationTime, tGals_aux2, obspar.max_zenith,obspar.Location)
                        if (visible):
                            visiMask = altaz.alt.value > 90 - (obspar.max_zenith + obspar.FOV)
                            visiGals2 = tGals_aux2[visiMask]
                            visiGals2 = ModifyCatalogue(prob,visiGals2, obspar.FOV, sum_dP_dV,nside)
                            
                            mask, minz = FulfillsRequirement(visiGals2, obspar.max_zenith,obspar.FOV,obspar.FulFillReq_Percentage,UsePix=False)

                            finalGals2 = visiGals2[mask]
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,obspar.Location,finalGals2,False,visiGals2,tGals_aux2,sum_dP_dV, alreadysumipixarray2,nside, minz,obspar.max_zenith, obspar.FOV, counter,name,dirName,obspar.doplot)

                            RAarray.append(np.float('{:3.4f}'.format(np.float(finalGals2['RAJ2000'][:1]))))
                            DECarray.append(np.float('{:3.4f}'.format(np.float(finalGals2['DEJ2000'][:1]))))
                            Round.append(2)

                    else:
                        #print("\n=================================")
                        #print("TARGET COORDINATES AND DETAILS...")
                        #print("=================================")
                        #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                        p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob,ObservationTime,obspar.Location,finalGals,False, visiGals,tGals_aux, sum_dP_dV,alreadysumipixarray1,nside, minz,obspar.max_zenith, obspar.FOV, counter,name, dirName, obspar.doplot)
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
                    print("NOT passing the cut on dp_dV_FOV > ",obspar.probCut)
        else:
            break

    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,P_GALarray,Round], names=['Observation Time UTC','RA[deg]','DEC[deg]','PGW','Pgal','Round'])
    return SuggestedPointings,cat

def PGalinFoV_PixRegion(filename,ObservationTime0,PointingFile,galFile, parameters,dirName):

    # Main parameters from config

    obspar = ObservationParameters()
    obspar.from_configfile(parameters)
    ###############################

    # load galaxy catalog from local file
    # this could be done at the beginning of the night to save time
    # galFile='./GLADE_2clean.txt'

    if not obspar.Mangrove:
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
    if not obspar.Mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D, obspar.MinimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.MinimumProbCutForCatalogue)

    alreadysumipixarray1 = []
    alreadysumipixarray2 = []
    #########################

    if (PointingFile=='False'):
       tGals = tGals0
       print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray1 = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,obspar.FOV,prob,nside)
        MaxRuns = obspar.MaxRuns - len(ra)
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
    print('MaxRuns: ', MaxRuns, 'MinimumProbCutForCatalogue: ', obspar.MinimumProbCutForCatalogue)

    if(obspar.UseGreytime):
        NightDarkRuns = NightDarkObservationwithGreyTime(ObservationTime0, obspar)
    else:
        NightDarkRuns = NightDarkObservation(ObservationTime0,obspar)

    totalProb = 0.
    n = 0

    ###############################

    # Get the RA & DEC of pixles of the pixels in an enclosed probability region (% precised by PercentCoverage).
    # Reduce these RA DEC to angles in maps with smaller resolution (ReducedNside)


    pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob, obspar.PercentCoverage, obspar.ReducedNside)



    ##############################
    counter=0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            if (nextround):
                ObservationTime = NightDarkRuns[j - 1]
                nextround = False
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, obspar.max_zenith,
                                                      obspar.Location)

            if (visible):

                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (obspar.max_zenith + obspar.FOV)
                visiGals = tGals_aux[visiMask]

                mask, minz = FulfillsRequirement(visiGals, obspar.max_zenith, obspar.FOV, obspar.FulFillReq_Percentage, UsePix=True)

                finalGals = visiGals[mask]
                visiPix = ModifyCataloguePIX(pix_ra1, pix_dec1, ObservationTime, obspar.max_zenith, prob, finalGals, obspar.FOV,
                                             sum_dP_dV, nside, obspar.NewNside, minz,obspar.Location)

                if (visiPix['PIXFOVPROB'][:1] > obspar.probCut):
                    n = n + 1
                    # final galaxies within the FoV

                    # print("\n=================================")
                    # print("TARGET COORDINATES AND DETAILS...")
                    # print("=================================")
                    # print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                    p_gal, p_gw, tGals_aux, alreadysumipixarray1 = ComputeProbPGALIntegrateFoV(prob, ObservationTime, obspar.Location,
                                                                                               visiPix, True, visiGals,
                                                                                               tGals_aux, sum_dP_dV,
                                                                                               alreadysumipixarray1,
                                                                                               nside, minz, obspar.max_zenith,
                                                                                               obspar.FOV, counter, name,dirName, obspar.doplot)
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
                    print("NOT passing the cut on dp_dV_FOV > ", obspar.probCut)

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
    print("Name", name, "Total GW probability covered: ", sum(P_GWarray), "Total Gal probability covered: ", sum(P_GALarray),
    "Number of runs that fulfill darkness condition  :", len(NightDarkRuns), "Number of effective pointings: ",
    len(ObservationTimearray))
    return SuggestedPointings,cat,area

    #print("===========================================================================================")
    #print()
    #print()


def PGWonFoV_WindowsfromIRFs(filename, InputChar, TC, parameters, dirName):

    UseObs = InputChar['Observatory']
    run = InputChar['run']
    mergerID = InputChar['MergerID']
    zenith = InputChar['Zenith']
    ObservationTime0 = datetime.datetime.strptime(InputChar['Time'], '%Y-%m-%d %H:%M:%S.%f')

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
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob, obspar.ReducedNside, power=-2)

    nside = obspar.ReducedNside

    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(prob, obspar.PercentCoverage, obspar.ReducedNside)
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
    tstar, tobs = ObtainObservingTimes(totalTime, total_followupDelay, run, mergerID, observatory, obspar.zenith)

    counter = 0
    # Loop where probabilities are computed comes here

    for j in range(0, len(tobs), 1):
        if (j > obspar.MaxRuns):
            break
        else:
            ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar[j])
            # print(ObservationTime)

            ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, obspar.MinProbCut, obspar.max_zenith, observatory.Location,
                                            False)
            # print(ObsBool)
            if ObsBool:
                # Round 1
                # print('Round ',j,'counter',counter)
                P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, obspar.ReducedNside, obspar.HRnside,
                                                                     obspar.MinProbCut, ObservationTime, observatory.Location,
                                                                     obspar.max_zenith, obspar.FOV, name, pixlist, ipixlistHR,
                                                                     counter,
                                                                     dirName, False, obspar.doplot)
                print('P_GW', P_GW)
                if ((P_GW >= obspar.MinProbCut)):
                    P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                    RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                    DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
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


def PGWonFoV_WindowOptimisation(filename, InputChar, TC, parameters, datasetDir, outDir):

    UseObs = InputChar['Observatory']
    run = InputChar['run']
    mergerID = InputChar['MergerID']
    # zenith=InputChar['Zenith']
    
    ObservationTime0 = datetime.datetime.strptime(InputChar['Time'], '%Y-%m-%d %H:%M:%S')
    ObservationTime0 = pytz.utc.localize(ObservationTime0)

    # Main parameters from config
    
    obspar = ObservationParameters()
    obspar.from_configfile(parameters)


    ##################
    path = outDir + '/PointingPlotting/' + run + '_' + mergerID + '/EvolutionPlot/'
    # Observatory
    if UseObs == 'South':
        print('Observed form the', UseObs)
        obspar.name = UseObs
        obspar.Lon = CTASouthObservatory().Lon
        obspar.Lat = CTASouthObservatory().Lat
        obspar.Height = CTASouthObservatory().Height
        obspar.Location = CTASouthObservatory().Location
    else:
        print('Observed from the', UseObs)
        obspar.name = UseObs
        obspar.Lon = CTANorthObservatory().Lon
        obspar.Lat = CTANorthObservatory().Lat
        obspar.Height = CTANorthObservatory().Height
        obspar.Location = CTANorthObservatory().Location

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
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob, obspar.ReducedNside, power=-2)

    nside = obspar.ReducedNside

    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(prob, obspar.PercentCoverage, obspar.ReducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='fk5', unit=(u.deg, u.deg))
    has3D = True

    if (len(distnorm) == 0):
        print("Found a generic map without 3D information")
        # flag the event for special treatment
        has3D = False
    else:
        print("Found a 3D reconstruction")

    # print('----------   NEW FOLLOW-UP ATTEMPT   ----------')

    # Set of delays and slewing times in seconds
    totalTime = 172800  # 48h
    followupDelay = 30  # Minimum delay to start observaiton
    SlewingTime = 210  # Slewing to first position
    interObsSlew = 20  # Slewing time between observations
    total_followupDelay = followupDelay + SlewingTime

    # From the injection time, look for the next window. Time is the time of the first observation
    ObservationTime = ObservationTime0 + datetime.timedelta(seconds=total_followupDelay)

    print("Main info of this scheduling:", totalTime, total_followupDelay, run, mergerID, TC, obspar.Location)

    TotalNights = 2
    counter = 0
    counterTotalPossible = 0
    MaxRuns = 20
    # GrbsensMax = 16384
    AuxTimeNextTry = 900  # 15 minutes.  Allows to check if 15 minutes later there is a better opportunity
    # MinProbCut = 0.005

    # case 1 = TstartNight == False
    # case 2 = "The source is not on the temporal FoV of the instrument"
    for nights in range(0, TotalNights):
        TstartNight = NextWindowTools.NextObservationWindow(ObservationTime, obspar)
        print(ObservationTime.tzinfo)
        print(ObservationTime0.tzinfo.utcoffset)
        print("TstartNight", TstartNight)
        if (TstartNight == False):
            TendNight = False
            print("Night number", nights, " window starts at", TstartNight, "finishes at", TendNight)
            ObsCase = 'NoDarknessFound'
            #print("===== RESULTS ========")
            print('++++ No Darkness Time Found +++++')
            ObservationTimearray.append(ObservationTime)
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
            TendNight = NextWindowTools.EndObservationWindow(TstartNight, obspar)
            print("Night number", nights, " window starts at", TstartNight, "finishes at", TendNight)

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
                ObsBool, yprob = ZenithAngleCut(prob, nside, checkTime, obspar.MinProbCut, obspar.max_zenith, obspar.Location,
                                                usegreytime=False)
                print(checkTime, ObsBool)
                if (ObsBool == True):
                    StartObsTime = checkTime
                    print("The first time the GW is observed is", StartObsTime, "MinProbCut", obspar.MinProbCut, "max_zenith",
                          obspar.max_zenith)
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
                    if (tt > MaxRuns):
                        break
                    if ((TendNight - datetime.timedelta(seconds=interObsSlew) - StartObsTime).total_seconds() <= 0):
                        print("---- NIGHT IS OVER ----")
                        break
                    TotalExposure = int((TendNight - StartObsTime).total_seconds())  # Total Exposure for the night
                    DelayObs = int((StartObsTime - ObservationTime0).total_seconds())

                    print("ObservationTime0", ObservationTime0)
                    print("TotalExposure", TotalExposure)
                    print("DelayObs", DelayObs)
                    print("----------------------------")
                    # TO-DO: Need to have a class to set all this parameters for an observation, so there are less arguments passed.
                    P_GW, TC, ObsExp, ZenIni, ZenEnd, ObsCase, pixlist, ipixlistHR = ComputeProbability2D_SelectClusters(
                        prob, highres, radecs,  TotalExposure, StartObsTime, DelayObs, interObsSlew, obspar, run, mergerID, pixlist, ipixlistHR, counter, datasetDir, outDir, False, False)
                    print("=============")
                    print("P_GW, ObsExp, ZenIni, ZenEnd, ObsCase")
                    print(P_GW, ObsExp, ZenIni, ZenEnd, ObsCase)
                    print("=============")
                    counterTotalPossible = counterTotalPossible + 1
                    # print("PGW", P_GW,'COORDINATES:', TC,'ZENITH CHANGE', ZenIni,'->', ZenEnd)
                    # print("TotalExposure: ",TotalExposure,"DelayObs:",DelayObs, "Observation Number:",counter, "Exposure:", ObsExp)

                    ObservationTimearray.append(str(StartObsTime).split('.')[0])

                    if (ObsCase == 'TimeNotEnoughIte' or ObsCase == 'TimeNotEnough'):
                        StartObsTime = StartObsTime + datetime.timedelta(seconds=ObsExp + interObsSlew + AuxTimeNextTry)
                    else:
                        StartObsTime = StartObsTime + datetime.timedelta(seconds=ObsExp + interObsSlew)

                    # PreDefWindow.append(predefWind[j])
                    # ObservationTime = ObservationTime0 + datetime.timedelta(seconds=tstar)
                    if ((P_GW >= obspar.MinProbCut)):
                        #print("===== RESULTS ========")
                        print('++++ SCHEDULING OBS +++++')
                        P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                        RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                        DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                        Obsarray.append(obspar.name)
                        ZenIniarray.append(np.float('{:1.4f}'.format(np.float(ZenIni))))
                        ZenEndarray.append(np.float('{:1.4f}'.format(np.float(ZenEnd))))
                        ObsBoolarray.append('True')
                        Exposure.append(ObsExp)
                        Delay.append(DelayObs)
                        counter = counter + 1
                    elif ((P_GW <= obspar.MinProbCut) and (P_GW > 0)):  # Although OBS AND NIGHT WAS TRUE
                        #print("===== RESULTS ========")
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
                        #print("===== RESULTS ========")
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
                #print("===== RESULTS ========")
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
        # nights = nights + 1
    print()
    print("===========================================================================================")
    print()
    # List of suggested pointings
    print("TOTAL POSSIBLE: ", counterTotalPossible, "DONE: ", counter, "Probability: ",np.sum(P_GWarray))
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

    print('GWGalaxyProbabilityIntegrated - Parameters:', max_zenith, MaxNights, FOV, MaxRuns, probCut,
          MinimumProbCutForCatalogue, doplot,
          Duration, MinDuration, SecondRound, FulFillReq_Percentage, dirName, UseGreytime)

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

    name = filename.split('.')[0].split('/')[-1]
    # correlate GW map with galaxy catalog, retrieve ordered list
    tGals, sum_dP_dV = GiveProbToGalaxy(prob, cat, distance, Edistance_max, Edistance_min, MinimumProbCutForCatalogue)
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
    print('MaxRuns: ', MaxRuns, 'MinimumProbCutForCatalogue: ', MinimumProbCutForCatalogue)
    # SlewingTime=datetime.timedelta(seconds=210)
    ObservationTime = ObservationTime0
    time = NextWindowTools.NextObservationWindow(ObservationTime, observatory)
    WindowDurations = InputObservationList['Interval']
    # WindowDurations = [15, 17, 20, 23, 27, 33, 40, 50, 64, 85,119,178,296,595,1905]
    NightDarkRuns = NextWindowTools.CheckWindowCreateArray(time, observatory, WindowDurations)

    # print('EffectiveRunsTime',len(NightDarkRuns),'being',NightDarkRuns)

    predefWind = InputObservationList['pointingNumber']
    totalProb = 0.
    counter = 0
    for j in range(0, len(NightDarkRuns)):
        if (len(ObservationTimearray) < MaxRuns):
            ObservationTime = NightDarkRuns[j]
            visible, altaz, tGals_aux = VisibleAtTime(ObservationTime, tGals_aux, max_zenith, observatory.Location)

            if (visible):

                # select galaxies within the slightly enlarged visiblity window
                visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                visiGals = tGals_aux[visiMask]
                visiGals = ModifyCatalogue(prob, visiGals, FOV, sum_dP_dV, nside)

                mask, minz = FulfillsRequirement(visiGals, max_zenith, FOV, FulFillReq_Percentage, UsePix=False)

                finalGals = visiGals[mask]

                if (finalGals['dp_dV_FOV'][:1] > probCut):
                    # final galaxies within the FoV
                    if ((finalGals['dp_dV_FOV'][:1] < (2 * probCut)) and (
                            sum(P_GWarray) > 0.40) and SecondRound):  # This notes LIGOVirgo type of signal
                        print('probability', finalGals['dp_dV_FOV'][:1])
                        visible, altaz, tGals_aux2 = VisibleAtTime(ObservationTime, tGals_aux2, max_zenith,
                                                                   observatory.Location)
                        if (visible):
                            visiMask = altaz.alt.value > 90 - (max_zenith + FOV)
                            visiGals2 = tGals_aux2[visiMask]
                            visiGals2 = ModifyCatalogue(prob, visiGals2, FOV, sum_dP_dV, nside)

                            mask, minz = FulfillsRequirement(visiGals2, max_zenith, FOV, FulFillReq_Percentage,
                                                             UsePix=False)

                            finalGals2 = visiGals2[mask]
                            p_gal, p_gw, tGals_aux2, alreadysumipixarray2 = ComputeProbPGALIntegrateFoV(prob,
                                                                                                        ObservationTime,
                                                                                                        observatory.Location,
                                                                                                        finalGals2,
                                                                                                        False,
                                                                                                        visiGals2,
                                                                                                        tGals_aux2,
                                                                                                        sum_dP_dV,
                                                                                                        alreadysumipixarray2,
                                                                                                        nside, minz,
                                                                                                        max_zenith, FOV,
                                                                                                        counter, name,
                                                                                                        dirName, doplot)

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
                                                                                                   observatory.Location,
                                                                                                   finalGals, False,
                                                                                                   visiGals, tGals_aux,
                                                                                                   sum_dP_dV,
                                                                                                   alreadysumipixarray1,
                                                                                                   nside, minz,
                                                                                                   max_zenith, FOV,
                                                                                                   counter, name,
                                                                                                   dirName, doplot)
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
