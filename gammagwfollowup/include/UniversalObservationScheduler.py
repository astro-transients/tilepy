from .PointingTools import (NightDarkObservation,
                                  NightDarkObservationwithGreyTime,
                                  LoadHealpixMap,LoadGalaxies,CorrelateGalaxies_LVC,
                                  Get90RegionPixReduced,
                                  SubstractPointings2D,ZenithAngleCut,ComputeProbability2D,
                                  VisibleAtTime, FulfillsRequirement, SimpleGWprob,ComputeProbBCFOV,
                                  Tools,LoadGalaxies_SteMgal, CorrelateGalaxies_LVC_SteMass, SubstractPointings,
                                  ModifyCatalogue,FulfillsRequirementGreyObservations,ComputeProbPGALIntegrateFoV,
                                  ModifyCataloguePIX, ObservationParameters,NextWindowTools)
import random
import numpy as np
from astropy.table import Table
import astropy.coordinates as co
from astropy import units as u
import healpy as hp
import datetime
from six.moves import configparser
import six
if six.PY2:
  ConfigParser = configparser.SafeConfigParser
else:
  ConfigParser = configparser.ConfigParser
############################################

#              General definitions              #

############################################

def PGWinFoV_NObs(filename, ObservationTime0, PointingsFile, parameters, dirName, ObsArray):

    #Finding the start time for each observatory and checking if it's now    
    FirstDark = np.full(len(ObsArray), False, dtype=bool)
    FirstDark_Flag = np.full(len(ObsArray), False, dtype=bool)
    obs_time = ObservationTime0
    ObsFirstTime = []
    ObsParameters = []

    j = 0
    for obspar1 in ObsArray:

        globals()[obspar1] = ObservationParameters.from_configfile(parameters[j])
        ObsParameters.append(globals()[obspar1])

        dark_at_start =False
        if ObsParameters[j].UseGreytime:
          dark_at_start = Tools.IsDarkness(obs_time, ObsParameters[j])
        if not ObsParameters[j].UseGreytime:
          dark_at_start = Tools.IsGreyness(obs_time, ObsParameters[j])
        FirstDark[j] = dark_at_start

        if FirstDark[j] == True:
          FirstDark_Flag[j] = True
          ObsFirstTime[j] = FirstDark[j]
        else:
          ObsFirstTime1 = NextWindowTools.NextObservationWindow(time = obs_time,obsSite=ObsParameters[j])
          ObsFirstTime.append(ObsFirstTime1)
          if ObsFirstTime1 < (obs_time + datetime.timedelta(hours=24)):
            FirstDark_Flag[j] = True

        j+=1


    #Checking which observatories are availabe for observations and saving their start time
    ActiveObsStart = []
    ActiveObs = []
    SameNight = np.full(len(ObsArray), False, dtype=bool)

    j = 0
    for obspar in ObsArray:
      if FirstDark_Flag[j]:
        ActiveObsStart.append(ObsFirstTime[j])
        ActiveObs.append(ObsParameters[j])
        SameNight[j] = True
      j+=1
    print(ActiveObs[0].name, ActiveObs[1].name, ActiveObs[2].name)


    #Sorting observatories according to their first obsevation time available
    NewActiveObsStart = np.sort(ActiveObsStart)
    NewActiveObs = ActiveObs
    ind = np.argsort(ActiveObsStart)
    ind = np.array(ind)
    NewActiveObs = np.take(ActiveObs, ind)

    #START
#################################################################################################################################################
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
    ObsName = []
#################################################################################################################################################
    obspar = ObsParameters[0]
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob,obspar.ReducedNside,power=-2)
    nside = obspar.ReducedNside
    highres=hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix,areapix=Get90RegionPixReduced(prob,obspar.PercentCoverage,obspar.ReducedNside)
    radecs= co.SkyCoord(rapix,decpix, frame='fk5', unit=(u.deg, u.deg))
    #Add observed pixels to pixlist
    if (PointingsFile=='False'):
       print('No pointings were given to be substracted')
    else:
        pixlist,P_GW=SubstractPointings2D(PointingsFile,prob,obspar.ReducedNside,obspar.FOV,pixlist)
        print('Already observed probability =',P_GW)
#################################################################################################################################################
    ITERATION_OBS = 0
    TIME_MIN_ALL = []
    TIME_MIN = obs_time + datetime.timedelta(hours=12)
    NewActiveObsTime = NewActiveObsStart
    NUMBER_OBS = np.zeros(len(NewActiveObs))
#################################################################################################################################################
    counter=0
    #print(SameNight)
    #print(NewActiveObs[0].name, NewActiveObs[1].name, NewActiveObs[2].name)
    #print(NewActiveObsTime)
    i = 0
    while (i < 500) & any(SameNight):
      for j in range(len(NewActiveObs)):
        obspar = NewActiveObs[j]
        #print(j)
        #print(NewActiveObs[0].name)
        #print(obspar.name)
        ObservationTime = NewActiveObsTime[j]
        if ITERATION_OBS == len(ObsArray):
          TIME_MIN_ALL = []
          ITERATION_OBS = 0

        ITERATION_OBS += 1

        if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
          ObsBool,yprob=ZenithAngleCut(prob,nside,ObservationTime,obspar.MinProbCut,obspar.max_zenith,obspar.Location,obspar.UseGreytime)
          if ObsBool:
            # Round 1
            P_GW,TC,pixlist,ipixlistHR = ComputeProbability2D(prob,highres,radecs,obspar.ReducedNside,obspar.HRnside,obspar.MinProbCut,ObservationTime,obspar.Location,obspar.max_zenith,obspar.FOV,name,pixlist,ipixlistHR,counter,dirName,obspar.UseGreytime,obspar.doplot)
            #print(P_GW, obspar.name)
            if ((P_GW <= obspar.MinProbCut) and obspar.SecondRound):
                #Try Round 2
                #print('The minimum probability cut being', MinProbCut * 100, '% is, unfortunately, not reached.')
                yprob1=highres
                P_GW, TC, pixlist1,ipixlistHR1 = ComputeProbability2D(prob,yprob1,radecs, nside,obspar.ReducedNside,obspar.HRnside,obspar.PercentCoverage, ObservationTime,obspar.Location, obspar.max_zenith,obspar.FOV, name, pixlist1,ipixlistHR1, counter,dirName,obspar.UseGreytime,obspar.doplot)
                if ((P_GW <= obspar.MinProbCut)):
                  print('Fail')
                else:
                  Round.append(2)
                  P_GWarray.append(P_GW)
                  RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                  DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                  ObservationTimearray.append(ObservationTime)
                  ObsName.append(obspar.name)
                  counter = counter + 1
                  

            elif(P_GW >= obspar.MinProbCut):
              Round.append(1)
              P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
              RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
              DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
              ObservationTimearray.append(ObservationTime)
              ObsName.append(obspar.name)
              counter=counter+1
              


          #HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
          NewActiveObsTime[j] = NewActiveObsTime[j] + datetime.timedelta(minutes=30)


          #HERE WE DETERMINE IF WE ARE STILL IN THE SAME NIGHT FOR THIS OBSERVATORY
          if (NewActiveObsTime[j] > Tools.NextSunrise(NewActiveObsStart[j], NewActiveObs[j])) | (NewActiveObsStart[j] > Tools.NextMoonrise(obs_time, NewActiveObs[j])):
            SameNight[j] = False

          NUMBER_OBS[j] += 1

          if SameNight[j]:
            TIME_MIN = NewActiveObsTime[j]
            TIME_MIN_ALL.append(TIME_MIN)
            TIME_MIN = np.min(TIME_MIN_ALL)
          else: 
            TIME_MIN = TIME_MIN + datetime.timedelta(hours=12)
      
      i+=1

    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,Round, ObsName], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','Round','ObsName'])


    return SuggestedPointings, ObsParameters


def PGalinFoV_NObs(filename,ObservationTime0,PointingFile,galFile, parameters,dirName, ObsArray):
    
    #Finding the start time for each observatory and checking if it's now
    FirstDark = np.full(len(ObsArray), False, dtype=bool)
    FirstDark_Flag = np.full(len(ObsArray), False, dtype=bool)
    obs_time = ObservationTime0
    ObsFirstTime = []
    ObsParameters = []

    j = 0
    for obspar1 in ObsArray:
        globals()[obspar1] = ObservationParameters.from_configfile(parameters[j])
        ObsParameters.append(globals()[obspar1])

        dark_at_start =False
        if ObsParameters[j].UseGreytime:
          dark_at_start = Tools.IsDarkness(obs_time, ObsParameters[j])
        if not ObsParameters[j].UseGreytime:
          dark_at_start = Tools.IsGreyness(obs_time, ObsParameters[j])
        FirstDark[j] = dark_at_start
        if FirstDark[j] == True:
          FirstDark_Flag[j] = True
          ObsFirstTime[j] = FirstDark[j]
        else:
          ObsFirstTime1 = NextWindowTools.NextObservationWindow(time = obs_time,obsSite=ObsParameters[j])
          ObsFirstTime.append(ObsFirstTime1)
          if ObsFirstTime1 < (obs_time + datetime.timedelta(hours=24)):
            FirstDark_Flag[j] = True
        j+=1


    #Checking which observatories are availabe for observations and saving their start time
    ActiveObsStart = []
    ActiveObs = []
    SameNight = np.full(len(ObsArray), False, dtype=bool)

    j = 0
    for obspar in ObsArray:
      if FirstDark_Flag[j]:
        ActiveObsStart.append(ObsFirstTime[j])
        ActiveObs.append(ObsParameters[j])
        SameNight[j] = True
      j+=1
    print(ActiveObs[0].name, ActiveObs[1].name, ActiveObs[2].name)


    #Sorting observatories according to their first obsevation time available
    NewActiveObsStart = np.sort(ActiveObsStart)
    NewActiveObs = ActiveObs
    ind = np.argsort(ActiveObsStart)
    ind = np.array(ind)
    NewActiveObs = np.take(ActiveObs, ind)

    #START
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
    obspar = ObsParameters[0]
    print('Loading map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    #load galaxu catalogue
    if not obspar.Mangrove:
        cat = LoadGalaxies(galFile)
    else:
        cat = LoadGalaxies_SteMgal(galFile)
    print('done loading galaxies')

    # correlate GW map with galaxy catalog, retrieve ordered list
    if not obspar.Mangrove:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D,obspar.MinimumProbCutForCatalogue)
    else:
        tGals0, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D, obspar.MinimumProbCutForCatalogue)



    #Add observed pixels to pixlist
    if (PointingFile=='False'):
       tGals = tGals0
       print('No pointings were given to be substracted')
    else:
        # tGals_aux = tGals
        ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal,alreadysumipixarray1 = SubstractPointings(PointingFile, tGals0,alreadysumipixarray1,sum_dP_dV,obspar.FOV,prob,nside)
        MaxRuns = obspar.MaxRuns - len(ra)
        sumPGW = sum(AlreadyObservedPgw)
        sumPGAL = sum(AlreadyObservedPgal)

    tGals_aux = tGals
    tGals_aux2 = tGals
#################################################################################################################################################
    ITERATION_OBS = 0
    TIME_MIN_ALL = []
    TIME_MIN = obs_time + datetime.timedelta(hours=12)
    NewActiveObsTime = NewActiveObsStart
    NUMBER_OBS = np.zeros(len(NewActiveObs))
#################################################################################################################################################

    totalProb=0.
    counter=0
    #print(SameNight)
    #print(NewActiveObs[0].name, NewActiveObs[1].name, NewActiveObs[2].name)
    #print(NewActiveObsTime)
    i = 0
    while (i < 500) & any(SameNight):
      for j in range(len(NewActiveObs)):
        obspar = NewActiveObs[j]
        ObservationTime = NewActiveObsTime[j]
        if ITERATION_OBS == len(ObsArray):
          TIME_MIN_ALL = []
          ITERATION_OBS = 0

        ITERATION_OBS += 1

        if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
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
              
              if(finalGals['dp_dV_FOV'][:1] > obspar.MinProbCut):
                  # final galaxies within the FoV
                  if ((finalGals['dp_dV_FOV'][:1] < (2 * obspar.MinProbCut)) and (sum(P_GWarray) > 0.40) and obspar.SecondRound):  # This notes LIGOVirgo type of signal
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
                  ObsName.append(obspar.name)
                  counter = counter + 1
                  #ObservationTimearrayNamibia.append(Tools.UTCtoNamibia(ObservationTime))
          
              else:
                  #print("Optimal pointing position is: ")
                  #print(finalGals['RAJ2000', 'DEJ2000', 'Bmag', 'Dist', 'Alt', 'dp_dV','dp_dV_FOV'][:1])
                  print(finalGals['dp_dV_FOV'][:1])
                  print("NOT passing the cut on dp_dV_FOV > ",obspar.MinProbCut)
            

          #HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
          NewActiveObsTime[j] = NewActiveObsTime[j] + datetime.timedelta(minutes=30)


          #HERE WE DETERMINE IF WE ARE STILL IN THE SAME NIGHT FOR THIS OBSERVATORY
          if (NewActiveObsTime[j] > Tools.NextSunrise(NewActiveObsStart[j], NewActiveObs[j])) | (NewActiveObsStart[j] > Tools.NextMoonrise(obs_time, NewActiveObs[j])):
            SameNight[j] = False

          NUMBER_OBS[j] += 1

          if SameNight[j]:
            TIME_MIN = NewActiveObsTime[j]
            TIME_MIN_ALL.append(TIME_MIN)
            TIME_MIN = np.min(TIME_MIN_ALL)
          else: 
            TIME_MIN = TIME_MIN + datetime.timedelta(hours=12)
      
      i+=1

    SuggestedPointings = Table([ObservationTimearray,RAarray,DECarray,P_GWarray,P_GALarray,Round, ObsName], names=['Observation Time UTC','RA(deg)','DEC(deg)','PGW','Pgal','Round', 'ObsName'])
    return SuggestedPointings, cat, ObsParameters


def PGWinFoV_NObs_Simulation(filename, ObservationTime0, PointingsFile, parameters, dirName, ObsArray):
    
    #TODO: Modify this.. NOT WORKING
    # Finding the start time for each observatory and checking if it's now
    FirstDark = np.full(len(ObsArray), False, dtype=bool)
    FirstDark_Flag = np.full(len(ObsArray), False, dtype=bool)
    obs_time = ObservationTime0
    ObsFirstTime = []
    ObsParameters = []
    
    j = 0
    for obspar1 in ObsArray:
        
        globals()[obspar1] = ObservationParameters.from_configfile(parameters[j])
        ObsParameters.append(globals()[obspar1])
        
        dark_at_start = False
        if ObsParameters[j].UseGreytime:
            dark_at_start = Tools.IsDarkness(obs_time, ObsParameters[j])
        if not ObsParameters[j].UseGreytime:
            dark_at_start = Tools.IsGreyness(obs_time, ObsParameters[j])
        FirstDark[j] = dark_at_start
        
        if FirstDark[j] == True:
            FirstDark_Flag[j] = True
            ObsFirstTime[j] = FirstDark[j]
        else:
            ObsFirstTime1 = NextWindowTools.NextObservationWindow(time=obs_time, obsSite=ObsParameters[j])
            ObsFirstTime.append(ObsFirstTime1)
            if ObsFirstTime1 < (obs_time + datetime.timedelta(hours=24)):
                FirstDark_Flag[j] = True
        
        j += 1
    
    # Checking which observatories are availabe for observations and saving their start time
    ActiveObsStart = []
    ActiveObs = []
    SameNight = np.full(len(ObsArray), False, dtype=bool)
    
    j = 0
    for obspar in ObsArray:
        if FirstDark_Flag[j]:
            ActiveObsStart.append(ObsFirstTime[j])
            ActiveObs.append(ObsParameters[j])
            SameNight[j] = True
        j += 1
    print(ActiveObs[0].name, ActiveObs[1].name, ActiveObs[2].name)
    
    # Sorting observatories according to their first obsevation time available
    NewActiveObsStart = np.sort(ActiveObsStart)
    NewActiveObs = ActiveObs
    ind = np.argsort(ActiveObsStart)
    ind = np.array(ind)
    NewActiveObs = np.take(ActiveObs, ind)

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
    obspar = ObsParameters[0]
    tprob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    prob = hp.pixelfunc.ud_grade(tprob, obspar.ReducedNside, power=-2)
    nside = obspar.ReducedNside
    highres = hp.pixelfunc.ud_grade(prob, obspar.HRnside, power=-2)
    # Create table for 2D probability at 90% containment
    rapix, decpix, areapix = Get90RegionPixReduced(prob, obspar.PercentCoverage, obspar.ReducedNside)
    radecs = co.SkyCoord(rapix, decpix, frame='fk5', unit=(u.deg, u.deg))
    # Add observed pixels to pixlist
    if (PointingsFile == 'False'):
        print('No pointings were given to be substracted')
    else:
        pixlist, P_GW = SubstractPointings2D(PointingsFile, prob, obspar.ReducedNside, obspar.FOV, pixlist)
        print('Already observed probability =', P_GW)
    #################################################################################################################################################
    ITERATION_OBS = 0
    TIME_MIN_ALL = []
    TIME_MIN = obs_time + datetime.timedelta(hours=12)
    NewActiveObsTime = NewActiveObsStart
    NUMBER_OBS = np.zeros(len(NewActiveObs))
    #################################################################################################################################################
    counter = 0
    # print(SameNight)
    # print(NewActiveObs[0].name, NewActiveObs[1].name, NewActiveObs[2].name)
    # print(NewActiveObsTime)
    i = 0
    while (i < 500) & any(SameNight):
        for j in range(len(NewActiveObs)):
            obspar = NewActiveObs[j]
            # print(j)
            # print(NewActiveObs[0].name)
            # print(obspar.name)
            ObservationTime = NewActiveObsTime[j]
            if ITERATION_OBS == len(ObsArray):
                TIME_MIN_ALL = []
                ITERATION_OBS = 0
            
            ITERATION_OBS += 1
            
            if (TIME_MIN >= NewActiveObsTime[j]) & SameNight[j]:
                ObsBool, yprob = ZenithAngleCut(prob, nside, ObservationTime, obspar.MinProbCut, obspar.max_zenith,
                                                obspar.Location, obspar.UseGreytime)
                if ObsBool:
                    # Round 1
                    P_GW, TC, pixlist, ipixlistHR = ComputeProbability2D(prob, highres, radecs, obspar.ReducedNside,
                                                                         obspar.HRnside, obspar.MinProbCut,
                                                                         ObservationTime, obspar.Location,
                                                                         obspar.max_zenith, obspar.FOV, name, pixlist,
                                                                         ipixlistHR, counter, dirName,
                                                                         obspar.UseGreytime, obspar.doplot)
                    # print(P_GW, obspar.name)
                    if ((P_GW <= obspar.MinProbCut) and obspar.SecondRound):
                        # Try Round 2
                        # print('The minimum probability cut being', MinProbCut * 100, '% is, unfortunately, not reached.')
                        yprob1 = highres
                        P_GW, TC, pixlist1, ipixlistHR1 = ComputeProbability2D(prob, yprob1, radecs, nside,
                                                                               obspar.ReducedNside, obspar.HRnside,
                                                                               obspar.PercentCoverage, ObservationTime,
                                                                               obspar.Location, obspar.max_zenith,
                                                                               obspar.FOV, name, pixlist1, ipixlistHR1,
                                                                               counter, dirName, obspar.UseGreytime,
                                                                               obspar.doplot)
                        if ((P_GW <= obspar.MinProbCut)):
                            print('Fail')
                        else:
                            Round.append(2)
                            P_GWarray.append(P_GW)
                            RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                            DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                            ObservationTimearray.append(ObservationTime)
                            ObsName.append(obspar.name)
                            counter = counter + 1
                    
                    
                    elif (P_GW >= obspar.MinProbCut):
                        Round.append(1)
                        P_GWarray.append(np.float('{:1.4f}'.format(np.float(P_GW))))
                        RAarray.append(np.float('{:3.4f}'.format(np.float(TC.ra.deg))))
                        DECarray.append(np.float('{:3.4f}'.format(np.float(TC.dec.deg))))
                        ObservationTimearray.append(ObservationTime)
                        ObsName.append(obspar.name)
                        counter = counter + 1
                
                # HERE WE DETERMINE THE OBSERVATION DURATION ... FOR NOW WE USE 30 MINS FOR ALL
                NewActiveObsTime[j] = NewActiveObsTime[j] + datetime.timedelta(minutes=30)
                
                # HERE WE DETERMINE IF WE ARE STILL IN THE SAME NIGHT FOR THIS OBSERVATORY
                if (NewActiveObsTime[j] > Tools.NextSunrise(NewActiveObsStart[j], NewActiveObs[j])) | (
                        NewActiveObsStart[j] > Tools.NextMoonrise(obs_time, NewActiveObs[j])):
                    SameNight[j] = False
                
                NUMBER_OBS[j] += 1
                
                if SameNight[j]:
                    TIME_MIN = NewActiveObsTime[j]
                    TIME_MIN_ALL.append(TIME_MIN)
                    TIME_MIN = np.min(TIME_MIN_ALL)
                else:
                    TIME_MIN = TIME_MIN + datetime.timedelta(hours=12)
        
        i += 1
    
    SuggestedPointings = Table([ObservationTimearray, RAarray, DECarray, P_GWarray, Round, ObsName],
                               names=['Observation Time UTC', 'RA(deg)', 'DEC(deg)', 'PGW', 'Round', 'ObsName'])
    
    return SuggestedPointings, ObsParameters

