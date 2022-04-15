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

    #Sorting observatories according to their first obsevation time available
    NewActiveObsStart = np.sort(ActiveObsStart)
    NewActiveObs = ActiveObs
    ind = np.argsort(ActiveObsStart)
    for i in range(len(ActiveObs)): #THIS IS NOT OPTIMAL AND NEEDS CHANGING
      NewActiveObs[i] = ActiveObs[ind[i]]

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
    while (i < 50) & any(SameNight):
      for j in range(len(NewActiveObs)):
        obspar = NewActiveObs[j]
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
            if ((P_GW <= obspar.MinProbCut)and obspar.SecondRound):
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
