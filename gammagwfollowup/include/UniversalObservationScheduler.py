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
    for obspar in ObsArray:

        globals()[obspar] = ObservationParameters.from_configfile(parameters[j])
        ObsParameters.append(globals()[obspar])

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
      print(FirstDark_Flag[j])
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




    return None
