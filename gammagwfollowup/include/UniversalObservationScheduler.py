from .PointingTools import (NightDarkObservation,
                                  NightDarkObservationwithGreyTime,
                                  LoadHealpixMap,LoadGalaxies,CorrelateGalaxies_LVC,
                                  Get90RegionPixReduced,
                                  SubstractPointings2D,ZenithAngleCut,ComputeProbability2D,
                                  VisibleAtTime, FulfillsRequirement, SimpleGWprob,ComputeProbBCFOV,
                                  Tools,LoadGalaxies_SteMgal, CorrelateGalaxies_LVC_SteMass, SubstractPointings,
                                  ModifyCatalogue,FulfillsRequirementGreyObservations,ComputeProbPGALIntegrateFoV,
                                  ModifyCataloguePIX, ObservationParameters)
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

def PGWinFoV_NObs(filename, ObservationTime0, PointingsFile, parameters, dirName, ObsArray):
    
    FirstDark = np.full(len(ObsArray), False, dtype=bool)
    FirstDark_Flag = np.full(len(ObsArray), False, dtype=bool)
    obs_time = ObservationTime0
    obs = np.zeros(len(ObsArray))

    j = 0
    for obspar in ObsArray:
        globals()[obspar] = ObservationParameters.from_configfile(parameters[j])

        dark_at_start =False
        if globals()[obspar].UseGreytime:
          dark_at_start = Tools.IsDarkness(obs_time, globals()[obspar])
        if not globals()[obspar].UseGreytime:
          dark_at_start = Tools.IsGreyness(obs_time, globals()[obspar])
        FirstDark[j] = dark_at_start

        j+=1

    return None
