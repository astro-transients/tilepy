from tilepy.tools.VisualizationTools import Pretty_Plot, LocateSource
from astropy.coordinates import SkyCoord
from astropy import units as u
from tilepy.include.PointingTools import (LoadHealpixMap, LoadGalaxies, CorrelateGalaxies_LVC, ObservationParameters)
import numpy as np


#adjust the following as desired
alertType = 'gw'
configDir = '/Users/hashkar/Desktop/UniversalScheduler/config'
datasetDir =  "/Users/hashkar/Desktop/Universa/Users/hashkar/Desktop/UniversalScheduler/call_vistool.pylScheduler/dataset/"
outDir = '/Users/hashkar/Desktop/UniversalScheduler/output'
galcatName = "converted_GLADE.h5"
pointingsFile = None
locCut = 500

url = 'https://dcc.ligo.org/public/0119/P1500071/007/927563_lalinference.fits.gz'
ObsArray = ['CTAN', 'CTAS']
obsTime = '2023-03-15 10:30:10'
parameters = []
obsparameters = []

name = 'P1500071_CTAO'
PointingsFile1 = '/Users/hashkar/Desktop/UniversalScheduler/output/P1500071/PGalinFoV_NObs/SuggestedPointings_GWOptimisation.txt'
dirName = 'test_paper'




for i in ObsArray:
    parameters.append("./config/FollowupParameters_%s.ini" % i)
print("===========================================================================================")
print('parameters', parameters)
obsparameters = []

for j in range(len(parameters)):
    obspar = ObservationParameters()
    obspar.add_parsed_args(url, obsTime, datasetDir, galcatName, outDir, pointingsFile, alertType, locCut)
    obspar.from_configfile(parameters[j])
    obsparameters.append(obspar)
print(obspar)


galaxies = obspar.datasetDir + obspar.galcatName


cat = LoadGalaxies(galaxies)
prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(url)
tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, True, obspar.minimumProbCutForCatalogue)



center = SkyCoord(195, 15, unit='deg', frame='icrs')
radius = '20 deg'

Pretty_Plot(url, name, PointingsFile1, dirName, obspar, tGals0, center, radius)



##################################
ra = ['08:25:08.857']
dec = ['-24:43:45.91']

filename = 'Bilby.offline0.multiorder-2.fits'

coordinates = SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
LocateSource(filename, coordinates.ra, coordinates.dec, PercentCov=90)
