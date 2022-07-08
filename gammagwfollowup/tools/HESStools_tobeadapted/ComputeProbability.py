from gammagwfollowup.include.PointingTools import (LoadHealpixMap,CorrelateGalaxies_LVC,SubstractPointings2D,
                                  SubstractPointings, LoadGalaxies)
#from PointingPlotting import *
import argparse
import healpy as hp


#Example python tools/CompareTwoSchedules.py G297595_LALInference.fits -file1 SuggestedPointings_GWOptimisation.txt  -file2 SuggestedPointings_GWOptimisation1.txt

parser = argparse.ArgumentParser(description='Plot together two pointing schedules')
parser.add_argument('GWname', metavar='GWmap',
                    help='the FITS file with the GW localization, e.g. Bayestar.fits.gz')
parser.add_argument('-file', metavar='PointingFile',
                    help='list with 1st subset of pointings\
                                       (RAarray DECarray )')
parser.add_argument('-galaxies', metavar='GLADE.txt', default="../GLADE.txt", help='galaxy catalog (default: %(default)s)')



args = parser.parse_args()


GWFile=args.GWname
PointingsFile=args.file
galFile=args.galaxies

print("===========================================================================================")
print("Starting the computation of covered probability from the following files\n")
print('Loading map from ', GWFile)
print("Filename : ",PointingsFile)


prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(GWFile)
npix = len(prob)
nside = hp.npix2nside(npix)

cat = LoadGalaxies(galFile)

has3D = True
if (len(distnorm) == 0):
    print("Found a generic map without 3D information")
    # flag the event for special treatment
    has3D = False
else:
    print("Found a 3D reconstruction")

alreadysumipixarray = []
FOV = 2.0
MinimumProbCutForCatalogue = 0.01

tGals0, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D, MinimumProbCutForCatalogue)
alreadysumipixarray,AlreadyObservedPgw = SubstractPointings2D(PointingsFile,prob, nside, FOV,alreadysumipixarray)
#ra, dec, tGals, AlreadyObservedPgw, AlreadyObservedPgal, alreadysumipixarray = SubstractPointings(PointingsFile, tGals0, alreadysumipixarray, sum_dP_dV, FOV, prob,nside)

#sumPGW = sum(AlreadyObservedPgw)
#sumPGAL = sum(AlreadyObservedPgal)

print('Total PGW covered: ', AlreadyObservedPgw) #sumPGW)
#print('Total PGAL covered: ',sumPGAL)