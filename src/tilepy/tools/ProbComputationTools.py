from tilepy.include.PointingTools import (LoadHealpixMap, CorrelateGalaxies_LVC, SubstractPointings2D,
                                          SubstractPointings, LoadGalaxies)

import healpy as hp


def ComputeCoveredProbability(GWFile, PointingsFile, galFile="../dataset/GLADEplus.h5", FOV=2.0, MinimumProbCutForCatalogue=0.01):

    print("===========================================================================================")
    print("Starting the computation of covered probability from the following files\n")
    print('Loading map from ', GWFile)
    print("Filename : ", PointingsFile)
    print("===========================================================================================")

    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(
        GWFile)
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

    tGals0, sum_dP_dV = CorrelateGalaxies_LVC(
        prob, distmu, distsigma, distnorm, cat, has3D, MinimumProbCutForCatalogue)
    # alreadysumipixarray,AlreadyObservedPgw = SubstractPointings2D(PointingsFile,prob, nside, FOV,alreadysumipixarray)
    ra, dec, tGals, alreadyObservedPgw, alreadyObservedPgal, alreadysumipixarray = SubstractPointings(
        PointingsFile, tGals0, alreadysumipixarray, sum_dP_dV, FOV, prob, nside)

    print('Total PGW covered: ', sum(alreadyObservedPgw))
    print('Total PGAL covered: ', sum(alreadyObservedPgal))
