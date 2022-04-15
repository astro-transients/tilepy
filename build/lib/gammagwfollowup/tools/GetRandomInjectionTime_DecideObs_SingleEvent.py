import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import healpy as hp
import random
from gammagwfollowup.include.PointingTools import randomDate, TableImportCTA
from astropy.time import Time
from astropy.table import Table
import astropy.coordinates as co

import numpy as np
from astropy.io import ascii
from astropy import units as u



if __name__ == "__main__":

    InputFileName = './../dataset/BNS-GW-new.txt'
    InputList = TableImportCTA(InputFileName)


    #Get Random time through the year
    random.seed()

    trial = []
    runA = []
    IDA= []
    ObsTime = []
    Observatory = []

    # Simple, first order approximations from the hypothesis that it should be observe with the most optimistic Zenith angle
    # Next step: who may be able to observe first? (instead of under a better angle)
    # The questions if there is OK visibility/ NIGHT?  comes with the scheduling
    #Table with: injection times, duration of the window , time of observation of the source


    # Obtain 100 different random times for a single event

    maxRun = 20
    i = 10247

    for j in range(0,maxRun):
        print('------------------------------------------')
        print('Iteration number',j,'out of',maxRun)
        AuxObservationTime = randomDate("2016-1-1 0:00:00", "2016-12-31 23:59:60", random.random())
        print(AuxObservationTime)

        ObservationTime0 = Time(AuxObservationTime)
        ObsTime.append(AuxObservationTime)
        # Maximum probability pixel
        # Choose who is following (not if it is successful)
        # CUT: dec cut

        InputGW = "/Users/mseglar/Documents/GitLab/dataset/skymaps/" + InputList['run'][i] + '_' + \
                  InputList['MergerID'][i] + "_skymap.fits.gz"

        print(InputGW)

        trial.append(j)
        runA.append(InputList['run'][j])
        IDA.append(InputList['MergerID'][j])

        # Look the maximum probability pixel
        tprob = hp.read_map(InputGW, field=range(1))
        prob = hp.pixelfunc.ud_grade(tprob, 512, power=-2)
        npix = len(prob)
        nside = hp.npix2nside(npix)
        ipix = np.argmax(prob)
        tt, pp = hp.pix2ang(nside, ipix)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)
        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        print(skycoord)
        if(dec2>0):
            Observatory.append('North')
            print("The random time is: ", ObservationTime0,'Observatory is North')

        else:
            Observatory.append('South')
            print("The random time is: ", ObservationTime0,'Observatory is South')

    ObservationTable = Table([trial, runA, IDA, ObsTime, Observatory],names=['Trial','run','MergerID','Time','Observatory'])
    outfilename= 'ObservationTable_' + InputList['run'][i] + '_' +InputList['MergerID'][i]+ '_fromGWHotspot.txt'
    ascii.write(ObservationTable, outfilename,overwrite=True)


