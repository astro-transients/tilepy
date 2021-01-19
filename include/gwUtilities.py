#!/usr/bin/env python
# Script name: Read.py
# Version: 1.0 (April 2018)
# Author(s): M. Razzano (massimiliano.razzano@pi.infn.it)
#
import sys,os
from astropy.table import Table
import numpy as np

#General Utilities

def TableImport(inputFileName):
    run,MergerID,RA,Dec,distance,z,theta,ndet,SNR,A90,A50 = np.genfromtxt(inputFileName,usecols=(0, 1, 2, 3,4,5,6,7,8,9,10),skip_header=3, unpack=True,dtype='str')
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    z=z.astype(np.float)
    distance = distance.astype(np.float)
    theta=theta.astype(np.float)
    ndet=ndet.astype(np.float)
    SNR=SNR.astype(np.float)
    A90=A90.astype(np.float)
    A50=A50.astype(np.float)
    OutputTable= Table([run,MergerID,RA,Dec,distance,z,theta,ndet, SNR,A90,A50],names=('run','MergerID','RA','Dec','Distance','redshift','theta','ndet','SNR','A90','A50'))
    return OutputTable

