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

def ParseConfiguration(inputFileName,Separator):
    "Parse a configuration file"
    print('Parsing file',inputFileName)
    datalines=open(inputFileName).readlines()
    OutputDict={}
    for di in range(len(datalines)):
        if (datalines[di][0]!='#'):
            data=datalines[di].strip('\n').split(Separator)         
            OutputDict[data[0]]=data[1]

    #Now, substitute the environment variables, in order to have a consisten enviornment
    Followup_dir=os.environ['FOLLOWUP_DIR']
    for key, val in OutputDict.iteritems():
        OutputDict[key]=val.replace('$datadir',os.path.join(Followup_dir,OutputDict['datadir']))
    for key, val in OutputDict.iteritems():
        OutputDict[key]=val.replace('$bayestardir',os.path.join(Followup_dir,OutputDict['bayestardir']))
    for key, val in OutputDict.iteritems():
        OutputDict[key]=val.replace('$maindir',Followup_dir)

    print('**Parsed variables:')
    for key, val in OutputDict.iteritems():
        print(key,'--',val)
    return OutputDict

def ParseCatalog(catFileName,Separator="\t"):
    "Parse a catalog of sources"
    print('Parsing catalog file',catFileName)
    datalist=open(catFileName).readlines()
    header=(datalist[1][1:-1]).split(Separator)
    print('Read',len(header),'columns:')
    print(header)
    SourceList=[]
    for di in range(3,len(datalist)):
        Source={}
        line=datalist[di].strip('\n').split(Separator)
        if (line[-1]==''):
            line=line[:-1]
        for li in range(len(line)):
            Source[header[li]]=line[li].strip('')
        SourceList.append(Source)    

    print('Read',len(SourceList),'sources')
    return SourceList
