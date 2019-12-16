#Packages
import numpy as np
from optparse import OptionParser
from gwUtilities import *
import sys
from astropy.coordinates import SkyCoord
sys.path.append('/Users/mseglar/Documents/GitHub/GWfollowup')
from GWHESSPointingTools import *
sys.path.append('/Users/mseglar/Documents/GitHub/CTASimulationsGW')
from GWCTAPointingTools import IsSourceInside,TableImportCTA
from GWCTAObservationScheduler import *

#Characteristics
FOV=2.5
nside=1024
# Input
j = sys.argv[1]
j = int(j)

random.seed()

auxObservationTime0 = randomDate("2016-1-1 0:00:00", "2016-12-31 23:59:60", random.random())
ObservationTime0 = datetime.datetime.strptime(auxObservationTime0, '%Y-%m-%d %H:%M:%S')
print(ObservationTime0)
#Get RA Dec
InputFileName='/Users/mseglar/Documents/GitHub/CTASimulationsGW/BNS-GW-new.txt'
InputList = TableImportCTA(InputFileName)
#HESS_SourcesTable=InputList[(InputList['Dec']<30.)&(InputList['Dec']>-60.)]
#HESS_SourcesDec=HESS_SourcesTable['Dec']
#HESS_SourcesRA=HESS_SourcesTable['RA']
Source=SkyCoord(InputList['RA'][j],InputList['Dec'][j], frame='fk5', unit=(u.deg, u.deg))
#outfilename='BNS-GW-HESS.txt'
#print(outfilename)
#ascii.write(HESS_SourcesTable, outfilename,overwrite=True)


#Load Pointing Coordinates
print(InputList['MergerID'][j])
InputGW="/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/skymaps/"+InputList['run'][j]+'_'+InputList['MergerID'][j]+"_skymap.fits.gz"

print('===============  2D algorithm seach ===========')
Pointings,Tmerge=PGWonFoV(InputGW,ObservationTime0)
GW_Pointings=SkyCoord(Pointings['RA(deg)'],Pointings['DEC(deg)'], frame='fk5', unit=(u.deg, u.deg))
Found, Npoiting=IsSourceInside(GW_Pointings,Source,FOV,nside)

#Get time to first observation
#InputList['Npointing']=Npoiting
#InputList['NTotal']=len(Pointings)

print('Time merger',Tmerge)
if(Npoiting==-1):
    InputList['T-To']=-1000
else:
    print('Time Difference',Pointings['Observation Time UTC'][Npoiting],'-',Tmerge)
    print('Time Difference',(Pointings['Observation Time UTC'][Npoiting]-Tmerge).total_seconds())
    InputList['T-To']=(Pointings['Observation Time UTC'][Npoiting]-Tmerge).total_seconds()


if Found:
    outfilename='2DAlgo_FindTheSource/Entry'+str(j)+'.txt'
    #print(outfilename)
    f=open(outfilename,'w')
    f.write(str(np.sum(Pointings['PGW']))+' '+str(Npoiting)+' '+str(len(Pointings))+'\n')
    #ascii.write(InputList[j], outfilename,overwrite=True)

print('===============  3D algorithm seach ===========')
#Galaxy catalog

galFile= "/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/GalaxyCatalogs/" + InputList['run'][j] +'_galaxylist.txt'

name = InputGW.split('.')[0].split('/')[-1]
prepath='/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/'
dirName=prepath+'%s/PGalonFoV' % name
import os
#if not os.path.exists(dirName):
#    os.makedirs(dirName)

#print(dirName)
Pointings,cat,nside=PGalonFoV(InputGW,galFile, InputList['Distance'][j],InputList['DistMax'][j],InputList['DistMin'][j], ObservationTime0,dirName)

GAL_Pointings=SkyCoord(Pointings['RA(deg)'],Pointings['DEC(deg)'], frame='fk5', unit=(u.deg, u.deg))
Found, Npoiting=IsSourceInside(GAL_Pointings,Source,FOV,nside)

if Found:
    outfilename='3DAlgo_FindTheSource/Entry'+str(j)+'.txt'
    #print(outfilename)
    f=open(outfilename,'w')
    f.write(str(np.sum(Pointings['PGW']))+' '+str(Npoiting)+' '+str(len(Pointings))+'\n')
    #ascii.write(InputList[j], outfilename,overwrite=True)
