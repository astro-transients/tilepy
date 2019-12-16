import sys
import argparse
sys.path.append('../../GW-Followup/include')
from PointingPlotting import *
from RankingObservationTimes import *

sys.path.append('./include')
from GWCTAObservationScheduler import *
from AnalysisTools import *
from SimulationTools import *


start = time.time()

###########################
#####    Parameters  ######
###########################

parser = argparse.ArgumentParser(description='Start the CTA GW pointing simulation and analysis (PgwinFov)')
parser.add_argument('index', metavar='GWmap',
                    help='the FITS file with the GW localization, e.g. Bayestar.fits.gz')
parser.add_argument('-params', metavar='parameters.ini', default="./configs/PGWinFoV.ini",
                    help='file holding the main parameters (default: %(default)s)')
args = parser.parse_args()

j=int(args.index)
parameters=args.params



###########################
#####      Files     ######
###########################

InputFileName = '../dataset/BNS-GW_onAxis5deg.txt'
InputList = TableImportCTA(InputFileName)
Source = SkyCoord(InputList['RA'][j], InputList['Dec'][j], frame='fk5', unit=(u.deg, u.deg))
#print(Source)
#Source = SkyCoord(345,-27, frame='fk5', unit=(u.deg, u.deg))
#t = 0.5 * np.pi - HESS_Sources.dec.rad
#p = HESS_Sources.ra.rad
#ipixSource = hp.ang2pix(512, t, p)


# GW file
GWFile = "../dataset/skymaps/" + InputList['run'][j] + '_' + InputList['MergerID'][j] + "_skymap.fits.gz"
name =  InputList['run'][j] + '_' + InputList['MergerID'][j]


#galFile= "../dataset/GalaxyCatalogs/" + InputList['run'][j] +'_galaxylist.txt'
galFile= "../dataset/GalaxyCatalogs/run0002_galaxylist.txt"

InjectTimeFile = '../dataset/BNS-GW-Time_onAxis5deg.txt'
InputTimeList = TableImportCTA_Time(InjectTimeFile)

#ObservationTime0 = datetime.datetime(2017, 7, 17, 20, 15, 13)
ObservationTime0 = datetime.datetime.strptime(InputTimeList['Time'][j], '%Y-%m-%d %H:%M:%S.%f')

InputObservationFile='../dataset/Selected_ObsTimes/'+InputList['run'][j] + '_' + InputList['MergerID'][j].split('r')[-1] +'_obstime.txt'



import os
#if not os.path.exists('./output/'+name):
#    os.makedirs('./output/'+name)
#prepath='/Users/mseglar/Documents/GitHub/CTASimulationsGW/'
dirName='./output/%s/PGWonFoV' % name
if not os.path.exists(dirName):
    os.makedirs(dirName)
print(' ')
print("===========================================================================================")
print("Starting the CTA pointing calculation with the following parameters\n")
print('GW file with merger ID ',InputList['MergerID'][j],'and run ',InputList['run'][j])
print('Input times: ',InjectTimeFile)
print("Date: ",InputTimeList['Time'][j])
print("Galaxy Catalog: ",galFile)
print("Parameters: ",parameters)
print("Output: ",dirName)
print("===========================================================================================")
print(' ')
try:
    InputObservationList = TableImportCTA_Obs(InputObservationFile)
except:
    filepath='../dataset/GammaCatalogV1.0/'+str(InputList['run'][j]) + '_' + str(InputList['MergerID'][j].split('r')[-1]) + ".fits"
    fitsfile = fits.open(filepath)
    luminosity = fitsfile[0].header['EISO']
    outfilename = dirName+'/SimuGW' + str("{:03d}".format(j)) + '.txt'
    f = open(outfilename, 'w')
    f.write(str(InputList['run'][j])+ ' ' +InputList['MergerID'][j].split('r')[-1]+ ' ' + str(InputList['Distance'][j])+ ' '+str(InputList['theta'][j])+ ' ' +str(InputList['A90'][j])+ ' ' +str(luminosity) + ' ' + str(0) + ' ' + str(0) + ' ' + str(-1) +' ' + str(0) +' '+str(-1)+' '+str(-1)+'\n')



#InputObservationList = TableImportCTA_Obs(InputObservationFile)
# ObservationTime0 = datetime.datetime.now(pytz.timezone("UTC"))
# ObservationTime0 = ObservationTime0.replace(tzinfo=None)


# First option: Obtain the Pointing by running the gravitational wave follow-up algorithm
SuggestedPointings,obs0,FOV,nside=PGWonFoV(GWFile,InputObservationList,InputTimeList['Observatory'][j], ObservationTime0,parameters,dirName)
end = time.time()
print('Execution time: ', end - start)

pointingsFile = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
ascii.write(SuggestedPointings, pointingsFile, overwrite=True, fast_writer=False)

'''
# Second option (mainly for testing the gammapy simulation+analysis): Obtain the Pointing from a table
nside=256
FOV=3.5 #High value in order to cover the source
pointingsFile = '%s/SuggestedPointings_GWOptimisation.txt' % dirName

RAarray, DECarray, P_GWarray= np.genfromtxt(pointingsFile, usecols=(2, 3,4), skip_header=1, unpack=True)  # ra, dec in degrees
preDefWind= np.genfromtxt(pointingsFile, usecols=(5), skip_header=1, unpack=True,dtype=int)  # ra, dec in degrees

SuggestedPointings = Table([RAarray, DECarray, P_GWarray,preDefWind], names=('RA(deg)','DEC(deg)','PGW','preDefWind'))
'''


print(SuggestedPointings)

print("===========================================================================================")
print()


if(len(SuggestedPointings)!=0):
    #Is the source inside?

    Pointings = SkyCoord(SuggestedPointings['RA(deg)'], SuggestedPointings['DEC(deg)'], frame='fk5',
                         unit=(u.deg, u.deg))
    # Add the predefined windows to the table
    predefWind = SuggestedPointings['preDefWind']
    totalPGW = sum(SuggestedPointings['PGW'])

    #Check if the source is covered by the scheduled pointings
    Found,nP=IsSourceInside(Pointings,Source,FOV,nside)
    #ProduceSummaryFile(InputList['run'][j],InputList['MergerID'][j], InputObservationList, predefWind,nP,j)

    if Found==False:
        nP = -1
        ProduceSummaryFile(InputList, InputObservationList,predefWind, nP, j,'GW',totalPGW,dirName)
        print('Source not covered')
    else:
        ProduceSummaryFile(InputList, InputObservationList,predefWind, nP, j,'GW',totalPGW,dirName)
        print('Found in scheduled observation:', nP)
        if type(nP) != int:
            nP = nP[0]
        print('Observation corresponds to the predefined window number:', np.int(predefWind[nP]))
        print(InputObservationList[nP])
        #outfilename = '%s/SuggestedPointings_GALOptimisation.txt' % dirName
        #ascii.write(SuggestedPointings, outfilename, overwrite=True,fast_writer=False)
        print()
        print("===========================================================================================")
        print("Starting the simulation of the CTA observations\n")
        IRFfile="../dataset/IRF-prod3b-v2/bcf/"+InputTimeList['Observatory'][j]+'_z'+str(InputTimeList['MeanZen'][j])+'_0.5h/irf_file.fits'
        print('Using IRF file:',IRFfile)
        gcat=LoadGalaxiesSimulation(galFile)

        print('Plotting the observations')
        PointingPlotting(GWFile,gcat,name,dirName,pointingsFile)


        ##RankingTimes(ObservationTime0, GWFile, gcat,parameters, dirName, pointingsFile)

        # Simulation of a single observation with single spectrum (not really used)
        #CubeSimulationSingleObservation_singleSpectrum(InputList['run'][j],InputList['MergerID'][j],Pointings,nP,Source,IRFfile,InputObservationList[np.int(predefWind[nP])],dirName)

        # Simulation of a single observation with LC evolution
        CubeSimulationSingleObservation_dynamicSpectrum(InputList['run'][j],InputList['MergerID'][j],Pointings,nP,Source,IRFfile,InputObservationList[np.int(predefWind[nP])],dirName)

        # Simulation of a several observations with LC evolution. Slow in case of lot of short observations..
        #CubeSimulationMultiObservation_dynamicSpectrum(InputList['run'][j],InputList['MergerID'][j],Pointings,nP,Source,IRFfile,InputObservationList[np.int(predefWind[nP])],dirName)

        # Analysis of the simulations in 3D
        LikelihoodFit_Analysis_3DCube(dirName)

        # Analysis of the simulations in 4D
        ##LikelihoodFit_Analysis_4DCube(dirName)
else:
    print('No observations are scheduled')
    filepath='../dataset/GammaCatalogV1.0/'+str(InputList['run'][j]) + '_' + str(InputList['MergerID'][j].split('r')[-1]) + ".fits"
    fitsfile = fits.open(filepath)
    luminosity = fitsfile[0].header['EISO']
    outfilename = dirName+'/SimuGW' + str("{:03d}".format(j)) + '.txt'
    f = open(outfilename, 'w')
    f.write(
        'RunList' + ' ' + 'MergerID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Luminosity' + ' ' + 'TotalInputObs' + ' ' + 'Obs' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'WindowinInputList' + ' ' + 'TotalProb' + '\n')
    f.write(str(InputList['run'][j])+ ' ' +InputList['MergerID'][j].split('r')[-1]+ ' ' + str(InputList['Distance'][j])+ ' '+str(InputList['theta'][j])+ ' ' +str(InputList['A90'][j])+ ' ' +str(luminosity) + ' ' + str(len(InputObservationList)) + ' ' + str(0) + ' ' + str(-1) +' ' + str(0) +' '+str(-1)+' '+str(-1)+'\n')


