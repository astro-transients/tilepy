import sys
sys.path.append('./include')
from PointingPlotting import *
import argparse


#Example python tools/CompareTwoSchedules.py G297595_LALInference.fits -file1 SuggestedPointings_GWOptimisation.txt  -file2 SuggestedPointings_GWOptimisation1.txt

parser = argparse.ArgumentParser(description='Plot together two pointing schedules')
parser.add_argument('GWname', metavar='GWmap',
                    help='the FITS file with the GW localization, e.g. Bayestar.fits.gz')
parser.add_argument('-file1', metavar='PointingFile1',default="False",
                    help='list with 1st subset of pointings\
                                       (YYYY-MM-DD hh:mm:ss RAarray DECarray P_GWarray P_Galarray Round)\
                                       (default: %(default)s)')
parser.add_argument('-file2', metavar='PointingFile2',default="False",
                    help='list with 2nd subset of pointings\
                                       (YYYY-MM-DD hh:mm:ss RAarray DECarray P_GWarray P_Galarray Round)\
                                       (default: %(default)s)')

args = parser.parse_args()


GWFile=args.GWname
PointingsFile1=args.file1
PointingsFile2=args.file2


print("===========================================================================================")
print("Starting the pointing plotting from the following files\n")

print('Loading map from ', GWFile)
print("Filename 1: ",PointingsFile1)
print("Filename 2: ",PointingsFile2)


prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(GWFile)
npix = len(prob)
nside = hp.npix2nside(npix)

if (PointingsFile1=='False') or (PointingsFile2=='False'):
    print('At least one of the pointings has not been given, try again')
    #Just do a plotting
    hp.mollview(prob,coord='C')
    hp.graticule()
    plt.show()

else:
    print('Loading pointings')

    ObservationTimearray1, Coordinates1, Pgw1, Pgal1 = LoadPointingsGAL(PointingsFile1)
    ObservationTimearray2, Coordinates2, Pgw2, Pgal2 = LoadPointingsGAL(PointingsFile2)


    print('Summary of 1st file: sum(PW)=', sum(Pgw1),'sum(PGAL)=',sum(Pgal1),'total pointings',len(Pgal1))
    print('Summary of 2nd file: sum(PW)=', sum(Pgw2),'sum(PGAL)=',sum(Pgal1),'total pointings',len(Pgal2))
    print("===========================================================================================")

    name1='File1'
    name2='File2'
    FOV=2
    PlotPointingsTogether(prob,ObservationTimearray1,Coordinates1,name1,Coordinates2,name2, nside, FOV, doplot=True)
    plt.show()


