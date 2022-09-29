from TilingDetermination import *
from RankingObservationTimes import *
from PointingPlotting import *
import os


# filename=sys.argv[1]
# name=filename.split('.')[0]

# filename=sys.argv[1]
# name=filename.split('.')[0]

# Input
j = sys.argv[1]
j = int(j)

# Round
r = sys.argv[2]
r = int(r)

# GW filename

GWfilename = '/Users/mseglar/Documents/CurrentPhD/HESS/GW/going-the-distance-o2-skymaps/bayestar.txt'
name= np.genfromtxt(GWfilename, usecols=(0), unpack=True, dtype='str')
filename='/Users/mseglar/Documents/CurrentPhD/HESS/GW/going-the-distance-o2-skymaps/bayestar/'+name[j]
name = filename.split('.')[0].split('/')[-1]

# Select Injection Time
#ObservationTime0 = datetime.datetime(2017, 8, 17, 10, 30, 13)
# ObservationTime0 = datetime.datetime.now(pytz.timezone("UTC"))
# ObservationTime0 = ObservationTime0.replace(tzinfo=None)

Timefilename = 'RandomTimes.txt'
AuxObservationTime1,AuxObservationTime2 =np.genfromtxt(Timefilename, usecols=(0,1), unpack=True, dtype='str')
AuxObservationTime=AuxObservationTime1[j+250*r]+' '+AuxObservationTime2[j+250*r]
print(AuxObservationTime)
ObservationTime0 = datetime.datetime.strptime(AuxObservationTime, '%Y-%m-%d %H:%M:%S')

start = time.time()
SuggestedPointings, cat= PGalonFoV(filename, ObservationTime0)
end = time.time()


print("===========================================================================================")
print()

print('Execution time: ', end - start)



#Include the other parameters in the .txt
prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)


print(SuggestedPointings)
nPoint=len(SuggestedPointings)
P_GAL=np.sum(SuggestedPointings['Pgal'])
P_GW=np.sum(SuggestedPointings['PGW'])


outfilename = './ResultsPGalOnFoV/GWInjection' + str(j) +'Round_'+ str(r) +'.txt'
f = open(outfilename, "w+")
concat = str(nPoint) + ' ' +str(end - start)+' '+ str(P_GW) + ' ' + str(P_GAL) + ' ' + str(detectors) + ' '+ str(thisDistance) +' ' + str(thisDistanceErr) +'\n'
f.write(concat)