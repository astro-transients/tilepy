from TilingDetermination import *
from RankingObservationTimes import *
from PointingPlotting import *
import os

random.seed()
outfilename = 'RandomTimes.txt'
f = open(outfilename, "w+")
for i in range(0,2500):
    AuxObservationTime = randomDate("2016-1-1 0:00:00", "2016-12-31 23:59:60", random.random())
    f.write(AuxObservationTime+'\n')
