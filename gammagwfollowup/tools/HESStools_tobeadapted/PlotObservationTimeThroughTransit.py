import sys
sys.path.append('./../include')
from RankingObservationTimes import *

### What is this script doing?
# input ?
# output?


if __name__ == "__main__":

    AltitudeCut=30.
    MaxNights=1

    #Get RA Dec
    InputFileName='/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/BNS-GW-HESS.txt'
    InputList=TableImport(InputFileName)
    HESS_Sources=SkyCoord(InputList['RA'],InputList['Dec'], frame='fk5', unit=(u.deg, u.deg))

    #Get Random time through the year
    random.seed()
    #ObservationTime0 = datetime.datetime.strptime(AuxObservationTime, '%Y-%m-%d %H:%M:%S')

    Vtimes=[]
    Vschhours=[]
    Vschdays=[]
    VFirstObsNight = []
    #Table with: injection times, duration of the window , time of observation of the source
    for i in range(0,2):
        AuxObservationTime = randomDate("2016-1-1 0:00:00", "2016-12-31 23:59:60", random.random())
        print(AuxObservationTime)
        ObservationTime0 = Time(AuxObservationTime)
        print(ObservationTime0)
        schhours,schdays,FirstObsNight=GetObservationPeriod(HESS_Sources[i], ObservationTime0, AltitudeCut,MaxNights,i)
        Vtimes.append(AuxObservationTime)
        Vschhours.append(schhours)
        Vschdays.append(schdays)
        VFirstObsNight.append(FirstObsNight)
        print('Iteration number',i)


    ObservationTable = Table([Vtimes,Vschhours,Vschdays,VFirstObsNight],names=['Time UTC', 'Hours observed', 'Nights observed','First observation night'])
    outfilename='ObservationTable_%g_%i.txt'% (AltitudeCut,MaxNights)
    ascii.write(ObservationTable, outfilename,overwrite=True)


