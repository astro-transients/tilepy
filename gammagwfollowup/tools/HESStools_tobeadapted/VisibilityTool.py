import sys
sys.path.append('./../include')
from RankingObservationTimes import *

if __name__ == "__main__":

    AltitudeCut=10.
    MaxNights=1

    #Get RA Dec
    Source=SkyCoord(270.3846454,14.0727500, frame='fk5', unit=(u.deg, u.deg))
    #SkyCoord('', frame='hmsdms', unit=(u.deg, u.deg))
    Source=SkyCoord('11 04 19 +38 11 41', unit=(u.hourangle, u.deg))
    #print('In sexagesimal:',Source.to_string('hmsdms'))
    #Get Random time through the year
    #random.seed()
    #ObservationTime0 = datetime.datetime.strptime(AuxObservationTime, '%Y-%m-%d %H:%M:%S')

    AuxObservationTime=datetime.datetime.now()
    print(AuxObservationTime)
    ObservationTime0 = Time(AuxObservationTime)
    print(ObservationTime0)
    GetObservationPeriod(ObservationTime0,Source, AltitudeCut,MaxNights,1,'.',True)



