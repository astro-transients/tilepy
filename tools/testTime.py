from GWCTAPointingTools import *
from GWHESSPointingTools import *

random.seed()

auxObservationTime0 = randomDate("2016-1-1 0:00:00", "2016-12-31 23:59:60", random.random())
ObservationTime0 = datetime.datetime.strptime(auxObservationTime0, '%Y-%m-%d %H:%M:%S')
time = NextWindowTools.NextObservationWindow(ObservationTime0,CTANorthObservatory())
WindowDurations=[15,17,20,23,27,33,40,50,64,85]
arraytest=NextWindowTools.CheckWindowCreateArray(time,CTANorthObservatory(),WindowDurations)
print(arraytest)