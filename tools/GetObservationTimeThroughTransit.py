#Packages
from optparse import OptionParser
from gwUtilities import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from astropy.table import Table
import sys
sys.path.append('/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/')
from gwUtilities import *
#sys.path.append('/Users/mseglar/Documents/GitHub/GWfollowup')
from GWHESSPointingTools import *
from GWHESSObservationScheduler import *




def GetObservationPeriod(msource,inputtime,AltitudeCut,nights,plotnumber):

    observatory=EarthLocation(lat=-23.271778*u.deg, lon=16.50022*u.deg, height=1835 *u.m)
    initialframe= AltAz(obstime=inputtime,location=observatory)

    ##############################################################################

    #If the alert arrives during night, the time is reduced to 12 hours so we dont observe the night after
    suninitial= get_sun(inputtime).transform_to(initialframe)

    if(suninitial.alt< -18.*u.deg):
        hoursinDay=12
    else:
        hoursinDay =24
    delta_day = np.linspace(0, hoursinDay+24*(nights-1), 1000*nights)*u.hour
    interval=(hoursinDay+24*(nights-1))/(1000.*nights)
    #print(interval)
    x=np.arange(int(hoursinDay/interval), dtype=int)
    firstN=np.full_like(x,1)
    ratio2=24./interval
    otherN=[]
    for i in range(2,nights+1):
        otherN.extend(np.full_like(np.arange(int(ratio2)),i))
    NightsCounter = []
    NightsCounter.extend(firstN)
    NightsCounter.extend(otherN)
    #print(len(NightsCounter), 'vs', len(delta_day))
    while len(NightsCounter)!=len(delta_day):
        NightsCounter.extend([nights])
    #print('After',len(NightsCounter),'vs',len(delta_day))
    times = inputtime + delta_day
    frame= AltAz(obstime=times,location=observatory)

    ##############################################################################
    #SUN
    sunaltazs_July12_to_13 = get_sun(times).transform_to(frame)
    #MOON
    moon_July12_to_13 = get_moon(times)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame)
    #SOURCE
    #print('msource',moon_July12_to_13)
    #print('moonaltazs_July12_to_13',moonaltazs_July12_to_13)
    msourcealtazs_July12_to_13 = msource.transform_to(frame)
    #print('msource',msource)
    #print('msourcealtazs_July12_to_13',msourcealtazs_July12_to_13)
    ##############################################################################



    Altitudes = Table([times,msourcealtazs_July12_to_13.alt, sunaltazs_July12_to_13.alt, moonaltazs_July12_to_13.alt,NightsCounter],
                           names=['Time UTC', 'Alt Source', 'Alt Sun', 'AltMoon','NightsCounter'])
    Times=Altitudes['Time UTC']
    mask=[(Altitudes['Alt Sun'] < -18.) & (Altitudes['Alt Source']> AltitudeCut) & (Altitudes['AltMoon'] < -0.5)]
    ScheduledTimes=Times[mask]
    print(ScheduledTimes[0],ScheduledTimes[-1])
    #outfilename='ObservationTable_Check_%g.txt'%plotnumber
    #ascii.write(Altitudes[mask], outfilename,overwrite=True)
    Nights=Altitudes[mask]['NightsCounter']
    nightcounter=0
    for i in range(1,nights+1):
        if i in Nights:
            nightcounter=nightcounter+1
    #Night when the source first observed
    FirstObsNight = 0
    for j in range(1,nights+1):
        if j in Nights:
            FirstObsNight = j
            break

    doplot=True
    if doplot:
        plt.figure(figsize=(200, 12))
        plt.plot(delta_day, sunaltazs_July12_to_13.alt, color='r', label='Sun')
        plt.plot(delta_day, moonaltazs_July12_to_13.alt, color=[0.75]*3, ls='--', label='Moon')
        plt.scatter(delta_day, msourcealtazs_July12_to_13.alt,c=msourcealtazs_July12_to_13.az, label='Point source', lw=0, s=8,cmap='viridis')
        plt.fill_between(delta_day.to('hr').value, 0, 90,sunaltazs_July12_to_13.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_day.to('hr').value, 0, 90,(sunaltazs_July12_to_13.alt < -18*u.deg)&(moonaltazs_July12_to_13.alt < -0.5*u.deg), color='k', zorder=0)
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        plt.xlim(0,hoursinDay+24*(nights-1))
        #plt.xticks(np.arange(13)*2 -12)
        plt.ylim(0, 90)
        plt.xlabel('Hours after injection')
        plt.ylabel('Altitude [deg]')
        plt.savefig('./Plots/Source%g.png'% plotnumber)
    if len(ScheduledTimes)==0:
        #print('No observaiton scheduled')
        return 0,nightcounter,FirstObsNight
    elif (len(ScheduledTimes)==1)&(MaxNights!=1):
        hours = (delta_day[1] - delta_day[0]) * len(ScheduledTimes)
        return float(str(hours).split(' ')[0]),nightcounter,FirstObsNight
    else:
        hours=(delta_day[1]-delta_day[0])*len(ScheduledTimes)
        #print(hours)
        #print(str(hours))
        if MaxNights==1:
            nights=1
            return float(str(hours).split(' ')[0]),nightcounter,FirstObsNight
        else:
            nights=(ScheduledTimes[-1].tt.datetime-ScheduledTimes[0].tt.datetime)
            #print(ScheduledTimes)
            #print(nights)
            try:
                retnights=float(str(nights).split(' ')[0])
                return float(str(hours).split(' ')[0]), nightcounter,FirstObsNight
            except ValueError:
                return float(str(hours).split(' ')[0]),nightcounter,FirstObsNight


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


