import sys
sys.path.append('./include')
from PointingPlotting import *
import argparse
from datetime import timezone
def getdate(x):
    if isinstance(x, datetime.datetime):
        return x
    elif isinstance(x, str):
        return datetime.datetime.strptime(x,'%Y-%m-%d %H:%M:%S')
    else:
        print("ERROR: something is wrong with the format of the date: ", x)
        return None

from astropy.utils.iers import conf as iers_conf
iers_conf.iers_auto_url
iers_conf.iers_auto_url = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
iers_conf.iers_auto_url_mirror = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'

parser = argparse.ArgumentParser(description='Plot together two pointing schedules')
parser.add_argument('GWname', metavar='GWmap',
                    help='the FITS file with the GW localization, e.g. Bayestar.fits.gz')

parser.add_argument('-date', metavar='\"YYYY-MM-DD HH:MM:SS\"', default=datetime.datetime.now(timezone.utc),
                    help='optional: date and time of the event (default: NOW, i.e. %(default)s)')

parser.add_argument('-zenith', metavar='ZenithCut',default='60',
                    help='Zenith angle cut of the observatory FoV. Usually set to 45 or to 60 degs')

args = parser.parse_args()

GWFile=args.GWname
time=getdate(args.date)
zenithcut=int(args.zenith)


#CTASouthObservatory()

print("===========================================================================================")

print('Loading map from ', GWFile)
print("Input time: ",time)

#Get darkness period
maxNights=1
maxDuration=28
minDuration=10
observatory=HESSObservatory()
NightDarkRuns =NightDarkObservation(time,observatory,maxNights,maxDuration,minDuration)
NightGreyRuns=NightDarkObservationwithGreyTime(time,observatory,maxNights,maxDuration,minDuration)
print("-----------------------------------------------------------------")
print('Start of Darkness:',NightDarkRuns[0])
print(NightDarkRuns)
print('Last run:',NightDarkRuns[-1])
print('Total number of runs fulfilling darkness condition:',len(NightDarkRuns))
print("-----------------------------------------------------------------")
#NightGreyRuns = NightDarkObservationwithGreyTime(time,observatory,maxNights,maxDuration,minDuration)
print('Start of Greyness:',NightGreyRuns[0])
print(NightGreyRuns)
print('Last run:',NightGreyRuns[-1])
print('Total number of runs fulfilling greyness condition:',len(NightGreyRuns))
print("-----------------------------------------------------------------")

prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(GWFile)
npix = len(prob)
nside = hp.npix2nside(npix)

TITLE_HESS = 'H.E.S.S.(red) & LST(blue) visbility at '+str(time)
#Just do a plotting
hp.mollview(prob,coord='C', title = TITLE_HESS) #Celestial=Equatorial
#frame = co.AltAz(obstime=time, location=observatory.Location)

#inputTime
altcoord = np.empty(4000)
altcoord.fill(90 - zenithcut)
azcoord = np.random.rand(4000) * 360
observatory=HESSObservatory()
RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory.Location)
RandomCoord_radec = RandomCoord.transform_to('fk5')
hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'red', lonlat=True, coord="C", linewidth=0.4)
'''
NightGreyRuns[0]
RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightDarkRuns[0],location=observatory.Location)
RandomCoord_radec = RandomCoord.transform_to('fk5')
hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'orange', lonlat=True, coord="C", linewidth=0.4)

NightGreyRuns[-1]
RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightDarkRuns[-1],location=observatory.Location)
RandomCoord_radec = RandomCoord.transform_to('fk5')
hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'orange', lonlat=True, coord="C", linewidth=0.4)
'''
altcoord = np.empty(4000)
altcoord.fill(90 - 75)
azcoord = np.random.rand(4000) * 360
observatory=CTANorthObservatory()
#GreyTime[0]
RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory.Location)
RandomCoord_radec = RandomCoord.transform_to('fk5')
hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'blue', lonlat=True, coord="C", linewidth=0.4)

#GreyTime[0]
#RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightGreyRuns[-1],location=observatory.Location)
#RandomCoord_radec = RandomCoord.transform_to('fk5')
#hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'blue', lonlat=True, coord="C", linewidth=0.4)

hp.graticule()
plt.savefig("Visibility_Overview.png")
print("\n plot created: Visibility_Overview.png")
#plt.show()

