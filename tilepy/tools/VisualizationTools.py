import os
import sys
from ..include.PointingTools import NightDarkObservation, NightDarkObservationwithGreyTime, ObservationParameters, getdate, Get90RegionPixReduced, TransformRADec
from ..include.Observatories import CTANorthObservatory, CTASouthObservatory, HESSObservatory
from ..include.PointingPlotting import LoadPointingsGAL, PlotPointingsTogether, PlotPointings_Pretty

import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
import datetime
from datetime import timezone
import matplotlib.pyplot as plt
from astropy import units as u
import astropy.coordinates as co
from astropy.utils import iers
import ligo.skymap.io.fits as lf

# from astropy.utils.iers import conf as iers_conf
# iers_conf.iers_auto_url = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
# iers_conf.iers_auto_url_mirror = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
iers_file = os.path.join(os.path.abspath(
    os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)


def LocateSource(filename, ra, dec, PercentCov=90):

    # prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]
    npix = len(prob)
    nside = hp.npix2nside(npix)

    pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob, PercentCov, nside)

    reduced_nside = 512
    coordinates = TransformRADec(ra, dec)
    ra1 = coordinates.ra.rad
    dec1 = coordinates.dec.rad

    radecs = co.SkyCoord(pix_ra1, pix_dec1, frame='fk5', unit=(u.deg, u.deg))
    pix_ra11 = radecs.ra.rad
    pix_dec11 = radecs.dec.rad

    Igrb = hp.ang2pix(reduced_nside, ra1, dec1, lonlat=True)
    Imap = hp.ang2pix(reduced_nside, pix_ra11, pix_dec11, lonlat=True)

    print(np.isin(Igrb, Imap))

    hp.mollview(prob, title="Locate source")
    hp.visufunc.projplot(coordinates.ra, coordinates.dec, 'r.', lonlat=True)
    hp.graticule()
    plt.show()


def VisibilityOverview_forZenithCut(filename, date=datetime.datetime.now(timezone.utc), zenithcut=60):
    time = getdate(date)

    print("===========================================================================================")

    print('Loading map from ', filename)
    print("Input time: ", time)
    print("Zenith angle cut for plots: ", zenithcut)

    print("===========================================================================================")

    #################################  Other parameters ##############################

    observatory = HESSObservatory()
    observatory2 = CTANorthObservatory()

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]

    # prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)

    titlePlot = observatory.Name + \
        '(red) &' + observatory2.Name + '(blue) visibility at ' + str(time)
    # hp.newvisufunc.projview(prob, unit="Probability", graticule=True, graticule_labels=True, projection_type="mollweide",
    #    xlabel="RA", ylabel="Dec", phi_convention="clockwise", longitude_grid_spacing=30,title=titlePlot,format="%#.3g")
    # Just do a plotting
    hp.mollview(prob, coord='C', title=titlePlot,
                unit="Probability")  # Celestial=Equatorial
    # frame = co.AltAz(obstime=time, location=observatory.Location)

    altcoord = np.empty(4000)
    altcoord.fill(90 - zenithcut)
    azcoord = np.random.rand(4000) * 360

    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                           location=observatory.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'red', lonlat=True, coord="C", linewidth=0.4)

    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                           location=observatory2.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'blue', lonlat=True, coord="C", linewidth=0.4)

    hp.graticule()
    plt.savefig("Visibility_Overview.png")
    print("-----------------------------------------------------------------")
    print('Observatory 1 is : ', observatory.Name,
          "and observatory 2 is: ", observatory2.Name)
    print("\n plot created: Visibility_Overview.png")


def Time_DarkTime_GreyTime(filename, cfgFile, date=datetime.datetime.now(timezone.utc), zenithcut=60, maxNights=1, maxDuration=28,
                           minDuration=10):
    time = getdate(date)

    ######################## Read parameters from config ##############################
    obspar = ObservationParameters()
    obspar.from_configfile(cfgFile)
    ######################### Other parameters ##############################

    observatory = HESSObservatory()

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]
    npix = len(prob)
    nside = hp.npix2nside(npix)

    altcoord = np.empty(4000)
    altcoord.fill(90 - zenithcut)
    azcoord = np.random.rand(4000) * 360

    NightDarkRuns = NightDarkObservation(time, obspar)
    NightGreyRuns = NightDarkObservationwithGreyTime(time, obspar)

    print("-----------------------------------------------------------------")
    print('Observatory is: ', observatory.Name)
    print('Start of Darkness:', NightDarkRuns[0])
    print(NightDarkRuns)
    print('Last run:', NightDarkRuns[-1])
    print('Total number of runs fulfilling darkness condition:', len(NightDarkRuns))
    print("-----------------------------------------------------------------")

    print('Start of Greyness:', NightGreyRuns[0])
    print(NightGreyRuns)
    print('Last run:', NightGreyRuns[-1])
    print('Total number of runs fulfilling greyness condition:', len(NightGreyRuns))
    print("-----------------------------------------------------------------")

    # InputTime
    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                           location=observatory.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'red', lonlat=True, coord="C", linewidth=0.4)

    # NightDarkRuns[0]
    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightDarkRuns[0],
                           location=observatory.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'orange', lonlat=True, coord="C", linewidth=0.4)

    # NightDarkRuns[-1]
    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightDarkRuns[-1],
                           location=observatory.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'orange', lonlat=True, coord="C", linewidth=0.4)

    # NightGreyRuns[0]
    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightGreyRuns[0],
                           location=observatory.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'blue', lonlat=True, coord="C", linewidth=0.4)

    # NightGreyRuns[-1]
    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=NightGreyRuns[-1],
                           location=observatory.Location)
    RandomCoord_radec = RandomCoord.transform_to('fk5')
    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                         'blue', lonlat=True, coord="C", linewidth=0.4)
    hp.graticule()
    plt.savefig("VisibilityTime_DarkTime_GreyTime_Overview" +
                observatory.Name + '.png')

    print('Plot created for observatory : ', observatory.Name)
    print("-----------------------------------------------------------------")


def CompareTwoTilings(filename, PointingsFile1=False, PointingsFile2=False, FOV=2):
    # Format of Pointing Files should be YYYY-MM-DD hh:mm:ss RAarray DECarray P_GWarray P_Galarray Round)

    print("===========================================================================================")
    print("Starting the pointing plotting from the following files\n")

    print('Loading map from ', filename)
    print("Filename 1: ", PointingsFile1)
    print("Filename 2: ", PointingsFile2)

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]
    npix = len(prob)
    nside = hp.npix2nside(npix)

    if (PointingsFile1 == 'False') or (PointingsFile2 == 'False'):
        print('At least one of the pointings has not been given, try again')
        # Just do a plotting
        hp.mollview(prob, coord='C')
        hp.graticule()
        plt.show()

    else:
        print('Loading pointings')

        ObservationTimearray1, Coordinates1, Pgw1, Pgal1 = LoadPointingsGAL(
            PointingsFile1)
        ObservationTimearray2, Coordinates2, Pgw2, Pgal2 = LoadPointingsGAL(
            PointingsFile2)

        print('Summary of 1st file: sum(PW)=', sum(Pgw1),
              'sum(PGAL)=', sum(Pgal1), 'total pointings', len(Pgal1))
        print('Summary of 2nd file: sum(PW)=', sum(Pgw2),
              'sum(PGAL)=', sum(Pgal1), 'total pointings', len(Pgal2))
        print("===========================================================================================")

        name1 = 'File1'
        name2 = 'File2'

        PlotPointingsTogether(prob, ObservationTimearray1, Coordinates1, name1, Coordinates2, name2, nside, FOV,
                              doplot=True)
    plt.show()


def Pretty_Plot(filename, name, PointingsFile1, dirName):
    PlotPointings_Pretty(filename, name, PointingsFile1, dirName)
