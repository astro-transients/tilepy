import datetime
from datetime import timezone

import astropy.coordinates as co
import healpy as hp
import ligo.skymap.io.fits as lf
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from tilepy.include.Observatories import CTANorthObservatory, HESSObservatory
from tilepy.include.PointingPlotting import (
    LoadPointings,
    PlotPointings_Pretty,
    PlotPointingsTogether,
)
from tilepy.include.PointingTools import (
    Get90RegionPixReduced,
    NightDarkObservation,
    NightDarkObservationwithGreyTime,
    TransformRADec,
    getdate,
)

from tilepy.include.CampaignDefinition import ObservationParameters


def LocateSource(filename, ra, dec, PercentCov=90):

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]
    npix = len(prob)
    nside = hp.npix2nside(npix)

    pix_ra1, pix_dec1, _ = Get90RegionPixReduced(prob, PercentCov, nside)

    reduced_nside = 512
    coordinates = TransformRADec(ra, dec)
    ra1 = coordinates.ra.rad
    dec1 = coordinates.dec.rad

    radecs = co.SkyCoord(pix_ra1, pix_dec1, frame="icrs", unit=(u.deg, u.deg))
    pix_ra11 = radecs.ra.rad
    pix_dec11 = radecs.dec.rad

    Igrb = hp.ang2pix(reduced_nside, ra1, dec1, lonlat=True)
    Imap = hp.ang2pix(reduced_nside, pix_ra11, pix_dec11, lonlat=True)

    print(np.isin(Igrb, Imap))

    hp.mollview(prob, title="Locate source")
    hp.visufunc.projplot(coordinates.ra, coordinates.dec, "r.", lonlat=True)
    hp.graticule()
    plt.show()


def VisibilityOverview_forZenithCut(
    filename, date=datetime.datetime.now(timezone.utc), zenithcut=60
):
    time = getdate(date)

    print(
        "==========================================================================================="
    )

    print("Loading map from ", filename)
    print("Input time: ", time)
    print("Zenith angle cut for plots: ", zenithcut)

    print(
        "==========================================================================================="
    )

    observatory = HESSObservatory()
    observatory2 = CTANorthObservatory()

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]

    titlePlot = (
        observatory.Name
        + "(red) &"
        + observatory2.Name
        + "(blue) visibility at "
        + str(time)
    )
    # Just do a plotting
    hp.mollview(
        prob, coord="C", title=titlePlot, unit="Probability"
    )  # Celestial=Equatorial

    altcoord = np.empty(4000)
    altcoord.fill(90 - zenithcut)
    azcoord = np.random.rand(4000) * 360

    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=time,
        location=observatory.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "red",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )

    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=time,
        location=observatory2.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "blue",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )

    hp.graticule()
    plt.savefig("Visibility_Overview.png")
    print("-----------------------------------------------------------------")
    print(
        "Observatory 1 is : ",
        observatory.Name,
        "and observatory 2 is: ",
        observatory2.Name,
    )
    print("\n plot created: Visibility_Overview.png")


def Time_DarkTime_GreyTime(
    filename,
    cfgFile,
    date=datetime.datetime.now(timezone.utc),
    zenithcut=60,
    maxNights=1,
    maxDuration=28,
    minDuration=10,
):
    time = getdate(date)

    obspar = ObservationParameters()
    obspar.from_configfile(cfgFile)

    altcoord = np.empty(4000)
    altcoord.fill(90 - zenithcut)
    azcoord = np.random.rand(4000) * 360

    NightDarkRuns = NightDarkObservation(time, obspar)
    NightGreyRuns = NightDarkObservationwithGreyTime(time, obspar)

    print("-----------------------------------------------------------------")
    print("Observatory is: ", obspar.obs_name)
    print("Start of Darkness:", NightDarkRuns[0])
    print(NightDarkRuns)
    print("Last run:", NightDarkRuns[-1])
    print("Total number of runs fulfilling darkness condition:", len(NightDarkRuns))
    print("-----------------------------------------------------------------")

    print("Start of Greyness:", NightGreyRuns[0])
    print(NightGreyRuns)
    print("Last run:", NightGreyRuns[-1])
    print("Total number of runs fulfilling greyness condition:", len(NightGreyRuns))
    print("-----------------------------------------------------------------")

    # InputTime
    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=time,
        location=obspar.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "red",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )

    # NightDarkRuns[0]
    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=NightDarkRuns[0],
        location=obspar.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "orange",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )

    # NightDarkRuns[-1]
    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=NightDarkRuns[-1],
        location=obspar.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "orange",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )

    # NightGreyRuns[0]
    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=NightGreyRuns[0],
        location=obspar.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "blue",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )

    # NightGreyRuns[-1]
    RandomCoord = SkyCoord(
        azcoord,
        altcoord,
        frame="altaz",
        unit=(u.deg, u.deg),
        obstime=NightGreyRuns[-1],
        location=obspar.location,
    )
    RandomCoord_radec = RandomCoord.transform_to("icrs")
    hp.visufunc.projplot(
        RandomCoord_radec.ra,
        RandomCoord_radec.dec,
        "blue",
        lonlat=True,
        coord="C",
        linewidth=0.4,
    )
    hp.graticule()
    plt.savefig("VisibilityTime_DarkTime_GreyTime_Overview" + obspar.obs_name + ".png")

    print("Plot created for observatory : ", obspar.obs_name)
    print("-----------------------------------------------------------------")


def CompareTwoTilings(
    filename, PointingsFile1=False, PointingsFile2=False, FOV=2, plotType="mollweide"
):
    # Format of Pointing Files should be YYYY-MM-DD hh:mm:ss RAarray DECarray P_GWarray P_Galarray Round)

    print(
        "==========================================================================================="
    )
    print("Starting the pointing plotting from the following files\n")

    print("Loading map from ", filename)
    print("Filename 1: ", PointingsFile1)
    print("Filename 2: ", PointingsFile2)

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]

    if (PointingsFile1 == "False") or (PointingsFile2 == "False"):
        print("At least one of the pointings has not been given, try again")
        # Just do a plotting
        hp.mollview(prob, coord="C")
        hp.graticule()
        plt.show()

    else:
        print("Loading pointings")
        df1 = LoadPointings(PointingsFile1)
        df2 = LoadPointings(PointingsFile2)

    if "Pgal" in df2.columns:
        print(
            "Summary of 1st file: sum(PW)=",
            sum(df1["PGW"]),
            "sum(PGAL)=",
            sum(df1["Pgal"]),
            "total pointings",
            len(df1["PGW"]),
        )
        print(
            "Summary of 2st file: sum(PW)=",
            sum(df2["PGW"]),
            "sum(PGAL)=",
            sum(df2["Pgal"]),
            "total pointings",
            len(df2["PGW"]),
        )
        print(
            "==========================================================================================="
        )
    else:
        print(
            "Summary of 1st file: sum(PW)=",
            sum(df1["PGW"]),
            "total pointings",
            len(df1["PGW"]),
        )
        print(
            "Summary of 2st file: sum(PW)=",
            sum(df2["PGW"]),
            "total pointings",
            len(df2["PGW"]),
        )
        print(
            "==========================================================================================="
        )

        Coordinates1 = co.SkyCoord(
            df1["RA(deg)"].tolist(),
            df1["DEC(deg)"].tolist(),
            frame="icrs",
            unit=(u.deg, u.deg),
        )
        Coordinates2 = co.SkyCoord(
            df2["RA(deg)"].tolist(),
            df2["DEC(deg)"].tolist(),
            frame="icrs",
            unit=(u.deg, u.deg),
        )
        PlotPointingsTogether(
            prob,
            Coordinates1,
            Coordinates2,
            FOV,
            plotType,
            doPlot=True,
        )
    plt.show()


def Pretty_Plot(
    filename, name, PointingsFile1, dirName, obspar, gal, center, radius, colors
):
    PlotPointings_Pretty(
        filename, name, PointingsFile1, dirName, obspar, gal, center, radius, colors
    )
