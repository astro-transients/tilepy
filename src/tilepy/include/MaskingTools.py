#
# Copyright (C) 2016-2024  tilepy developers (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
##################################################################################################
#                           Low-level tools to mask skymaps                                     #
##################################################################################################


import datetime

#####################################################################
# Packages
import os

import astropy.coordinates as co
import ephem
import healpy as hp
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pytz
import six
import tables
from astropy import units as u
from astropy.coordinates import AltAz, get_body
from astropy.table import Table
from astropy.time import Time
from matplotlib.path import Path
from pytz import timezone
from six.moves import configparser


if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser

__all__ = [
    "ZenithAngleCut",
    "VisibleAtTime",
    "FulfillsRequirement",
    "FulfillsRequirementGreyObservations",
    "GetEarthOccultedPix",
    "GetMoonOccultedPix",
    "GetSunOccultedPix",
    "OccultationCut",
    "SAA_Times",
    "GetBestGridPos2D",
    "GetBestGridPos3D",
]


def ZenithAngleCut(prob, time, obspar):
    """
    Mask in the pixels with zenith angle larger than the max Zenith cut
    """
    nside = obspar.reducedNside
    minProbcut = obspar.minProbcut
    maxZenith = obspar.maxZenith
    observatory = obspar.location
    minMoonSourceSeparation = obspar.minMoonSourceSeparation
    useGreytime = obspar.useGreytime

    frame = co.AltAz(obstime=time, location=observatory)
    pprob = prob

    mzenith = hp.ma(pprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame="fk5", unit=(u.deg, u.deg))
    altaz_map = targetCoord_map.transform_to(frame)
    maskzenith[altaz_map.alt.value < 90 - maxZenith] = 1
    mzenith.mask = maskzenith

    yprob = ma.masked_array(pprob, mzenith.mask)

    if np.sum(yprob) < minProbcut:
        ObsBool = False
    else:
        ObsBool = True

    if useGreytime and ObsBool:
        moonaltazs = get_body("moon", Time(time, scale="utc")).transform_to(
            AltAz(obstime=Time(time, scale="utc"), location=observatory)
        )
        separations = altaz_map.separation(moonaltazs)
        mask_moonDistance = np.zeros(hp.nside2npix(nside), dtype=bool)
        mask_moonDistance[separations < minMoonSourceSeparation * u.deg] = 1
        mzenith = hp.ma(pprob)
        mzenith.mask = mask_moonDistance
        yprob = ma.masked_array(pprob, mzenith.mask)

        if np.sum(yprob) < minProbcut:
            ObsBool = False
        else:
            ObsBool = True

    return ObsBool, yprob


def VisibleAtTime(test_time, galaxies, obspar):
    """
    Determine if prompt or afterglow follow-up is possible by checking if there are galaxies
    with non-negligible probability of hosting the NSM in the FoV.

    Process:
        1. Check if any galaxy is visible; if not, follow-up is 'AFTERGLOW'.
        2. Loop over zenith angle and select subsets of galaxies.
        3. Stop if maximum p-value of this subset is smaller than 75% of the previous subset.
        4. Otherwise, apply a stricter zenith cut and repeat.
        5. Select the galaxy with highest p-value fulfilling both criteria as the target.

    Returns
    -------
    is_vis : bool
        True if a galaxy is visible now, False otherwise.
    alt_az : numpy.ndarray
        Altitude and azimuth locations of galaxies.

    """

    # print()
    # print("Check visibility at time {0}".format(test_time))

    # observatory time and location to look up visibility of objects

    # observatory = co.EarthLocation(lat=-23.271333 * u.deg,lon=16.5 * u.deg, height=1800 * u.m)
    maxz = obspar.maxZenith
    observatory = obspar.location
    frame = co.AltAz(obstime=test_time, location=observatory)
    # print('galaxies',galaxies)
    # print('galaxies',len(galaxies['RAJ2000']))
    radecs = co.SkyCoord(
        galaxies["RAJ2000"], galaxies["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )
    if len(radecs) > 0:
        thisaltaz = radecs.transform_to(frame)

        # add altitude to topGals array
        # already sorted by descending probability value

        galaxies["Alt"] = thisaltaz.alt.value

        # check if any galaxy is visible at the moment
        thismask = thisaltaz.alt.value > 90 - maxz

        nGals = len(galaxies[thismask])

        # print('nGals',nGals)

        if nGals == 0:
            # print("No galaxies visible within {0} deg zenith angle --> AFTERGLOW".format(maxz))

            return False, thisaltaz, galaxies
        else:
            # print("{0} galaxies are visible within {1} deg zenith angle ""--> Test for prompt follow up".format(nGals, maxz))

            return True, thisaltaz, galaxies
    else:
        thisaltaz = []
        return False, thisaltaz, galaxies


def FulfillsRequirement(theseGals, obspar, UsePix):
    """
    Apply filter criteria to visible galaxy sample and compares them to get the best option of zenith angle
    """

    maxz = obspar.maxZenith
    FOV = obspar.FOV
    zenithWeighting = obspar.zenithWeighting

    maxp = 1
    mask = 0
    thisminz = 0

    alt = theseGals["Alt"]

    if obspar.strategy == "targeted":
        testedProba = "dp_dV"
    if obspar.strategy == "integrated":
        testedProba = "dp_dV_FOV"

    for minz_aux in range(maxz, 5, -5):

        tmpmask = alt > 90 - (minz_aux)
        tmpGals = theseGals.copy()

        if len(tmpGals[tmpmask]) > 0:
            cur_maxp = (
                tmpGals[tmpmask][testedProba].max() / theseGals[testedProba].max()
            )
            if maxz == minz_aux:
                maxp = cur_maxp
                mask = tmpmask
                thisminz = minz_aux
            if cur_maxp > zenithWeighting * maxp:
                mask = tmpmask
                thisminz = minz_aux
            else:
                thisminz = minz_aux + 5
                break
    if UsePix:
        mask = alt > 90 - (thisminz + FOV)
    return mask, thisminz


def FulfillsRequirementGreyObservations(Ktime, theseGals, obspar):
    observatory = obspar.location
    minMoonSourceSeparation = obspar.minMoonSourceSeparation
    targetCoord = co.SkyCoord(
        theseGals["RAJ2000"], theseGals["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )
    frame = co.AltAz(obstime=Ktime, location=observatory)
    moonaltazs = get_body("moon", Time(Ktime, scale="utc")).transform_to(frame)

    altaz_map = targetCoord.transform_to(frame)
    separations = altaz_map.separation(moonaltazs)

    # Mask
    greymask = separations > minMoonSourceSeparation * u.deg
    return greymask


def GetEarthOccultedPix(
    nside, time, earth_radius, earth_sep, satellite_position, satellite_location
):
    # for equatorial frame

    distance_to_satellite = np.linalg.norm(satellite_position)
    earth_altitude = np.arcsin(
        -satellite_position[2] / distance_to_satellite
    )  # Altitude in radians
    earth_azimuth = np.arctan2(
        -satellite_position[1], -satellite_position[0]
    )  # Azimuth in radians

    angle_of_occlusion = np.arcsin(earth_radius / distance_to_satellite)

    altaz_frame = AltAz(obstime=time, location=satellite_location)
    earthCoord = SkyCoord(
        earth_azimuth * u.rad, earth_altitude * u.rad, frame=altaz_frame
    )
    earthCoord_equatorial = earthCoord.transform_to("icrs")

    earth_phipix = np.deg2rad(earthCoord_equatorial.ra.deg)
    earth_thetapix = 0.5 * np.pi - np.deg2rad(earthCoord_equatorial.dec.deg)
    earth_xyzpix = hp.ang2vec(earth_thetapix, earth_phipix)
    earth_occulted_pix = hp.query_disc(
        nside, earth_xyzpix, np.deg2rad(np.rad2deg(angle_of_occlusion) + earth_sep)
    )
    return earth_occulted_pix, earthCoord_equatorial


def GetMoonOccultedPix(nside, moon_sep, time):
    time = Time(time)
    # for equatorial frame
    MoonCoord_equatorial = get_body("moon", time, location=EarthLocation(0, 0, 0))
    MoonCoord_equatorial = MoonCoord_equatorial.transform_to("icrs")
    phipix_moon = np.deg2rad(MoonCoord_equatorial.ra.deg)
    thetapix_moon = 0.5 * np.pi - np.deg2rad(MoonCoord_equatorial.dec.deg)
    moon_xyzpix = hp.ang2vec(thetapix_moon, phipix_moon)
    moon_occulted_pix = hp.query_disc(nside, moon_xyzpix, np.deg2rad(moon_sep))
    return moon_occulted_pix, MoonCoord_equatorial


def GetSunOccultedPix(nside, sun_sep, time):
    time = Time(time)
    # for equatorial frame
    SunCoord_equatorial = get_body("sun", time, location=EarthLocation(0, 0, 0))
    SunCoord_equatorial = SunCoord_equatorial.transform_to("icrs")
    phipix_sun = np.deg2rad(SunCoord_equatorial.ra.deg)
    thetapix_sun = 0.5 * np.pi - np.deg2rad(SunCoord_equatorial.dec.deg)
    sun_xyzpix = hp.ang2vec(thetapix_sun, phipix_sun)
    sun_occulted_pix = hp.query_disc(nside, sun_xyzpix, np.deg2rad(sun_sep))
    return sun_occulted_pix, SunCoord_equatorial


def OccultationCut(
    prob,
    nside,
    time,
    minProbcut,
    satellite_position,
    observatory,
    sun_sep,
    moon_sep,
    earth_sep,
):
    """
    Mask in the pixels that are occulted by Earth, Sun and Moon
    """
    pixlist = []
    pprob = prob

    mOcc = hp.ma(pprob)
    maskOcc = np.zeros(hp.nside2npix(nside), dtype=bool)
    mpixels = []

    mEarth, posEarth = GetEarthOccultedPix(
        nside, time, 6371, earth_sep, satellite_position, observatory
    )
    mpixels.extend(mEarth)

    mSun, posSun = GetSunOccultedPix(nside, sun_sep, time)
    mpixels.extend(mSun)

    mMoon, poSMoon = GetMoonOccultedPix(nside, moon_sep, time)
    mpixels.extend(mMoon)

    pixlist.extend(mpixels)

    maskOcc[mSun] = 1
    maskOcc[mMoon] = 1
    maskOcc[mEarth] = 1

    mOcc.mask = maskOcc

    yprob = ma.masked_array(pprob, mOcc.mask)

    if np.sum(yprob) < minProbcut:
        ObsBool = False
    else:
        ObsBool = True

    return ObsBool, yprob, pixlist


def SAA_Times(
    duration,
    start_time,
    current_time,
    SatelliteName,
    saa,
    SatTimes,
    step,
    doPlot,
    dirName,
    datasetDir,
    saathershold,
):
    SatTimes = []
    i = 0
    saa = []
    while current_time <= start_time + datetime.timedelta(minutes=duration):
        SatelliteTime = GetSatelliteTime(SatelliteName, current_time)
        # satellite_position, satellite_location = GetSatellitePositions(
        #    SatelliteName, SatelliteTime
        # )
        # if Tools.is_in_saa(satellite_location.lat.deg, satellite_location.lon.deg):
        #    saa.append(True)
        # else:
        #    saa.append(False)
        saa.append(
            Tools.is_in_saa_opt(SatelliteName, SatelliteTime, saathershold, datasetDir)
        )
        SatTimes.append(current_time)

        current_time += step
        i += 1

    saa_numeric = [1 if value else 0 for value in saa]
    # Plot
    if doPlot:
        path = dirName + "/SAAPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)
        plt.figure(figsize=(12, 6))
        plt.plot(
            SatTimes, saa_numeric, drawstyle="steps-post", label="SAA (True=1, False=0)"
        )
        # plt.title("SAA Times")
        plt.xlabel("Time")
        plt.ylabel("SAA Status")
        plt.ylim(-0.1, 1.1)  # Set limits to make binary values clear
        plt.legend()
        plt.grid(True)
        plt.savefig("%s/SAA_Times.png" % (path))
    return SatTimes, saa


def GetBestGridPos2D(
    prob,
    highres,
    HRnside,
    reducedNside,
    newpix,
    radius,
    maxRuns,
    Occultedpixels,
    doPlot,
    dirName,
    n_sides,
    ipixlistHR,
    minProbcut,
):

    dp_dV_FOV = []
    ra = []
    dec = []
    newpixfinal = []
    sum_PGW = []
    prob1 = prob[newpix]
    newpix = newpix[np.argsort(prob1)[::-1]]
    rotation = 0

    for i in range(0, len(newpix)):
        if n_sides == 0:
            xyzpix = hp.pix2vec(reducedNside, newpix[i])
            # xyzpix = np.column_stack(xyzpix1)
            ipix_discfull = hp.query_disc(HRnside, xyzpix, np.deg2rad(radius))
        elif n_sides > 0:
            theta, phi = hp.pix2ang(reducedNside, newpix[i])
            ra_center = np.rad2deg(phi)
            dec_center = 90 - np.rad2deg(theta)
            vertices = Tools.get_regular_polygon_vertices(
                ra_center, dec_center, radius, n_sides, rotation
            )
            ipix_discfull = hp.query_polygon(HRnside, vertices, inclusive=True)
        else:
            raise ValueError("Shape must be 'circle' or 'polygon'.")

        if len(ipixlistHR) == 0:
            # No mask needed
            HRprob = highres[ipix_discfull].sum()
            HRprobf = highres[ipix_discfull].sum()
            ipixlistHR.extend(ipix_discfull)
        else:
            # Mask the ipix_discfull with the pixels that are already observed. I think the problem is here
            maskComputeProb = np.isin(ipix_discfull, ipixlistHR, invert=True)
            # Obtain list of pixel ID after the mask what has been observed already
            m_ipix_discfull = ma.compressed(
                ma.masked_array(ipix_discfull, mask=np.logical_not(maskComputeProb))
            )
            HRprob = highres[m_ipix_discfull].sum()
            HRprobf = highres[ipix_discfull].sum()
            ipixlistHR.extend(m_ipix_discfull)

        if HRprob > minProbcut:
            sum_PGW.append(HRprob)
            dp_dV_FOV.append(HRprobf)
            newpixfinal.append(newpix[i])
            theta, phi = hp.pix2ang(reducedNside, newpix[i])
            ra.append(np.degrees(phi))  # RA in degrees
            dec.append(90 - np.degrees(theta))

    if len(dp_dV_FOV) > 0:
        print("sum_PGW", sum(sum_PGW))
        cat_pix = Table(
            [newpixfinal, ra, dec, dp_dV_FOV],
            names=("PIX", "PIXRA", "PIXDEC", "PIXFOVPROB"),
        )

        sortcat1 = cat_pix[np.flipud(np.argsort(cat_pix["PIXFOVPROB"]))]
        first_values = sortcat1[:maxRuns]

    else:
        raise ValueError("No pointing were found with current minProbCut")

    if doPlot:
        hp.gnomview(
            highres,
            rot=(first_values["PIXRA"][0], first_values["PIXDEC"][0]),
            xsize=500,
            ysize=500,
        )

        path = dirName + "/GridPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)

        if n_sides > 0:
            for ra1, dec1 in zip(first_values["PIXRA"], first_values["PIXDEC"]):
                # Get unit vector vertices (assumed in Cartesian coords)
                vertices_radec = Tools.get_regular_polygon_vertices(
                    ra1, dec1, radius, n_sides, rotation
                )
                # Convert Cartesian to angular (theta, phi)
                theta, phi = hp.vec2ang(vertices_radec)  # in radians
                ra_deg = np.degrees(phi)
                dec_deg = 90 - np.degrees(theta)

                hp.projplot(
                    ra_deg,
                    dec_deg,
                    "r",
                    markersize=4,
                    lonlat=True,
                    coord="C",  # RA/Dec mode
                )
                hp.projplot(
                    np.append(ra_deg, ra_deg[0]),
                    np.append(dec_deg, dec_deg[0]),
                    "r-",
                    linewidth=1,
                    lonlat=True,
                    coord="C",
                )

        hp.graticule()
        try:
            tt, pp = hp.pix2ang(reducedNside, Occultedpixels)
            ra2 = np.rad2deg(pp)
            dec2 = np.rad2deg(0.5 * np.pi - tt)
            skycoord = co.SkyCoord(ra2, dec2, frame="fk5", unit=(u.deg, u.deg))
            hp.visufunc.projplot(
                skycoord.ra.deg,
                skycoord.dec.deg,
                "g.",
                lonlat=True,
                coord="C",
                linewidth=0.1,
            )
        except Exception:
            print("No occulted pix")

        hp.visufunc.projplot(
            first_values["PIXRA"],
            first_values["PIXDEC"],
            "b.",
            lonlat=True,
            coord="C",
            linewidth=0.1,
        )

        plt.savefig("%s/Grid_Pointing.png" % (path), bbox_inches="tight")
        plt.close()

    return first_values


def GetBestGridPos3D(
    prob,
    cat,
    galpix,
    newpix,
    FOV,
    totaldPdV,
    HRnside,
    n_sides,
    maxRuns,
    doPlot,
    dirName,
    reducedNside,
    Occultedpixels,
    minProbCut,
):

    prob1 = prob[newpix]
    galpix = galpix[np.argsort(prob1)[::-1]]
    SelectedGals = galpix
    dp_dV_FOV = []
    BestGalsRA = []
    BestGalsDec = []
    # galaxx = []
    cat0 = cat
    for element in range(0, len(SelectedGals)):
        if element < len(SelectedGals):
            dp_dV_FOV1, galax = ComputePGalinFOV(
                prob,
                cat,
                SelectedGals[element],
                FOV,
                totaldPdV,
                n_sides,
                UsePix=True,
            )
        if dp_dV_FOV1 > minProbCut:
            dp_dV_FOV2, galax2 = ComputePGalinFOV(
                prob, cat0, SelectedGals[element], FOV, totaldPdV, n_sides, UsePix=True
            )
            dp_dV_FOV.append(dp_dV_FOV2)
            BestGalsRA.append(SelectedGals[element].ra.deg)
            BestGalsDec.append(SelectedGals[element].dec.deg)
            cat = galax

    if len(dp_dV_FOV) > 0:
        print("sum(dp_dV_FOV)", sum(dp_dV_FOV))
        cat_pix = Table(
            [BestGalsRA, BestGalsDec, dp_dV_FOV],
            names=("PIXRA", "PIXDEC", "PIXFOVPROB"),
        )

        sortcat = cat_pix[np.flipud(np.argsort(cat_pix["PIXFOVPROB"]))]
        first_values = sortcat[:maxRuns]

    else:
        raise ValueError("No pointing were found with current minProbCut")

    if doPlot:
        path = dirName + "/GridPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)

        # mpl.rcParams.update({'font.size':14})
        hp.gnomview(
            prob,
            rot=(first_values["PIXRA"][0], first_values["PIXDEC"][0]),
            xsize=500,
            ysize=500,
        )
        hp.graticule()

        # Filter out rows with NaN in dp_dV
        mask = ~np.isnan(cat["dp_dV"])
        cat_clean = cat[mask]

        # Sort descending by dp_dV
        cat_sorted = cat_clean.copy()
        cat_sorted.sort("dp_dV", reverse=True)
        # Select top 100 galaxies with highest dp_dV
        top100 = cat_sorted[:500]

        # Plot with projplot
        hp.visufunc.projplot(
            top100["RAJ2000"],
            top100["DEJ2000"],
            "r.",
            lonlat=True,
            coord="C",
            linewidth=0.1,
        )

        try:
            tt, pp = hp.pix2ang(reducedNside, Occultedpixels)
            ra2 = np.rad2deg(pp)
            dec2 = np.rad2deg(0.5 * np.pi - tt)
            skycoord = co.SkyCoord(ra2, dec2, frame="fk5", unit=(u.deg, u.deg))
            hp.visufunc.projplot(
                skycoord.ra.deg,
                skycoord.dec.deg,
                "g.",
                lonlat=True,
                coord="C",
                linewidth=0.1,
            )

        except Exception:
            print("No occulted pix")

        hp.visufunc.projplot(
            first_values["PIXRA"],
            first_values["PIXDEC"],
            "b.",
            lonlat=True,
            coord="C",
            linewidth=0.1,
        )
        plt.savefig("%s/Grid_Pointing.png" % (path))
        plt.close()

    return first_values
