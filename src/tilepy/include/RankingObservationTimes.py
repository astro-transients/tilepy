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
#           Tools to rank the observations from the largest to the lowest probability covered,
#           adding the observability window of each of them, which gives a comprehensive view
##################################################################################################

import datetime

import astropy.coordinates as co
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six
from astropy import units as u
from astropy.coordinates import AltAz, SkyCoord, get_body
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from six.moves import configparser
from sklearn.cluster import AgglomerativeClustering

from .PointingTools import FilterGalaxies, Tools

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


# iers_file = os.path.join(os.path.abspath(
#    os.path.dirname(__file__)), '../dataset/finals2000A.all')
# iers.IERS.iers_table = iers.IERS_A.open(iers_file)

__all__ = [
    "load_pointingFile",
    "VisibilityWindow",
    "GetObservationPeriod",
    "GetVisibility",
    "ProbabilitiesinPointings3D",
    "PGGPGalinFOV",
    "ProbabilitiesinPointings2D",
    "PGinFOV",
    "Sortingby",
    "EvolutionPlot",
    "RankingTimes",
    "RankingTimes_2D",
]


def load_pointingFile(tpointingFile):
    # Read PointingsFile

    print("Loading pointings from " + tpointingFile)
    (
        time1,
        time2,
        ra,
        dec,
    ) = np.genfromtxt(
        tpointingFile,
        usecols=(0, 1, 2, 3),
        dtype="str",
        skip_header=1,
        delimiter=" ",
        unpack=True,
    )  # ra, dec in degrees
    time1 = np.atleast_1d(time1)
    time2 = np.atleast_1d(time2)

    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    time = []

    for i in range(len(time1)):
        time.append(
            (
                time1[i] + " " + time2[i].split(":")[0] + ":" + time2[i].split(":")[1]
            ).split('"')[1]
        )

    ra = ra.astype(float)
    dec = dec.astype(float)

    index_list = list(range(len(ra)))
    Pointings = Table(
        [index_list, time, ra, dec], names=("Pointing", "Time", "RAJ2000", "DEJ2000")
    )

    return Pointings


def VisibilityWindow(ObservationTime, Pointing, obspar, dirName):

    source = SkyCoord(
        Pointing["RAJ2000"], Pointing["DEJ2000"], frame="icrs", unit=(u.deg, u.deg)
    )
    WINDOW = []
    ZENITH = []
    SZENITH = []

    try:
        auxtime = datetime.datetime.strptime(
            Pointing["Time"][0], "%Y-%m-%d %H:%M:%S.%f"
        )
    except ValueError:
        try:
            auxtime = datetime.datetime.strptime(
                Pointing["Time"][0], "%Y-%m-%d %H:%M:%S"
            )
        except ValueError:
            auxtime = datetime.datetime.strptime(Pointing["Time"][0], "%Y-%m-%d %H:%M")

    # frame = co.AltAz(obstime=auxtime, location=observatory)
    timeInitial = auxtime - datetime.timedelta(minutes=30)
    for i in range(0, len(source)):
        NonValidwindow, Stepzenith = GetVisibility(
            Pointing["Time"], source[i], obspar.maxZenith, obspar.location
        )
        window, zenith = GetObservationPeriod(
            timeInitial, source[i], obspar, i, dirName, False
        )
        WINDOW.append(window)
        ZENITH.append(zenith)
        SZENITH.append(Stepzenith)

        # At input ObservationTime the night is over, the scheduling has been computed for the next night (with the condition <24h holding)
        if not Tools.IsGreyness(ObservationTime, obspar):
            window, zenith = GetObservationPeriod(
                ObservationTime + datetime.timedelta(hours=12),
                source[i],
                obspar,
                i,
                dirName,
                True,
            )
        else:
            window, zenith = GetObservationPeriod(
                ObservationTime, source[i], obspar, i, dirName, True
            )

    Pointing["Observation window"] = WINDOW
    Pointing["Array of zenith angles"] = ZENITH
    Pointing["Zenith angles in steps"] = SZENITH

    return Pointing


def GetObservationPeriod(inputtime0, msource, obspar, plotnumber, dirName, doPlot):

    AltitudeCut = 90 - obspar.maxZenith
    nights = obspar.maxNights
    useGreytime = obspar.useGreytime
    moonGrey = obspar.moonGrey
    moonDown = obspar.moonDown
    moonPhase = obspar.moonPhase
    sunDown = obspar.sunDown
    minMoonSourceSeparation = obspar.minMoonSourceSeparation
    maxMoonSourceSeparation = obspar.maxMoonSourceSeparation

    inputtime = Time(inputtime0)
    initialframe = AltAz(obstime=inputtime, location=obspar.location)

    ##############################################################################
    suninitial = get_body("sun", inputtime).transform_to(initialframe)

    if suninitial.alt < -18.0 * u.deg:
        hoursinDay = 12
    else:
        hoursinDay = 24
    delta_day = np.linspace(0, hoursinDay + 24 * (nights - 1), 1000 * nights) * u.hour
    interval = (hoursinDay + 24 * (nights - 1)) / (1000.0 * nights)

    x = np.arange(int(hoursinDay / interval), dtype=int)
    firstN = np.full_like(x, 1)
    ratio2 = 24.0 / interval
    otherN = []
    for i in range(2, nights + 1):
        otherN.extend(np.full_like(np.arange(int(ratio2)), i))
    NightsCounter = []
    NightsCounter.extend(firstN)
    NightsCounter.extend(otherN)

    while len(NightsCounter) != len(delta_day):
        NightsCounter.extend([nights])

    times = inputtime + delta_day
    frame = AltAz(obstime=times, location=obspar.location)

    ##############################################################################
    # SUN
    sunaltazs = get_body("sun", times).transform_to(frame)

    # MOON
    moon = get_body("moon", times)
    moonaltazs = moon.transform_to(frame)
    msourcealtazs = msource.transform_to(frame)

    # Add Moon phase
    moonPhase = np.full(len(msourcealtazs), Tools.MoonPhase(inputtime0, obspar))

    MoonDistance = msourcealtazs.separation(moonaltazs)
    ##############################################################################
    if useGreytime:
        Altitudes = Table(
            [
                times,
                msourcealtazs.alt,
                sunaltazs.alt,
                moonaltazs.alt,
                moonPhase,
                MoonDistance,
                NightsCounter,
            ],
            names=[
                "Time UTC",
                "Alt Source",
                "Alt Sun",
                "AltMoon",
                "moonPhase",
                "MoonDistance",
                "NightsCounter",
            ],
        )
        # selectedTimes=Altitudes['Time UTC']
        selection = (
            (Altitudes["Alt Sun"] < -18.0)
            & (Altitudes["Alt Source"] > AltitudeCut)
            & (Altitudes["AltMoon"] < -0.5)
        )
        DTaltitudes = Altitudes[selection]
        newtimes = []
        newtimes.extend(DTaltitudes["Time UTC"].mjd)
        selectionGreyness = (
            (Altitudes["AltMoon"] < moonGrey)
            & (Altitudes["AltMoon"] > moonDown)
            & (Altitudes["moonPhase"] < moonPhase)
            & (Altitudes["Alt Sun"] < sunDown)
            & (Altitudes["MoonDistance"] > minMoonSourceSeparation)
            & (Altitudes["MoonDistance"] < maxMoonSourceSeparation)
            & (Altitudes["Alt Source"] > AltitudeCut)
        )
        GTaltitudes = Altitudes[selectionGreyness]
        newtimes.extend(GTaltitudes["Time UTC"].mjd)
        newtimes = sorted(newtimes)
        ScheduledTimes = Time(newtimes, format="mjd").iso

    else:
        Altitudes = Table(
            [times, msourcealtazs.alt, sunaltazs.alt, moonaltazs.alt, NightsCounter],
            names=["Time UTC", "Alt Source", "Alt Sun", "AltMoon", "NightsCounter"],
        )
        Times = Altitudes["Time UTC"]
        selection = (
            (Altitudes["Alt Sun"] < -18.0)
            & (Altitudes["Alt Source"] > AltitudeCut)
            & (Altitudes["AltMoon"] < -0.5)
        )
        ScheduledTimes = Time(Times[selection], format="mjd").iso

    try:
        return (
            str(ScheduledTimes[0]).split(".")[0]
            + "-->"
            + str(ScheduledTimes[-1]).split(".")[0]
        ), msourcealtazs.alt
    except Exception:
        ScheduledTimesUni = str(ScheduledTimes).split(".")
        ScheduledTimes1 = ScheduledTimesUni[0]
        ScheduledTimes2 = ScheduledTimesUni[-1]
        return (str(ScheduledTimes1) + "-->" + str(ScheduledTimes2)), msourcealtazs.alt


def GetVisibility(time, radecs, maxZenith, obsLoc):

    visibility = []
    altitude = []

    for i in range(0, len(time)):
        try:
            auxtime = datetime.datetime.strptime(time[i], "%Y-%m-%d %H:%M:%S.%f")
        except ValueError:
            try:
                auxtime = datetime.datetime.strptime(time[i], "%Y-%m-%d %H:%M:%S")
            except ValueError:
                auxtime = datetime.datetime.strptime(time[i], "%Y-%m-%d %H:%M")
        frame = co.AltAz(obstime=auxtime, location=obsLoc)
        thisaltaz = radecs.transform_to(frame)
        visible = thisaltaz.alt.value > (90 - maxZenith)

        if visible:
            visibility.append(auxtime)
            altitude.append(thisaltaz.alt.value)
        else:
            #    visibility.append(auxtime)
            altitude.append(thisaltaz.alt.value)
    lasttime = auxtime + datetime.timedelta(minutes=30)
    frame = co.AltAz(obstime=lasttime, location=obsLoc)
    thisaltaz = radecs.transform_to(frame)
    visible = thisaltaz.alt.value > (90 - maxZenith)

    if visible:
        visibility.append(auxtime)
        altitude.append(thisaltaz.alt.value)
    else:
        #    visibility.append(auxtime)
        altitude.append(thisaltaz.alt.value)

    window = (
        visibility[0].strftime("%H:%M:%S") + "-" + visibility[-1].strftime("%H:%M:%S")
    )
    return window, altitude


def ProbabilitiesinPointings3D(cat, galPointing, FOV, totaldPdV, prob, nside):

    ra = galPointing["RAJ2000"]
    dec = galPointing["DEJ2000"]
    PGW = []
    PGAL = []

    # bucle
    for i in range(0, len(ra)):
        pgwcircle, pgalcircle = PGGPGalinFOV(
            cat, ra[i], dec[i], prob, totaldPdV, FOV, nside
        )
        PGW.append(float("{:1.4f}".format(pgwcircle)))
        PGAL.append(float("{:1.4f}".format(pgalcircle)))

    galPointing["Pgw"] = PGW
    galPointing["Pgal"] = PGAL

    return galPointing


def PGGPGalinFOV(cat, ra, dec, prob, totaldPdV, FOV, nside):

    targetCoordcat = co.SkyCoord(
        cat["RAJ2000"], cat["DEJ2000"], frame="icrs", unit=(u.deg, u.deg)
    )
    targetCoordpointing = co.SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))
    dp_dV = cat["dp_dV"]

    # Array of indices of pixels inside circle of FoV

    radius = FOV
    t = 0.5 * np.pi - targetCoordpointing.dec.rad
    p = targetCoordpointing.ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)
    xyz = hp.ang2vec(t, p)

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    P_GW = prob[ipix_disc].sum()
    Pgal_inFoV = (
        dp_dV[targetCoordcat.separation(targetCoordpointing).deg <= radius].sum()
        / totaldPdV
    )

    return P_GW, Pgal_inFoV


def ProbabilitiesinPointings2D(Pointing, FOV, prob, nside):

    ra = Pointing["RAJ2000"]
    dec = Pointing["DEJ2000"]
    PGW = []
    PGAL = []
    for i in range(0, len(ra)):
        pgwcircle = PGinFOV(ra[i], dec[i], prob, FOV, nside)
        PGW.append(float("{:1.4f}".format(pgwcircle)))
        PGAL.append(float("{:1.4f}".format(0)))

    Pointing["Pgw"] = PGW
    Pointing["Pgal"] = PGAL

    return Pointing


def PGinFOV(ra, dec, prob, radius, nside):

    targetCoordpointing = co.SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))

    # Array of indices of pixels inside circle of FoV

    t = 0.5 * np.pi - targetCoordpointing.dec.rad
    p = targetCoordpointing.ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)
    xyz = hp.ang2vec(t, p)

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    P_GW = prob[ipix_disc].sum()

    return P_GW


def Sortingby(galPointing, name, exposure):

    gggalPointing = galPointing[np.flipud(np.argsort(galPointing["Pgal"]))]
    prioritygal = list(range(len(galPointing["Pgal"])))
    ra = gggalPointing["RAJ2000"]
    dec = gggalPointing["DEJ2000"]
    coord = SkyCoord(ra, dec, unit="deg")
    # print(coord.to_string('hmsdms'))
    gggalPointing["RA(HH:MM:SS) Dec (DD:MM:SS)"] = coord.to_string("hmsdms")
    gggalPointing["PriorityGal"] = prioritygal
    gggalPointing.remove_column("Array of zenith angles")
    gggalPointing.remove_column("Zenith angles in steps")

    # Prepare filename which is going to be complete
    gwgalPointing = gggalPointing[np.flipud(np.argsort(gggalPointing["Pgw"]))]
    prioritygw = list(range(len(galPointing["Pgw"])))
    gwgalPointing["PriorityGW"] = prioritygw

    # gwgalPointing.remove_column('Array of zenith angles')
    # gwgalPointing.remove_column('Zenith angles in steps')
    # print(gwgalPointing)
    outfilename = "%s/RankingObservationTimes_Complete.txt" % name
    ascii.write(
        gwgalPointing[np.argsort(gwgalPointing["Pointing"])],
        outfilename,
        overwrite=True,
    )

    gwgalPointing.remove_column("Pgal")
    gwgalPointing.remove_column("Pgw")
    gwgalPointing.remove_column("PriorityGW")
    gwgalPointing.rename_column("PriorityGal", "Priority")
    gwgalPointing.remove_column("RA(HH:MM:SS) Dec (DD:MM:SS)")
    outfilename = "%s/RankingObservationTimes_forShifters.txt" % name
    ascii.write(
        gwgalPointing[np.argsort(gwgalPointing["Pointing"])],
        outfilename,
        overwrite=True,
    )

    gwgalPointing.remove_column("Observation window")
    gwgalPointing.remove_column("Priority")
    print(name)

    target = [
        (name.split("/")[2] + "_{0}").format(element)
        for element in gwgalPointing["Pointing"]
    ]
    gwgalPointing["Target"] = target
    gwgalPointing.rename_column("Pointing", "Id")
    gwgalPointing["Duration"] = exposure

    gwgalPointing_TH = gwgalPointing[
        "Target", "Id", "RAJ2000", "DEJ2000", "Time", "Duration"
    ]
    # new_order = ['Target', 'Id', 'RAJ2000','DEJ2000']  # List or tuple
    # gwgalPointing_TH = gwgalPointing[new_order]
    outfilename = "%s/RankingObservationTimes_forAlerter.txt" % name
    ascii.write(
        gwgalPointing_TH[np.argsort(gwgalPointing_TH["Id"])],
        outfilename,
        overwrite=True,
    )


def EvolutionPlot(galPointing, tname, ObsArray):

    fig = plt.figure(figsize=(18, 10))
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.8])
    ra = galPointing["RAJ2000"]
    dec = galPointing["DEJ2000"]
    pgw = galPointing["Pgw"]
    pgal = galPointing["Pgal"]
    time = galPointing["Time"]
    NUM_COLORS = len(time)
    hour = []
    for j in range(0, len(time)):
        selecttime = time[j].split(" ")
        hour.append(selecttime[1].split(".")[0])
    try:
        lasttime = datetime.datetime.strptime(
            time[len(time) - 1], "%Y-%m-%d %H:%M"
        ) + datetime.timedelta(minutes=30)
    except ValueError:
        lasttime = datetime.datetime.strptime(
            time[len(time) - 1], "%Y-%m-%d %H:%M"
        ) + datetime.timedelta(minutes=30)

    hour.append(lasttime.strftime("%H:%M"))
    GWordered = galPointing[np.flipud(np.argsort(galPointing["Pgw"]))]

    ax.set_prop_cycle(plt.cycler("color", plt.cm.Accent(np.linspace(0, 1, NUM_COLORS))))
    for i in range(0, len(ra)):
        # x = np.arange(0, len(ra), 1)
        ZENITH = GWordered["Zenith angles in steps"][i]
        x = np.arange(0, len(ZENITH), 1)
        ax.plot(
            x,
            ZENITH,
            label="ra:%.2f dec:%.2f- Pgw:%.3f - Pgal:%.3f "
            % (ra[i], dec[i], 100 * pgw[i], 100 * pgal[i]),
        )

    ax.set_xticks(x)
    ax.set_xticklabels(hour)
    ax.set_ylabel("Altitude (deg)", fontsize=14)
    ax.set_xlabel("Time", fontsize=14)
    ax.grid()
    ax.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
    plt.savefig("%s/AltitudevsTime_%s.png" % (tname, ObsArray))


def RankingTimes(ObservationTime, skymap, cat, obspar, dirName, PointingFile, ObsArray):

    point = load_pointingFile(PointingFile)

    ################################################################

    print()
    print(
        "---------  RANKING THE OBSERVATIONS AND PRODUCING THE OUTPUT FILES   ----------"
    )
    print()

    nside = obspar.HRnside
    prob = skymap.getMap("prob", obspar.HRnside)

    # correlate GW map with galaxy catalog, retrieve ordered list
    cat = skymap.computeGalaxyProbability(cat)
    tGals = FilterGalaxies(cat, obspar.minimumProbCutForCatalogue)
    sum_dP_dV = cat["dp_dV"].sum()
    point = ProbabilitiesinPointings3D(tGals, point, obspar.FOV, sum_dP_dV, prob, nside)
    point = VisibilityWindow(ObservationTime, point, obspar, dirName)
    EvolutionPlot(point, dirName, ObsArray)
    Sortingby(point, dirName, obspar.duration)


def RankingTimes_2D(ObservationTime, prob, obspar, dirName, PointingFile, ObsArray):

    point = load_pointingFile(PointingFile)

    ################################################################

    print()
    print(
        "---------  RANKING THE OBSERVATIONS AND PRODUCING THE OUTPUT FILES   ----------"
    )
    print()

    npix = len(prob)
    nside = hp.npix2nside(npix)

    # In this function, the 2D probability is computed and the 3D probability is set to zero
    point = ProbabilitiesinPointings2D(point, obspar.FOV, prob, nside)
    point = VisibilityWindow(ObservationTime, point, obspar, dirName)

    EvolutionPlot(point, dirName, ObsArray)
    Sortingby(point, dirName, obspar.duration)


# Function to compute 2D distance between two rows
def distance(entry1, entry2):
    ra1, dec1 = entry1["RA(deg)"], entry1["DEC(deg)"]
    ra2, dec2 = entry2["RA(deg)"], entry2["DEC(deg)"]

    # Handle circular distance for RA
    delta_ra = min(abs(ra1 - ra2), 360 - abs(ra1 - ra2))
    delta_dec = abs(dec1 - dec2)

    # Euclidean distance in 2D
    return np.sqrt(delta_ra**2 + delta_dec**2)


# Ranking function
def Ranking_Space(dirName, PointingFile):
    # Read the data from the pointing file
    file_path = f"{PointingFile}"
    data = pd.read_csv(file_path, delim_whitespace=True)

    # Sort by PGW in descending order
    try:
        data = data.sort_values(by="PGW", ascending=False).reset_index(drop=True)
    except Exception:
        data = data.sort_values(by="PGal", ascending=False).reset_index(drop=True)

    # Initialize ranked list with the first (highest PGW) entry
    ranked = [data.iloc[0]]
    data = data.iloc[1:].reset_index(drop=True)  # Exclude the first entry

    # Iteratively find the closest entry
    while not data.empty:
        last_entry = ranked[-1]
        # Compute distances to the last entry
        data["distance"] = data.apply(lambda row: distance(last_entry, row), axis=1)
        # Find the closest entry
        closest_idx = data["distance"].idxmin()
        closest_entry = data.loc[closest_idx]
        ranked.append(closest_entry)
        # Remove the closest entry from the dataset
        data = data.drop(index=closest_idx).reset_index(drop=True)

    # Output the ranked list
    print("Ranked List:")
    for idx, entry in enumerate(ranked, start=1):
        print(f"Rank {idx}: {entry.to_dict()}")

    # Save the ranked list to a file
    output_file = "%s/RankingObservations_Space.txt" % dirName
    pd.DataFrame(ranked).to_csv(output_file, index=False, sep="\t")
    print(f"Ranked file saved to {output_file}")


def Ranking_Space_AI(dirName, PointingFile):
    # Convert to DataFrame for easier handling
    file_path = f"{PointingFile}"
    data = pd.read_csv(file_path, delim_whitespace=True)
    df = pd.DataFrame(data)

    # Extract RA and DEC for clustering
    coordinates = df[["RA(deg)", "DEC(deg)"]].to_numpy()

    # Clustering with Agglomerative Clustering
    clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=1.0
    )  # Distance can be adjusted
    df["Cluster"] = clustering.fit_predict(coordinates)

    # Sort within each cluster by PGW
    ranked_data = []
    for cluster_id in sorted(df["Cluster"].unique()):
        try:
            cluster_data = df[df["Cluster"] == cluster_id].sort_values(
                by="PGW", ascending=False
            )
        except Exception:
            cluster_data = df[df["Cluster"] == cluster_id].sort_values(
                by="PGal", ascending=False
            )
        ranked_data.append(cluster_data)

    # Combine ranked clusters
    final_ranked = pd.concat(ranked_data)

    # Save the ranked list to a file
    output_file = "%s/RankingObservations_AI_Space.txt" % dirName
    pd.DataFrame(final_ranked).to_csv(output_file, index=False, sep="\t")
    print(f"Ranked file saved to {output_file}")
