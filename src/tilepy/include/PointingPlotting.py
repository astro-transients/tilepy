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
#                        Various plotting tools to obtain plots with the scheduling              #
##################################################################################################

import datetime
import os

import astropy.coordinates as co
import healpy as hp
import ligo.skymap.io.fits as lf
import ligo.skymap.plot
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six
from astropy import units as u
from astropy.coordinates import SkyCoord
from six.moves import configparser

from .Observatories import CTANorthObservatory, CTASouthObservatory
from .PointingTools import TransformRADec

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser

from matplotlib import colors
from matplotlib.patches import Circle
import matplotlib.colors as mcolors

# iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

# iers_file = os.path.join(os.path.abspath(
#    os.path.dirname(__file__)), '../dataset/finals2000A.all')
# iers.IERS.iers_table = iers.IERS_A.open(iers_file)

# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))

__all__ = [
    "LoadPointingsGW",
    "LoadPointingsGAL",
    "LoadPointings",
    "PointingPlotting",
    "PlotPointings",
    "PlotPointingsTogether",
    "PointingPlottingGWCTA",
    "PlotPointings_Pretty",
]


Colors = ['b', 'm', 'y', 'c', 'g', 'w', 'k', 'c', 'b', 'c', 'm', 'b', 'g', 'y', 'b', 'c', 'm', 'b', 'g', 'y',
          'b', 'c', 'm', 'b','b', 'm', 'y', 'c', 'g', 'w', 'k', 'c', 'b', 'c', 'm', 'b', 'g', 'y', 'b', 'c', 'm', 'b', 'g', 'y',
          'b', 'c', 'm', 'b','b', 'm', 'y', 'c', 'g', 'w', 'k', 'c', 'b', 'c', 'm', 'b', 'g', 'y', 'b', 'c', 'm', 'b', 'g', 'y',
          'b', 'c', 'm', 'b','b', 'm', 'y', 'c', 'g', 'w', 'k', 'c', 'b', 'c', 'm', 'b', 'g', 'y', 'b', 'c', 'm', 'b', 'g', 'y','b', 'm', 'y', 'c', 'g', 'w', 'k', 'c', 'b', 'c', 'm', 'b', 'g', 'y', 'b', 'c', 'm', 'b', 'g', 'y',]


def LoadPointingsGW(tpointingFile):

    print("Loading pointings from " + tpointingFile)

    time1, time2, ra, dec = np.genfromtxt(tpointingFile, usecols=(0, 1, 2, 3), dtype="str", skip_header=1,
                                          delimiter=' ',
                                          unpack=True)  # ra, dec in degrees

    time1 = np.atleast_1d(time1)
    time2 = np.atleast_1d(time2)
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    time = []
    for i, time1 in enumerate(time1):
        time.append(time1.split('"')[1] + ' ' + time2[i].split('"')[0])

    ra = ra.astype(float)
    dec = dec.astype(float)
    coordinates = co.SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg))
    # pgw = Pgw.astype(float)
    pgw = np.genfromtxt(tpointingFile, usecols=4,
                        skip_header=1, delimiter=' ', unpack=True)
    return time, coordinates, pgw


def LoadPointingsGAL(tpointingFile):

    print("Loading pointings from " + tpointingFile)
    time1, time2, ra, dec, Pgw, Pgal = np.genfromtxt(tpointingFile, usecols=(0, 1, 2, 3, 4, 5), dtype="str", skip_header=1, delimiter=' ',
                                                     unpack=True)  # ra, dec in degrees
    time1 = np.atleast_1d(time1)
    time2 = np.atleast_1d(time2)
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    time = []
    for i, time1 in enumerate(time1):
        try:
            time.append((time1[i] + ' ' + time2[i]).split('"')[1])
        except IndexError:
            time.append((time1 + ' ' + time2).split('"')[1])
            break
    coordinates = TransformRADec(ra, dec)
    Pgw = Pgw.astype(float)
    Pgal = Pgal.astype(float)
    return time, coordinates, Pgw, Pgal


def LoadPointings(tpointingFile):
    print("Loading pointings from " + tpointingFile)
    # Read the first line of the file to determine column names
    with open(tpointingFile, 'r') as f:
        header_line = f.readline().strip()
    # Read the data into a DataFrame
    data = pd.read_csv(tpointingFile, delimiter=' ', header=None, names=header_line.split(), skiprows=1)
    return data


def PointingPlotting(prob, obspar, name, dirName, PointingsFile1, ObsArray, gal):
    npix = len(prob)
    nside = hp.npix2nside(npix)

    ObservationTimearray1, Coordinates1, Probarray1 = LoadPointingsGW(
        PointingsFile1)

    Probarray1 = np.atleast_1d(Probarray1)

    print('----------   PLOTTING THE SCHEDULING   ----------')
    print('Total covered probability with the scheduled tiles is PGW= {0:.5f}'.format(
        sum(Probarray1)))
    converted_time1 = []
    for i, time1 in enumerate(ObservationTimearray1):
        try:
            converted_time1.append(datetime.datetime.strptime(
                time1, '%Y-%m-%d %H:%M:%S.%f'))
        except ValueError:
            try:
                converted_time1.append(datetime.datetime.strptime(
                    time1, '%Y-%m-%d %H:%M:%S.%f '))
            except ValueError:
                try:
                    converted_time1.append(
                        datetime.datetime.strptime(time1, '%Y-%m-%d %H:%M:%S'))
                except ValueError:
                    converted_time1.append(time1)

    PlotPointings(prob, converted_time1, Coordinates1, sum(Probarray1), nside, obspar, name, dirName, ObsArray)
    PlotPointings_Pretty(prob, name, PointingsFile1, dirName, obspar, gal)


def PlotPointings(prob, time, targetCoord, Totalprob, nside, obspar, name, dirName, ObsArray):
    FOV = obspar.FOV
    maxzenith = obspar.maxZenith
    doPlot = obspar.doPlot

    t = 0.5 * np.pi - targetCoord[0].dec.rad
    p = targetCoord[0].ra.rad

    xyz = hp.ang2vec(t, p)

    # translate pixel indices to coordinates
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))

    if (doPlot):

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='icrs', unit=(u.deg, u.deg))
        observatory = obspar.location
        frame = co.AltAz(obstime=time[0], location=observatory)
        altaz_all = skycoord.transform_to(frame)

        dirName = '%s/Pointing_Plotting_%s' % (dirName, ObsArray)
        if not os.path.exists(dirName):
            os.makedirs(dirName)

        hp.mollview(prob, rot=[180, 0], coord='C', title="GW prob map (Equatorial) + %s %g  %s/%s/%s %s:%s:%s UTC" %
                                                         (name, Totalprob * 100, time[0].day, time[0].month,
                                                          time[0].year,
                                                          time[0].hour, time[0].minute, time[0].second))
        hp.graticule()

        theta = np.random.rand(400) * 360
        tarcoordra = np.empty(400)
        tarcoorddec = np.empty(400)
        Fov_array = np.empty(400)
        Fov_array.fill(FOV)

        # Position of the GRB: Add here your coordinates if you want to add them to the plot
        # RA_GRB = 0.526
        # DEC_GRB = -2.917
        # skycoordGRB = co.SkyCoord(RA_GRB, DEC_GRB, frame='icrs', unit=(u.deg, u.deg))

        for j in range(0, len(targetCoord.ra)):
            tarcoordra.fill(targetCoord[j].ra.deg)
            tarcoorddec.fill(targetCoord[j].dec.deg)
            racoord = tarcoordra + Fov_array * np.cos(theta)
            deccoord = tarcoorddec + Fov_array * np.sin(theta)
            hp.visufunc.projscatter(
                racoord, deccoord, lonlat=True, marker='.', color=Colors[j], coord='C')

            # Plotting the visibility of the instrument as specified in the config
            altcoord = np.empty(1000)
            altcoord.fill(90 - maxzenith)
            azcoord = np.random.rand(1000) * 360
            # print(time)
            RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time[j],
                                   location=observatory)
            RandomCoord_radec = RandomCoord.transform_to('icrs')
            hp.visufunc.projplot(RandomCoord_radec.ra,
                                 RandomCoord_radec.dec, 'b.', lonlat=True)
            plt.savefig('%s/Pointings%s.png' % (dirName, j))


def PlotPointingsTogether(prob, time, targetCoord1, n1, targetCoord2, n2, nside, FOV, plotType, doPlot=True):

    t = 0.5 * np.pi - targetCoord1[0].dec.rad
    p = targetCoord1[0].ra.rad

    xyz = hp.ang2vec(t, p)
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))

    if (doPlot):

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        # skycoord = co.SkyCoord(ra2, dec2, frame='icrs', unit=(u.deg, u.deg))
        # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

        # frame = co.AltAz(obstime=time, location=observatory)
        # altaz_all = skycoord.transform_to(frame)

        # path = os.path.dirname(os.path.realpath(__file__)) + '/Pointing_Plotting'
        # if not os.path.exists(path):
        #    os.mkdir(path, 493)
        if(plotType=='gnomonic'): 
            hp.gnomview(prob, xsize=1000, ysize=1000, rot=[targetCoord1[0].ra.deg, targetCoord1[0].dec.deg], reso=5.0)
        if(plotType=='mollweide'): 
            hp.mollview(prob,title="GW prob map (Equatorial)",coord='C')
        hp.graticule()

        hp.visufunc.projscatter(
            targetCoord1.ra.deg, targetCoord1.dec.deg, lonlat=True, marker='.', color=Colors[4])
        hp.visufunc.projscatter(
            targetCoord2.ra.deg, targetCoord2.dec.deg, lonlat=True, marker='.',color=Colors[5])

        npoints = 400
        Fov_array = np.empty(npoints)
        Fov_array.fill(FOV)

        theta = np.random.rand(npoints) * 360
        tarcoordra1 = np.empty(npoints)

        tarcoorddec1 = np.empty(npoints)
        for j in range(0, len(targetCoord1.ra)):
            tarcoordra1.fill(targetCoord1[j].ra.deg)
            tarcoorddec1.fill(targetCoord1[j].dec.deg)
            racoord1 = tarcoordra1 + Fov_array * np.cos(theta)
            deccoord1 = tarcoorddec1 + Fov_array * np.sin(theta)
            hp.visufunc.projscatter(racoord1, deccoord1, lonlat=True, marker='.', color=Colors[2])   
            hp.projtext(targetCoord1[j].ra, targetCoord1[j].dec, str(j), lonlat=True, color=Colors[1])
        # plt.savefig("Pointing_Plotting_%s/G274296Pointing_GW_%g_%g:%g.png" % (ObsArray, time.day, time.hour,time.minute))

        tarcoordra2 = np.empty(npoints)

        tarcoorddec2 = np.empty(npoints)
        for j in range(0, len(targetCoord2.ra)):
            tarcoordra2.fill(targetCoord2[j].ra.deg)
            tarcoorddec2.fill(targetCoord2[j].dec.deg)
            racoord2 = tarcoordra2 + Fov_array * np.cos(theta)
            deccoord2 = tarcoorddec2 + Fov_array * np.sin(theta)
            hp.visufunc.projscatter(racoord2, deccoord2, lonlat=True, marker='.', color=Colors[1])
            hp.projtext(targetCoord2[j].ra, targetCoord2[j].dec, str(j), lonlat=True, color=Colors[2])


def PointingPlottingGWCTA(filename, ID, outDir, SuggestedPointings, obspar):

    print()
    print('-------------------   PLOTTING SCHEDULE   --------------------')
    print()

    UseObs = obspar.obs_name
    FOV = obspar.FOV 
    # Mask table if necesary
    maskClean = (SuggestedPointings['ObsInfo'] == 'True')
    SuggestedPointingsC = SuggestedPointings[maskClean]
    SuggestedPointingsC.remove_column('ObsInfo')
   
    # Observatory
    if UseObs == 'South':
        observatory = CTASouthObservatory()
    else:
        observatory = CTANorthObservatory()

    skymap_OD = lf.read_sky_map(filename)
    prob = skymap_OD[0]
    npix = len(prob)
    nside = hp.npix2nside(npix)

    # ObservationTimearray, Coordinates, Probarray = LoadPointingsGW(SuggestedPointings)
    ObservationTimearray = SuggestedPointingsC["Time[UTC]"]
    Coordinates = SkyCoord(SuggestedPointingsC['RA[deg]'], SuggestedPointingsC['DEC[deg]'], frame='icrs',
                           unit=(u.deg, u.deg))
    Probarray = SuggestedPointingsC["PGW"]
    # making sure to have an array on which we can call the sum() function
    Probarray = np.atleast_1d(Probarray)

    print('----------   BUILDING A MAP   ----------')
    print('Total probability of map 1 that maximises PGW= {0:.5f}'.format(
        sum(Probarray)))

    converted_time = []
    for i, time in enumerate(ObservationTimearray):
        try:
            converted_time.append(datetime.datetime.strptime(
                time, '%Y-%m-%d %H:%M:%S.%f'))
        except ValueError:
            try:
                converted_time.append(datetime.datetime.strptime(
                    time, '%Y-%m-%d %H:%M:%S.%f '))
            except ValueError:
                converted_time.append(
                    datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S'))

    # PlotPointings(prob,cat,converted_time,Coordinates,sum(Probarray), nside, FOV, name, dirName, doPlot=True)

    t = 0.5 * np.pi - Coordinates[0].dec.rad
    p = Coordinates[0].ra.rad
    xyz = hp.ang2vec(t, p)

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))

    tt, pp = hp.pix2ang(nside, ipix_disc)
    ra2 = np.rad2deg(pp)
    dec2 = np.rad2deg(0.5 * np.pi - tt)

    skycoord = co.SkyCoord(ra2, dec2, frame='icrs', unit=(u.deg, u.deg))

    frame = co.AltAz(obstime=converted_time[0], location=observatory.location)
    altaz_all = skycoord.transform_to(frame)

    dirName = '%s/Pointing_Plotting_%s/%s' % (outDir, UseObs, ID)
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    # path = os.path.dirname(os.path.realpath(__file__)) + '/Pointing_Plotting'
    # if not os.path.exists(path):
    #    os.mkdir(path, 493)

    hp.mollview(prob, rot=[180, 0], coord='C', title="GW prob map (Equatorial) + %s %g  %s/%s/%s %s:%s:%s UTC" %
                (str(ID), sum(Probarray) * 100, converted_time[0].day, converted_time[0].month, converted_time[0].year,
                 converted_time[0].hour, converted_time[0].minute, converted_time[0].second))
    hp.graticule()
    # plt.show()
    theta = np.random.rand(400) * 360
    tarcoordra = np.empty(400)
    tarcoorddec = np.empty(400)
    Fov_array = np.empty(400)
    Fov_array.fill(FOV)

    for j in range(0, len(Coordinates.ra)):
        tarcoordra.fill(Coordinates[j].ra.deg)
        tarcoorddec.fill(Coordinates[j].dec.deg)
        racoord = tarcoordra + Fov_array * np.cos(theta)
        deccoord = tarcoorddec + Fov_array * np.sin(theta)
        hp.visufunc.projscatter(
            racoord, deccoord, lonlat=True, marker='.', color=Colors[1], coord='C')
        # plt.show()
    plt.savefig('%s/Pointings.png' % dirName)

    # dist = cat['Dist']
    # hp.visufunc.projscatter(cat['RAJ2000'][dist < 200], cat['DEJ2000'][dist < 200], lonlat=True, marker='.',color='g', linewidth=0.1, coord='C')


def PlotPointings_Pretty(filename, name, PointingsFile1, dirName, obspar, gal, centerMap = None, radiusMap = None, colorspar = None):
    try:
        ragal = gal['RAJ2000']
        decgal = gal['DEJ2000']
        galprob = gal['dp_dV']
        print("Plotting galaxies")
    except:
        print("No galaxies given")
    # Read the pointings file
    tpointingFile = PointingsFile1
    # tpointingFile = '/Users/mseglar/Documents/GitLab/lst_gwfollowup/output/bn180720598/PGWinFoV/RankingObservationTimes_Complete.txt'
    time = []
    try:
        time1, time2, ra, dec, pgw, pgal, Round, nametel, duration, fov = np.genfromtxt(tpointingFile, usecols=(
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9), skip_header=1, delimiter=' ', unpack=True, dtype='str')  # ra, dec in degrees
    except:
        try:
            time1, time2, ra, dec, pgw, Round, nametel, duration, fov = np.genfromtxt(tpointingFile, usecols=(
                0, 1, 2, 3, 4, 5, 6, 7, 8), skip_header=1, delimiter=' ', unpack=True, dtype='str')  # ra, dec in degrees
            pgal = pgw
        except:
            try:
                time1, time2, ra, dec, pgw, pgal, Round, nametel, duration, fov  = np.genfromtxt(tpointingFile, usecols=(
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9), skip_header=1, delimiter=' ', unpack=True, dtype='str')  # ra, dec in degrees
            except:
                time1, time2, ra, dec, pgw, Round, nametel, duration, fov  = np.genfromtxt(tpointingFile, usecols=(
                    0, 1, 2, 3, 4, 5, 6, 7, 8), skip_header=1, delimiter=' ', unpack=True, dtype='str')  # ra, dec in degrees
                pgal = pgw

    # print(time1, time2)

    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    fov = np.atleast_1d(fov)
    nametel = np.atleast_1d(nametel)
    ra = ra.astype(float)
    dec = dec.astype(float)
    pgw = pgw.astype(float)
    pgal = pgal.astype(float)
    fov = fov.astype(float)

    coordinates = SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg))

    # Prepare the color map
    cdict_coolheat = {
        'red':  ((0., 0., 0.), (0.25, 0., 0.), (0.5, 1., 1.), (0.75, 1.0, 1.0),  (1., 1., 1.)),
        'green':  ((0., 0., 0.), (0.25, 0., 0.), (0.5, 0., 0.), (0.75, 1.0, 1.0),  (1., 1., 1.)),
        'blue':  ((0., 0., 0.), (0.25, 1., 1.), (0.5, 0., 0.), (0.75, 0.0, 0.0),  (1., 1., 1.))
    }
    coolheat = colors.LinearSegmentedColormap('coolheat', cdict_coolheat, 1024)

    # center = SkyCoord(12, -25, unit='deg', frame='icrs')
    if centerMap == None:    
        center = SkyCoord(ra[0], dec[0], unit='deg', frame='icrs')
    else: 
        center = centerMap

    if radiusMap == None:    
        radius = '20 deg'
    else: 
        radius = radiusMap

    # center = SkyCoord(195, 15, unit='deg', frame='icrs')
    # center = SkyCoord(30, 20, unit='deg', frame='icrs')
    center_str = '%fd %fd' % (center.ra.deg, center.dec.deg)

    # start preparing figure and inset figure
    fig = plt.figure(figsize=(9, 6))

    ax = plt.axes([0.1, 0.1, 0.4, 0.4],
                  projection='astro degrees globe',
                  center=center_str)

    ax_inset = plt.axes([0.5, 0.1, 0.45, 0.45],
                        projection='astro degrees zoom',
                        center=center_str, radius=radius)

    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(True)
        ax_inset.coords[key].set_ticks_visible(True)

    ax.grid()
    ax_inset.grid()
    ax.mark_inset_axes(ax_inset)
    ax_inset.set_xlabel('RA (J2000)')
    ax_inset.set_ylabel('Dec (J2000)')

    # ax.connect_inset_axes(ax_inset, loc='upper left')
    # ax.connect_inset_axes(ax_inset, loc='lower left')

    ax.imshow_hpx(filename, cmap='cylon')

    try:
        vmin_value = np.min(galprob)
        vmax_value = np.max(galprob)
    except:
        vmin_value = 0.000001
        vmax_value = 0.00001

    try:
        norm = mcolors.Normalize(vmin=vmin_value, vmax=vmax_value)
        sc_inset = ax_inset.scatter(ragal, decgal, c = galprob, transform=ax_inset.get_transform('icrs'), alpha=0.5, s = 0.1, norm = norm)
        cbar_inset = plt.colorbar(sc_inset, ax=ax_inset)
        cbar_inset.set_label('Galaxy probability density')
    except:
        print('No galaxies given for plot 2')

    unique_obs_names = np.unique(nametel) 
    if colorspar == None:
        colors1 = plt.cm.get_cmap('tab10', len(unique_obs_names))
        obs_name_to_color = {name: colors1(i) for i, name in enumerate(unique_obs_names)}
    else: 
        colors1 = colorspar
        obs_name_to_color = {name: colors1[i] for i, name in enumerate(unique_obs_names)}

    # color_list = ['blue', 'green']
    # obs_name_to_color = {name: color_list[i % len(color_list)] for i, name in enumerate(unique_obs_names)}

    pos = ax_inset.imshow_hpx(filename, cmap='cylon')

    for i in range(0, len(ra)):
        obs_name = nametel[i]
        color1 = obs_name_to_color[obs_name]

        # ax_inset.scatter(ra[i], dec[i], s = fov[i], edgecolor=color1, facecolor="None", transform=ax_inset.get_transform('icrs'), alpha=1)
        circle_patch = Circle((ra[i], dec[i]), fov[i], edgecolor=color1, facecolor="None", transform=ax_inset.get_transform('icrs'), alpha=1)
        ax_inset.add_patch(circle_patch)
        # c = Circle((ra[i], dec[i]), fov_plot, edgecolor=color_map[name], facecolor="None", transform=ax_inset.get_transform('icrs'), alpha=1)

        legend_drawn_flag = True

    # Step 6: Create a legend with unique labels
    unique_labels = {}
    for obs_name, color in obs_name_to_color.items():
        unique_labels[color] = obs_name

    handles = []
    for color, label in unique_labels.items():
        handles.append(mpatches.Patch(color=color, label=label))
        # mpatches.Patch(color=color, label=label)

    plt.legend(handles=handles, loc='best')

    ax.mark_inset_axes(ax_inset)
    ax.connect_inset_axes(ax_inset, loc = 'lower left')
    ax.connect_inset_axes(ax_inset, loc = 'upper left')

    cbar = fig.colorbar(pos, ax=ax, location="left", fraction=0.046, pad=0.05)
    cbar.formatter.set_powerlimits((0, 0))
    # to get 10^3 instead of 1e3
    cbar.formatter.set_useMathText(True)
    cbar.set_label("Map probability density",  color='black', fontsize=9)
    plt.savefig("%s/Plot_PrettyMap_%s.png" % (dirName, name), dpi=300, bbox_inches='tight')
    plt.close()
