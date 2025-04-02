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
#                           Low-level tools to schedule tiled observations                     #
##################################################################################################


import datetime

#####################################################################
# Packages
import os

import astropy.coordinates as co
import ephem
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pytz
import six
import tables
from astropy import units as u
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord, get_body
from astropy.table import Table
from astropy.time import Time
from gdpyc import DustMap
from pytz import timezone
from six.moves import configparser
from skyfield import almanac
from skyfield.api import E, N, load, wgs84

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser

__all__ = [
    "Tools",
    "Observer",
    "NextWindowTools",
    "getdate",
    "UNIQSkymap_toNested",
    "get_lvk_uniq_maps",
    "NightDarkObservation",
    "NightDarkObservationwithGreyTime",
    "ZenithAngleCut",
    "ComputeProbability2D",
    "SubstractPointings2D",
    "TransformRADec",
    "LoadGalaxies",
    "LoadGalaxies_SteMgal",
    "VisibleAtTime",
    "FulfillsRequirement",
    "FulfillsRequirementGreyObservations",
    "ObtainHighestProbabilityCoordinates",
    "ComputeProbGalTargeted",
    "SubstractPointings",
    "SubstractGalaxiesCircle",
    "ComputePGalinFOV",
    "ModifyCatalogue",
    "ComputeProbPGALIntegrateFoV",
    "Get90RegionPixReduced",
    "Get90RegionPixGal",
    "IsSourceInside",
    "FillSummary",
    "GetSatelliteName",
    "GetSatelliteTime",
    "GetSatellitePositions",
    "OccultationCut",
    "GetEarthOccultedPix",
    "GetMoonOccultedPix",
    "GetSunOccultedPix",
    "SAA_Times",
    "PlotSpaceOcc",
    "GetBestNSIDE",
    "FillSummary",
]


class Tools:
    """
    class with different visibility check functions and other setting and rising of the sun and the moon functions
    """

    @classmethod
    def IsDarkness(cls, obsTime, obspar):
        sunAlt = Tools.SunAlt(obsTime, obspar)
        moonAlt = Tools.MoonAlt(obsTime, obspar)
        SunDown = obspar.sunDown
        MoonDown = obspar.moonDown

        if sunAlt > SunDown:
            return False
        if moonAlt > MoonDown:
            return False

        return True

    @classmethod
    def IsGreyness(cls, obsTime, obspar):
        # SUN altitude
        sunAlt = Tools.SunAlt(obsTime, obspar)
        # MOON altitude
        moonAlt = Tools.MoonAlt(obsTime, obspar)
        # MOON azimuth
        # moonAz = Tools.MoonAz(obsTime,obsSite)
        # MOON phase
        moonPhase = Tools.MoonPhase(obsTime, obspar)
        SunDown = obspar.sunDown
        MoonDown = obspar.moonDown
        MoonGrey = obspar.moonGrey
        MoonPhase = obspar.moonPhase
        if sunAlt > SunDown:
            return False
        if moonAlt > MoonGrey:
            return False
        if moonPhase > MoonPhase and moonAlt > MoonDown:
            return False
        return True

    @classmethod
    def MoonPhase(cls, obsTime, obspar):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        moon.compute(obs)

        # print("Phase of the moon = %s percent" % moon.phase)

        return moon.phase

    @classmethod
    def SunAlt(cls, obsTime, obspar):
        sun = get_body("sun", Time(obsTime, scale="utc")).transform_to(
            AltAz(obstime=Time(obsTime, scale="utc"), location=obspar.location)
        )
        # print(Time(obsTime),obsSite.location)
        # print(get_sun(Time(obsTime)))
        # print(sun.alt/u.deg)
        return sun.alt / u.deg

    @classmethod
    def MoonAlt(cls, obsTime, obspar):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        # print(obs)
        moon.compute(obs)
        # print('Altitude of the moon = ',moon.alt * 180. / np.pi)
        return moon.alt * 180.0 / np.pi

    @classmethod
    def NextSunrise(cls, obsTime, obspar):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        # obs.horizon = obspar.horizonSun
        obs.horizon = Angle(obspar.sunDown, u.deg).to_string(unit=u.degree, sep=":")
        sun.compute(obs)
        nextSunrise = (
            obs.next_rising(sun, use_center=True).datetime().replace(tzinfo=pytz.utc)
        )
        return nextSunrise

    @classmethod
    def PreviousSunrise(cls, obsTime, obspar):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        # obs.horizon = obspar.horizonSun
        obs.horizon = Angle(obspar.sunDown, u.deg).to_string(unit=u.degree, sep=":")
        sun.compute(obs)
        previousSunrise = (
            obs.previous_rising(sun, use_center=True)
            .datetime()
            .replace(tzinfo=pytz.utc)
        )
        return previousSunrise

    @classmethod
    def NextSunset(cls, obsTime, obspar):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        # obs.horizon = obspar.horizonSun
        obs.horizon = Angle(obspar.sunDown, u.deg).to_string(unit=u.degree, sep=":")
        sun.compute(obs)
        nextSunset = (
            obs.next_setting(sun, use_center=True).datetime().replace(tzinfo=pytz.utc)
        )
        return nextSunset

    @classmethod
    def PreviousMoonset(cls, obsTime, obspar):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        # obs.horizon = obspar.HorizonMoon
        obs.horizon = Angle(obspar.moonDown, u.deg).to_string(unit=u.degree, sep=":")
        moon.compute()
        previousMoonset = (
            obs.previous_setting(moon, use_center=True)
            .datetime()
            .replace(tzinfo=pytz.utc)
        )
        return previousMoonset

    @classmethod
    def NextMoonset(cls, obsTime, obspar):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        obs.horizon = Angle(obspar.moonDown, u.deg).to_string(unit=u.degree, sep=":")
        moon.compute()
        nextMoonset = (
            obs.next_setting(moon, use_center=True).datetime().replace(tzinfo=pytz.utc)
        )
        return nextMoonset

    @classmethod
    def TrustingDarknessSun(cls, obsTime, obspar):
        DarkObsTime = obsTime
        referencetime = obsTime
        while not Tools.IsDarkness(DarkObsTime, obspar) and (
            (
                DarkObsTime.hour >= referencetime.hour
                and DarkObsTime.day == referencetime.day
            )
            or (
                DarkObsTime.hour <= Tools.NextSunrise(referencetime, obspar).hour
                and DarkObsTime.day == Tools.NextSunrise(referencetime, obspar).day
            )
        ):
            DarkObsTime = DarkObsTime + datetime.timedelta(minutes=1)
        return DarkObsTime

    @classmethod
    def TrustingGreynessSun(cls, obsTime, obspar):
        GreyObsTime = obsTime
        referencetime = obsTime
        while not Tools.IsGreyness(GreyObsTime, obspar) and (
            (
                GreyObsTime.hour >= referencetime.hour
                and GreyObsTime.day == referencetime.day
            )
            or (
                GreyObsTime.hour <= Tools.NextSunrise(referencetime, obspar).hour
                and GreyObsTime.day == Tools.NextSunrise(referencetime, obspar).day
            )
        ):
            GreyObsTime = GreyObsTime + datetime.timedelta(minutes=1)
            # print(Tools.IsDarkness(DarkObsTime,obsSite))
        return GreyObsTime

    @classmethod
    def UTCtoNamibia(cls, UTCtime):
        TimezonesDifference = datetime.timedelta(hours=2)
        NamibianTime = UTCtime + TimezonesDifference
        return NamibianTime

    @classmethod
    def CheckWindow(cls, time, obspar):
        MinimalWindowDuration = datetime.timedelta(minutes=obspar.minDuration)
        if (Tools.IsDarkness(time, obspar) is True) and (
            Tools.IsDarkness(time + MinimalWindowDuration, obspar) is True
        ):
            Observe = True
        else:
            Observe = False

        return Observe

    @classmethod
    def CheckWindowGrey(cls, time, obspar):
        MinimalWindowDuration = datetime.timedelta(minutes=obspar.minDuration)
        if (Tools.IsGreyness(time, obspar) is True) and (
            Tools.IsGreyness(time + MinimalWindowDuration, obspar) is True
        ):
            Observe = True
        else:
            Observe = False
        return Observe

    @classmethod
    def GalacticPlaneBorder(cls, coords):
        lon = coords.galactic.l.value  # x-coordinate
        lat = coords.galactic.b.value  # y-coordinate
        # print(lon)
        # print(lat)
        YouAreInside = False
        n = 20
        if lat <= 10 and lat >= 0 and lon <= 130:
            n = lat - (1.0 / 13) * lon
        elif lat <= 10 and lat >= 0 and lon >= 240:
            n = lat - (1.0 / 12) * lon + 20
        elif lat >= -10 and lat <= 0 and lon <= 130:
            n = lat + (1.0 / 13) * lon
        elif lat >= -10 and lat <= 0 and lon >= 240:
            n = lat + (1.0 / 12) * lon - 20
        if np.absolute(n) <= 10:
            YouAreInside = True
        return YouAreInside

    @classmethod
    def GetGalacticExtinction(cls, coords, dustmap="SFD", filters="SDSS_r"):
        # Extinction = DustMap.ebv(coords)
        extinction = DustMap.extinction(coords, dustmap="SFD", filters="SDSS_r")
        # GasMap.plot_map('HI4PI')
        return extinction

    @classmethod
    def is_in_saa(cls, latitude, longitude):
        saa_lat_min = -40.0  # Minimum latitude for the SAA
        saa_lat_max = 0.0  # Maximum latitude for the SAA
        saa_lon_min = -50.0  # Minimum longitude for the SAA
        saa_lon_max = -30.0  # Maximum longitude for the SAA

        # Check if the satellite's position falls within the SAA region
        if (saa_lat_min <= latitude <= saa_lat_max) and (
            saa_lon_min <= longitude <= saa_lon_max
        ):
            return True
        else:
            return False

    @classmethod
    def query_square(nside, center, side_length_rad):
        # Convert side length to radians

        # Calculate corner offsets from the center point (assuming a small angle approximation)
        dx = side_length_rad / np.sqrt(2)

        # Get four corners in the form of xyz offsets
        corners = [
            center + np.array([dx, dx, 0]),
            center + np.array([-dx, dx, 0]),
            center + np.array([dx, -dx, 0]),
            center + np.array([-dx, -dx, 0]),
        ]

        # Query discs at each corner point
        pixels = set()
        for corner in corners:
            pix_ids = hp.query_disc(nside, corner, side_length_rad, inclusive=True)
            pixels.update(pix_ids)

        return list(pixels)

    @classmethod
    def hexagon_vertices(center, radius):
        """Calculate hexagon vertices around a center point on the sphere."""
        theta_c, phi_c = center
        vertices = []

        # Angle step for each vertex (60 degrees apart)
        angle_step = 2 * np.pi / 6

        for i in range(6):
            angle = i * angle_step
            theta_v = theta_c + radius * np.cos(angle)
            phi_v = phi_c + radius * np.sin(angle)
            vertices.append((theta_v, phi_v))

        return vertices


class Observer:
    """Class to store information and handle operation related to the observatory used for the observations."""

    def __init__(
        self,
        longitude,
        latitude,
        elevation,
        run_duration,
        minimal_run_duration,
        max_sun_altitude,
        max_moon_altitude,
        max_moon_phase,
    ):
        """
        Initialize the class with all the needed parameters.

        Args:
        longitude (float): The longitude of the observatory (degree).
        latitude (float): The latitude of the observatory (degree).
        elevation (float): The elevation of the observatory (m).
        run_duration (datetime.timedelta): The duration of each run.
        minimal_run_duration (datetime.timedelta): The minimum duration of each run.
        max_sun_altitude (float): The maximum altitude of the sun (degree).
        max_moon_altitude (float): The maximum altitude of the moon (degree).
        max_moon_phase (float): The maximum phase of the moon (illumination fraction).
        """

        self.eph = load("de440s.bsp")
        self.observatory_location = wgs84.latlon(latitude * N, longitude * E, elevation)

        self.max_sun_altitude = max_sun_altitude
        self.max_moon_altitude = max_moon_altitude
        self.max_moon_phase = max_moon_phase

        self.run_duration = run_duration
        self.minimal_run_duration = minimal_run_duration

        # Setup time converter for skyfield
        self.timescale_converter = load.timescale()

    def get_time_window(self, start_time, nb_observation_night):
        """
        Calculate the time window for observations.

        Args:
        start_time (datetime): Earliest start time to start observations, if no timezone, assume UTC
        nb_observation_night (int): The number of observation nights.

        Returns:
        list: The start times for each run within the valid time range.
        """

        # Compute time interval
        if start_time.tzinfo is None:
            zone = timezone("UTC")
            start_time = zone.localize(start_time)
        delta_time_obs = datetime.timedelta(days=nb_observation_night + 1)
        stop_time = start_time + delta_time_obs

        time_interval_sun = self.get_sun_constraint_time_interval(
            start_time, stop_time, nb_observation_night
        )
        time_interval_moon = self.get_moon_constraint_time_interval(
            start_time, stop_time
        )
        valid_time_interval = self.compute_interval_intersection(
            time_interval_sun, time_interval_moon
        )
        return self.compute_run_start_time(valid_time_interval)

    def compute_interval_intersection(self, time_range_1, time_range_2):
        """
        Compute the intersection of two time ranges.

        Args:
        time_range_1 (list): The first time range.
        time_range_2 (list): The second time range.

        Returns:
        list: The intersection of the two time ranges.
        """
        # Initialisation
        i = j = 0
        n = len(time_range_1)
        m = len(time_range_2)
        intersection_time_intervals = []

        # Loop through all intervals unless one
        # of the interval gets exhausted
        while i < n and j < m:
            # Determine if the intersection of the two currently selected intervals is valid
            left = max(time_range_1[i][0], time_range_2[j][0])
            right = min(time_range_1[i][1], time_range_2[j][1])
            if left <= right:
                intersection_time_intervals.append([left, right])

            # Move to the next interval for the youngest one
            if time_range_1[i][1] < time_range_2[j][1]:
                i += 1
            else:
                j += 1

        return intersection_time_intervals

    def compute_run_start_time(self, valid_time_range):
        """
        Compute the start times for each run within a valid time range.

        Args:
        valid_time_range (list): The valid time range.

        Returns:
        list: The start times for each run.
        """
        run_start_time = []
        for i in range(len(valid_time_range)):
            observation_time_available = valid_time_range[i][1] - valid_time_range[i][0]
            nb_observation_run = int(
                np.rint(observation_time_available // self.run_duration)
            )
            remaining_observation_time = observation_time_available % self.run_duration
            if remaining_observation_time > self.minimal_run_duration:
                nb_observation_run += 1
            run_start_time += list(
                valid_time_range[i][0]
                + np.arange(nb_observation_run) * self.run_duration
            )
        return run_start_time

    def get_sun_constraint_time_interval(
        self, start_time, stop_time, nb_observation_night
    ):
        """
        Get the time interval for the sun constraint.

        Args:
        start_time (datetime): The start time of observations.
        stop_time (datetime): The stop time of observations.
        nb_observation_night (int): The number of observation nights.

        Returns:
        list: The time interval for the sun constraint.
        """
        rise_time, set_time = self.get_risings_and_settings(
            "sun", self.max_sun_altitude, start_time, stop_time
        )

        time_interval_sun = []
        for i in range(min(nb_observation_night, len(set_time))):
            if set_time[i] < rise_time[i]:
                time_interval_sun.append([set_time[i], rise_time[i]])
            elif i < (len(set_time) - 1):
                time_interval_sun.append([set_time[i], rise_time[i + 1]])
            else:
                continue

        return time_interval_sun

    def get_moon_constraint_time_interval(self, start_time, stop_time):
        """
        Get the time interval for the moon constraint.

        Args:
        start_time (datetime): The start time of observations.
        stop_time (datetime): The stop time of observations.

        Returns:
        list: The time interval for the moon constraint.
        """
        rise_time, set_time = self.get_risings_and_settings(
            "moon", self.max_moon_altitude, start_time, stop_time
        )

        # Initialise data for time intervals
        time_interval_moon = []
        tmp_time_interval_moon = []
        i, j = 0, 0
        time_cursor_start_window = start_time
        time_cursor_end_window = start_time

        # Loop over the time to determine valid time window
        while time_cursor_end_window < stop_time:
            # Determine interval range
            time_cursor_start_window = min(rise_time[i], set_time[j])
            time_cursor_end_window = max(rise_time[i], set_time[j])

            # Determine the interval is valid
            valid_interval = False
            if rise_time[i] > set_time[j]:
                valid_interval = True
                j += 1
            else:
                time_middle_window = rise_time[i] + (set_time[j] - rise_time[i]) / 2.0
                valid_interval = (
                    self.get_moon_phase(time_middle_window) < self.max_moon_phase
                )
                i += 1

            # Apply action on the wider interval based on the validity results
            if valid_interval and len(tmp_time_interval_moon) == 0:
                tmp_time_interval_moon.append(time_cursor_start_window)
            elif not valid_interval and len(tmp_time_interval_moon) == 1:
                tmp_time_interval_moon.append(time_cursor_start_window)
                time_interval_moon.append(tmp_time_interval_moon)
                tmp_time_interval_moon = []

        # If the current time interval is open after iteration, close it
        if len(tmp_time_interval_moon) == 1:
            tmp_time_interval_moon.append(time_cursor_end_window)
            time_interval_moon.append(tmp_time_interval_moon)

        return time_interval_moon

    def get_moon_phase(self, observation_time):
        """
        Get the moon phase at a given observation time.

        Args:
        observation_time (datetime): The time of observation.

        Returns:
        float: The moon phase at the given observation time.
        """
        sun, moon, earth = self.eph["sun"], self.eph["moon"], self.eph["earth"]

        return (
            earth.at(self.timescale_converter.from_datetime(observation_time))
            .observe(moon)
            .apparent()
            .fraction_illuminated(sun)
        )

    def get_risings_and_settings(self, celestial_body, horizon, start_time, stop_time):
        """
        Get the rise and set times of a celestial body within a given time range.

        Args:
        celestial_body (str): The celestial body.
        horizon (float): The horizon to consider as rised or set (degree).
        start_time (datetime): The start time.
        stop_time (datetime): The stop time.

        Returns:
        list, list: The rise and set times of the celestial body.
        """
        f = almanac.risings_and_settings(
            self.eph,
            self.eph[celestial_body],
            self.observatory_location,
            horizon_degrees=horizon,
        )
        time, rising_indicator = almanac.find_discrete(
            self.timescale_converter.from_datetime(start_time),
            self.timescale_converter.from_datetime(stop_time),
            f,
        )
        if len(time) == 0:
            alt, az, distance = (
                (self.observatory_location + self.eph["earth"])
                .at(self.timescale_converter.from_datetime(start_time))
                .observe(self.eph[celestial_body])
                .apparent()
                .altaz()
            )
            if alt.degrees < horizon:
                rise_time = [
                    stop_time,
                ]
                set_time = [
                    start_time,
                ]
            else:
                rise_time = [
                    start_time,
                ]
                set_time = [
                    stop_time,
                ]
        else:
            rise_time = list(time[rising_indicator == 1].utc_datetime())
            set_time = list(time[rising_indicator == 0].utc_datetime())

            # Set the start time and end time as either rise of set time to fully cover the time range
            if rising_indicator[0] == 0:
                rise_time = [
                    start_time,
                ] + rise_time
            else:
                set_time = [
                    start_time,
                ] + set_time
            if rising_indicator[-1] == 1:
                set_time = set_time + [
                    stop_time,
                ]
            else:
                rise_time = rise_time + [
                    stop_time,
                ]

        return rise_time, set_time


######################################################

# Functions related to the Skymap handling

######################################################


def getdate(x):
    """
    Bottom-level function that takes a date and prints it in ISO format

    :param x: the date to be formatted
    :type x: datetime

    :return: none
    rtype: none
    """

    if isinstance(x, datetime.datetime):
        return x
    elif isinstance(x, str):
        return datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
    else:
        print("ERROR: something is wrong with the format of the date: ", x)
        return None


def UNIQSkymap_toNested(skymap_fname):
    """
    Bottom-level function that takes a skymap and computes from it the uniq map

    :param skymap_fname: Healpix skymap
    :type skymap_fname: Table

    :return: healpix_skymaps_dict
    rtype: dict
    """

    sky_tab = Table.read(skymap_fname)
    healpix_skymaps_dict = get_lvk_uniq_maps(sky_tab, "max")
    # prob = healpix_skymaps_dict['PROBDENSITY']
    return healpix_skymaps_dict


def get_lvk_uniq_maps(sky_map, Order, map_names="all"):
    un_inds = sky_map["UNIQ"]

    order = (np.log2(un_inds / 4).astype(int) / 2).astype(int)
    inds = (un_inds - 4 * (np.power(4, order))).astype(int)

    if Order == "max":
        Order = np.max(order)

    Nside = int(2**Order)
    Npix = hp.nside2npix(Nside)

    if map_names == "all":
        keys = ["PROB", "DISTMU", "DISTSIGMA", "DISTNORM"]
    else:
        keys = map_names
    maps = {}

    for k in keys:
        maps[k] = np.zeros(Npix)

    # print np.min(order), np.max(order)

    for ii in range(np.max(order), np.min(order) - 1, -1):

        nside = 2**ii
        npix = hp.nside2npix(nside)
        bl = order == ii

        for k in maps.keys():
            a = hp.UNSEEN * np.ones(npix)
            if k == "PROB":
                a[inds[bl]] = sky_map["PROBDENSITY"][bl]

            else:
                a[inds[bl]] = sky_map[k][bl]
            if ii == Order:
                bl_ = a != hp.UNSEEN
                maps[k][bl_] += a[bl_]
                del a
            else:
                a_ = hp.ud_grade(
                    a, nside_out=Nside, order_in="Nested", order_out="Nested"
                )
                bl_ = a_ != hp.UNSEEN
                maps[k][bl_] += a_[bl_]
                del a, a_

    maps["PROB"] = (
        maps["PROB"] * (np.pi / 180) ** 2 * hp.nside2pixarea(Nside, degrees=True)
    )
    # print('Total probability is:', maps['PROB'].sum())
    return maps


def NightDarkObservation(time, obspar):
    """
    Function that searches for an array of observation times that fulfilled darkness condition and window

    """
    obs = Observer(
        longitude=obspar.lon.to_value(u.deg),
        latitude=obspar.lat.to_value(u.deg),
        elevation=obspar.height.to_value(u.m),
        run_duration=datetime.timedelta(minutes=obspar.duration),
        minimal_run_duration=datetime.timedelta(minutes=obspar.minDuration),
        max_sun_altitude=obspar.sunDown,
        max_moon_altitude=obspar.moonDown,
        max_moon_phase=-1.0,
    )
    return obs.get_time_window(start_time=time, nb_observation_night=obspar.maxNights)


def NightDarkObservationwithGreyTime(time, obspar):
    """
    Function that searches for an array of observation times that fulfilled darkness condition and window

    """
    obs = Observer(
        longitude=obspar.lon.to_value(u.deg),
        latitude=obspar.lat.to_value(u.deg),
        elevation=obspar.height.to_value(u.m),
        run_duration=datetime.timedelta(minutes=obspar.duration),
        minimal_run_duration=datetime.timedelta(minutes=obspar.minDuration),
        max_sun_altitude=obspar.sunDown,
        max_moon_altitude=obspar.moonGrey,
        max_moon_phase=obspar.moonPhase / 100.0,
    )
    return obs.get_time_window(start_time=time, nb_observation_night=obspar.maxNights)


def ZenithAngleCut(
    prob,
    nside,
    time,
    minProbcut,
    maxZenith,
    observatory,
    minMoonSourceSeparation,
    useGreytime,
):
    """
    Mask in the pixels with zenith angle larger than 45
    """
    # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)
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
    # hp.mollview(mzenith)
    # plt.show()
    # plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_%g.png")

    yprob = ma.masked_array(pprob, mzenith.mask)
    # hp.mollview(yprob)
    # plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_prob_%g.png")

    # print('Integrated probability of the masked map', np.sum(yprob))

    if np.sum(yprob) < minProbcut:
        ObsBool = False
    else:
        ObsBool = True

    if useGreytime and ObsBool:
        # Get Alt/Az of the Moon
        moonaltazs = get_body("moon", Time(time, scale="utc")).transform_to(
            AltAz(obstime=Time(time, scale="utc"), location=observatory)
        )
        separations = altaz_map.separation(moonaltazs)
        mask_moonDistance = np.zeros(hp.nside2npix(nside), dtype=bool)
        mask_moonDistance[separations < minMoonSourceSeparation * u.deg] = 1
        mzenith = hp.ma(pprob)
        mzenith.mask = mask_moonDistance
        yprob = ma.masked_array(pprob, mzenith.mask)
        # hp.mollview(pprob)
        # hp.mollview(yprob)
        # plt.show()
        if np.sum(yprob) < minProbcut:
            ObsBool = False
        else:
            ObsBool = True
        # print('Integrated probability of the masked map', np.sum(yprob))
        # hp.mollview(mzenith)
        # plt.show()
        # Get the mask that does a radius around, of 30 degs
        # Plot to check
        # Return a bool if there is any observable region

    return ObsBool, yprob


def GetSatelliteName(satellitename, stationsurl):
    stations_url = stationsurl
    satellites = load.tle_file(stations_url)
    # print('Loaded', len(satellites), 'satellites')
    by_name = {sat.name: sat for sat in satellites}
    satellite_name = by_name.get(satellitename)
    return satellite_name


def GetSatelliteTime(satellite_name, t):
    dt = t  # Get the underlying datetime object
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second + dt.microsecond / 1e6  # Include microseconds
    skyfield_time = load.timescale().utc(year, month, day, hour, minute, second)
    return skyfield_time


def GetSatellitePositions(satellite_name, t):
    geocentric = satellite_name.at(t)
    # Print position in latitude, longitude, altitude
    subpoint = geocentric.subpoint()

    # Get satellite's current position in astronomical units (AU)
    satellite_position = geocentric.position.km
    satellite_location = EarthLocation(
        lat=subpoint.latitude.degrees,
        lon=subpoint.longitude.degrees,
        height=subpoint.elevation.m,
    )
    return satellite_position, satellite_location


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


def OccultationCut(
    prob, nside, time, minProbcut, satellite_position, observatory, sun_sep, moon_sep
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
        nside, time, 6371, 1, satellite_position, observatory
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
):
    SatTimes = []
    i = 0
    saa = []
    while current_time <= start_time + datetime.timedelta(minutes=duration):
        SatelliteTime = GetSatelliteTime(SatelliteName, current_time)
        satellite_position, satellite_location = GetSatellitePositions(
            SatelliteName, SatelliteTime
        )

        if Tools.is_in_saa(satellite_location.lat.deg, satellite_location.lon.deg):
            saa.append(True)
        else:
            saa.append(False)

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
        plt.title("SAA Times")
        plt.xlabel("Time")
        plt.ylabel("SAA Status")
        plt.ylim(-0.1, 1.1)  # Set limits to make binary values clear
        plt.legend()
        plt.grid(True)
        plt.savefig("%s/SAA_Times.png" % (path))
    return SatTimes, saa


def GetBestNSIDE(ReducedNSIDE, HRnside, fov):

    if isinstance(HRnside, int) and HRnside > 0 and (HRnside & (HRnside - 1)) == 0:
        max_nside = HRnside
    else:
        max_nside = 512

    if (
        isinstance(ReducedNSIDE, int)
        and ReducedNSIDE > 0
        and (ReducedNSIDE & (ReducedNSIDE - 1)) == 0
    ):
        best_nside = ReducedNSIDE
        print("The NSIDE is already given. No optimization...")

    else:
        nside_values = [2**i for i in range(1, 13)]  # From NSIDE=2 to NSIDE=4096
        nside_values = [nside for nside in nside_values if nside <= max_nside]
        pixel_sizes = {nside: (180.0 / (np.sqrt(3) * nside)) for nside in nside_values}
        valid_nsides = [nside for nside, size in pixel_sizes.items() if size <= fov]

        if not valid_nsides:
            best_nside = max_nside  # Default to max_nside if no valid NSIDE is found
        else:
            best_nside = max(valid_nsides)  # Choose the best NSIDE in range

        print("NO REDUCED NSIDE GIVEN. Optimizing...")
        print(
            f"Best NSIDE for FoV of {fov}° (Min NSIDE {ReducedNSIDE}, Max NSIDE {max_nside}): {best_nside} (Pixel Size ≈ {pixel_sizes[best_nside]:.3f}°)"
        )

    print("best_nside", best_nside)
    return max_nside, best_nside


def ComputeProbability2D(
    obspar,
    prob,
    highres,
    radecs,
    time,
    ipixlist,
    ipixlistHR,
    counter,
    dirName,
    ipixlistOcc=None,
):
    """
    Compute probability in 2D by taking the highest probability in FoV value
    """

    reducedNside = obspar.reducedNside
    HRnside = obspar.HRnside
    minProbcut = obspar.minProbcut
    observatory = obspar.location
    maxZenith = obspar.maxZenith
    radius = obspar.FOV
    useGreytime = obspar.useGreytime
    plot = obspar.doPlot

    frame = co.AltAz(obstime=time, location=observatory)
    thisaltaz = radecs.transform_to(frame)

    if useGreytime:
        moonaltazs = get_body("moon", Time(time, scale="utc")).transform_to(
            AltAz(obstime=Time(time, scale="utc"), location=observatory)
        )
        # Zenith and Moon angular distance mask
        pix_ra = radecs.ra.value[
            (thisaltaz.alt.value > 90 - maxZenith)
            & (thisaltaz.separation(moonaltazs) > (90 - maxZenith) * u.deg)
        ]
        pix_dec = radecs.dec.value[
            (thisaltaz.alt.value > 90 - maxZenith)
            & (thisaltaz.separation(moonaltazs) > (90 - maxZenith) * u.deg)
        ]

    else:
        # Zenith angle mask
        pix_ra = radecs.ra.value[(thisaltaz.alt.value > 90 - maxZenith)]
        pix_dec = radecs.dec.value[thisaltaz.alt.value > 90 - maxZenith]

    phipix = np.deg2rad(pix_ra)
    thetapix = 0.5 * np.pi - np.deg2rad(pix_dec)

    ipix = hp.ang2pix(reducedNside, thetapix, phipix)

    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)

    cat_pix = Table(
        [ipix, pix_ra, pix_dec, dp_Pix_Fov],
        names=("PIX", "PIXRA", "PIXDEC", "PIXFOVPROB"),
    )

    dp_dV_FOV = []

    xyzpix = hp.ang2vec(thetapix, phipix)

    # Grid-scheme and the connection between HR and LR
    for i in range(0, len(cat_pix)):
        # Pixels associated to a disk of radius centered in xyzpix[i] for HR NSIDE
        ipix_discfull = hp.query_disc(HRnside, xyzpix[i], np.deg2rad(radius))
        if len(ipixlistHR) == 0:
            # No mask needed
            HRprob = highres[ipix_discfull].sum()
        else:
            # Mask the ipix_discfull with the pixels that are already observed. I think the problem is here
            maskComputeProb = np.isin(ipix_discfull, ipixlistHR, invert=True)
            # Obtain list of pixel ID after the mask what has been observed already
            m_ipix_discfull = ma.compressed(
                ma.masked_array(ipix_discfull, mask=np.logical_not(maskComputeProb))
            )
            HRprob = highres[m_ipix_discfull].sum()
            # HRprob = 0
            # for j in ipix_discfullNotCovered:
            #    HRprob = HRprob+highres[j]
            # print('Length of list of pixels:', m_ipix_discfull, 'vs', ipix_discfull, 'vs', ipixlistHR)
            # print('Comparison to see if mask is considered: ',HRprob, 'vs',highres[ipix_discfull].sum())
        dp_dV_FOV.append(HRprob)
    cat_pix["PIXFOVPROB"] = dp_dV_FOV

    # Mask already observed pixels
    mask = np.isin(cat_pix["PIX"], ipixlist, invert=True)
    if all(np.isin(cat_pix["PIX"], ipixlist, invert=False)):
        maskcat_pix = cat_pix
    else:
        maskcat_pix = cat_pix[mask]
    # Sort table
    sortcat1 = maskcat_pix[np.flipud(np.argsort(maskcat_pix["PIXFOVPROB"]))]

    # Mask occulted pixels for this round without affecting ipixlist and ipixlistHR
    mask2 = np.isin(sortcat1["PIX"], ipixlistOcc, invert=True)
    if all(np.isin(sortcat1["PIX"], ipixlistOcc, invert=False)):
        sortcat2 = sortcat1[mask2]
    else:
        sortcat2 = sortcat1[mask2]
    # Sort table
    sortcat = sortcat2[np.flipud(np.argsort(sortcat2["PIXFOVPROB"]))]

    # Chose highest
    targetCoord = co.SkyCoord(
        sortcat["PIXRA"][:1], sortcat["PIXDEC"][:1], frame="fk5", unit=(u.deg, u.deg)
    )

    P_GW = sortcat["PIXFOVPROB"][:1]

    if P_GW >= minProbcut:
        phip = float(np.deg2rad(targetCoord.ra.deg))
        thetap = float(0.5 * np.pi - np.deg2rad(targetCoord.dec.deg))
        xyz = hp.ang2vec(thetap, phip)

        ipixlistHR.extend(hp.query_disc(HRnside, xyz, np.deg2rad(radius)))
        ipix_disc = hp.query_disc(reducedNside, xyz, np.deg2rad(radius))
        ipixlist.extend(ipix_disc)

        ##################################
        # PLOT THE RESULTS
        if plot:
            path = dirName + "/EvolutionPlot"
            if not os.path.exists(path):
                os.mkdir(path, 493)
            # nside = 1024

            # hp.mollview(highres,title="With FoV circle")

            hp.mollview(prob, title=str(time))

            hp.graticule()

            ipix_discplot = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
            tt, pp = hp.pix2ang(HRnside, ipix_discplot)
            ra2 = np.rad2deg(pp)
            dec2 = np.rad2deg(0.5 * np.pi - tt)
            skycoord = co.SkyCoord(ra2, dec2, frame="fk5", unit=(u.deg, u.deg))

            # hp.visufunc.projplot(skycoord.ra, skycoord.dec, 'y.', lonlat=True, coord="C")
            # plt.show()
            # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

            hp.visufunc.projplot(
                sortcat["PIXRA"][:1],
                sortcat["PIXDEC"][:1],
                "r.",
                lonlat=True,
                coord="C",
            )
            MaxCoord = SkyCoord(
                sortcat["PIXRA"][:1],
                sortcat["PIXDEC"][:1],
                frame="fk5",
                unit=(u.deg, u.deg),
            )
            separations = skycoord.separation(MaxCoord)
            tempmask = separations < (radius + 0.05 * radius) * u.deg
            tempmask2 = separations > (radius - 0.05 * radius) * u.deg
            hp.visufunc.projplot(
                skycoord[tempmask & tempmask2].ra,
                skycoord[tempmask & tempmask2].dec,
                "r.",
                lonlat=True,
                coord="C",
                linewidth=0.1,
            )
            hp.visufunc.projplot(
                skycoord[tempmask & tempmask2].ra,
                skycoord[tempmask & tempmask2].dec,
                "r.",
                lonlat=True,
                coord="C",
                linewidth=0.1,
            )

            if ipixlistOcc is not None:
                try:
                    tt, pp = hp.pix2ang(reducedNside, ipixlistOcc)
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
                    print("No occulted pixel")

            plt.savefig("%s/Zoom_Pointing_%g.png" % (path, counter))
            # for i in range(0,1):
            #    altcoord.fill(90-(maxZenith-5*i))
            #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
            #    RandomCoord_radec = RandomCoord.transform_to('fk5')
            #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
            # plt.show()
            # plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))
            plt.close()

    return P_GW, targetCoord, ipixlist, ipixlistHR


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
):

    xyzpix1 = hp.pix2vec(reducedNside, newpix)
    xyzpix = np.column_stack(xyzpix1)

    dp_dV_FOV = []

    for i in range(0, len(newpix)):
        ipix_discfull = hp.query_disc(HRnside, xyzpix[i], np.deg2rad(radius))
        HRprob = highres[ipix_discfull].sum()
        dp_dV_FOV.append(HRprob)

    theta, phi = hp.pix2ang(reducedNside, newpix)
    ra = np.degrees(phi)  # RA in degrees
    dec = 90 - np.degrees(theta)

    cat_pix = Table(
        [newpix, ra, dec, dp_dV_FOV], names=("PIX", "PIXRA", "PIXDEC", "PIXFOVPROB")
    )

    cat_pix["PIXFOVPROB"] = dp_dV_FOV

    sortcat1 = cat_pix[np.flipud(np.argsort(cat_pix["PIXFOVPROB"]))]
    first_values = sortcat1[:maxRuns]

    # mapsize = 200
    # centerRA = 314
    # centerDEC = 10

    if doPlot:
        path = dirName + "/GridPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)

        # mpl.rcParams.update({'font.size':14})
        hp.mollview(prob)
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
            print("No pcculted pix")

        hp.visufunc.projplot(
            first_values["PIXRA"],
            first_values["PIXDEC"],
            "b.",
            lonlat=True,
            coord="C",
            linewidth=0.1,
        )
        plt.savefig("%s/Occ_Pointing.png" % (path))
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
    UsePix,
    maxRuns,
    doPlot,
    dirName,
    reducedNside,
    Occultedpixels,
):
    SelectedGals = galpix
    dp_dV_FOV = []
    for element in range(0, len(SelectedGals)):
        if element < len(SelectedGals):
            dp_dV_FOV.append(
                ComputePGalinFOV(
                    prob,
                    cat,
                    SelectedGals[element],
                    FOV,
                    totaldPdV,
                    reducedNside,
                    UsePix=True,
                )
            )
        else:
            dp_dV_FOV.append(0)

    print("Length of SelectedGals:", len(SelectedGals))
    print("Length of dp_dV_FOV[:lengthSG]:", len(dp_dV_FOV))

    cat_pix = Table(
        [SelectedGals.ra.deg, SelectedGals.dec.deg, dp_dV_FOV],
        names=("PIXRA", "PIXDEC", "PIXFOVPROB"),
    )

    sortcat = cat_pix[np.flipud(np.argsort(cat_pix["PIXFOVPROB"]))]
    first_values = sortcat[:maxRuns]

    if doPlot:
        path = dirName + "/GridPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)

        # mpl.rcParams.update({'font.size':14})
        hp.mollview(prob)
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
        plt.savefig("%s/Occ_Pointing.png" % (path))
        plt.close()

    return first_values


def PlotSpaceOcc(prob, dirName, reducedNside, Occultedpixels, first_values):
    path = dirName + "/Occ_Space_Obs"
    if not os.path.exists(path):
        os.mkdir(path, 493)

    # mpl.rcParams.update({'font.size':14})
    hp.mollview(prob)
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
        print("No pcculted pix")

    hp.visufunc.projplot(
        first_values["PIXRA"],
        first_values["PIXDEC"],
        "b.",
        lonlat=True,
        coord="C",
        linewidth=0.1,
    )
    plt.savefig("%s/Occ_Pointing.png" % (path))
    plt.close()


def SubstractPointings2D(tpointingFile, prob, obspar, pixlist, pixlistHR):
    nside = obspar.reducedNside
    radius = obspar.FOV

    print("Subtracting pointings from " + tpointingFile)
    raPointing, decPointing = np.genfromtxt(
        tpointingFile,
        usecols=(2, 3),
        dtype="str",
        skip_header=1,
        delimiter=" ",
        unpack=True,
    )  # ra, dec in degrees
    raPointing = np.atleast_1d(raPointing)
    decPointing = np.atleast_1d(decPointing)

    coordinates = TransformRADec(raPointing, decPointing)
    P_GW = []
    for i, valuei in enumerate(raPointing):
        t = 0.5 * np.pi - coordinates[i].dec.rad
        p = coordinates[i].ra.rad
        # Get the pixels for the ipix_disc (low res)
        xyz = hp.ang2vec(t, p)
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
        effectiveipix_disc = []
        for j, valuej in enumerate(ipix_disc):
            if valuej not in pixlist:
                effectiveipix_disc.append(valuej)
            pixlist.append(valuej)
        P_GW.append(prob[effectiveipix_disc].sum())

        print(
            "Coordinates ra:",
            raPointing[i],
            "dec:",
            decPointing[i],
            "Pgw:",
            P_GW[i],
            "vs",
            prob[ipix_disc].sum(),
        )
        # Save the ipixels in HR
        ipix_discHR = hp.query_disc(obspar.HRnside, xyz, np.deg2rad(radius))
        for k, valuek in enumerate(ipix_discHR):
            if valuek not in pixlistHR:
                pixlistHR.append(valuek)

    return pixlist, pixlistHR, np.sum(P_GW), len(raPointing)


def TransformRADec(vra, vdec):
    if "h" in vra[0]:
        ra = []
        dec = []
        for i in range(0, len(vra)):
            coord = SkyCoord(vra[i].split('"')[1], vdec[i].split('"')[0], frame="fk5")
            # print(coord)
            ra.append(coord.ra.deg)
            dec.append(coord.dec.deg)
    else:
        ra = vra.astype(float)
        dec = vdec.astype(float)

    coordinates = co.SkyCoord(ra, dec, frame="fk5", unit=(u.deg, u.deg))
    return coordinates


def LoadGalaxies(tgalFile):
    """
    Load galaxy catalog as an Astropy Table
    """

    print("Loading galaxy catalogue from " + tgalFile)

    # Load data
    h5file = tables.open_file(tgalFile, mode="r")
    tcat = Table(h5file.root.catalog.read()[["no_GLADE", "RA", "Dec", "d_L", "B_mag"]])
    h5file.close()

    # Rename column to match naming scheme
    tcat.rename_columns(
        ["RA", "Dec", "d_L", "B_mag"], ["RAJ2000", "DEJ2000", "Dist", "Bmag"]
    )

    return tcat


def LoadGalaxies_SteMgal(tgalFile):
    """
    Load galaxy catalog as an Astropy Table
    """

    print("Loading galaxy catalogue from " + tgalFile)

    # Load data
    h5file = tables.open_file(tgalFile, mode="r")
    tcat = Table(
        h5file.root.catalog.read()[["no_GLADE", "RA", "Dec", "d_L", "B_mag", "mass"]]
    )
    h5file.close()

    # Rename column to match naming scheme
    tcat.rename_columns(
        ["RA", "Dec", "d_L", "B_mag", "mass"],
        ["RAJ2000", "DEJ2000", "Dist", "Bmag", "SteMgal"],
    )

    return tcat


def FilterGalaxies(catalog, MinimumProbCutForCatalogue):
    """
    Filters galaxies to only keep the highest probabilities ones
    """
    min_prob_cut = catalog["dp_dV"] > MinimumProbCutForCatalogue * max(catalog["dp_dV"])
    Gals = catalog[min_prob_cut]
    # return array with list of Galaxies passing cuts, ordered by p-value

    tGals = Gals[np.flipud(np.argsort(Gals["dp_dV"]))]
    # ascii.write(tGals, '/Users/hashkar/Desktop/GWfollowup/GW-Followup/tGals_noM.txt', names = ['RAJ2000','DEJ2000','Dist','Bmag','dp_dV'],overwrite=True)
    return tGals


def MangroveGalaxiesProbabilities(catalog):
    """
    Computes new probabilities for each galaxy based on the Mangrove method
    """

    beta = 1
    alpha = 0
    Mgal1 = catalog["SteMgal"]
    Pgal_pos = catalog["dp_dV"]

    Mgal1 = np.nan_to_num(Mgal1)
    Mgal = 10 ** (Mgal1)
    Pgal_pos = np.nan_to_num(Pgal_pos)

    Gmass = Mgal / (np.sum(Mgal))
    alpha = (Pgal_pos).sum() / (Pgal_pos * Gmass).sum()
    catalog["dp_dV"] = (Pgal_pos) + (Pgal_pos * (alpha * beta * Gmass))

    return catalog


def VisibleAtTime(test_time, galaxies, maxz, observatory):
    """Determine if prompt or afterglow follow-up is possible by knowing if there are galaxies with non-negligible probability of hosting the NSM in the FoV
     1) check if any galaxy is visible, if not --> AFTERGLOW

    2) loop over zenith angle and select subsets of galaxies

    3) stop if maximum p-value of this subset is smaller than 75% of the previous subset

    4) else: stricter cut on zenith and repeat

    5) take galaxy with highest p-value fulfilling both criteria as target

    RETURNS:
    --------
    bool `is_vis` : is visible now?
    np.ndarray `alt_az` : alt_az location of galaxies
    """

    # print()
    # print("Check visibility at time {0}".format(test_time))

    # observatory time and location to look up visibility of objects

    # observatory = co.EarthLocation(lat=-23.271333 * u.deg,lon=16.5 * u.deg, height=1800 * u.m)

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


def FulfillsRequirementGreyObservations(
    Ktime, theseGals, observatory, minMoonSourceSeparation
):
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


def ObtainHighestProbabilityCoordinates(filename):
    hpx = UNIQSkymap_toNested(filename)
    hpx_prob = hpx["PROB"]
    print(sum(hpx_prob))
    ipix_max = np.argmax(hpx_prob)
    nside = hp.npix2nside(len(hpx_prob))
    theta, phi = hp.pix2ang(nside, ipix_max)
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    return ra, dec


def ComputeProbGalTargeted(
    prob,
    time,
    finalGals,
    visiGals,
    allGals,
    tsum_dP_dV,
    talreadysumipixarray,
    nside,
    thisminz,
    obspar,
    counter,
    dirName,
):
    """Computes probability Pgal and Pgw in FoV but it takes into account a list of pixels to avoid recounting already observed zones.
    Returns saved circle too (is it really needed? )
    bool doPlot when  = True is used to plot the maps

    RETURNS:

    --------

        P_Gal: Probability of galaxies within FoV in the LIGO signal region
        P_GW: Total probability within  FoV of the Ligo signal.
        noncircleGal: Table of galaxies that are outside the circle(s) and inside the LIGO signal region

    """
    observatory = obspar.location
    maxZenith = obspar.maxZenith
    FOV = obspar.FOV
    doPlot = obspar.doPlot

    targetCoord = co.SkyCoord(
        finalGals["RAJ2000"][:1],
        finalGals["DEJ2000"][:1],
        frame="fk5",
        unit=(u.deg, u.deg),
    )

    targetCoord2 = co.SkyCoord(
        visiGals["RAJ2000"], visiGals["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )

    targetCoord3 = co.SkyCoord(
        allGals["RAJ2000"], allGals["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )

    dp_dVfinal = visiGals["dp_dV"]

    # Array of indices of pixels inside circle of  FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoord[0].dec.rad
    p = targetCoord[0].ra.rad
    xyz = hp.ang2vec(t, p)

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    P_GW = prob[effectiveipix_disc].sum()

    P_Gal = (
        dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius].sum() / tsum_dP_dV
    )
    noncircleGal = allGals[targetCoord3.separation(targetCoord).deg > radius]

    if doPlot:

        path = dirName + "/EvolutionPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame="fk5", unit=(u.deg, u.deg))

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        hp.gnomview(
            prob,
            xsize=500,
            ysize=500,
            rot=[targetCoord.ra.deg, targetCoord.dec.deg],
            reso=5.0,
        )

        hp.graticule()
        # plt.savefig("%s/ExampleGW_%g.png" % (tname,j))

        # draw all galaxies within zenith-angle cut
        # hp.visufunc.projscatter(allGalsaftercuts['RAJ2000'], allGalsaftercuts['DEJ2000'], lonlat=True, marker='.',color='g', linewidth=0.1)
        # plt.savefig("%s/ExampleGW_Galaxies_%g.png" % (tname,j))

        # If I want to plot all gals, plot also the ones that are out of the circle
        # hp.visufunc.projscatter(noncircleGal['RAJ2000'], noncircleGal['DEJ2000'], lonlat=True, marker='*', color='g')

        # draw observation position, which is equivalent to galaxy with highest
        # probability
        # hp.visufunc.projscatter(finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], lonlat=True, marker='.', color='r',linewidth=0.1)

        # draw circle of FoV around best fit position

        hp.visufunc.projplot(
            skycoord[tempmask & tempmask2].ra,
            skycoord[tempmask & tempmask2].dec,
            "r.",
            lonlat=True,
            coord="C",
        )

        # Draw H.E.S.S. visibility

        altcoord = np.empty(4000)

        altcoord.fill(90 - maxZenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(
            azcoord,
            altcoord,
            frame="altaz",
            unit=(u.deg, u.deg),
            obstime=time,
            location=observatory,
        )

        RandomCoord_radec = RandomCoord.transform_to("fk5")

        hp.visufunc.projplot(
            RandomCoord_radec.ra, RandomCoord_radec.dec, "b.", lonlat=True, coord="C"
        )
        # MOON

        # Draw MinZ area

        altcoordmin = np.empty(4000)
        altcoordmin.fill(90 - thisminz)

        plt.savefig("%s/Zoom_Pointing_%g.png" % (path, counter))
        plt.close()

    return P_Gal, P_GW, noncircleGal, talreadysumipixarray


def SubstractPointings(
    tpointingFile, galaxies, talreadysumipixarray, tsum_dP_dV, prob, obspar, nside
):
    FOV = obspar.FOV

    # Read PointingsFile

    print("Subtracting pointings from " + tpointingFile)
    (
        rap,
        decP,
    ) = np.genfromtxt(
        tpointingFile,
        usecols=(2, 3),
        dtype="str",
        skip_header=1,
        delimiter=" ",
        unpack=True,
    )  # ra, dec in degrees

    coordinates = TransformRADec(rap, decP)
    ra = coordinates.ra.deg
    dec = coordinates.dec.deg

    PGW = []
    PGAL = []
    updatedGalaxies = galaxies
    if np.isscalar(ra):
        updatedGalaxies, pgwcircle, pgalcircle, talreadysumipixarray = (
            SubstractGalaxiesCircle(
                updatedGalaxies,
                ra,
                dec,
                talreadysumipixarray,
                tsum_dP_dV,
                FOV,
                prob,
                nside,
            )
        )
        PGW.append(pgwcircle)
        PGAL.append(pgalcircle)
        print(
            "Coordinates ra:", ra, "dec:", dec, "Pgw:", pgwcircle, "PGAL:", pgalcircle
        )
    else:
        for i, coord in enumerate(coordinates):
            ra = coord.ra.deg
            dec = coord.dec.deg
            updatedGalaxies, pgwcircle, pgalcircle, talreadysumipixarray = (
                SubstractGalaxiesCircle(
                    updatedGalaxies,
                    ra,
                    dec,
                    talreadysumipixarray,
                    tsum_dP_dV,
                    FOV,
                    prob,
                    nside,
                )
            )
            PGW.append(pgwcircle)
            PGAL.append(pgalcircle)
            print(
                "Coordinates ra:",
                ra,
                "dec:",
                dec,
                "Pgw:",
                pgwcircle,
                "PGAL:",
                pgalcircle,
            )
    return (
        ra,
        dec,
        updatedGalaxies,
        PGW,
        PGAL,
        talreadysumipixarray,
        len(np.atleast_1d(ra)),
    )


def SubstractGalaxiesCircle(
    galaux, ra, dec, talreadysumipixarray, tsum_dP_dV, FOV, prob, nside
):
    radius = FOV
    coordinates = co.SkyCoord(ra, dec, frame="fk5", unit=(u.deg, u.deg))

    targetCoord = co.SkyCoord(
        galaux["RAJ2000"], galaux["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )
    dp_dVfinal = galaux["dp_dV"]

    t = 0.5 * np.pi - coordinates.dec.rad
    p = coordinates.ra.rad
    xyz = hp.ang2vec(t, p)
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    P_GW = prob[effectiveipix_disc].sum()
    P_Gal = (
        dp_dVfinal[targetCoord.separation(coordinates).deg < radius].sum() / tsum_dP_dV
    )

    print("PGW", P_GW, "P_GAL", P_Gal)

    newgalaxies = galaux[targetCoord.separation(coordinates).deg > radius]

    return newgalaxies, P_GW, P_Gal, talreadysumipixarray


def ComputePGalinFOV(prob, cat, galpix, FOV, totaldPdV, nside, UsePix):
    """
    Computes probability Pgal in FoV
    """
    if UsePix:
        try:
            targetCoord = co.SkyCoord(
                galpix["PIXRA"], galpix["PIXDEC"], frame="fk5", unit=(u.deg, u.deg)
            )
        except Exception:
            targetCoord = galpix
    else:
        targetCoord = co.SkyCoord(
            galpix["RAJ2000"], galpix["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
        )

    targetCoord2 = co.SkyCoord(
        cat["RAJ2000"], cat["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )

    dp_dV = cat["dp_dV"]

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV

    # translate pixel indices to coordinates

    Pgal_inFoV = (
        dp_dV[targetCoord2.separation(targetCoord).deg <= radius].sum() / totaldPdV
    )

    return Pgal_inFoV


def ModifyCatalogue(prob, cat, FOV, totaldPdV, nside):
    """
    Computes the integrated Pgal in FoV for a list of calues using Pgal in FoV and sorts the catalog
    using that quantity as a criteria
    """
    lengthSG = 100
    SelectedGals = cat[:lengthSG]
    dp_dV_FOV = []
    for element in range(0, len(cat["dp_dV"])):
        if element < len(SelectedGals["dp_dV"]):
            dp_dV_FOV.append(
                ComputePGalinFOV(
                    prob,
                    cat,
                    SelectedGals[element],
                    FOV,
                    totaldPdV,
                    nside,
                    UsePix=False,
                )
            )
        else:
            dp_dV_FOV.append(0)

    cat["dp_dV_FOV"] = dp_dV_FOV

    tcat = cat[np.flipud(np.argsort(cat["dp_dV_FOV"]))]

    return tcat


def ComputeProbPGALIntegrateFoV(
    prob,
    time,
    observatory,
    centerPoint,
    UsePix,
    visiGals,
    allGalsaftercuts,
    tsum_dP_dV,
    talreadysumipixarray,
    nside,
    thisminz,
    obspar,
    counter,
    tname,
    dirName,
    doPlot,
):
    """
    Same as ComputeProbGalTargetted but it does not return circle coordinates.
    """
    maxZenith = obspar.maxZenith
    FOV = obspar.FOV
    if UsePix:
        try:
            targetCoord = co.SkyCoord(
                centerPoint["PIXRA"][:1],
                centerPoint["PIXDEC"][:1],
                frame="fk5",
                unit=(u.deg, u.deg),
            )
        except Exception:
            targetCoord = centerPoint

    else:
        targetCoord = co.SkyCoord(
            centerPoint["RAJ2000"][:1],
            centerPoint["DEJ2000"][:1],
            frame="fk5",
            unit=(u.deg, u.deg),
        )

    targetCoord2 = co.SkyCoord(
        visiGals["RAJ2000"], visiGals["DEJ2000"], frame="fk5", unit=(u.deg, u.deg)
    )

    targetCoord3 = co.SkyCoord(
        allGalsaftercuts["RAJ2000"],
        allGalsaftercuts["DEJ2000"],
        frame="fk5",
        unit=(u.deg, u.deg),
    )

    dp_dVfinal = visiGals["dp_dV"]

    # Array of indices of pixels inside circle of FoV

    radius = FOV
    t = 0.5 * np.pi - targetCoord[0].dec.rad

    p = targetCoord[0].ra.rad

    xyz = hp.ang2vec(t, p)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    P_GW = prob[effectiveipix_disc].sum()

    P_Gal = (
        dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius].sum() / tsum_dP_dV
    )

    # all galaxies inside the current observation circle

    # all galaxies outside the current observation circle, no visibility selection

    noncircleGal = allGalsaftercuts[targetCoord3.separation(targetCoord).deg > radius]

    if doPlot:

        path = dirName + "/EvolutionPlot"
        if not os.path.exists(path):
            os.mkdir(path, 493)

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame="fk5", unit=(u.deg, u.deg))

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        hp.gnomview(
            prob,
            xsize=500,
            ysize=500,
            rot=[targetCoord.ra.deg, targetCoord.dec.deg],
            reso=5.0,
        )

        hp.graticule()
        # draw circle of FoV around best fit position

        hp.visufunc.projplot(
            skycoord[tempmask & tempmask2].ra,
            skycoord[tempmask & tempmask2].dec,
            "r.",
            lonlat=True,
            coord="C",
        )

        # Draw H.E.S.S. visibility

        altcoord = np.empty(4000)

        altcoord.fill(90 - maxZenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(
            azcoord,
            altcoord,
            frame="altaz",
            unit=(u.deg, u.deg),
            obstime=time,
            location=observatory,
        )

        RandomCoord_radec = RandomCoord.transform_to("fk5")

        hp.visufunc.projplot(
            RandomCoord_radec.ra, RandomCoord_radec.dec, "b.", lonlat=True, coord="C"
        )
        # MOON

        altcoordmin = np.empty(4000)

        altcoordmin.fill(90 - thisminz)

        plt.savefig("%s/Zoom_Pointing_%g.png" % (path, counter))
        plt.close()

    return P_Gal, P_GW, noncircleGal, talreadysumipixarray


def Get90RegionPixReduced(hpxx, percentage, Nnside):
    nside = Nnside  # size of map used for contour determination
    hpx = hp.ud_grade(
        hpxx, nside_out=nside, power=-2, order_in="Nested", order_out="Nested"
    )

    sort = sorted(hpx, reverse=True)
    cumsum = np.cumsum(sort)
    index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

    # finding ipix indices confined in a given percentage
    index_hpx = range(0, len(hpx))
    hpx_index = np.c_[hpx, index_hpx]

    sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
    value_contour = sort_2array[0:index]

    j = 1
    table_ipix_contour = []

    for i in range(0, len(value_contour)):
        ipix_contour = int(value_contour[i][j])
        table_ipix_contour.append(ipix_contour)
    # from index to polar coordinates
    theta1, phi1 = hp.pix2ang(nside, table_ipix_contour)
    area = len(table_ipix_contour) * hp.nside2pixarea(nside, True)

    # reducing resolution to et a faser execution
    # list of pixel indices in the new map
    R_ipix = hp.ang2pix(Nnside, theta1, phi1)
    R_ipix = list(set(R_ipix))  # Removing/keeping 1 duplicate from list)

    # from index to polar coordinates
    theta, phi = hp.pix2ang(Nnside, R_ipix)

    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    return ra, dec, area


def Get90RegionPixGal(hpxx, percentage, Nside):
    hpx = hpxx
    sort = sorted(hpx, reverse=True)
    cumsum = np.cumsum(sort)
    index, _ = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

    # finding ipix indices confined in a given percentage
    index_hpx = range(0, len(hpx))
    hpx_index = np.c_[hpx, index_hpx]

    sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
    value_contour = sort_2array[0:index]

    j = 1
    table_ipix_contour = []

    for i in range(0, len(value_contour)):
        ipix_contour = int(value_contour[i][j])
        table_ipix_contour.append(ipix_contour)
    return table_ipix_contour


def IsSourceInside(Pointings, Sources, FOV, nside):
    tt = 0.5 * np.pi - Sources.dec.rad
    tp = Sources.ra.rad
    txyz = hp.ang2pix(nside, tt, tp)
    Npoiting = ""
    Found = False
    try:
        for i in range(0, len(Pointings)):
            t = 0.5 * np.pi - Pointings[i].dec.rad
            p = Pointings[i].ra.rad
            xyz = hp.ang2vec(t, p)
            try:
                ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))
            except Exception:
                ipix_disc = hp.query_disc(nside, xyz[0], np.deg2rad(FOV))
            if txyz in ipix_disc:
                print("Found in pointing number", i)
                # Npoiting.append(i)
                Npoiting = Npoiting + str(i) + ","
                Found = True
        if not Found:
            print("Source not covered!")
    except TypeError:
        t = 0.5 * np.pi - Pointings.dec.rad
        p = Pointings.ra.rad
        xyz = hp.ang2vec(t, p)
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))
        if txyz in ipix_disc:
            Npoiting = "0,"
            Found = True
            print("Found in pointing number 0")
        else:
            print("Source not covered!")
    # Reformat output
    if Found:
        Npoiting = Npoiting[:-1]
    return Found, Npoiting


def FillSummary(
    outfilename,
    ID,
    doneObservations,
    totalPoswindow,
    foundFirst,
    nP,
    totalProb,
    ObsInfo,
):
    f = open(outfilename, "w")
    f.write(
        "ID"
        + " "
        + "TotalObservations"
        + " "
        + "TotalPossible"
        + " "
        + "FirstCovered"
        + " "
        + "TimesFound"
        + " "
        + "TotalProb"
        + " "
        + "ObsInfo"
        + "\n"
    )
    f.write(
        str(ID)
        + " "
        + str(doneObservations)
        + " "
        + str(totalPoswindow)
        + " "
        + str(foundFirst)
        + " "
        + str(nP)
        + " "
        + str(totalProb)
        + " "
        + str(ObsInfo)
        + "\n"
    )


class NextWindowTools:

    @classmethod
    def CheckWindowCreateArray(cls, time, obsSite, WindowDurations):
        FullWindow = datetime.timedelta(seconds=np.float64(WindowDurations[-1]))
        if (Tools.IsDarkness(time, obsSite) is True) and (
            Tools.IsDarkness(time + FullWindow, obsSite) is True
        ):
            LastItem = len(WindowDurations)
        else:
            print("Window is smaller")
            for i in range(0, len(WindowDurations)):
                if Tools.IsDarkness(
                    time + datetime.timedelta(minutes=np.float64(WindowDurations[-i])),
                    obsSite,
                ):
                    LastItem = len(WindowDurations) - i
        cumsumWindow = np.cumsum(WindowDurations)
        arr = np.array(
            [
                time + datetime.timedelta(seconds=np.float64(cumsumWindow[j]))
                for j in range(LastItem)
            ]
        )
        return arr

    @classmethod
    def NextObservationWindow(cls, time, obspar):
        if (
            Tools.NextSunset(time, obspar).hour
            >= time.hour
            >= Tools.PreviousSunrise(time, obspar).hour
            and time.day == Tools.NextSunset(time, obspar).day
        ):
            time = Tools.NextSunset(time, obspar)
            time = Tools.TrustingDarknessSun(time, obspar)
        if Tools.IsDarkness(time, obspar) is True:
            return time
        elif (Tools.IsDarkness(time, obspar) is False) and (
            Tools.IsDarkness(Tools.NextMoonset(time, obspar), obspar) is True
        ):
            time = Tools.NextMoonset(time, obspar)
            return time
        else:
            print("No window is found")
            return False

    @classmethod
    def NextObservationWindowGrey(cls, time, obspar):
        if (
            Tools.NextSunset(time, obspar).hour
            >= time.hour
            >= Tools.PreviousSunrise(time, obspar).hour
            and time.day == Tools.NextSunset(time, obspar).day
        ):
            time = Tools.NextSunset(time, obspar)
            # print('Sunset', time)
            time = Tools.TrustingGreynessSun(time, obspar)
            # print('Trusted', time)
        if Tools.IsGreyness(time, obspar) is True:
            return time
        elif Tools.IsGreyness(time, obspar) is False:
            time = Tools.TrustingGreynessSun(time, obspar)
            # time=Tools.NextMoonset(time, obsSite)
            return time
        else:
            print("No window is found")
            return False

    @classmethod
    def EndObservationWindow(cls, time, obsSite):
        time = Tools.NextSunrise(time, obsSite)
        # Check if the night ends before due to the moon.
        if Tools.IsDarkness(time, obsSite) is False:
            time = Tools.PreviousMoonset(time, obsSite)
        return time

    @classmethod
    def AddRunDuration_to_StartingTime(cls, obspar):
        previousTime = np.genfromtxt(
            obspar.pointingsFile, usecols=(0), skip_header=1, unpack=True, dtype="str"
        )
        obspar.obsTime = datetime.datetime.strptime(
            str(previousTime), "%Y-%m-%dT%H:%M:%S"
        ) + datetime.timedelta(minutes=np.float64(obspar.duration))
