import datetime
import numpy as np
from pytz import timezone
from skyfield import almanac
from skyfield.api import wgs84, N, E, load


class Observatory:
    """Class to store information and handle operation related to the observatory used for the observations."""

    def __init__(self, longitude, latitude, elevation, run_duration,
                 minimal_run_duration, max_sun_altitude, max_moon_altitude,
                 max_moon_phase):
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

        self.eph = load('de440s.bsp')
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
            zone = timezone('UTC')
            start_time = zone.localize(start_time)
        delta_time_obs = datetime.timedelta(days=nb_observation_night + 1)
        stop_time = start_time + delta_time_obs

        time_interval_sun = self.get_sun_constraint_time_interval(start_time, stop_time, nb_observation_night)
        time_interval_moon = self.get_moon_constraint_time_interval(start_time, stop_time)
        valid_time_interval = self.compute_interval_intersection(time_interval_sun, time_interval_moon)
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
            nb_observation_run = int(np.rint(observation_time_available // self.run_duration))
            remaining_observation_time = observation_time_available % self.run_duration
            if remaining_observation_time > self.minimal_run_duration:
                nb_observation_run += 1
            run_start_time += list(valid_time_range[i][0] + np.arange(nb_observation_run) * self.run_duration)
        return run_start_time

    def get_sun_constraint_time_interval(self, start_time, stop_time, nb_observation_night):
        """
        Get the time interval for the sun constraint.

        Args:
        start_time (datetime): The start time of observations.
        stop_time (datetime): The stop time of observations.
        nb_observation_night (int): The number of observation nights.

        Returns:
        list: The time interval for the sun constraint.
        """
        rise_time, set_time = self.get_risings_and_settings('sun', self.max_sun_altitude, start_time, stop_time)

        time_interval_sun = []
        for i in range(nb_observation_night):
            if set_time[i] < rise_time[i]:
                time_interval_sun.append([set_time[i], rise_time[i]])
            else:
                time_interval_sun.append([set_time[i], rise_time[i + 1]])

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
        rise_time, set_time = self.get_risings_and_settings('moon', self.max_moon_altitude, start_time, stop_time)

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
                time_middle_window = rise_time[i] + (set_time[j] - rise_time[i]) / 2.
                valid_interval = self.get_moon_phase(time_middle_window) < self.max_moon_phase
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
        sun, moon, earth = self.eph['sun'], self.eph['moon'], self.eph['earth']

        return earth.at(self.timescale_converter.from_datetime(observation_time)).observe(
            moon).apparent().fraction_illuminated(sun)

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
        f = almanac.risings_and_settings(self.eph, self.eph[celestial_body], self.observatory_location,
                                         horizon_degrees=horizon)
        time, rising_indicator = almanac.find_discrete(self.timescale_converter.from_datetime(start_time),
                                                       self.timescale_converter.from_datetime(stop_time),
                                                       f)
        rise_time = list(time[rising_indicator == 1].utc_datetime())
        set_time = list(time[rising_indicator == 0].utc_datetime())

        # Set the start time and end time as either rise of set time to fully cover the time range
        if rising_indicator[0] == 0:
            rise_time = [start_time, ] + rise_time
        else:
            set_time = [start_time, ] + set_time
        if rising_indicator[-1] == 1:
            set_time = set_time + [stop_time, ]
        else:
            rise_time = rise_time + [stop_time, ]

        return rise_time, set_time
