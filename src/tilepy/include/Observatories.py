from astropy import units as u
from astropy.coordinates import EarthLocation

__all__ = [
    "HESSObservatory",
    "LST",
    "CTASouthObservatory",
    "CTANorthObservatory",
]


class HESSObservatory:

    def __init__(self):
        self.Name = "HESS"
        self.Lat = -23.271778 * u.deg
        self.Lon = 16.50022 * u.deg
        self.Height = 1835 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class LST:

    def __init__(self):
        self.Name = "LST"
        self.Lat = 28.75 * u.deg
        self.Lon = -17.5 * u.deg
        self.Height = 2200 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class CTASouthObservatory:
    def __init__(self):
        self.Name = "South"
        self.Lat = -24.5 * u.deg
        self.Lon = -70.17 * u.deg
        self.Height = 2635 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class CTANorthObservatory:

    def __init__(self):
        self.Name = "North"
        self.Lat = 28.75 * u.deg
        self.Lon = -17.5 * u.deg
        self.Height = 2200 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)
