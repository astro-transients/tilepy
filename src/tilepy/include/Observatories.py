from astropy import units as u
from astropy.coordinates import EarthLocation

__all__ = [
    "HESSObservatory",
    "LST",
    "CTASouthObservatory",
    "CTANorthObservatory",
]


class HESSObservatory:
    r"""
    Coordinates and location for the HESS Observatory.

    The **High Energy Stereoscopic System (H.E.S.S.)** is an array of five Imaging Atmospheric Cherenkov Telescopes (IACTs)
    located in Namibia, dedicated to very high-energy gamma-ray astronomy. The array consists of four 12 m telescopes
    arranged in a 120 m square for stereoscopic imaging, and a central 28 m telescope (operational since 2012).
    H.E.S.S. has provided full service since 2022, enabling major discoveries of galactic and extragalactic sources
    of gamma rays :footcite:`2018A&A...612A...1H,2022A&A...666A.124A`.

    For more information, see the `H.E.S.S. website <https://www.mpi-hd.mpg.de/HESS/>`_.

    References
    ----------
    .. footbibliography::
    """

    def __init__(self):
        self.Name = "HESS"
        self.Lat = -23.271778 * u.deg
        self.Lon = 16.50022 * u.deg
        self.Height = 1835 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class LST:
    """Coordinates and location for the LST Observatory."""

    def __init__(self):
        self.Name = "LST"
        self.Lat = 28.75 * u.deg
        self.Lon = -17.5 * u.deg
        self.Height = 2200 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class CTASouthObservatory:
    """Coordinates and location for the CTA South Observatory."""

    def __init__(self):
        self.Name = "South"
        self.Lat = -24.5 * u.deg
        self.Lon = -70.17 * u.deg
        self.Height = 2635 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class CTANorthObservatory:
    """Coordinates and location for the CTA North Observatory."""

    def __init__(self):
        self.Name = "North"
        self.Lat = 28.75 * u.deg
        self.Lon = -17.5 * u.deg
        self.Height = 2200 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)
