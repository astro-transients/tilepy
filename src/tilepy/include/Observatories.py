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
    Coordinates and basic information for the H.E.S.S. Observatory.

    The **High Energy Stereoscopic System (H.E.S.S.)** is an array of five Imaging Atmospheric Cherenkov Telescopes (IACTs)
    located in Namibia, dedicated to very high-energy gamma-ray astronomy. The array consists of four 12 m telescopes
    arranged in a 120 m square for stereoscopic imaging, and a central 28 m telescope (operational since 2012).
    H.E.S.S. has provided full service since 2022, enabling major discoveries of galactic and extragalactic sources
    of gamma rays :footcite:`2018A&A...612A...1H`.

    - **Latitude**: -23.271778°
    - **Longitude**: 16.50022°
    - **Altitude**: 1835 m

    More information: `H.E.S.S. website <https://www.mpi-hd.mpg.de/HESS/>`_

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
    r"""
    Coordinates and description for the Large-Sized Telescope (LST), CTAO.

    The **Large-Sized Telescope (LST)** is the flagship instrument of the Cherenkov Telescope Array Observatory (CTAO),
    optimized for the detection of gamma rays at the lowest energies accessible from the ground (20-150 GeV).
    Each LST features a 23-meter segmented mirror (collecting area : :math:`400 m^2`) and a field of view of 4.3°, allowing
    high sensitivity to faint Cherenkov flashes :footcite:`2016APh....72...76A`.

    Four LSTs are being installed at the center of the CTAO northern site, enabling rapid and sensitive follow-up
    of transient events. Thanks to their lightweight structure, LSTs can repoint in less than 20 seconds, which is
    crucial for observing phenomena such as gamma-ray bursts.

    Each LST is equipped with a high-speed camera containing 1855 photomultiplier tubes (PMTs), able to record fast
    Cherenkov signals produced when gamma rays interact with the atmosphere. These data provide crucial insights
    into extreme cosmic sources.

    The LST, together with the Medium-Sized and Small-Sized Telescopes, forms the CTAO, offering broad energy
    coverage from 20 GeV to 300 TeV :footcite:`2023arXiv230512888H`.

    More information:
    `LST official page <https://www.ctao.org/emission-to-discovery/telescopes/lst/>`_

    .. figure:: ../_static/LST.jpg
        :alt: The Large-Sized Telescope (LST) at the CTAO site (Credit: Otger Ballester, IFAE).
        :width: 80%
        :align: center

        The Large-Sized Telescope (LST) at the CTAO site (Credit: Otger Ballester, IFAE).

    References
    ----------
    .. footbibliography::

    """

    def __init__(self):
        self.Name = "LST"
        self.Lat = 28.75 * u.deg
        self.Lon = -17.5 * u.deg
        self.Height = 2200 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class CTASouthObservatory:
    r"""
    Coordinates and site information for the CTAO-South Observatory.

    The **CTAO-South Observatory** is located in the Atacama Desert in Chile, about 10 km southeast of the ESO Paranal Observatory.
    This remote desert site is one of the driest places on Earth, providing excellent observing conditions for very high-energy gamma-ray astronomy.

    - **Latitude**: -24.5°
    - **Longitude**: -70.3°
    - **Altitude**: 2653 m

    The southern array is designed to detect gamma rays from 150 GeV up to 300 TeV, covering the highest energies accessible
    from the ground and focusing mainly on Galactic sources .

    The initial array—called the “Alpha Configuration”—includes:
    - 14 Medium-Sized Telescopes (MSTs), covering energies from 150 GeV to 5 TeV,
    - 37 Small-Sized Telescopes (SSTs), which extend sensitivity above 5 TeV,
    and the full array covers about :math:`3 km^{2}`.

    No Large-Sized Telescopes (LSTs) are installed at CTAO-South for now, but the site is prepared for future LST additions.

    `The SSTs <https://www.ctao.org/emission-to-discovery/telescopes/sst/>`_ use a dual-mirror design with silicon photomultiplier (SiPM) cameras
    and a wide field of view (about 8.8 degrees). `The MSTs <https://www.ctao.org/emission-to-discovery/telescopes/mst/>`_ have 12 m mirrors and fast PMT-based cameras (FlashCam), with a field of view of about 8 degrees.

    `The CTAO-South site <https://www.ctao.org/emission-to-discovery/array-sites/ctao-south/>`_  is built to explore some of the most powerful particle accelerators in our Galaxy,
    such as supernova remnants, pulsar wind nebulae, and massive star clusters :footcite:`The_CTA_Consortium_2019`.

    .. figure:: ../_static/CTAO_south.jpg
        :alt: An artistic illustration of the proposed CTA. Image credit: Gabriel Pérez Diaz, IAC / Marc-André Besel, CTAO.
        :width: 80%
        :align: center

        An artistic illustration of the proposed CTA (credit: Gabriel Pérez Diaz, IAC / Marc-André Besel.

    References
    ----------
    .. footbibliography::
    """

    def __init__(self):
        self.Name = "South"
        self.Lat = -24.5 * u.deg
        self.Lon = -70.17 * u.deg
        self.Height = 2635 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)


class CTANorthObservatory:
    r"""
    Coordinates and location for the CTA North Observatory.

    The **CTAO-North Observatory** is located at the Roque de los Muchachos Observatory on La Palma (Canary Islands, Spain), at :math:`\approx 2200` m altitude.
    This site is optimized for low- and mid-energy gamma rays (20 GeV to 5 TeV), with a layout focused on extragalactic sources.

    - **Latitude**: 28.762° N
    - **Longitude**: 17.891° W

    The “Alpha Configuration” includes:
      - 4 Large-Sized Telescopes (LSTs)
      - 9 Medium-Sized Telescopes (MSTs)
    covering about :math:`0.5 km^2`.

    As of mid‑2025:
      - **LST‑1**: installed in October 2018 and under commissioning
      - **LST‑2 / LST‑3 / LST‑4**: under construction, expected to be in place before end 2025
      - **9 MSTs**: being deployed on-site

    CTAO–North is ideal for studying low-energy cosmic transients, distant active galactic nuclei, and extragalactic gamma-ray sources.

    More information: `CTAO North site <https://www.ctao.org/emission-to-discovery/array-sites/ctao-north/>`_

    References
    ----------
    .. footbibliography::
    """

    def __init__(self):
        self.Name = "North"
        self.Lat = 28.75 * u.deg
        self.Lon = -17.5 * u.deg
        self.Height = 2200 * u.m
        self.location = EarthLocation(lat=self.Lat, lon=self.Lon, height=self.Height)
