# Copyright (C) 2016-2025  tilepy developers
# (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

##################################################################################################
#                            Classes and tools for campaign definition                           #
##################################################################################################

#####################################################################
# Packages
import six
from astropy import units as u
from astropy.coordinates import EarthLocation
from six.moves import configparser

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


__all__ = ["ObservationParameters", "set_gaussian_source"]


def set_gaussian_source(obspar, ra, dec, sigma, name="gaussian_event"):
    """
    Configure an ObservationParameters instance to use a simulated Gaussian probability map.

    Parameters
    ----------
    obspar : ObservationParameters
        The instance to configure.
    ra : float
        Right Ascension of the source in degrees.
    dec : float
        Declination of the source in degrees.
    sigma : float
        Standard deviation (1-sigma) of the Gaussian in degrees.
    name : str, optional
        Event name to assign if not already set (default is "gaussian_event").

    Returns
    -------
    None
        This function modifies `obspar` in place and returns nothing.

    Example
    -------
    >>> obspar = ObservationParameters()
    >>> set_gaussian_source(obspar, ra=180.0, dec=30.0, sigma=2.5)
    >>> print(obspar.mode)
    gaussian

    """

    obspar.raSource = ra
    obspar.decSource = dec
    obspar.sigmaSource = sigma
    obspar.mode = "gaussian"
    if not hasattr(obspar, "event_name") or obspar.event_name is None:
        obspar.event_name = name


class ObservationParameters(object):
    """
    Stores all the configuration parameters from the .ini file

    This class collects observatory, scheduling, and event parameters used for observation planning.
    Attributes are typically loaded via the `from_configfile()` method.

    See source code for the complete list of available attributes.
    """

    # Observatory

    def __init__(
        self,
        obs_name=None,
        event_name=None,
        lat=0,
        lon=0,
        height=0,
        wobbleOffset=0,
        sunDown=None,
        moonDown=None,
        EarthDown=None,
        moonGrey=None,
        moonPhase=None,
        minMoonSourceSeparation=None,
        maxMoonSourceSeparation=None,
        SAAThreshold=0,
        maxZenith=None,
        FOV=None,
        maxRuns=None,
        maxNights=None,
        duration=None,
        minDuration=None,
        useGreytime=None,
        minSlewing=0,
        locCut=None,
        minimumProbCutForCatalogue=None,
        minProbcut=None,
        distCut=None,
        doPlot=False,
        secondRound=None,
        zenithWeighting=None,
        percentageMOC=None,
        reducedNside=None,
        HRnside=None,
        mangrove=None,
        skymap=None,
        mode=None,
        obsTime=None,
        datasetDir=None,
        galcatName=None,
        outDir=None,
        pointingsFile=None,
        countPrevious=False,
        MO=False,
        algorithm=None,
        strategy=None,
        doRank=False,
        downloadMaxRetry=0,
        downloadWaitPeriodRetry=20,
        shape=None,
        numberSides=None,
        igrfcoeffs=None,
        FoVRotation=None,
        alphaR=None,
        betaR=None,
    ):
        self.obs_name = obs_name
        self.event_name = event_name
        self.lat = lat
        self.lon = lon
        self.height = height
        self.location = EarthLocation(lat=self.lat, lon=self.lon, height=self.height)
        self.wobbleOffset = wobbleOffset * u.deg

        # Visibility
        self.sunDown = sunDown
        self.moonDown = moonDown
        self.EarthDown = EarthDown
        self.moonGrey = moonGrey
        self.moonPhase = moonPhase
        self.minMoonSourceSeparation = minMoonSourceSeparation
        self.maxMoonSourceSeparation = maxMoonSourceSeparation
        self.SAAThreshold = SAAThreshold
        self.igrfcoeffs = igrfcoeffs

        # Operations
        self.maxZenith = maxZenith
        self.FOV = FOV
        self.FoVRotation = FoVRotation
        self.maxRuns = maxRuns
        self.maxNights = maxNights
        self.duration = duration
        self.minDuration = minDuration
        self.useGreytime = useGreytime
        self.minSlewing = minSlewing
        self.shape = (shape,)
        self.numberSides = (numberSides,)

        # Tiling
        self.locCut = locCut
        self.minimumProbCutForCatalogue = minimumProbCutForCatalogue
        self.minProbcut = minProbcut
        self.distCut = distCut
        self.doPlot = doPlot
        self.secondRound = secondRound
        self.zenithWeighting = zenithWeighting
        self.percentageMOC = percentageMOC
        self.reducedNside = reducedNside
        self.HRnside = HRnside
        self.mangrove = mangrove
        self.algorithm = algorithm
        self.strategy = strategy
        self.doRank = doRank
        self.countPrevious = countPrevious
        self.alphaR = (alphaR,)
        self.betaR = (betaR,)

        # Parsed args
        self.skymap = skymap
        self.obsTime = obsTime
        self.datasetDir = datasetDir
        self.galcatName = galcatName
        self.outDir = outDir
        self.pointingsFile = pointingsFile

        # Download arguments
        self.downloadMaxRetry = downloadMaxRetry
        self.downloadWaitPeriodRetry = downloadWaitPeriodRetry

        # Characterstics of the event
        self.MO = MO

        self.mode = mode or None

        # Source localization parameters (used in "gaussian" mode)
        self.raSource = None
        self.decSource = None
        self.sigmaSource = None

    def __str__(self):
        return "\n".join(
            [
                "============== Observation Parameters ======================",
                f"Observatory Name: {self.obs_name}",
                f"Event Name: {self.event_name}",
                f"obsTime: {self.obsTime}",
                "---------------------- Strategy ----------------------",
                f"Algorithm = {self.algorithm}, Strategy = {self.strategy},  Mangrove = {self.mangrove}",
                f"Do Plot = {self.doPlot}, Do Rank = {self.doRank}, Count Previous= {self.countPrevious}, Second Round= {self.secondRound}, Use Grey Time= {self.useGreytime}",
                "--------------------- Observatory ---------------------",
                f"Observatory Location: {self.lat}, {self.lon}, {self.height}",
                f"Wobble Offset: {self.wobbleOffset}",
                f"FOV: {self.FOV}, Duration: {self.duration}, Min Duration: {self.minDuration}, Min Slewing: {self.minSlewing}",
                f"Max Runs: {self.maxRuns}, Max Nights: {self.maxNights}",
                f"Visibility: {self.sunDown}, {self.moonDown}, {self.moonGrey}, {self.moonPhase}, {self.EarthDown}",
                f"Min Moon Source Separation: {self.minMoonSourceSeparation}",
                f"Max Moon Source Separation: {self.maxMoonSourceSeparation}",
                f"Geomagnetic Threshold for SAA: {self.SAAThreshold}",
                f"Max Zenith: {self.maxZenith}, Zenith Weighting: {self.zenithWeighting}",
                f"FoV number of sides: {self.numberSides}, "
                f"FoV rotation: {self.FoVRotation},"
                f"Priority for FoV proximity and Probability: {self.alphaR}, Zenith Weighting: {self.betaR}",
                "--------------------- Skymap considerations ----------------",
                f"Skymap: {self.skymap}",
                f"Cuts: MinProbcut {self.minProbcut}, Dist Cut: {self.distCut}, Minimum Prob Cut for Catalogue: {self.minimumProbCutForCatalogue}",
                f"Percentage MOC: {self.percentageMOC}",
                f"NSIDE: HR = {self.HRnside}, reduced = {self.reducedNside}",
                "--------------------- Directories and files ----------------",
                f"DatasetDir: {self.datasetDir}",
                f"Galaxy Catalog Name: {self.galcatName}",
                f"Geomagnetic Coefficient Data Name: {self.igrfcoeffs}",
                f"Output Directory: {self.outDir}",
                f"Pointings File: {self.pointingsFile}",
                "============================================================",
            ]
        )

    def add_parsed_args(
        self,
        skymap,
        obsTime,
        datasetDir,
        galcatName,
        outDir,
        pointingsFile,
        igrfcoeffs=None,
        eventName=None,
        mode="healpix",
        ra=None,
        dec=None,
        sigma=None,
        nside=None,
    ):
        """Update instance attributes from parsed command-line arguments."""

        # Parsed args in command line
        self.skymap = skymap
        self.obsTime = obsTime
        self.datasetDir = datasetDir
        self.galcatName = galcatName
        self.igrfcoeffs = igrfcoeffs
        self.outDir = outDir
        self.pointingsFile = pointingsFile
        self.event_name = self.event_name if eventName is None else eventName
        self.mode = mode
        self.raSource = ra
        self.decSource = dec
        self.sigmaSource = sigma
        self.nside = nside

    def from_configfile(self, filepath):
        """Update instance attributes using parsed command-line arguments."""
        ##################
        cfg = filepath
        parser = ConfigParser()
        parser.read(cfg)
        parser.sections()
        section = "observatory"
        self.obs_name = str(parser.get(section, "name", fallback=None))
        self.lat = float(parser.get(section, "lat", fallback=0)) * u.deg
        self.lon = float(parser.get(section, "lon", fallback=0)) * u.deg
        self.height = float(parser.get(section, "height", fallback=0)) * u.m
        self.location = EarthLocation(lat=self.lat, lon=self.lon, height=self.height)
        self.base = str(parser.get(section, "base", fallback=None))
        self.stationsurl = str(parser.get(section, "stationsurl", fallback=None))
        self.wobbleOffset = (
            float(parser.get(section, "wobbleoffset", fallback=0)) * u.deg
        )

        section = "visibility"
        self.sunDown = int(parser.get(section, "sundown", fallback=0))
        self.moonDown = float(parser.get(section, "moondown", fallback=0))
        self.EarthDown = float(parser.get(section, "earthdown", fallback=0))
        # Altitude in degrees
        self.moonGrey = int(parser.get(section, "moongrey", fallback=0))
        self.moonPhase = int(
            parser.get(section, "gmoonphase", fallback=0)
        )  # Phase in %
        self.minMoonSourceSeparation = int(
            parser.get(section, "minmoonsourceseparation", fallback=0)
        )  # Separation in degrees
        self.maxMoonSourceSeparation = int(
            parser.get(section, "maxmoonsourceseparation", fallback=0)
        )  # Max separation in degrees
        self.SAAThreshold = int(parser.get(section, "SAAThreshold", fallback=0))

        section = "operations"
        self.maxZenith = int(parser.get(section, "maxzenith", fallback=0))
        self.FOV = float(parser.get(section, "fov", fallback=0))
        self.maxRuns = int(parser.get(section, "maxRuns", fallback=0))
        self.maxNights = int(parser.get(section, "maxNights", fallback=0))
        self.duration = float(parser.get(section, "duration", fallback=0))
        self.minDuration = float(parser.get(section, "minduration", fallback=0))
        self.useGreytime = parser.getboolean(section, "useGreytime", fallback=0)
        self.minSlewing = float(parser.get(section, "minSlewing", fallback=0))
        self.shape = str(parser.get(section, "shape", fallback=None))
        self.numberSides = int(parser.get(section, "numberSides", fallback=0))
        self.FoVRotation = int(parser.get(section, "FoVRotation", fallback=0))

        section = "tiling"
        self.locCut = float(parser.get(section, "locCut", fallback=99999))
        self.minimumProbCutForCatalogue = float(
            parser.get(section, "minimumprobcutforcatalogue", fallback=0)
        )
        self.minProbcut = float(parser.get(section, "minProbcut", fallback=0))
        self.distCut = float(parser.get(section, "distcut", fallback=0))
        self.doPlot = parser.getboolean(section, "doPlot", fallback=None)
        self.secondRound = parser.getboolean(section, "secondRound", fallback=None)
        self.zenithWeighting = float(parser.get(section, "zenithWeighting", fallback=0))
        self.percentageMOC = float(parser.get(section, "percentageMOC", fallback=0.90))
        try:
            self.reducedNside = int(parser.get(section, "reducedNside", fallback=0))
        except Exception:
            self.reducedNside = parser.getboolean(
                section, "reducedNside", fallback=None
            )

        try:
            self.HRnside = int(parser.get(section, "hrnside", fallback=0))
        except Exception:
            self.HRnside = parser.getboolean(section, "hrnside", fallback=None)
        self.mangrove = parser.getboolean(section, "mangrove", fallback=None)
        self.algorithm = str(parser.get(section, "algorithm", fallback=None))
        self.strategy = str(parser.get(section, "strategy", fallback=None))
        self.doRank = parser.getboolean(section, "doRank", fallback=None)
        self.countPrevious = parser.getboolean(section, "countPrevious", fallback=None)
        self.alphaR = float(parser.get(section, "alphaR", fallback=0))
        self.betaR = float(parser.get(section, "betaR", fallback=0))

        section = "general"
        self.downloadMaxRetry = int(parser.get(section, "downloadMaxRetry", fallback=0))
        self.downloadWaitPeriodRetry = float(
            parser.get(section, "downloadWaitPeriodRetry", fallback=0)
        )
        if parser.has_option(section, "eventName"):
            self.event_name = parser.get(section, "eventName")

    def from_args(
        self,
        obsName,
        eventName,
        lat,
        lon,
        height,
        sunDown,
        moonDown,
        EarthDown,
        moonGrey,
        moonPhase,
        minMoonSourceSeparation,
        maxMoonSourceSeparation,
        maxZenith,
        FOV,
        maxRuns,
        maxNights,
        duration,
        minDuration,
        useGreytime,
        minSlewing,
        minimumProbCutForCatalogue,
        minProbcut,
        distCut,
        doPlot,
        secondRound,
        zenithWeighting,
        percentageMOC,
        reducedNside,
        HRnside,
        mangrove,
    ):
        """Set instance attributes using direct function arguments."""

        self.obs_name = obsName
        self.event_name = eventName
        self.lat = lat * u.deg
        self.lon = lon * u.deg
        self.height = height * u.m
        self.location = EarthLocation(lat=self.lat, lon=self.lon, height=self.height)

        # Visibility
        self.sunDown = sunDown
        self.moonDown = moonDown
        self.EarthDown = EarthDown
        self.moonGrey = moonGrey
        self.moonPhase = moonPhase
        self.minMoonSourceSeparation = minMoonSourceSeparation
        self.maxMoonSourceSeparation = maxMoonSourceSeparation

        # Operations
        self.maxZenith = maxZenith
        self.FOV = FOV
        self.maxRuns = maxRuns
        self.maxNights = maxNights
        self.duration = duration
        self.minDuration = minDuration
        self.useGreytime = useGreytime
        self.minSlewing = minSlewing

        # Tiling
        self.minimumProbCutForCatalogue = minimumProbCutForCatalogue
        self.minProbcut = minProbcut
        self.distCut = distCut
        self.doPlot = doPlot
        self.secondRound = secondRound
        self.zenithWeighting = zenithWeighting
        self.percentageMOC = percentageMOC
        self.reducedNside = reducedNside
        self.HRnside = HRnside
        self.mangrove = mangrove
