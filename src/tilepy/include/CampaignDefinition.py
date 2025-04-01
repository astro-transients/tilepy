#####################################################################
# Packages
from astropy import units as u
from astropy.coordinates import EarthLocation
from six.moves import configparser
import six

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


__all__ = ["ObservationParameters"]


class ObservationParameters(object):
    """Stores all the parameters in the .ini file"""

    # Observatory

    def __init__(
        self,
        obs_name=None,
        event_name=None,
        lat=0,
        lon=0,
        height=0,
        sunDown=None,
        moonDown=None,
        moonGrey=None,
        moonPhase=None,
        minMoonSourceSeparation=None,
        maxMoonSourceSeparation=None,
        maxZenith=None,
        FOV=None,
        maxRuns=None,
        maxNights=None,
        duration=None,
        minDuration=None,
        useGreytime=None,
        minSlewing=None,
        locCut90=None,
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
    ):
        self.obs_name = obs_name
        self.event_name = event_name
        self.lat = lat
        self.lon = lon
        self.height = height
        self.location = EarthLocation(lat=self.lat, lon=self.lon, height=self.height)

        # Visibility
        self.sunDown = sunDown
        self.moonDown = moonDown
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
        self.locCut90 = locCut90
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

    def __str__(self):
        txt = ""
        txt += "============== Main parsed observation parameters ==============  \n".format()
        txt += "Observatory Name: {}\n".format(self.obs_name)
        txt += "Event Name: {}\n".format(self.event_name)
        txt += "Observatory: {}\n".format(self.lat)
        txt += "Observatory: {}\n".format(self.lon)
        txt += "Observatory: {}\n".format(self.height)
        txt += "Max zenith: {}\n".format(self.maxZenith)
        txt += "Using Greytime is: {}\n".format(self.useGreytime)
        txt += "FOV: {}\n".format(self.FOV)
        txt += "Max runs: {}\n".format(self.maxRuns)
        txt += "Duration: {}\n".format(self.duration)
        txt += "High Resolution NSIDE: {}\n".format(self.HRnside)
        txt += "Low Resolution NSIDE: {}\n".format(self.reducedNside)
        txt += (
            "The strategy is ({algorithm}, {strategy}, mangrove={mangrove})\n".format(
                algorithm=self.algorithm, strategy=self.strategy, mangrove=self.mangrove
            )
        )
        txt += "The level of details is (doPlot={doPlot}, doRank = {doRank})\n".format(
            doPlot=self.doPlot, doRank=self.doRank
        )

        # txt += '----------------------------------------------------------------------\n'.format()
        return txt

    def add_parsed_args(
        self,
        skymap,
        obsTime,
        datasetDir,
        galcatName,
        outDir,
        pointingsFile,
        eventName=None,
    ):
        # Parsed args in command line
        self.skymap = skymap
        self.obsTime = obsTime
        self.datasetDir = datasetDir
        self.galcatName = galcatName
        self.outDir = outDir
        self.pointingsFile = pointingsFile
        self.event_name = self.event_name if eventName is None else eventName

    def from_configfile(self, filepath):
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

        section = "visibility"
        self.sunDown = int(parser.get(section, "sundown", fallback=0))
        self.moonDown = float(parser.get(section, "moondown", fallback=0))
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

        section = "operations"
        self.maxZenith = int(parser.get(section, "maxzenith", fallback=0))
        self.FOV = float(parser.get(section, "fov", fallback=0))
        self.maxRuns = int(parser.get(section, "maxRuns", fallback=0))
        self.maxNights = int(parser.get(section, "maxNights", fallback=0))
        self.duration = int(parser.get(section, "duration", fallback=0))
        self.minDuration = int(parser.get(section, "minduration", fallback=0))
        self.useGreytime = parser.getboolean(section, "useGreytime", fallback=0)
        self.minSlewing = float(parser.get(section, "minslewing", fallback=0))

        section = "tiling"
        self.locCut90 = float(parser.get(section, "locCut90", fallback=99999))
        self.minimumProbCutForCatalogue = float(
            parser.get(section, "minimumprobcutforcatalogue", fallback=0)
        )
        self.minProbcut = float(parser.get(section, "minProbcut", fallback=0))
        self.distCut = float(parser.get(section, "distcut", fallback=0))
        self.doPlot = parser.getboolean(section, "doPlot", fallback=None)
        self.secondRound = parser.getboolean(section, "secondRound", fallback=None)
        self.zenithWeighting = float(parser.get(section, "zenithWeighting", fallback=0))
        self.percentageMOC = float(parser.get(section, "percentageMOC", fallback=0))
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

        section = "general"
        self.downloadMaxRetry = int(
            parser.getboolean(section, "downloadMaxRetry", fallback=0)
        )
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
        self.obs_name = obsName
        self.event_name = eventName
        self.lat = lat * u.deg
        self.lon = lon * u.deg
        self.height = height * u.m
        self.location = EarthLocation(lat=self.lat, lon=self.lon, height=self.height)

        # Visibility
        self.sunDown = sunDown
        self.moonDown = moonDown
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
