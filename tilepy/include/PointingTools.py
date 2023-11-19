#####################################################################
# Author: Monica Seglar-Arroyo
# Contributors: Halim Ashkar,  Fabian Schussler, Mathieu de Bony
# All the tools that are needed to follow-up a GW with an IACT (HESS)
# are described and implemented in the following.
#####################################################################
# Packages
import datetime
import os
import time

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
from astropy.coordinates import EarthLocation, get_sun
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import get_moon
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers
from astropy.coordinates import Angle
from gdpyc import DustMap
from mocpy import MOC
from scipy.stats import norm
from six.moves import configparser
import ligo.skymap.io.fits as lf
from .gwobserve import Sensitivity, GRB
from .observatory import Observatory
import pandas as pd
from astropy.table import QTable
import astropy_healpix as ah

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser


iers_file = os.path.join(os.path.abspath(
    os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)
# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))


##                      Classes                     ##

######################################################


class Tools:
    '''
        class with different visibility check functions and other setting and rising of the sun and the moon functions
        '''

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
        sun = get_sun(Time(obsTime, scale='utc')).transform_to(AltAz(obstime=Time(obsTime, scale='utc'),
                                                                     location=obspar.location))
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
        return moon.alt * 180. / np.pi

    @classmethod
    def NextSunrise(cls, obsTime, obspar):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        #obs.horizon = obspar.horizonSun
        obs.horizon = Angle(obspar.sunDown, u.deg).to_string(unit=u.degree, sep=':')
        sun.compute(obs)
        nextSunrise = obs.next_rising(
            sun, use_center=True).datetime().replace(tzinfo=pytz.utc)
        return nextSunrise

    @classmethod
    def PreviousSunrise(cls, obsTime, obspar):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        #obs.horizon = obspar.horizonSun
        obs.horizon = Angle(obspar.sunDown, u.deg).to_string(unit=u.degree, sep=':')
        sun.compute(obs)
        previousSunrise = obs.previous_rising(
            sun, use_center=True).datetime().replace(tzinfo=pytz.utc)
        return previousSunrise

    @classmethod
    def NextSunset(cls, obsTime, obspar):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        #obs.horizon = obspar.horizonSun
        obs.horizon = Angle(obspar.sunDown, u.deg).to_string(unit=u.degree, sep=':')
        sun.compute(obs)
        nextSunset = obs.next_setting(
            sun, use_center=True).datetime().replace(tzinfo=pytz.utc)
        return nextSunset

    @classmethod
    def PreviousMoonset(cls, obsTime, obspar):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        #obs.horizon = obspar.HorizonMoon
        obs.horizon = Angle(obspar.moonDown, u.deg).to_string(unit=u.degree, sep=':')
        moon.compute()
        previousMoonset = obs.previous_setting(
            moon, use_center=True).datetime().replace(tzinfo=pytz.utc)
        return previousMoonset

    @classmethod
    def NextMoonset(cls, obsTime, obspar):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obspar.lon / u.deg)
        obs.lat = str(obspar.lat / u.deg)
        obs.elev = obspar.height / u.m
        obs.date = obsTime  # Requires time in UTC
        #obs.horizon = obspar.HorizonMoon
        obs.horizon = Angle(obspar.moonDown, u.deg).to_string(unit=u.degree, sep=':')
        moon.compute()
        nextMoonset = obs.next_setting(
            moon, use_center=True).datetime().replace(tzinfo=pytz.utc)
        # print('NextMoonset',nextMoonset)
        previousMoonset = obs.previous_setting(
            moon, use_center=True).datetime().replace(tzinfo=pytz.utc)
        # print(previousMoonset)
        return nextMoonset

    @classmethod
    def TrustingDarknessSun(cls, obsTime, obspar):
        DarkObsTime = obsTime
        referencetime = obsTime
        while (Tools.IsDarkness(DarkObsTime, obspar) == False and ((DarkObsTime.hour >= referencetime.hour and DarkObsTime.day == referencetime.day) or (
                DarkObsTime.hour <= Tools.NextSunrise(referencetime, obspar).hour and DarkObsTime.day == Tools.NextSunrise(
                referencetime, obspar).day))):
            DarkObsTime = DarkObsTime + datetime.timedelta(minutes=1)
        return DarkObsTime

    @classmethod
    def TrustingGreynessSun(cls, obsTime, obspar):
        GreyObsTime = obsTime
        referencetime = obsTime
        while (Tools.IsGreyness(GreyObsTime, obspar) == False and ((GreyObsTime.hour >= referencetime.hour and GreyObsTime.day == referencetime.day) or (
                GreyObsTime.hour <= Tools.NextSunrise(referencetime, obspar).hour and GreyObsTime.day == Tools.NextSunrise(
                referencetime, obspar).day))):
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
        MinimalWindowDuration = datetime.timedelta(minutes=10)
        if (Tools.IsDarkness(time, obspar) is True) and (Tools.IsDarkness(time + MinimalWindowDuration, obspar) is True):
            Observe = True
        else:
            print('No window found')
            Observe = False

        return Observe

    @classmethod
    def CheckWindowGrey(cls, time, obspar):
        MinimalWindowDuration = datetime.timedelta(minutes=10)
        if (Tools.IsGreyness(time, obspar) is True) and (Tools.IsGreyness(time + MinimalWindowDuration, obspar) is True):
            Observe = True
        else:
            print('No window found')
            Observe = False

        return Observe

    @classmethod
    def GalacticPlaneBorder(cls, coords):
        lon = coords.galactic.l.value  # x-coordinate
        lat = coords.galactic.b.value  # y-coordinate
        #print(lon)
        #print(lat)
        YouAreInside = False
        n = 20
        if (lat <= 10 and lat >= 0 and lon <= 130):
            n = lat-(1.0/13)*lon
        elif (lat <= 10 and lat >= 0 and lon >= 240):
            n = lat - (1.0/12) * lon + 20
        elif (lat >= -10 and lat <= 0 and lon <= 130):
            n = lat + (1.0/13) * lon
        elif (lat >= -10 and lat <= 0 and lon >= 240):
            n = lat + (1.0/12) * lon - 20
        if np.absolute(n) <= 10:
            YouAreInside = True
        return YouAreInside

    @classmethod
    def GetGalacticExtinction(cls, coords, dustmap='SFD', filters='SDSS_r'):
        # Extinction = DustMap.ebv(coords)
        extinction = DustMap.extinction(
            coords, dustmap='SFD', filters='SDSS_r')
        # GasMap.plot_map('HI4PI')
        return extinction

#####################################################

# Parse Observation Parameters

#####################################################


class ObservationParameters(object):
    """Stores all the parameters in the .ini file"""
    # Observatory

    def __init__(self, name=None, lat=0, lon=0, height=0, sunDown=None, moonDown=None,
                 moonGrey=None, moonPhase=None, minMoonSourceSeparation=None,
                 maxMoonSourceSeparation=None, maxZenith=None, FOV=None, maxRuns=None, maxNights=None,
                 duration=None, minDuration=None, useGreytime=None, minSlewing=None, online=False,
                 minimumProbCutForCatalogue=None, minProbcut=None, distCut=None, doPlot=False, secondRound=None,
                 zenithWeighting=None, percentageMOC=None, reducedNside=None, HRnside=None,
                 mangrove=None, url=None,obsTime=None,datasetDir=None,galcatName=None,outDir=None,pointingsFile=None,alertType=None, locCut=None, MO=False, algorithm=None, strategy=None, doRank=False):

        self.name = name
        self.lat = lat
        self.lon = lon
        self.height = height
        self.location = EarthLocation(lat=self.lat, lon=self.lon,
                                      height=self.height)

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
        self.online = online
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

        # Parsed args
        self.url = url 
        self.obsTime = obsTime
        self.datasetDir = datasetDir
        self.galcatName = galcatName
        self.outDir = outDir
        self.pointingsFile = pointingsFile
        self.alertType = alertType
        self.locCut = locCut

        #Characterstics of the event
        self.MO = MO

    def __str__(self):
        txt = ''
        txt += '----------------- Main parsed observation parameters ----------------- \n'.format()
        txt += 'Observatory: {}\n'.format(self.lat)
        txt += 'Observatory: {}\n'.format(self.lon)
        txt += 'Observatory: {}\n'.format(self.height)
        txt += 'Name: {}\n'.format(self.name)
        txt += 'Max zenith: {}\n'.format(self.maxZenith)
        txt += 'FOV: {}\n'.format(self.FOV)
        txt += 'Max runs: {}\n'.format(self.maxRuns)
        txt += 'Duration: {}\n'.format(self.duration)
        txt += 'High Resolution NSIDE: {}\n'.format(self.HRnside)
        txt += 'Low Resolution NSIDE: {}\n'.format(self.reducedNside)
        txt += 'The strategy is ({algorithm},{strategy}, mangrove={mangrove})\n'.format(algorithm = self.algorithm, strategy = self.strategy,mangrove = self.mangrove)
        txt += 'The level of details is (doPlot={doPlot}, doRank = {doRank})\n'.format(doPlot = self.doPlot, doRank = self.doRank)

        # txt += '----------------------------------------------------------------------\n'.format()
        return txt

    def add_parsed_args(self, url,obsTime,datasetDir,galcatName,outDir,pointingsFile,alertType,locCut):
        # Parsed args in command line
        self.url = url 
        self.obsTime = obsTime 
        self.datasetDir = datasetDir
        self.galcatName = galcatName
        self.outDir = outDir
        self.pointingsFile = pointingsFile
        self.alertType = alertType
        self.locCut = locCut

    def from_configfile(self, filepath):

        ##################
        cfg = filepath
        parser = ConfigParser()
        parser.read(cfg)
        parser.sections()
        section = 'observatory'
        self.name = str(parser.get(section, 'name', fallback=None))
        self.lat = float(parser.get(section, 'lat', fallback=0))*u.deg
        self.lon = float(parser.get(section, 'lon', fallback=0))*u.deg
        self.height = float(parser.get(section, 'height', fallback=0))*u.m
        self.location = EarthLocation(lat=self.lat, lon=self.lon,
                                      height=self.height)

        section = 'visibility'
        self.sunDown = int(parser.get(section, 'sundown', fallback=0))
        self.moonDown = float(parser.get(section, 'moondown', fallback=0))
        # Altitude in degrees
        self.moonGrey = int(parser.get(section, 'moongrey', fallback=0))
        self.moonPhase = int(parser.get(
            section, 'gmoonphase', fallback=0))  # Phase in %
        self.minMoonSourceSeparation = int(parser.get(
            section, 'minmoonsourceseparation', fallback=0))  # Separation in degrees
        self.maxMoonSourceSeparation = int(parser.get(
            section, 'maxmoonsourceseparation', fallback=0))  # Max separation in degrees

        section = 'operations'
        self.maxZenith = int(parser.get(section, 'maxzenith', fallback=0))
        self.FOV = float(parser.get(section, 'fov', fallback=0))
        self.maxRuns = int(parser.get(section, 'maxRuns', fallback=0))
        self.maxNights = int(parser.get(section, 'maxNights', fallback=0))
        self.duration = int(parser.get(section, 'duration', fallback=0))
        self.minDuration = int(parser.get(section, 'minduration', fallback=0))
        self.useGreytime = (parser.getboolean(
            section, 'useGreytime', fallback=0))
        self.minSlewing = float(parser.get(section, 'minslewing', fallback=0))

        section = 'tiling'
        self.online = (parser.getboolean(section, 'online', fallback=None))
        self.minimumProbCutForCatalogue = float(parser.get(
            section, 'minimumprobcutforcatalogue', fallback=0))
        self.minProbcut = float(parser.get(section, 'minProbcut', fallback=0))
        self.distCut = float(parser.get(section, 'distcut', fallback=0))
        self.doPlot = (parser.getboolean(section, 'doPlot', fallback=None))
        self.secondRound = (parser.getboolean(
            section, 'secondRound', fallback=None))
        self.zenithWeighting = float(parser.get(
            section, 'zenithWeighting', fallback=0))
        self.percentageMOC = float(parser.get(
            section, 'percentageMOC', fallback=0))
        self.reducedNside = int(parser.get(
            section, 'reducedNside', fallback=0))
        self.HRnside = int(parser.get(section, 'hrnside', fallback=0))
        self.mangrove = (parser.getboolean(section, 'mangrove', fallback=None))
        self.algorithm = str(parser.get(section, 'algorithm', fallback=None))
        self.strategy = str(parser.get(section, 'strategy', fallback=None))
        self.doRank = (parser.getboolean(section, 'doRank', fallback=None))

    def from_args(self, name, lat, lon, height, sunDown, moonDown,
                  moonGrey, moonPhase, minMoonSourceSeparation,
                  maxMoonSourceSeparation, maxZenith, FOV, maxRuns, maxNights,
                  duration, minDuration, useGreytime, minSlewing, online,
                  minimumProbCutForCatalogue, minProbcut,distCut, doPlot, secondRound,
                  zenithWeighting, percentageMOC, reducedNside, HRnside,
                  mangrove):

        self.name = name
        self.lat = lat * u.deg
        self.lon = lon * u.deg
        self.height = height * u.m
        self.location = EarthLocation(lat=self.lat,
                                      lon=self.lon,
                                      height=self.height)

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
        self.online = online
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
        return datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    else:
        print("ERROR: something is wrong with the format of the date: ", x)
        return None


def GetGBMMap(URL):
    """
    Bottom-level function that takes a url searches for the localisatios maps from the Fermi-GBM database, or waits until it is uplaoded. 
    
    :param URL: the URL of the map
    :type URL: str

    :return: fitsfile, filename
    rtype: fits, str
    """
    
    filename = URL.split("/")[-1]
    filename = filename.split(".")[0]
    filename = "./maps/" + filename + ".fit"
    # filename = "glg_healpix_all_bn211130636_2.fits"
    print("The GBM filename is ", filename)
    try:
        fits_map_url_intial = URL
        fits_map_url1 = fits_map_url_intial.split("/")
        fits_map_url2 = fits_map_url_intial.split("_")[-1]
        # fits_map_url1[-1] = ""
        i = 0
        fits_map_url = ""
        for i in range(len(fits_map_url1) - 1):
            fits_map_url += fits_map_url1[i] + "/"
        fits_map_url += "glg_healpix_all" + "_" + \
            fits_map_url2.split(".")[0] + ".fit"

        # command = 'curl %s -o %s' % (fits_map_url, filename)
        # print(command)
        # os.system(command)
        # should fix this issue to save the fits file ... then comment the following line.
        filename = fits_map_url
    except:
        warn = "Caught exception: "
        print(warn)
        pass

    max_delay = 20
    delay = 0
    d = 0
    while delay == 0:
        delay = 1
        d = d + 1
        try:
            fitsfile = fits.open(filename)
        except:
            print('map is not uploaded yet... Waiting for minute:', d)
            time.sleep(60)
            delay = 0
            if d > max_delay:
                print(
                    f"Waited for {max_delay} minutes... can't wait anymore... I'm leaving")
                fitsfile = None
                filename = None
                break

    return fitsfile, filename


def GetGWMap_Flat(URL):
    """
    Bottom-level function that takes a url searches for the localisatios maps from the GW database, or waits until it is uplaoded. 
    
    :param URL: the URL of the map
    :type URL: str

    :return: fitsfile, filename
    rtype: fits, str
    """

    filename = URL.split("/")[-1]
    print("The filename is ", filename)
    fits_map_url = URL
    # fits_map_url = self.What["GW_SKYMAP"]['skymap_fits']['value']
    if 'multiorder.' in filename:
        fits_map_url = str(fits_map_url).replace('multiorder.', '')
        fits_map_url = str(fits_map_url).replace('fits', 'fits.gz')
        print('The GW map is in the right multiorder format')
    else:
        print('The GW map is not in multiorder format, we will try the .fits.gz format, you are welcome')

    newFilename = filename + "_"+str(int(time.time() * 1e6))
    print("internal filename: ", newFilename)

    try:
        command = f'curl {fits_map_url} -o {newFilename}'
        print(command)
        os.system(command)
    except x:
        print('Problem with downloading map from url, it was not multiorder or fits.gz')
        warn = "Caught exception: %s" % x
        print(warn)
        pass

    fitsfile = fits.open(newFilename)

    return fitsfile, newFilename


def GetGWMap(URL):
    """
    Bottom-level function that takes a url searches for the localisation maps in multi-order format from the GW database, or waits until it is uplaoded. 
    
    :param URL: the URL of the map
    :type URL: str

    :return: fitsfile, filename
    rtype: fits, str
    """
    filename = URL.split("/")[-1]
    print("The filename is ", filename)
    fits_map_url = URL
    try:
        command = f'curl {fits_map_url} -o {filename}'
        print(command)
        os.system(command)

    except x:
        print('Problem with downloading map from url, it was not multiorder or fits.gz')
        warn = "Caught exception: %s" % x
        print(warn)
        pass

    fitsfile = fits.open(filename)
    
    return fitsfile, filename

def UNIQSkymap_toNested(skymap_fname):
    """
    Bottom-level function that takes a skymap and computes from it the uniq map 
    
    :param skymap_fname: Healpix skymap
    :type skymap_fname: Table

    :return: healpix_skymaps_dict
    rtype: dict
    """

    sky_tab = Table.read(skymap_fname)
    healpix_skymaps_dict = get_lvk_uniq_maps(sky_tab, 'max')
    # prob = healpix_skymaps_dict['PROBDENSITY']
    return healpix_skymaps_dict


def get_lvk_uniq_maps(sky_map, Order, map_names='all'):

    un_inds = sky_map['UNIQ']

    order = (np.log2(un_inds / 4).astype(int) /
             2).astype(int)
    inds = (un_inds - 4 * (np.power(4, order))).astype(int)

    if Order == 'max':
        Order = np.max(order)

    Nside = int(2 ** Order)
    Npix = hp.nside2npix(Nside)

    if map_names == 'all':
        keys = ['PROB', 'DISTMU', 'DISTSIGMA', 'DISTNORM']
    else:
        keys = map_names
    maps = {}

    for k in keys:
        maps[k] = np.zeros(Npix)

    # print np.min(order), np.max(order)

    for ii in range(np.max(order), np.min(order) - 1, -1):

        nside = 2 ** ii
        npix = hp.nside2npix(nside)
        bl = (order == ii)

        for k in maps.keys():
            a = hp.UNSEEN * np.ones(npix)
            if k == 'PROB':
                a[inds[bl]] = sky_map['PROBDENSITY'][bl]

            else:
                a[inds[bl]] = sky_map[k][bl]
            if ii == Order:
                bl_ = (a != hp.UNSEEN)
                maps[k][bl_] += a[bl_]
                del a
            else:
                a_ = hp.ud_grade(a, nside_out=Nside,
                                 order_in='Nested', order_out='Nested')
                bl_ = (a_ != hp.UNSEEN)
                maps[k][bl_] += a_[bl_]
                del a, a_

    maps['PROB'] = maps['PROB'] * \
        (np.pi / 180) ** 2 * hp.nside2pixarea(Nside, degrees=True)
    # print('Total probability is:', maps['PROB'].sum())
    return maps


def uniq2order_ind(uniq):
    order = (np.log2(uniq / 4).astype(int) / 2).astype(int)
    inds = (uniq - 4 * (np.power(4, order))).astype(int)
    return order, inds


def order_inds2uniq(order, inds):
    uniq = 4 * (np.power(4, order)).astype(int) + inds
    return uniq


def LoadHealpixMap(thisfilename):
    """
    Bottom-level function that downloads aLIGO HEALpix map and keep in cache. 

    :param thisfilename: name of the fits file containing the localisation map
    :type thisfilename: str

    :return tprob, tdistmu, tdistsigma, distnorm, detectors, event_id, distmean, disterr
    :rtype: array, array, array, array, array, str, float, foat, float, 
    """
    '''
   :return tprob : array of p-values as a function of sky position
    :return tdistmu : array of distance estimate
    :return tdistsigma : array of error on distance estimates
    :return distnorm : array of distance normalisations
    :return detectors: which interferometers triggered
    :return event_id: ID of the event
    :return distmean: mean distance from the header
    :return disterr: error on distance from the header
    :rtype: array
    '''
    PrintFileName = "Loading LVC HEALPix map from file: " + thisfilename
    print(PrintFileName)
    fitsfile = fits.open(thisfilename)

    tevent_id = "Non specified"
    tdetectors = ""
    tdistmean = 0
    tdisterr = 0
    tdistmu = []
    tdistsigma = []
    tdistnorm = []

    if 'OBJECT' in fitsfile[1].header:
        tevent_id = fitsfile[1].header['OBJECT']
    else:
        tevent_id = "Non specified"

    if 'INSTRUME' in fitsfile[1].header:
        tdetectors = fitsfile[1].header['INSTRUME']
    else:
        tdetectors = "Non specified"

    if (fitsfile[1].header['TFIELDS'] <= 2):
        skymap = lf.read_sky_map(thisfilename)  
        tprob = skymap[0]
    else:
        skymap = lf.read_sky_map(thisfilename, distances = True) 
        tprob = skymap[0][0]
        tdistmu = skymap[0][1]
        tdistsigma = skymap[0][2]
        tdistnorm  = skymap[0][3]
        tdistmean = fitsfile[1].header['DISTMEAN']
        tdisterr = fitsfile[1].header['DISTSTD']
        print('Event has triggered ', tdetectors, ' => distance = {0:.2f}'.format(
            tdistmean), ' +- {0:.2f}'.format(tdisterr), ' Mpc') 

    fitsfile.close()

    return tprob, tdistmu, tdistsigma, tdistnorm, tdetectors, tevent_id, tdistmean, tdisterr


def LoadHealpixMap_Flat(thisfilename):
    """
    Bottom-level function that downloads aLIGO HEALpix map and keep in cache. 

    :param thisfilename: name of the fits file containing the localisation map
    :type thisfilename: str

    :return tprob, tdistmu, tdistsigma, distnorm, detectors, event_id, distmean, disterr
    :rtype: array, array, array, array, array, str, float, foat, float, 
    """
    '''
   :return tprob : array of p-values as a function of sky position
    :return tdistmu : array of distance estimate
    :return tdistsigma : array of error on distance estimates
    :return distnorm : array of distance normalisations
    :return detectors: which interferometers triggered
    :return event_id: ID of the event
    :return distmean: mean distance from the header
    :return disterr: error on distance from the header
    :rtype: array
    '''
    PrintFileName = "Loading LVC HEALPix map from file: " + thisfilename
    print(PrintFileName)
    fitsfile = fits.open(thisfilename)

    tevent_id = "Non specified"
    tdetectors = ""
    tdistmean = 0
    tdisterr = 0
    tdistmu = []
    tdistsigma = []
    tdistnorm = []

    if 'OBJECT' in fitsfile[1].header:
        tevent_id = fitsfile[1].header['OBJECT']
    else:
        tevent_id = "Non specified"

    if 'INSTRUME' in fitsfile[1].header:
        tdetectors = fitsfile[1].header['INSTRUME']
    else:
        tdetectors = "Non specified"

    if (fitsfile[1].header['TFIELDS'] == 4):
        tprob, tdistmu, tdistsigma, tdistnorm = hp.read_map(
            thisfilename, field=range(4))
        tdistmean = fitsfile[1].header['DISTMEAN']
        tdisterr = fitsfile[1].header['DISTSTD']
        print('Event has triggered ', tdetectors, ' => distance = {0:.2f}'.format(
            tdistmean), ' +- {0:.2f}'.format(tdisterr), ' Mpc')
    else:
        tprob = hp.read_map(thisfilename, field=range(1))
    # raise

    fitsfile.close()

    return tprob, tdistmu, tdistsigma, tdistnorm, tdetectors, tevent_id, tdistmean, tdisterr


def LoadHealpixUNIQMap(thisfilename):
    '''Download aLIGO HEALpix map and keep in cache
        RETURNS:
        --------

        tprob : array of p-values as a function of sky position
        tdistmu : array of distance estimate
        tdistsigma : array of error on distance estimates
        distnorm : array of distance normalisations
        detectors: which interferometers triggered
        event_id: ID of the event
        distmean: mean distance from the header
        disterr: error on distance from the header
        '''

    tdistmu = []
    tdistsigma = []
    tdistnorm = []

    PrintFileName = "Loading LVC HEALPix UNIQ map from file: " + thisfilename
    print(PrintFileName)
    healpix_skymaps_dict = UNIQSkymap_toNested(thisfilename)

    tprob = healpix_skymaps_dict['PROB']
    tdistmu = healpix_skymaps_dict['DISTMU']
    tdistsigma = healpix_skymaps_dict['DISTSIGMA']
    tdistnorm = healpix_skymaps_dict['DISTNORM']

    return tprob, tdistmu, tdistsigma, tdistnorm


def MOC_confidence_region_Flat(infile, percentage, short_name=' ', save2File=False):

    # reading skymap
    hpx = hp.read_map(infile, verbose=False)
    npix = len(hpx)
    nside = hp.npix2nside(npix)

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
    theta, phi = hp.pix2ang(nside, table_ipix_contour)
    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    from astropy.table import Table
    contour_ipix = Table([ra, dec], names=(
        'RA[deg]', 'DEC[deg]'), meta={'ipix': 'ipix table'})

    # setting MOC order
    from math import log
    moc_order = int(log(nside, 2))

    # creating a MOC map from the contour_ipix table
    moc = MOC.from_table(contour_ipix, 'RA[deg]', 'DEC[deg]', moc_order)

    # writing MOC file in fits
    if (save2File):
        moc.write(short_name + '_MOC_' + str(percentage), format='fits')
    return moc


def MOC_confidence_region2D_Flat(hpx, percentage, short_name=' ', save2File=False):
    """
        Multi-Order coverage map (MOC) of sky area enclosed within a contour plot
        at a given confidence level.

        Input:
        infile: healpix format
        LVC probability sky map
        percentage: float
        probability percentage of the enclosed area
        short_name: str
        output file name

        Output: fits format
        MOC map named "short_name"_"percentage"

        Remark: for json format change the statement
        "moc.write(short_name+'_MOC_'+str(percentage), format='fits' )" -->
        "moc.write(short_name+'_MOC_'+str(percentage), format='json' )"
        """

    # reading skymap
    npix = len(hpx)
    nside = hp.npix2nside(npix)

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
    theta, phi = hp.pix2ang(nside, table_ipix_contour)
    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    from astropy.table import Table
    contour_ipix = Table([ra, dec], names=(
        'RA[deg]', 'DEC[deg]'), meta={'ipix': 'ipix table'})

    # setting MOC order
    from math import log
    moc_order = int(log(nside, 2))

    # creating a MOC map from the contour_ipix table
    moc = MOC.from_table(contour_ipix, 'RA[deg]', 'DEC[deg]', moc_order)

    # writing MOC file in fits
    if (save2File):
        moc.write(short_name + '_MOC_' + str(percentage), format='fits')
    return moc


def IsMultiOrder(fields):
    isMO = True
    if fields==1:
        isMO = False
    if fields == 4:
        isMO = False
    return isMO

def Intersect2D(filename,intersectionThres,obspar):

    distCut = obspar.distCut
    distnorm = []
    tdistmean = 0
    tdiststd = 0
    fitsfile = fits.open(filename)
    has3D = True
    skymap = lf.read_sky_map(filename)  
    prob = skymap[0]
    
    #Check if the skymap has 3D information
    obspar.MO = IsMultiOrder(fitsfile[1].header['TFIELDS'])
    if (fitsfile[1].header['TFIELDS'] <= 2):
        has3D = False
    else:
        tdistmean = fitsfile[1].header['DISTMEAN']
        tdiststd= fitsfile[1].header['DISTSTD']
        # Check if the object is too far away to use a catalog
        if tdistmean+2*tdiststd > distCut:
            has3D = False
    
    nside = hp.npix2nside(len(prob))
    npix = hp.nside2npix(nside)

    theta, phi = hp.pix2ang(nside, np.arange(npix))

    totprob = 0
    ver = True

    for i in range(len(theta)):
        T = np.degrees(theta[i])
        P = np.degrees(phi[i])

        dec = np.radians(90 - T)
        ra =  np.radians(P)

        decngp = np.radians(27.13)
        rangp = np.radians(192.85)
        lncp = np.radians(122.93314)

        sinb = np.sin(decngp) * np.sin(dec) + np.cos(decngp)*np.cos(dec)*np.cos(ra - rangp)
        b = np.degrees(np.arcsin(sinb))
        l = np.degrees(-np.arcsin( np.cos(dec)*np.sin(ra - rangp) / np.cos(b) ) +lncp)

        if abs(b) <= 5:
            totprob += prob[i]

    if 1 >= totprob >= intersectionThres/100:
        ver = False
        return has3D, nside, totprob, ver #a big part of the GW falls behind the galactic plane, so we need to use the 3D method
    elif intersectionThres/100 > totprob >= 0:
        ver = True
        return has3D, nside, totprob, ver #the percentage of GW behind the galactic plane is still low enough to apply the 2D method

def Check2Dor3D(fitsfile, filename, obspar):
    
    distCut = obspar.distCut
    distnorm = []
    tdistmean = 0
    tdiststd = 0
    fitsfile = fits.open(filename)
    has3D = True
    skymap = lf.read_sky_map(filename)  
    prob = skymap[0]

    #Check if the skymap has 3D information
    obspar.MO = IsMultiOrder(fitsfile[1].header['TFIELDS'])
    if (fitsfile[1].header['TFIELDS'] <= 2):
        has3D = False
    else:
        tdistmean = fitsfile[1].header['DISTMEAN']
        tdiststd= fitsfile[1].header['DISTSTD']
        # Check if the object is too far away to use a catalog
        if tdistmean+2*tdiststd > distCut:
            has3D = False

    # Check if no galaxy catalog was given as has3D should be False
    if obspar.galcatName == None:
        has3D = False

    # Check if the algorithm type is declared in the config file 
    if obspar.algorithm != None: 
        if has3D == True:
            if obspar.algorithm == '2D':
                has3D = False 
            
    # Check if the hotspot is in the galactic plane
    npix = len(prob)
    NSide = hp.npix2nside(npix)
    MaxPix = np.argmax(prob)
    MaxTheta, MaxPhi = hp.pix2ang(NSide, MaxPix)
    raMax = np.rad2deg(MaxPhi)
    decMax = np.rad2deg(0.5 * np.pi - MaxTheta)
    c_icrs = SkyCoord(raMax, decMax, frame='fk5', unit=(u.deg, u.deg))

    InsidePlane = Tools.GalacticPlaneBorder(c_icrs)
    print('Is the hotspot in the galactic plane?',InsidePlane)
    if InsidePlane:
        has3D = False 
    fitsfile.close()
    return prob, has3D, NSide


def Check2Dor3D_Flat(fitsfile, filename, obspar):

    distCut = obspar.distCut
    distnorm = []
    tdistmean = 0
    tdiststd = 0
    fitsfile = fits.open(filename)
    if (fitsfile[1].header['TFIELDS'] == 4):
        prob, distmu, distsigma, distnorm = hp.read_map(filename,
                                                        field=range(4))
        tdistmean = fitsfile[1].header['DISTMEAN']
        tdiststd= fitsfile[1].header['DISTSTD']
    else:
        prob = hp.read_map(fitsfile, field=range(1))

    has3D = True
    if len(distnorm) == 0:
        has3D = False

    # Check if no galaxy catalog was given as has3D should be False
    if obspar.galcatName =='False':
        has3D = False
    
    # The distance fullfils the catalog cut
    if tdistmean+2*tdiststd > distCut:
        has3D = False

    # Check if the algorithm type is declared in the config file 
    if obspar.algorithm != None: 
        if has3D == True:
            if obspar.algorithm == '2D':
                has3D = False 
    npix = len(prob)
    NSide = hp.npix2nside(npix)
    if obspar.algorithm == None:  
        MaxPix = np.argmax(prob)
        MaxTheta, MaxPhi = hp.pix2ang(NSide, MaxPix)
        raMax = np.rad2deg(MaxPhi)
        decMax = np.rad2deg(0.5 * np.pi - MaxTheta)
        c_icrs = SkyCoord(raMax, decMax, frame='fk5', unit=(u.deg, u.deg))

        InsidePlane = Tools.GalacticPlaneBorder(c_icrs)
        print('Is the hotspot in the galactic plane?',InsidePlane)
        if InsidePlane:
            has3D = False
    return prob, has3D, NSide

######################################################

def NightDarkObservation(time, obspar):
    '''
    Function that searches for an array of observation times that fulfilled darkness condition and window

    '''
    obs = Observatory(longitude=obspar.lon.to_value(u.deg),
                      latitude=obspar.lat.to_value(u.deg),
                      elevation=obspar.height.to_value(u.m),
                      run_duration=datetime.timedelta(minutes=obspar.duration),
                      minimal_run_duration=datetime.timedelta(minutes=obspar.minDuration),
                      max_sun_altitude=obspar.sunDown,
                      max_moon_altitude=obspar.moonDown,
                      max_moon_phase=-1.0)
    return obs.get_time_window(start_time=time,
                               nb_observation_night=obspar.maxNights)


def NightDarkObservationwithGreyTime(time, obspar):
    '''
    Function that searches for an array of observation times that fulfilled darkness condition and window

    '''
    obs = Observatory(longitude=obspar.lon.to_value(u.deg),
                      latitude=obspar.lat.to_value(u.deg),
                      elevation=obspar.height.to_value(u.m),
                      run_duration=datetime.timedelta(minutes=obspar.duration),
                      minimal_run_duration=datetime.timedelta(minutes=obspar.minDuration),
                      max_sun_altitude=obspar.sunDown,
                      max_moon_altitude=obspar.moonDown,
                      max_moon_phase=obspar.moonPhase/100.)
    return obs.get_time_window(start_time=time,
                               nb_observation_night=obspar.maxNights)


def ZenithAngleCut(prob, nside, time, minProbcut, maxZenith, observatory, minMoonSourceSeparation, useGreytime):
    '''
    Mask in the pixels with zenith angle larger than 45
    '''
    # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)
    frame = co.AltAz(obstime=time, location=observatory)
    pprob = prob

    mzenith = hp.ma(pprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    altaz_map = targetCoord_map.transform_to(frame)
    maskzenith[altaz_map.alt.value < 90-maxZenith] = 1
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
        moonaltazs = get_moon(Time(time, scale='utc')).transform_to(AltAz(obstime=Time(time, scale='utc'),
                                                                          location=observatory))
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


def ComputeProbability2D(prob, highres, radecs, reducedNside, HRnside, minProbcut, time, observatory, maxZenith, FOV, tname, ipixlist, ipixlistHR, counter, dirName, useGreytime, plot):
    '''
    Compute probability in 2D by taking the highest probability in FoV value
    '''
    radius = FOV
    frame = co.AltAz(obstime=time, location=observatory)
    thisaltaz = radecs.transform_to(frame)
    # pix_alt1 = thisaltaz.alt.value

    if useGreytime:
        moonaltazs = get_moon(Time(time, scale='utc')).transform_to(
            AltAz(obstime=Time(time, scale='utc'), location=observatory))
        # Zenith and Moon angular distance mask
        pix_ra = radecs.ra.value[(thisaltaz.alt.value > 90-maxZenith)
                                 & (thisaltaz.separation(moonaltazs) > 30 * u.deg)]
        pix_dec = radecs.dec.value[(thisaltaz.alt.value > 90 - maxZenith)
                                   & (thisaltaz.separation(moonaltazs) > 30 * u.deg)]

    else:
        # Zenith angle mask
        pix_ra = radecs.ra.value[(thisaltaz.alt.value > 90-maxZenith)]
        pix_dec = radecs.dec.value[thisaltaz.alt.value > 90 - maxZenith]
        # pix_alt = pix_alt1[thisaltaz.alt.value > 90 - maxZenith]

    phipix = np.deg2rad(pix_ra)
    thetapix = 0.5 * np.pi - np.deg2rad(pix_dec)

    ipix = hp.ang2pix(reducedNside, thetapix, phipix)

    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)

    cat_pix = Table([ipix, pix_ra, pix_dec, dp_Pix_Fov],
                    names=('PIX', 'PIXRA', 'PIXDEC', 'PIXFOVPROB'))

    dp_dV_FOV = []
    # dp_dV_FOV = np.zero(len(dp_Pix_Fov))

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
            m_ipix_discfull = ma.compressed(ma.masked_array(
                ipix_discfull, mask=np.logical_not(maskComputeProb)))
            HRprob = highres[m_ipix_discfull].sum()
            # HRprob = 0
            # for j in ipix_discfullNotCovered:
            #    HRprob = HRprob+highres[j]
            # print('Length of list of pixels:', m_ipix_discfull, 'vs', ipix_discfull, 'vs', ipixlistHR)
            # print('Comparison to see if mask is considered: ',HRprob, 'vs',highres[ipix_discfull].sum())
        dp_dV_FOV.append(HRprob)
    cat_pix['PIXFOVPROB'] = dp_dV_FOV

    # Mask already observed pixels
    mask = np.isin(cat_pix['PIX'], ipixlist, invert=True)
    if all(np.isin(cat_pix['PIX'], ipixlist, invert=False)):
        maskcat_pix = cat_pix
    else:
        maskcat_pix = cat_pix[mask]

    # Sort table
    sortcat = maskcat_pix[np.flipud(np.argsort(maskcat_pix['PIXFOVPROB']))]
    # Chose highest

    targetCoord = co.SkyCoord(
        sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))

    P_GW = sortcat['PIXFOVPROB'][:1]

    # Include to the list of pixels already observed

    if (P_GW >= minProbcut):
        phip = float(np.deg2rad(targetCoord.ra.deg))
        thetap = float(0.5 * np.pi - np.deg2rad(targetCoord.dec.deg))
        xyz = hp.ang2vec(thetap, phip)

        ipixlistHR.extend(hp.query_disc(HRnside, xyz, np.deg2rad(radius)))
        ipix_disc = hp.query_disc(reducedNside, xyz, np.deg2rad(radius))
        ipixlist.extend(ipix_disc)

        ##################################
        # PLOT THE RESULTS
        if plot:
            path = dirName + '/EvolutionPlot'
            if not os.path.exists(path):
                os.mkdir(path, 493)
            # nside = 1024

            # hp.mollview(highres,title="With FoV circle")

            hp.gnomview(prob, xsize=500, ysize=500, rot=[
                        targetCoord.ra.deg, targetCoord.dec.deg], reso=8.0)
            hp.graticule()

            ipix_discplot = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
            tt, pp = hp.pix2ang(HRnside, ipix_discplot)
            ra2 = np.rad2deg(pp)
            dec2 = np.rad2deg(0.5 * np.pi - tt)
            skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))

            # hp.visufunc.projplot(skycoord.ra, skycoord.dec, 'y.', lonlat=True, coord="C")
            # plt.show()
            # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

            hp.visufunc.projplot(
                sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], 'r.', lonlat=True, coord="C")
            MaxCoord = SkyCoord(
                sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))
            separations = skycoord.separation(MaxCoord)
            tempmask = separations < (radius + 0.05 * radius) * u.deg
            tempmask2 = separations > (radius - 0.05 * radius) * u.deg
            hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask &
                                 tempmask2].dec, 'r.', lonlat=True, coord="C", linewidth=0.1)
            altcoord = np.empty(1000)
            azcoord = np.random.rand(1000) * 360

            plt.savefig('%s/Zoom_Pointing_%g.png' % (path, counter))
            # for i in range(0,1):
            #    altcoord.fill(90-(maxZenith-5*i))
            #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
            #    RandomCoord_radec = RandomCoord.transform_to('fk5')
            #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
            # plt.show()
            # plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))

    return P_GW, targetCoord, ipixlist, ipixlistHR


def SubstractPointings2D(tpointingFile, prob, nside, FOV, pixlist):
    radius = FOV
    print("Loading pointings from " + tpointingFile)
    ra, dec = np.genfromtxt(tpointingFile, usecols=(2, 3), dtype="str", skip_header=1,
                            delimiter=' ',
                            unpack=True)  # ra, dec in degrees
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    coordinates = TransformRADec(ra, dec)
    P_GW = []
    for i in range(0, len(ra)):
        t = 0.5 * np.pi - coordinates[i].dec.rad
        p = coordinates[i].ra.rad
        xyz = hp.ang2vec(t, p)
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
        effectiveipix_disc = []
        for j in range(0, len(ipix_disc)):
            if not (ipix_disc[j] in pixlist):
                effectiveipix_disc.append(ipix_disc[j])
            pixlist.append(ipix_disc[j])
        P_GW.append(prob[effectiveipix_disc].sum())
        print('Coordinates ra:', ra[i], 'dec:', dec[i],
              'Pgw:', P_GW[i], 'vs', prob[ipix_disc].sum())
    return pixlist, np.sum(P_GW)

######################################################


def TransformRADec(vra, vdec):
    if ('h' in vra[0]):
        ra = []
        dec = []
        for i in range(0, len(vra)):
            coord = SkyCoord(vra[i].split('"')[1],
                             vdec[i].split('"')[0], frame='fk5')
            # print(coord)
            ra.append(coord.ra.deg)
            dec.append(coord.dec.deg)
    else:
        # print(vra,vdec)
        ra = vra.astype(float)
        dec = vdec.astype(float)
        # float(vra)
        # dec = float(vdec)
    # print(ra,dec)
    coordinates = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    return coordinates

# Extra functions from BestCandidateon PGal ==> 3D

######################################################


def LoadGalaxies(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    # Load data
    h5file = tables.open_file(tgalFile, mode="r")
    tcat = Table(h5file.root.catalog.read()[
                 ['no_GLADE', 'RA', 'Dec', 'd_L', 'B_mag']])
    h5file.close()

    # Rename column to match naming scheme
    tcat.rename_columns(['RA', 'Dec', 'd_L', 'B_mag'], [
                        'RAJ2000', 'DEJ2000', 'Dist', 'Bmag'])

    return tcat


def LoadGalaxies_SteMgal(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    # Load data
    h5file = tables.open_file(tgalFile, mode="r")
    tcat = Table(h5file.root.catalog.read()[
                 ['no_GLADE', 'RA', 'Dec', 'd_L', 'B_mag', 'mass']])
    h5file.close()

    # Rename column to match naming scheme
    tcat.rename_columns(['RA', 'Dec', 'd_L', 'B_mag', 'mass'], [
                        'RAJ2000', 'DEJ2000', 'Dist', 'Bmag', 'SteMgal'])

    return tcat


def CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, Info3D_available, MinimumProbCutForCatalogue):
    '''
    Correlates galaxies with GW 3D information following Going the Distance, then sort catalog by that value
    In case there is no 3 extra layers, the probability in 2D is asigned.
    '''
    ra = cat['RAJ2000']
    dec = cat['DEJ2000']
    dist = cat['Dist']

    # Translate RA,Dec of galaxies into theta,phi angles
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)

    # Get corresponding healpix pixel IDs
    npix = len(prob)
    nside = hp.npix2nside(npix)
    ipix = hp.ang2pix(nside, theta, phi)

    # Calculate probability in the space volumes

    pixarea = hp.nside2pixarea(nside)

    if (Info3D_available):
        dp_dV = prob[ipix] * distnorm[ipix] * \
            norm(distmu[ipix], distsigma[ipix]).pdf(dist) / pixarea

    else:
        dp_dV = prob[ipix] / pixarea

    # Add dp_dV to catalogue
    cat['dp_dV'] = dp_dV

    # Select all values > 1% of the peak prob.

    total_dP_dV = dp_dV.sum()
    min_prob_cut = dp_dV > MinimumProbCutForCatalogue * max(dp_dV)
    Gals = cat[min_prob_cut]
    # return array with list of Galaxies passing cuts, ordered by p-value

    tGals = Gals[np.flipud(np.argsort(Gals['dp_dV']))]
    # ascii.write(tGals, '/Users/hashkar/Desktop/GWfollowup/GW-Followup/tGals_noM.txt', names = ['RAJ2000','DEJ2000','Dist','Bmag','dp_dV'],overwrite=True)
    return tGals, total_dP_dV


def CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, distmean, disterr, distnorm, cat, Info3D_available, MinimumProbCutForCatalogue):
    '''
    Correlates galaxies with GW 3D information following Going the Distance, then sort catalog by that value
    In case there is no 3 extra layers, the probability in 2D is asigned.
    '''
    beta = 1
    alpha = 0
    ra = cat['RAJ2000']
    dec = cat['DEJ2000']
    dist = cat['Dist']

    # Translate RA,Dec of galaxies into theta,phi angles

    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)

    # Get corresponding healpix pixel IDs
    npix = len(prob)
    nside = hp.npix2nside(npix)
    ipix = hp.ang2pix(nside, theta, phi)

    # get the pixles that are confined in the 90% area
    pixtab = Get90RegionPixGal(prob, 0.9, nside)

    # create new catalog with galaxies inside 90% region
    cat['pixtab'] = ipix
    pixtablist = np.in1d(ipix, pixtab)
    Gals = cat[pixtablist]

    # filter out galaxies outside the distance uncertaintity
    dist = Gals['Dist']
    distok = (dist < (distmean+2*disterr)) & (dist > (distmean-2*disterr))
    Gals = Gals[distok]

    # compute a new ipix
    ra = Gals['RAJ2000']
    dec = Gals['DEJ2000']
    dist = Gals['Dist']

    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)

    NewIpix = hp.ang2pix(nside, theta, phi)

    # Calculate probability in the space volumes
    pixarea = hp.nside2pixarea(nside)

    if (Info3D_available):
        dp_dV_pos = prob[NewIpix] * distnorm[NewIpix] * \
            norm(distmu[NewIpix], distsigma[NewIpix]).pdf(dist) / pixarea

    else:
        dp_dV_pos = prob[NewIpix] / pixarea

    # Add dp_dV to new catalogue
    Gals['dp_dV'] = dp_dV_pos

    # Select all values > 1% of the peak prob.

    total_dP_dV = dp_dV_pos.sum()
    min_prob_cut = dp_dV_pos > MinimumProbCutForCatalogue * max(dp_dV_pos)
    Gals = Gals[min_prob_cut]

    if (Info3D_available):
        Mgal1 = Gals['SteMgal']
        Pgal_pos = Gals['dp_dV']

        Mgal1 = np.nan_to_num(Mgal1)
        Mgal = 10**(Mgal1)
        Pgal_pos = np.nan_to_num(Pgal_pos)

        Gmass = Mgal/(np.sum(Mgal))
        alpha = (Pgal_pos).sum()/(Pgal_pos*Gmass).sum()
        dp_dV = (Pgal_pos)+(Pgal_pos*(alpha*beta*Gmass))

    Gals['dp_dV'] = dp_dV

    total_dP_dV = dp_dV.sum()
    # print(total_dP_dV)

    tGals = Gals[np.flipud(np.argsort(Gals['dp_dV']))]
    # ascii.write(tGals, '/Users/hashkar/Desktop/GWfollowup/GW-Followup/tGals.txt', names = ['RAJ2000','DEJ2000','Dist','Bmag','SteMgal', 'index' ,'dp_dV'],overwrite=True)

    return tGals, total_dP_dV


def VisibleAtTime(test_time, galaxies, maxz, observatory):
    '''Determine if prompt or afterglow follow-up is possible by knowing if there are galaxies with non-negligible probability of hosting the NSM in the FoV
     1) check if any galaxy is visible, if not --> AFTERGLOW

    2) loop over zenith angle and select subsets of galaxies

    3) stop if maximum p-value of this subset is smaller than 75% of the previous subset

    4) else: stricter cut on zenith and repeat

    5) take galaxy with highest p-value fulfilling both criteria as target

    RETURNS:
    --------
    bool `is_vis` : is visible now?
    np.ndarray `alt_az` : alt_az location of galaxies
    '''

    # print()
    # print("Check visibility at time {0}".format(test_time))

    # observatory time and location to look up visibility of objects

    # observatory = co.EarthLocation(lat=-23.271333 * u.deg,lon=16.5 * u.deg, height=1800 * u.m)

    frame = co.AltAz(obstime=test_time, location=observatory)
    # print('galaxies',galaxies)
    # print('galaxies',len(galaxies['RAJ2000']))
    radecs = co.SkyCoord(
        galaxies['RAJ2000'], galaxies['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    if (len(radecs) > 0):
        thisaltaz = radecs.transform_to(frame)

        # add altitude to topGals array
        # already sorted by descending probability value

        galaxies['Alt'] = thisaltaz.alt.value

        # check if any galaxy is visible at the moment
        thismask = thisaltaz.alt.value > 90 - maxz

        nGals = len(galaxies[thismask])

        # print('nGals',nGals)

        if (nGals == 0):
            # print("No galaxies visible within {0} deg zenith angle --> AFTERGLOW".format(maxz))

            return False, thisaltaz, galaxies
        else:
            # print("{0} galaxies are visible within {1} deg zenith angle ""--> Test for prompt follow up".format(nGals, maxz))

            return True, thisaltaz, galaxies
    else:
        thisaltaz = []
        return False, thisaltaz, galaxies


def FulfillsRequirement(theseGals, maxz, FOV, zenithWeighting, UsePix):
    '''
    Apply filter criteria to visible galaxy sample and compares them to get the best option of zenith angle

    '''

    # print("Check if galaxies with minimum p-value are among candidates...")

    # Initialise maximum p-value in map to 1

    maxp = 1
    mask = 0
    thisminz = 0

    alt = theseGals['Alt']

    for minz_aux in range(maxz, 5, -5):

        tmpmask = alt > 90 - (minz_aux)
        # print('TheseGals', theseGals[tmpmask])
        tmpGals = theseGals.copy()
        # print('len(tmpGals[tmpmask])', len(tmpGals[tmpmask]), 'without mask', len(tmpGals))
        # cut on zenith angle and select most probable galaxy

        if (len(tmpGals[tmpmask]) > 0):

            cur_maxp = tmpGals[tmpmask]['dp_dV'].max() / \
                theseGals['dp_dV'].max()

            # print("{0} galaxies visible with zen < {1} deg - maximum p-value {2:0.3f}"

            #     .format(len(tmpGals[tmpmask]), minz_aux, cur_maxp))

           # print("Maximum probability of {0:0.3f} of global maximum of {1:3f}"

            #      .format(cur_maxp, theseGals['dp_dV'].max()))

            # define final mask

            if (maxz == minz_aux):
                maxp = cur_maxp
                mask = tmpmask
                thisminz = minz_aux

            if (cur_maxp > zenithWeighting * maxp):
                mask = tmpmask
                thisminz = minz_aux
            else:
                thisminz = minz_aux + 5
                break
    if UsePix:
        mask = alt > 90-(thisminz+FOV)
    return mask, thisminz


def FulfillsRequirementGreyObservations(Ktime, theseGals, observatory, minMoonSourceSeparation):

    targetCoord = co.SkyCoord(
        theseGals['RAJ2000'], theseGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    frame = co.AltAz(obstime=Ktime, location=observatory)
    moonaltazs = get_moon(Time(Ktime, scale='utc')).transform_to(frame)

    altaz_map = targetCoord.transform_to(frame)
    separations = altaz_map.separation(moonaltazs)

    # Mask
    greymask = separations > minMoonSourceSeparation*u.deg
    return greymask


def FulfillsRequirement_MinProb(thisGals_aux, maxz):
    ''' Same as FulfillsRequirements but at the end a supplementary cut is performed.
    This algorithm comes from the first version of the code but in the newer versions
    the algorithms separates this two options.

    '''

    print("Check if galaxies with minimum p-value are among candidates...")

    # Initialise maximum p-value in map to 1

    maxp = 1
    mask = 0

    for minz_aux in range(maxz, 5, -5):

        tmpmask = altaz.alt.value > 90 - minz_aux
        tmpGals = thisGals_aux.copy()

        # cut on zenith angle and select most probable galaxy

        if (len(tmpGals[tmpmask]) > 0):

            cur_maxp = tmpGals[tmpmask]['dp_dV'].max() / tGals['dp_dV'].max()

            print("{0} galaxies visible with zen < {1} deg - maximum p-value {2:0.3f}".format(len(tmpGals[tmpmask]),
                                                                                              minz_aux, cur_maxp))
            print("Maximum probability of {0:0.3f} of global maximum of {1:3f}".format(
                cur_maxp, tGals['dp_dV'].max()))

            if (maxz == minz_aux):
                maxp = cur_maxp
                mask = tmpmask
                minz = minz_aux

            if (cur_maxp > 0.75 * maxp):
                mask = tmpmask
                minz = minz_aux
            else:
                minz = minz_aux + 5
                break

    if (maxp < 0.02 * thisGals_aux['dp_dV'].max()):

        print("Probability too low, postpone observations --> AFTERGLOW")
        return False, mask, minz

    else:
        print('This minz= ', minz)
        return True, mask, minz


def Afterglow():
    print('Afterglow!')


def ObtainHighestProbabilityCoordinates(filename):
    hpx = UNIQSkymap_toNested(filename)
    hpx_prob = hpx['PROB']
    print(sum(hpx_prob))
    ipix_max = np.argmax(hpx_prob)
    nside = hp.npix2nside(len(hpx_prob))
    theta, phi = hp.pix2ang(nside, ipix_max)
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    return ra, dec


def SelectObservatory_fromHotspot(filename):
    ra, dec = ObtainHighestProbabilityCoordinates(filename)
    if dec > 0:
        UseObs = 'North'
    else:
        UseObs = 'South'
    return UseObs


def ComputeProbBCFOVSimple(prob, time, observatory, visiGals, allGals, tsum_dP_dV, nside, thisminz, maxZenith, FOV, tname,
                           tsavedcircle, dirName, doPlot):
    '''Computes probability pgal and pgw in FoV and draws everything

    bool doPlot when  = True is used to plot the maps

    RETURNS:

    --------

        P_Gal: Probability of galaxies within FoV in the LIGO signal region
        P_GW: Total probability within FoV of the Ligo signal.
        noncircleGal: Table of galaxies that are outside the circle(s) and inside the LIGO signal region


    '''

    targetCoord = co.SkyCoord(
        visiGals['RAJ2000'][:1], visiGals['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(
        visiGals['RAJ2000'], visiGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord3 = co.SkyCoord(
        allGals['RAJ2000'], allGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    dp_dVfinal = visiGals['dp_dV']
    # dp_dV = tGals['dp_dV']

    # Array of indices of pixels inside circle of FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoord[0].dec.rad

    p = targetCoord[0].ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg',t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg)

    xyz = hp.ang2vec(t, p)

    # print(xyz)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    # print('ipix_disc',ipix_disc)
    P_GW = prob[ipix_disc].sum()

    P_Gal = dp_dVfinal[targetCoord2.separation(
        targetCoord).deg < radius].sum() / tsum_dP_dV

    # print("Total probability within H.E.S.S. FoV: {0}".format(P_GW))

    # print("Probability of galaxies within H.E.S.S. FoV in the LIGO signal region:{0}".format(P_Gal))

    # all galaxies inside the current observation circle

    circleGal = visiGals[targetCoord2.separation(targetCoord).deg < radius]
    # print('Galaxies within the FoV: ', len(circleGal['RAJ2000']))

    # all galaxies outside the current observation circle, no visibility selection

    noncircleGal = allGals[targetCoord3.separation(targetCoord).deg > radius]

    if (doPlot):
        path = dirName + '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)
        hp.mollview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (
            time.day, time.month, time.year, time.hour, time.minute))
        hp.graticule()
        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)
        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))

        frame = co.AltAz(obstime=time, location=observatory.location)
        altaz_all = skycoord.transform_to(frame)
        tmask = altaz_all.alt.value > 90 - thisminz

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        # hp.gnomview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (time.day, time.month, time.year, time.hour, time.minute),xsize = 4000,ysize=6000,rot = [90,-50],reso=0.8)
        # plt.show()
        # plt.savefig("Figures/ExampleGW_%g.png" % (j))

        # draw all galaxies within zenith-angle cut
        # hp.visufunc.projscatter(finalGals['RAJ2000'], finalGals['DEJ2000'], lonlat=True, marker='*', color='g')
        # plt.savefig("Figures/ExampleGW_Galaxies_%g.png" % (j))

        # If I want to plot all gals, plot also the ones that are out of the circle
        # hp.visufunc.projscatter(noncircleGal['RAJ2000'], noncircleGal['DEJ2000'], lonlat=True, marker='*', color='g')

        # draw observation position, which is equivalent to galaxy with highest
        # probability
        hp.visufunc.projscatter(visiGals['RAJ2000'][:1], visiGals['DEJ2000'][:1], lonlat=True, marker='.', color='r',
                                linewidth=0.1)

        # draw circle of FoV around best fit position
        # hp.visufunc.projscatter(allGals['RAJ2000'], allGals['DEJ2000'], lonlat=True, marker='.', color='g',linewidth=0.1)
        # plt.show()
        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'k.', lonlat=True,
                             coord="C", linewidth=0.1)

        # hp.visufunc.projplot(tsavedcircle.ra, tsavedcircle.dec, 'r.', lonlat=True, coord="C",linewidth=0.1)

        # Draw H.E.S.S. visibility

        # altcoord= [np.random.randint(-90,90-thisminz) for _ in range(4000)]

        altcoord = np.empty(4000)
        altcoord.fill(90 - maxZenith)
        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(
            u.deg, u.deg), obstime=time, location=observatory)
        RandomCoord_radec = RandomCoord.transform_to('fk5')
        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec,
                             'b.', lonlat=True, coord="C", linewidth=0.1)
        # plt.show()
        # Draw MinZ area

        #print('Min Zenith= ', thisminz)

        altcoordmin = np.empty(4000)
        azcoordmin = np.random.rand(4000) * 360

        RandomCoordmin = SkyCoord(azcoordmin, altcoordmin, frame='altaz', unit=(
            u.deg, u.deg), obstime=time, location=observatory)
        RandomCoordmin_radec = RandomCoordmin.transform_to('fk5')

        # hp.visufunc.projplot(RandomCoordmin_radec.ra, RandomCoordmin_radec.dec, 'y.', lonlat=True, coord="C", marker='.', markersize = 8 )

        # plt.show()
        plt.savefig("%s/Pointing_%g.png" % (path, len(ObservationTimearray)))

    return P_Gal, P_GW, talreadysumipixarray2


def ComputeProbGalTargetted(prob, time, finalGals, visiGals, allGals, tsum_dP_dV, talreadysumipixarray, nside, thisminz,obspar,counter,dirName):
    '''Computes probability Pgal and Pgw in FoV but it takes into account a list of pixels to avoid recounting already observed zones.
    Returns saved circle too (is it really needed? )
    bool doPlot when  = True is used to plot the maps

    RETURNS:

    --------

        P_Gal: Probability of galaxies within FoV in the LIGO signal region
        P_GW: Total probability within  FoV of the Ligo signal.
        noncircleGal: Table of galaxies that are outside the circle(s) and inside the LIGO signal region

    '''
    observatory = obspar.location
    maxZenith = obspar.maxZenith
    FOV = obspar.FOV
    doPlot = obspar.doPlot

    targetCoord = co.SkyCoord(
        finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(
        visiGals['RAJ2000'], visiGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord3 = co.SkyCoord(
        allGals['RAJ2000'], allGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    dp_dVfinal = visiGals['dp_dV']

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

    P_Gal = dp_dVfinal[targetCoord2.separation(
        targetCoord).deg < radius].sum() / tsum_dP_dV
    circleGal = visiGals[targetCoord2.separation(targetCoord).deg < radius]
    noncircleGal = allGals[targetCoord3.separation(targetCoord).deg > radius]

    if (doPlot):

        path = dirName + '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))

        frame = co.AltAz(obstime=time, location=observatory)
        altaz_all = skycoord.transform_to(frame)
        tmask = altaz_all.alt.value > 90 - thisminz

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        # path = os.path.dirname(os.path.realpath(__file__)) + tname
        # if not os.path.exists(path):
        #    os.mkdir(path, 493)

        # hp.mollview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (time.day, time.month, time.year, time.hour, time.minute), xsize=2000)
        hp.gnomview(prob, xsize=500, ysize=500, rot=[
                    targetCoord.ra.deg, targetCoord.dec.deg], reso=5.0)

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

        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'r.', lonlat=True,
                                coord="C")

        # Draw H.E.S.S. visibility

        # altcoord= [np.random.randint(-90,90-thisminz) for _ in range(4000)]

        altcoord = np.empty(4000)

        altcoord.fill(90 - maxZenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                                location=observatory)

        RandomCoord_radec = RandomCoord.transform_to('fk5')

        hp.visufunc.projplot(
            RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # MOON

        # hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # Draw MinZ area

        #print('Min Zenith= ', thisminz)

        altcoordmin = np.empty(4000)
        altcoordmin.fill(90 - thisminz)
        azcoordmin = np.random.rand(4000) * 360

        RandomCoordmin = SkyCoord(azcoordmin, altcoordmin, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                                    location=observatory)
        
        #RandomCoordmin_radec = RandomCoordmin.transform_to('fk5')
        # hp.visufunc.projplot(RandomCoordmin_radec.ra, RandomCoordmin_radec.dec, 'y.', lonlat=True, coord="C")
        plt.savefig("%s/Zoom_Pointing_%g.png" % (path, counter))

    return P_Gal, P_GW, noncircleGal, talreadysumipixarray


def SimpleGWprob(prob, finalGals, talreadysumipixarray, FOV, nside):
    '''Computes probability Pgw in FoV but it takes into account a list of pixels to avoid recounting already observed zones.
    bool doPlot when  = True is used to plot the maps
    '''
    targetCoord = co.SkyCoord(
        finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))
    radius = FOV
    t = 0.5 * np.pi - targetCoord[0].dec.rad
    p = targetCoord[0].ra.rad
    xyz = hp.ang2vec(t, p)
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])

    probability = prob[effectiveipix_disc].sum()
    return probability


def SubstractPointings(tpointingFile, galaxies, talreadysumipixarray, tsum_dP_dV, FOV, prob, nside):

    # targetCoord = co.SkyCoord(galaxies['RAJ2000'], galaxies['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    # Read PointingsFile

    print("Loading pointings from " + tpointingFile)
    rap, decP, = np.genfromtxt(tpointingFile, usecols=(2, 3), dtype="str", skip_header=1,
                               delimiter=' ',
                               unpack=True)  # ra, dec in degrees

    coordinates = TransformRADec(rap, decP)
    ra = coordinates.ra.deg
    dec = coordinates.dec.deg

    PGW = []
    PGAL = []
    updatedGalaxies = galaxies
    if np.isscalar(ra):
        updatedGalaxies, pgwcircle, pgalcircle, talreadysumipixarray = SubstractGalaxiesCircle(
            updatedGalaxies, ra, dec, talreadysumipixarray, tsum_dP_dV, FOV, prob, nside)
        PGW.append(pgwcircle)
        PGAL.append(pgalcircle)
        print('Coordinates ra:', ra, 'dec:', dec,
              'Pgw:', pgwcircle, 'PGAL:', pgalcircle)
    else:
        for i, coord in enumerate(coordinates):
            ra = coord.ra.deg
            dec = coord.dec.deg
            updatedGalaxies, pgwcircle, pgalcircle, talreadysumipixarray = SubstractGalaxiesCircle(
                updatedGalaxies, ra, dec, talreadysumipixarray, tsum_dP_dV, FOV, prob, nside)
            PGW.append(pgwcircle)
            PGAL.append(pgalcircle)
            print('Coordinates ra:', ra, 'dec:', dec,
                  'Pgw:', pgwcircle, 'PGAL:', pgalcircle)
    return ra, dec, updatedGalaxies, PGW, PGAL, talreadysumipixarray


def SubstractGalaxiesCircle(galaux, ra, dec, talreadysumipixarray, tsum_dP_dV, FOV, prob, nside):

    radius = FOV
    coordinates = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))

    targetCoord = co.SkyCoord(
        galaux['RAJ2000'], galaux['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    dp_dVfinal = galaux['dp_dV']

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
    P_Gal = dp_dVfinal[targetCoord.separation(
        coordinates).deg < radius].sum() / tsum_dP_dV

    print('PGW', P_GW, 'P_GAL', P_Gal)

    newgalaxies = galaux[targetCoord.separation(coordinates).deg > radius]

    return newgalaxies, P_GW, P_Gal, talreadysumipixarray
#####################################################

## Extra functions from PGalinFoVOptimised ==> 3D ##

#####################################################


def ComputePGalinFOV(prob, cat, galpix, FOV, totaldPdV, nside, UsePix):
    '''
        Computes probability Pgal in FoV
    '''
    if UsePix:
        targetCoord = co.SkyCoord(
            galpix['PIXRA'], galpix['PIXDEC'], frame='fk5', unit=(u.deg, u.deg))
    else:
        targetCoord = co.SkyCoord(
            galpix['RAJ2000'], galpix['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(
        cat['RAJ2000'], cat['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    # dp_dVfinal = finalGals['dp_dV']
    dp_dV = cat['dp_dV']

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoord.dec.rad

    p = targetCoord.ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)

    xyz = hp.ang2vec(t, p)

    # print(xyz)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    P_GW = prob[ipix_disc].sum()

    Pgal_inFoV = dp_dV[targetCoord2.separation(
        targetCoord).deg <= radius].sum() / totaldPdV

    return Pgal_inFoV


def ModifyCatalogue(prob, cat, FOV, totaldPdV, nside):
    '''
     Computes the integrated Pgal in FoV for a list of calues using Pgal in FoV and sorts the catalog
     using that quantity as a criteria
    '''
    # lengthSG=0.02*len(cat)
    lengthSG = 100
    SelectedGals = cat[:lengthSG]
    dp_dV_FOV = []
    # print('len(cat[RAJ2000])', len(cat['RAJ2000']))
    # print('len(SelectedGals[RAJ2000])',len(SelectedGals['RAJ2000']))
    for l in range(0, len(cat['dp_dV'])):
        if (l < len(SelectedGals['dp_dV'])):
            dp_dV_FOV.append(ComputePGalinFOV(
                prob, cat, SelectedGals[l], FOV, totaldPdV, nside, UsePix=False))
        else:
            dp_dV_FOV.append(0)

    cat['dp_dV_FOV'] = dp_dV_FOV

    tcat = cat[np.flipud(np.argsort(cat['dp_dV_FOV']))]

    return tcat


def ComputeProbPGALIntegrateFoV(prob, time, observatory, centerPoint, UsePix, visiGals, allGalsaftercuts, tsum_dP_dV, talreadysumipixarray, nside,
                                thisminz, obspar, counter, tname, dirName, doPlot):
    '''
        Same as ComputeProbGalTargetted but it does not return circle coordinates.
    '''
    maxZenith = obspar.maxZenith 
    FOV = obspar.FOV
    if UsePix:
        targetCoord = co.SkyCoord(
            centerPoint['PIXRA'][:1], centerPoint['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))

    else:
        targetCoord = co.SkyCoord(
            centerPoint['RAJ2000'][:1], centerPoint['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(
        visiGals['RAJ2000'], visiGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord3 = co.SkyCoord(allGalsaftercuts['RAJ2000'], allGalsaftercuts['DEJ2000'], frame='fk5',
                               unit=(u.deg, u.deg))

    dp_dVfinal = visiGals['dp_dV']

    # Array of indices of pixels inside circle of FoV

    radius = FOV
    # print('FOV',FOV)
    t = 0.5 * np.pi - targetCoord[0].dec.rad

    p = targetCoord[0].ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg)

    xyz = hp.ang2vec(t, p)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    P_GW = prob[effectiveipix_disc].sum()

    P_Gal = dp_dVfinal[targetCoord2.separation(
        targetCoord).deg < radius].sum() / tsum_dP_dV
    # TotalGal_FOV = len(dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius])

    # print("Total probability within H.E.S.S. FoV: {0}".format(P_GW))

    # print("Probability of galaxies within H.E.S.S. FoV in the LIGO signal region:{0}".format(P_Gal))

    # all galaxies inside the current observation circle

    circleGal = visiGals[targetCoord2.separation(targetCoord).deg < radius]
    # print('Galaxies within the FoV: ', len(circleGal['RAJ2000']))

    # all galaxies outside the current observation circle, no visibility selection

    noncircleGal = allGalsaftercuts[targetCoord3.separation(
        targetCoord).deg > radius]

    if (doPlot):

        path = dirName + '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))

        frame = co.AltAz(obstime=time, location=observatory)
        altaz_all = skycoord.transform_to(frame)
        tmask = altaz_all.alt.value > 90 - thisminz

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        # path = os.path.dirname(os.path.realpath(__file__)) + tname
        # if not os.path.exists(path):
        #    os.mkdir(path, 493)

        # hp.mollview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (time.day, time.month, time.year, time.hour, time.minute), xsize=2000)
        hp.gnomview(prob, xsize=500, ysize=500, rot=[
                    targetCoord.ra.deg, targetCoord.dec.deg], reso=5.0)

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

        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'r.', lonlat=True,
                             coord="C")

        # Draw H.E.S.S. visibility

        # altcoord= [np.random.randint(-90,90-thisminz) for _ in range(4000)]

        altcoord = np.empty(4000)

        altcoord.fill(90 - maxZenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                               location=observatory)

        RandomCoord_radec = RandomCoord.transform_to('fk5')

        hp.visufunc.projplot(
            RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # MOON

        # hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # Draw MinZ area

        #print('Min Zenith= ', thisminz)

        altcoordmin = np.empty(4000)

        altcoordmin.fill(90 - thisminz)

        azcoordmin = np.random.rand(4000) * 360

        RandomCoordmin = SkyCoord(azcoordmin, altcoordmin, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                                  location=observatory)

        RandomCoordmin_radec = RandomCoordmin.transform_to('fk5')

        # hp.visufunc.projplot(RandomCoordmin_radec.ra, RandomCoordmin_radec.dec, 'y.', lonlat=True, coord="C")

        # plt.show()
        plt.savefig("%s/Zoom_Pointing_%g.png" % (path, counter))

    return P_Gal, P_GW, noncircleGal, talreadysumipixarray


def randomDate(start, end, prop):
    return strTimeProp(start, end, '%Y-%m-%d %H:%M:%S', prop)


def strTimeProp(start, end, format, prop):
    """Get a time at a proportion of a range of two formatted times.

        start and end should be strings specifying times formated in the
        given format (strftime-style), giving an interval [start, end].
        prop specifies how a proportion of the interval to be taken after
        start.  The returned time will be in the specified format.
        """

    stime = time.mktime(time.strptime(start, format))
    etime = time.mktime(time.strptime(end, format))
    # print etime

    ptime = stime + prop * (etime - stime)
    # print time.strftime(format, time.localtime(ptime))
    # print time.strftime(format, time.gmtime(ptime))

    return time.strftime(format, time.localtime(ptime))

######################################################

## Extra functions from the use of Center of Pixels ##
# as center of pointings, in a 3D treatment

######################################################


def ModifyCataloguePIX(pix_ra1, pix_dec1, test_time, maxz, prob, cat, FOV, totaldPdV, nside, NewNside, minz,
                       observatory):
    # To do:
    # for a faster time:
    # subtract the summed pixels
    # chose only visbile pixels / done
    # adjust the number of interations and pixel

    #####################
    # pprob = hp.pixelfunc.ud_grade(prob, 64)#power = -2
    # pixel_theta, pixel_phi = hp.pix2ang((hp.npix2nside(len(pprob))), np.arange(len(pprob)))

    # pix_ra1 = np.rad2deg(pixel_phi)
    # pix_dec1 = np.rad2deg(0.5 * np.pi - pixel_theta)
    ##################

    # Cuts on azimuth angle  (in the probability region)

    frame = co.AltAz(obstime=test_time, location=observatory)

    radecs = co.SkyCoord(pix_ra1, pix_dec1, frame='fk5', unit=(u.deg, u.deg))
    thisaltaz = radecs.transform_to(frame)

    # pix_alt1 = thisaltaz.alt.value

    pix_ra = radecs.ra.value[thisaltaz.alt.value > 90 - (minz)]
    pix_dec = radecs.dec.value[thisaltaz.alt.value > 90 - (minz)]
    # pix_alt = pix_alt1[thisaltaz.alt.value > 90 - (minz)]

    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)

    cat_pix = Table([pix_ra, pix_dec, dp_Pix_Fov],
                    names=('PIXRA', 'PIXDEC', 'PIXFOVPROB'))

    ##############################################################
    # Possible:  select the pixels that only have a prob > certain value
    #       To do so attribute for each pix its prob: start with nside initial, put in table prob, reduce resolution, make cut...
    # Note : maybe highest PROBFOV pix is not visible ? ? ? check fullfills requirements
    ###############################################################

    dp_dV_FOV = []

    # iteration on chosen pixel to calculate the probability on their field of view using galaxies
    for l in range(0, len(cat_pix)):
        dp_dV_FOV.append(ComputePGalinFOV(
            prob, cat, cat_pix[l], FOV, totaldPdV, nside, UsePix=True))

    cat_pix['PIXFOVPROB'] = dp_dV_FOV

    ttcat = cat_pix[np.flipud(np.argsort(cat_pix['PIXFOVPROB']))]
    return ttcat

def GetAreaSkymap5090(filename):
    skymap = QTable.read(filename)
    skymap.sort('PROBDENSITY', reverse=True)
    level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * skymap['PROBDENSITY']
    cumprob = np.cumsum(prob)
    
    i = cumprob.searchsorted(0.9)
    area_90 = pixel_area[:i].sum()
    area_90_deg =area_90.to_value(u.deg**2)
    
    j = cumprob.searchsorted(0.5)
    area_50 = pixel_area[:j].sum()
    area_50_deg = area_50.to_value(u.deg**2)
    
    return area_50_deg, area_90_deg

def GetAreaSkymap5090_Flat(filename):
    hpx = hp.read_map(f'{filename}')
    npix = len(hpx)
    nside = hp.npix2nside(npix)

    i = np.flipud(np.argsort(hpx))
    sorted_credible_levels = np.cumsum(hpx[i])
    credible_levels = np.empty_like(sorted_credible_levels)
    credible_levels[i] = sorted_credible_levels

    area_50 = np.sum(credible_levels <= 0.5) * hp.nside2pixarea(nside, degrees=True)
    print(f"50% area: {area_50} deg2")

    area_90 = np.sum(credible_levels <= 0.9) * hp.nside2pixarea(nside, degrees=True)
    print(f"90% area: {area_90} deg2")

    return area_50, area_90


def Get90RegionPixReduced(hpxx, percentage, Nnside):

    nside = Nnside  # size of map used for contour determination
    hpx = hp.ud_grade(hpxx, nside_out = nside, power=-2, order_in='Nested', order_out='Nested')

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
    # area = len(table_ipix_contour)*hp.nside2pixarea(nside, True)
    # from index to polar coordinates
    theta1, phi1 = hp.pix2ang(nside, table_ipix_contour)
    area = len(table_ipix_contour)*hp.nside2pixarea(nside, True)

    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    # contour_ipix = Table([ra, dec], names=('RA[deg]', 'DEC[deg]'), meta={'ipix': 'ipix table'})

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

    nside = Nside  # size of map used for contour determination
    # hpx = hp.ud_grade(hpxx, nside, power=-2)
    hpx = hpxx
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
    return table_ipix_contour


######################################################
#####     Functions used in CTA Simulations       ####
######################################################

def GiveProbToGalaxy(prob, cat, distance, Edistance_max, Edistance_min, MinimumProbCutForCatalogue):
    # Cut cat to a cube in distance!

    dist = cat['Dist']
    mask1 = dist < Edistance_max
    mask2 = dist > Edistance_min
    subcat = cat[mask1 & mask2]
    # print(len(subcat),len(cat),Edistance_max,Edistance_min)
    ra = subcat['RAJ2000']
    dec = subcat['DEJ2000']
    # Translate RA,Dec of galaxies into theta,phi angles

    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)

    # Get corresponding healpix pixel IDs

    npix = len(prob)
    nside = hp.npix2nside(npix)
    ipix = hp.ang2pix(nside, theta, phi)

    # Calculate probability in the space volumes

    pixarea = hp.nside2pixarea(nside)

    # Give probability to galaxy
    dp_dV = prob[ipix] / pixarea

    # plt.hist(dp_dV/dp_dV.sum(),histtype='step', stacked=True,  fill=False,cumulative=True,density=True)
    # plt.savefig('/Users/mseglar/Documents/GitHub/CTASimulationsGW/3DComparison_plots/BarbaraConvolution.png')
    TotalProb = dp_dV.sum()
    subcat['dp_dV'] = dp_dV
    min_prob_cut = dp_dV > MinimumProbCutForCatalogue * max(dp_dV)
    Gals = subcat[min_prob_cut]
    # return array with list of Galaxies passing cuts, ordered by p-value

    tGals = Gals[np.flipud(np.argsort(Gals['dp_dV']))]

    return tGals, TotalProb


def LoadGalaxies_GladeCTASimu(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    name, dist, z, ra, dec, flag = np.genfromtxt(tgalFile, usecols=(0, 1, 2, 3, 4, 5), skip_header=3, unpack=True,
                                                 dtype='str')  # ra, dec in degrees
    ra = ra.astype(float)
    dec = dec.astype(float)
    z = z.astype(float)
    dist = dist.astype(float) / 1000  # change to Mpc!
    tcat = Table([name, ra, dec, dist, z, flag], names=(
        'Galaxy', 'RAJ2000', 'DEJ2000', 'Dist', 'z', 'flag'))
    # print(tcat)
    return tcat


def LoadGalaxiesSimulation(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    dist, z, ra, dec = np.genfromtxt(tgalFile, usecols=(
        1, 2, 3, 4), skip_header=3, unpack=True)  # ra, dec in degrees

    tcat = Table([ra, dec, dist], names=('RAJ2000', 'DEJ2000', 'Dist'))
    return tcat


def PointingFileReadCTA(pointingFile):
    time1, time2, RA, Dec, Observatory, ZenIni, ZenEnd, Duration = np.genfromtxt(pointingFile,
                                                                                 usecols=(
                                                                                     0, 1, 2, 3, 4, 6, 7),
                                                                                 skip_header=1, unpack=True,
                                                                                 dtype='str')
    time = []
    RA = RA.astype(float)
    Dec = Dec.astype(float)
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([time, RA, Dec, Observatory, ZenIni, ZenEnd, Duration],
                        names=['Observation Time UTC', 'RA[deg]', 'DEC[deg]', 'Observatory', 'ZenIni[deg]', 'ZenEnd[deg]',
                               'Duration[s]'])
    return OutputTable


def TableImportCTA_simple(inputFileName):
    run, MergerID, RA, Dec, distance, z, theta, ndet, SNR, A90, A50 = np.genfromtxt(inputFileName, usecols=(
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), skip_header=3, unpack=True, dtype='str')
    RA = RA.astype(float)
    Dec = Dec.astype(float)
    z = z.astype(float)
    distance = distance.astype(float)
    theta = theta.astype(float)
    ndet = ndet.astype(float)
    SNR = SNR.astype(float)
    A90 = A90.astype(float)
    A50 = A50.astype(float)
    OutputTable = Table([run, MergerID, RA, Dec, distance, z, theta, ndet, SNR, A90, A50], names=(
        'run', 'MergerID', 'RA', 'Dec', 'Distance', 'redshift', 'theta', 'ndet', 'SNR', 'A90', 'A50'))
    return OutputTable


def TableImportCTA(tgalFile):
    run, MergerID, RA, Dec, distance, distMin, distMax, z, theta, ndet, SNR, A90, A50 = np.genfromtxt(tgalFile,
                                                                                                      usecols=(
                                                                                                          0, 1, 2, 3, 4, 5,
                                                                                                          6, 7, 8, 9, 10,
                                                                                                          11, 12),
                                                                                                      skip_header=1,
                                                                                                      unpack=True,
                                                                                                      dtype='str')
    RA = RA.astype(float)
    Dec = Dec.astype(float)
    z = z.astype(float)
    distance = distance.astype(float) / 1000  # to Mpc!
    distMax = distMax.astype(float) / 1000  # to Mpc!
    distMin = distMin.astype(float) / 1000  # to Mpc!
    theta = theta.astype(float)
    ndet = ndet.astype(float)
    SNR = SNR.astype(float)
    A90 = A90.astype(float)
    A50 = A50.astype(float)
    OutputTable = Table([run, MergerID, RA, Dec, distance, distMin, distMax, z, theta, ndet, SNR, A90, A50], names=(
        'run', 'MergerID', 'RA', 'Dec', 'Distance', 'DistMin', 'DistMax', 'redshift', 'theta', 'ndet', 'SNR', 'A90', 'A50'))
    # print(OutputTable)
    return OutputTable


def TableImportCTA_Glade(tgalFile):
    run, gal, MergerID, RA, Dec, distance, z, theta, ndet, SNR, A90, A50 = np.genfromtxt(tgalFile, usecols=(
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), skip_header=2, unpack=True, dtype='str')
    RA = RA.astype(float)
    Dec = Dec.astype(float)
    z = z.astype(float)
    distance = distance.astype(float) / 1000  # to Mpc!
    theta = theta.astype(float)
    ndet = ndet.astype(float)
    SNR = SNR.astype(float)
    A90 = A90.astype(float)
    A50 = A50.astype(float)
    OutputTable = Table([run, gal, MergerID, RA, Dec, distance, z, theta, ndet, SNR, A90, A50], names=(
        'run', 'Galaxy', 'MergerID', 'RA', 'Dec', 'Distance', 'redshift', 'theta', 'ndet', 'SNR', 'A90', 'A50'))
    # print(OutputTable)
    return OutputTable


def TableImportCTA_TimeNoZenith(ttimeFile):
    run, MergerID, time1, time2, Observatory = np.genfromtxt(ttimeFile, usecols=(0, 1, 2, 3, 4), skip_header=1,
                                                             unpack=True, dtype='str')
    time = []
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([run, MergerID, time, Observatory], names=[
                        'run', 'MergerID', 'Time', 'Observatory'])

    return OutputTable


def TableImportCTA_Petrov(ttimeFile):
    ID, time1, time2 = np.genfromtxt(ttimeFile, usecols=(
        0, 1, 2), skip_header=1, unpack=True, dtype='str')
    time = []
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([ID, time], names=['ID', 'Time'])

    return OutputTable


def TableImportCTA_SetOfTimes(ttimeFile):
    trial, run, MergerID, time1, time2, Observatory = np.genfromtxt(ttimeFile, usecols=(0, 1, 2, 3, 4, 5),
                                                                    skip_header=1, unpack=True, dtype='str')
    time = []
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([run, MergerID, trial, time, Observatory],
                        names=['run', 'MergerID', 'trial', 'Time', 'Observatory'])

    return OutputTable


def TableImportCTA_Time(ttimeFile):
    run, MergerID, time1, time2, MeanAlt, Observatory = np.genfromtxt(ttimeFile, usecols=(0, 1, 2, 3, 4, 5),
                                                                      skip_header=1, unpack=True, dtype='str')
    MeanAlt = MeanAlt.astype(int)
    time = []
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([run, MergerID, time, MeanAlt, Observatory],
                        names=['run', 'MergerID', 'Time', 'Zenith', 'Observatory'])

    return OutputTable


def TableImportCTA_Obs(tobsFile):
    pointing, tI, tF, interval = np.genfromtxt(tobsFile, usecols=(
        0, 1, 2, 3), skip_header=3, unpack=True, dtype='int')
    OutputTable = Table([pointing, tI, tF, interval], names=[
                        'pointingNumber', 'tI', 'tF', 'Interval'])
    return OutputTable


def TableImportCTA_LS(tgalFile):
    eventid, RA, Dec, distance = np.genfromtxt(tgalFile, usecols=(
        0, 3, 4, 8), skip_header=1, unpack=True, dtype='str')
    RA = RA.astype(float)
    Dec = Dec.astype(float)
    distance = distance.astype(float) / 1000  # to Mpc!
    OutputTable = Table([eventid, RA, Dec, distance],
                        names=('MergerID', 'RA', 'Dec', 'Distance'))
    return OutputTable


def IsSourceInside(Pointings, Sources, FOV, nside):
    tt = 0.5 * np.pi - Sources.dec.rad
    tp = Sources.ra.rad
    txyz = hp.ang2pix(nside, tt, tp)
    Npoiting = ''
    Found = False
    try:
        for i in range(0, len(Pointings)):
            # t = 0.5 * np.pi - Pointings[i].dec.rad[0]
            # p = Pointings[i].ra.rad[0]
            t = 0.5 * np.pi - Pointings[i].dec.rad
            p = Pointings[i].ra.rad
            # print('targetCoord1[0].ra.deg, targetCoord1[0].dec.deg',HESS_Sources.ra.deg, HESS_Sources.dec.deg, Pointings[i].ra.deg, Pointings[i].dec.deg)
            xyz = hp.ang2vec(t, p)
            try:
                # print(xyz)
                ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))
            except:
                ipix_disc = hp.query_disc(nside, xyz[0], np.deg2rad(FOV))
                # print('Smething went wrong')
            # print(txyz)
            # print(ipix_disc)
            # print(txyz in ipix_disc)
            if (txyz in ipix_disc):
                print('Found in pointing number', i)
                #Npoiting.append(i)
                Npoiting = Npoiting+str(i)+','
                Found = True
        if Found == False:
            print('Source not covered!')
    except TypeError:
        t = 0.5 * np.pi - Pointings.dec.rad
        p = Pointings.ra.rad
        # print('targetCoord1[0].ra.deg, targetCoord1[0].dec.deg',HESS_Sources.ra.deg, HESS_Sources.dec.deg, Pointings[i].ra.deg, Pointings[i].dec.deg)
        xyz = hp.ang2vec(t, p)
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))
        # print(txyz)
        # print(ipix_disc)
        # print(txyz in ipix_disc)
        if (txyz in ipix_disc):
            Npoiting = '0,'
            Found = True
            print('Found in pointing number 0')
        else:
            print('Source not covered!')
    #Reformat output
    if Found == True: 
        Npoiting = Npoiting[:-1]
    return Found, Npoiting


def ProduceSummaryFileOld(Found, InputList, InputObservationList, allPossiblePoint, foundIn, j, typeSimu, totalProb, datasetDir,
                          outDir, name):

    filepath = datasetDir + '/GammaCatalogV2.0/' + str(InputList['run'][j]) + '_' + str(
        InputList['MergerID'][j].split('r')[-1]) + ".fits"
    fitsfile = fits.open(filepath)
    luminosity = fitsfile[0].header['EISO']
    dirNameFile = outDir + '/SummaryFile'
    if not os.path.exists(dirNameFile):
        os.makedirs(dirNameFile)

    # Obtain the luminosity
    if foundIn == -1:
        outfilename = outDir + '/SummaryFile/' + name + \
            '_SimuS' + typeSimu + str("{:03d}".format(j)) + '.txt'
        f = open(outfilename, 'w')
        f.write(
            'ID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Luminosity' + ' ' + 'TotalObservations' + ' ' + 'TotalPossible' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'TotalProb' + ' ' + 'ObsInfo' + '\n')
        f.write(InputList['ID'][j] + ' ' + str(InputList['Distance'][j]) + ' ' + str(InputList['theta'][j]) + ' ' + str(InputList['A90'][j]) + ' ' + str(
            luminosity) + ' ' + str(len(InputObservationList)) + ' ' + str(allPossiblePoint) + ' ' + str(
            foundIn) + ' ' + str(0) + ' ' + str(totalProb) + ' ' + 'True' + '\n')
    else:
        if type(foundIn) != int:
            foundFirst = foundIn[0]
            foundTimes = len(foundIn)  # Has it been observed several times?
        else:
            foundFirst = foundIn
            foundTimes = 1

        duration = InputObservationList['Duration[s]']

        outfilename = outDir + '/SummaryFile/' + name + \
            '_SimuSF' + typeSimu + str("{:03d}".format(j)) + '.txt'
        # print(outfilename)
        f = open(outfilename, 'w')
        f.write(
            'ID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Luminosity' + ' ' + 'TotalObservations' + ' ' + 'TotalPossible' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'TotalProb' + ' ' + 'ObsInfo' + '\n')
        f.write(str(InputList['ID'][j]) + ' ' + str(
            InputList['Distance'][j]) + ' ' + str(InputList['theta'][j]) + ' ' + str(InputList['A90'][j]) + ' ' + str(
            luminosity) + ' ' + str(len(InputObservationList)) + ' ' + str(allPossiblePoint) + ' ' + str(
            foundFirst) + ' ' + str(foundTimes) + ' ' + str(totalProb) + ' ' + 'True' + '\n')


def FillSummary(outfilename, ID, doneObservations, totalPoswindow, foundFirst, nP, totalProb, ObsInfo):
    f = open(outfilename, 'w')
    f.write(
        'ID' + ' ' + 'TotalObservations' + ' ' + 'TotalPossible' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'TotalProb' + ' ' + 'ObsInfo' + '\n')
    f.write(str(ID) + ' ' + str(doneObservations) + ' ' + str(totalPoswindow) + ' ' + str(
        foundFirst) + ' ' + str(nP) + ' ' + str(totalProb) + ' ' + str(ObsInfo) + '\n')


def ProduceSummaryFile(Source, SuggestedPointings, totalPoswindow, ID, obspar, typeSimu, datasetDir, outDir):

    # Where to save results
    dirNameFile = outDir + '/SummaryFile/'
    print(dirNameFile)
    if not os.path.exists(dirNameFile):
        os.makedirs(dirNameFile)

    dirNameSch = outDir + '/ScheduledObs'
    if not os.path.exists(dirNameSch):
        os.makedirs(dirNameSch)

    dirNameSchCommas = outDir + '/ScheduledObsCommas'
    if not os.path.exists(dirNameSchCommas):
        os.makedirs(dirNameSchCommas)

    # grbFilename = datasetDir +'GRB-GW_TeV_catO5/catO5_'+ str(ID) + '.fits'
    # fitsfile = fits.open(grbFilename)

    # Interesting parameters f
    # luminosity = fitsfile[0].header['EISO']

    # Is the source inside?
    maskClean = (SuggestedPointings['ObsInfo'] == 'True')
    SuggestedPointingsC = SuggestedPointings[maskClean]
    SuggestedPointingsC.remove_column('ObsInfo')

    print('Source coordinates', Source)
    print(SuggestedPointingsC)

    Pointings = SkyCoord(SuggestedPointingsC['RA[deg]'], SuggestedPointingsC['DEC[deg]'], frame='fk5',
                         unit=(u.deg, u.deg))

    totalPGW = float('{:1.4f}'.format(float(sum(SuggestedPointingsC['PGW']))))

    # Check if the source is covered by the scheduled pointings
    Found, nP = IsSourceInside(
        Pointings, Source, obspar.FOV, obspar.reducedNside)
    foundFirst = -1
    #if len(nP) == 0:
    #    nP = 0
    if 'True' in SuggestedPointings['ObsInfo'] and Found == True:
        print('Found in scheduled observation:', nP)
        FoundFirst = nP[0]

        # --- Writting down the results ---
        pointingsFileC = '%s/%s_cov.txt' % (dirNameSch, ID)
        ascii.write(SuggestedPointingsC, pointingsFileC,
                    overwrite=True, fast_writer=False)

        pointingsFileCommas = '%s/%s_cov.txt' % (dirNameSchCommas, ID)
        ascii.write(SuggestedPointingsC, pointingsFileCommas,
                    format='csv', overwrite=True, fast_writer=False)

        outfilename = dirNameFile + str(ID) + '_SimuSF_' + typeSimu + '.txt'
        FillSummary(outfilename, ID, len(SuggestedPointingsC), totalPoswindow,
                    FoundFirst, nP, np.sum(SuggestedPointings['PGW']), str(2))

    if 'True' in SuggestedPointings['ObsInfo'] and Found == False:
        print('Source not covered')

        # print('Plotting the observations')
        # --- Writting down the results ---
        pointingsFileC = '%s/%s_NOTcov.txt' % (dirNameSch, ID)
        ascii.write(SuggestedPointingsC, pointingsFileC,
                    overwrite=True, fast_writer=False)

        pointingsFileCommas = '%s/%s_NOTcov.txt' % (dirNameSchCommas, ID)
        ascii.write(SuggestedPointingsC, pointingsFileCommas,
                    format='csv', overwrite=True, fast_writer=False)

        outfilename = dirNameFile + str(ID) + '_SimuS_' + typeSimu + '.txt'
        FillSummary(outfilename, ID, len(SuggestedPointingsC), totalPoswindow,
                    foundFirst, nP, np.sum(SuggestedPointings['PGW']), str(1))

    if 'True' not in SuggestedPointings['ObsInfo']:
        print('No observations are scheduled, lets write it down')
        # --- Writting down the results ---
        totalPGW = 0
        outfilename = dirNameFile + str(ID) + '_Simu_' + typeSimu + '.txt'
        FillSummary(outfilename, ID, 0, totalPoswindow,
                    foundFirst, nP, totalPGW, str(0))


def ProducePandasSummaryFile(Source, SuggestedPointings, totalPoswindow, ID, obspar, typeSimu, datasetDir, outDir, configID):

    # Where to save results
    dirNameFile = outDir + '/PandasSummaryFile/'
    #print(dirNameFile)
    
    if not os.path.exists(dirNameFile):
        os.makedirs(dirNameFile)
    # What should be in the pandas file is the following: 
    # 
    # Is the source inside?
    maskClean = (SuggestedPointings['ObsInfo'] == 'True')
    SuggestedPointingsC = SuggestedPointings[maskClean]

    Pointings = SkyCoord(SuggestedPointingsC['RA[deg]'], SuggestedPointingsC['DEC[deg]'], frame='fk5',
                         unit=(u.deg, u.deg))

    # Check if the source is covered by the scheduled pointings
    Found, nP = IsSourceInside(Pointings, Source, obspar.FOV, obspar.reducedNside)

    totalPGW = float('{:1.4f}'.format(float(sum(SuggestedPointingsC['PGW']))))
    
    # Rename columns 
    SuggestedPointingsC.rename_column('Observation Time UTC', 'obs_time_utc')
    SuggestedPointingsC.rename_column('DEC[deg]', 'dec')
    SuggestedPointingsC.rename_column('RA[deg]', 'ra')
    SuggestedPointingsC.rename_column('Observatory', 'observatory')
    SuggestedPointingsC.rename_column('ZenIni[deg]', 'zenith_init')
    SuggestedPointingsC.rename_column('ZenEnd[deg]', 'zenith_end')
    SuggestedPointingsC.rename_column('Duration[s]','duration')
    SuggestedPointingsC.rename_column('Delay[s]','delay')

    SuggestedPointingsC.remove_column('ObsInfo')
    SuggestedPointingsC.remove_column('PGW')


    if 'True' not in SuggestedPointings['ObsInfo']:
        print('No observations are scheduled')
        totalPGW = 0
        totalObservations = 0
        pointings = []
    else: 
        totalPGW = np.sum(SuggestedPointings['PGW'])
        totalObservations = len(SuggestedPointingsC['ra'])
        pointings = SuggestedPointingsC.to_pandas().to_dict(orient='records')
        # Convert all values to strings in the list of dictionaries
        # for my_dict in pointings:
        #    for key in my_dict:
        #        my_dict[key] = str(my_dict[key])
        # Convert all values to strings
        # pointings_as_strings = {key: str(value) for key, value in pointings.items()}

        #print(pointings)
    data = {
        'obs_run': int(ID),
        'config_file': str(configID),
        'total_observations': str(totalObservations),
        'total_prob': str(totalPGW),
        'pointings_found': str(nP),
        'found': str(Found),
        'n_found': str(len(nP)), 
        'pointings': list([pointings]), 
    }
       
    obs_run = configID.split('_')[0]
    layout = configID.split('_')[1] 
    has_ebl = False
    if 'EBL' in configID:
        has_ebl = True

    data = {
            'event_id': int(ID),
            'cta_layout': str(layout),
            'ebl': bool(has_ebl),
            'obs_run': str(obs_run),
            'config_file': str(configID),
            'total_observations': int(totalObservations),
            'total_prob': float(totalPGW),
            'pointings_found': list([nP]),
            'found': bool(Found),
            'n_found': int(len(nP)), 
            'pointings': list([pointings]), 
        }



    #print(data)
    # Create the DataFrame
    df = pd.DataFrame(data)

    # Display the DataFrame
    print(df.iloc[0].pointings)
    df.to_parquet(dirNameFile+str(ID)+'_'+configID+'.parquet')
    #import pickle
    #with open(dirNameFile+str(ID)+'_'+configID+'.pkl', "wb") as f:
    #    pickle.dump(df, f)

def ReadSummaryFile(summaryFile):
    print('Note that this function needs to be adapted to the output')
    nP, Found = np.genfromtxt(summaryFile, usecols=(
        8, 9), skip_header=1, unpack=True, dtype='str')
    return nP, Found


class NextWindowTools:

    @classmethod
    def CheckWindowCreateArray(cls, time, obsSite, WindowDurations):
        FullWindow = datetime.timedelta(
            seconds=np.float64(WindowDurations[-1]))
        # print('time',time)
        # print('FullWindow',FullWindow)
        # print(len(WindowDurations))
        if (Tools.IsDarkness(time, obsSite) is True) and (Tools.IsDarkness(time + FullWindow, obsSite) is True):
            LastItem = len(WindowDurations)
        else:
            print('Window is smaller')
            for i in range(0, len(WindowDurations)):
                if Tools.IsDarkness(time + datetime.timedelta(minutes=np.float64(WindowDurations[-i])), obsSite):
                    # print('Found!')
                    LastItem = len(WindowDurations) - i
        cumsumWindow = np.cumsum(WindowDurations)
        # print(cumsumWindow)
        # print(datetime.timedelta(seconds=np.float64(cumsumWindow[j])) for j in range(LastItem))
        arr = np.array([time + datetime.timedelta(seconds=np.float64(cumsumWindow[j]))
                       for j in range(LastItem)])
        return arr

    @classmethod
    def NextObservationWindow(cls, time, obspar):
        if (Tools.NextSunset(time, obspar).hour >= time.hour >= Tools.PreviousSunrise(time, obspar).hour and time.day == Tools.NextSunset(time, obspar).day):
            time = Tools.NextSunset(time, obspar)
            time = Tools.TrustingDarknessSun(time, obspar)
        if (Tools.IsDarkness(time, obspar) is True):
            return time
        elif ((Tools.IsDarkness(time, obspar) is False) and (
                Tools.IsDarkness(Tools.NextMoonset(time, obspar), obspar) is True)):
            time = Tools.NextMoonset(time, obspar)
            return time
        else:
            print('No window is found')
            return False

    @classmethod
    def NextObservationWindowGrey(cls, time, obspar):
        if (Tools.NextSunset(time, obspar).hour >= time.hour >= Tools.PreviousSunrise(time, obspar).hour and time.day == Tools.NextSunset(time, obspar).day):
            time = Tools.NextSunset(time, obspar)
            # print('Sunset', time)
            time = Tools.TrustingGreynessSun(time, obspar)
            # print('Trusted', time)
        if (Tools.IsGreyness(time, obspar) is True):
            return time
        elif ((Tools.IsGreyness(time, obspar) is False)):
            time = Tools.TrustingGreynessSun(time, obspar)
            # time=Tools.NextMoonset(time, obsSite)
            return time
        else:
            print('No window is found')
            return False

    @classmethod
    def EndObservationWindow(cls, time, obsSite):
        time = Tools.NextSunrise(time, obsSite)
        # Check if the night ends before due to the moon.
        if (Tools.IsDarkness(time, obsSite) is False):
            time = Tools.PreviousMoonset(time, obsSite)
        return time

    @classmethod
    def AddRunDuration_to_StartingTime(cls, obspar):  
        previousTime = np.genfromtxt(obspar.pointingsFile, usecols=(0),skip_header=1, unpack=True, dtype='str')   
        obspar.obsTime = datetime.datetime.strptime(str(previousTime), '%Y-%m-%dT%H:%M:%S')+ datetime.timedelta(minutes=np.float64(obspar.duration))

def ZenithAngleCut_TwoTimes(prob, nside, time, time1, minProbcut, maxZenith, observatory):
    '''
    Mask in the pixels with zenith angle larger than maxZenith
    '''

    # Initial time
    frame = co.AltAz(obstime=time, location=observatory)
    pprob = prob

    mzenith = hp.ma(pprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
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

    # End time
    frame = co.AltAz(obstime=time1, location=observatory)
    ppprob = pprob

    mzenith = hp.ma(ppprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    altaz_map = targetCoord_map.transform_to(frame)
    maskzenith[altaz_map.alt.value < 90 - maxZenith] = 1
    mzenith.mask = maskzenith

    # plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_%g.png")

    xprob = ma.masked_array(yprob, mzenith.mask)

    if np.sum(yprob) < minProbcut:
        ObsBool = False
    else:
        ObsBool = True

    return ObsBool, yprob


def ComputeProbability2D_SelectClusters(prob, highres, radecs, conf, time, DelayObs, interObsSlew, obspar, ID, ipixlist, ipixlistHR, counter, datasetDir, outDir, useGreytime, plot):
    '''
    Compute probability in 2D by taking the highest value pixel
    '''

    reducedNside = obspar.reducedNside
    HRnside = obspar.HRnside
    minProbcut = obspar.minProbcut
    maxZenith = obspar.maxZenith
    radius = obspar.FOV

    frame = co.AltAz(obstime=time, location=obspar.location)
    thisaltaz = radecs.transform_to(frame)
    pix_alt1 = thisaltaz.alt.value

    grbFilename = datasetDir + 'GRB-GW_TeV_catO5/catO5_' + ID + '.fits'

    if useGreytime:
        moonaltazs = get_moon(Time(time, scale='utc')).transform_to(
            AltAz(obstime=Time(time), location=obspar.location))
        # Zenith and Moon angular distance mask
        pix_ra = radecs.ra.value[
            (thisaltaz.alt.value > 90 - maxZenith) & (thisaltaz.separation(moonaltazs) > 30 * u.deg)]
        pix_dec = radecs.dec.value[
            (thisaltaz.alt.value > 90 - maxZenith) & (thisaltaz.separation(moonaltazs) > 30 * u.deg)]

    else:
        # Zenith angle mask
        pix_ra = radecs.ra.value[thisaltaz.alt.value > (90 - maxZenith)]
        pix_dec = radecs.dec.value[thisaltaz.alt.value > (90 - maxZenith)]
        zenith_ini = 90 - pix_alt1[thisaltaz.alt.value > (90 - maxZenith)]

    phipix = np.deg2rad(pix_ra)
    thetapix = 0.5 * np.pi - np.deg2rad(pix_dec)

    ipix = hp.ang2pix(reducedNside, thetapix, phipix)

    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)
    exposure = np.empty(len(pix_ra), dtype=object)
    zenith_end = np.empty(len(pix_ra), dtype=object)
    # print(len(phipix),len(thetapix),len(ipix),len(pix_ra),len(pix_dec), len(dp_Pix_Fov), len(pix_alt), len(zenith_end), len(exposure))

    cat_pix = Table([ipix, pix_ra, pix_dec, dp_Pix_Fov, zenith_ini, zenith_end, exposure],
                    names=('PIX', 'PIXRA', 'PIXDEC', 'PIXFOVPROB', 'ZENITH_INI', 'ZENITH_END', 'EXPOSURE'))

    dp_dV_FOV = []

    xyzpix = hp.ang2vec(thetapix, phipix)

    for i in range(0, len(cat_pix)):
        ipix_discfull = hp.query_disc(HRnside, xyzpix[i], np.deg2rad(radius))
        maskComputeProb = np.isin(ipix_discfull, ipixlistHR, invert=True)
        dp_dV_FOV.append(highres[ipix_discfull[maskComputeProb]].sum())

    cat_pix['PIXFOVPROB'] = dp_dV_FOV

    # Mask already observed pixels

    mask = np.isin(cat_pix['PIX'], ipixlist, invert=True)

    if all(np.isin(cat_pix['PIX'], ipixlist, invert=False)):
        maskcat_pix = cat_pix
    else:
        maskcat_pix = cat_pix[mask]

    # Sort table
    sortcat = maskcat_pix[np.flipud(np.argsort(maskcat_pix['PIXFOVPROB']))]
    ObsCase = 'SourceOutFoV'  # Default case, the source is not in CTA FoV

    # Fill a the column EXPOSURE column. Corresponds to the time that one needs to observe to get 5sigma for the highest of the list
    # Three cases depending on the IRFs that should be used (60,40,20)
    grbSensPath = '/grbsens_output_v3_Sep_2022/'+ conf +'_configuration_EBL/grbsens-5.0sigma_t1s-t16384s_irf-'
    print(datasetDir + grbSensPath +obspar.name+ "_z60_0.5h.txt")
    if (np.any(sortcat['ZENITH_INI'] > 55)):
        # ObsCase, texp60 = ObtainSingleObservingTimes(TotalExposure, DelayObs, interObsSlew, ID, obspar,datasetDir, zenith=60)
        # if observatory.name == "North":
        # "sensitivity-5sigma_irf-North_z20_0.5.txt"
        grbSensFile = datasetDir + grbSensPath +obspar.name+ "_z60_0.5h.txt"
        grb_result = GetExposureForDetection(
            grbSensFile, grbFilename, DelayObs)
        print(grb_result)
        if (grb_result['obs_time'] == -1):
            ObsCase = 'TimeNotEnough'
            sortcat['EXPOSURE'][sortcat['ZENITH_INI'] > 55] = False
        else:
            texp60 = grb_result['obs_time']
            print("ObsCase60", ObsCase, 'time =', texp60)
            # Cat60 = sortcat[sortcat['ZENITH_INI'] >55]
            # print("ObsCase60", ObsCase)
            sortcat['EXPOSURE'][sortcat['ZENITH_INI'] > 55] = texp60
            frame = co.AltAz(
                obstime=time + datetime.timedelta(seconds=texp60), location=obspar.location)
            catCoord60 = co.SkyCoord(sortcat['PIXRA'][sortcat['ZENITH_INI'] > 55],
                                     sortcat['PIXDEC'][sortcat['ZENITH_INI'] > 55], frame='fk5', unit=(u.deg, u.deg))
            thisaltaz60 = catCoord60.transform_to(frame)
            pix_zen60 = 90 - thisaltaz60.alt.value
            sortcat['ZENITH_END'][sortcat['ZENITH_INI'] > 55] = pix_zen60

    mask1 = sortcat['ZENITH_INI'] >= 30
    mask2 = sortcat['ZENITH_INI'] <= 55

    if (sortcat['ZENITH_INI'][mask1 & mask2].any()):
        # ObsCase, texp40 = ObtainSingleObservingTimes(TotalExposure, DelayObs, interObsSlew, ID, obspar,datasetDir, zenith=40)

        # "sensitivity-5sigma_irf-North_z20_0.5.txt"
        grbSensFile = datasetDir + grbSensPath + obspar.name + "_z40_0.5h.txt"
        grb_result = GetExposureForDetection(
            grbSensFile, grbFilename, DelayObs)
        print(grb_result)
        if (grb_result['obs_time'] == -1):
            ObsCase = 'TimeNotEnough'
            sortcat['EXPOSURE'][(30 <= sortcat['ZENITH_INI']) & (
                sortcat['ZENITH_INI'] <= 55)] = False
        else:
            texp40 = grb_result['obs_time']
            print("ObsCase40", ObsCase, 'time=', texp40)
            sortcat['EXPOSURE'][(30 <= sortcat['ZENITH_INI']) & (
                sortcat['ZENITH_INI'] <= 55)] = texp40
            # Cat40 = sortcat[(30 < sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] < 55)]
            # print(Cat40)
            frame = co.AltAz(
                obstime=time + datetime.timedelta(seconds=texp40), location=obspar.location)
            catCoord40 = co.SkyCoord(sortcat['PIXRA'][(30 <= sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] <= 55)],
                                     sortcat['PIXDEC'][(30 <= sortcat['ZENITH_INI']) & (
                                         sortcat['ZENITH_INI'] <= 55)],
                                     frame='fk5', unit=(u.deg, u.deg))
            # print("radecs",radecs)
            thisaltaz40 = catCoord40.transform_to(frame)
            pix_zen40 = 90 - thisaltaz40.alt.value
            sortcat["ZENITH_END"][(30 <= sortcat['ZENITH_INI']) & (
                sortcat['ZENITH_INI'] <= 55)] = pix_zen40

    if (np.any(sortcat['ZENITH_INI'] < 30)):
        # ObsCase, texp20 = ObtainSingleObservingTimes(TotalExposure, DelayObs, interObsSlew, ID, obspar, datasetDir, zenith=20)

        # "sensitivity-5sigma_irf-North_z20_0.5.txt"
        grbSensFile = datasetDir + grbSensPath + obspar.name + "_z20_0.5h.txt"
        grb_result = GetExposureForDetection(
            grbSensFile, grbFilename, DelayObs)
        if (grb_result['obs_time'] == -1):
            ObsCase = 'TimeNotEnough'
            sortcat['EXPOSURE'][sortcat['ZENITH_INI'] < 30] = False
        else:
            texp20 = grb_result['obs_time']
            print("ObsCase20", ObsCase, 'time=', texp20)
            sortcat['EXPOSURE'][sortcat['ZENITH_INI'] < 30] = texp20
            # Cat20 = sortcat[sortcat['ZENITH_INI'] < 30]
            # print(Cat30)
            # print("time + datetime.timedelta(seconds=texp20)",time + datetime.timedelta(seconds=texp20))
            # print("observatory.location",observatory.location)
            frame = co.AltAz(
                obstime=time + datetime.timedelta(seconds=texp20), location=obspar.location)
            catCoord20 = co.SkyCoord(sortcat['PIXRA'][sortcat['ZENITH_INI'] < 30],
                                     sortcat['PIXDEC'][sortcat['ZENITH_INI'] < 30], frame='fk5', unit=(u.deg, u.deg))
            thisaltaz20 = catCoord20.transform_to(frame)
            pix_zen20 = 90 - thisaltaz20.alt.value
            sortcat["ZENITH_END"][sortcat['ZENITH_INI'] < 30] = pix_zen20

    # Check how many false values there are in the array
    # false_count = sum(not i for i in sortcat['EXPOSURE'])
    
    # print("Number of False values:", false_count,'vs. the number of entries being', len(sortcat['EXPOSURE']))
    # Mask the catalog from the entries that are actually not feaseable from exposure value
    if False in sortcat['EXPOSURE']:
        maskExposure = (sortcat['EXPOSURE'] != False)
        print(maskExposure)
        sortcat = sortcat[maskExposure]

    # Mask if it is not visible at the end of the window
    mask = np.isin(sortcat['ZENITH_END'], 66, invert=True)
    maskcat_zen = sortcat[mask]
    if (len(sortcat[mask]) == 0):
        P_GW = 0
        targetCoord = False
        ObsExp = False
        ZenIni = False
        ZenEnd = False
        print('Obscase', ObsCase)
    else:
        sortcat = maskcat_zen[np.flipud(np.argsort(maskcat_zen['PIXFOVPROB']))]
        targetCoord = co.SkyCoord(
            sortcat['PIXRA'][:1][0], sortcat['PIXDEC'][:1][0], frame='fk5', unit=(u.deg, u.deg))
        # print('targetCoord',targetCoord)

        P_GW = sortcat['PIXFOVPROB'][:1][0]
        ObsExp = sortcat['EXPOSURE'][:1][0]
        ZenIni = sortcat['ZENITH_INI'][:1][0]
        ZenEnd = sortcat['ZENITH_END'][:1][0]

        # Include to the list of pixels already observed

        if (P_GW >= minProbcut):
            # Chose highest
            phip = float(np.deg2rad(targetCoord.ra.deg))
            thetap = float(0.5 * np.pi - np.deg2rad(targetCoord.dec.deg))
            xyz = hp.ang2vec(thetap, phip)

            # ipix_discComputeProb = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
            # maskComputeProb=[np.isin(ipix_discComputeProb, ipixlist,invert=True)]
            # print('maskedP_GW',highres[ipix_discComputeProb[maskComputeProb]].sum())
            ipixlistHR.extend(hp.query_disc(HRnside, xyz, np.deg2rad(radius)))
            ipix_disc = hp.query_disc(reducedNside, xyz, np.deg2rad(radius))
            ipixlist.extend(ipix_disc)

            ######################################

            # PLOT THE RESULTS
            if (plot):
                # path = outDir + '/EvolutionPlot'
                path = '%s/Pointing_Plotting/%s/EvolutionPlot' % (outDir, ID)
                if not os.path.exists(path):
                    os.makedirs(path)
                # nside = 1024

                hp.mollview(prob, title="With FoV circle")

                # hp.gnomview(prob, xsize=500, ysize=500, rot=[targetCoord.ra.deg, targetCoord.dec.deg], reso=8.0)
                hp.graticule()
                # print('This skymap has nside equals to',hp.npix2nside(len(highres)))
                # plt.savefig('%s/Pointing-prob_%g.png' % (path, counter))

                ipix_discplot = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
                tt, pp = hp.pix2ang(HRnside, ipix_discplot)
                ra2 = np.rad2deg(pp)
                dec2 = np.rad2deg(0.5 * np.pi - tt)
                skycoord = co.SkyCoord(
                    ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
                # hp.visufunc.projplot(skycoord.ra, skycoord.dec, 'y.', lonlat=True, coord="C")
                # plt.show()
                # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

                hp.visufunc.projplot(
                    sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], 'r.', lonlat=True, coord="C")
                MaxCoord = SkyCoord(
                    sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))
                separations = skycoord.separation(MaxCoord)
                tempmask = separations < (radius + 0.01 * radius) * u.deg
                tempmask2 = separations > (radius - 0.01 * radius) * u.deg
                hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'g.',
                                     lonlat=True,
                                     coord="C", linewidth=0.1)
                plt.savefig('%s/Pointing-prob-FoV_%g.png' % (path, counter))
                altcoord = np.empty(100)
                azcoord = np.random.rand(100) * 360
                plt.savefig('%s/Pointing-prob-FoV_%g.png' % (path, counter))
            # for i in range(0,1):
            #    altcoord.fill(90-(maxZenith-5*i))
            #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
            #    RandomCoord_radec = RandomCoord.transform_to('fk5')
            #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
            # plt.show()
            # plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))
        altcoord = np.empty(100)
        azcoord = np.random.rand(100) * 360
        # for i in range(0,1):
        #    altcoord.fill(90-(maxZenith-5*i))
        #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
        #    RandomCoord_radec = RandomCoord.transform_to('fk5')
        #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # plt.show()
        # plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))

    return P_GW, targetCoord, ObsExp, ZenIni, ZenEnd, ObsCase, ipixlist, ipixlistHR


def GetExposureForDetection(grbSensFile, grbFilename, DelayObs):
    sens = Sensitivity(grbSensFile, min_energy=0.3, max_energy=10000,)
    grb = GRB(grbFilename)
    return grb.observe(sensitivity=sens, start_time=DelayObs, target_precision=0.1)