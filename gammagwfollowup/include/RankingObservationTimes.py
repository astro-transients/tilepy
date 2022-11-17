import datetime
import os
import astropy.coordinates as co
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun,get_moon
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.data import download_file
import copy
from .PointingTools import (LoadHealpixMap, Tools, CorrelateGalaxies_LVC,
                                  CorrelateGalaxies_LVC_SteMass, ObservationParameters)

from six.moves import configparser
import six
if six.PY2:
  ConfigParser = configparser.SafeConfigParser
else:
  ConfigParser = configparser.ConfigParser


iers_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)

'''
    gloably defined darkenss criteria:
    max sun and moon altitude in degrees.
    '''
'''#################################################
#   Global parameters about darkness criteria   #
#################################################
# max sun and moon altitude in degrees.
cfg = "./configs/FollowupParameters.ini"
parser = ConfigParser()
parser.read(cfg)
parser.sections()
section = 'visibility'

try:
    gSunDown = int(parser.get(section, 'gSunDown'))
    HorizonSun = parser.get(section, 'HorizonSun')
    gMoonDown = float(parser.get(section, 'gMoonDown'))
    HorizonMoon = (parser.get(section, 'HorizonMoon'))
    gMoonGrey = int(parser.get(section, 'gMoonGrey')) # Altitude in degrees
    gMoonPhase = int(parser.get(section, 'gMoonPhase'))  # Phase in %
    MoonSourceSeparation = int(parser.get(section, 'MoonSourceSeparation')) # Separation in degrees
    MaxMoonSourceSeparation = int(parser.get(section, 'MaxMoonSourceSeparation'))  # Max separation in degrees

except Exception as x:
    print(x)
'''


def load_healpix_map(filename):
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
    PrintFileName = "Loading LVC HEALPix map from file: " + filename
    #print(PrintFileName)
    fitsfile = fits.open(filename)

    tevent_id = "Non specified"
    tdetectors=""
    tdistmean=0
    tdisterr=0
    tdistmu=[]
    tdistsigma=[]
    tdistnorm=[]

    if 'OBJECT' in fitsfile[1].header:
        tevent_id = fitsfile[1].header['OBJECT']
    else:
        tevent_id = "Non specified"

    if 'INSTRUME' in fitsfile[1].header:
        tdetectors = fitsfile[1].header['INSTRUME']
    else:
        tdetectors = "Non specified"

    if(fitsfile[1].header['TFIELDS']==4):
        tprob, tdistmu, tdistsigma, tdistnorm = hp.read_map(filename, field=range(4))
        tdistmean=fitsfile[1].header['DISTMEAN']
        tdisterr=fitsfile[1].header['DISTSTD']
        print('Event has triggered ',tdetectors,' => distance = ',tdistmean,' +- ',tdisterr, ' Mpc')
    else:
        tprob = hp.read_map(filename, field=range(1))
    #raise

    fitsfile.close()

    return tprob, tdistmu, tdistsigma, tdistnorm, tdetectors, tevent_id, tdistmean, tdisterr


def load_pointingFile(tpointingFile):
    #Read PointingsFile

    print("Loading pointings from " + tpointingFile)
    time1,time2, ra, dec, = np.genfromtxt(tpointingFile, usecols=(0, 1, 2,3), dtype="str", skip_header=1,
                                             delimiter=' ',
                                             unpack=True)  # ra, dec in degrees
    time1 = np.atleast_1d(time1)
    time2 = np.atleast_1d(time2)

    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    time=[]

    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i].split(':')[0]+ ':'+time2[i].split(':')[1]).split('"')[1])

    ra = ra.astype(np.float)
    dec = dec.astype(np.float)

    l = list(range(len(ra)))
    Pointings = Table([l,time,ra,dec],names=('Pointing','Time','RAJ2000', 'DEJ2000'))

    return Pointings

def VisibilityWindow(ObservationTime,Pointing,AltitudeCut,nights,UseGreytime,dirName,max_zenith,observatory, gMoonGrey, gMoonDown, gMoonPhase, gSunDown, MoonSourceSeparation, MaxMoonSourceSeparation):
    source = SkyCoord(Pointing['RAJ2000'],Pointing['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    WINDOW=[]
    ZENITH=[]
    SZENITH=[]

    try:
        auxtime = datetime.datetime.strptime(Pointing['Time'][0], '%Y-%m-%d %H:%M:%S.%f')
    except ValueError:
        try:
            auxtime = datetime.datetime.strptime(Pointing['Time'][0], '%Y-%m-%d %H:%M:%S')
        except ValueError:
            auxtime = datetime.datetime.strptime(Pointing['Time'][0], '%Y-%m-%d %H:%M')

    #frame = co.AltAz(obstime=auxtime, location=observatory)
    timeInitial=auxtime-datetime.timedelta(minutes=30)
    for i in range(0,len(source)):
        NonValidwindow,Stepzenith = GetVisibility(Pointing['Time'],source[i],max_zenith, observatory)
        window, zenith = GetObservationPeriod(timeInitial, source[i],AltitudeCut,observatory,nights,i,dirName,UseGreytime,False, gMoonGrey, gMoonDown, gMoonPhase, gSunDown, MoonSourceSeparation, MaxMoonSourceSeparation)
        WINDOW.append(window)
        ZENITH.append(zenith)
        SZENITH.append(Stepzenith)
        window, zenith = GetObservationPeriod(ObservationTime, source[i],AltitudeCut,observatory,nights,i,dirName,UseGreytime,True, gMoonGrey, gMoonDown, gMoonPhase, gSunDown, MoonSourceSeparation, MaxMoonSourceSeparation)

    Pointing['Observation window'] = WINDOW
    Pointing['Array of zenith angles']=ZENITH
    Pointing['Zenith angles in steps'] = SZENITH

    return Pointing


def GetObservationPeriod(inputtime0,msource,AltitudeCut,observatory,nights,plotnumber,dirName,UseGreytime,doplot, gMoonGrey, gMoonDown, gMoonPhase, gSunDown, MoonSourceSeparation, MaxMoonSourceSeparation):


    inputtime = Time(inputtime0)
    initialframe = AltAz(obstime=inputtime,location=observatory.Location)

    ##############################################################################
    suninitial= get_sun(inputtime).transform_to(initialframe)

    if(suninitial.alt< -18.*u.deg):
        hoursinDay=12
    else:
        hoursinDay =24
    delta_day = np.linspace(0, hoursinDay+24*(nights-1), 1000*nights)*u.hour
    interval=(hoursinDay+24*(nights-1))/(1000.*nights)

    x=np.arange(int(hoursinDay/interval), dtype=int)
    firstN=np.full_like(x,1)
    ratio2=24./interval
    otherN=[]
    for i in range(2,nights+1):
        otherN.extend(np.full_like(np.arange(int(ratio2)),i))
    NightsCounter = []
    NightsCounter.extend(firstN)
    NightsCounter.extend(otherN)

    while len(NightsCounter)!=len(delta_day):
        NightsCounter.extend([nights])

    times = inputtime + delta_day
    frame = AltAz(obstime=times,location=observatory.Location)

    ##############################################################################
    #SUN
    sunaltazs = get_sun(times).transform_to(frame)

    #MOON
    moon = get_moon(times)
    moonaltazs = moon.transform_to(frame)
    msourcealtazs = msource.transform_to(frame)

    #Add Moon phase
    moonPhase = np.full(len(msourcealtazs),Tools.MoonPhase(inputtime0,observatory))


    MoonDistance=msourcealtazs.separation(moonaltazs)
    ##############################################################################
    if UseGreytime:
        Altitudes = Table([times,msourcealtazs.alt, sunaltazs.alt, moonaltazs.alt,moonPhase,MoonDistance,NightsCounter],
                           names=['Time UTC', 'Alt Source', 'Alt Sun', 'AltMoon','moonPhase','MoonDistance','NightsCounter'])
        #selectedTimes=Altitudes['Time UTC']
        selection=(Altitudes['Alt Sun'] < -18.) & (Altitudes['Alt Source']> AltitudeCut) & (Altitudes['AltMoon'] < -0.5)
        DTaltitudes=Altitudes[selection]
        newtimes=[]
        newtimes.extend(DTaltitudes['Time UTC'].mjd)
        selectionGreyness = (Altitudes['AltMoon'] < gMoonGrey)&(Altitudes['AltMoon'] > gMoonDown)&(Altitudes['moonPhase'] <gMoonPhase)&(Altitudes['Alt Sun'] < gSunDown) & (Altitudes['MoonDistance']>MoonSourceSeparation)&(Altitudes['MoonDistance']<MaxMoonSourceSeparation)&(Altitudes['Alt Source']> AltitudeCut)
        GTaltitudes = Altitudes[selectionGreyness]
        newtimes.extend(GTaltitudes['Time UTC'].mjd)
        newtimes = sorted(newtimes)
        ScheduledTimes = Time(newtimes,format='mjd').iso
    else:
        Altitudes = Table([times,msourcealtazs.alt, sunaltazs.alt, moonaltazs.alt,NightsCounter],
                           names=['Time UTC', 'Alt Source', 'Alt Sun', 'AltMoon','NightsCounter'])
        #print(Altitudes)
        Times=Altitudes['Time UTC']
        selection=(Altitudes['Alt Sun'] < -18.) & (Altitudes['Alt Source']> AltitudeCut) & (Altitudes['AltMoon'] < -0.5)
        ScheduledTimes=Times[selection]

    '''if doplot:
        plotDir = '%s/TransitPlots' % dirName
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        #print('delta_day',delta_day)
        #print('delta_day.to(hr)',delta_day.to('hr').value)
        #print('sunaltazs.alt.value',sunaltazs.alt)
        #print((sunaltazs.alt < -18 * u.deg))


        plt.figure(figsize=(20, 12))
        plt.plot(delta_day, sunaltazs.alt, color='r', label='Sun')
        plt.plot(delta_day, moonaltazs.alt, color=[0.75] * 3, ls='--', label='Moon')
        plt.plot(delta_day, msourcealtazs.alt,color = 'b', label='Point source')
        # ToDo: Debug these lines
        #plt.scatter(delta_day, msourcealtazs.alt, c=msourcealtazs.az, label='Point source', lw=0, s=8, cmap='viridis')
        plt.fill_between(delta_day.to('hr').value, 0, 90, sunaltazs.alt < -0 * u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_day.to('hr').value, 0, 90, (sunaltazs.alt < -18 * u.deg) & (moonaltazs.alt < -0.5 * u.deg), color='k', zorder=0)
        #plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        plt.ylim(0, 90)
        #plt.xlabel('Hours after injections')
        plt.ylabel('Altitude [deg]')
        plt.grid()
        plt.savefig('%s/Source%g.png' % (plotDir, plotnumber))'''

    return (str(ScheduledTimes[0]).split('.')[0]+'-->'+str(ScheduledTimes[-1]).split('.')[0]),msourcealtazs.alt

def GetVisibility(time,radecs,max_zenith, observatory):


    visibility = []
    altitude = []

    for i in range(0,len(time)):
        try:
            auxtime = datetime.datetime.strptime(time[i], '%Y-%m-%d %H:%M:%S.%f')
        except ValueError:
            try:
                auxtime = datetime.datetime.strptime(time[i], '%Y-%m-%d %H:%M:%S')
            except ValueError:
                auxtime = datetime.datetime.strptime(time[i], '%Y-%m-%d %H:%M')
        frame = co.AltAz(obstime=auxtime, location=observatory.Location)
        thisaltaz = radecs.transform_to(frame)
        visible = thisaltaz.alt.value > (90-max_zenith)

        if(visible):
            visibility.append(auxtime)
            altitude.append(thisaltaz.alt.value)
        else:
        #    visibility.append(auxtime)
             altitude.append(thisaltaz.alt.value)
    lasttime=auxtime+datetime.timedelta(minutes=30)
    frame = co.AltAz(obstime=lasttime, location=observatory.Location)
    thisaltaz = radecs.transform_to(frame)
    visible = thisaltaz.alt.value > (90 - max_zenith)

    if (visible):
        visibility.append(auxtime)
        altitude.append(thisaltaz.alt.value)
    else:
        #    visibility.append(auxtime)
        altitude.append(thisaltaz.alt.value)

    window = visibility[0].strftime('%H:%M:%S') +'-'+ visibility[-1].strftime('%H:%M:%S')
    return window,altitude

def ProbabilitiesinPointings3D(cat,galPointing,FOV, totaldPdV,prob,nside):

    ra= galPointing['RAJ2000']
    dec=galPointing['DEJ2000']
    PGW= []
    PGAL=[]

    #bucle
    for i in range(0,len(ra)):
        pgwcircle,pgalcircle =PGGPGalinFOV(cat,ra[i], dec[i],prob,totaldPdV,FOV,nside)
        PGW.append(np.float('{:1.4f}'.format(pgwcircle)))
        PGAL.append(np.float('{:1.4f}'.format(pgalcircle)))

    galPointing['Pgw'] = PGW
    galPointing['Pgal'] = PGAL

    return galPointing


def PGGPGalinFOV(cat,ra,dec,prob,totaldPdV,FOV,nside):

    targetCoordcat = co.SkyCoord(cat['RAJ2000'], cat['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    targetCoordpointing = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    dp_dV = cat['dp_dV']

    # Array of indices of pixels inside circle of FoV

    radius = FOV
    t = 0.5 * np.pi - targetCoordpointing.dec.rad
    p = targetCoordpointing.ra.rad

    #print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)
    xyz = hp.ang2vec(t, p)


    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    P_GW = prob[ipix_disc].sum()
    Pgal_inFoV = dp_dV[targetCoordcat.separation(targetCoordpointing).deg <= radius].sum() / totaldPdV

    return P_GW,Pgal_inFoV

def ProbabilitiesinPointings2D(Pointing,FOV,prob,nside):

    ra = Pointing['RAJ2000']
    dec = Pointing['DEJ2000']
    PGW= []
    PGAL=[]
    for i in range(0,len(ra)):
        pgwcircle = PGinFOV(ra[i], dec[i],prob,FOV,nside)
        PGW.append(np.float('{:1.4f}'.format(pgwcircle)))
        PGAL.append(np.float('{:1.4f}'.format(0)))

    Pointing['Pgw'] = PGW
    Pointing['Pgal'] = PGAL

    return Pointing
def PGinFOV(ra,dec,prob,radius,nside):

    targetCoordpointing = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))

    # Array of indices of pixels inside circle of FoV

    t = 0.5 * np.pi - targetCoordpointing.dec.rad
    p = targetCoordpointing.ra.rad

    #print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)
    xyz = hp.ang2vec(t, p)

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    P_GW = prob[ipix_disc].sum()

    return P_GW



def Sortingby(galPointing,targetType, name, exposure):

    gggalPointing = galPointing[np.flipud(np.argsort(galPointing['Pgal']))]
    prioritygal = list(range(len(galPointing['Pgal'])))
    ra=gggalPointing['RAJ2000']
    dec=gggalPointing['DEJ2000']
    coord=SkyCoord(ra,dec,unit='deg')
    #print(coord.to_string('hmsdms'))
    gggalPointing['RA(HH:MM:SS) Dec (DD:MM:SS)'] = coord.to_string('hmsdms')
    gggalPointing['PriorityGal'] = prioritygal
    gggalPointing.remove_column('Array of zenith angles')
    gggalPointing.remove_column('Zenith angles in steps')

    #Prepare filename which is going to be complete
    gwgalPointing = gggalPointing[np.flipud(np.argsort(gggalPointing['Pgw']))]
    prioritygw = list(range(len(galPointing['Pgw'])))
    gwgalPointing['PriorityGW'] = prioritygw

    #gwgalPointing.remove_column('Array of zenith angles')
    #gwgalPointing.remove_column('Zenith angles in steps')
    #print(gwgalPointing)
    outfilename='%s/RankingObservationTimes_Complete.txt' % name
    ascii.write(gwgalPointing[np.argsort(gwgalPointing['Pointing'])], outfilename, overwrite=True)

    gwgalPointing.remove_column('Pgal')
    gwgalPointing.remove_column('Pgw')
    gwgalPointing.remove_column('PriorityGW')
    gwgalPointing.rename_column('PriorityGal','Priority')
    gwgalPointing.remove_column('RA(HH:MM:SS) Dec (DD:MM:SS)')
    outfilename='%s/RankingObservationTimes_forShifters.txt' % name
    ascii.write(gwgalPointing[np.argsort(gwgalPointing['Pointing'])], outfilename, overwrite=True)

    gwgalPointing.remove_column('Observation window')
    gwgalPointing.remove_column('Priority')


    target = [(targetType+ '_' + name.split('/')[2] +'_{0}').format(element) for element in gwgalPointing['Pointing']]
    gwgalPointing['Target'] = target
    gwgalPointing.rename_column('Pointing','Id')
    gwgalPointing['Duration'] = exposure

    gwgalPointing_TH = gwgalPointing['Target', 'Id', 'RAJ2000','DEJ2000','Time','Duration']
    #new_order = ['Target', 'Id', 'RAJ2000','DEJ2000']  # List or tuple
    #gwgalPointing_TH = gwgalPointing[new_order]
    outfilename='%s/RankingObservationTimes_forAlerter.txt' % name
    ascii.write(gwgalPointing_TH[np.argsort(gwgalPointing_TH['Id'])], outfilename, overwrite=True)


    #gggalPointing.remove_column('DEJ2000')
    #gggalPointing.remove_column('RAJ2000')

    #gggalPointing['RA(HH:MM:SS) Dec (DD:MM:SS)'] = coord.to_string('hmsdms')
    #outfilename='%s/RankingObservationTimes_Complete.txt' % name
    #ascii.write(gggalPointing[np.argsort(gggalPointing['Pointing'])], outfilename,overwrite=True)



def EvolutionPlot(galPointing,tname, ObsArray):
    cm = plt.get_cmap('gist_rainbow')

    fig = plt.figure(figsize=(18, 10))
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.8])
    ra=galPointing['RAJ2000']
    dec=galPointing['DEJ2000']
    pgw=galPointing['Pgw']
    pgal=galPointing['Pgal']
    time = galPointing['Time']
    NUM_COLORS = len(time)
    hour=[]
    for j in range(0,len(time)):
        selecttime=time[j].split(" ")
        hour.append(selecttime[1].split('.')[0])
        #print(time[j])
    try:
        lasttime = datetime.datetime.strptime(time[len(time) - 1], '%Y-%m-%d %H:%M') + datetime.timedelta(minutes=30)
    except ValueError:
        lasttime = datetime.datetime.strptime(time[len(time) - 1], '%Y-%m-%d %H:%M') + datetime.timedelta(minutes=30)

    hour.append(lasttime.strftime("%H:%M"))
    GWordered=galPointing[np.flipud(np.argsort(galPointing['Pgw']))]

    ax.set_prop_cycle(plt.cycler('color', plt.cm.Accent(np.linspace(0, 1, NUM_COLORS))))
    for i in range(0,len(ra)):
        #x = np.arange(0, len(ra), 1)
        ZENITH = GWordered['Zenith angles in steps'][i]
        #print('ZNITHHHHH',len(galPointing['Array of zenith angles'][i]))
        #print('x',len(x))
        x=np.arange(0, len(ZENITH), 1)
        ax.plot(x,ZENITH,label='ra:%.2f dec:%.2f- Pgw:%.3f - Pgal:%.3f ' % (ra[i],dec[i],100*pgw[i],100*pgal[i]))

    ax.set_xticks(x)
    ax.set_xticklabels(hour)
    ax.set_ylabel('Altitude (deg)',fontsize=14)
    ax.set_xlabel('Time',fontsize=14)
    ax.grid()
    ax.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
    plt.savefig("%s/AltitudevsTime_%s.png" %(tname, ObsArray))


def RankingTimes(ObservationTime, filename, cat, obspar, targetType, dirName, PointingFile, ObsArray):


    point = load_pointingFile(PointingFile)

    ################################################################

    print()
    print('---------  RANKING THE OBSERVATIONS AND PRODUCING THE OUTPUT FILES   ----------')
    print()

    print('Loading map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    has3D = True
    if (len(distnorm) == 0):
        has3D = False
    # correlate GW map with galaxy catalog, retrieve ordered list
    tGals, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D, obspar.MinimumProbCutForCatalogue)
    point = ProbabilitiesinPointings3D(tGals, point, obspar.FOV, sum_dP_dV, prob, nside)
    point = VisibilityWindow(ObservationTime, point, 90 - obspar.max_zenith, obspar.MaxNights, obspar.UseGreytime, dirName, obspar.max_zenith,
                             obspar, obspar.gMoonGrey, obspar.gMoonDown, obspar.gMoonPhase, obspar.gSunDown, obspar.MoonSourceSeparation, obspar.MaxMoonSourceSeparation)

    EvolutionPlot(point, dirName, ObsArray)
    Sortingby(point, targetType, dirName, obspar.Duration)

def RankingTimes_SkyMapInput_2D(ObservationTime, prob, obspar, targetType, dirName, PointingFile, ObsArray):

    point = load_pointingFile(PointingFile)

    ################################################################

    print()
    print('---------  RANKING THE OBSERVATIONS AND PRODUCING THE OUTPUT FILES   ----------')
    print()

    npix = len(prob)
    nside = hp.npix2nside(npix)

    point = ProbabilitiesinPointings2D(point, obspar.FOV, prob, nside)
    point = VisibilityWindow(ObservationTime, point, 90 - obspar.max_zenith, obspar.MaxNights, obspar.UseGreytime, dirName, obspar.max_zenith,
                             obspar, obspar.gMoonGrey, obspar.gMoonDown, obspar.gMoonPhase, obspar.gSunDown, obspar.MoonSourceSeparation, obspar.MaxMoonSourceSeparation)

    EvolutionPlot(point, dirName, ObsArray)
    Sortingby(point, targetType, dirName, obspar.Duration)



