import datetime
import os
import astropy.coordinates as co
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun,get_moon
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.data import download_file
import copy
from .GWHESSPointingTools import (LoadHealpixMap, Tools, CorrelateGalaxies_LVC,
                                  CorrelateGalaxies_LVC_SteMass)

from six.moves import configparser
import six
if six.PY2:
  ConfigParser = configparser.SafeConfigParser
else:
  ConfigParser = configparser.ConfigParser

# iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))

iers_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)

'''
    gloably defined darkenss criteria:
    max sun and moon altitude in degrees.
    '''
gSunDown = -18
#HorizonSun = '-17:43:48'
HorizonSun = '-18:00:00'
gMoonDown = -0.5
HorizonMoon = '-00:30:00'
gMoonPhase=60
MoonSourceSeparation=30

class HESSObservatory:
    def __init__(self):
        self.Lat = -23.271778 * u.deg
        self.Lon = 16.50022 * u.deg
        self.Height = 1835 * u.m
        self.Location = EarthLocation(lat=self.Lat, lon=self.Lon,
                                      height=self.Height)


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
    print(ra)
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    time=[]
    raclean=[]
    decclean=[]
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i].split(':')[0]+ ':'+time2[i].split(':')[1]).split('"')[1])


    ra = ra.astype(np.float)
    dec = dec.astype(np.float)

    l = list(range(len(ra)))
    Pointings = Table([l,time,ra,dec],names=('Pointing','Time','RAJ2000', 'DEJ2000'))

    #print(Pointings)
    return Pointings


def VisibilityWindow(ObservationTime,galPointing,AltitudeCut,nights,UseGreytime,dirName):
    source = SkyCoord(galPointing['RAJ2000'],galPointing['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    WINDOW=[]
    ZENITH=[]
    SZENITH=[]
    try:
        auxtime = datetime.datetime.strptime(galPointing['Time'][0], '%Y-%m-%d %H:%M')
    except ValueError:
        auxtime = datetime.datetime.strptime(galPointing['Time'][0], '%Y-%m-%d %H:%M')
    #frame = co.AltAz(obstime=auxtime, location=observatory)
    timeInitial=auxtime-datetime.timedelta(minutes=30)
    for i in range(0,len(source)):
        NonValidwindow,Stepzenith = GetVisibility(galPointing['Time'],source[i])
        window, zenith = GetObservationPeriod(timeInitial, source[i],AltitudeCut,nights,i,dirName,UseGreytime,False)
        print(window)
        WINDOW.append(window)
        ZENITH.append(zenith)
        SZENITH.append(Stepzenith)
        window, zenith = GetObservationPeriod(ObservationTime, source[i],AltitudeCut,nights,i,dirName,UseGreytime,True)

    galPointing['Observation window'] = WINDOW
    galPointing['Array of zenith angles']=ZENITH
    galPointing['Zenith angles in steps'] = SZENITH

    return galPointing

def GetVisibility(time,radecs):
    #Se tiene que ver el FoV entero, no vale solo algo mas del centro MEJORAR ESTE CRITERIOOOOOOOOOOO
    observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1835 * u.m)
    visibility = []
    altitude = []

    for i in range(0,len(time)):
        try:
            auxtime = datetime.datetime.strptime(time[i], '%Y-%m-%d %H:%M')
        except ValueError:
            auxtime = datetime.datetime.strptime(time[i], '%Y-%m-%d %H:%M')
        frame = co.AltAz(obstime=auxtime, location=observatory)
        thisaltaz = radecs.transform_to(frame)
        visible = thisaltaz.alt.value > (90-60)

        if(visible):
            visibility.append(auxtime)
            altitude.append(thisaltaz.alt.value)
        else:
        #    visibility.append(auxtime)
             altitude.append(thisaltaz.alt.value)
    lasttime=auxtime+datetime.timedelta(minutes=30)
    frame = co.AltAz(obstime=lasttime, location=observatory)
    thisaltaz = radecs.transform_to(frame)
    visible = thisaltaz.alt.value > (90 - 60)

    if (visible):
        visibility.append(auxtime)
        altitude.append(thisaltaz.alt.value)
    else:
        #    visibility.append(auxtime)
        altitude.append(thisaltaz.alt.value)

    window = (visibility[0].strftime('%H:%M:%S')).split(':')[0]+':'+(visibility[0].strftime('%H:%M:%S')).split(':')[1] +'-'+ (visibility[-1].strftime('%H:%M:%S')).split(':')[0] +':'+(visibility[-1].strftime('%H:%M:%S')).split(':')[1]
    return window,altitude

def ProbabilitiesinPointings(cat,galPointing,FOV, totaldPdV,prob,nside):
    #targetCoord2 = co.SkyCoord(cat['RAJ2000'], cat['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

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

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoordpointing.dec.rad

    p = targetCoordpointing.ra.rad

    #print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)

    xyz = hp.ang2vec(t, p)


    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    P_GW = prob[ipix_disc].sum()

    Pgal_inFoV = dp_dV[targetCoordcat.separation(targetCoordpointing).deg <= radius].sum() / totaldPdV

    return P_GW,Pgal_inFoV

def Sortingby(galPointing,name):

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
    print(gwgalPointing)
    outfilename='%s/RankingObservationTimes_Complete.txt' % name
    ascii.write(gwgalPointing[np.argsort(gwgalPointing['Pointing'])], outfilename, overwrite=True)
    
    gwgalPointing.remove_column('Pgal')
    gwgalPointing.remove_column('Pgw')
    gwgalPointing.remove_column('RA(HH:MM:SS) Dec (DD:MM:SS)')

    #file to send to shifter in case PGW is used
    gwgalPointingGW = copy.deepcopy(gwgalPointing)
    gwgalPointingGW.remove_column('PriorityGal')
    gwgalPointingGW.rename_column('PriorityGW', 'Priority')
    outfilename='%s/RankingObservationTimesPGW_forShifters.txt' % name
    ascii.write(gwgalPointingGW[np.argsort(gwgalPointingGW['Pointing'])], outfilename, overwrite=True)
    
    #file to send to shifters in case PGal is used
    gwgalPointing.remove_column('PriorityGW')
    gwgalPointing.rename_column('PriorityGal','Priority')
    outfilename='%s/RankingObservationTimesPGal_forShifters.txt' % name
    ascii.write(gwgalPointing[np.argsort(gwgalPointing['Pointing'])], outfilename, overwrite=True)

    #gggalPointing.remove_column('DEJ2000')
    #gggalPointing.remove_column('RAJ2000')

    #gggalPointing['RA(HH:MM:SS) Dec (DD:MM:SS)'] = coord.to_string('hmsdms')
    #outfilename='%s/RankingObservationTimes_Complete.txt' % name
    #ascii.write(gggalPointing[np.argsort(gggalPointing['Pointing'])], outfilename,overwrite=True)



def GetObservationPeriod(inputtime0,msource,AltitudeCut,nights,plotnumber,dirName,UseGreytime,doplot):

    observatory=EarthLocation(lat=-23.271778*u.deg, lon=16.50022*u.deg, height=1835 *u.m)
    inputtime=Time(inputtime0)
    initialframe= AltAz(obstime=inputtime,location=observatory)

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
    #print('After',len(NightsCounter),'vs',len(delta_day))
    times = inputtime + delta_day
    frame= AltAz(obstime=times,location=observatory)

    ##############################################################################
    #SUN
    sunaltazs = get_sun(times).transform_to(frame)
    #MOON
    moon = get_moon(times)
    moonaltazs = moon.transform_to(frame)
    msourcealtazs = msource.transform_to(frame)

    #Add Moon phase
    moonPhase = np.full(len(msourcealtazs),Tools.MoonPhase(inputtime0))
    #print(moonPhase)
    #Moon Distance

    MoonDistance=msourcealtazs.separation(moonaltazs)
    ##############################################################################
    if UseGreytime:
        Altitudes = Table([times,msourcealtazs.alt, sunaltazs.alt, moonaltazs.alt,moonPhase,MoonDistance,NightsCounter],
                           names=['Time UTC', 'Alt Source', 'Alt Sun', 'AltMoon','moonPhase','MoonDistance','NightsCounter'])
        #selectedTimes=Altitudes['Time UTC']
        selection=[(Altitudes['Alt Sun'] < -18.) & (Altitudes['Alt Source']> AltitudeCut) & (Altitudes['AltMoon'] < -0.5)]
        DTaltitudes=Altitudes[selection]
        newtimes=[]
        newtimes.extend(DTaltitudes['Time UTC'].mjd)
        selectionGreyness= [(Altitudes['AltMoon'] < 50)&(Altitudes['AltMoon'] > -0.5)&(Altitudes['moonPhase'] <gMoonPhase)&(Altitudes['Alt Sun'] < -18.) & (Altitudes['MoonDistance']>MoonSourceSeparation)&(Altitudes['Alt Source']> AltitudeCut)]
        GTaltitudes=Altitudes[selectionGreyness]
        newtimes.extend(GTaltitudes['Time UTC'].mjd)
        newtimes=sorted(newtimes)
        ScheduledTimes=Time(newtimes,format='mjd').iso
    else:
        Altitudes = Table([times,msourcealtazs.alt, sunaltazs.alt, moonaltazs.alt,NightsCounter],
                           names=['Time UTC', 'Alt Source', 'Alt Sun', 'AltMoon','NightsCounter'])
        Times=Altitudes['Time UTC']
        selection=[(Altitudes['Alt Sun'] < -18.) & (Altitudes['Alt Source']> AltitudeCut) & (Altitudes['AltMoon'] < -0.5)]
        ScheduledTimes=Times[selection]

    #print(ScheduledTimes)
    #doplot=True
    if doplot:
        plotDir = '%s/TransitPlots' % dirName
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        plt.figure(figsize=(20, 12))
        plt.plot(delta_day, sunaltazs.alt, color='r', label='Sun')
        plt.plot(delta_day, moonaltazs.alt, color=[0.75]*3, ls='--', label='Moon')
        plt.scatter(delta_day, msourcealtazs.alt,c=msourcealtazs.az, label='Point source', lw=0, s=8,cmap='viridis')
        plt.fill_between(delta_day.to('hr').value, 0, 90,sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_day.to('hr').value, 0, 90,(sunaltazs.alt < -18*u.deg)&(moonaltazs.alt < -0.5*u.deg), color='k', zorder=0)
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        #x=np.arange(0, len(ZENITH), 1)
        #ax.set_xticks(x)
        #ax.set_xticklabels(hour)
        #print('hour axis',inputtime.tt.datetime.hour+delta_day.to('hr').value)
        #plt.xlim(inputtime.tt.datetime.hour+delta_day.to('hr').value)
        #plt.xlim(0,hoursinDay+24*(nights-1))
        #plt.xticks(ScheduledTimes['Times UTC'])
        plt.ylim(0, 90)
        #plt.xlabel('Hours after injections')
        plt.ylabel('Altitude [deg]')
        plt.grid()
        print("plotDir, plotnumber",plotDir, plotnumber)
        plt.savefig('%s/Source%g.png' % (plotDir, plotnumber))

    return (str(ScheduledTimes[0]).split('.')[0]+'-->'+str(ScheduledTimes[-1]).split('.')[0]),msourcealtazs.alt


def EvolutionPlot(galPointing,tname):
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
        print(time[j])
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
    plt.savefig("%s/ZenithvsTime.png"%tname)

def RankingTimes(ObservationTime,filename,cat,parameters,dirName,PointingFile):
    # Main parameters

    ##################
    cfg = parameters
    parser = ConfigParser()
    print(parser.read(cfg))
    print(parser.sections())
    section = 'GWBestGalaxyParameters'

    try:
        max_zenith = int(parser.get(section, 'max_zenith'))
        nights = int(parser.get(section, 'MaxNights'))
        FOV = float(parser.get(section, 'FOV'))
        MinimumProbCutForCatalogue = float(parser.get(section, 'MinimumProbCutForCatalogue'))
        UseGreytime = (parser.getboolean(section, 'UseGreytime'))
        Mangrove = (parser.getboolean(section, 'Mangrove'))

    except Exception as x:
        print(x)

    #print(UseGreytime)
    #########################

    # Load galaxy catalog
    #TODO: need variable location for this!!
    galFile = '/Users/mseglar/Documents/CurrentPhD/HESS/GW/GLADE/GLADE_2.3clean.txt'
    #cat = LoadGalaxies(galFile)
    print('done loading galaxies')

    #PointingFile = '/Users/mseglar/Documents/CurrentPhD/HESS/GW/Follow-up/ContinuePointing/G275404_SuggestedPointings_GALOptimisation.txt'
    #PointingFile = '/Users/mseglar/Documents/GitHub/GWfollowup/'
    point = load_pointingFile(PointingFile)

    #filename = '/Users/mseglar/Documents/CurrentPhD/HESS/GW/Follow-up/ContinuePointing/G275404_bayestar.fits'

    ################################################################

    #if ('G' in filename):
    #    names = filename.split("_")
    #    name = names[0]
    name = filename.split('.')[0].split('/')[-1]

    print()
    print('-------------------   NEW LVC EVENT   --------------------s')
    print()

    print('Loading GW map from ', filename)
    prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    has3D=True
    if(len(distnorm)==0):
        print("Found a generic map without 3D information")
    # flag the event for special treatment
        has3D=False
    else:
        print("Found a 3D reconstruction")
    #hp.mollview(prob, title="GW prob map (Ecliptic)")
        # plt.show
    
    # correlate GW map with galaxy catalog, retrieve ordered list
    if not Mangrove:
        tGals, sum_dP_dV = CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, has3D,MinimumProbCutForCatalogue)
    else:
        tGals, sum_dP_dV = CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, thisDistance, thisDistanceErr, distnorm, cat, has3D,MinimumProbCutForCatalogue)

    point = ProbabilitiesinPointings(tGals,point,FOV,sum_dP_dV,prob,nside)
    point = VisibilityWindow(ObservationTime,point,90-max_zenith,nights,UseGreytime,dirName)

    EvolutionPlot(point,dirName)
    Sortingby(point,dirName)


