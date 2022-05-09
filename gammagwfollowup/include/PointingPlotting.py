from .PointingTools import (LoadHealpixMap, TransformRADec,
                                  CTANorthObservatory,CTASouthObservatory,ObservationParameters)
import os
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.utils import iers
import astropy.coordinates as co
from astropy.coordinates import SkyCoord
import datetime
from six.moves import configparser
import six
if six.PY2:
  ConfigParser = configparser.SafeConfigParser
else:
  ConfigParser = configparser.ConfigParser


# iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

iers_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)

# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))

colors = ['b', 'm', 'y', 'c', 'g', 'w', 'k', 'c', 'b', 'c', 'm', 'b', 'g', 'y', 'b', 'c', 'm', 'b', 'g', 'y',
          'b', 'c', 'm', 'b']



def LoadPointingsGW(tpointingFile):

    print("Loading pointings from " + tpointingFile)

    time1,time2, ra, dec = np.genfromtxt(tpointingFile, usecols=(0, 1, 2,3), dtype="str", skip_header=1,
                                             delimiter=' ',
                                             unpack=True)  # ra, dec in degrees
    time1 = np.atleast_1d(time1)
    time2 = np.atleast_1d(time2)
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    time = []
    for i in range(len(time1)):
        try:
            time.append((time1[i] + ' ' + time2[i]).split('"')[1])
        except IndexError:
            time.append((time1 + ' ' + time2).split('"')[1])
            break

    ra = ra.astype(np.float)
    dec = dec.astype(np.float)
    coordinates = co.SkyCoord(ra,dec, frame='fk5', unit=(u.deg, u.deg))
    #pgw = Pgw.astype(np.float)
    pgw = np.genfromtxt(tpointingFile, usecols=5, skip_header=1,delimiter=' ',unpack=True)

    return time, coordinates , pgw


def LoadPointingsGAL(tpointingFile):


    print("Loading pointings from " + tpointingFile)
    time1, time2, ra ,dec,Pgw, Pgal = np.genfromtxt(tpointingFile, usecols=(0, 1, 2,3,4,5),dtype= "str" , skip_header=1,delimiter = ' ',
                                        unpack=True)  # ra, dec in degrees
    time1 = np.atleast_1d(time1)
    time2 = np.atleast_1d(time2)
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    time = []
    for i in range(len(time1)):
        try:
            time.append((time1[i] + ' ' + time2[i]).split('"')[1])
        except IndexError:
            time.append((time1 + ' ' + time2).split('"')[1])
            break
    coordinates=TransformRADec(ra, dec)
    print(Pgw)
    Pgw = Pgw.astype(np.float)
    Pgal = Pgal.astype(np.float)
    return time, coordinates ,Pgw, Pgal


def PointingPlotting(prob, parameters, name,dirName,PointingsFile1):

    # Main parameters
    obspar = ObservationParameters.from_configfile(parameters)

    npix = len(prob)
    nside = hp.npix2nside(npix)

    ObservationTimearray1, Coordinates1, Probarray1 = LoadPointingsGW(PointingsFile1)

    Probarray1 = np.atleast_1d(Probarray1)

    print('----------   PLOTTING THE SCHEDULING   ----------')
    print('Total probability of map 1 that maximises PGW= {0:.5f}'.format(sum(Probarray1)))
    converted_time1=[]
    for i in range(0,len(ObservationTimearray1)):
        time1 = ObservationTimearray1[i]
        try:
            converted_time1.append(datetime.datetime.strptime(time1, '%Y-%m-%d %H:%M:%S.%f'))
        except ValueError:
            try:
                converted_time1.append(datetime.datetime.strptime(time1, '%Y-%m-%d %H:%M:%S.%f '))
            except ValueError:
                converted_time1.append(datetime.datetime.strptime(time1, '%Y-%m-%d %H:%M:%S'))



    #PlotPointingsTogether(prob,converted_time1[0],Coordinates1,sum(Probarray1),name1,Coordinates2,sum(Probarray2),name2, nside, obspar.FOV, doplot=True)
    PlotPointings(prob,converted_time1,Coordinates1,sum(Probarray1), nside, obspar, name, dirName)


def PlotPointings(prob, time, targetCoord, Totalprob, nside, obspar, name, dirName):

    FOV = obspar.FOV
    maxzenith= obspar.max_zenith
    doplot = obspar.doplot
    radius = FOV

    t = 0.5 * np.pi - targetCoord[0].dec.rad
    p = targetCoord[0].ra.rad

    xyz = hp.ang2vec(t, p)


    # translate pixel indices to coordinates
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    if (doplot):

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        observatory = obspar.Location
        frame = co.AltAz(obstime=time[0], location=observatory)
        altaz_all = skycoord.transform_to(frame)

        dirName = '%s/Pointing_Plotting' % dirName
        if not os.path.exists(dirName):
            os.makedirs(dirName)


        hp.mollview(prob, rot=[180, 0], coord='C', title="GW prob map (Ecliptic) + %s %g  %s/%s/%s %s:%s:%s UTC" %
                                                         (name, Totalprob * 100, time[0].day, time[0].month,
                                                          time[0].year,
                                                          time[0].hour, time[0].minute, time[0].second))
        hp.graticule()
        # plt.show()
        theta = np.random.rand(400) * 360
        tarcoordra = np.empty(400)
        tarcoorddec = np.empty(400)
        Fov_array = np.empty(400)
        Fov_array.fill(FOV)

        # Position of the GRB: Add here your coordinates if you want to add them to the plot
        #RA_GRB = 0.526
        #DEC_GRB = -2.917
        #skycoordGRB = co.SkyCoord(RA_GRB, DEC_GRB, frame='fk5', unit=(u.deg, u.deg))

        for j in range(0, len(targetCoord.ra)):
            tarcoordra.fill(targetCoord[j].ra.deg)
            tarcoorddec.fill(targetCoord[j].dec.deg)
            racoord = tarcoordra + Fov_array * np.cos(theta)
            deccoord = tarcoorddec + Fov_array * np.sin(theta)
            hp.visufunc.projscatter(racoord, deccoord, lonlat=True, marker='.', color=colors[j], coord='C')
            #hp.visufunc.projscatter(RA_GRB, DEC_GRB, lonlat=True, marker='+', color=colors[1], coord='C')

            # Plotting the visibility of the instrument as specified in the config

            altcoord = np.empty(1000)
            altcoord.fill(90 - maxzenith)
            azcoord = np.random.rand(1000) * 360
            # print(time)
            RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time[j],
                                   location=observatory)
            RandomCoord_radec = RandomCoord.transform_to('fk5')
            hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True)
            #plt.show()
            plt.savefig('%s/Pointings%s.png' % (dirName, j))

        #dist = cat['Dist']
        #hp.visufunc.projscatter(cat['RAJ2000'][dist < 200], cat['DEJ2000'][dist < 200], lonlat=True, marker='.',
        #                        color='g', linewidth=0.1, coord='C')
        # draw all galaxies within zenith-angle cut

        # draw observation position, which is equivalent to galaxy with highest
        # probability

        # hp.visufunc.projscatter(targetCoord.ra.deg, targetCoord.dec.deg, lonlat=True, marker='.', color=colors[1],coord=['E','G'])

        # plt.savefig("Pointing_Plotting/G274296Pointing_Comparison_%g_%g:%g.png" % (time.day, time.hour, time.minute))
        # draw circle of HESS-I FoV around best fit position

        # hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'r.', lonlat=True, coord="E")

        # Draw H.E.S.S. FOV

        # dec = np.empty(500)
        # dec.fill(10)
        # ra = np.random.rand(500) * 360
        # hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        # dec = np.empty(500)
        # dec.fill(-10)
        # hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        # dec = np.random.rand(500) * 20 - 10
        # ra = np.empty(500)
        # ra.fill(130)
        # hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        # ra = np.empty(500)
        # ra.fill(-100)
        # hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        '''
        ra = np.random.rand(500) * 130
        dec = ra * (-10 / 130) + 10
        hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        ra = np.random.rand(500) * 120 - 120
        dec = ra * (10 / 120) + 10
        hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        ra = np.random.rand(500) * 130
        dec = ra * (10 / 130) - 10
        hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')

        ra = np.random.rand(500) * 120 - 120
        dec = ra * (-10 / 120) - 10
        hp.visufunc.projscatter(ra, dec, lonlat=True, marker='.', color='w', coord='G')



        altcoord = np.empty(1000)
        altcoord.fill(45)
        azcoord = np.random.rand(1000) * 360
        #print(time)
        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time[0],
                               location=observatory)
        RandomCoord_radec = RandomCoord.transform_to('fk5')
        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord=['E','G'])

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time[-1],
                               location=observatory)
        RandomCoord_radec = RandomCoord.transform_to('fk5')
        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord=['E','G'])

        #plt.show()

        #LOAD galaxies
        c_icrs= SkyCoord(ra=cat['RAJ2000'][dist<200]*u.degree, dec=cat['DEJ2000'][dist<200]*u.degree, frame='icrs')
        subset0=c_icrs[(np.absolute(c_icrs.galactic.l.value-180)>0)&(np.absolute(c_icrs.galactic.l.value-180)<30)]
        subset1 = c_icrs[(np.absolute(c_icrs.galactic.l.value-180)>30)&(np.absolute(c_icrs.galactic.l.value-180)<60)]
        subset2 = c_icrs[(np.absolute(c_icrs.galactic.l.value-180)>60)&(np.absolute(c_icrs.galactic.l.value-180)<90)]
        subset3 = c_icrs[(np.absolute(c_icrs.galactic.l.value-180)>90)&(np.absolute(c_icrs.galactic.l.value-180)<120)]
        subset4 = c_icrs[(np.absolute(c_icrs.galactic.l.value-180)>120)&(np.absolute(c_icrs.galactic.l.value-180)<150)]
        subset5 = c_icrs[(np.absolute(c_icrs.galactic.l.value-180) > 150) & (np.absolute(c_icrs.galactic.l.value-180) < 180)]

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        #plt.hist(c_icrs.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g')
        #plt.hist(c_icrs.galactic.b.value,90,fill=False,facecolor='b',stacked=True,histtype='step')

        #plt.hist(subset0.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g')
        plt.hist(subset0.galactic.b.value,40,fill=False,facecolor='b',stacked=True,histtype='step',label='|l-180|<30')

        #plt.hist(subset1.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g')
        plt.hist(subset1.galactic.b.value,40,fill=False,facecolor='b',stacked=True,histtype='step',label='30<|l-180|<60')

        #plt.hist(subset2.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g')
        plt.hist(subset2.galactic.b.value,40,fill=False,facecolor='b',stacked=True,histtype='step',label='60<|l-180|<90')
        plt.hist(subset3.galactic.b.value,40,fill=False,facecolor='b',stacked=True,histtype='step',label='90<|l-180|<120')
        plt.hist(subset4.galactic.b.value,40,fill=False,facecolor='b',stacked=True,histtype='step',label='120<|l-180|<150')
        plt.hist(subset5.galactic.b.value,40,fill=False,facecolor='b',stacked=True,histtype='step',label='150<|l-180|<180')

        # Major ticks every 20, minor ticks every 5
        #major_ticks = np.arange(-200, 200, 50)
        major_ticks = np.arange(-50, 50, 10)
        minor_ticks = np.arange(-50, 50, 5)

        #ax.set_xticks(major_ticks)
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)

        #ax.set_yticks(major_ticks)
        #ax.set_yticks(minor_ticks, minor=True)
        ax.legend(loc='lower left')
        ax.set_xlabel('Galactic latitude')
        ax.set_ylabel('#')
        #ax.set_ylim(0, 4000)
        ax.set_xlim(-50, 50)
        ax.set_yscale("log")
        # And a corresponding grid
        ax.grid(which='major', alpha=0.2)
        ax.grid(which='minor', alpha=0.2)
        #ax.grid(which='major', alpha=0.5)
        plt.savefig("Latitude-Longitude.png")

        #include them in the plot
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.hist(c_icrs.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='l 0-30')
        plt.hist(subset0.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='Lat 0-30')
        plt.hist(subset1.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='Lat 30-60')
        plt.hist(subset2.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='Lat 60-90')
        plt.hist(subset3.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='Lat 90-120')
        plt.hist(subset4.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='Lat 120-150')
        plt.hist(subset5.galactic.l.value-180,90,histtype='step', stacked=True, fill=False,facecolor='g',label='Lat 120-150')

        ax.legend()

        #plt.show()

        plt.savefig("%s/Pointing_Galaxies.png"% dirName)
        '''


def PlotPointingsTogether(prob, time, targetCoord1, n1, targetCoord2, n2, nside, FOV, doplot=True):
    radius = FOV

    t = 0.5 * np.pi - targetCoord1[0].dec.rad
    p = targetCoord1[0].ra.rad

    # print('t, p, targetCoord1[0].ra.deg, targetCoord1[0].dec.deg',t, p, targetCoord1[0].ra.deg, targetCoord1[0].dec.deg)

    xyz = hp.ang2vec(t, p)

    # print(xyz)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    if (doplot):

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        # skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

        # frame = co.AltAz(obstime=time, location=observatory)
        # altaz_all = skycoord.transform_to(frame)

        # path = os.path.dirname(os.path.realpath(__file__)) + '/Pointing_Plotting'
        # if not os.path.exists(path):
        #    os.mkdir(path, 493)
        hp.gnomview(prob, xsize=500, ysize=500, rot=[targetCoord1[0].ra.deg, targetCoord1[0].dec.deg], reso=5.0)
        # hp.mollview(prob,title="GW prob map (Ecliptic)",coord='C')
        hp.graticule()

        hp.visufunc.projscatter(targetCoord1.ra.deg, targetCoord1.dec.deg, lonlat=True, marker='.', color=colors[4])
        hp.visufunc.projscatter(targetCoord2.ra.deg, targetCoord2.dec.deg, lonlat=True, marker='.', color=colors[5])

        # plt.savefig("Pointing_Plotting/G274296Pointing_Comparison_%g_%g:%g.png" % (time.day, time.hour, time.minute))
        # draw circle of HESS-I FoV around best fit position

        # hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'r.', lonlat=True, coord="E")

        # Draw H.E.S.S. FOV
        Fov_array = np.empty(400)
        Fov_array.fill(FOV)

        theta = np.random.rand(400) * 360
        tarcoordra1 = np.empty(400)

        tarcoorddec1 = np.empty(400)
        for j in range(0, len(targetCoord1.ra)):
            tarcoordra1.fill(targetCoord1[j].ra.deg)
            tarcoorddec1.fill(targetCoord1[j].dec.deg)
            racoord1 = tarcoordra1 + Fov_array * np.cos(theta)
            deccoord1 = tarcoorddec1 + Fov_array * np.sin(theta)
            hp.visufunc.projscatter(racoord1, deccoord1, lonlat=True, marker='.', color=colors[2])
            hp.projtext(targetCoord1[j].ra, targetCoord1[j].dec, str(j), lonlat=True, color=colors[2])
        # plt.savefig("Pointing_Plotting/G274296Pointing_GW_%g_%g:%g.png" % (time.day, time.hour,time.minute))

        tarcoordra2 = np.empty(400)

        tarcoorddec2 = np.empty(400)
        for j in range(0, len(targetCoord2.ra)):
            tarcoordra2.fill(targetCoord2[j].ra.deg)
            tarcoorddec2.fill(targetCoord2[j].dec.deg)
            racoord2 = tarcoordra2 + Fov_array * np.cos(theta)
            deccoord2 = tarcoorddec2 + Fov_array * np.sin(theta)
            hp.visufunc.projscatter(racoord2, deccoord2, lonlat=True, marker='.', color=colors[1])
            hp.projtext(targetCoord1[j].ra, targetCoord1[j].dec, str(j), lonlat=True, color=colors[1])
        # plt.show()
        # plt.savefig("Pointing_Plotting/PointingFOV_Comparison.png")


def PointingPlottingGWCTA(filename,name, dirName, PointingsFile,FOV,UseObs):

    print()
    print('-------------------   PLOTTING SCHEDULE   --------------------')
    print()

    # Observatory
    if UseObs == 'South':
        print('Observed form the', UseObs)
        observatory = CTASouthObservatory()
    else:
        print('Observed from the', UseObs)
        observatory = CTANorthObservatory()

    if ('G' in filename):
        names = filename.split("_")
        name = names[0]


    print('Loading map from ', filename)
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

    ObservationTimearray, Coordinates, Probarray = LoadPointingsGW(PointingsFile)

    Probarray = np.atleast_1d(Probarray) #making sure to have an array on which we can call the sum() function

    print('----------   BUILDING A MAP   ----------')
    print('Total probability of map 1 that maximises PGW= {0:.5f}'.format(sum(Probarray)))

    converted_time=[]
    for i in range(0,len(ObservationTimearray)):
        time = ObservationTimearray[i]
        try:
            converted_time.append(datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f'))
        except ValueError:
            try:
                converted_time.append(datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f '))
            except ValueError:
                converted_time.append(datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S'))


    '''
    converted_time2=[]
    for i in range(0,len(ObservationTimearray2)):
        time2 = ObservationTimearray2[i]
        try:
            converted_time2.append(datetime.datetime.strptime(time2, '%Y-%m-%d %H:%M:%S.%f '))
        except ValueError:
            converted_time2.append(datetime.datetime.strptime(time2, '%Y-%m-%d %H:%M:%S '))
    '''

    #PlotPointings(prob,cat,converted_time,Coordinates,sum(Probarray), nside, FOV, name, dirName, doplot=True)
    radius = FOV

    t = 0.5 * np.pi - Coordinates[0].dec.rad
    p = Coordinates[0].ra.rad
    xyz = hp.ang2vec(t, p)

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    tt, pp = hp.pix2ang(nside, ipix_disc)
    ra2 = np.rad2deg(pp)
    dec2 = np.rad2deg(0.5 * np.pi - tt)

    skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))


    frame = co.AltAz(obstime=converted_time[0], location=observatory.Location)
    altaz_all = skycoord.transform_to(frame)

    dirName = '%s/Pointing_Plotting/%s' % (dirName,name)
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    #path = os.path.dirname(os.path.realpath(__file__)) + '/Pointing_Plotting'
    #if not os.path.exists(path):
    #    os.mkdir(path, 493)

    hp.mollview(prob, rot=[180,0],coord='C', title="GW prob map (Ecliptic) + %s %g  %s/%s/%s %s:%s:%s UTC" %
                                                  (name, sum(Probarray) * 100, converted_time[0].day, converted_time[0].month, converted_time[0].year,
                                                   converted_time[0].hour, converted_time[0].minute, converted_time[0].second))
    hp.graticule()
    #plt.show()
    theta = np.random.rand(400) * 360
    tarcoordra = np.empty(400)
    tarcoorddec = np.empty(400)
    Fov_array = np.empty(400)
    Fov_array.fill(FOV)


    for j in range(0, len(Coordinates.ra)):
        tarcoordra.fill(Coordinates[j].ra.deg)
        tarcoorddec.fill(Coordinates[j].dec.deg)
        racoord = tarcoordra + Fov_array * np.cos(theta)
        deccoord = tarcoorddec + Fov_array * np.sin(theta)
        hp.visufunc.projscatter(racoord, deccoord, lonlat=True, marker='.', color=colors[1],coord='C')
        #plt.show()
    plt.savefig('%s/Pointings.png'%dirName)

    #dist = cat['Dist']
    #hp.visufunc.projscatter(cat['RAJ2000'][dist < 200], cat['DEJ2000'][dist < 200], lonlat=True, marker='.',color='g', linewidth=0.1, coord='C')


def PointingPlottingGW_ZenithSteps(filename,name,dirName,FOV,InputTimeObs):

    print()
    print('-------------------   PLOTTING SCHEDULE   --------------------')
    print()

    UseObs = InputTimeObs['Observatory']
    time = InputTimeObs['Time']

    # Observatory
    if UseObs == 'South':
        print('Observed form the', UseObs)
        observatory = CTASouthObservatory()
    else:
        print('Observed from the', UseObs)
        observatory = CTANorthObservatory()

    if ('G' in filename):
        names = filename.split("_")
        name = names[0]


    print('Loading map from ', filename)
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

    print('----------   BUILDING A MAP   ----------')

    try:
        converted_time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f')
    except ValueError:
        try:
            converted_time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f ')
        except ValueError:
            converted_time = datetime.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')

    radius = FOV

    #t = 0.5 * np.pi - Coordinates[0].dec.rad
    #p = Coordinates[0].ra.rad
    #xyz = hp.ang2vec(t, p)

    #ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    #tt, pp = hp.pix2ang(nside, ipix_disc)
    #ra2 = np.rad2deg(pp)
    #dec2 = np.rad2deg(0.5 * np.pi - tt)

    #skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))


    #frame = co.AltAz(obstime=converted_time, location=observatory.Location)
    #altaz_all = skycoord.transform_to(frame)

    dirName = '%s/Pointing_Plotting/%s' % (dirName,name)
    if not os.path.exists(dirName):
        os.makedirs(dirName)

    max_zenith=45
    hp.mollview(prob)
    for i in range(0,6):
        altcoord = np.empty(2500-200*i)
        azcoord = np.random.rand(2500-200*i) * 360
        altcoord.fill(max_zenith+5*i)
        print(max_zenith+5*i)
        RandomCoord = co.SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory.Location)
        RandomCoord_radec = RandomCoord.transform_to('fk5')
        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, lonlat=True)
    #plt.show()
    hp.graticule()
    plt.savefig('%s/Pointings.png'%dirName)