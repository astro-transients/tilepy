import os
from .GWHESSPointingTools import Tools
from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel, PowerLaw
from .ObservingTimes import ObtainSingleObservingTimes
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.table import Table
from astropy import units as u
from astropy.utils import iers
import astropy.coordinates as co
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import get_moon
from astropy.io import fits
import datetime
import numpy.ma as ma

#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))
# iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iersf/finals2000A.all'
# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

# iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))
# iers_file = './gammagwfollowup/finals2000A.all'


iers_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)



def GiveProbToGalaxy(prob,cat,distance,Edistance_max,Edistance_min,MinimumProbCutForCatalogue):
    #Cut cat to a cube in distance!

    dist=cat['Dist']
    mask1=dist<Edistance_max
    mask2=dist>Edistance_min
    subcat=cat[mask1&mask2]
    #print(len(subcat),len(cat),Edistance_max,Edistance_min)
    ra=subcat['RAJ2000']
    dec=subcat['DEJ2000']
    # Translate RA,Dec of galaxies into theta,phi angles
    
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    
    # Get corresponding healpix pixel IDs
    
    npix = len(prob)
    nside = hp.npix2nside(npix)
    ipix = hp.ang2pix(nside, theta, phi)
    
    # Calculate probability in the space volumes
    
    pixarea = hp.nside2pixarea(nside)
    
    #Give probability to galaxy
    dp_dV = prob[ipix] / pixarea

    #plt.hist(dp_dV/dp_dV.sum(),histtype='step', stacked=True,  fill=False,cumulative=True,density=True)
    #plt.savefig('/Users/mseglar/Documents/GitHub/CTASimulationsGW/3DComparison_plots/BarbaraConvolution.png')
    TotalProb=dp_dV.sum()
    subcat['dp_dV'] = dp_dV
    min_prob_cut = dp_dV > MinimumProbCutForCatalogue * max(dp_dV)
    Gals = subcat[min_prob_cut]
    # return array with list of Galaxies passing cuts, ordered by p-value

    tGals = Gals[np.flipud(np.argsort(Gals['dp_dV']))]
    
    return tGals,TotalProb


def LoadGalaxies_GladeCTASimu(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    name, dist, z, ra, dec, flag = np.genfromtxt(tgalFile, usecols=(0, 1, 2, 3, 4, 5), skip_header=3, unpack=True,dtype='str')  # ra, dec in degrees
    ra=ra.astype(np.float)
    dec=dec.astype(np.float)
    z=z.astype(np.float)
    dist = dist.astype(np.float)/1000  # change to Mpc!
    tcat = Table([name, ra, dec, dist, z, flag], names=('Galaxy','RAJ2000', 'DEJ2000', 'Dist', 'z', 'flag'))
    #print(tcat)
    return tcat
def LoadGalaxiesSimulation(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    dist, z, ra, dec = np.genfromtxt(tgalFile, usecols=(1, 2, 3,4), skip_header=3, unpack=True)  # ra, dec in degrees

    tcat = Table([ra, dec, dist], names=('RAJ2000', 'DEJ2000', 'Dist'))
    return tcat


def PointingFileReadCTA(pointingFile):

    time1,time2,RA, Dec, Observatory, ZenIni, ZenEnd, Duration = np.genfromtxt(pointingFile,usecols=(0, 1, 2, 3, 4, 6, 7),skip_header=1, unpack=True,dtype='str')
    time=[]
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([time, RA, Dec, Observatory, ZenIni, ZenEnd, Duration],names=['Observation Time UTC','RA(deg)','DEC(deg)','Observatory','ZenIni','ZenEnd','Duration'])
    return OutputTable

def TableImportCTA_simple(inputFileName):
    run,MergerID,RA,Dec,distance,z,theta,ndet,SNR,A90,A50 = np.genfromtxt(inputFileName,usecols=(0, 1, 2, 3,4,5,6,7,8,9,10),skip_header=3, unpack=True,dtype='str')
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    z=z.astype(np.float)
    distance = distance.astype(np.float)
    theta=theta.astype(np.float)
    ndet=ndet.astype(np.float)
    SNR=SNR.astype(np.float)
    A90=A90.astype(np.float)
    A50=A50.astype(np.float)
    OutputTable= Table([run,MergerID,RA,Dec,distance,z,theta,ndet, SNR,A90,A50],names=('run','MergerID','RA','Dec','Distance','redshift','theta','ndet','SNR','A90','A50'))
    return OutputTable


def TableImportCTA(tgalFile):
    run,MergerID,RA,Dec,distance,distMin,distMax,z,theta,ndet,SNR,A90,A50 = np.genfromtxt(tgalFile,usecols=(0, 1, 2, 3,4,5,6,7,8,9,10,11,12),skip_header=1, unpack=True,dtype='str')
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    z=z.astype(np.float)
    distance = distance.astype(np.float)/1000 #to Mpc!
    distMax = distMax.astype(np.float)/1000 #to Mpc!
    distMin= distMin.astype(np.float)/1000 #to Mpc!
    theta=theta.astype(np.float)
    ndet=ndet.astype(np.float)
    SNR=SNR.astype(np.float)
    A90=A90.astype(np.float)
    A50=A50.astype(np.float)
    OutputTable= Table([run,MergerID,RA,Dec,distance,distMin,distMax,z,theta,ndet, SNR,A90,A50],names=('run','MergerID','RA','Dec','Distance','DistMin','DistMax','redshift','theta','ndet','SNR','A90','A50'))
    #print(OutputTable)
    return OutputTable

def TableImportCTA_Glade(tgalFile):
    run,gal,MergerID,RA,Dec,distance,z,theta,ndet,SNR,A90,A50 = np.genfromtxt(tgalFile,usecols=(0, 1, 2, 3,4,5,6,7,8,9,10,11),skip_header=2, unpack=True,dtype='str')
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    z=z.astype(np.float)
    distance = distance.astype(np.float)/1000 #to Mpc!
    theta=theta.astype(np.float)
    ndet=ndet.astype(np.float)
    SNR=SNR.astype(np.float)
    A90=A90.astype(np.float)
    A50=A50.astype(np.float)
    OutputTable= Table([run,gal,MergerID,RA,Dec,distance,z,theta,ndet, SNR,A90,A50],names=('run','Galaxy','MergerID','RA','Dec','Distance','redshift','theta','ndet','SNR','A90','A50'))
    #print(OutputTable)
    return OutputTable

def TableImportCTA_TimeNoZenith(ttimeFile):

    run,MergerID,time1,time2,Observatory = np.genfromtxt(ttimeFile,usecols=(0, 1, 2, 3,4),skip_header=1, unpack=True,dtype='str')
    time=[]
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([run,MergerID,time,Observatory],names=['run', 'MergerID', 'Time', 'Observatory'])

    return OutputTable

def TableImportCTA_SetOfTimes(ttimeFile):
    trial,run, MergerID,time1,time2,Observatory = np.genfromtxt(ttimeFile,usecols=(0, 1, 2, 3,4,5),skip_header=1, unpack=True,dtype='str')
    time=[]
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([run,MergerID,trial,time,Observatory],names=['run', 'MergerID','trial', 'Time', 'Observatory'])

    return OutputTable


def TableImportCTA_Time(ttimeFile):

    run,MergerID,time1,time2,MeanAlt,Observatory = np.genfromtxt(ttimeFile,usecols=(0, 1, 2, 3,4,5),skip_header=1, unpack=True,dtype='str')
    MeanAlt=MeanAlt.astype(np.int)
    time=[]
    for i in range(len(time1)):
        time.append((time1[i] + ' ' + time2[i]).split('"')[1])
    OutputTable = Table([run,MergerID,time,MeanAlt,Observatory],names=['run', 'MergerID', 'Time', 'Zenith', 'Observatory'])

    return OutputTable

def TableImportCTA_Obs(tobsFile):
    pointing,tI,tF,interval = np.genfromtxt(tobsFile,usecols=(0, 1, 2,3),skip_header=3, unpack=True,dtype='int')
    OutputTable = Table([pointing,tI,tF,interval],names=['pointingNumber', 'tI', 'tF', 'Interval'])
    return OutputTable


def TableImportCTA_LS(tgalFile):
    eventid,RA,Dec,distance = np.genfromtxt(tgalFile,usecols=(0,3,4,8),skip_header=1, unpack=True,dtype='str')
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    distance = distance.astype(np.float)/1000 #to Mpc!
    OutputTable= Table([eventid,RA,Dec,distance],names=('MergerID','RA','Dec','Distance'))
    return OutputTable

def IsSourceInside(Pointings,HESS_Sources,FOV,nside):
    tt = 0.5 * np.pi - HESS_Sources.dec.rad
    tp = HESS_Sources.ra.rad
    # print('t, p, targetCoord1[0].ra.deg, targetCoord1[0].dec.deg',t, p, targetCoord1[0].ra.deg, targetCoord1[0].dec.deg)
    txyz = hp.ang2pix(nside, tt, tp)
    # hipix = hp.query_disc(nside, txyz, np.deg2rad(FOV))
    #Npoiting = -1
    Npoiting=[]
    Found = False
    for i in range(0, len(Pointings)):
        #t = 0.5 * np.pi - Pointings[i].dec.rad[0]
        #p = Pointings[i].ra.rad[0]
        t = 0.5 * np.pi - Pointings[i].dec.rad
        p = Pointings[i].ra.rad
        #print('targetCoord1[0].ra.deg, targetCoord1[0].dec.deg',HESS_Sources.ra.deg, HESS_Sources.dec.deg, Pointings[i].ra.deg, Pointings[i].dec.deg)
        xyz = hp.ang2vec(t, p)
        try:
            #print(xyz)
            ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(FOV))
        except:
            ipix_disc = hp.query_disc(nside, xyz[0], np.deg2rad(FOV))
            #print('Smething went wrong')
        #print(txyz)
        #print(ipix_disc)
        #print(txyz in ipix_disc)
        if (txyz in ipix_disc):
            print('Found in pointing number', i)
            Npoiting.append(i)
            Found=True
    return Found, Npoiting

def ProduceSummaryFile(InputList,InputObservationList,allPossiblePoint,foundIn,j,typeSimu,totalProb,datasetDir,outDir,name):
    #print(InputList['run'][j])
    #print(InputList['MergerID'][j].split('r')[-1])
    filepath= datasetDir +'/GammaCatalogV2.0/'+str(InputList['run'][j]) + '_' + str(InputList['MergerID'][j].split('r')[-1]) + ".fits"
    fitsfile = fits.open(filepath)
    luminosity = fitsfile[0].header['EISO']
    dirNameFile = outDir +'/SummaryFile'
    if not os.path.exists(dirNameFile):
        os.makedirs(dirNameFile)
    # Obtain the luminosity
    if foundIn ==-1:
        outfilename = outDir+'/SummaryFile/'+name+'_SimuS' +typeSimu+ str("{:03d}".format(j)) + '.txt'
        f = open(outfilename, 'w')
        f.write('RunList' + ' ' + 'MergerID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Luminosity' + ' ' + 'TotalObservations' + ' ' + 'Obs' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'WindowinInputList' + ' ' + 'TotalProb'+ ' ' +'ObsInfo' + '\n')
        f.write(str(InputList['run'][j])+ ' ' +InputList['MergerID'][j].split('r')[-1]+ ' ' + str(InputList['Distance'][j])+ ' '+str(InputList['theta'][j])+ ' ' +str(InputList['A90'][j])+ ' ' +str(luminosity) + ' ' + str(len(InputObservationList)) + ' ' + str(allPossiblePoint) + ' ' + str(foundIn) +' ' + str(0) +' '+str(-1)+' '+str(totalProb)+' '+'True'+'\n')
    else:
        if type(foundIn) != int:
            foundFirst = foundIn[0]
            foundTimes=len(foundIn) # Has it been observed several times?
        else:
            foundFirst=foundIn
            foundTimes=1

        duration=InputObservationList['Duration']

        outfilename = outDir+'/SummaryFile/'+name+'_SimuSF' +typeSimu+ str("{:03d}".format(j)) + '.txt'
        # print(outfilename)
        f = open(outfilename, 'w')
        f.write('RunList' + ' ' + 'MergerID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Lum' + ' ' + 'Duration'+ '  '+'TotObs' + ' ' + 'TotalPossible' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'TotalProb'+ ' ' + 'ObsInfo'+'\n')
        f.write(str(InputList['run'][j])+ ' ' +InputList['MergerID'][j].split('r')[-1]+ ' ' + str(InputList['Distance'][j])+ ' '+str(InputList['theta'][j])+ ' ' +str(InputList['A90'][j])+ ' '  + str(luminosity) + ' ' + str(duration[foundFirst]) + ' '+ str(len(InputObservationList)) + ' ' + str(allPossiblePoint) + ' ' + str(foundFirst) + ' ' + str(foundTimes) + ' ' +str(totalProb)+ ' '+ 'True' +'\n')


def ReadSummaryFile(summaryFile):
    print('Note that this function needs to be adapted to the output')
    nP, Found = np.genfromtxt(summaryFile,usecols=(8,9),skip_header=1, unpack=True,dtype='str')
    return nP, Found


class NextWindowTools:

    @classmethod
    def CheckWindowCreateArray(cls,time, obsSite,WindowDurations):
        FullWindow = datetime.timedelta(seconds=np.float64(WindowDurations[-1]))
        #print('time',time)
        #print('FullWindow',FullWindow)
        #print(len(WindowDurations))
        if (Tools.IsDarkness(time, obsSite) is True) and (Tools.IsDarkness(time + FullWindow, obsSite) is True):
            LastItem=len(WindowDurations)
        else:
            print('Window is smaller')
            for i in range(0,len(WindowDurations)):
                if Tools.IsDarkness(time + datetime.timedelta(minutes=np.float64(WindowDurations[-i])), obsSite):
                    #print('Found!')
                    LastItem=len(WindowDurations)-i
        cumsumWindow=np.cumsum(WindowDurations)
        #print(cumsumWindow)
        #print(datetime.timedelta(seconds=np.float64(cumsumWindow[j])) for j in range(LastItem))
        arr = np.array([time + datetime.timedelta(seconds=np.float64(cumsumWindow[j])) for j in range(LastItem)])
        return arr
    @classmethod
    def NextObservationWindow(cls,time, obsSite):
        if(Tools.NextSunset(time, obsSite).hour >= time.hour >= Tools.PreviousSunrise(time,obsSite).hour and time.day == Tools.NextSunset(time, obsSite).day):
            time = Tools.NextSunset(time, obsSite)
            # print('Sunset', time)
            time = Tools.TrustingDarknessSun(time, obsSite)
            # print('Trusted', time)
        if(Tools.IsDarkness(time, obsSite) is True):
            return time
        elif ((Tools.IsDarkness(time, obsSite) is False) and (Tools.IsDarkness(Tools.NextMoonset(time, obsSite), obsSite)is True)):
            time=Tools.NextMoonset(time, obsSite)
            return time
        else:
            print('No window is found')
            return False

    @classmethod
    def EndObservationWindow(cls,time, obsSite):
        time = Tools.NextSunrise(time,obsSite)
        # Check if the night ends before due to the moon.
        if(Tools.IsDarkness(time, obsSite) is False):
            time = Tools.PreviousMoonset(time,obsSite)
        return time
class GRB(object):
    """
        Class to store GRB properties.

        Simulations are done with appropriate functions
        """

    def __init__(self, filepath=None,
                 name=None,
                 z=None,
                 time_interval=None,
                 spectral_model=None,
                 energy_interval=None):
        # Path of the GRB properties/model
        self.filepath = filepath

        # GRB properties
        self.name = name
        self.z = z

        # Time intervals
        self.time_interval = time_interval
        # Gammapy models
        self.spectral_model = spectral_model
        self.energy_interval = energy_interval

    def __str__(self):
        txt = ''
        txt += 'GRB summary\n'.format()
        txt += 'Name: {}\n'.format(self.name)
        txt += 'Redshift: {}\n'.format(self.z)
        txt += 'Times:\n'.format(self.time_interval)
        #txt += 'Energy steps:\n'.format(self.energy_interval)
        #txt += 'Flux:\n'.format(self.spectral_model)
        #for t in self.time_interval:
        #    txt += '{} -- {}\n'.format(t[0], t[1])

        #if self.stack_obs is not None:
        #    txt += str(self.stack_obs.total_stats_safe_range)

        return txt
    @classmethod
    def from_fitsfile(cls,filepath,absorption):
        data = fits.open(filepath)
        energy_interval=data[1].data.field(0)*u.GeV
        time_interval=data[2].data.field(0)* u.s
        flux=data[3].data
        z = 0.0 # No EBL correction for the moment, it is already 0.01 corrected
        name = filepath.split('/')[-1].split('.')[0]
        spectral_model = []
        for interval in range(0,len(time_interval)):
            flux_t = flux[interval]* u.Unit('1 / (cm2 s GeV)')
            newflux = flux_t #Use Factor 1000000.0  if needed
            #print(np.log(energy_interval.value))
            #print(np.log(flux_t.value))
            table_model = TableModel(energy=energy_interval,values=newflux,values_scale='log')
            table_model.plot(energy_range=(0.001, 10) * u.TeV)
            #name='/Users/mseglar/Documents/GitHub/CTASimulationsGW/run0017_MergerID000132_skymap/PGalonFoVanalysis_3d/plot'+str(interval)
            #plt.show()
            #plt.savefig('SpectralModel.png')
            #spectral_model.append(SpectralModel(spectral_model=table_model))
            spectral_model.append(AbsorbedSpectralModel(spectral_model=table_model,absorption=absorption,parameter=z,parameter_name='redshift'))
            #spectral_model.append(PowerLaw(index=2.2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV"))

        return cls(
            filepath=filepath,
            name=name,
            z=z,
            time_interval=time_interval,
            spectral_model=spectral_model,
            energy_interval=energy_interval,
        )

def ZenithAngleCut_TwoTimes(prob, nside, time, time1, MinProbCut,max_zenith,observatory):
    '''
    Mask in the pixels with zenith angle larger than max_zenith
    '''

    # Initial time
    frame = co.AltAz(obstime=time, location=observatory)
    pprob = prob

    mzenith = hp.ma(pprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=np.bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    altaz_map = targetCoord_map.transform_to(frame)
    maskzenith[altaz_map.alt.value < 90-max_zenith] = 1
    mzenith.mask = maskzenith
    #hp.mollview(mzenith)
    #plt.show()
    #plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_%g.png")

    yprob = ma.masked_array(pprob, mzenith.mask)
    #hp.mollview(yprob)
    #plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_prob_%g.png")

    #print('Integrated probability of the masked map', np.sum(yprob))

    # End time
    frame = co.AltAz(obstime=time1, location=observatory)
    ppprob = pprob

    mzenith = hp.ma(ppprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=np.bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    altaz_map = targetCoord_map.transform_to(frame)
    maskzenith[altaz_map.alt.value < 90-max_zenith] = 1
    mzenith.mask = maskzenith

    #plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_%g.png")

    xprob = ma.masked_array(yprob, mzenith.mask)

    if np.sum(yprob) < MinProbCut:
        ObsBool = False
    else:
        ObsBool = True


    return ObsBool, yprob


def ComputeProbability2D_SelectClusters(prob, highres, radecs, ReducedNside, HRnside, MinProbCut, TotalExposure, time,DelayObs,interObsSlew, observatory, max_zenith, FOV,
                         run, mergerID, ipixlist, ipixlistHR, counter, datasetDir,outDir, usegreytime, plot):
    '''
    Compute probability in 2D by taking the highest value pixel
    '''
    radius = FOV


    frame = co.AltAz(obstime=time, location=observatory.Location)
    thisaltaz = radecs.transform_to(frame)
    pix_alt1 = thisaltaz.alt.value

    if usegreytime:
        moonaltazs = get_moon(Time(time)).transform_to(AltAz(obstime=Time(time), location=observatory.Location))
        # Zenith and Moon angular distance mask
        pix_ra = radecs.ra.value[
            (thisaltaz.alt.value > 90 - max_zenith) & (thisaltaz.separation(moonaltazs) > 30 * u.deg)]
        pix_dec = radecs.dec.value[
            (thisaltaz.alt.value > 90 - max_zenith) & (thisaltaz.separation(moonaltazs) > 30 * u.deg)]

    else:
        # Zenith angle mask
        pix_ra = radecs.ra.value[thisaltaz.alt.value > (90 - max_zenith)]
        pix_dec = radecs.dec.value[thisaltaz.alt.value > (90 - max_zenith)]
        zenith_ini = 90 - pix_alt1[thisaltaz.alt.value > (90 - max_zenith)]

    phipix = np.deg2rad(pix_ra)
    thetapix = 0.5 * np.pi - np.deg2rad(pix_dec)

    ipix = hp.ang2pix(ReducedNside, thetapix, phipix)

    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)
    exposure = np.empty(len(pix_ra), dtype=object)
    zenith_end = np.empty(len(pix_ra), dtype=object)
    #print(len(phipix),len(thetapix),len(ipix),len(pix_ra),len(pix_dec), len(dp_Pix_Fov), len(pix_alt), len(zenith_end), len(exposure))

    cat_pix = Table([ipix, pix_ra, pix_dec, dp_Pix_Fov, zenith_ini, zenith_end, exposure], names=('PIX', 'PIXRA', 'PIXDEC', 'PIXFOVPROB','ZENITH_INI','ZENITH_END','EXPOSURE'))

    dp_dV_FOV = []

    xyzpix = hp.ang2vec(thetapix, phipix)

    for i in range(0, len(cat_pix)):
        ipix_discfull = hp.query_disc(HRnside, xyzpix[i], np.deg2rad(radius))
        maskComputeProb = [np.isin(ipix_discfull, ipixlistHR, invert=True)]
        dp_dV_FOV.append(highres[ipix_discfull[maskComputeProb]].sum())

    cat_pix['PIXFOVPROB'] = dp_dV_FOV

    # Mask already observed pixels

    mask = [np.isin(cat_pix['PIX'], ipixlist, invert=True)]

    if all(np.isin(cat_pix['PIX'], ipixlist, invert=False)):
        maskcat_pix = cat_pix
    else:
        maskcat_pix = cat_pix[mask]

    # Sort table
    sortcat = maskcat_pix[np.flipud(np.argsort(maskcat_pix['PIXFOVPROB']))]
    ObsCase = 'SourceOutFoV' # Default case, the source is not in CTA FoV

    # Fill column for time that one needs to observe them to get 5sigma for the highest of the list

    if(np.any(sortcat['ZENITH_INI']>55)):
        ObsCase, texp60 = ObtainSingleObservingTimes(TotalExposure, DelayObs, interObsSlew, run, mergerID, observatory,datasetDir, zenith=60)
        print("Lets have a look on what was found for Zenith= 60 ", texp60)
        # Cat60 = sortcat[sortcat['ZENITH_INI'] >55]
        # print("ObsCase60", ObsCase)
        sortcat['EXPOSURE'][sortcat['ZENITH_INI'] > 55] = texp60
        frame = co.AltAz(obstime= time + datetime.timedelta(seconds=texp60), location=observatory.Location)
        catCoord60 = co.SkyCoord(sortcat['PIXRA'][sortcat['ZENITH_INI'] > 55], sortcat['PIXDEC'][sortcat['ZENITH_INI'] > 55], frame='fk5', unit=(u.deg, u.deg))
        thisaltaz60 = catCoord60.transform_to(frame)
        pix_zen60 = 90 - thisaltaz60.alt.value
        sortcat['ZENITH_END'][sortcat['ZENITH_INI'] > 55] = pix_zen60

    mask1 = sortcat['ZENITH_INI'] >= 30
    mask2 = sortcat['ZENITH_INI'] <= 55

    if (sortcat['ZENITH_INI'][mask1&mask2].any()):
        ObsCase, texp40 = ObtainSingleObservingTimes(TotalExposure, DelayObs,interObsSlew, run, mergerID, observatory, datasetDir, zenith=40)
        print("ObsCase40", ObsCase)
        sortcat['EXPOSURE'][(30 <= sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] <= 55)] = texp40
        # Cat40 = sortcat[(30 < sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] < 55)]
        # print(Cat40)
        frame = co.AltAz(obstime=time + datetime.timedelta(seconds=texp40), location=observatory.Location)
        catCoord40 = co.SkyCoord(sortcat['PIXRA'][(30 <= sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] <= 55)], sortcat['PIXDEC'][(30 <= sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] <= 55)], frame='fk5', unit=(u.deg, u.deg))
        # print("radecs",radecs)
        thisaltaz40 = catCoord40.transform_to(frame)
        pix_zen40 = 90 - thisaltaz40.alt.value
        sortcat["ZENITH_END"][(30 <= sortcat['ZENITH_INI']) & (sortcat['ZENITH_INI'] <= 55)] = pix_zen40


    if(np.any(sortcat['ZENITH_INI']<30)):
        ObsCase, texp20 = ObtainSingleObservingTimes(TotalExposure, DelayObs, interObsSlew, run, mergerID, observatory, datasetDir, zenith=20)
        print("ObsCase20", ObsCase)
        sortcat['EXPOSURE'][sortcat['ZENITH_INI'] < 30] = texp20
        # Cat20 = sortcat[sortcat['ZENITH_INI'] < 30]
        # print(Cat30)
        # print("time + datetime.timedelta(seconds=texp20)",time + datetime.timedelta(seconds=texp20))
        # print("observatory.Location",observatory.Location)
        frame = co.AltAz(obstime=time + datetime.timedelta(seconds=texp20), location=observatory.Location)
        catCoord20 = co.SkyCoord(sortcat['PIXRA'][sortcat['ZENITH_INI'] < 30], sortcat['PIXDEC'][sortcat['ZENITH_INI'] < 30], frame='fk5', unit=(u.deg, u.deg))
        thisaltaz20 = catCoord20.transform_to(frame)
        pix_zen20 = 90 - thisaltaz20.alt.value
        sortcat["ZENITH_END"][sortcat['ZENITH_INI'] < 30] = pix_zen20


    # Mask if texp is False due to the fact that no observation window is achieved.
    # print('HERE', len(sortcat))
    # print(sortcat['EXPOSURE'])
    if False in sortcat['EXPOSURE']:
        maskExposure = [sortcat['EXPOSURE'] != False]
        sortcat = sortcat[maskExposure]

    # Mask if it is not visible at the end of the window
    mask = [np.isin(sortcat['ZENITH_END'], 66, invert=True)]
    maskcat_zen = sortcat[mask]
    if(len(sortcat[mask]) == 0):
        P_GW = 0
        targetCoord = False
        ObsExp = False
        ZenIni = False
        ZenEnd = False
    else:
        sortcat = maskcat_zen[np.flipud(np.argsort(maskcat_zen['PIXFOVPROB']))]
        targetCoord = co.SkyCoord(sortcat['PIXRA'][:1][0], sortcat['PIXDEC'][:1][0], frame='fk5', unit=(u.deg, u.deg))
        # print('targetCoord',targetCoord)

        P_GW = sortcat['PIXFOVPROB'][:1][0]
        ObsExp = sortcat['EXPOSURE'][:1][0]
        ZenIni = sortcat['ZENITH_INI'][:1][0]
        ZenEnd = sortcat['ZENITH_END'][:1][0]


        # Include to the list of pixels already observed

        if (P_GW >= MinProbCut):
            # Chose highest
            phip = float(np.deg2rad(targetCoord.ra.deg))
            thetap = float(0.5 * np.pi - np.deg2rad(targetCoord.dec.deg))
            xyz = hp.ang2vec(thetap, phip)

            # ipix_discComputeProb = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
            # maskComputeProb=[np.isin(ipix_discComputeProb, ipixlist,invert=True)]
            # print('maskedP_GW',highres[ipix_discComputeProb[maskComputeProb]].sum())
            ipixlistHR.extend(hp.query_disc(HRnside, xyz, np.deg2rad(radius)))
            ipix_disc = hp.query_disc(ReducedNside, xyz, np.deg2rad(radius))
            ipixlist.extend(ipix_disc)

            ######################################

            # PLOT THE RESULTS
            if (plot):
                #path = outDir + '/EvolutionPlot'
                path = '%s/Pointing_Plotting/%s_%s/EvolutionPlot' % (outDir,run,mergerID)
                if not os.path.exists(path):
                    os.makedirs(path)
                # nside = 1024

                hp.mollview(prob,title="With FoV circle")

                #hp.gnomview(prob, xsize=500, ysize=500, rot=[targetCoord.ra.deg, targetCoord.dec.deg], reso=8.0)
                hp.graticule()
                # print('This skymap has nside equals to',hp.npix2nside(len(highres)))
                # plt.savefig('%s/Pointing-prob_%g.png' % (path, counter))

                ipix_discplot = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
                tt, pp = hp.pix2ang(HRnside, ipix_discplot)
                ra2 = np.rad2deg(pp)
                dec2 = np.rad2deg(0.5 * np.pi - tt)
                skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
                # hp.visufunc.projplot(skycoord.ra, skycoord.dec, 'y.', lonlat=True, coord="C")
                # plt.show()
                # observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

                hp.visufunc.projplot(sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], 'r.', lonlat=True, coord="C")
                MaxCoord = SkyCoord(sortcat['PIXRA'][:1], sortcat['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))
                separations = skycoord.separation(MaxCoord)
                tempmask = separations < (radius + 0.01 * radius) * u.deg
                tempmask2 = separations > (radius - 0.01 * radius) * u.deg
                hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'g.', lonlat=True,
                             coord="C", linewidth=0.1)
                plt.savefig('%s/Pointing-prob-FoV_%g.png' % (path, counter))
                altcoord = np.empty(100)
                azcoord = np.random.rand(100) * 360
                plt.savefig('%s/Pointing-prob-FoV_%g.png' % (path, counter))
		# for i in range(0,1):
                #    altcoord.fill(90-(max_zenith-5*i))
                #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
                #    RandomCoord_radec = RandomCoord.transform_to('fk5')
                #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
                # plt.show()
                # plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))
        altcoord = np.empty(100)
        azcoord = np.random.rand(100) * 360
        # for i in range(0,1):
        #    altcoord.fill(90-(max_zenith-5*i))
        #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
        #    RandomCoord_radec = RandomCoord.transform_to('fk5')
        #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # plt.show()
        # plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))

    return P_GW, targetCoord,ObsExp, ZenIni, ZenEnd, ObsCase, ipixlist, ipixlistHR
