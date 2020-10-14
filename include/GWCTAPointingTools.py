import sys
sys.path.append('../../GW-Followup/include')
from GWHESSPointingTools import *
from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel
#from gammapy.spectrum.models import PowerLaw

iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))



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

def LoadGalaxiesSimulation(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    dist, z, ra, dec = np.genfromtxt(tgalFile, usecols=(1, 2, 3,4), skip_header=3, unpack=True)  # ra, dec in degrees

    tcat = Table([ra, dec, dist], names=('RAJ2000', 'DEJ2000', 'Dist'))
    return tcat

def TableImportCTA(tgalFile):
    run,MergerID,RA,Dec,distance,distMin,distMax,z,theta,ndet,SNR,A90,A50 = np.genfromtxt(tgalFile,usecols=(0, 1, 2, 3,4,5,6,7,8,9,10,11,12),skip_header=1, unpack=True,dtype='str')
    RA=RA.astype(np.float)
    Dec=Dec.astype(np.float)
    z=z.astype(np.float)
    distance = distance.astype(np.float)
    distMax = distMax.astype(np.float)
    distMin= distMin.astype(np.float)
    theta=theta.astype(np.float)
    ndet=ndet.astype(np.float)
    SNR=SNR.astype(np.float)
    A90=A90.astype(np.float)
    A50=A50.astype(np.float)
    OutputTable= Table([run,MergerID,RA,Dec,distance,distMin,distMax,z,theta,ndet, SNR,A90,A50],names=('run','MergerID','RA','Dec','Distance','DistMin','DistMax','redshift','theta','ndet','SNR','A90','A50'))
    #print(OutputTable)
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
    distance = distance.astype(np.float)
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
            print(xyz)
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

def ProduceSummaryFile(InputList,InputObservationList,allPossiblePoint,foundIn,j,typeSimu,totalProb,dirName,name):
    #print(InputList['run'][j])
    #print(InputList['MergerID'][j].split('r')[-1])
    filepath='../dataset/GammaCatalogV1.0/'+str(InputList['run'][j]) + '_' + str(InputList['MergerID'][j].split('r')[-1]) + ".fits"
    fitsfile = fits.open(filepath)
    luminosity = fitsfile[0].header['EISO']
    dirNameFile = dirName+'/SummaryFile'
    if not os.path.exists(dirNameFile):
        os.makedirs(dirNameFile)
    # Obtain the luminosity
    if foundIn ==-1:
        outfilename = dirNameFile+'/'+name+'_Simu' +typeSimu+ str("{:03d}".format(j)) + '.txt'
        f = open(outfilename, 'w')
        f.write('RunList' + ' ' + 'MergerID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Luminosity' + ' ' + 'TotalObservations' + ' ' + 'Obs' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'WindowinInputList' + ' ' + 'TotalProb' + '\n')
        f.write(str(InputList['run'][j])+ ' ' +InputList['MergerID'][j].split('r')[-1]+ ' ' + str(InputList['Distance'][j])+ ' '+str(InputList['theta'][j])+ ' ' +str(InputList['A90'][j])+ ' ' +str(luminosity) + ' ' + str(len(InputObservationList)) + ' ' + str(allPossiblePoint) + ' ' + str(foundIn) +' ' + str(0) +' '+str(-1)+' '+str(totalProb)+'\n')
    else:
        if type(foundIn) != int:
            foundFirst = foundIn[0]
            foundTimes=len(foundIn) # Has it been observed several times?
        else:
            foundFirst=foundIn
            foundTimes=1

        duration=InputObservationList['Duration']

        outfilename = dirNameFile+'/'+name+'_Simu' +typeSimu+ str("{:03d}".format(j)) + '.txt'
        # print(outfilename)
        f = open(outfilename, 'w')
        f.write('RunList' + ' ' + 'MergerID' + ' ' + 'Distance' + ' ' + 'Theta' + ' ' + 'A90' + ' ' + 'Lum' + ' ' + 'Duration'+ '  '+'TotObs' + ' ' + 'TotalPossible' + ' ' + 'FirstCovered' + ' ' + 'TimesFound' + ' ' + 'TotalProb' + '\n')
        f.write(str(InputList['run'][j])+ ' ' +InputList['MergerID'][j].split('r')[-1]+ ' ' + str(InputList['Distance'][j])+ ' '+str(InputList['theta'][j])+ ' ' +str(InputList['A90'][j])+ ' '  + str(luminosity) + ' ' + str(duration[foundFirst]) + ' '+ str(len(InputObservationList)) + ' ' + str(allPossiblePoint) + ' ' + str(foundFirst) + ' ' + str(foundTimes) + ' ' +str(totalProb)+'\n')

class NextWindowTools:

    @classmethod
    def CheckWindowCreateArray(cls,time, obsSite,WindowDurations):
        FullWindow = datetime.timedelta(seconds=np.float64(WindowDurations[-1]))
        #print('time',time)
        #print('FullWindow',FullWindow)
        print(len(WindowDurations))
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
        #print('ARRRRRRRRR',arr)
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
        z=0.0 # No EBL correction for the moment, it is already 0.01 corrected
        name=filepath.split('/')[-1].split('.')[0]
        spectral_model=[]
        for interval in range(0,len(time_interval)):
            flux_t=flux[interval]* u.Unit('1 / (cm2 s GeV)')
            newflux=flux_t #Use Factor 1000000.0  if needed
            #print(np.log(energy_interval.value))
            #print(np.log(flux_t.value))
            table_model = TableModel(energy=energy_interval,values=newflux,norm=1.,values_scale='log')
            table_model.plot(energy_range=(0.001, 10) * u.TeV)
            #name='/Users/mseglar/Documents/GitHub/CTASimulationsGW/run0017_MergerID000132_skymap/PGalonFoVanalysis_3d/plot'+str(interval)
            #plt.show()
            #plt.savefig('Muere.png')
            #spectral_model.append(SpectralModel(spectral_model=table_model))
            spectral_model.append(AbsorbedSpectralModel(spectral_model=table_model,absorption=absorption,parameter=z,parameter_name='redshift'))
            #spectral_model.append(PowerLaw(index=2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV"))

        return cls(
            filepath=filepath,
            name=name,
            z=z,
            time_interval=time_interval,
            spectral_model=spectral_model,
            energy_interval=energy_interval,
        )
