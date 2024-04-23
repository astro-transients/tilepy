import numpy as np
import astropy.units as u
from gammapy.extern.pathlib import Path
from gammapy.maps import WcsGeom, MapAxis, Map, WcsNDMap
from gammapy.cube import  MapEvaluator, PSFKernel
from gammapy.cube.models import SkyModel
from gammapy.image.models import SkyGaussian
from gammapy.irf import load_cta_irfs
from astropy.coordinates import SkyCoord, Angle
from gammapy.cube import make_map_exposure_true_energy, make_map_background_irf
from gammapy.spectrum.models import Absorption, PowerLaw
import matplotlib.pyplot as plt
from .PointingTools import GRB
from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel, PowerLaw


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
        return txt

    @classmethod
    def from_fitsfile(cls, filepath, absorption):
        data = fits.open(filepath)
        energy_interval = data[1].data.field(0) * u.GeV
        time_interval = data[2].data.field(0) * u.s
        flux = data[3].data
        z = 0.0  # No EBL correction for the moment, needs to be added
        name = filepath.split('/')[-1].split('.')[0]
        spectral_model = []
        for interval in range(0, len(time_interval)):
            flux_t = flux[interval] * u.Unit('1 / (cm2 s GeV)')
            newflux = flux_t  # Use Factor 1000000.0  if needed
            # print(np.log(energy_interval.value))
            # print(np.log(flux_t.value))
            table_model = TableModel(energy=energy_interval, values=newflux, values_scale='log')
            table_model.plot(energy_range=(0.001, 10) * u.TeV)
            # name='/Users/mseglar/Documents/GitHub/CTASimulationsGW/run0017_MergerID000132_skymap/PGalonFoVanalysis_3d/plot'+str(interval)
            # plt.show()
            # plt.savefig('SpectralModel.png')
            # spectral_model.append(SpectralModel(spectral_model=table_model))
            spectral_model.append(AbsorbedSpectralModel(spectral_model=table_model, absorption=absorption, parameter=z,
                                                        parameter_name='redshift'))
            # spectral_model.append(PowerLaw(index=2.2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV"))

        return cls(
            filepath=filepath,
            name=name,
            z=z,
            time_interval=time_interval,
            spectral_model=spectral_model,
            energy_interval=energy_interval,
        )


def PointingSimulation(sky_model,axis,geom,pointing,livetime,offset,irfs):
    exposure = make_map_exposure_true_energy(pointing=pointing, livetime=livetime, aeff=irfs["aeff"], geom=geom)
    background = make_map_background_irf(pointing=pointing, ontime=livetime, bkg=irfs["bkg"], geom=geom)
    psf = irfs["psf"].to_energy_dependent_table_psf(theta=offset)
    psf_kernel = PSFKernel.from_table_psf(psf, geom, max_radius=0.2 * u.deg)
    energy = axis.edges * axis.unit
    edisp = irfs["edisp"].to_energy_dispersion(offset, e_reco=energy, e_true=energy)

    #Evaluate and compute things
    evaluator = MapEvaluator(model=sky_model,exposure=exposure,background=background,psf=psf_kernel,edisp=edisp,)
    npred = evaluator.compute_npred()
    rng = np.random.RandomState(seed=None)
    npred = npred.clip(min=0)
    counts = rng.poisson(npred)
    return npred,counts, background,exposure

def DynamicPointingSimulation(slew,step,TimeAxis,efLivetime,spatial_model,grb,axis,geom,geom3D,pointing,LCintervals,offset,irfs,Cube4DMap,Background4DMap,Exposure4DMap):
# ToDo: need to update this function
    LC_intTimes=[]
    boolmask=False

    for i in range(len(TimeAxis)):
        boolmask=step<TimeAxis[i]
        if boolmask:
            #print('First interval found is',TimeAxis[i],'from',AT)
            LC_intTimes.append(TimeAxis[i]-step)
            break

    EndWindow=efLivetime+slew
    boolmask=False
    for j in range(i,len(TimeAxis)):
        boolmask=EndWindow<TimeAxis[j+1]
        if boolmask:
            #print('End of this window happens at',EndWindow,'with the last interval being',TimeAxis[j])
            LC_intTimes.append(EndWindow-TimeAxis[j])
            break
        else:
            LC_intTimes.append(TimeAxis[j+1]-TimeAxis[j])

    for k in range(0,len(LC_intTimes)):
        livetime1=LC_intTimes[k]* u.second
        mask=[i for i in LCintervals.value if i <=(step+LC_intTimes[k])]
        sky_model=SkyModel(spatial_model=spatial_model, spectral_model=grb.spectral_model[len(mask)-1])
        npred,c, b,e=PointingSimulation(sky_model,axis,geom,pointing,livetime1,offset,irfs)
        print('index of the data',i+k)
        # toDo: THE INDEX IS NOT GOOOOOD!!!!!!!
        # toDo: when availableMapDataset
        Cube4DMap.data[i+k] = WcsNDMap(geom,c).reproject(geom3D).data
        Background4DMap.data[i+k] = b.reproject(geom3D).data
        Exposure4DMap.data[i+k] = e.reproject(geom3D).data

    return Cube4DMap,Background4DMap,Exposure4DMap

def CubeSimulationMultiObservations(InputListrun,InputListMergerID,pointing,nP,Source,IRFfilename,InputObservationTime, datasetDir, outDir):
    # ToDo: need to update this function
    absorption = Absorption.read_builtin('dominguez')
    filepath = datasetDir +'/GammaCatalogV2.0/'+InputListrun + '_' + InputListMergerID.split('r')[-1] + ".fits"

    #filepath='/Users/mseglar/Documents/GitHub/CTASimulationsGW/GRBtemplates/SGRB001.fits'
    grb = GRB.from_fitsfile(filepath=filepath, absorption=absorption)
    #filenameP='ExamplePointing.txt'
    #RAP=Pointings['RA(deg)']
    #DecP=Pointings['DEC(deg)']

    RAP=pointing.ra.deg
    DecP=pointing.dec.deg

    #filenameS='ExampleSource.txt'
    #RAS,DecS=np.genfromtxt(filenameS,skip_header=1)

    spatial_model = SkyGaussian(lon_0=Source.ra.deg*u.deg, lat_0=Source.dec.deg*u.deg, sigma="0.1 deg")

    print()
    print('Source location is,',Source.ra.deg*u.deg,Source.dec.deg*u.deg)
    print('Pointing coordinates is,',RAP[nP], DecP[nP])
    print()

    for i in range(0,len(grb.spectral_model)):
        grb.spectral_model[i].plot(energy_range=(0.001, 10) * u.TeV)
    #plt.grid()
    #plt.savefig('./CubeSimulation_Figures/Fits_Spectra.png')

    # Get Spectral model from fits!
    #spectral_model = PowerLaw(index=2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV")
    #plt.savefig('./CubeSimulation_Figures/LGRB_Evolution.png')

    #IRFfilename = ("/Users/mseglar/Documents/CurrentPhD/CTA/grb_paper-master/IRFs_prod3b/bcf/North_NSBx05_z20_N_LST_30m/irf_file.fits")
    irfs = load_cta_irfs(IRFfilename)

    ###########################################
    #Cube Definition
    axis = MapAxis.from_edges(np.logspace(-2.0, 1.0, 20), unit="TeV", name="energy", interp="log")

    ###########################################
    # TimeRun => Pointing time relative to MERGER!
    # This matches the format that is it given in the LC evolution, so a merge is posible.
    livetimes=[InputObservationTime['tI'],InputObservationTime['tF']]


    #livetimes=np.array([0,20,59,101,250,595,1905,3505,4905,6505,9905])#12000,15000
    print(livetimes)
    TimeAxis=np.unique(np.concatenate((livetimes,grb.time_interval.value)))

    #offset_max = 2.5 * u.deg
    offset = Angle("2.5 deg")

    energy_axis = MapAxis.from_edges(np.logspace(-2.0, 10, 10), unit="TeV", name="energy", interp="log")
    time_axis = MapAxis.from_edges(TimeAxis, unit='second', name="time", interp="log")
    geom3D = WcsGeom.create(skydir=(RAP[0],DecP[0]), binsz=0.06, width=(30, 24), coordsys="CEL", axes=[energy_axis])
    geom4D = WcsGeom.create(skydir=(RAP[0],DecP[0]), binsz=0.06, width=(30, 24), coordsys="CEL", axes=[energy_axis, time_axis])
    Cube4DMap = Map.from_geom(geom4D)
    Background4DMap = Map.from_geom(geom4D)
    Exposure4DMap = Map.from_geom(geom4D)

    timestep=0
    slew=0
    counterPointings=0
    for i in range(0,len(livetimes)):
        efLivetimes=livetimes

        if efLivetimes[i]>0:
            geom = WcsGeom.create(skydir=(RAP[counterPointings], DecP[counterPointings]), binsz=0.04, width=(10, 10), coordsys="CEL", axes=[axis])
            pointing = SkyCoord(RAP[counterPointings], DecP[counterPointings], unit="deg", frame="fk5")
            Cube4DMap, Background4DMap, Exposure4DMap=DynamicPointingSimulation(slew,slew+timestep,TimeAxis,efLivetimes[i],spatial_model,grb,axis,geom,geom3D,pointing,grb.time_interval,offset,irfs,Cube4DMap,Background4DMap,Exposure4DMap)
            timestep=efLivetimes[i]
            counterPointings=counterPointings+1
    print('========= RESULTS ========')
    print(Cube4DMap.data.shape)
    print(Background4DMap.data.shape)
    print(Exposure4DMap.data.shape)

    maps = {
        "counts": Cube4DMap,
        "background": Background4DMap,
        "exposure": Exposure4DMap,
    }

    print(maps)

    path = "analysis_3d"
    print(path)
    print(outDir)
    fullpath=Path(outDir+path)

    fullpath.mkdir(exist_ok=True)

    # write maps
    maps["counts"].write(str(path / "counts_multiObs.fits"), overwrite=True)
    maps["background"].write(str(path / "background_multiObs.fits"), overwrite=True)
    maps["exposure"].write(str(path / "exposure_multiObs.fits"), overwrite=True)

def CubeSimulationMultiObservations_dynamicSpectrum(InputListrun,InputListMergerID,pointing,nP,Source,IRFfilename,InputObservationTime,datasetDir,outDir):
    # ToDo: need to update this function
    absorption = Absorption.read_builtin('dominguez')
    filepath= datasetDir +'/GammaCatalogV2.0/'+InputListrun + '_' + InputListMergerID.split('r')[-1] + ".fits"

    #filepath='/Users/mseglar/Documents/GitHub/CTASimulationsGW/GRBtemplates/SGRB001.fits'
    grb = GRB.from_fitsfile(filepath=filepath, absorption=absorption)

    RAP=pointing.ra.deg
    DecP=pointing.dec.deg

    #filenameS='ExampleSource.txt'
    #RAS,DecS=np.genfromtxt(filenameS,skip_header=1)

    spatial_model = SkyGaussian(lon_0=Source.ra.deg*u.deg, lat_0=Source.dec.deg*u.deg, sigma="0.1 deg")
    print()
    print('Source location is,',Source.ra.deg*u.deg,Source.dec.deg*u.deg)
    print('Pointing coordinates is,',RAP[nP], DecP[nP])
    print()

    for i in range(0,len(grb.spectral_model)):
        #print('iteration',i)
        grb.spectral_model[i].plot(energy_range=(0.001, 10) * u.TeV)
    #plt.show()
    #plt.grid()
    #plt.savefig('./CubeSimulation_Figures/Fits_Spectra.png')

    # Get Spectral model from fits!
    #spectral_model = PowerLaw(index=2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV")
    #plt.savefig('./CubeSimulation_Figures/LGRB_Evolution.png')

    #Select the IRFs depending on observation condition and site
    irfs = load_cta_irfs(IRFfilename)

    ###########################################
    #Cube Definition
    #axis = MapAxis.from_edges(np.logspace(-2.0, 1.0, 20), unit="TeV", name="energy", interp="log")

    ###########################################
    # TimeRun => Pointing time relative to MERGER!
    # This matches the format that is it given in the LC evolution, so a merge is posible.
    #livetimes=InputObservationTime['Interval']
    #startingObsTimes=[InputObservationTime['tI'],InputObservationTime['tF']]

    #print('GRB time intervals',grb.time_interval.value)
    #print('startingObsTimes',startingObsTimes)
    #print('GRB simulation should take into account this:',InputObservationTime)
    #print('Livetime of the observation is:',livetimes[0])
    #TimeAxis=np.unique(np.concatenate((startingObsTimes,grb.time_interval.value)))
    #selection1 = TimeAxis>=startingObsTimes[0]
    #selection2 = TimeAxis<=startingObsTimes[1]

    #mask1 = grb.time_interval.value<=startingObsTimes[0]
    #mask2 = grb.time_interval.value>=startingObsTimes[1]
    #print(len(grb.time_interval.value))
    #print(len(grb.time_interval.value[mask1]))
    #print(len(grb.time_interval.value[mask2])+len(startingObsTimes))
    #print(TimeAxis[selection1&selection2])
    #SimulationInterval=TimeAxis[selection1&selection2]
    offset_max = 2.5 * u.deg
    offset = Angle("2.5 deg")

    energy_axis = MapAxis.from_edges(np.logspace(-2.0, 1.0, 20), unit="TeV", name="energy", interp="log")
    time_axis = MapAxis.from_edges(np.arange(len(pointing)+1), unit='second', name="time")
    print(time_axis)
    print(energy_axis)
    geom3D = WcsGeom.create(skydir=(RAP[3],DecP[3]), binsz=0.06, width=(52, 48), coordsys="CEL", axes=[energy_axis])
    Cube3DMap = Map.from_geom(geom3D)
    Background3DMap = Map.from_geom(geom3D)
    Exposure3DMap = Map.from_geom(geom3D)

    geom4D = WcsGeom.create(skydir=(RAP[3],DecP[3]), binsz=0.06, width=(52, 48), coordsys="CEL", axes=[energy_axis, time_axis])
    Cube4DMap = Map.from_geom(geom4D)
    Background4DMap = Map.from_geom(geom4D)
    Exposure4DMap = Map.from_geom(geom4D)

    #counterPointings=0
    #geom = WcsGeom.create(skydir=(RAP[nP], DecP[nP]), binsz=0.04, width=(10, 10), coordsys="CEL", axes=[axis])

    #FirstSpectralModel=len(grb.time_interval.value[mask1])
    for i in range(0,len(pointing)):
        pointing = SkyCoord(RAP[i], DecP[i], unit="deg", frame="fk5")
        geom = WcsGeom.create(skydir=(RAP[i], DecP[i]), binsz=0.04, width=(12, 12),coordsys="CEL", axes=[energy_axis])
        sky_model = SkyModel(spatial_model=spatial_model, spectral_model=grb.spectral_model[10])
        npred, counts, background, exposure= PointingSimulation(sky_model,energy_axis,geom,pointing,15*u.second,offset,irfs)

        Cube4DMap.data[i] = WcsNDMap(geom,counts).reproject(geom3D).data
        Background4DMap.data[i] = background.reproject(geom3D).data
        Exposure4DMap.data[i] = exposure.reproject(geom3D).data


    #Cube3DMap.data=counts
    #Background3DMap.data=background.data
    #Exposure3DMap.data=exposure.data

    print('========= RESULTS ========')
    print(Cube4DMap.data.shape)
    print(Background4DMap.data.shape)
    print(Exposure4DMap.data.shape)

    maps = {
        "counts": Cube4DMap,
        "background": Background4DMap,
        "exposure": Exposure4DMap,
    }

    print(maps)

    path = "analysis_3d"
    print(path)
    print(outDir)
    fullpath=Path(outDir+path)

    fullpath.mkdir(exist_ok=True)

    # write maps
    maps["counts"].write(str(fullpath / "counts_multipleObs.fits"), overwrite=True)
    maps["background"].write(str(fullpath / "background_multipleObs.fits"), overwrite=True)
    maps["exposure"].write(str(fullpath / "exposure_multipleObs.fits"), overwrite=True)


def CubeSimulationSingleObservation_dynamicSpectrum(InputListrun,InputListMergerID,nP,Source,IRFfilename,InputObservationTime,pointing,datasetDir,outDir):
    # ToDo: need to update this function
    absorption = Absorption.read_builtin('dominguez')
    filepath= datasetDir +'/GammaCatalogV2.0/'+InputListrun + '_' + InputListMergerID.split('r')[-1] + ".fits"

    #filepath='/Users/mseglar/Documents/GitHub/CTASimulationsGW/GRBtemplates/SGRB001.fits'
    grb = GRB.from_fitsfile(filepath=filepath, absorption=absorption)

    PointingsCoord = SkyCoord(Pointings['RA(deg)'], Pointings['DEC(deg)'], frame='fk5', unit=(u.deg, u.deg))

    RAP = PointingsCoord.ra.deg
    DecP = PointingsCoord.dec.deg

    #filenameS='ExampleSource.txt'
    #RAS,DecS=np.genfromtxt(filenameS,skip_header=1)

    spatial_model = SkyGaussian(lon_0=Source.ra.deg*u.deg, lat_0=Source.dec.deg*u.deg, sigma="0.1 deg")
    print()
    print('Source location is,',Source.ra.deg*u.deg,Source.dec.deg*u.deg)
    print('Pointing coordinates is,',RAP[nP], DecP[nP])
    print()

    for i in range(0,len(grb.spectral_model)):
        grb.spectral_model[i].plot(energy_range=(0.001, 10) * u.TeV)

    #plt.show()
    #plt.grid()
    #plt.savefig('./CubeSimulation_Figures/Fits_Spectra.png')

    # Get Spectral model from fits!
    #spectral_model = PowerLaw(index=2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV")
    #plt.savefig('./CubeSimulation_Figures/LGRB_Evolution.png')

    #Select the IRFs depending on observation condition and site
    irfs = load_cta_irfs(IRFfilename)

    ###########################################
    #Cube Definition
    axis = MapAxis.from_edges(np.logspace(-2.0, 1.0, 20), unit="TeV", name="energy", interp="log")

    ###########################################
    # TimeRun => Pointing time relative to MERGER!
    # This matches the format that is it given in the LC evolution, so a merge is posible.
    #livetimes=InputObservationTime['Interval']
    #startingObsTimes=[InputObservationTime['tI'],InputObservationTime['tF']]

    #print('GRB time intervals',grb.time_interval.value)
    print('startingObsTimes', startingObsTimes)
    #print('startingObsTimes',InputObservationTime['tI'].data)
    #print('GRB simulation should take into account this:',InputObservationTime)
    #print('Livetime of the observation is:',livetimes[0])
    TimeAxis=np.unique(np.concatenate((startingObsTimes,grb.time_interval.value)))
    selection1 = TimeAxis>=startingObsTimes[0]
    selection2 = TimeAxis<=startingObsTimes[1]

    mask1 = grb.time_interval.value<=startingObsTimes[0]
    mask2 = grb.time_interval.value>=startingObsTimes[1]
    #print(len(grb.time_interval.value))
    #print(len(grb.time_interval.value[mask1]))
    #print(len(grb.time_interval.value[mask2])+len(startingObsTimes))
    #print(TimeAxis[selection1&selection2])
    SimulationInterval=TimeAxis[selection1&selection2]
    offset_max = 2.5 * u.deg
    offset = Angle("2.5 deg")

    energy_axis = MapAxis.from_edges(np.logspace(-2.0, 1.0, 20), unit="TeV", name="energy", interp="log")
    time_axis = MapAxis.from_edges(SimulationInterval, unit='second', name="time", interp="log")
    geom3D = WcsGeom.create(skydir=(RAP[nP],DecP[nP]), binsz=0.04, width=(16, 14), coordsys="CEL", axes=[energy_axis])
    Cube3DMap = Map.from_geom(geom3D)
    Background3DMap = Map.from_geom(geom3D)
    Exposure3DMap = Map.from_geom(geom3D)

    geom4D = WcsGeom.create(skydir=(RAP[nP],DecP[nP]), binsz=0.04, width=(16, 14), coordsys="CEL", axes=[energy_axis, time_axis])
    Cube4DMap = Map.from_geom(geom4D)
    Background4DMap = Map.from_geom(geom4D)
    Exposure4DMap = Map.from_geom(geom4D)

    #counterPointings=0
    #geom = WcsGeom.create(skydir=(RAP[nP], DecP[nP]), binsz=0.04, width=(10, 10), coordsys="CEL", axes=[axis])
    pointing = SkyCoord(RAP[nP],DecP[nP], unit="deg", frame="fk5")
    FirstSpectralModel=len(grb.time_interval.value[mask1])
    for i in range(0,len(SimulationInterval)-1):
        efLivetimes=SimulationInterval[i+1]-SimulationInterval[i]
        print(efLivetimes)
        geom = WcsGeom.create(skydir=(RAP[nP], DecP[nP]), binsz=0.04, width=(16, 14),coordsys="CEL", axes=[axis])
        sky_model = SkyModel(spatial_model=spatial_model, spectral_model=grb.spectral_model[FirstSpectralModel+i])
        npred, counts, background, exposure= PointingSimulation(sky_model,axis,geom,pointing,efLivetimes*u.second,offset,irfs)

        Cube4DMap.data[i] = WcsNDMap(geom,counts).reproject(geom3D).data
        Background4DMap.data[i] = background.reproject(geom3D).data
        Exposure4DMap.data[i] = exposure.reproject(geom3D).data


    #Cube3DMap.data=counts
    #Background3DMap.data=background.data
    #Exposure3DMap.data=exposure.data

    print('========= RESULTS ========')
    print(Cube4DMap.data.shape)
    print(Background4DMap.data.shape)
    print(Exposure4DMap.data.shape)

    maps = {
        "counts": Cube4DMap,
        "background": Background4DMap,
        "exposure": Exposure4DMap,
    }

    print(maps)

    path = "/analysis_3d"
    print(path)
    print(outDir)
    fullpath=Path(outDir+path)

    fullpath.mkdir(exist_ok=True)

    # write maps
    maps["counts"].write(str(fullpath / "counts_singleObs_dynamic.fits"), overwrite=True)
    maps["background"].write(str(fullpath / "background_singleObs_dynamic.fits"), overwrite=True)
    maps["exposure"].write(str(fullpath / "exposure_singleObs_dynamic.fits"), overwrite=True)

def BinnedZenithAngle(zenithAngle):
    if(zenithAngle<55):
        zenIRF = 60
    elif(zenithAngle>55 and zenithAngle<30):
        zenIRF = 40
    elif(zenithAngle>30):
        zenIRF=20
    return zenIRF


def CubeSimulationSingleObservation_singleSpectrum(InputListrun,InputListMergerID,nP,Source,IRFfilename,zenithAngle,Pointings,datasetDir,outDir):

    path = '/' + InputListrun + '_' + InputListMergerID + '_' + str(nP)
    fullpath=Path(outDir+path)
    fullpath.mkdir(exist_ok=True)

    absorption = Absorption.read_builtin('dominguez')
    filepath= datasetDir+'/GammaCatalogV2.0/'+InputListrun + '_' + InputListMergerID.split('r')[-1] + ".fits"

    #filepath='/Users/mseglar/Documents/GitHub/CTASimulationsGW/GRBtemplates/SGRB001.fits'
    grb = GRB.from_fitsfile(filepath=filepath, absorption=absorption)

    PointingsCoord = SkyCoord(Pointings['RA(deg)'], Pointings['DEC(deg)'], frame='fk5', unit=(u.deg, u.deg))
    Source = SkyCoord(Source.ra.deg, Source.dec.deg, frame='fk5', unit=(u.deg, u.deg))

    RAP = PointingsCoord.ra.deg
    DecP = PointingsCoord.dec.deg

    #filenameS='ExampleSource.txt'


    spatial_model = SkyGaussian(lon_0=Source.ra.deg*u.deg, lat_0=Source.dec.deg*u.deg, sigma="0.2 deg")
    # Testing results if the source was it the center of the FoV
    #spatial_model = SkyGaussian(lon_0=PointingsCoord.ra.deg, lat_0=PointingsCoord.ra.deg * u.deg, sigma="0.2 deg")
    print()
    print('Source location is,',Source.ra.deg*u.deg,Source.dec.deg*u.deg)
    print('Pointing coordinates is,',RAP, DecP)
    print('GRB:',filepath)
    print()

    #for i in range(0,len(grb.spectral_model)):
    #    grb.spectral_model[i].plot(energy_range=(0.001, 10) * u.TeV)
    #plt.show()
    #plt.grid()
    #plt.savefig(str(fullpath)+'/Fits_Spectra.png')

    # Get Spectral model from fits!

    #plt.savefig('./CubeSimulation_Figures/LGRB_Evolution.png')

    #Select the IRFs depending on observation condition and site
    irfs = load_cta_irfs(IRFfilename)

    ###########################################
    #Cube Definition
    axis = MapAxis.from_edges(np.logspace(-3.0, 1.0, 20), unit="TeV", name="energy", interp="log")

    ###########################################

    livetime = 1800*u.second

    offset = Angle("2.5 deg")

    energy_axis = MapAxis.from_edges(np.logspace(-2.0, 1.0, 20), unit="TeV", name="energy", interp="log")
    geom3D = WcsGeom.create(skydir=(RAP,DecP), binsz=0.04, width=(16, 14), coordsys="CEL", axes=[energy_axis])
    Cube3DMap = Map.from_geom(geom3D)
    Background3DMap = Map.from_geom(geom3D)
    Exposure3DMap = Map.from_geom(geom3D)

    pointing = SkyCoord(RAP,DecP, unit="deg", frame="fk5")

    print(grb.spectral_model[0])

    #ToDo: Select the grb that corresponds to time of the observation (now, set to 6 for testing)
    spectral_model = grb.spectral_model[6]
    #spectral_model = PowerLaw(index=2.2, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV")

    sky_model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model)
    npred, counts, background, exposure= PointingSimulation(sky_model,axis,geom3D,pointing,livetime,offset,irfs)
    Cube3DMap.data = counts
    Background3DMap.data=background.data
    Exposure3DMap.data=exposure.data

    #print('========= RESULTS ========')
    #print(Cube3DMap.data.shape)
    #print(Background3DMap.data.shape)
    #print(Exposure3DMap.data.shape)

    maps = {
        "counts": Cube3DMap,
        "background": Background3DMap,
        "exposure": Exposure3DMap,
    }

    #print(maps)

    # write maps
    maps["counts"].write(str(fullpath / "counts_singleObs.fits"), overwrite=True)
    maps["background"].write(str(fullpath / "background_singleObs.fits"), overwrite=True)
    maps["exposure"].write(str(fullpath / "exposure_singleObs.fits"), overwrite=True)

    return fullpath