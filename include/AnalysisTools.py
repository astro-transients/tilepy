from gammapy.maps import WcsGeom, MapAxis, Map, WcsNDMap
from gammapy.extern.pathlib import Path
from astropy.convolution import Gaussian2DKernel
from gammapy.detect import TSMapEstimator, find_peaks
from astropy.io import ascii, fits
import matplotlib.pyplot as plt
from astropy.convolution import Tophat2DKernel
from gammapy.detect import compute_lima_image
from gammapy.image.models import SkyGaussian, SkyPointSource
from gammapy.detect import TSMapEstimator
from gammapy.scripts import SpectrumAnalysisIACT
import astropy.units as u
from regions import CircleSkyRegion
from astropy.coordinates import SkyCoord
import numpy as np
from gammapy.spectrum.models import PowerLaw
from gammapy.cube.models import SkyModel
from gammapy.cube import MapMaker, MapFit

def TSEstimation(maps):
    kernel = Gaussian2DKernel(2.5, mode="oversample")
    estimator = TSMapEstimator()
    images = estimator.run(maps, kernel)
    return images


def LikelihoodFit_Analysis_3DCube(dirName):
    path = "LikelihoodFit_Analysis_ResultPlots"
    fullpath=Path(dirName+'/analysis_3d/'+path)
    fullpath.mkdir(exist_ok=True)

    counter=0
    maps = {
            "counts": Map.read(dirName+'/analysis_3d/'+ "counts_singleObs_dynamic.fits").sum_over_axes(),
            "background": Map.read(dirName+'/analysis_3d/' + "background_singleObs_dynamic.fits").sum_over_axes(),
            "exposure": Map.read(dirName+'/analysis_3d/' + "exposure_singleObs_dynamic.fits").sum_over_axes()
    }


    # Using LIMA !!
    #kernel = Tophat2DKernel(2.5)
    #images = compute_lima_image(maps["counts"].sum_over_axes(), maps["background"].sum_over_axes(), kernel)

    #hdul = fits.open(str(path / "counts.fits"))
    #h1 = hdul[1]
    ##print(h1.columns)
    #print(h1.data['energy'])
    #print(h1.data['TIME'])
    #fig, ax = plt.subplots()
    #plt.plot(h1.data['TIME'],h1.data['energy'],'+')
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #plt.show()

    # print(images.meta['runtime'])

    #print(images["excess"])
    #print(images["significance"])

    print('------------------ TIME-AVERAGED TS COMPUTATION  ----------------------- ')

    images = TSEstimation(maps)

    print('Plotting results')

    plt.figure(figsize=(15, 5))
    maps["background"].sum_over_axes().plot(stretch="sqrt",add_cbar=True)
    plt.savefig(str(fullpath)+'/Background.png')

    plt.figure(figsize=(15, 5))
    maps["exposure"].sum_over_axes().plot(stretch="sqrt",add_cbar=True)
    plt.savefig(str(fullpath)+'/Exposure.png')


    plt.figure(figsize=(15, 5))
    maps["counts"].sum_over_axes().plot(stretch="sqrt",add_cbar=True)
    plt.savefig(str(fullpath)+'/Counts.png')
    plt.figure(figsize=(15, 5))

    #images["significance"].sum_over_axes().plot(stretch="sqrt",add_cbar=True)
    #sources = find_peaks(images["significance"].sum_over_axes(), threshold=4) #FOR LIMA
    images["sqrt_ts"].sum_over_axes().plot(stretch="sqrt",add_cbar=True)

    print('------------------ Find hotspot in TS ----------------------- ')

    sources = find_peaks(images["sqrt_ts"], threshold=3) # FOR TS

    print(sources)
    plt.gca().scatter(sources["ra"],sources["dec"],transform=plt.gca().get_transform("icrs"),color="none",edgecolor="white",marker="o",s=600,lw=1.5,);
    plt.savefig(str(fullpath)+'/Significance.png')

    sources = find_peaks(images["sqrt_ts"], threshold=0) # FOR TS
    plt.figure(figsize=(10, 5))

    mu, sigma = 0, 1  # mean and standard deviation
    count, bins, ignored =plt.hist(sources['value'].data,100, histtype='step',fill=False,density=True, stacked=True)
    plt.plot(bins, 1 / (sigma * np.sqrt(2 * np.pi)) *np.exp(- (bins - mu) ** 2 / (2 * sigma ** 2)),linewidth = 1, color = 'r')
    plt.savefig(str(fullpath)+'/histogram_SQRT_TS.png')


    '''
    print('====== SPECTRAL ANALYSIS ? =======')

    spatial_model = SkyPointSource(lon_0="339.29 deg", lat_0="-29.623 deg")
    spectral_model = PowerLaw(index=2.6, amplitude="5e-11 cm-2 s-1 TeV-1", reference="1 TeV")
    model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model)

    fit = MapFit(
                model=model,
                counts=maps["counts"],
                exposure=maps["exposure"],
                background=maps["background"],
                psf=psf_kernel,
                )
    result = fit.run()
    print(result)
    print(result.model.parameters.to_table())

    '''

def LikelihoodFit_Analysis_4DCube(dirname):
    #ToDo This function may not work. Gammapy doesnt support temporal analysis and anyways it may not be justified (too slow)

    path = dirname
    counter = 0
    mapsCounts = Map.read(str(path / "counts.fits"))
    for i in range(0, mapsCounts.data.shape[0]):
        counter = counter + 1
        print(i, np.nansum(mapsCounts.slice_by_idx({"TIME": i}).data))

    # Check a cut of the 4D cube and define a new geometry
    # Here, the slices have been chosen but they should be found

    newTime = mapsCounts.geom.get_axis_by_name('TIME').slice(slice(71, 110))
    newgeom = mapsCounts.geom.slice_by_idx({"TIME": slice(71, 110)})
    TS3DMap = Map.from_geom(newgeom)

    print('NEW GEOMETRY')
    print(newgeom)
    mapscounter = 0

    print('Number of temporal bins', mapsCounts.geom.get_axis_by_name("Time").nbin)
    '''
    #for i in range(71,mapsCounts.geom.get_axis_by_name("Time").nbin):
    print('------------------ TS CUBE COMPUTATION ----------------------- ')
    for i in range(71, 110):
        maps = {
            "counts": Map.read(str(path / "counts.fits")).slice_by_idx({"TIME":i}).sum_over_axes(),
            "background": Map.read(str(path / "background.fits")).slice_by_idx({"TIME":i}).sum_over_axes(),
            "exposure": Map.read(str(path / "exposure.fits")).slice_by_idx({"TIME":i}).sum_over_axes()
        }
        try:
            images = TSEstimation(maps)
            images["sqrt_ts"].sum_over_axes().plot(stretch="sqrt", add_cbar=True)
            TS3DMap.data[mapscounter]=images["sqrt_ts"].data
            mapscounter=mapscounter+1
            print('All good over here')
        except ValueError:
            print('Something went wrong')

    plt.figure(figsize=(15, 5))
    TS3DMap.sum_over_axes().plot(stretch="sqrt",add_cbar=True)
    sources = find_peaks(TS3DMap.sum_over_axes(), threshold=4) # FOR TS
    print(sources)
    plt.savefig('./CubeSimulation_Figures/TS_test.png')
    '''
    # Using LIMA !!
    # kernel = Tophat2DKernel(2.5)
    # images = compute_lima_image(maps["counts"].sum_over_axes(), maps["background"].sum_over_axes(), kernel)

    # hdul = fits.open(str(path / "counts.fits"))
    # h1 = hdul[1]
    ##print(h1.columns)
    # print(h1.data['energy'])
    # print(h1.data['TIME'])
    # fig, ax = plt.subplots()
    # plt.plot(h1.data['TIME'],h1.data['energy'],'+')
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    # plt.show()

    # print(images.meta['runtime'])

    # print(images["excess"])
    # print(images["significance"])

    print('------------------ TIME-AVERAGED TS COMPUTATION  ----------------------- ')

    maps = {
        "counts": Map.read(str(path / "counts.fits")).sum_over_axes(),
        "background": Map.read(str(path / "background.fits")).sum_over_axes(),
        "exposure": Map.read(str(path / "exposure.fits")).sum_over_axes()
    }
    plt.figure(figsize=(15, 5))
    maps["background"].sum_over_axes().plot(stretch="sqrt", add_cbar=True)
    plt.savefig('./CubeSimulation_Figures/Background.png')

    plt.figure(figsize=(15, 5))
    maps["exposure"].sum_over_axes().plot(stretch="sqrt", add_cbar=True)
    plt.savefig('./CubeSimulation_Figures/Exposure.png')

    images = TSEstimation(maps)

    maps["counts"].sum_over_axes().plot(stretch="sqrt", add_cbar=True)
    plt.savefig('./CubeSimulation_Figures/Counts.png')
    plt.figure(figsize=(15, 5))
    # images["significance"].sum_over_axes().plot(stretch="sqrt",add_cbar=True)
    # sources = find_peaks(images["significance"].sum_over_axes(), threshold=4) #FOR LIMA
    images["sqrt_ts"].sum_over_axes().plot(stretch="sqrt", add_cbar=True)
    sources = find_peaks(images["sqrt_ts"], threshold=4)  # FOR TS

    print(sources)
    plt.gca().scatter(sources["ra"], sources["dec"], transform=plt.gca().get_transform("icrs"), color="none",
                      edgecolor="white", marker="o", s=600, lw=1.5, );
    plt.savefig('./CubeSimulation_Figures/Significance.png')

    sources = find_peaks(images["sqrt_ts"], threshold=-3)  # FOR TS
    plt.figure(figsize=(10, 5))

    mu, sigma = 0, 1  # mean and standard deviation
    count, bins, ignored = plt.hist(sources['value'].data, 100, histtype='step', fill=False, density=True, stacked=True)
    plt.plot(bins, 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(- (bins - mu) ** 2 / (2 * sigma ** 2)), linewidth=1,
             color='r')
    plt.savefig(str(fullpath) + '/histogram_SQRT_TS.png')

    '''
    print('====== SPECTRAL ANALYSIS ? =======')

    spatial_model = SkyPointSource(lon_0="339.29 deg", lat_0="-29.623 deg")
    spectral_model = PowerLaw(index=2.6, amplitude="5e-11 cm-2 s-1 TeV-1", reference="1 TeV")
    model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model)

    fit = MapFit(
        model=model,
        counts=maps["counts"],
        exposure=maps["exposure"],
        background=maps["background"],
        psf=psf_kernel,
    )
    result = fit.run()
    print(result)
    print(result.model.parameters.to_table())

    '''