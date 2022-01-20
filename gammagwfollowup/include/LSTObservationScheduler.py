############################################################################
#       Authors: Monica Seglar-Arroyo, Halim Ashkar,  Fabian Schussler     #
#           LST observation scheduler of GBM alerts and GW events          #
############################################################################


from .GWHESSObservationScheduler import PGWinFoV, PGalinFOV
from .RankingObservationTimes import RankingTimes
from .PointingPlotting import PointingPlotting
from astropy.coordinates import SkyCoord
from .GWHESSPointingTools import Tools, LoadGalaxies
from astropy.io import fits, ascii
import time
import healpy as hp
import numpy as np
from astropy import units as u
import datetime
import os


def getdate(x):
    if isinstance(x, datetime.datetime):
        return x
    elif isinstance(x, str):
        return datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    else:
        print("ERROR: something is wrong with the format of the date: ", x)
        return None


def GetSchedule_GW(URL, date, observatory,datasetDir,outDir):
    FOLLOWUP = True
    targetType = 'GW_Pointing'

    filename = URL.split("/")[-1]
    print("The filename is ", filename)
    try:
        fits_map_url = URL
        ##fits_map_url = self.What["GW_SKYMAP"]['skymap_fits']['value']
        command = 'curl %s -o %s' % (fits_map_url, filename)
        print(command)
        os.system(command)
    except x:
        warn = "Caught exeption: %s" % x
        print(warn)
        pass

    prob = []
    distmu = []
    distsigma = []
    distnorm = []
    tdistmean = 0
    fitsfile = fits.open(filename)
    if (fitsfile[1].header['TFIELDS'] == 4):
        prob, distmu, distsigma, distnorm = hp.read_map(filename,
                                                        field=range(4))
        tdistmean = fitsfile[1].header['DISTMEAN']
    else:
        prob = hp.read_map(filename, field=range(1))

    has3D = True
    if len(distnorm) == 0:
        has3D = False

    # Halim: Check if max pix is in the Galactic Plane
    npix = len(prob)
    NSide = hp.npix2nside(npix)
    MaxPix = np.argmax(prob)
    MaxTheta, MaxPhi = hp.pix2ang(NSide, MaxPix)
    raMax = np.rad2deg(MaxPhi)
    decMax = np.rad2deg(0.5 * np.pi - MaxTheta)

    c_icrs = SkyCoord(raMax, decMax, frame='fk5', unit=(u.deg, u.deg))

    InsidePlane = False
    InsidePlane = Tools.GalacticPlaneBorder(c_icrs)
    if InsidePlane:
        has3D = False

    #print("the distance of the event is:", tdistmean)
    if tdistmean > 150:
        has3D = False

    name = URL.split('/')[-3]

    print("===========================================================================================")

    if has3D:
        # date = datetime.datetime.now()
        # date = datetime.datetime.now(timezone.utc)
        # date = datetime.datetime.utcnow()
        # date = datetime.datetime(2019,7,28,20,20,20)
        done = "False"
        galaxies = datasetDir + "/GLADE.txt"
        params = "./configs/PGalinFoV.ini"


        #ObservationTime = getdate(date)
        ObservationTime = date
        PointingsFile = done
        galFile = galaxies
        parameters = params
        outputDir = "%s/%s" % (outDir, name)

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        dirName = '%s/PGallinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Starting the LST PGalinFoV pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galFile)
        print("Parameters: ", parameters)
        print("Output: ", outputDir)

        SuggestedPointings, cat = PGalinFoV(filename, ObservationTime, PointingsFile, galFile, parameters, dirName,observatory)

    else:
        # date = datetime.datetime.now()
        # date = datetime.datetime.now(timezone.utc)
        # date = datetime.datetime.utcnow()
        # date = datetime.datetime(2021,6,8,14,20,20)
        done = "False"
        galaxies = "GLADE.txt"
        params = "./configs/PGWinFoV.ini"

        ObservationTime = date
        PointingsFile = done
        galFile = galaxies
        parameters = params
        outputDir = "%s/%s" % (outDir, name)

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        dirName = '%s/PGWinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Starting the LST PGwinFoV pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galFile)
        print("Parameters: ", parameters)
        print("Output: ", outputDir)

        SuggestedPointings, t0 = PGWinFoV(filename, ObservationTime, PointingsFile, galFile, parameters, dirName, observatory)

    print(SuggestedPointings)
    print("===========================================================================================")
    print()

    if (len(SuggestedPointings) != 0):
        FOLLOWUP = True
        outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
        ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
        print()
        cat = LoadGalaxies(galFile)
        RankingTimes(ObservationTime, filename, cat, parameters, observatory, targetType, dirName,
                     '%s/SuggestedPointings_GWOptimisation.txt' % dirName)
        PointingPlotting(filename, cat,observatory,params, name, dirName, '%s/SuggestedPointings_GWOptimisation.txt' % dirName)
    else:
        FOLLOWUP = False
        print('No observations are scheduled')

    return SuggestedPointings, FOLLOWUP


def GetSchedule_GBM(URL, date, observatory,datasetDir,outDir):
    FOLLOWUP = True
    targetType = 'GBM_Pointing'


    filename = URL.split("/")[-1]
    filename = filename.split(".")[0]
    filename = "./" + filename + ".fits"
    # filename = "glg_healpix_all_bn211130636_2.fits"
    print("The filename is ", filename)
    try:
        fits_map_url_intial = URL
        fits_map_url1 = fits_map_url_intial.split("/")
        fits_map_url2 = fits_map_url_intial.split("_")[-1]
        # fits_map_url1[-1] = ""
        i = 0
        fits_map_url = ""
        for i in range(len(fits_map_url1) - 1):
            fits_map_url += fits_map_url1[i] + "/"
        fits_map_url += "glg_healpix_all" + "_" + fits_map_url2.split(".")[0] + ".fit"

        # command = 'curl %s -o %s' % (fits_map_url, filename)
        # print(command)
        # os.system(command)
        # should fix this issue to save the fits file ... then comment the following line.
        filename = fits_map_url
    except:
        warn = "Caught exception: "
        print(warn)
        pass

    prob = []
    distmu = []
    distsigma = []
    distnorm = []
    tdistmean = 0
    print('Filename is ', filename)

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
            if d > 20:
                print("Waited for 20 minutes... can't wait anymore... I'm leaving")
                break

    if (fitsfile[1].header['TFIELDS'] == 4):
        prob, distmu, distsigma, distnorm = hp.read_map(filename,
                                                        field=range(4))
        tdistmean = fitsfile[1].header['DISTMEAN']
    else:
        prob = hp.read_map(filename, field=range(1))

    has3D = True
    if len(distnorm) == 0:
        has3D = False

    # Check if max pix is in the Galactic Plane
    npix = len(prob)
    NSide = hp.npix2nside(npix)
    MaxPix = np.argmax(prob)
    MaxTheta, MaxPhi = hp.pix2ang(NSide, MaxPix)
    raMax = np.rad2deg(MaxPhi)
    decMax = np.rad2deg(0.5 * np.pi - MaxTheta)
    c_icrs = SkyCoord(raMax, decMax, frame='fk5', unit=(u.deg, u.deg))

    InsidePlane = False
    InsidePlane = Tools.GalacticPlaneBorder(c_icrs)
    if InsidePlane:
        has3D = False

    #print(" The distance of the event is:", tdistmean)
    if tdistmean > 150:
        has3D = False

    # filename=args.name
    name = URL.split('/')[-3]

    if has3D:
        # date = datetime.datetime.now()
        # date = datetime.datetime.now(timezone.utc)
        # date = datetime.datetime.utcnow()
        # date = datetime.datetime(2021,3,3,14,20,20)
        done = "False"
        galaxies = datasetDir + "/GLADE.txt"
        params = "./configs/PGalinFoV.ini"

        ObservationTime = date
        PointingsFile = done
        galFile = galaxies
        parameters = params
        outputDir = "%s/%s" % (outDir, name)

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        dirName = '%s/PGallinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Starting the LST PGalinFoV pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galFile)
        print("Parameters: ", parameters)
        print("Output: ", outputDir)

        SuggestedPointings, cat = PGalinFoV(filename, ObservationTime, PointingsFile, galFile, parameters, dirName, observatory)

    else:
        # date = datetime.datetime.now()
        # date = datetime.datetime.now(timezone.utc)
        # date = datetime.datetime.utcnow()
        # date = datetime.datetime(2021, 11, 30, 15, 15, 15)
        done = "False"
        galaxies = datasetDir + "/GLADE.txt"
        params = "./configs/PGWinFoV.ini"

        ObservationTime = date
        PointingsFile = done
        galFile = galaxies
        parameters = params
        outputDir = "%s/%s" % (outDir, name)

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        dirName = '%s/PGWinFoV' % outputDir

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        print("===========================================================================================")
        print("Starting the LST PGwinFoV pointing calculation with the following parameters\n")
        print("Filename: ", name)
        print("Date: ", ObservationTime)
        print("Previous pointings: ", PointingsFile)
        print("Catalog: ", galFile)
        print("Parameters: ", parameters)
        print("Output: ", outputDir)

        SuggestedPointings, t0 = PGWinFoV(filename, ObservationTime, PointingsFile, galFile, parameters, dirName,observatory)

        print(SuggestedPointings)
        print("===========================================================================================")
        print()

    if (len(SuggestedPointings) != 0):
        FOLLOWUP = True
        outfilename = '%s/SuggestedPointings_GWOptimisation.txt' % dirName
        ascii.write(SuggestedPointings, outfilename, overwrite=True, fast_writer=False)
        print()
        cat = LoadGalaxies(galFile)
        RankingTimes(ObservationTime, filename, cat, parameters,observatory, targetType, dirName,
                     '%s/SuggestedPointings_GWOptimisation.txt' % dirName)
        PointingPlotting(filename, cat, observatory, params, name, dirName, '%s/SuggestedPointings_GWOptimisation.txt' % dirName)
    else:
        FOLLOWUP = False
        print('No observations are scheduled')

    #print ("followup", FOLLOWUP)
    return SuggestedPointings, FOLLOWUP


