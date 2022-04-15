#!/usr/bin/env python
# Script name: ObservingTimes.py
# Version: 7.0 (August 2020)
# Author(s): B. Patricelli (barbara.patricelli@pi.infn.it)

# Modified by: Monica Seglar-Arroyo
#######################################################
# Import
import os
import sys
from astropy.io import fits
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import interp1d


#######################################################

#Get Script name
ScriptName = os.path.split(sys.argv[0])[1].split('.')[0]


def ParseCatalog(catFileName,Separator="\t"):
    #"Parse a catalog of sources"
    #print('Parsing catalog file',catFileName)
    datalist=open(catFileName).readlines()
    header=(datalist[1][1:-1]).split(Separator)
    #print('Read',len(header),'columns:')
    #print(header)
    SourceList=[]
    for di in range(3,len(datalist)):
        Source={}
        line=datalist[di].strip('\n').split(Separator)
        if (line[-1]==''):
            line=line[:-1]
        for li in range(len(line)):
            Source[header[li]]=line[li].strip('')
        SourceList.append(Source)

    #print('Read',len(SourceList),'sources')
    return SourceList
#
# def TableMoonSun(startTime,observatory):
#
#     initialframe = AltAz(obstime=inputtime, location=observatory)
#     suninitial = get_sun(inputtime).transform_to(initialframe)


# Defining the function to fit the sensitivity vs observing time
def fs(x, a, b, c, d):
    return a * (x ** b + x ** c + x ** d)


def spectrum(x):
    return (x/1)**(-2.2)


def ObtainObservingTimes(totaltime,delayAlert,run,id,observatory,dirName,zenith):
  dt=10
  x1=2
  x2=10000
  #tslew=30 #talert should include initial slewing
  tslewbtwn=10

  #frame = co.AltAz(obstime=tAlert, location=observatory.Location)
  #targetCoord_map = co.SkyCoord(targetCoor.ra, targetCoor.dec, frame='fk5', unit=(u.deg, u.deg))
  #altaz = targetCoord_map.transform_to(frame)
  #altitude=altaz.alt.deg
  #print(zenith)

  ########################################################################
  # Defining which sensitivity to use for each event (CTA site and zenith)
  ########################################################################
  #site=site[j].strip()

  if observatory.Name=="North":
        
    if int(zenith)==20:
        # sensitivity CTA North, zenit=20 deg
        sensifilenorthz20 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-North_z20_0.5h.txt"  # "sensitivity-5sigma_irf-North_z20_0.5.txt"

        SensiListNorthz20 = ParseCatalog(sensifilenorthz20)
        ndatasnz20 = len(SensiListNorthz20)

        xn20 = []
        yn20 = []
        for i in range(0, ndatasnz20):
            obstimen20 = float(SensiListNorthz20[i]['Obs time'])
            pfluxn20 = float(SensiListNorthz20[i]['photon_flux'])
            efluxn20 = float(SensiListNorthz20[i]['energy_flux'])

            xn20.append(obstimen20)
            yn20.append(pfluxn20)

        #parn20, pcovn20 = curve_fit(fs, xn20, yn20)
        fsn20=interp1d(xn20,yn20)
        def sensi(x):
            return fsn20(x)
            #return fs(x,*parn20)
    elif int(zenith)==40:
        ###############
        # sensitivity CTA North, zenit=40 deg
        sensifilenorthz40 = dirName +"/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-North_z40_0.5h.txt"

        SensiListNorthz40 = ParseCatalog(sensifilenorthz40)
        ndatasnz40 = len(SensiListNorthz40)

        xn40 = []
        yn40 = []
        for i in range(0, ndatasnz40):
            obstimen40 = float(SensiListNorthz40[i]['Obs time'])
            pfluxn40 = float(SensiListNorthz40[i]['photon_flux'])
            efluxn40 = float(SensiListNorthz40[i]['energy_flux'])

            xn40.append(obstimen40)
            yn40.append(pfluxn40)

        fsn40=interp1d(xn40,yn40)
        #parn40, pcovn40 = curve_fit(fs, xn40, yn40)

        def sensi(x):
            return fsn40(x)
            #return fs(x,*parn40)
    else:
        ###############
        # sensitivity CTA North, zenit=60 deg
        sensifilenorthz60 = dirName +"/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-North_z60_0.5h.txt"

        SensiListNorthz60 = ParseCatalog(sensifilenorthz60)
        ndatasnz60 = len(SensiListNorthz60)

        xn60 = []
        yn60 = []
        for i in range(0, ndatasnz60):
            obstimen60 = float(SensiListNorthz60[i]['Obs time'])
            pfluxn60 = float(SensiListNorthz60[i]['photon_flux'])
            efluxn60 = float(SensiListNorthz60[i]['energy_flux'])

            xn60.append(obstimen60)
            yn60.append(pfluxn60)

        #parn60, pcovn60 = curve_fit(fs, xn60, yn60)
        fsn60=interp1d(xn60,yn60)
        
        def sensi(x):
            return fsn60(x)
            #return fs(x,*parn60)
  else:
    if int(zenith)==20:
        ###############
        # sensitivity CTA South, zenit=20 deg

        sensifilesouthz20 = dirName +"/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-South_z20_0.5h.txt"

        SensiListSouthz20 = ParseCatalog(sensifilesouthz20)
        ndatassz20 = len(SensiListSouthz20)

        xs20 = []
        ys20 = []
        for i in range(0, ndatassz20):
            obstimes20 = float(SensiListSouthz20[i]['Obs time'])
            pfluxs20 = float(SensiListSouthz20[i]['photon_flux'])
            efluxs20 = float(SensiListSouthz20[i]['energy_flux'])

            xs20.append(obstimes20)
            ys20.append(pfluxs20)

        #pars20, pcovs20 = curve_fit(fs, xs20, ys20)
        fss20=interp1d(xs20,ys20)
        def sensi(x):
            return fss20(x)
            #return fs(x,*pars20)
    elif int(zenith)==40:
        ###############
        # sensitivity CTA South, zenit=40 deg
        sensifilesouthz40 = dirName+"/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-South_z40_0.5h.txt"

        SensiListSouthz40 = ParseCatalog(sensifilesouthz40)
        ndatassz40 = len(SensiListSouthz40)

        xs40 = []
        ys40 = []
        for i in range(0, ndatassz40):
            obstimes40 = float(SensiListSouthz40[i]['Obs time'])
            pfluxs40 = float(SensiListSouthz40[i]['photon_flux'])
            efluxs40 = float(SensiListSouthz40[i]['energy_flux'])

            xs40.append(obstimes40)
            ys40.append(pfluxs40)

        fss40=interp1d(xs40,ys40)
        #pars40, pcovs40 = curve_fit(fs, xs40, ys40)
        def sensi(x):
            return fss40(x)
            #return fs(x,*pars40)
    else:
        ###############
        # sensitivity CTA South, zenit=40 deg
        sensifilesouthz60 = dirName +"/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-South_z60_0.5h.txt"

        SensiListSouthz60 = ParseCatalog(sensifilesouthz60)
        ndatassz60 = len(SensiListSouthz60)

        xs60 = []
        ys60 = []
        for i in range(0, ndatassz60):
            obstimes60 = float(SensiListSouthz60[i]['Obs time'])
            pfluxs60 = float(SensiListSouthz60[i]['photon_flux'])
            efluxs60 = float(SensiListSouthz60[i]['energy_flux'])

            xs60.append(obstimes60)
            ys60.append(pfluxs60)

        fss60=interp1d(xs60,ys60)
        #pars60, pcovs60 = curve_fit(fs, xs60, ys60)
        def sensi(x):
            return fss60(x)
            #return fs(x,*pars60)

#######################################################
# Associating a GRB to the BNS mergers
#######################################################

  InputGRB= dirName +"/GammaCatalogV2.0/%s_%s.fits"%(run,id.split('ger')[1])
  hdu_list = fits.open(InputGRB)
# hdu_list.info()

  datalc = hdu_list[3].data
  datatime = hdu_list[2].data
  dataenergy = hdu_list[1].data


  lc=datalc.field(0)
  time=datatime.field(0)
  energy=dataenergy.field(0)
  spec=datalc[0]
  #Norm=spec[0]

  # Spectrum
  # Please note that we are assuming here that the spectral index does not evolve with time,
  # and that it is equal to the one of GRB 090510 (this is true for the first version of the GRB catalog)


  # Integral of the spectrum
  intl,errl=integrate.quad(lambda x: spectrum(x),x1,x2) #GeV

  # Interpolation of the flux with time
  flux = interpolate.interp1d(time,lc)

  #################################################
  # Starting the procedure
  #################################################


# starting time of the observations
  tstart=delayAlert#+tslew # (s)
    
  # Loop over consecutive pointings
  #for n in range(1,totaltime):
  Tstart=[]
  Obst=[]
  # increasing the observing time to get the GRB detection
  if tstart >=totaltime:
    return(tstart,-1)

  # Loop over consecutive pointings
  for n in range(0, totaltime-1):
    if (n >= 1):  # Accounting for slewing time
        tstart = tstart + 10
        t = t + 10
    for m in range(0,totaltime,1):
        t=tstart+m+dt
        obst=t-tstart

        fluencen, errorn = integrate.quad(lambda x: flux(x), tstart, t) #ph/cm2/GeV
        averagefluxn=fluencen*intl/obst #ph/cm2/s
            
        if averagefluxn >sensi(obst):
                
            if tstart+obst < totaltime:
                Tstart.append(tstart)
                Obst.append(obst)
                tstart = tstart + obst
                break
                #outfile2n.write("%d\t%d\t%d\t%d\n"%(n,tstart,tend,obst))
        if tstart+obst > totaltime:
            return(Tstart, Obst)
            #return (tstar,-1)

        #tstart=tstart+obst+tslewbtwn


def ObtainSingleObservingTimes(TotalExposure, tstart, interObsSlew , run, id, observatory,dirName, zenith):

    E1 = 2
    E2 = 100000
    steps = 10  # Check every 10 seconds

    # tslew=30

    # frame = co.AltAz(obstime=tAlert, location=observatory.Location)
    # targetCoord_map = co.SkyCoord(targetCoor.ra, targetCoor.dec, frame='fk5', unit=(u.deg, u.deg))
    # altaz = targetCoord_map.transform_to(frame)
    # altitude=altaz.alt.deg
    # print(zenith)

    ########################################################################
    # Defining which sensitivity to use for each event (CTA site and zenith)
    ########################################################################
    # site=site[j].strip()

    if observatory.Name == "North":

        if int(zenith) == 20:
            # sensitivity CTA North, zenit=20 deg
            sensifilenorthz20 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-North_z20_0.5h.txt"  # "sensitivity-5sigma_irf-North_z20_0.5.txt"

            SensiListNorthz20 = ParseCatalog(sensifilenorthz20)
            ndatasnz20 = len(SensiListNorthz20)

            xn20 = []
            yn20 = []
            for i in range(0, ndatasnz20):
                obstimen20 = float(SensiListNorthz20[i]['Obs time'])
                pfluxn20 = float(SensiListNorthz20[i]['photon_flux'])
                efluxn20 = float(SensiListNorthz20[i]['energy_flux'])

                xn20.append(obstimen20)
                yn20.append(pfluxn20)

            # parn20, pcovn20 = curve_fit(fs, xn20, yn20)
            fsn20 = interp1d(xn20, yn20)

            def sensi(x):
                return fsn20(x)
                # return fs(x,*parn20)
        elif int(zenith) == 40:
            ###############
            # sensitivity CTA North, zenit=40 deg
            sensifilenorthz40 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-North_z40_0.5h.txt"

            SensiListNorthz40 = ParseCatalog(sensifilenorthz40)
            ndatasnz40 = len(SensiListNorthz40)

            xn40 = []
            yn40 = []
            for i in range(0, ndatasnz40):
                obstimen40 = float(SensiListNorthz40[i]['Obs time'])
                pfluxn40 = float(SensiListNorthz40[i]['photon_flux'])
                efluxn40 = float(SensiListNorthz40[i]['energy_flux'])

                xn40.append(obstimen40)
                yn40.append(pfluxn40)

            fsn40 = interp1d(xn40, yn40)

            # parn40, pcovn40 = curve_fit(fs, xn40, yn40)

            def sensi(x):
                return fsn40(x)
                # return fs(x,*parn40)
        else:
            ###############
            # sensitivity CTA North, zenit=60 deg
            sensifilenorthz60 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-North_z60_0.5h.txt"

            SensiListNorthz60 = ParseCatalog(sensifilenorthz60)
            ndatasnz60 = len(SensiListNorthz60)

            xn60 = []
            yn60 = []
            for i in range(0, ndatasnz60):
                obstimen60 = float(SensiListNorthz60[i]['Obs time'])
                pfluxn60 = float(SensiListNorthz60[i]['photon_flux'])
                efluxn60 = float(SensiListNorthz60[i]['energy_flux'])

                xn60.append(obstimen60)
                yn60.append(pfluxn60)
            # parn60, pcovn60 = curve_fit(fs, xn60, yn60)
            fsn60 = interp1d(xn60, yn60)

            def sensi(x):
                return fsn60(x)
                # return fs(x,*parn60)
    else:
        if int(zenith) == 20:
            ###############
            # sensitivity CTA South, zenit=20 deg

            sensifilesouthz20 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-South_z20_0.5h.txt"

            SensiListSouthz20 = ParseCatalog(sensifilesouthz20)
            ndatassz20 = len(SensiListSouthz20)

            xs20 = []
            ys20 = []
            for i in range(0, ndatassz20):
                obstimes20 = float(SensiListSouthz20[i]['Obs time'])
                pfluxs20 = float(SensiListSouthz20[i]['photon_flux'])
                efluxs20 = float(SensiListSouthz20[i]['energy_flux'])

                xs20.append(obstimes20)
                ys20.append(pfluxs20)

            # pars20, pcovs20 = curve_fit(fs, xs20, ys20)
            fss20 = interp1d(xs20, ys20)

            def sensi(x):
                return fss20(x)
                # return fs(x,*pars20)
        elif int(zenith) == 40:
            ###############
            # sensitivity CTA South, zenit=40 deg
            sensifilesouthz40 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-South_z40_0.5h.txt"

            SensiListSouthz40 = ParseCatalog(sensifilesouthz40)
            ndatassz40 = len(SensiListSouthz40)

            xs40 = []
            ys40 = []
            for i in range(0, ndatassz40):
                obstimes40 = float(SensiListSouthz40[i]['Obs time'])
                pfluxs40 = float(SensiListSouthz40[i]['photon_flux'])
                efluxs40 = float(SensiListSouthz40[i]['energy_flux'])

                xs40.append(obstimes40)
                ys40.append(pfluxs40)

            fss40 = interp1d(xs40, ys40)

            # pars40, pcovs40 = curve_fit(fs, xs40, ys40)
            def sensi(x):
                return fss40(x)
                # return fs(x,*pars40)
        else:
            ###############
            # sensitivity CTA South, zenit=40 deg
            sensifilesouthz60 = dirName + "/grbsens_output_v2_May_2021/grbsens-5.0sigma_t1s-t16384s_irf-South_z60_0.5h.txt"

            SensiListSouthz60 = ParseCatalog(sensifilesouthz60)
            ndatassz60 = len(SensiListSouthz60)

            xs60 = []
            ys60 = []
            for i in range(0, ndatassz60):
                obstimes60 = float(SensiListSouthz60[i]['Obs time'])
                pfluxs60 = float(SensiListSouthz60[i]['photon_flux'])
                efluxs60 = float(SensiListSouthz60[i]['energy_flux'])

                xs60.append(obstimes60)
                ys60.append(pfluxs60)

            fss60 = interp1d(xs60, ys60)

            # pars60, pcovs60 = curve_fit(fs, xs60, ys60)
            def sensi(x):
                return fss60(x)
                # return fs(x,*pars60)

    #######################################################
    # Associating a GRB to the BNS mergers
    #######################################################

    InputGRB = dirName + "/GammaCatalogV2.0/%s_%s.fits" % (run, id.split('ger')[1])
    hdu_list = fits.open(InputGRB)
    # hdu_list.info()

    datalc = hdu_list[3].data
    datatime = hdu_list[2].data
    dataenergy = hdu_list[1].data

    lc = datalc.field(0)
    time = datatime.field(0)
    energy = dataenergy.field(0)
    spec = datalc[0]
    # Norm=spec[0]

    # Spectrum
    # Please note that we are assuming here that the spectral index does not evolve with time,
    # and that it is equal to the one of GRB 090510 (this is true for the first version of the GRB catalog)

    # Integral of the spectrum
    intl, errl = integrate.quad(lambda x: spectrum(x), E1, E2)  # GeV

    # Interpolation of the flux with time
    flux = interpolate.interp1d(time, lc)

    #################################################
    # Starting the procedure
    #################################################

    # starting time of the observations
    # in the case of i!=0, this time should be the already performed observing time

    # increasing the observing time to get the GRB detection

    totaltime = TotalExposure + tstart            # The total time is defined for a night
    maxcounter = totaltime-(tstart+interObsSlew)
    MaxObservationTime = 16380      # Max. observation time from GRBsens files
    for j in range(1, maxcounter):
        #t = tstart + m + interObsSlew
        t = tstart + j*steps + interObsSlew
        print("INCREASE from tstart = ", tstart," to t = ",t, "Total time = ",totaltime)

        fluencen, errorn = integrate.quad(lambda x: flux(x), tstart, t)  # ph/cm2/GeV
        obst = t - (tstart+interObsSlew)
        #print("Tobs mi reina",obst)
        averagefluxn = fluencen * intl / obst  # ph/cm2/s
        #print(' Averageflux achieved =', averagefluxn)
        #print(" CTA sensi(obst) = ", sensi(obst))

        if averagefluxn > sensi(obst): # The detection is possible
            if tstart + interObsSlew + obst < totaltime:
                obsCase = 'True'
                return obsCase, obst
            else:
                obsCase = 'TimeNotEnough'
                obst = False
                return obsCase, obst
        if (tstart + interObsSlew + obst > totaltime): # There is not enough time to obtain a detection
            print("Flux didnt reach the minimum value and the total time is OVER")
            obst = False
            obsCase = 'TimeNotEnough'
            return obsCase, obst
        if (obst+interObsSlew+steps>MaxObservationTime): # There is not enough time to obtain a detection
            print("Flux didnt reach the minimum value for this iteration")
            obst = False
            obsCase = 'TimeNotEnoughIte'
            return obsCase, obst

