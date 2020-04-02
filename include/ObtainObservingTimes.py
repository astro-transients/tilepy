#!/usr/bin/env python
# Script name: ObservingTimes.py
# Version: 6.0 (March 2020)
# Author(s): B. Patricelli (barbara.patricelli@pi.infn.it)


#######################################################
# Import
import os, sys
sys.path.append('../')
import matplotlib.pyplot as plt
from optparse import OptionParser
from astropy.io import fits
from gwUtilities import *

from scipy import interpolate
from scipy import integrate
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import math
import numpy as np

#######################################################

#Get Script name
ScriptName = os.path.split(sys.argv[0])[1].split('.')[0]

#######################################################
# Main 
#######################################################

def ObtainObservingTimes(totaltime,talert):
  dt=10
  x1=20
  x2=10000
  tslew=30
  tslewbtwn=10

  '''
  usg = "\033[1;31m%prog [ options] inputFile\033[1;m \n"

  desc = "\033[34mThis is a sample script\033[0m \n"

  parser=OptionParser(description=desc,usage=usg)
  
  parser.add_option("-i","--InputBNS",default="BNS-GW-Time_onAxis5deg-test.txt",help="File with the BNS mergers, zenith and CTA site")
  parser.add_option("-d","--deltat",type='float',default=10.,help="time bin width (s)")
  parser.add_option("-t","--totaltime",type='int',default=10210,help="maximum integration time")
  
  parser.add_option("-l","--lowerenergy",type=float,default=20,help="lower energy limit (GeV)")
  parser.add_option("-m","--higherenergy",type=float,default=10000,help="higher energy limit (GeV)")
 
  parser.add_option("-s","--slewtime",type='int',default=30,help="initial slewing time (s)")
  parser.add_option("-b","--slewtimebtwn",type='int',default=10,help="slewing time between pointings (s)")
  parser.add_option("-a","--alerttime",type='int',default=180,help="latency for the GW alert (s)")


  (options,args) = parser.parse_args()
  if (len(args)!=0):
      parser.error("incorrect number of arguments. Use -h for help")

  InputData = options.InputBNS
  dt=options.deltat
  totaltime=options.totaltime
  x1=options.lowerenergy
  x2=options.higherenergy
  tslew=options.slewtime
  tslewbtwn=options.slewtimebtwn
  talert=options.alerttime
  '''

  #####################################################

  print('************************************')
  print('** '+ScriptName)
  print('************************************')

  print('** Options:')

  #List parameters and args
  print('\n**Input parameters:')
  for key, val in parser.values.__dict__.iteritems():
    print(key,':', val)
    if (val==None):
        print('\nERROR!',key,'is a required option! Exit.\nType ',ScriptName+'.py -h for help\n')
        sys.exit(1)

  
#################################################################
# Defining the function to fit the sensitivity vs observing time
#################################################################

  def fs(x,a,b,c,d):
    return a*(x**b+x**c+x**d)
  
#################################################
# Reading the file with the sensitivity and fit
#################################################

###############
# sensitivity CTA North, zenit=20 deg
  sensifilenorthz20= "sensitivity_CO30000_North_0.5h.txt" #"sensitivity-5sigma_irf-North_z20_0.5.txt"

  SensiListNorthz20=ParseCatalog(sensifilenorthz20)
  ndatasnz20=len(SensiListNorthz20)

  xn20=[]
  yn20=[]
  for i in range(0,ndatasnz20):
      obstimen20=float(SensiListNorthz20[i]['Obs time'])
      pfluxn20=float(SensiListNorthz20[i]['photon_flux'])
      efluxn20=float(SensiListNorthz20[i]['energy_flux'])
      
      xn20.append(obstimen20)
      yn20.append(pfluxn20)

  parn20, pcovn20 = curve_fit(fs, xn20, yn20)


###############
# sensitivity CTA North, zenit=40 deg
  sensifilenorthz40="sensitivity_CO30000_South_0.5h.txt"

  SensiListNorthz40=ParseCatalog(sensifilenorthz40)
  ndatasnz40=len(SensiListNorthz40)

  xn40=[]
  yn40=[]
  for i in range(0,ndatasnz40):
      obstimen40=float(SensiListNorthz40[i]['Obs time'])
      pfluxn40=float(SensiListNorthz40[i]['photon_flux'])
      efluxn40=float(SensiListNorthz40[i]['energy_flux'])
      
      xn40.append(obstimen40)
      yn40.append(pfluxn40)

  parn40, pcovn40 = curve_fit(fs, xn40, yn40)

###############
# sensitivity CTA South, zenit=20 deg
  sensifilesouthz20="sensitivity_CO30000_North_0.5h.txt"

  SensiListSouthz20=ParseCatalog(sensifilesouthz20)
  ndatassz20=len(SensiListSouthz20)

  xs20=[]
  ys20=[]
  for i in range(0,ndatassz20):
      obstimes20=float(SensiListSouthz20[i]['Obs time'])
      pfluxs20=float(SensiListSouthz20[i]['photon_flux'])
      efluxs20=float(SensiListSouthz20[i]['energy_flux'])
      
      xs20.append(obstimes20)
      ys20.append(pfluxs20)

  pars20, pcovs20 = curve_fit(fs, xs20, ys20)


###############
# sensitivity CTA South, zenit=40 deg
  sensifilesouthz40="sensitivity_CO30000_South_0.5h.txt"

  SensiListSouthz40=ParseCatalog(sensifilesouthz40)
  ndatassz40=len(SensiListSouthz40)

  xs40=[]
  ys40=[]
  for i in range(0,ndatassz40):
      obstimes40=float(SensiListSouthz40[i]['Obs time'])
      pfluxs40=float(SensiListSouthz40[i]['photon_flux'])
      efluxs40=float(SensiListSouthz40[i]['energy_flux'])
      
      xs40.append(obstimes40)
      ys40.append(pfluxs40)

  pars40, pcovs40 = curve_fit(fs, xs40, ys40)


#######################################################################

  
#########################################
# Creating the output directories with the obs time
  os.system('mkdir ObsTimes')
#########################################

#######################################################
# Reading the file with the BNS mergers id numbers,
# together with zenith angle and CTA site
#######################################################

  infile=open(InputData)
  run,id,zenith,site=np.genfromtxt(InputData,usecols=(0, 1, 2, 3),unpack=True,dtype='str')
  print(zenith)
  '''    
  if infile:
      for line in infile:
        print(line)
        print(line.split("\t")[:-1])
        a,b,c,d=[str(t) for t in line.split("\t")]
        run.append(a)
        id.append(b)
        zenith.append(c)
        site.append(d)


  '''
  ndata=len(run)

  # Loop over the BNS mergers
  for j in range(0,ndata):

########################################################################
# Defining which sensitivity to use for each event (CTA site and zenith)
########################################################################
    ctasite=site[j].strip()

    if ctasite=="North":
        
        if int(zenith[j])==20:
            def sensi(x):
                return fs(x,*parn20)
        else:
            def sensi(x):
                return fs(x,*parn40)

    else:
        if int(zenith[j])==20:
            def sensi(x):
                return fs(x,*pars20)
        else:
            def sensi(x):
                return fs(x,*pars40)

#######################################################
# Associating a GRB to the BNS mergers
#######################################################

    InputGRB="../../../../dataset/GammaCatalogV1.0/%s_%s.fits"%(run[j],id[j].split('ger')[1])
    hdu_list = fits.open(InputGRB)
#    hdu_list.info()

    datalc = hdu_list[3].data
    datatime = hdu_list[2].data
    dataenergy = hdu_list[1].data


    lc=datalc.field(0)
    time=datatime.field(0)
    energy=dataenergy.field(0)
    spec=datalc[0]
#    Norm=spec[0]

# Spectrum
# Please note that we are assuming here that the spectral index does not evolve with time,
# and that it is equal to the one of GRB 090510 (this is true for the first version of the GRB catalog)

    def spectrum(x):
        return (x/1)**(-2.1)

# Integral of the spectrum

    intl,errl=integrate.quad(lambda x: spectrum(x),x1,x2) #GeV


# Interpolation of the flux with time
    flux = interpolate.interp1d(time,lc)


#################################################
# Starting the procedure
#################################################

# Output file
    outtsn="ObsTimes/%s_%s_obstime.txt"%(run[j],id[j])
    
    outfile2n = open(outtsn,"w")
    outfile2n.write('#\n')
    outfile2n.write('#N Pointing\ttstart\ttend\tObs time\n')
    outfile2n.write('#\t\ts\n')


# starting time of the observations
    tstart=talert+tslew # (s)
    
    # Loop over consecutive pointings
    for n in range(1,totaltime):
        
        # increasing the observing time to get the GRB detection
        for m in range(1,totaltime,1):
            t=tstart+m*dt
            obst=t-tstart
            
            fluencen, errorn = integrate.quad(lambda x: flux(x), tstart, t) #ph/cm2/GeV
            averagefluxn=fluencen*intl/obst #ph/cm2/s
            
            if averagefluxn >sensi(obst):
                
                if tstart+obst < totaltime:
                    tend=tstart+obst
                    outfile2n.write("%d\t%d\t%d\t%d\n"%(n,tstart,tend,obst))
                
                tstart=tstart+obst+tslewbtwn
                
                break


        if tstart >=totaltime:
                    break
