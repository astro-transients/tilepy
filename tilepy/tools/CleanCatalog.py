import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import astropy.coordinates as co
from astropy import units as u
import pandas as pd
import collections

def Clean_catalog_from_NaNs_5var(catalog, name):
    
    diste=catalog['col9']
    dist = np.where(np.isnan(diste), -99, diste)
    print('max(dist),min(dist)',max(dist),min(dist))
    countere=0
    for n,i in enumerate(dist):
        if i==-99:
            countere=countere+1
    print(countere,'objects are nans')
    
    #bmag = np.nan_to_num(bmag)
    bmage=catalog['col12']
    bmag = np.where(np.isnan(bmage), 99, bmage)
    counterbmag=0
    for n,i in enumerate(bmag):
        if i==99:
            counterbmag=counterbmag+1
    print(counterbmag,'objects are nans')
    
    
    ze=catalog['col11']
    z = np.where(np.isnan(ze), 99, bmage)
    counterz=0
    for n,i in enumerate(z):
        if i==99:
            counterbmag=counterbmag+1
    print(counterbmag,'objects are nans')
    
    
    print('bmag',bmag.sum())
    print('maximin',max(bmag),min(bmag))
    plt.hist(dist,50, facecolor='green', alpha=0.75)
    
    for n,i in enumerate(dist):
       if i==0:
         dist[n]=-1
    
    for n,i in enumerate(bmag):
       if i==0:
         bmag[n]=99
    
    for n,i in enumerate(z):
        if i==0:
            z[n]=-1
    
    ascii.write([ra, dec,dist,bmag,z], name+'_noNans.txt',names = ['RAJ2000','DEJ2000','Dist','Bmag','z'])

def Plot_Catalog(tcat):
    plt.hist(tcat['Dist'], 100, facecolor='blue', alpha=0.75, log=True)
    plt.xlabel('Distance(Mpc)')
    plt.ylabel('#')
    plt.show()
    
    distCut = tcat['Dist']
    largedistCUT = distCut < 1000
    shortCat = tcat[largedistCUT]
    print('If we make a cut on distance<1000,lenght=', len(shortCat['Dist']))

    # plt.hist(tcat['Dist'],50, facecolor='blue', alpha=0.75,log=True)
    plt.hist(shortCat['Dist'], 50, facecolor='blue', alpha=0.75, log=True)
    plt.xlabel('Distance(Mpc)')
    plt.ylabel('#')
    plt.show()

def PlotValueandError(value, error, string,stringError):
    fig = plt.figure()
    plt.hist(value, 100, facecolor='blue', alpha=0.75, log=True)
    plt.xlabel(string)
    plt.ylabel('#')
    plt.savefig(string+'.png')
    
    fig = plt.figure()
    plt.hist(error, 100, facecolor='blue', alpha=0.75)
    plt.xlabel(stringError)
    plt.ylabel('#')
    plt.savefig(stringError+'.png')
    
    fig = plt.figure()
    plt.plot(value, error,'.')
    plt.savefig(string+'_'+stringError+'.png')
    
dirName = '/Users/hashkar/Desktop/GWFollowup_lib/gw-follow-up-simulations/tilepy/tools/'

# Cat file should be .cvs
catName = dirName+'GLADE+.cvs'

#Method1 returns a astropy Table
cats = ascii.read(catName,delimiter=' ', format='csv', guess=False,
                 fast_reader={'chunk_size': 100 * 1000000})

# Follow https://docs.astropy.org/en/stable/io/ascii/read.html#chunk-reading to read catalog



#ra_C = []
#dec_C = []
#dist_C = []
#mass_C = []

outcats = []

for cat in cats:
    cutCatalog =  cats['col34'] < 500
    noZero = cats['col34'] != 'null'
    
    #mergerRate = np.array(cat['col39'])
    #EmergerRate = np.array(cat['col40'])
    
    #PlotValueandError(dist,Edist,'Distance [Mpc]','ErrorDistance[Mpc]')
    #PlotValueandError(mass,Emass,'Solar Mass','Error Solar Mass')
    #PlotValueandError(mergerRate,EmergerRate, 'Merger Rate', 'Error Merger Rate')
    #mass = cat['M*']
    #mergerRate = cat['Merger rate']
    
    #Method 2 uses numpy but one needs to handle NaNs
    #ra,dec, dist,Edist, mass,Emass, mergerRate, EmergerRate = np.loadtxt(catName, usecols = (9,10,33,34,36,37,39,40), dtype=str, unpack = True)
    
    # NSBH RANGE IS 300–330 Mpc FOR O4 (aLIGO), so I used 500 Mpc as a reasonable cut
    if np.count_nonzero(cutCatalog and noZero):
        outcats.append(cat[cutCatalog and noZero])
        #ra_C.append(ra[cutCatalog and noZero])
        #dec_C.append(dec[cutCatalog and noZero])
        #dist_C.append(dist[cutCatalog and noZero])
        #mass_C.append(mass[cutCatalog and noZero])
outcat = vstack(outcats)

#ra_out = vstack(ra_C)
#dec_out = vstack(dec_C)
#dist_out = vstack(dist_C)
#mass_out = vstack(mass_C)

print(len(cat),'vs',len(outcat),'when cleaned')
#print('There are',len(dist[(dist)=='null')]),'NaNs in Distance')
#print('There are',len(mass[(np.isnan(mass)==True)]),'NaNs in mass')
#print('There are',len(mergerRate[(np.isnan(mergerRate)==True)]),'NaNs in mergerRate')
# print(cat)


tcat = Table([outcat['col9'],outcat['col10'],outcat['col33'],outcat['col36']], names=('RAJ2000', 'DEJ2000', 'Dist','Mass'))

#print('If we make a cut on nans in distance,lenght=',len(tcat['Dist']))
ascii.write(tcat, dirName+'GLADE+clean.txt',overwrite=True)


#GLADE = Table.read('GLADE_2.2.txt', format='ascii')
