import numpy as np
import sys
sys.path.append('/Users/mseglar/Documents/GitHub/GWfollowup/include')
from GWHESSPointingTools import *
from GWCTAPointingTools import *


#InputFileName = 'BNS-GW-new.txt'
#InputList = TableImportCTA(InputFileName)
MinimumProbCutForCatalogue=0.01

i = sys.argv[1]
i = int(i)
InputFileNameLS = 'BNS_LeoSinger.txt'
InputListLS= TableImportCTA_LS(InputFileNameLS)

#GWFile_LS = '/Users/mseglar/Documents/CurrentPhD/HESS/GW/going-the-distance-o2-skymaps/bayestar/288172_bayestar.fits.gz'
GWFile_LS='/Users/mseglar/Documents/CurrentPhD/HESS/GW/going-the-distance-o2-skymaps/bayestar/' + InputListLS['MergerID'][i] +'_bayestar.fits.gz'
#GWFile_LS='/Users/mseglar/Documents/CurrentPhD/HESS/GW/2016_fits/' + InputListLS['MergerID'][i] +'/bayestar.fits.gz'

GalCat='/Users/mseglar/Documents/CurrentPhD/HESS/GW/GLADE/GLADE_2.3clean.txt'
prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(GWFile_LS)

cat_GLADE=LoadGalaxies(GalCat)
Info3D_available=True
#CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat_GLADE, Info3D_available,MinimumProbCutForCatalogue)
ra = cat_GLADE['RAJ2000']
dec = cat_GLADE['DEJ2000']
dist = cat_GLADE['Dist']

# Translate RA,Dec of galaxies into theta,phi angles

theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)

# Get corresponding healpix pixel IDs

npix = len(prob)
nside = hp.npix2nside(npix)
ipix = hp.ang2pix(nside, theta, phi)

# Calculate probability in the space volumes
pixarea = hp.nside2pixarea(nside)
dp_dV = prob[ipix] * distnorm[ipix] * norm(distmu[ipix], distsigma[ipix]).pdf(dist) / pixarea
fig, ax = plt.subplots()
#x=np.full(len(dp_dV),hp.nside2pixarea(nside,degrees=True))
#plt.plot(np.cumsum(x),np.cumsum(np.flipud(np.sort(dp_dV))/dp_dV.sum()),color='orange',label='LS')
plt.plot(np.sort(prob)/prob.sum(),1-np.cumsum((np.sort(prob))/prob.sum()),color='blue',label='2D GtD')
plt.plot(np.sort(dp_dV)/dp_dV.sum(),1-np.cumsum((np.sort(dp_dV))/dp_dV.sum()),color='orange',label='3D GtD')
ax.set_xscale('log')
#plt.hist(dp_dV / dp_dV.sum(), histtype='step', stacked=True, fill=False, cumulative=True, density=True)
#plt.savefig('/Users/mseglar/Documents/GitHub/CTASimulationsGW/3DComparison_plots/LeoConvolution.png')
print('========================Pixel where the source is found - LEO ========================')
#Get the source and the probability value
raSourceLS=InputListLS['RA'][i]
decSourceLS=InputListLS['Dec'][i]
thetaSourceLS = 0.5 * np.pi - np.deg2rad(decSourceLS)
phiSourceLS = np.deg2rad(raSourceLS)
npix = len(prob)
nside = hp.npix2nside(npix)
ipixSourceLS = hp.ang2pix(nside, thetaSourceLS, phiSourceLS)
pixarea = hp.nside2pixarea(nside)
dp_dVipixSourceLS =prob[ipixSourceLS] * distnorm[ipixSourceLS] * norm(distmu[ipixSourceLS], distsigma[ipixSourceLS]).pdf(InputListLS['Distance'][i]) / pixarea
cumcut=dp_dV< dp_dVipixSourceLS
cdfvalue=np.sum(dp_dV[cumcut])/np.sum(dp_dV)

plt.plot([prob[ipixSourceLS]/np.sum(prob)],[1-(np.sum(prob[prob<prob[ipixSourceLS]])/np.sum(prob))],'+',markersize=10,color='blue')
plt.plot([dp_dVipixSourceLS/dp_dV.sum()],[1-cdfvalue],'+',markersize=10,color='orange')

print('The LS-source is found in: PGW',1-(np.sum(prob[prob<prob[ipixSourceLS]])/np.sum(prob)))
print('The LS-source is found in: PGAL',1-cdfvalue)

prob2DS_ls=1-(np.sum(prob[prob<prob[ipixSourceLS]])/np.sum(prob))
prob3DS_ls=1-cdfvalue
print('================================================================================')
InputFileName = 'BNS-GW-new.txt'
InputList = TableImportCTA(InputFileName)

GWFile_B="/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/skymaps/" + InputList['run'][i] + '_' + InputList['MergerID'][i] + "_skymap.fits.gz"
GalCat_B='/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/GalaxyCatalogs/' + InputList['run'][i] +'_galaxylist.txt'

cat_B=LoadGalaxiesSimulation(GalCat_B)
prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(GWFile_B)
#GiveProbToGalaxy(prob,cat_B,InputList['Distance'][j],InputList['DistMax'][j],InputList['DistMin'][j],MinimumProbCutForCatalogue)
dist = cat_B['Dist']
mask1 = dist < InputList['DistMax'][i]
mask2 = dist > InputList['DistMin'][i]
subcat = cat_B[mask1 & mask2]
# print(len(subcat),len(cat),Edistance_max,Edistance_min)
ra = subcat['RAJ2000']
dec = subcat['DEJ2000']
# Translate RA,Dec of galaxies into theta,phi angles
theta = 0.5 * np.pi - np.deg2rad(dec)
phi = np.deg2rad(ra)

# Get corresponding healpix pixel IDs

npix = len(prob)
nside = hp.npix2nside(npix)
ipix = hp.ang2pix(nside, theta, phi)

# Calculate probability in the space volumes
pixarea = hp.nside2pixarea(nside)
dp_dV = prob[ipix] / pixarea


print('========================Pixel where the source is found - BARBARA ========================')
#Get the source and the probability value
raSourceB=InputList['RA'][i]
decSourceB=InputList['Dec'][i]
thetaSourceB = 0.5 * np.pi - np.deg2rad(decSourceB)
phiSourceB = np.deg2rad(raSourceB)
npix = len(prob)
nside = hp.npix2nside(npix)
ipixSourceB = hp.ang2pix(nside, thetaSourceB, phiSourceB)
pixarea = hp.nside2pixarea(nside)
dp_dVipixSourceB = prob[ipixSourceB] / pixarea
cumcut=dp_dV< dp_dVipixSourceB
cdfvalue=np.sum(dp_dV[cumcut])/np.sum(dp_dV)
print('The LS-source is found in: PGW',1-(np.sum(prob[prob<prob[ipixSourceB]])/np.sum(prob)))
print('The LS-source is found in: PGAL',1-cdfvalue)
print('======================== Comparison  ========================')

#x=np.full(len(dp_dV),hp.nside2pixarea(nside,degrees=True))
#plt.plot(np.cumsum(x),np.cumsum(np.flipud(np.sort(dp_dV))/dp_dV.sum()),color='blue',label='GWCosmos')

plt.plot(np.sort(prob)/prob.sum(),1-np.cumsum((np.sort(prob))/prob.sum()),color='green',label='2D GWCosmos')

plt.plot(np.sort(dp_dV)/dp_dV.sum(),1-np.cumsum((np.sort(dp_dV))/dp_dV.sum()),color='red',label='3D GWCosmos')
plt.plot([prob[ipixSourceB]/np.sum(prob)],[1-(np.sum(prob[prob<prob[ipixSourceB]])/np.sum(prob))],'+', markersize=10,color='green')
plt.plot([dp_dVipixSourceB/dp_dV.sum()],[1-cdfvalue],'+',markersize=10,color='red')
ax.legend()
ax.set_xlim((0.0000001,1))
ax.set_xlabel('Probabability')
ax.set_ylabel('1-c.d.f.')
ax.set_xscale('log')
plt.grid()
#plt.hist(dp_dV / dp_dV.sum(), histtype='step', stacked=True, fill=False, cumulative=True, density=True)
plt.savefig('/Users/mseglar/Documents/GitHub/CTASimulationsGW/3DComparison_plots/PorbabilityCDF.png')

prob2DS_b=1-(np.sum(prob[prob<prob[ipixSourceB]])/np.sum(prob))
prob3DS_b=1-cdfvalue

outfilename='3DComparison_plots/Entry'+str(i)+'.txt'
#print(outfilename)
f=open(outfilename,'w')
f.write(str(prob2DS_ls)+' '+str(prob3DS_ls)+' '+str(prob2DS_b)+' '+str(prob3DS_b)+'\n')
