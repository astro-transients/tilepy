

import sys
sys.path.append('./include')

from PointingTools import *

#Takes in argument  ra , dec (of the g grb/frb), the percentage of the localization area desired and the name of the GW map.
#  It will print if the GRB/FRB is inside this area (True or False) and will display the map with the given point

NewNside = 512

ra = sys.argv[1]
dec = sys.argv[2]
PercentCov = float(sys.argv[3])

filename = sys.argv[4]
#filename = '/Users/hashkar/Desktop/GWfollowup_master_alert/G297595_LALInference.fits'

prob, distmu, distsigma, distnorm, detectors, fits_id, thisDistance, thisDistanceErr = LoadHealpixMap(filename)
npix = len(prob)
nside = hp.npix2nside(npix)

pix_ra1, pix_dec1, area = Get90RegionPixReduced(prob, PercentCov, NewNside)


coordinates=TransformRADec(ra,dec)
ra1= coordinates.ra.rad
dec1 = coordinates.dec.rad


radecs = co.SkyCoord(pix_ra1, pix_dec1, frame='fk5', unit=(u.deg, u.deg))
pix_ra11 = radecs.ra.rad
pix_dec11 = radecs.dec.rad



Igrb = hp.ang2pix(NewNside, ra1, dec1, lonlat=True)
Imap = hp.ang2pix(NewNside, pix_ra11, pix_dec11, lonlat=True)


print(np.isin(Igrb, Imap))

hp.mollview(prob,title="Inside Region")
hp.visufunc.projplot(coordinates.ra, coordinates.dec, 'r.', lonlat=True)
hp.graticule()
plt.show()
