#############################################################
# Author: Halim Ashkar                                      #
# Modified by: Monica Seglar-Arroyo                         #
# WARNING: This file is included as a part of the package   #
# but doesnt run in the same env, since it needs updated    #
# versions of matplotlib, etc                               #
# check the installation in the webside of ligo.skymap.plot #
#############################################################

from math import log

import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib import pyplot as plt
from mocpy import MOC

# Prepare the color map
cdict_coolheat = {
    "red": (
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 1.0, 1.0),
        (0.75, 1.0, 1.0),
        (1.0, 1.0, 1.0),
    ),
    "green": (
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 0.0, 0.0),
        (0.75, 1.0, 1.0),
        (1.0, 1.0, 1.0),
    ),
    "blue": (
        (0.0, 0.0, 0.0),
        (0.25, 1.0, 1.0),
        (0.5, 0.0, 0.0),
        (0.75, 0.0, 0.0),
        (1.0, 1.0, 1.0),
    ),
}
coolheat = colors.LinearSegmentedColormap("coolheat", cdict_coolheat, 1024)

# load the GW map
center = SkyCoord(250.31, -26.61, unit="deg", frame="icrs")
filename = "/Users/hashkar/Desktop/GW_HESS_PUBLICATION/GW_O2_O3_MAPS/maps/O3/GW190512_180714_PublicationSamples_flattened.fits.gz,0"

# load the results map
hdu = fits.open(
    "/Users/hashkar/Desktop/GW_HESS_PUBLICATION/GW_ANALYSIS_MAPS/S190512at_IntULMap_Stereo1Std_E1-10.fits"
)[1]
wcs = WCS(hdu.header)
data = hdu.data
data[data == 0] = np.nan

VMAX = np.nanmax(data)
VMIN = np.nanmin(data)
ticks = np.linspace(VMIN, VMAX, 6)


# start preparing figure and inset figure
fig = plt.figure(figsize=(14, 6))

ax = plt.axes([0.04, 0.05, 0.6, 0.6], projection="astro globe", center=center)

ax_inset = plt.axes([0.35, 0.2, 0.6, 0.6], projection=wcs)

for key in ["ra", "dec"]:
    ax_inset.coords[key].set_ticklabel_visible(True)
    ax_inset.coords[key].set_ticks_visible(True)


ax.grid()
ax.mark_inset_axes(ax_inset)
ax.connect_inset_axes(ax_inset, loc="upper left")
ax.connect_inset_axes(ax_inset, loc="lower left")

ax.imshow_hpx(filename, cmap="cylon")


percentage = 0.5

hpx = hp.read_map(filename, verbose=False)
npix = len(hpx)
nside = hp.npix2nside(npix)

sort = sorted(hpx, reverse=True)
cumsum = np.cumsum(sort)
index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

# finding ipix indices confined in a given percentage
index_hpx = range(0, len(hpx))
hpx_index = np.c_[hpx, index_hpx]

sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
value_contour = sort_2array[0:index]

j = 1
table_ipix_contour = []

for i in range(0, len(value_contour)):
    ipix_contour = int(value_contour[i][j])
    table_ipix_contour.append(ipix_contour)

# from index to polar coordinates
theta, phi = hp.pix2ang(nside, table_ipix_contour)
# converting these to right ascension and declination in degrees
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)

radecs = SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))

# creating an astropy.table with RA[deg] and DEC[deg] ipix positions


contour_ipix = Table(
    [ra, dec], names=("RA[deg]", "DEC[deg]"), meta={"ipix": "ipix table"}
)

# setting MOC order

moc_order = int(log(nside, 2))

# transforming to MOC
moc1 = MOC.from_lonlat(radecs.ra, radecs.dec, max_norder=moc_order)


percentage = 0.9
hpx = hp.read_map(filename, verbose=False)
npix = len(hpx)
nside = hp.npix2nside(npix)

sort = sorted(hpx, reverse=True)
cumsum = np.cumsum(sort)
index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

# finding ipix indices confined in a given percentage
index_hpx = range(0, len(hpx))
hpx_index = np.c_[hpx, index_hpx]

sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
value_contour = sort_2array[0:index]

j = 1
table_ipix_contour = []

for i in range(0, len(value_contour)):
    ipix_contour = int(value_contour[i][j])
    table_ipix_contour.append(ipix_contour)

# from index to polar coordinates
theta, phi = hp.pix2ang(nside, table_ipix_contour)
# converting these to right ascension and declination in degrees
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)

radecs = SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))

# creating an astropy.table with RA[deg] and DEC[deg] ipix positions

contour_ipix = Table(
    [ra, dec], names=("RA[deg]", "DEC[deg]"), meta={"ipix": "ipix table"}
)

# setting MOC order

moc_order = int(log(nside, 2))
moc2 = MOC.from_lonlat(radecs.ra, radecs.dec, max_norder=moc_order)


print(VMIN)

moc1.border(ax=ax_inset, wcs=wcs, alpha=3, color="cyan")
moc2.border(ax=ax_inset, wcs=wcs, alpha=3, color="magenta")

pos = ax_inset.imshow(data, cmap=coolheat, origin="lower", vmin=VMIN, vmax=VMAX)
ax_inset.set_xlim(-0.5, hdu.data.shape[1] - 0.5)
ax_inset.set_ylim(-0.5, hdu.data.shape[0] - 0.5)

plt.grid(color="black", linestyle="--", linewidth=0.5)

lon = ax_inset.coords[0]
lat = ax_inset.coords[1]

ax_inset.coords[0].set_ticks(spacing=5 * u.degree)
ax_inset.coords[1].set_ticks(spacing=5 * u.degree)
lon.set_ticklabel(color="black", size=12)
lat.set_ticklabel(color="black", size=12)


cbar = fig.colorbar(pos, ax=ax_inset, ticks=ticks)
cbar.set_label("Integral Upper Limit ($m^{-2}.s^{-1}$)", color="black", fontsize=9)


"""
ax_inset.plot(
              center.ra.deg, center.dec.deg,
              transform=ax_inset.get_transform('world'),
              marker=ligo.skymap.plot.reticle(),
              markersize=30,
              markeredgewidth=3)
"""
plt.savefig("S190512at_1-10_PrettyMap.png", bbox_inches="tight")
plt.show()
