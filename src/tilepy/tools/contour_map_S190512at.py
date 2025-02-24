# Imports
from math import log

import astropy.coordinates as co
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
from astropy.utils.data import download_file
from matplotlib.patches import Circle
from mocpy import MOC, World2ScreenMPL

fig = plt.figure(111, figsize=(5, 5))
percentage = 0.9

# Download and read sky map Initial S190512at.
url = "https://gracedb.ligo.org/api/superevents/S190512at/files/LALInference.fits.gz,0"
filename = download_file(url, cache=True)


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

radecs = co.SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))


# creating an astropy.table with RA[deg] and DEC[deg] ipix positions
contour_ipix = Table(
    [ra, dec], names=("RA[deg]", "DEC[deg]"), meta={"ipix": "ipix table"}
)

# setting MOC order
moc_order = int(log(nside, 2))

moc1 = MOC.from_lonlat(radecs.ra, radecs.dec, max_norder=moc_order)

#########################################################################################################

# Do the same for a 50% interval
percentage = 0.5

# Download and read sky map Update S190728q.
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

radecs = co.SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))


# creating an astropy.table with RA[deg] and DEC[deg] ipix positions
contour_ipix = Table(
    [ra, dec], names=("RA[deg]", "DEC[deg]"), meta={"ipix": "ipix table"}
)

# setting MOC order
moc_order = int(log(nside, 2))

moc2 = MOC.from_lonlat(radecs.ra, radecs.dec, max_norder=moc_order)

#########################################################################################################

# Download and read sky map Update S190512at.
percentage = 0.9
url = "https://gracedb.ligo.org/api/superevents/S190512at/files/GW190512_180714_PublicationSamples_flattened.fits.gz,0"
filename = download_file(url, cache=True)
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

radecs = co.SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))


# creating an astropy.table with RA[deg] and DEC[deg] ipix positions
contour_ipix = Table(
    [ra, dec], names=("RA[deg]", "DEC[deg]"), meta={"ipix": "ipix table"}
)

# setting MOC order
moc_order = int(log(nside, 2))

moc3 = MOC.from_lonlat(radecs.ra, radecs.dec, max_norder=moc_order)

#########################################################################################################
# redoo the same for a 50% contour
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

radecs = co.SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))


# creating an astropy.table with RA[deg] and DEC[deg] ipix positions
contour_ipix = Table(
    [ra, dec], names=("RA[deg]", "DEC[deg]"), meta={"ipix": "ipix table"}
)

# setting MOC order
moc_order = int(log(nside, 2))

moc4 = MOC.from_lonlat(radecs.ra, radecs.dec, max_norder=moc_order)


#########################################################################################################
"""
#add the pointings and cinstruct the figure

RA_OBS = [312.8027, 313.09, 317.109, 314.385, 312.979, 318.691, 316.143]
DEC_OBS = [7.0304, 8.16, 15.018, 10.807, 5.679, 17.192, 12.941]
RADIUS_OBS = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]
#RADECS_OBS = co.SkyCoord(RA_OBS, DEC_OBS, frame='icrs', unit=(u.deg, u.deg))
"""

# construct the figure
with World2ScreenMPL(
    fig,
    fov=20 * u.deg,
    center=SkyCoord(252.5, -27.5, unit="deg", frame="icrs"),
    coordsys="icrs",
    rotation=Angle(0, u.degree),
    projection="AIT",
) as wcs:

    # add the moc contours
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    # moc1.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green")
    moc1.border(ax=ax, wcs=wcs, alpha=0.5, color="blue", antialiased=True)

    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    # moc2.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red")
    moc2.border(ax=ax, wcs=wcs, alpha=0.5, color="blue", antialiased=True)

    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    # moc2.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red")
    moc3.border(ax=ax, wcs=wcs, alpha=0.5, color="red", antialiased=True)

    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    # moc2.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red")
    moc4.border(ax=ax, wcs=wcs, alpha=0.5, color="red", antialiased=True)

    ###############################
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    # add the schedule FoVs
    c = Circle(
        (251.72, -25.28),
        1.5,
        edgecolor="grey",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=0.8,
    )
    ax.add_patch(c)
    c = Circle(
        (248.91, -27.95),
        1.5,
        edgecolor="grey",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=0.8,
    )
    ax.add_patch(c)
    c = Circle(
        (253.12, -26.61),
        1.5,
        edgecolor="grey",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=0.8,
    )
    ax.add_patch(c)
    c = Circle(
        (255.94, -29.31),
        1.5,
        edgecolor="grey",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=0.8,
    )
    ax.add_patch(c)

    # add the observed FoVs
    c = Circle(
        (250.313, -26.61),
        1.5,
        edgecolor="black",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=15,
    )
    ax.add_patch(c)
    ax.text(250.313, -26.61, "1", transform=ax.get_transform("icrs"))
    c = Circle(
        (251.719, -27.953),
        1.5,
        edgecolor="black",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=15,
    )
    ax.add_patch(c)
    ax.text(251.719, -27.953, "2", transform=ax.get_transform("icrs"))
    c = Circle(
        (248.906, -25.283),
        1.5,
        edgecolor="black",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=15,
    )
    ax.add_patch(c)
    ax.text(248.906, -25.283, "3", transform=ax.get_transform("icrs"))
    c = Circle(
        (254.531, -27.953),
        1.5,
        edgecolor="black",
        facecolor="none",
        transform=ax.get_transform("icrs"),
        alpha=15,
    )
    ax.add_patch(c)
    ax.text(254.531, -27.953, "4", transform=ax.get_transform("icrs"))
    ###########################

plt.xlabel("Right ascension")
plt.ylabel("Declination")
plt.title("H.E.S.S. coverage of S190512at GW event")
plt.grid(color="black", linestyle="dotted")
plt.show()
