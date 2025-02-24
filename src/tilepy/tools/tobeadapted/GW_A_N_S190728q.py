#############################################################
# Author: Halim Ashkar                                      #
# Modified by: Monica Seglar-Arroyo                         #
# WARNING: This file is included as a part of the package   #
# but doesnt run in the same env, since it needs updated    #
# versions of matplotlib, etc                               #
# check the installation in the webside of ligo.skymap.plot #
#############################################################

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import colors
from matplotlib.patches import Circle

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
RA_GRB = 217
DEC_GRB = 34
# url = 'https://gracedb.ligo.org/api/superevents/S190728q/files/GW190728_064510_PublicationSamples_flattened.fits.gz,0'
center = SkyCoord(RA_GRB, DEC_GRB, unit="deg", frame="icrs")
center_str = "%fd %fd" % (center.ra.deg, center.dec.deg)
# center = SkyCoord(232.418, 28.071, unit='deg', frame='icrs')
# filename = download_file(url, cache=True)
# filename = '/Users/hashkar/Desktop/GW_HESS_PUBLICATION/GW_O2_O3_MAPS/maps/O3/GW190512_180714_PublicationSamples_flattened.fits.gz,0'
# filename = '/Users/mseglar/Documents/GitLab/dataset/skymaps/run0378_MergerID000234_skymap.fits.gz'
filename = "/Users/hashkar/Desktop/lst_gwfollowup/bayestar.fits-2.gz,1"


# start preparing figure and inset figure
fig = plt.figure(figsize=(9, 6))

ax = plt.axes([0.1, 0.1, 0.4, 0.4], projection="astro globe", center=center_str)

ax_inset = plt.axes(
    [0.4, 0.2, 0.45, 0.45], projection="astro zoom", center=center_str, radius="10 deg"
)

for key in ["ra", "dec"]:
    ax_inset.coords[key].set_ticklabel_visible(True)
    ax_inset.coords[key].set_ticks_visible(True)


ax.grid()
ax.mark_inset_axes(ax_inset)
ax.connect_inset_axes(ax_inset, loc="upper left")
ax.connect_inset_axes(ax_inset, loc="lower left")
# ax.connect_inset_axes(ax_inset, loc = 'center')


# ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()
# ax_inset.compass(0.9, 0.1, 0.2)

ax.imshow_hpx(filename, cmap="cylon")


# BNS coordinates
# Source = SkyCoord(InputList['RA'][j], InputList['Dec'][j], frame='icrs', unit=(u.deg, u.deg))
source = Circle(
    (RA_GRB, DEC_GRB),
    0.1,
    edgecolor="green",
    facecolor="none",
    transform=ax_inset.get_transform("icrs"),
)
# source = Circle((232.418, 32.271), 0.1, edgecolor='green', facecolor='none', transform=ax_inset.get_transform('icrs'))
ax_inset.add_patch(source)

# neutrino fro S190728q
# c = Circle((312.87, 5.85), 4.81, edgecolor='green', facecolor='none', transform=ax_inset.get_transform('icrs'), linestyle='dashed')
# ax_inset.add_patch(c)

"""
c = Circle((312.8027, 7.0304), 1.5, edgecolor='black', facecolor='none', transform=ax_inset.get_transform('icrs'),alpha=15)#this is the first pos that was actually observed
ax_inset.add_patch(c)
ax_inset.text(312.8027, 7.0304, "1\n26 min\n54 deg",transform=ax_inset.get_transform('icrs'))
"""
# read coordinates from file

tpointingFile = "/Users/hashkar/Desktop/lst_gwfollowup/output/MS221117n/PGallinFoV/SuggestedPointings_GWOptimisation.txt"
# tpointingFile = '/Users/mseglar/Documents/GitLab/lst_gwfollowup/output/bn180720598/PGWinFoV/RankingObservationTimes_Complete.txt'
time = []
time1, time2, ra, dec, pgw, pgal, Round = np.genfromtxt(
    tpointingFile,
    usecols=(0, 1, 2, 3, 4, 5, 6),
    skip_header=1,
    delimiter=" ",
    unpack=True,
    dtype="str",
)  # ra, dec in degrees
print(time1, time2)
for i in range(len(time1)):
    time.append((time1[i] + " " + time2[i]).split('"')[1])
ra = np.atleast_1d(ra)
dec = np.atleast_1d(dec)
ra = ra.astype(float)
dec = dec.astype(float)
pgw = pgw.astype(float)
pgal = pgal.astype(float)
coordinates = SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))
print(ra, dec, pgw, pgal)

for i in range(0, len(ra)):
    print(ra[i])
    c = Circle(
        (ra[i], dec[i]),
        2.0,
        edgecolor="black",
        facecolor="none",
        transform=ax_inset.get_transform("icrs"),
        alpha=15,
    )
    ax_inset.add_patch(c)
    ax_inset.text(
        ra[i] - 2.5,
        dec[i] - 1,
        "%d\n%s \n %d%% %d deg " % (i, time[i], 100 * pgw[i], pgal[i]),
        transform=ax_inset.get_transform("icrs"),
        color="k",
        rotation=-15,
        fontsize=8,
    )


pos = ax_inset.imshow_hpx(filename, cmap="cylon")


cbar = fig.colorbar(pos, ax=ax, location="left", fraction=0.046, pad=0.05)
cbar.set_label("Probability", color="black", fontsize=7)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("bottom", size="1%", pad=0.05)

"""
ax_inset.plot(
              center.ra.deg, center.dec.deg,
              transform=ax_inset.get_transform('world'),
              marker=ligo.skymap.plot.reticle(),
              markersize=30,
              markeredgewidth=3)
"""
# plt.savefig("S190728q_OBS_PrettyMap.png", bbox_inches = 'tight', pad_inches = 0.1)
plt.show()
