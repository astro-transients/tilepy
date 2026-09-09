"""Regression tests for the subtraction of previously observed pointings.

See issue #137: a pointing lying inside the ``percentageMOC`` region was
refused ("... as it is outside of the 90.0% area") because the containment
test compared the distance to the region pixel centres against a tolerance
derived from ``HRnside`` while the region is defined at ``reducedNside``. The
pointing was then scheduled a second time.
"""

import astropy.units as u
import healpy as hp
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from tilepy.include.CampaignDefinition import (
    ObservationParameters,
    set_gaussian_source,
)
from tilepy.include.MapManagement.MapReader import create_map_reader
from tilepy.include.MapManagement.SkyMap import SkyMap
from tilepy.include.PointingTools import SubtractPointings, SubtractPointings2D

SOURCE_RA, SOURCE_DEC = 240.0, 20.0

CONFIG = """[observatory]
name = CTAO-N
lat = 28.75
lon = -17.5
height = 2200
[visibility]
sunDown = -18
moonDown = -0.5
moonGrey = 65
moonPhase = 60
minMoonSourceSeparation = 30
maxMoonSourceSeparation = 145
[operations]
maxZenith = 70
FOV = 2.0
maxRuns = 5
maxNights = 1
duration = 15
minDuration = 10
useGreytime = False
[tiling]
minimumProbCutforCatalogue = 0.01
minProbcut = 0.002
distCut = 500
doPlot = False
secondRound = False
zenithWeighting = 0.75
percentageMOC = 0.90
reducedNside = 128
HRnside = 512
mangrove = False
algorithm = 3D
strategy = integrated
doRank = False
countPrevious = False
countSubtractedPointingsOutside = False
[general]
downloadMaxRetry = 3
downloadWaitPeriodRetry = 20
"""


@pytest.fixture
def setup(tmp_path):
    """A synthetic Gaussian map (no network) and its percentageMOC pixel list."""
    cfg = tmp_path / "config.ini"
    cfg.write_text(CONFIG)

    obspar = ObservationParameters()
    set_gaussian_source(obspar, ra=SOURCE_RA, dec=SOURCE_DEC, sigma=5.0)
    obspar.from_configfile(cfg)

    skymap = SkyMap(obspar, create_map_reader(obspar))
    regionPix = skymap.getPixIdArea(obspar.percentageMOC, obspar.reducedNside)

    return obspar, skymap, regionPix


def write_pointings(path, coord):
    """Write a pointings file in the format produced by tilepy."""
    path.write_text(
        "Time[UTC] RA[deg] DEC[deg] PGW Round ObsName Duration FoV\n"
        f'"2024-06-16 02:12:17" {coord.ra.deg:.4f} {coord.dec.deg:.4f} '
        "0.1 1 CTAO-N 15 2.0\n"
    )
    return path


def offset_from_pixel_centre(regionPix, nside, is_nested, offset_deg=0.25):
    """A coordinate inside the region but away from any region pixel centre.

    ``offset_deg`` is larger than ``nside2resol(HRnside)`` (0.115 deg for
    HRnside = 512), which is the tolerance the buggy containment test used.
    """
    lon, lat = hp.pix2ang(nside, regionPix[0], nest=is_nested, lonlat=True)
    centre = SkyCoord(lon * u.deg, lat * u.deg)
    return centre.directional_offset_by(0.0 * u.deg, offset_deg * u.deg)


def make_galaxies(coord):
    """A minimal galaxy catalogue with one galaxy at ``coord``."""
    return Table(
        {
            "RAJ2000": [coord.ra.deg],
            "DEJ2000": [coord.dec.deg],
            "dp_dV": [1.0],
        }
    )


def test_3D_subtracts_pointing_between_region_pixel_centres(setup, tmp_path):
    """A pointing inside the region must be subtracted (issue #137)."""
    obspar, skymap, regionPix = setup
    coord = offset_from_pixel_centre(regionPix, obspar.reducedNside, skymap.is_nested)
    pointingFile = write_pointings(tmp_path / "Pointings.txt", coord)

    prob = skymap.getMap("prob", obspar.HRnside)
    galaxies = make_galaxies(coord)

    *_, PGW, PGAL, _, pointings_subtracted = SubtractPointings(
        str(pointingFile),
        galaxies,
        [],
        galaxies["dp_dV"].sum(),
        prob,
        skymap.is_nested,
        obspar,
        obspar.HRnside,
        regionPix,
    )

    assert pointings_subtracted == 1
    assert sum(PGW) > 0
    assert sum(PGAL) > 0


def test_2D_subtracts_pointing_between_region_pixel_centres(setup, tmp_path):
    obspar, skymap, regionPix = setup
    coord = offset_from_pixel_centre(regionPix, obspar.reducedNside, skymap.is_nested)
    pointingFile = write_pointings(tmp_path / "Pointings.txt", coord)

    prob = skymap.getMap("prob", obspar.reducedNside)

    _, _, sumPGW, pointings_subtracted = SubtractPointings2D(
        str(pointingFile), prob, skymap.is_nested, obspar, [], [], regionPix
    )

    assert pointings_subtracted == 1
    assert sumPGW > 0


def test_pointing_far_from_region_is_not_subtracted(setup, tmp_path):
    """A pointing whose FoV does not reach the region must still be skipped."""
    obspar, skymap, regionPix = setup
    coord = SkyCoord(SOURCE_RA * u.deg, SOURCE_DEC * u.deg).directional_offset_by(
        0.0 * u.deg, 90.0 * u.deg
    )
    pointingFile = write_pointings(tmp_path / "Pointings.txt", coord)

    prob = skymap.getMap("prob", obspar.reducedNside)

    _, _, sumPGW, pointings_subtracted = SubtractPointings2D(
        str(pointingFile), prob, skymap.is_nested, obspar, [], [], regionPix
    )

    assert pointings_subtracted == 0
    assert sumPGW == 0.0


def test_pointing_centred_outside_but_overlapping_is_subtracted(setup, tmp_path):
    """PGalinFoV centres pointings on galaxies, which may sit just outside the
    contour; their FoV still covers part of the region, so they must count."""
    obspar, skymap, regionPix = setup

    # Walk outwards from the source until just outside the region, then check
    # that the FoV of a pointing there still overlaps it.
    inRegion = np.isin(
        hp.ang2pix(
            obspar.reducedNside,
            SOURCE_RA,
            SOURCE_DEC,
            nest=skymap.is_nested,
            lonlat=True,
        ),
        regionPix,
    )
    assert inRegion, "sanity check: the source must be inside the region"

    source = SkyCoord(SOURCE_RA * u.deg, SOURCE_DEC * u.deg)
    coord = None
    for offset in np.arange(0.5, 30.0, 0.25):
        candidate = source.directional_offset_by(0.0 * u.deg, offset * u.deg)
        pix = hp.ang2pix(
            obspar.reducedNside,
            candidate.ra.deg,
            candidate.dec.deg,
            nest=skymap.is_nested,
            lonlat=True,
        )
        if pix not in regionPix:
            coord = candidate
            break
    assert coord is not None, "could not find a coordinate outside the region"

    pointingFile = write_pointings(tmp_path / "Pointings.txt", coord)
    prob = skymap.getMap("prob", obspar.reducedNside)

    _, _, _, pointings_subtracted = SubtractPointings2D(
        str(pointingFile), prob, skymap.is_nested, obspar, [], [], regionPix
    )

    assert pointings_subtracted == 1
