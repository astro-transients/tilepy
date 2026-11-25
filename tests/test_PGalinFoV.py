"""Regression tests for the 3D (galaxy-targeted) schedulers on NumPy >= 2.4.

PGalinFoV and PGalinFoV_NObs stored each pointing as
``float(finalGals["RAJ2000"][:1])``. Converting a 1-element array to a Python
scalar has been deprecated since NumPy 1.25 and raises ``TypeError: only
0-dimensional arrays can be converted to Python scalars`` from NumPy 2.4 on,
so every 3D schedule failed as soon as it found its first pointing.

The tests run both schedulers on a synthetic 3D sky map and a four-galaxy
catalogue, so they need no network access and no galaxy catalogue download.
"""

import datetime

import healpy as hp
import numpy as np
import pytest
import tables

from tilepy.include.CampaignDefinition import ObservationParameters
from tilepy.include.MapManagement.MapReader import create_map_reader
from tilepy.include.MapManagement.SkyMap import SkyMap
from tilepy.include.TilingDetermination import PGalinFoV, PGalinFoV_NObs

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
maxRuns = 4
maxNights = 1
duration = 30
minDuration = 10
useGreytime = False
[tiling]
minimumProbCutforCatalogue = 0.001
minProbcut = 0.03
distCut = 500
doPlot = False
secondRound = {second_round}
zenithWeighting = 0.75
percentageMOC = 0.90
reducedNside = 128
HRnside = 256
mangrove = False
algorithm = 3D
strategy = {strategy}
doRank = False
countPrevious = False
countSubtractedPointingsOutside = False
[general]
downloadMaxRetry = 0
downloadWaitPeriodRetry = 1
"""

# The source is a Gaussian (sigma 1.5 deg) at RA 240, Dec 20, which is dark
# and well above the horizon from CTAO-N at 02:30 UTC on 2024-06-16.
SOURCE_RA, SOURCE_DEC, SOURCE_SIGMA = 240.0, 20.0, 1.5
NSIDE = 128
OBS_TIME = datetime.datetime(2024, 6, 16, 2, 30, 0)

# Three galaxies near the centre of the map and one 3 deg away. A pointing on
# the central group covers PGW ~ 0.56, more than the 0.40 needed for a second
# round. The outlier carries Pgal ~ 0.045, between minProbcut and
# 2 * minProbcut, so with secondRound on the next two pointings are round 2.
# Once round 2 has used up its galaxies, the fourth pointing falls back to
# round 1. With secondRound off, the schedule has one pointing per group.
# Coordinates have more than 4 decimals to check that the schedule rounds the
# coordinates of the selected galaxy.
GALAXIES = [
    (240.123456, 20.012345),
    (240.412121, 20.287654),
    (239.876543, 19.765432),
    (240.054321, 23.012345),
]


@pytest.fixture(scope="module")
def input_files(tmp_path_factory):
    """A 3D HEALPix map and a galaxy catalogue in the GLADE HDF5 layout."""
    directory = tmp_path_factory.mktemp("pgalinfov")

    npix = hp.nside2npix(NSIDE)
    centre = hp.ang2vec(SOURCE_RA, SOURCE_DEC, lonlat=True)
    pixels = np.array(hp.pix2vec(NSIDE, np.arange(npix), nest=True))
    separation = np.degrees(np.arccos(np.clip(centre @ pixels, -1.0, 1.0)))
    prob = np.exp(-0.5 * (separation / SOURCE_SIGMA) ** 2)
    prob /= prob.sum()
    skymap_file = directory / "synthetic_skymap.fits"
    hp.write_map(
        str(skymap_file),
        [prob, np.full(npix, 100.0), np.full(npix, 20.0), np.full(npix, 1e-4)],
        nest=True,
        column_names=["PROB", "DISTMU", "DISTSIGMA", "DISTNORM"],
        column_units=["pix-1", "Mpc", "Mpc", "Mpc-2"],
        extra_header=[("DISTMEAN", 100.0), ("DISTSTD", 20.0), ("OBJECT", "SYNTH")],
        dtype=np.float64,
    )

    rows = np.array(
        [(i, ra, dec, 100.0, 15.0) for i, (ra, dec) in enumerate(GALAXIES)],
        dtype=[
            ("no_GLADE", "i8"),
            ("RA", "f8"),
            ("Dec", "f8"),
            ("d_L", "f8"),
            ("B_mag", "f8"),
        ],
    )
    galaxy_file = directory / "synthetic_galaxies.h5"
    with tables.open_file(str(galaxy_file), mode="w") as h5file:
        h5file.create_table("/", "catalog", obj=rows)

    return directory, skymap_file, galaxy_file


def make_obspar(tmp_path, input_files, strategy, second_round, constraints=None):
    directory, skymap_file, galaxy_file = input_files
    cfg = tmp_path / "config.ini"
    cfg.write_text(CONFIG.format(strategy=strategy, second_round=second_round))

    obspar = ObservationParameters()
    obspar.from_configfile(cfg)
    obspar.add_parsed_args(
        str(skymap_file),
        OBS_TIME,
        str(directory),
        galaxy_file.name,
        str(tmp_path),
        None,
    )
    obspar.userConstraints = list(constraints or [])
    return obspar


def run_pgalinfov(tmp_path, input_files, strategy, second_round, constraints=None):
    obspar = make_obspar(tmp_path, input_files, strategy, second_round, constraints)
    skymap = SkyMap(obspar, create_map_reader(obspar))
    assert skymap.is3D
    pointings, _ = PGalinFoV(
        skymap, "SYNTH", str(input_files[2]), obspar, str(tmp_path)
    )
    return pointings


def run_pgalinfov_nobs(tmp_path, input_files, strategy, second_round, constraints=None):
    obspar = make_obspar(tmp_path, input_files, strategy, second_round, constraints)
    skymap = SkyMap(obspar, create_map_reader(obspar))
    assert skymap.is3D
    pointings, _, _ = PGalinFoV_NObs(
        skymap,
        "SYNTH",
        obspar.obsTime,
        None,
        str(input_files[2]),
        [obspar],
        str(tmp_path),
    )
    return pointings


@pytest.mark.parametrize("second_round", [False, True])
@pytest.mark.parametrize("strategy", ["integrated", "targeted"])
@pytest.mark.parametrize(
    "scheduler",
    [run_pgalinfov, run_pgalinfov_nobs],
    ids=["PGalinFoV", "PGalinFoV_NObs"],
)
def test_schedule_points_at_catalogue_galaxies(
    tmp_path, input_files, scheduler, strategy, second_round
):
    pointings = scheduler(tmp_path, input_files, strategy, second_round)

    assert len(pointings) >= 2

    # Every pointing is centred on a catalogue galaxy, rounded to 4 decimals.
    galaxies = {(float(f"{ra:3.4f}"), float(f"{dec:3.4f}")) for ra, dec in GALAXIES}
    for row in pointings:
        assert (float(row["RA[deg]"]), float(row["DEC[deg]"])) in galaxies

    # Round 2 pointings come from separate code paths, so make sure the
    # second-round cases really reach them.
    if second_round:
        assert 2 in pointings["Round"]
    else:
        assert set(pointings["Round"]) == {1}


# ---------------------------------------------------------------------------
# Optional user-defined constraints (issue #147)
# ---------------------------------------------------------------------------


def outside_central_group(time, coords, obspar):
    """Forbid the three galaxies at the centre of the map."""
    separation = hp.rotator.angdist(
        [SOURCE_RA, SOURCE_DEC], [coords.ra.deg, coords.dec.deg], lonlat=True
    )
    return np.degrees(separation) > 1.0


@pytest.mark.parametrize("strategy", ["integrated", "targeted"])
@pytest.mark.parametrize(
    "scheduler",
    [run_pgalinfov, run_pgalinfov_nobs],
    ids=["PGalinFoV", "PGalinFoV_NObs"],
)
def test_user_constraint_excludes_galaxies(tmp_path, input_files, scheduler, strategy):
    # Without constraint the schedule starts with the central group (see above)
    pointings = scheduler(
        tmp_path, input_files, strategy, False, [outside_central_group]
    )

    # Only the outlying galaxy is left: it is scheduled, the others never
    outlier = tuple(float(f"{x:3.4f}") for x in GALAXIES[3])
    assert len(pointings) >= 1
    for row in pointings:
        assert (float(row["RA[deg]"]), float(row["DEC[deg]"])) == outlier


@pytest.mark.parametrize("strategy", ["integrated", "targeted"])
@pytest.mark.parametrize(
    "scheduler",
    [run_pgalinfov, run_pgalinfov_nobs],
    ids=["PGalinFoV", "PGalinFoV_NObs"],
)
def test_user_constraint_vetoing_everything(tmp_path, input_files, scheduler, strategy):
    pointings = scheduler(
        tmp_path, input_files, strategy, False, [lambda t, c, o: False]
    )

    assert len(pointings) == 0
