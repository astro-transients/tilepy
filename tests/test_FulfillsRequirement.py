"""Tests of the galaxy selection in FulfillsRequirement (3D scheduling),
including the optional user-defined constraints (issue #147).
"""

from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.table import Table

from tilepy.include.MaskingTools import FulfillsRequirement
from tilepy.include.PointingTools import getdate

OBSTIME = getdate("2024-06-16 02:12:17")  # naive datetime, as in the CLI


def make_obspar(strategy, constraints=None):
    return SimpleNamespace(
        maxZenith=60,
        FOV=2.0,
        zenithWeighting=0.75,
        strategy=strategy,
        userConstraints=constraints,
        location=EarthLocation(lat=28.76 * u.deg, lon=-17.89 * u.deg, height=2200),
        obs_name="LST",
    )


@pytest.fixture
def galaxies():
    """Galaxies sorted by decreasing probability, as in the scheduling loops."""
    return Table(
        {
            "RAJ2000": [10.0, 20.0, 30.0, 40.0, 50.0],
            "DEJ2000": [5.0, 15.0, 25.0, 35.0, 45.0],
            "Alt": [80.0, 70.0, 45.0, 35.0, 29.0],
            "dp_dV": [0.5, 0.2, 0.15, 0.1, 0.05],
            "dp_dV_FOV": [0.5, 0.2, 0.15, 0.1, 0.05],
        }
    )


def exclude_ra(*ras):
    def constraint(time, coords, obspar):
        return ~np.isin(np.round(coords.ra.deg, 3), ras)

    return constraint


@pytest.mark.parametrize("strategy", ["targeted", "integrated"])
def test_time_without_constraints_matches_legacy(galaxies, strategy):
    obspar = make_obspar(strategy)
    mask0, minz0 = FulfillsRequirement(galaxies, obspar, UsePix=False)
    mask1, minz1 = FulfillsRequirement(galaxies, obspar, UsePix=False, time=OBSTIME)

    assert mask0.dtype == bool
    np.testing.assert_array_equal(mask0, mask1)
    assert minz0 == minz1
    # The most probable galaxy is selected first
    assert galaxies[mask0]["RAJ2000"][0] == 10.0


def test_constraints_ignored_without_time(galaxies):
    obspar = make_obspar("targeted", [exclude_ra(10.0)])
    mask, _ = FulfillsRequirement(galaxies, obspar, UsePix=False)
    assert galaxies[mask]["RAJ2000"][0] == 10.0


@pytest.mark.parametrize("strategy", ["targeted", "integrated"])
def test_vetoed_galaxy_is_not_selected(galaxies, strategy):
    obspar = make_obspar(strategy, [exclude_ra(10.0)])
    mask, _ = FulfillsRequirement(galaxies, obspar, UsePix=False, time=OBSTIME)

    selected = galaxies[mask]
    assert 10.0 not in selected["RAJ2000"]
    # The next most probable galaxy becomes the pointing
    assert selected["RAJ2000"][0] == 20.0


def test_zenith_band_is_chosen_among_allowed_galaxies(galaxies):
    # Without the two galaxies at small zenith angle, the remaining ones are
    # only found at large zenith angle: they must still be selected.
    obspar = make_obspar("targeted", [exclude_ra(10.0, 20.0)])
    mask, _ = FulfillsRequirement(galaxies, obspar, UsePix=False, time=OBSTIME)

    assert galaxies[mask]["RAJ2000"][0] == 30.0


def test_use_pix_applies_constraints(galaxies):
    obspar = make_obspar("targeted", [exclude_ra(10.0)])
    mask, _ = FulfillsRequirement(galaxies, obspar, UsePix=True, time=OBSTIME)
    assert not mask[0]
    assert mask[1]


@pytest.mark.parametrize("constraints", [None, [lambda t, c, o: False]])
def test_empty_selection_is_a_boolean_mask(galaxies, constraints):
    """No galaxy left (zenith cut and/or constraints): empty table, no crash."""
    if constraints is None:
        galaxies["Alt"] = 10.0  # all below the zenith cut
    obspar = make_obspar("targeted", constraints)
    mask, _ = FulfillsRequirement(galaxies, obspar, UsePix=False, time=OBSTIME)

    assert mask.dtype == bool
    assert mask.shape == (len(galaxies),)
    assert not mask.any()
    assert len(galaxies[mask]) == 0
