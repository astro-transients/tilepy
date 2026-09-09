"""Regression tests for GetRegionPixReduced.

The credible region is decoded from the map returned by ``hp.ud_grade``, which
is called with ``order_out="NESTED"``. Deriving the ``nest`` flag from the
scheme of the *input* map instead returned coordinates scattered over the whole
sky for NUNIQ (multi-order LVK maps) and RING maps, which corrupted the grid of
candidate pointings used by the 2D algorithm.
"""

import astropy.units as u
import healpy as hp
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from tilepy.include.PointingTools import GetRegionPixReduced

NSIDE = 128
PERCENTAGE = 0.9
SOURCE = SkyCoord(240.0 * u.deg, 20.0 * u.deg)
SIGMA = 5.0


@pytest.fixture
def prob():
    """A Gaussian localisation, as a NESTED probability map."""
    npix = hp.nside2npix(NSIDE)
    lon, lat = hp.pix2ang(NSIDE, np.arange(npix), nest=True, lonlat=True)
    sep = SkyCoord(lon * u.deg, lat * u.deg).separation(SOURCE).deg
    p = np.exp(-0.5 * (sep / SIGMA) ** 2)
    return p / p.sum()


@pytest.mark.parametrize("scheme", ["NESTED", "NUNIQ"])
def test_region_covers_the_localisation(prob, scheme):
    """The region must sit on the source, not somewhere else on the sky."""
    ra, dec, area = GetRegionPixReduced(prob, PERCENTAGE, NSIDE, scheme)

    separations = SkyCoord(ra * u.deg, dec * u.deg).separation(SOURCE)
    assert separations.max() < 4 * SIGMA * u.deg
    assert area > 0


def test_nested_and_nuniq_agree(prob):
    """NUNIQ maps are nested; both must give the same region."""
    nested = GetRegionPixReduced(prob, PERCENTAGE, NSIDE, "NESTED")
    nuniq = GetRegionPixReduced(prob, PERCENTAGE, NSIDE, "NUNIQ")

    np.testing.assert_array_equal(nested[0], nuniq[0])
    np.testing.assert_array_equal(nested[1], nuniq[1])


def test_ring_input_gives_the_same_region(prob):
    """A RING map of the same localisation must give the same region."""
    nested = GetRegionPixReduced(prob, PERCENTAGE, NSIDE, "NESTED")
    ring = GetRegionPixReduced(hp.reorder(prob, n2r=True), PERCENTAGE, NSIDE, "RING")

    np.testing.assert_allclose(np.sort(nested[0]), np.sort(ring[0]))
    np.testing.assert_allclose(np.sort(nested[1]), np.sort(ring[1]))
