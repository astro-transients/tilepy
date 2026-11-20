"""Regression tests for polygonal fields of view centred near a celestial pole.

``Tools.get_regular_polygon_vertices`` adds the radius offsets directly to RA
and Dec. For a centre within one radius of a pole, a vertex got
``|dec| > 90 deg`` and ``SkyCoord`` raised "Latitude angle(s) must be within
-90 deg <= angle <= 90 deg". With ``strategy = integrated``, ``ModifyCatalogue``
calls ``ComputePGalinFOV`` for the visible galaxies, so a single galaxy near a
pole made the whole schedule fail (for example S190425z seen from the H.E.S.S.
site, with a GLADE galaxy at dec -89.09 deg and a 2 deg FoV).
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from tilepy.include.PointingTools import ComputePGalinFOV, ModifyCatalogue, Tools

RA_CENTER = 123.4
RADIUS = 2.0


def vertices_to_radec(vertices):
    """Convert (N, 3) unit vectors to RA and Dec in degrees."""
    coords = SkyCoord(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        representation_type="cartesian",
    ).represent_as("spherical")
    return coords.lon.deg, coords.lat.deg


def flat_offset_vertices(ra_center, dec_center, radius, n_sides, rotation):
    """The unclipped RA/Dec vertices, as the helper builds them."""
    angles = np.radians(np.linspace(0, 360, n_sides, endpoint=False) + rotation)
    ra = (ra_center + radius * np.cos(angles) + 360) % 360
    dec = dec_center + radius * np.sin(angles)
    return ra, dec


@pytest.mark.parametrize("dec_center", [-89.09, 89.5, -90.0, 90.0])
@pytest.mark.parametrize("n_sides, rotation", [(4, 0), (6, 15)])
def test_vertices_beyond_a_pole_are_clipped_on_their_meridian(
    dec_center, n_sides, rotation
):
    vertices = Tools.get_regular_polygon_vertices(
        RA_CENTER, dec_center, RADIUS, n_sides, rotation
    )

    assert vertices.shape == (n_sides, 3)
    assert np.all(np.isfinite(vertices))
    np.testing.assert_allclose(np.linalg.norm(vertices, axis=1), 1.0)

    ra, dec = vertices_to_radec(vertices)
    ra_flat, dec_flat = flat_offset_vertices(
        RA_CENTER, dec_center, RADIUS, n_sides, rotation
    )
    beyond = np.abs(dec_flat) > 90 - 1e-6
    assert beyond.any()

    # A vertex beyond the pole stays on its own meridian, just short of the
    # pole, instead of being reflected to RA + 180 deg.
    np.testing.assert_allclose(ra, ra_flat, atol=1e-6)
    np.testing.assert_allclose(np.abs(dec[beyond]), 90, atol=1e-5)
    assert np.all(np.sign(dec[beyond]) == np.sign(dec_flat[beyond]))
    # The other vertices keep their offsets.
    np.testing.assert_allclose(dec[~beyond], dec_flat[~beyond], atol=1e-9)


@pytest.mark.parametrize(
    "ra_center, dec_center, n_sides, rotation",
    [
        (0.0, 0.0, 4, 0),
        (359.5, 45.0, 4, 0),
        (200.0, -60.0, 6, 15),
        # Touches dec -89.5 deg, still clear of the pole.
        (10.0, -87.5, 4, 0),
        (300.0, 87.9, 8, 22.5),
    ],
)
def test_vertices_clear_of_the_poles_are_unchanged(
    ra_center, dec_center, n_sides, rotation
):
    vertices = Tools.get_regular_polygon_vertices(
        ra_center, dec_center, RADIUS, n_sides, rotation
    )

    ra_flat, dec_flat = flat_offset_vertices(
        ra_center, dec_center, RADIUS, n_sides, rotation
    )
    expected = SkyCoord(ra=ra_flat * u.deg, dec=dec_flat * u.deg).cartesian
    np.testing.assert_allclose(vertices, expected.xyz.value.T, rtol=0, atol=1e-12)


def galaxy_catalogue():
    """Three galaxies: one near the south pole, a neighbour inside the same
    2 deg FoV, and one far away."""
    return Table(
        {
            "RAJ2000": [RA_CENTER, RA_CENTER + 0.5, 40.0],
            "DEJ2000": [-89.09, -88.5, 10.0],
            "dp_dV": [0.5, 0.3, 0.2],
        }
    )


def test_pgal_in_polygonal_fov_near_the_pole():
    cat = galaxy_catalogue()

    pgal, remaining = ComputePGalinFOV(
        None, cat, cat[0], RADIUS, 1.0, n_sides=4, UsePix=False
    )

    assert pgal == pytest.approx(0.8)
    assert list(remaining["DEJ2000"]) == [10.0]


def test_integrated_strategy_catalogue_near_the_pole():
    # ModifyCatalogue is what TilingDetermination.PGalinFoV calls for
    # strategy == "integrated"; it passes the HEALPix nside as n_sides.
    cat = galaxy_catalogue()

    tcat = ModifyCatalogue(None, cat, RADIUS, 1.0, 512)

    pgal_by_dec = dict(zip(tcat["DEJ2000"], tcat["dp_dV_FOV"]))
    assert pgal_by_dec == pytest.approx({-89.09: 0.8, -88.5: 0.8, 10.0: 0.2})
    assert tcat["DEJ2000"][-1] == 10.0
