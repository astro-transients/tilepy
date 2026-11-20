"""Regression tests for polygonal fields of view centred near a celestial pole.

``Tools.get_regular_polygon_vertices`` used to add the radius offsets directly
to RA and Dec. For a centre within one radius of a pole, a vertex got
``|dec| > 90 deg`` and ``SkyCoord`` raised "Latitude angle(s) must be within
-90 deg <= angle <= 90 deg". With ``strategy = integrated``, ``ModifyCatalogue``
called ``ComputePGalinFOV`` for the visible galaxies, so a single galaxy near a
pole made the whole schedule fail (for example S190425z seen from the H.E.S.S.
site, with a GLADE galaxy at dec -89.09 deg and a 2 deg FoV). The vertices are
now built on the sphere, so a polygon around a pole simply extends over it.
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from tilepy.include.PointingTools import ComputePGalinFOV, ModifyCatalogue, Tools

RA_CENTER = 123.4
RADIUS = 2.0


def vertices_to_skycoord(vertices):
    """Convert (N, 3) unit vectors to a SkyCoord."""
    return SkyCoord(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        representation_type="cartesian",
    )


@pytest.mark.parametrize("dec_center", [-89.09, 89.5, -90.0, 90.0])
@pytest.mark.parametrize("n_sides, rotation", [(4, 0), (6, 15)])
def test_vertices_near_a_pole_are_on_the_sphere(dec_center, n_sides, rotation):
    vertices = Tools.get_regular_polygon_vertices(
        RA_CENTER, dec_center, RADIUS, n_sides, rotation
    )

    assert vertices.shape == (n_sides, 3)
    assert np.all(np.isfinite(vertices))
    np.testing.assert_allclose(np.linalg.norm(vertices, axis=1), 1.0)

    # Every vertex is one radius away from the centre, also when the polygon
    # extends over the pole.
    center = SkyCoord(RA_CENTER * u.deg, dec_center * u.deg)
    separations = vertices_to_skycoord(vertices).separation(center).deg
    np.testing.assert_allclose(separations, RADIUS, atol=1e-9)


def galaxy_catalogue():
    """Four galaxies: one near the south pole, a neighbour on the same side,
    one across the pole (1.41 deg from the first one), and one far away."""
    return Table(
        {
            "RAJ2000": [RA_CENTER, RA_CENTER + 0.5, RA_CENTER + 180, 40.0],
            "DEJ2000": [-89.09, -88.8, -89.5, 10.0],
            "dp_dV": [0.4, 0.3, 0.2, 0.1],
        }
    )


@pytest.mark.parametrize("n_sides", [0, 4, 6])
def test_pgal_in_fov_near_the_pole(n_sides):
    cat = galaxy_catalogue()

    pgal, remaining = ComputePGalinFOV(
        None, cat, cat[0], RADIUS, 1.0, n_sides=n_sides, UsePix=False
    )

    # The galaxy across the pole is inside the FoV as well.
    assert pgal == pytest.approx(0.9)
    assert list(remaining["DEJ2000"]) == [10.0]


@pytest.mark.parametrize("n_sides", [0, 4])
def test_integrated_strategy_catalogue_near_the_pole(n_sides):
    # ModifyCatalogue is what TilingDetermination.PGalinFoV calls for
    # strategy == "integrated".
    cat = galaxy_catalogue()

    tcat = ModifyCatalogue(None, cat, RADIUS, 1.0, n_sides)

    pgal_by_dec = dict(zip(tcat["DEJ2000"], tcat["dp_dV_FOV"]))
    assert pgal_by_dec == pytest.approx(
        {-89.09: 0.9, -88.8: 0.9, -89.5: 0.9, 10.0: 0.1}
    )
    assert tcat["DEJ2000"][-1] == pytest.approx(10.0)
