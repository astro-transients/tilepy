"""Tests for the geometry of fields of view on the sphere.

Polygonal FoVs used to be built and tested in flat RA/Dec: the vertices were
offset in RA without the cos(dec) factor, and ``ComputePGalinFOV`` tested
containment with a ``matplotlib.path.Path`` in RA/Dec, which breaks across
RA 0/360. ``ComputePGalinFOV`` also always used 4 vertices, and
``ModifyCatalogue`` passed the HEALPix nside as the number of sides, so the
integrated strategy ranked its galaxies with a square instead of the circle
that it then counts and subtracts.
"""

import astropy.units as u
import healpy as hp
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from tilepy.include.MaskingTools import GetBestGridPos2D
from tilepy.include.PointingTools import ComputePGalinFOV, ModifyCatalogue, Tools

RADIUS = 2.0

CENTRES = [
    (10.0, 0.0),
    (359.5, -30.0),
    (0.2, 45.0),
    (200.0, 60.0),
    (123.4, -89.09),
    (300.0, 90.0),
]


def unit_vectors(coords):
    return coords.cartesian.xyz.value.T


def offset(ra, dec, position_angle, separation):
    """SkyCoord at the given position angle (deg, east of north) and
    separation (deg) from (ra, dec)."""
    return SkyCoord(ra * u.deg, dec * u.deg).directional_offset_by(
        position_angle * u.deg, separation * u.deg
    )


@pytest.mark.parametrize("ra_center, dec_center", CENTRES)
@pytest.mark.parametrize(
    "n_sides, rotation", [(3, 0), (4, 0), (4, 30), (6, 15), (8, 0)]
)
def test_vertices_match_directional_offset_by(ra_center, dec_center, n_sides, rotation):
    vertices = Tools.get_regular_polygon_vertices(
        ra_center, dec_center, RADIUS, n_sides, rotation
    )

    # Vertex k is at the angle 360 * k / n_sides + rotation from east towards
    # north, which is a position angle of 90 deg minus that angle.
    angles = np.linspace(0, 360, n_sides, endpoint=False) + rotation
    expected = offset(ra_center, dec_center, 90 - angles, RADIUS)
    np.testing.assert_allclose(vertices, unit_vectors(expected), atol=1e-12)


def test_vertices_on_the_equator_keep_their_layout():
    # Without rotation the first vertex is due east, then north, west, south.
    vertices = Tools.get_regular_polygon_vertices(10.0, 0.0, RADIUS, 4, 0)

    expected = SkyCoord([12.0, 10.0, 8.0, 10.0] * u.deg, [0.0, 2.0, 0.0, -2.0] * u.deg)
    np.testing.assert_allclose(vertices, unit_vectors(expected), atol=1e-12)


@pytest.mark.parametrize("dec_center", [0.0, 60.0, -75.0, 89.0])
def test_polygon_has_the_same_size_on_the_sky_at_any_declination(dec_center):
    vertices = Tools.get_regular_polygon_vertices(50.0, dec_center, RADIUS, 4, 0)
    coords = SkyCoord(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        representation_type="cartesian",
    )

    center = SkyCoord(50.0 * u.deg, dec_center * u.deg)
    np.testing.assert_allclose(coords.separation(center).deg, RADIUS, atol=1e-9)
    sides = coords.separation(np.roll(coords, -1)).deg
    np.testing.assert_allclose(sides, sides[0], atol=1e-9)


@pytest.mark.parametrize(
    "n_sides, radius", [(1, RADIUS), (2, RADIUS), (4, 0.0), (4, 90.0), (4, -1.0)]
)
def test_invalid_polygons_raise(n_sides, radius):
    with pytest.raises(ValueError):
        Tools.get_regular_polygon_vertices(10.0, 0.0, radius, n_sides, 0)


@pytest.mark.parametrize("ra_center, dec_center", CENTRES)
@pytest.mark.parametrize("n_sides", [3, 4, 6])
@pytest.mark.parametrize("reverse", [False, True])
def test_points_in_spherical_polygon(ra_center, dec_center, n_sides, reverse):
    vertices = Tools.get_regular_polygon_vertices(
        ra_center, dec_center, RADIUS, n_sides, 0
    )
    if reverse:
        vertices = vertices[::-1]
    angles = np.linspace(0, 360, n_sides, endpoint=False)
    inradius = np.degrees(
        np.arctan(np.tan(np.radians(RADIUS)) * np.cos(np.pi / n_sides))
    )
    edge_angles = angles + 180 / n_sides

    center = SkyCoord(ra_center * u.deg, dec_center * u.deg)
    inside = [
        center,
        offset(ra_center, dec_center, 90 - angles, 0.999 * RADIUS),
        offset(ra_center, dec_center, 90 - edge_angles, 0.999 * inradius),
    ]
    on_the_edge = [
        offset(ra_center, dec_center, 90 - angles, RADIUS),
        offset(ra_center, dec_center, 90 - edge_angles, inradius),
    ]
    outside = [
        offset(ra_center, dec_center, 90 - angles, 1.001 * RADIUS),
        offset(ra_center, dec_center, 90 - edge_angles, 1.001 * inradius),
        # The antipode of the centre, and the other side of the sky.
        SkyCoord((ra_center + 180) * u.deg, -dec_center * u.deg),
        offset(ra_center, dec_center, 0, 90.0),
    ]

    for coords in inside + on_the_edge:
        mask = Tools.points_in_spherical_polygon(vertices, unit_vectors(coords))
        assert np.all(mask)
    for coords in outside:
        mask = Tools.points_in_spherical_polygon(vertices, unit_vectors(coords))
        assert not np.any(mask)


@pytest.mark.parametrize("ra_center, dec_center", CENTRES)
@pytest.mark.parametrize("n_sides", [4, 6])
def test_containment_agrees_with_healpy_query_polygon(ra_center, dec_center, n_sides):
    # Without inclusive, healpy returns the pixels whose centres lie inside the
    # polygon, which is the same test as the one used for galaxies.
    nside = 256
    vertices = Tools.get_regular_polygon_vertices(
        ra_center, dec_center, RADIUS, n_sides, 0
    )
    healpy_pix = hp.query_polygon(nside, vertices)

    candidates = hp.query_disc(
        nside,
        unit_vectors(SkyCoord(ra_center * u.deg, dec_center * u.deg)),
        np.radians(RADIUS + 1),
    )
    centres = np.column_stack(hp.pix2vec(nside, candidates))
    mask = Tools.points_in_spherical_polygon(vertices, centres)

    assert len(healpy_pix) > 0
    np.testing.assert_array_equal(np.sort(candidates[mask]), np.sort(healpy_pix))


def catalogue(coords, dp_dV):
    return Table({"RAJ2000": coords.ra.deg, "DEJ2000": coords.dec.deg, "dp_dV": dp_dV})


@pytest.mark.parametrize("n_sides", [0, 4])
def test_pgal_in_fov_across_ra_zero(n_sides):
    # The FoV centred at RA 359.5 deg used to contain the galaxy at RA 180 deg
    # and to miss its own centre and the galaxy at RA 0.5 deg.
    coords = SkyCoord([359.5, 0.5, 358.0, 180.0, 2.5] * u.deg, [-30.0] * 5 * u.deg)
    cat = catalogue(coords, [0.4, 0.3, 0.15, 0.1, 0.05])

    pgal, remaining = ComputePGalinFOV(
        None, cat, cat[0], RADIUS, 1.0, n_sides=n_sides, UsePix=False
    )

    assert pgal == pytest.approx(0.85)
    np.testing.assert_allclose(remaining["RAJ2000"], [180.0, 2.5])


@pytest.mark.parametrize("n_sides", [0, 4])
def test_pgal_in_fov_has_the_right_size_at_high_declination(n_sides):
    # 3 deg of RA at dec 60 deg is 1.5 deg on the sky: inside a 2 deg FoV.
    # With the flat RA/Dec polygon, the FoV was only 1 deg wide in that
    # direction.
    coords = SkyCoord([200.0, 203.0, 205.0] * u.deg, [60.0] * 3 * u.deg)
    cat = catalogue(coords, [0.5, 0.3, 0.2])

    pgal, remaining = ComputePGalinFOV(
        None, cat, cat[0], RADIUS, 1.0, n_sides=n_sides, UsePix=False
    )

    assert pgal == pytest.approx(0.8)
    np.testing.assert_allclose(remaining["RAJ2000"], [205.0])


def north_east_catalogue():
    """A galaxy and a neighbour 1.6 deg to its north east: inside a 2 deg
    circle or hexagon, outside a square with its vertices at 2 deg (whose
    edge passes 1.41 deg from the centre in that direction)."""
    centre = SkyCoord(100.0 * u.deg, 20.0 * u.deg)
    neighbour = offset(100.0, 20.0, 45.0, 1.6)
    coords = SkyCoord([centre.ra, neighbour.ra], [centre.dec, neighbour.dec])
    return catalogue(coords, [0.6, 0.4])


@pytest.mark.parametrize("n_sides, expected", [(0, 1.0), (4, 0.6), (6, 1.0)])
def test_pgal_in_fov_uses_the_number_of_sides(n_sides, expected):
    cat = north_east_catalogue()

    pgal, _ = ComputePGalinFOV(
        None, cat, cat[0], RADIUS, 1.0, n_sides=n_sides, UsePix=False
    )

    assert pgal == pytest.approx(expected)


def test_integrated_strategy_ranks_with_the_circle_by_default():
    # The integrated strategy counts and subtracts the galaxies within a
    # circle of radius FOV around the chosen pointing
    # (ComputeProbPGALIntegrateFoV), so ModifyCatalogue ranks with that
    # circle unless it is given a number of sides.
    tcat = ModifyCatalogue(None, north_east_catalogue(), RADIUS, 1.0)
    np.testing.assert_allclose(tcat["dp_dV_FOV"], [1.0, 1.0])

    tcat = ModifyCatalogue(None, north_east_catalogue(), RADIUS, 1.0, 4)
    np.testing.assert_allclose(sorted(tcat["dp_dV_FOV"]), [0.4, 0.6])


@pytest.mark.parametrize("dec_center", [-89.09, -90.0, 0.0, 89.5])
def test_grid_positions_2d_with_a_polygonal_fov_near_a_pole(dec_center):
    # The 2D path for polygonal FoVs passes the vertices to hp.query_polygon.
    # Near a pole, the flat RA/Dec polygon was either rejected (vertex beyond
    # the pole) or, once clipped to the pole, covered only 10 to 13% of this
    # 3 deg disc. On the sphere, the best square covers about 40% of it, as
    # it does at the equator.
    reduced_nside, hr_nside = 32, 128
    centre = unit_vectors(SkyCoord(123.4 * u.deg, dec_center * u.deg))
    highres = np.zeros(hp.nside2npix(hr_nside))
    highres[hp.query_disc(hr_nside, centre, np.radians(3))] = 1.0
    highres /= highres.sum()
    prob = hp.ud_grade(highres, reduced_nside, power=-2)
    newpix = hp.query_disc(reduced_nside, centre, np.radians(3))

    pointings = GetBestGridPos2D(
        prob,
        False,
        highres,
        hr_nside,
        reduced_nside,
        newpix,
        RADIUS,
        3,
        [],
        False,
        None,
        4,
        [],
        0.0,
    )

    assert len(pointings) == 3
    assert pointings["PIXFOVPROB"][0] > 0.35
