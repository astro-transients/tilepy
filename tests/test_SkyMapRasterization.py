"""Regression tests for fixed-order rasterization of multi-order sky maps."""

from types import SimpleNamespace

import healpy as hp
import mhealpy as mh
import numpy as np
import pytest

from tilepy.include.MapManagement.SkyMap import SkyMap


def make_skymap(raw_map):
    obspar = SimpleNamespace(algorithm="2D", minimumProbCutForCatalogue=0.01)
    reader = SimpleNamespace(getMap=lambda map_type: raw_map)
    return SkyMap(obspar, reader)


@pytest.mark.parametrize("nside", [2, 4, 8, 16])
def test_multiorder_probability_rasterization(nside, monkeypatch):
    raw_map = mh.HealpixMap.moc_from_pixels(
        nside=8, pixels=[0, 1, 2], nest=True, density=True
    )
    raw_map.data[:] = np.arange(1, raw_map.npix + 1)
    raw_map.data[:] /= np.sum(raw_map.data * raw_map.pixarea().value)
    original_density = raw_map.data.copy()
    skymap = make_skymap(raw_map)

    expected_density = raw_map.rasterize(nside=nside, scheme="NESTED").data

    def reject_generic_rasterization(*args, **kwargs):
        raise AssertionError("Multi-order maps should use ligo.skymap rasterization")

    monkeypatch.setattr(raw_map, "rasterize", reject_generic_rasterization)
    density = skymap.getMap("prob_density", nside)
    probability = skymap.getMap("prob", nside)

    np.testing.assert_allclose(density, expected_density, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(probability, density * hp.nside2pixarea(nside))
    np.testing.assert_allclose(probability.sum(), 1.0)
    np.testing.assert_array_equal(raw_map.data, original_density)
    assert skymap.getMap("prob_density", nside) is density


@pytest.mark.parametrize("scheme", ["NESTED", "RING"])
def test_fixed_order_map_keeps_existing_rasterization(scheme):
    raw_map = mh.HealpixMap(
        data=np.arange(1, hp.nside2npix(4) + 1, dtype=float),
        scheme=scheme,
        density=True,
    )
    skymap = make_skymap(raw_map)

    np.testing.assert_array_equal(
        skymap.getMap("prob_density", 8),
        raw_map.rasterize(nside=8, scheme=scheme).data,
    )
