"""Tests of the observation windows computed for the ranking of the pointings.

All tests are network-free.
"""

import re

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time

from tilepy.include import RankingObservationTimes
from tilepy.include.CampaignDefinition import ObservationParameters
from tilepy.include.PointingTools import getdate
from tilepy.include.RankingObservationTimes import (
    EvolutionPlot,
    GetObservationPeriod,
    VisibilityWindow,
)

CONFIG = """[observatory]
name = LST
lat = 28.761944
lon = -17.89
height = 2200
[visibility]
sunDown = {sunDown}
moonDown = {moonDown}
moonGrey = 65
moonPhase = {moonPhase}
minMoonSourceSeparation = 30
maxMoonSourceSeparation = 145
[operations]
maxZenith = 70
FOV = 2.5
maxRuns = 4
maxNights = 1
duration = 20
minDuration = 10
useGreytime = {useGreytime}
"""

# Night of the 16th of June 2024 at La Palma: the Moon (70% illuminated) sets
# around 02:15 UTC, the astronomical twilight starts around 04:40 UTC.
START = getdate("2024-06-16 00:30:00")
SOURCE = SkyCoord(240.0, 20.0, unit="deg")
# Position still high in the sky at dawn
SOURCE_AT_DAWN = SkyCoord(300.0, 30.0, unit="deg")


@pytest.fixture
def make_obspar(tmp_path):
    def _make(sunDown=-18, moonDown=-0.5, moonPhase=60, useGreytime=False):
        cfg = tmp_path / "config.ini"
        cfg.write_text(
            CONFIG.format(
                sunDown=sunDown,
                moonDown=moonDown,
                moonPhase=moonPhase,
                useGreytime=useGreytime,
            )
        )
        obspar = ObservationParameters()
        obspar.from_configfile(cfg)
        return obspar

    return _make


def window(obspar, source=SOURCE):
    """Start and end of the observation period, as ISO strings."""
    period, _altitudes = GetObservationPeriod(START, source, obspar)
    start, end = period.split("-->")
    return start, end


def test_default_window(make_obspar):
    start, end = window(make_obspar())

    # From the Moon set to the start of the astronomical twilight
    assert "2024-06-16 02:10:00" < start < "2024-06-16 02:20:00"
    assert "2024-06-16 04:30:00" < end < "2024-06-16 04:45:00"


def test_sun_altitude_from_configuration(make_obspar):
    _start, end = window(make_obspar(), SOURCE_AT_DAWN)
    _start, end_nautical = window(make_obspar(sunDown=-12), SOURCE_AT_DAWN)

    # Observations are allowed further into the twilight
    assert end_nautical > end
    assert "2024-06-16 05:00:00" < end_nautical < "2024-06-16 05:20:00"


def test_moon_altitude_from_configuration(make_obspar):
    start, end = window(make_obspar())
    start_moon, end_moon = window(make_obspar(moonDown=90))

    # The Moon does not prevent observations anymore
    assert start_moon < start
    assert start_moon.startswith("2024-06-16 00:30")
    assert end_moon == end


def test_grey_time_extends_the_window(make_obspar):
    start_dark, end_dark = window(make_obspar())

    # Moon phase (70%) below the limit: observations start before the Moon set
    start, end = window(make_obspar(moonPhase=80, useGreytime=True))
    assert start < start_dark
    assert start.startswith("2024-06-16 00:30")
    assert end == end_dark

    # Moon phase above the limit: dark time only
    assert window(make_obspar(moonPhase=60, useGreytime=True)) == (
        start_dark,
        end_dark,
    )


def test_visibility_window_computes_each_period_once(make_obspar, monkeypatch):
    calls = []

    def fake_period(inputtime0, msource, obspar):
        calls.append((inputtime0, msource))
        return "start-->end", msource.dec

    monkeypatch.setattr(RankingObservationTimes, "GetObservationPeriod", fake_period)

    pointings = Table(
        {
            "Pointing": [0, 1],
            "Time": ["2024-06-16 02:30:00", "2024-06-16 02:50:00"],
            "RA[deg]": [240.0, 243.0],
            "DEC[deg]": [20.0, 22.0],
        }
    )
    obspar = make_obspar()
    result = VisibilityWindow(START, pointings, obspar, ".")

    assert len(calls) == len(pointings)
    # The periods start one run duration before the first pointing
    assert all(time == getdate("2024-06-16 02:10:00") for time, _ in calls)
    assert list(result["Observation window"]) == ["start-->end"] * 2
    assert u.allclose(result["Array of zenith angles"], [20.0, 22.0] * u.deg)
    assert "Zenith angles in steps" in result.colnames


def test_evolution_plot_time_axis_is_continuous(make_obspar, tmp_path):
    # The last pointing is not consecutive to the others
    pointings = Table(
        {
            "Pointing": [0, 1, 2],
            "Time": ["2024-06-16 02:30", "2024-06-16 02:50", "2024-06-16 04:10"],
            "RA[deg]": [240.0, 243.0, 300.0],
            "DEC[deg]": [20.0, 22.0, 30.0],
            "Pgw": [0.02, 0.05, 0.01],
            "Pgal": [0.03, 0.01, 0.02],
        }
    )
    obspar = make_obspar()
    fig = EvolutionPlot(pointings, tmp_path, "LST", obspar)

    assert (tmp_path / "AltitudevsTime_LST.png").exists()
    curves = [line for line in fig.axes[0].get_lines() if line.get_label()[0] != "_"]
    assert len(curves) == len(pointings)

    for line in curves:
        times = Time(list(line.get_xdata()))
        # From the start of the first pointing to the end of the last one,
        # without jumps over the gap between the pointings
        assert abs(times[0] - Time("2024-06-16 02:30")) < 1 * u.s
        assert abs(times[-1] - Time("2024-06-16 04:30")) < 1 * u.s
        assert np.diff(times.mjd).max() * u.day <= 1 * u.min

        # The curve is the one of the pointing given in its label
        ra, dec = (
            float(re.search(rf"{name}:(\S+)", line.get_label()).group(1))
            for name in ("ra", "dec")
        )
        expected = (
            SkyCoord(ra, dec, unit="deg")
            .transform_to(AltAz(obstime=times, location=obspar.location))
            .alt.deg
        )
        assert np.allclose(line.get_ydata(), expected)


def test_evolution_plot_distinguishes_many_pointings(make_obspar, tmp_path):
    n = 20
    pointings = Table(
        {
            "Pointing": range(n),
            "Time": [
                f"2024-06-16 {2 + i // 6:02d}:{10 * (i % 6):02d}" for i in range(n)
            ],
            "RA[deg]": np.linspace(230.0, 260.0, n),
            "DEC[deg]": np.linspace(10.0, 40.0, n),
            "Pgw": np.linspace(0.05, 0.01, n),
            "Pgal": np.full(n, 0.01),
        }
    )
    fig = EvolutionPlot(pointings, tmp_path, "LST", make_obspar())
    curves = [line for line in fig.axes[0].get_lines() if line.get_label()[0] != "_"]

    # The 8 most probable pointings have distinct colors and solid lines, the
    # next ones reuse the colors with another line style
    colors = [line.get_color() for line in curves]
    assert len(set(colors[:8])) == 8
    assert all(line.get_linestyle() == "-" for line in curves[:8])
    assert colors[8:16] == colors[:8]
    assert all(line.get_linestyle() != "-" for line in curves[8:])
