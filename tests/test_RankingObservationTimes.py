"""Tests of the observation windows computed for the ranking of the pointings.

All tests are network-free.
"""

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from tilepy.include import RankingObservationTimes
from tilepy.include.CampaignDefinition import ObservationParameters
from tilepy.include.PointingTools import getdate
from tilepy.include.RankingObservationTimes import (
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
