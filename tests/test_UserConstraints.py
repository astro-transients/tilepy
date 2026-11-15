"""Tests of the optional user-defined constraints layer (issue #147).

All tests are network-free.
"""

import datetime
import os
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

from tilepy.include.CampaignDefinition import ObservationParameters
from tilepy.include.PointingTools import getdate
from tilepy.include.UserConstraints import (
    EvaluateUserConstraints,
    LoadUserConstraints,
)

LOCATION = EarthLocation(lat=28.76 * u.deg, lon=-17.89 * u.deg, height=2200 * u.m)
OBSTIME = getdate("2024-06-16 02:12:17")  # naive datetime, as in the CLI
COORDS = SkyCoord(np.linspace(0.0, 350.0, 36), np.linspace(-60.0, 80.0, 36), unit="deg")


def make_obspar(constraints=None):
    return SimpleNamespace(
        userConstraints=constraints, location=LOCATION, obs_name="LST"
    )


def north_only(time, coords, obspar):
    return coords.dec.deg > 0


def ra_below_180(time, coords, obspar):
    return coords.ra.deg < 180


# --------------------------------------------------------------------------
# EvaluateUserConstraints
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "obspar",
    [
        SimpleNamespace(),  # no attribute at all
        make_obspar(None),
        make_obspar([]),
        ObservationParameters(),
    ],
)
def test_no_constraints_returns_none(obspar):
    assert EvaluateUserConstraints(OBSTIME, COORDS, obspar) is None


def test_callable_constraint():
    mask = EvaluateUserConstraints(OBSTIME, COORDS, make_obspar([north_only]))
    assert mask.dtype == bool
    np.testing.assert_array_equal(mask, COORDS.dec.deg > 0)


def test_single_constraint_without_list():
    mask = EvaluateUserConstraints(OBSTIME, COORDS, make_obspar(north_only))
    np.testing.assert_array_equal(mask, COORDS.dec.deg > 0)


def test_constraints_are_combined_with_and():
    obspar = make_obspar([north_only, ra_below_180])
    mask = EvaluateUserConstraints(OBSTIME, COORDS, obspar)
    np.testing.assert_array_equal(mask, (COORDS.dec.deg > 0) & (COORDS.ra.deg < 180))


def test_constraint_receives_utc_time_coords_and_obspar():
    seen = {}

    def spy(time, coords, obspar):
        seen["time"], seen["coords"], seen["obspar"] = time, coords, obspar
        return True

    obspar = make_obspar([spy])
    utc = datetime.timezone.utc
    plus2 = datetime.timezone(datetime.timedelta(hours=2))
    for time in (
        OBSTIME,  # naive, interpreted as UTC
        OBSTIME.replace(tzinfo=utc),
        (OBSTIME + datetime.timedelta(hours=2)).replace(tzinfo=plus2),
        Time(OBSTIME, scale="utc"),
        Time(OBSTIME, scale="utc").tt,
    ):
        EvaluateUserConstraints(time, COORDS, obspar)
        assert isinstance(seen["time"], Time)
        assert seen["time"].isscalar
        assert seen["time"].scale == "utc"
        assert abs(seen["time"] - Time(OBSTIME, scale="utc")) < 1 * u.ms
        assert len(seen["coords"]) == len(COORDS)
        assert seen["obspar"] is obspar


def test_scalar_result_is_broadcast():
    assert EvaluateUserConstraints(
        OBSTIME, COORDS, make_obspar([lambda t, c, o: True])
    ).all()
    assert not EvaluateUserConstraints(
        OBSTIME, COORDS, make_obspar([lambda t, c, o: False])
    ).any()


def test_scores_are_converted_to_booleans():
    scores = np.linspace(0.0, 1.0, len(COORDS))
    scores[5] = np.nan
    mask = EvaluateUserConstraints(
        OBSTIME, COORDS, make_obspar([lambda t, c, o: scores])
    )
    expected = scores > 0
    expected[5] = False
    np.testing.assert_array_equal(mask, expected)


def test_empty_coordinates():
    def fail(time, coords, obspar):
        raise AssertionError("must not be called without candidates")

    mask = EvaluateUserConstraints(OBSTIME, COORDS[:0], make_obspar([fail]))
    assert mask.shape == (0,)


def test_wrong_size_raises():
    obspar = make_obspar([lambda t, c, o: np.ones(3, dtype=bool)])
    with pytest.raises(ValueError, match="3 values for 36 sky positions"):
        EvaluateUserConstraints(OBSTIME, COORDS, obspar)


def test_class_instead_of_instance_raises():
    class MyConstraint:
        def compute_constraint(self, times, observer, targets):
            return True

    with pytest.raises(TypeError, match="pass an instance"):
        EvaluateUserConstraints(OBSTIME, COORDS, make_obspar([MyConstraint]))


def test_not_a_constraint_raises():
    with pytest.raises(TypeError, match="neither"):
        EvaluateUserConstraints(OBSTIME, COORDS, make_obspar(["north_only"]))


# --------------------------------------------------------------------------
# astroplan constraints
# --------------------------------------------------------------------------


def test_astroplan_builtin_constraint():
    astroplan = pytest.importorskip("astroplan")

    obspar = make_obspar([astroplan.AltitudeConstraint(min=30 * u.deg)])
    mask = EvaluateUserConstraints(OBSTIME, COORDS, obspar)

    altaz = COORDS.transform_to(
        AltAz(obstime=Time(OBSTIME, scale="utc"), location=LOCATION)
    )
    assert mask.shape == (len(COORDS),)
    np.testing.assert_array_equal(mask, altaz.alt.deg >= 30)
    assert 0 < mask.sum() < len(COORDS)


def test_astroplan_time_only_constraint():
    astroplan = pytest.importorskip("astroplan")

    night = astroplan.AtNightConstraint.twilight_astronomical()
    # 02:12 UTC is night-time at La Palma, 14:12 UTC is not
    assert EvaluateUserConstraints(OBSTIME, COORDS, make_obspar([night])).all()
    noon = OBSTIME + datetime.timedelta(hours=12)
    assert not EvaluateUserConstraints(noon, COORDS, make_obspar([night])).any()


def test_astroplan_user_defined_constraint():
    astroplan = pytest.importorskip("astroplan")

    class MoonAngleConstraint(astroplan.Constraint):
        """User-defined constraint following the astroplan tutorial."""

        def __init__(self, min_separation):
            self.min_separation = min_separation

        def compute_constraint(self, times, observer, targets):
            altaz = observer.altaz(times, targets)
            moon = observer.moon_altaz(times)
            return altaz.separation(moon) > self.min_separation

    constraint = MoonAngleConstraint(60 * u.deg)
    mask = EvaluateUserConstraints(
        OBSTIME, COORDS, make_obspar([constraint, north_only])
    )

    observer = astroplan.Observer(location=LOCATION)
    obstime = Time(OBSTIME, scale="utc")
    separation = observer.altaz(obstime, COORDS).separation(
        observer.moon_altaz(obstime)
    )
    expected = (separation > 60 * u.deg) & (COORDS.dec.deg > 0)
    np.testing.assert_array_equal(mask, expected)
    assert 0 < mask.sum() < len(COORDS)


# --------------------------------------------------------------------------
# LoadUserConstraints and configuration file
# --------------------------------------------------------------------------

USER_FILE = """
def north_only(time, coords, obspar):
    return coords.dec.deg > 0

def south_only(time, coords, obspar):
    return coords.dec.deg < 0

both = [north_only, south_only]

class Holder:
    cut = staticmethod(north_only)

class NotAnInstance:
    def compute_constraint(self, times, observer, targets):
        return True
"""


@pytest.fixture
def user_file(tmp_path):
    path = tmp_path / "cuts" / "my_cuts.py"
    path.parent.mkdir()
    path.write_text(USER_FILE)
    return path


def test_load_from_absolute_path(user_file):
    constraints = LoadUserConstraints(f"{user_file}:north_only")
    assert [c.__name__ for c in constraints] == ["north_only"]


def test_load_relative_to_base_dir(user_file, tmp_path, monkeypatch):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    constraints = LoadUserConstraints("cuts/my_cuts.py:north_only", baseDir=tmp_path)
    assert [c.__name__ for c in constraints] == ["north_only"]


def test_load_list_and_several_entries(user_file):
    spec = f"{user_file}:both,\n  {user_file}:Holder.cut\n"
    constraints = LoadUserConstraints(spec)
    assert [c.__name__ for c in constraints] == [
        "north_only",
        "south_only",
        "north_only",
    ]


def test_load_from_module():
    constraints = LoadUserConstraints("operator:truth")
    assert len(constraints) == 1
    assert callable(constraints[0])


def test_load_empty_spec():
    assert LoadUserConstraints("") == []
    assert LoadUserConstraints(None) == []


@pytest.mark.parametrize(
    "spec, error",
    [
        ("my_cuts.py", ValueError),  # no attribute
        ("{file}:", ValueError),
        ("{file}:not-an-identifier", ValueError),
        ("{file}:missing", ValueError),
        ("does_not_exist.py:cut", ValueError),
        ("{file}:NotAnInstance", TypeError),
        ("math:pi", TypeError),  # neither a constraint nor a callable
    ],
)
def test_load_invalid_reference(user_file, spec, error):
    with pytest.raises(error):
        LoadUserConstraints(spec.format(file=user_file), baseDir=user_file.parent)


CONFIG = """[observatory]
name = LST
lat = 28.76
lon = -17.89
height = 2200
[visibility]
sunDown = -18
{extra}
[operations]
maxZenith = 60
"""


def write_config(tmp_path, extra=""):
    cfg = tmp_path / "config.ini"
    cfg.write_text(CONFIG.format(extra=extra))
    return cfg


def test_default_is_no_constraint():
    assert ObservationParameters().userConstraints == []
    assert ObservationParameters(userConstraints=north_only).userConstraints == [
        north_only
    ]


def test_configfile_without_key_keeps_python_constraints(tmp_path):
    cfg = write_config(tmp_path)

    obspar = ObservationParameters()
    obspar.from_configfile(cfg)
    assert obspar.userConstraints == []

    # Both call orders of from_configfile / add_parsed_args keep constraints
    obspar = ObservationParameters(userConstraints=[north_only])
    obspar.from_configfile(cfg)
    obspar.add_parsed_args(None, OBSTIME, None, None, None, None)
    assert obspar.userConstraints == [north_only]

    obspar = ObservationParameters()
    obspar.userConstraints = [north_only]
    obspar.add_parsed_args(None, OBSTIME, None, None, None, None)
    obspar.from_configfile(cfg)
    assert obspar.userConstraints == [north_only]
    assert "User constraints: 1" in str(obspar)


def test_configfile_with_key(user_file, tmp_path, monkeypatch):
    # Relative to the configuration file, whatever the working directory
    cfg = write_config(tmp_path, "userConstraints = cuts/my_cuts.py:both")
    monkeypatch.chdir(os.path.dirname(user_file))

    obspar = ObservationParameters()
    obspar.from_configfile(cfg)
    assert [c.__name__ for c in obspar.userConstraints] == [
        "north_only",
        "south_only",
    ]
    # north_only and south_only together veto everything
    assert not EvaluateUserConstraints(OBSTIME, COORDS, obspar).any()
