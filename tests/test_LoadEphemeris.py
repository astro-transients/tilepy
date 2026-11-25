"""Tests of the search order and fallbacks of LoadEphemeris.

All tests are network-free: opening and downloading the ephemeris are replaced
by stand-ins, except for the bundled ephemeris of skyfield-data.
"""

import datetime
import logging
import sys
from pathlib import Path

import pytest

from tilepy.include import PointingTools
from tilepy.include.PointingTools import LoadEphemeris, Observer


class FakeKernel:
    """Stand-in for skyfield's SpiceKernel: only remembers the opened path."""

    def __init__(self, path):
        self.path = Path(path)
        OPENED.append(self.path)


OPENED = []


@pytest.fixture
def dirs(tmp_path, monkeypatch):
    """Empty working and home directories, no cache, no environment variable."""
    cwd = tmp_path / "cwd"
    home = tmp_path / "home"
    cwd.mkdir()
    home.mkdir()
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.delenv("TILEPY_EPHEMERIS", raising=False)
    monkeypatch.setattr(PointingTools, "_EPHEMERIS_CACHE", {})
    monkeypatch.setattr(PointingTools, "SpiceKernel", FakeKernel)
    OPENED.clear()
    return cwd, home


def forbid_download(monkeypatch):
    def loader(directory):
        raise AssertionError("no download expected")

    monkeypatch.setattr(PointingTools, "Loader", loader)


def fail_download(monkeypatch):
    """Downloads fail as they do without network access."""

    class FailingLoader:
        def __init__(self, directory):
            self.directory = directory

        def __call__(self, filename):
            raise OSError(f"cannot download {filename}")

    monkeypatch.setattr(PointingTools, "Loader", FailingLoader)


def test_environment_variable(dirs, tmp_path, monkeypatch):
    cwd, _home = dirs
    forbid_download(monkeypatch)
    kernel = tmp_path / "my_kernel.bsp"
    kernel.touch()
    (cwd / "de440s.bsp").touch()  # ignored: the variable takes precedence
    monkeypatch.setenv("TILEPY_EPHEMERIS", str(kernel))

    assert LoadEphemeris().path == kernel


def test_environment_variable_without_file(dirs, tmp_path, monkeypatch):
    monkeypatch.setenv("TILEPY_EPHEMERIS", str(tmp_path / "missing.bsp"))

    with pytest.raises(FileNotFoundError, match="TILEPY_EPHEMERIS"):
        LoadEphemeris()


def test_working_directory_before_home(dirs, monkeypatch):
    cwd, home = dirs
    forbid_download(monkeypatch)
    (cwd / "de440s.bsp").touch()
    (home / ".skyfield").mkdir()
    (home / ".skyfield" / "de440s.bsp").touch()

    assert LoadEphemeris().path == cwd / "de440s.bsp"


def test_home_directory(dirs, monkeypatch):
    _cwd, home = dirs
    forbid_download(monkeypatch)
    (home / ".skyfield").mkdir()
    (home / ".skyfield" / "de440s.bsp").touch()

    assert LoadEphemeris().path == home / ".skyfield" / "de440s.bsp"


def test_download_into_home_directory(dirs, monkeypatch):
    _cwd, home = dirs
    calls = []

    class FakeLoader:
        def __init__(self, directory):
            self.directory = directory

        def __call__(self, filename):
            calls.append((Path(self.directory), filename))
            return "downloaded"

    monkeypatch.setattr(PointingTools, "Loader", FakeLoader)

    assert LoadEphemeris() == "downloaded"
    assert calls == [(home / ".skyfield", "de440s.bsp")]


def test_opened_once(dirs, monkeypatch):
    cwd, _home = dirs
    forbid_download(monkeypatch)
    (cwd / "de440s.bsp").touch()

    assert LoadEphemeris() is LoadEphemeris()
    assert len(OPENED) == 1


def test_no_download_and_no_bundled_ephemeris(dirs, monkeypatch):
    fail_download(monkeypatch)
    monkeypatch.setitem(sys.modules, "skyfield_data", None)  # not installed

    with pytest.raises(OSError) as excinfo:
        LoadEphemeris()

    message = str(excinfo.value)
    assert ".skyfield" in message
    assert "TILEPY_EPHEMERIS" in message
    assert "skyfield-data" in message


def test_bundled_ephemeris_when_download_fails(dirs, monkeypatch, caplog):
    pytest.importorskip("skyfield_data")
    from skyfield.jpllib import SpiceKernel

    fail_download(monkeypatch)
    monkeypatch.setattr(PointingTools, "SpiceKernel", SpiceKernel)

    with caplog.at_level(logging.WARNING):
        ephemeris = LoadEphemeris()

    assert "DE421" in caplog.text
    for body in ("sun", "moon", "earth"):
        assert ephemeris[body] is not None


def test_observer_uses_load_ephemeris(monkeypatch):
    sentinel = object()
    monkeypatch.setattr(PointingTools, "LoadEphemeris", lambda: sentinel)

    observer = Observer(
        longitude=-17.89,
        latitude=28.76,
        elevation=2200,
        run_duration=datetime.timedelta(minutes=20),
        minimal_run_duration=datetime.timedelta(minutes=10),
        max_sun_altitude=-18,
        max_moon_altitude=-0.5,
        max_moon_phase=-1.0,
    )

    assert observer.eph is sentinel
