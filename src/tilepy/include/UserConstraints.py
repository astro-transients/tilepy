# Copyright (C) 2016-2025  tilepy developers
# (Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

#
##################################################################################################
#                   Optional layer for user-defined observability constraints                    #
##################################################################################################

import importlib
import importlib.util
import logging
import re
import sys
from pathlib import Path

import numpy as np
from astropy.time import Time

__all__ = [
    "EvaluateUserConstraints",
    "LoadUserConstraints",
]

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

# Modules loaded from a file path, keyed by (resolved path, modification time)
_LOADED_FILES = {}


def _ConstraintName(constraint):
    """Readable name of a constraint for log and error messages."""
    name = getattr(constraint, "__name__", None)
    return name or type(constraint).__name__


def _CheckConstraint(constraint):
    """Raise a TypeError if the object cannot be used as a constraint."""
    if isinstance(constraint, type):
        raise TypeError(
            f"User constraint '{constraint.__name__}' is a class: "
            "pass an instance instead, e.g. "
            f"{constraint.__name__}(...)."
        )
    if not hasattr(constraint, "compute_constraint") and not callable(constraint):
        raise TypeError(
            f"User constraint {constraint!r} is neither an astroplan-style "
            "Constraint (with a compute_constraint method) nor a callable "
            "f(time, coords, obspar)."
        )


def _ToUTCTime(time):
    """Convert a (naive or timezone-aware) datetime or a Time to a UTC Time."""
    if isinstance(time, Time):
        return time.utc
    return Time(time, scale="utc")


def _ToBoolMask(result, nCoords, constraint):
    """Convert the output of a constraint to a boolean array of length nCoords."""
    result = np.asarray(result)
    if result.dtype != bool:
        # Scores: a candidate is allowed if its score is positive; NaN vetoes
        result = np.nan_to_num(result.astype(float), nan=0.0) > 0
    if result.size == 1:
        return np.full(nCoords, bool(result.reshape(-1)[0]))
    if result.size != nCoords:
        raise ValueError(
            f"User constraint '{_ConstraintName(constraint)}' returned "
            f"{result.size} values for {nCoords} sky positions."
        )
    return result.reshape(nCoords)


def _BuildAstroplanObserver(obspar):
    """astroplan Observer at the location of the observatory."""
    try:
        from astroplan import Observer
    except ImportError as err:  # pragma: no cover
        raise ImportError(
            "astroplan is required to evaluate astroplan-style user constraints."
        ) from err
    return Observer(location=obspar.location, name=str(obspar.obs_name))


def EvaluateUserConstraints(time, coords, obspar):
    """
    Evaluate the user-defined constraints for a set of candidate pointings.

    This is an optional layer on top of the built-in visibility cuts. The
    constraints are taken from ``obspar.userConstraints`` and each of them can
    be either

    - an `astroplan <https://astroplan.readthedocs.io>`_ ``Constraint``
      instance (built-in or user-defined, i.e. any object implementing
      ``compute_constraint(times, observer, targets)``), or
    - a callable ``f(time, coords, obspar)`` returning an array of booleans.

    A candidate is kept only if all constraints are fulfilled.

    Parameters
    ----------
    time : datetime.datetime or astropy.time.Time
        Time of the observation. Naive datetimes are interpreted as UTC.
    coords : astropy.coordinates.SkyCoord
        Positions of the candidate pointings (centres of the field of view).
    obspar : ObservationParameters
        Observation parameters. Only ``userConstraints``, ``location`` and
        ``obs_name`` are used here; the full object is handed to callables.

    Returns
    -------
    mask : numpy.ndarray of bool or None
        ``True`` for the positions fulfilling all constraints. ``None`` if no
        user constraint is defined, in which case nothing is evaluated.

    Notes
    -----
    Constraints receive the time as a scalar `~astropy.time.Time` in UTC.
    Non-boolean outputs (e.g. astroplan scores in [0, 1]) are interpreted as
    allowed when larger than zero.
    """
    constraints = getattr(obspar, "userConstraints", None)
    if not constraints:
        return None
    if not isinstance(constraints, (list, tuple)):
        constraints = [constraints]

    coords = coords.reshape(-1)
    nCoords = coords.size
    mask = np.ones(nCoords, dtype=bool)
    if nCoords == 0:
        return mask

    obsTime = _ToUTCTime(time)
    # A new astroplan observer for each evaluation: astroplan caches the
    # alt/az of all targets on the observer and never releases them.
    observer = None
    for constraint in constraints:
        _CheckConstraint(constraint)
        if hasattr(constraint, "compute_constraint"):
            if observer is None:
                observer = _BuildAstroplanObserver(obspar)
            result = constraint(observer, coords, times=obsTime)
        else:
            result = constraint(obsTime, coords, obspar)
        mask &= _ToBoolMask(result, nCoords, constraint)

    logger.debug(
        f"User constraints at {obsTime.iso}: {nCoords - np.count_nonzero(mask)} "
        f"of {nCoords} candidates removed."
    )
    return mask


def _ImportFile(path):
    """Import a python file given by its path and return the module."""
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError(f"User constraints file not found: {path}")
    key = (str(path), path.stat().st_mtime_ns)
    if key not in _LOADED_FILES:
        name = f"_tilepy_user_constraints_{len(_LOADED_FILES)}"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        try:
            spec.loader.exec_module(module)
        except BaseException:
            del sys.modules[name]
            raise
        _LOADED_FILES[key] = module
    return _LOADED_FILES[key]


def LoadUserConstraints(spec, baseDir=None):
    """
    Load user-defined constraints from references to python objects.

    Parameters
    ----------
    spec : str
        One or several references separated by commas or new lines. Each
        reference is either ``path/to/file.py:attribute`` or
        ``package.module:attribute``, where the attribute is a constraint or a
        list of constraints (see `EvaluateUserConstraints`).
    baseDir : str or pathlib.Path, optional
        Directory used to resolve relative file paths (the directory of the
        configuration file when called from
        `ObservationParameters.from_configfile`). Default: current directory.

    Returns
    -------
    constraints : list
        The constraints, in the order they are given.

    Notes
    -----
    The referenced python file or module is imported, i.e. its code is
    executed. Only use references to code you trust.
    """
    constraints = []
    for entry in re.split(r"[,\n]", spec or ""):
        entry = entry.strip()
        if not entry:
            continue
        ref, sep, attr = entry.rpartition(":")
        ref, attr = ref.strip(), attr.strip()
        if not sep or not ref or not all(p.isidentifier() for p in attr.split(".")):
            raise ValueError(
                f"Invalid user constraint reference '{entry}': expected "
                "'path/to/file.py:attribute' or 'package.module:attribute'."
            )

        if ref.endswith(".py") or "/" in ref or "\\" in ref:
            path = Path(ref).expanduser()
            if not path.is_absolute():
                path = Path(baseDir or ".") / path
            module = _ImportFile(path)
        else:
            module = importlib.import_module(ref)

        obj = module
        for part in attr.split("."):
            try:
                obj = getattr(obj, part)
            except AttributeError:
                raise ValueError(
                    f"User constraint reference '{entry}': '{attr}' not found in '{ref}'."
                ) from None

        found = list(obj) if isinstance(obj, (list, tuple)) else [obj]
        for constraint in found:
            _CheckConstraint(constraint)
        constraints.extend(found)

    logger.info(
        f"Loaded {len(constraints)} user constraint(s): "
        f"{', '.join(_ConstraintName(c) for c in constraints)}"
    )
    return constraints
