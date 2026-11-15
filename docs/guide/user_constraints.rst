========================
User-defined constraints
========================

|tilepy| applies a set of built-in visibility cuts (Sun and Moon altitude, Moon
phase, Moon-source separation, zenith angle) that are configured in the
``[visibility]`` and ``[operations]`` sections of the configuration file.

On top of those, you can provide **your own constraints** without modifying
|tilepy|, for example a limit on the airmass, or a cut on an
instrument-specific quantity such as the expected night-sky background in the
camera, which depends on the pointing direction and on the position and phase
of the Moon.

This layer is optional: without user-defined constraints the schedules are
exactly the same as before, and nothing changes in the configuration files, in
the command line or in the output files.


How constraints are applied
===========================

At each time step of the scheduling, every candidate pointing that passed the
built-in cuts is tested against the user-defined constraints. Candidates that
do not fulfil **all** of them are discarded for that time step, and the best
remaining candidate is scheduled. A constraint therefore does not simply remove
pointings from the final schedule: the scheduler looks for the best pointing
that is allowed at that time.

A constraint is evaluated at the **centre of the field of view** of the
candidate pointing and at the **start time** of the observation.


Writing a constraint
====================

Two forms are accepted, and they can be mixed.

astroplan constraints
---------------------

Any `astroplan <https://astroplan.readthedocs.io>`_ ``Constraint`` can be used:
the ones shipped with astroplan (``AirmassConstraint``,
``MoonSeparationConstraint``, ``LocalTimeConstraint``, ...) as well as your own
subclasses, written as described in the
`astroplan documentation <https://astroplan.readthedocs.io/en/latest/tutorials/constraints.html#user-defined-constraints>`_.
|tilepy| builds the astroplan ``Observer`` from the location of the observatory
given in the configuration file.

.. code-block:: python

   # my_cuts.py
   import astropy.units as u
   from astroplan import AirmassConstraint, Constraint


   class MoonBrightnessCut(Constraint):
       """Toy proxy for the sky brightness induced by the Moon."""

       def __init__(self, max_proxy):
           self.max_proxy = max_proxy

       def compute_constraint(self, times, observer, targets):
           altaz = observer.altaz(times, targets)
           moon = observer.moon_altaz(times)
           illumination = observer.moon_illumination(times)
           proxy = illumination / altaz.separation(moon).to_value(u.deg)
           return proxy < self.max_proxy


   cuts = [MoonBrightnessCut(max_proxy=0.02), AirmassConstraint(max=2.0)]

Plain functions
---------------

A constraint can also be any callable with the signature
``f(time, coords, obspar)``, returning one boolean per candidate (``True`` if
the pointing is allowed):

.. code-block:: python

   def northern_sky_only(time, coords, obspar):
       return coords.dec.deg > 0

- ``time`` is a scalar `~astropy.time.Time` (UTC),
- ``coords`` is a `~astropy.coordinates.SkyCoord` array with the candidate
  pointings,
- ``obspar`` is the
  :class:`~tilepy.include.CampaignDefinition.ObservationParameters` of the
  observatory (``obspar.location``, ``obspar.maxZenith``, ...).

A single boolean is accepted for constraints that only depend on time. If a
constraint returns numbers instead of booleans (e.g. astroplan constraints
with ``boolean_constraint=False``), a candidate is allowed when the value is
larger than zero.


Using the constraints
=====================

From python
-----------

Set the ``userConstraints`` attribute of the observation parameters before
calling the scheduler:

.. code-block:: python

   from my_cuts import cuts

   obspar.userConstraints = cuts
   GetSchedule(obspar)

With :func:`~tilepy.include.ObservationScheduler.GetUniversalSchedule`, each
observatory has its own ``ObservationParameters`` and therefore its own
constraints.

From the configuration file
---------------------------

Alternatively, add the optional key ``userConstraints`` to the ``[visibility]``
section. This also works with the ``Tiling_Observations`` command line tool:

.. code-block:: ini

   [visibility]
   sunDown = -18
   moonDown = -0.5
   userConstraints = my_cuts.py:cuts

Each entry is ``reference:attribute``, where the reference is either the path
to a python file (relative paths are resolved with respect to the directory of
the configuration file) or the name of an installed module
(``mypackage.constraints:cuts``). The attribute is a constraint or a list of
constraints. Several entries can be given, separated by commas or new lines.

.. warning::

   The python file or module referenced by ``userConstraints`` is imported, so
   its code is executed. Only use configuration files you trust, and do not
   build this key from untrusted input, e.g. in a service exposed to external
   users.


Current limitations
===================

- User-defined constraints are applied to **ground-based** observatories
  (single and multi-observatory scheduling, 2D and 3D algorithms). They are not
  applied to space-based observatories.
- The observation windows reported by the ranking tools (``doRank``) only take
  the built-in cuts into account.
- With the 3D ``integrated`` strategy only the most probable galaxies are
  considered as pointing candidates at each time step. If all of them are
  excluded by a constraint, no pointing is scheduled for that time step.


Reference/API
=============

.. automodapi:: tilepy.include.UserConstraints
    :no-inheritance-diagram:
