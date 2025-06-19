Scheduling and Tiling Tools (`tilepy.include`)
==============================================

The `tilepy.include` package provides core modules for scheduling, tiling computation, and visualization
in the |tilepy| workflow. It enables rapid follow-up observations, campaign configuration, and plotting for
multi-messenger astronomy.

Observatories
=============

The **Observatories** subpackage defines the main ground- and space-based facilities supported in the |tilepy| workflow,
including location, field of view, and technical constraints for each telescope.

.. automodapi:: tilepy.include.Observatories

Observation Scheduling
=====================

The **ObservationScheduler** module provides algorithms and utilities to schedule rapid follow-up observations
of electromagnetic counterparts to gravitational wave transients. It computes optimal tiling schedules for telescopes,
generates visibility plots, and supports both single- and multi-observatory coordination.

.. automodapi:: tilepy.include.ObservationScheduler

Other Modules
=============

.. automodapi:: tilepy.include
