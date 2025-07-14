Scheduling and Tiling Tools (`tilepy.include`)
==============================================

The `tilepy.include` package provides core modules for scheduling, tiling computation, and visualization
in the |tilepy| workflow. It enables rapid follow-up observations, campaign configuration, and plotting for
multi-messenger astronomy.

Read Skymap
-----------

Provides classes and utilities for reading and interpreting localization skymaps.

.. automodapi:: tilepy.include.MapManagement.SkyMap
    :no-inheritance-diagram:


Observatories
-------------

Defines the main ground- and space-based facilities supported in the |tilepy| workflow,
including location, field of view, and technical constraints for each telescope.

.. automodapi:: tilepy.include.Observatories
    :no-inheritance-diagram:

Observatoire Configuration
--------------------------

.. automodapi:: tilepy.include.CampaignDefinition
    :no-inheritance-diagram:

Tiling Determination
--------------------

Functions to compute tiling observation schedules for one or multiple observatories,
using 2D/3D probability maps and considering constraints.

.. automodapi:: tilepy.include.TilingDetermination
     :no-inheritance-diagram:

Observation Scheduling
----------------------

Algorithms and utilities to schedule rapid follow-up observations, generate visibility plots, and support coordination.

.. automodapi:: tilepy.include.ObservationScheduler
     :no-inheritance-diagram:

Ranking Observation Times
-------------------------

Tools to rank observations by probability covered, adding the observability window for a comprehensive view.

.. automodapi:: tilepy.include.RankingObservationTimes
     :no-inheritance-diagram:

Plotting Tools
--------------

Utilities for plotting sky maps, pointings, and other visualization tasks.

.. automodapi:: tilepy.include.PointingTools
     :no-inheritance-diagram:
