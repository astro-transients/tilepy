==========
Quickstart
==========

.. button-link:: https://colab.research.google.com/github/astro-transients/tilepy/blob/main/docs/tutorials/quickstart_tilepy.ipynb
   :color: info
   :shadow:

   Open the notebook in Colab

.. note::

   This page presents the **full workflow** of |tilepy|, from the localization
   map to the observation schedule: display the sky map, extract the 90% credible
   region and its pixels, then run the scheduling. The companion notebook (button
   above) follows exactly these steps and can be run directly.

.. admonition:: What you will see
   :class: info

   #. What a **HEALPix sky map** looks like and how to display it.
   #. How to **extract the 90%** probability and **which pixels** make it up.
   #. How |tilepy| turns this region into a **schedule of pointings**.


The sky map
===========

A localization map is a **HEALPix** map: the sky is tiled into pixels of equal
area and each pixel ``i`` carries the probability that the source lies inside it.
|tilepy| reads these maps internally (module
:mod:`tilepy.include.MapManagement`), but to clearly visualize the starting point
one can display it with ``ligo.skymap``.

.. plot::
   :include-source: False

   import healpy as hp
   import matplotlib.pyplot as plt
   from ligo.skymap.io import read_sky_map
   import ligo.skymap.plot  # registers the "astro ..." projections

   url = ("https://gracedb.ligo.org/api/superevents/S190728q/files/"
          "GW190728_064510_PublicationSamples_flattened.fits.gz,0")

   prob, meta = read_sky_map(url)          # prob[i]: probability of pixel i
   nside = hp.npix2nside(len(prob))

   ax = plt.axes(projection="astro hours mollweide")
   ax.imshow_hpx((prob, "ICRS"), cmap="cylon")
   ax.grid()

.. tip::

   The figure above is produced by this snippet **when the docs are built** (it
   downloads only the small GW190728 skymap, a few MB). ``read_sky_map`` returns
   a **flattened** map (a single NSIDE); recent IGWN alerts are often
   multi-resolution (MOC): in that case use ``read_sky_map(url, moc=True)``.
   |tilepy| handles both formats automatically.


HEALPix resolution and order
============================

HEALPix is **hierarchical**: the sphere starts as 12 base pixels (order
:math:`\ell = 0`) and each pixel is recursively split into 4. The **order** (or
*level*) :math:`\ell` therefore sets the resolution:

.. math::

   \mathrm{NSIDE} = 2^{\ell}, \qquad
   N_\mathrm{pix} = 12\ \times \mathrm{NSIDE}^2 = 12 \cdot 4^{\ell}.

Every pixel covers the **same** solid angle:

.. math::

   \Omega_\mathrm{pix} = \frac{4\pi}{N_\mathrm{pix}}
       = \frac{41253\ \mathrm{deg}^2}{N_\mathrm{pix}}, \qquad
   \theta_\mathrm{pix} \approx \sqrt{\Omega_\mathrm{pix}}
       \approx \frac{58.6°}{\mathrm{NSIDE}}.

**Worked example** (:math:`\ell = 6`):
:math:`\mathrm{NSIDE} = 2^6 = 64`,
:math:`N_\mathrm{pix} = 12\cdot64^2 = 49152`,
:math:`\Omega_\mathrm{pix} = 41253/49152 \approx 0.84\ \mathrm{deg}^2`,
:math:`\theta_\mathrm{pix} \approx 58.6/64 \approx 0.92°`.

.. list-table:: Resolution as a function of the order
   :header-rows: 1
   :widths: 10 12 22 24 22

   * - :math:`\ell`
     - NSIDE
     - :math:`N_\mathrm{pix}`
     - Pixel area
     - Pixel size
   * - 6
     - 64
     - 49,152
     - :math:`\approx 0.84\ \mathrm{deg}^2`
     - :math:`\approx 0.92°`
   * - 8
     - 256
     - 786,432
     - :math:`\approx 0.052\ \mathrm{deg}^2`
     - :math:`\approx 0.23°`
   * - 9
     - 512
     - 3,145,728
     - :math:`\approx 0.013\ \mathrm{deg}^2`
     - :math:`\approx 0.11°`

Going from :math:`\ell` to :math:`\ell+1` multiplies the number of pixels by 4
and halves the pixel size. This is the multi-resolution "tree" used by IGWN MOC
maps: well-localized regions are stored at high :math:`\ell`, empty regions at
low :math:`\ell`. The figure below (rendered at build time) shows the same sphere
pixelized at two orders:

.. plot::

   import numpy as np, healpy as hp
   import matplotlib.pyplot as plt

   fig = plt.figure(figsize=(8, 3.2))
   for i, ns in enumerate([2, 8]):
       hp.mollview(np.arange(hp.nside2npix(ns)), sub=(1, 2, i + 1), cbar=False,
                   cmap="tab20",
                   title=f"l={int(np.log2(ns))}  (NSIDE={ns}, Npix={hp.nside2npix(ns)})")
       hp.graticule()

.. important::

   |tilepy| works with **two** orders, configurable in the ``.ini`` file:

   * ``HRnside`` (high :math:`\ell`) to integrate the probability inside a field
     of view precisely;
   * ``reducedNside`` (low :math:`\ell`) to scan the sky quickly.


Extracting the 90% region
==========================

**Principle (greedy algorithm).** Sort the pixels from most to least probable,
accumulate their probability, and keep the first ones until the desired fraction
is reached (here 0.9). The area is simply the number of retained pixels times the
area of one pixel.

.. tab-set::

   .. tab-item:: With healpy (pedagogical)

      .. code-block:: python

         import numpy as np

         order = np.argsort(prob)[::-1]      # decreasing sort
         cumul = np.cumsum(prob[order])      # cumulative sum
         pix_90 = order[cumul <= 0.90]       # pixels of the 90%

         pix_area = hp.nside2pixarea(nside, degrees=True)
         area_90 = len(pix_90) * pix_area
         print(f"{len(pix_90)} pixels  ->  {area_90:.1f} deg²")

   .. tab-item:: With tilepy

      Once the map is loaded into a
      :class:`~tilepy.include.MapManagement.SkyMap.SkyMap` object, the same
      computation is directly available:

      .. code-block:: python

         skymap.getArea(0.9)          # area of the 90% region (deg²)
         skymap.getPixIdArea(0.9)     # indices of the pixels in the region

      And to retrieve the (RA, Dec) coordinates on a reduced grid — this is the
      list of positions the algorithm then scans:

      .. code-block:: python

         from tilepy.include.PointingTools import GetRegionPixReduced
         ra, dec, area = GetRegionPixReduced(prob, 0.9, reducedNside)

.. seealso::

   The functions :func:`~tilepy.include.PointingTools.GetRegionPixReduced` and
   :meth:`~tilepy.include.MapManagement.SkyMap.SkyMap.getPixIdArea` implement this
   sort/accumulate step. It is the same "cumulative sum of sorted pixels" scheme
   found in most multi-messenger follow-up tools.


The scheduling process
=======================

|tilepy| is driven by a ``.ini`` **configuration file** (telescope, visibility
constraints, strategy) and an
:class:`~tilepy.include.CampaignDefinition.ObservationParameters` object. We stay
here on the simplest case: ``algorithm = 2D`` (probability integration) for a
ground-based telescope.

.. important::

   In **2D** mode tilepy does **not** read any galaxy catalog, so there is
   **nothing heavy to download** (the GLADE+ catalog ``Gladeplus.h5`` is only
   needed in 3D mode). Pass ``None`` for both the dataset directory and the
   catalog name.

.. dropdown:: Minimal configuration file (LST, 2D)

   .. code-block:: ini

      [observatory]
      name = LST
      lat = 28.75
      lon = -17.5
      height = 2200
      base = ground

      [visibility]
      sundown = -18
      moondown = -0.5
      earthdown = 0
      moongrey = -0.5
      gmoonphase = 60
      minmoonsourceseparation = 30
      maxmoonsourceseparation = 150
      SAAThreshold = 0

      [operations]
      maxzenith = 60
      fov = 2.0
      maxRuns = 20
      maxNights = 1
      duration = 30
      minduration = 10
      useGreytime = False
      minSlewing = 0
      shape = circle
      numberSides = 0
      FoVRotation = 0

      [tiling]
      locCut = 99999
      minimumprobcutforcatalogue = 0.01
      minProbcut = 0.01
      distcut = 0.5
      doPlot = True
      secondRound = False
      zenithWeighting = 0.75
      percentageMOC = 0.9
      reducedNside = 64
      hrnside = 512
      mangrove = False
      algorithm = 2D
      strategy = integrated
      doRank = False
      countPrevious = False
      alphaR = 0
      betaR = 0

      [general]
      downloadMaxRetry = 1
      downloadWaitPeriodRetry = 20

.. dropdown:: Run the computation
   :open:

   .. tab-set::

      .. tab-item:: From a skymap (small download)

         .. code-block:: python

            import datetime
            from tilepy.include.CampaignDefinition import ObservationParameters
            from tilepy.include.ObservationScheduler import GetSchedule

            obsTime = datetime.datetime.fromisoformat("2019-07-28 08:30:00")

            obspar = ObservationParameters()
            # datasetDir = None, galcatName = None: 2D needs no galaxy catalog
            obspar.add_parsed_args(url, obsTime, None, None,
                                   "./output", None, None, "GW190728")
            obspar.from_configfile("FollowupParameters_LST.ini")

            GetSchedule(obspar)

      .. tab-item:: Fully offline (no download)

         tilepy can build a **simulated Gaussian map** in memory with
         :func:`~tilepy.include.CampaignDefinition.set_gaussian_source` — same
         scheduling, zero network access.

         .. code-block:: python

            from tilepy.include.CampaignDefinition import (
                ObservationParameters, set_gaussian_source)
            from tilepy.include.ObservationScheduler import GetSchedule

            obspar = ObservationParameters()
            obspar.add_parsed_args(None, obsTime, None, None,
                                   "./output", None, None, "sim_event")
            set_gaussian_source(obspar, ra=240.0, dec=20.0, sigma=5.0)
            obspar.nside = 256
            obspar.from_configfile("FollowupParameters_LST.ini")

            GetSchedule(obspar)

.. admonition:: What does ``GetSchedule`` do in 2D mode?
   :class: note

   #. reads the map and computes the 90% region (previous section);
   #. splits the night into observable **dark-time windows** at the site;
   #. at each window, centers the field of view on the pixel with the **highest
      probability still visible**, adds up the probability contained in the FoV,
      then **masks** those pixels so they are not counted again;
   #. repeats until ``maxRuns`` pointings (or until the probability above
      ``minProbcut`` is exhausted).

   The result is a table of pointings
   ``SuggestedPointings_2DProbOptimisation.txt`` (columns ``Time[UTC]``,
   ``RA[deg]``, ``DEC[deg]``, ``PGW``, ...), together with the figures if
   ``doPlot = True``.


Reading the produced schedule
==============================

.. code-block:: python

   from astropy.table import Table
   res = Table.read(
       "output/GW190728/PGinFoV/SuggestedPointings_2DProbOptimisation.txt",
       format="ascii")
   print(res)
   print("Total covered probability:", res["PGW"].sum())


.. seealso::

   * :doc:`guide/observation_scheduler` for the full description of the modules.
   * The |tilepyGitHub| for examples and advanced configuration (3D mode, galaxy
     catalog, multi-telescope).

.. footbibliography::
