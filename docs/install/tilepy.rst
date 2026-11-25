============
Installation
============

This guide describes the step-by-step process to install **tilepy** and its dependencies.

.. note::

   **tilepy**  is a toolkit designed to simulate observations based on gravitational wave (GW) skymaps, and to schedule follow-up observations (gamma-ray bursts, neutrinos, and optical transients) with both space- and ground-based telescopes, such as HESS, LST, and other CTA telescopes.

   We strongly recommend installing tilepy inside a virtual environment (`conda` or `venv`) to avoid dependency conflicts.


.. dropdown:: Requirements

   You will need the following:

   .. button-link:: https://www.python.org/downloads/
      :color: info
      :shadow:

      Python >= 3.11

.. dropdown:: ⚡ Fastest & simplest install -- uv :bdg-success:`recommended`

   `uv <https://docs.astral.sh/uv/>`_ installs tilepy and all its dependencies
   in a single command, with no environment setup needed.

   .. tab-set::

      .. tab-item:: Linux / macOS

         .. code-block:: bash

            curl -LsSf https://astral.sh/uv/install.sh | sh
            git clone https://github.com/astro-transients/tilepy.git
            cd tilepy
            uv sync --python 3.11

      .. tab-item:: Windows

         .. code-block:: powershell

            powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
            git clone https://github.com/astro-transients/tilepy.git
            cd tilepy
            uv sync --python 3.11


   Once installed, run tilepy commands directly with ``uv run`` -- no need to activate any environment:

   .. tab-set::

      .. tab-item:: Using uv run (recommended)

         .. code-block:: bash

            uv run python -c "import tilepy; print(tilepy.__version__)"
            uv run Tiling_Observations ...

      .. tab-item:: Activating the environment (optional)

         .. code-block:: bash

            source .venv/bin/activate
            python -c "import tilepy; print(tilepy.__version__)"


.. dropdown:: Other installation methods

   Prefer another tool? Choose one of the following:

   .. tab-set::

      .. tab-item:: Conda

         **1. Create and activate the environment:**

         .. code-block:: bash

            conda create -n tilepy_env python=3.11
            conda activate tilepy_env

         **2. Install tilepy:**

         .. code-block:: bash

            pip install tilepy

      .. tab-item:: pip + venv

         **1. Create and activate the environment:**

         .. code-block:: bash

            python -m venv tilepy_venv
            source tilepy_venv/bin/activate  # Linux/macOS
            tilepy_venv\Scripts\activate     # Windows
            pip install --upgrade pip

         **2. Install tilepy:**

         .. code-block:: bash

            pip install tilepy

      .. tab-item:: GitHub (unreleased)

         .. admonition:: Unreleased version
            :class: danger

            Install the current master branch directly from GitHub:

            .. code-block:: console

               $ pip install git+https://github.com/astro-transients/tilepy.git@master

   .. warning:: Troubleshooting MOCpy

      If you have issues with ``mocpy`` under conda, install it separately:

      .. code-block:: console

         $ pip install mocpy


.. dropdown:: Developer install

   For contributors who want to modify the source code:

   .. tab-set::

      .. tab-item:: uv

         .. code-block:: bash

            git clone https://github.com/astro-transients/tilepy.git
            cd tilepy
            uv sync --python 3.11

      .. tab-item:: pip

         .. code-block:: bash

            git clone https://github.com/astro-transients/tilepy.git
            cd tilepy
            pip install -e ".[dev]"


.. dropdown:: Verify your installation

   .. tab-set::

      .. tab-item:: uv

         .. code-block:: bash

            uv run python -c "import tilepy; print(tilepy.__version__)"

      .. tab-item:: conda / pip

         .. code-block:: console

            $ python -c "import tilepy; print(tilepy.__version__)"

   If no error appears and the version number is printed, your installation is correct!


.. dropdown:: Ephemeris file (offline or restricted networks)

   The dark-time windows are computed with the JPL planetary ephemeris
   ``de440s.bsp`` (32 MB). Nothing has to be done on a machine with internet
   access: the file is downloaded once, into ``~/.skyfield/``, the first time
   a schedule is computed.

   **tilepy** uses the first of the following that is available:

   1. the ephemeris file (``.bsp``) given by the environment variable
      ``TILEPY_EPHEMERIS``;
   2. ``de440s.bsp`` in the current working directory;
   3. ``de440s.bsp`` in ``~/.skyfield/``;
   4. a download from the JPL server into ``~/.skyfield/``;
   5. the DE421 ephemeris shipped with the optional package ``skyfield-data``.

   On a machine that cannot reach the JPL server, either copy
   `de440s.bsp <https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp>`_ to
   ``~/.skyfield/``, or install the bundled ephemeris:

   .. code-block:: bash

      pip install "tilepy[offline]"

   A warning is logged when the DE421 ephemeris is used. It is valid until
   2053, and the Sun and Moon rising and setting times obtained with it agree
   with those of DE440 to about a millisecond.

   The Earth orientation parameters needed by Astropy (``finals2000A``) do
   not require any configuration: if they cannot be downloaded, Astropy uses
   the table shipped with the ``astropy-iers-data`` package. On a machine that
   stays offline, update that package from time to time
   (``pip install -U astropy-iers-data``) and disable the download attempts
   before importing **tilepy**:

   .. code-block:: python

      from astropy.utils import iers

      iers.conf.auto_download = False



.. note::

   - **Dependencies** are automatically managed if you use the provided environment.
   - For advanced usage (simulations, scripts, or examples), check the `examples/ <https://github.com/astro-transients/tilepy/tree/master/examples>`__ directory.


.. seealso::

   - Open an issue on the |tilepyGitHub|
   - Contact the dev team at: |tilepyEmail|
   - Join the forum: |Forum|
