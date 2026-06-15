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

      Python >= 3.9

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

            Install the current main branch directly from GitHub:

            .. code-block:: console

               $ pip install git+https://github.com/astro-transients/tilepy.git@main

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



.. note::

   - **Dependencies** are automatically managed if you use the provided environment.
   - For advanced usage (simulations, scripts, or examples), check the `examples/ <https://github.com/astro-transients/tilepy/tree/main/examples>`__ directory.


.. seealso::

   - Open an issue on the |tilepyGitHub|
   - Contact the dev team at: |tilepyEmail|
   - Join the forum: |Forum|
