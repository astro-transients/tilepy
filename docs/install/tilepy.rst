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

   .. button-link:: https://docs.conda.io/projects/conda/en/stable/
      :color: info
      :shadow:

      Conda

   **or**

   .. button-link:: https://docs.python.org/3/library/venv.html
      :color: info
      :shadow:

      virtualenv

   - (Optional) mamba for faster environment setup


.. dropdown:: Environment setup

   .. tab-set::

      .. tab-item:: Conda

         .. code-block:: bash

            conda create -n tilepy_env python=3.9
            conda activate tilepyenv

      .. tab-item:: Python

         .. code-block:: bash

            python -m venv tilepy_venv
            source tilepy_venv/bin/activate
            pip install --upgrade pip

.. dropdown:: Install tilepy

   .. tab-set::

      .. tab-item:: PyPI

         .. admonition:: Recommended
            :class: tip

            The recommended way to install tilepy is using pip:

            .. code-block:: console

               $ pip install tilepy

      .. tab-item:: Latest

         .. admonition:: Unreleased version
            :class: danger

            If you want to use the latest version of tilepy, you can install the current main branch directly from GitHub:

            .. code-block:: console

               $ pip install git+https://github.com/astro-transients/tilepy.git@main

      .. tab-item:: Development

         .. admonition:: Development install
            :class: info

            To install tilepy in development mode, run:

            .. code-block:: bash

               git clone https://github.com/astro-transients/tilepy.git
               cd tilepy
               conda env create -n tilepyenv -f environment.yml
               conda activate tilepyenv
               pip install -e .


   .. warning:: Troubleshooting MOCpy Installation (If Needed)

      If you have issues with the `mocpy` package when using conda, try installing it separately with pip after creating your environment:

      .. code-block:: console

         $ pip install mocpy


.. dropdown:: Verify your installation

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
