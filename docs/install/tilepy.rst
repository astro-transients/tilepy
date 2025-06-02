
TilePy Installation Guide
==========================

This is the process to install end to end TilePy.

1. Clone and Install TilePy
---------------------------

Create and Activate a Conda Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/astro-transients/tilepy.git
   cd tilepy
   conda env create -n tilepyenv -f environment.yml
   conda activate tilepyenv
   pip install -e .

If you use `mamba`, you can replace `conda` by `mamba` for faster environment creation.

Alternative: Install with Virtualenv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python -m venv tilepy_venv
   source tilepy_venv/bin/activate
   pip install --upgrade pip
   pip install -e .

*Note: TilePy requires Python >= 3.9.*

2. Troubleshooting MOCpy Installation (If Needed)
-------------------------------------------------

If you have issues with the `mocpy` package when using conda, try installing it separately with pip after creating your environment :

.. code-block:: bash

   pip install mocpy

3. Verify Installation
----------------------

.. code-block:: bash

   python -c "import tilepy; print(tilepy.__version__)"

If no error appears and the version number is printed, your installation is correct!

4. Additional Notes
-------------------

- **Dependencies** are automatically managed if you use the provided environment.
- For advanced usage (simulations, scripts, or examples), check the `examples/ <https://github.com/astro-transients/tilepy/tree/main/examples>`__ directory.


5. Need Help?
-------------

- Open an issue on the |TilepyGitHub|
- Contact the dev team at: |TilepyEmail|
- Join the forum: |Forum|
