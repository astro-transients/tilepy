# Tilepy

## Installation

We clone the repo, create an environment to work, activate the enviroment and install the package.

```python
git clone git@drf-gitlab.cea.fr:multimessenger-IRFU/general-mwl-mm-transients-analyses/tilepy.git
conda env create -n tilepyenv -f environment.yml
conda activate tilepyenv
python -m pip install -e .
```

If you prefer to avoid conda and use a virtual environment with your favorite python version, use the following sequence:
```
python -m venv tilepy_venv
source tilepy_venv/bin/activate
pip3 install --upgrade pip
python -m pip install -r requirements.txt
```

Requirements of the installation: 

- The current version of the package **only** runs with `python=3.9`. Careful as well with the versions of matplotlib and healpy, they should be the ones explicited in the requirements.yml, otherwise there will be conflicts between them when plotting skymaps.  
- Note that by creating the env from the environment.yml, the libraries and versions needed will be installed authomatically.
- Note that everytime we made changes to the package, you should update the installation of the package doing ```pip install .``` in the folder where the setup.py is. The changes will be only applied to the env where you are working. 

In the case you are working in CC-Lyon, the easiest solution is to do```ccenv conda ``` and then follow the instructions given above. 

## Description

Package including functions to perform GW follow-up scheduling and simulations in IACTS. The package can be found in the folder tilepy, which contains the following folders: 
 
- tilepy: Folder including the package
    - tilepy.tools: Includes several scripts that have been used so far for different aims related to visualization and catalog cleaning 
    - tilepy.include: The main functions used by the two main scripts are in this folder. It includes the Pointing Tools specifically for CTA (the others are in GWHESSPointing tool which is imported by GWCTAPointingTools), the CTA observation scheduler, simulation tools and analysis tools (both using gammapy)
    - tilepy.dataset: Some files that have been used in the past and in ongoing simulation efforts. The important file in this folder finals2000.all. To have the entire dataset, 

- relics: Old scripts that may be useful in the future

- examples: Example on how to use tilepy:
    - Tiling_Observations.py: example of a python script to run a tiling schedule. Default values works for the alertType = gw
    - config/ExampleConfig.ini includes the minimal required parameters to run a tile. These are 
       - [observatory]
          - name: name of the observatory (it is not critical, you can use any name) 
          - lat: lat coordinates of the observatory 
          - lon: lon coordinates of the observatory 
          - height: height of the observatory 

        - [visibility]
          - sunDown: altitude coordinate of the Sun in deg to define darkness conditions (for astronomic darkness sunDown= -18)
          - horizonSun: altitude coordinate of the Sun in hh:mm:ss to define darkness conditions (for astronomic darkness horizonSun= -18:00:00)
          - moonDown: altitude coordinate of the Moon in deg to define darkness conditions (for astronomic darkness moonDown= -0.5)
          - horizonMoon: altitude coordinate of the Moon in hh:mm:ss to define darkness conditions (for astronomic darkness moonDown= -00:30:00)
          - moonGrey: altitude coordinate of the Moon in deg to define greyness conditions
          - moonPhase: phase of the Moon to define greyness conditions
          - moonSourceSeparation: minimum separation Source-Moon in deg to define greyness conditions
          - maxMoonSourceSeparation: max separation Source-Moon in deg to define greyness conditions

        - [operations]
          - maxZenith: max zenith angle which will be considered as accesible sky
          - FOV: radious of the circular FoV defining the tiles
          - maxRuns: maximum number of tiles that could be scheduled
          - maxNights: total number of nights considered
          - duration: standard exposure per tile 
          - minDuration: minimal duration of tile if standard exposure is not allocable
          - useGreytime: flag to schedule greyness observations in addition to darkness

        - [tiling]
          - online: ??
          - minimumProbCutforCatalogue: minimal 3D(GWxgalaxyCat)probability covered per tile to schedule observation
          - minProbcut:  minimal 2D probability covered per tile to schedule observation
          - distCut: distance cut to define the mandatory use of 2D strategy (coming from galaxy catalogue completeness)
          - doPlot: produce detailed various plots of the tile
          - secondRound: consider two maps for scheduling 
          - zenithWeighting: weight on probability that would be applied to prioritize coordinates that have lower zenith angle values. Step size is 5 deg in zenith (0.75 is a resonable value)
          - percentageMOC: percentage of the sky localization region that would be considered to compute the MOC
          - reducedNside: nside of the low-resolution skymap used as a grid to speed up the computation
          - HRnside: nside of the high-resolution skymap used to compute the probability covered
          - mangrove: flag to use the mangrove method of weighting by the mass of the host galaxy

## Creation of the reduced galaxy file

For using the 3D algorithm tilepy need to have access to a galaxy catalog. Currently, the only supported catalog is GLADE+. You could found the download link on this webpage : https://glade.elte.hu.
Then to prepare it for utilisation by tilepy you could use the `ConvertGalaxyCatalog.py` script to convert it to a file format compatible with tilepy. It is located in the `tilepy/tools` repository.

Example of use to keep only galaxies that are below 500 Mpc (recommended). It read the downloaded `GLADE+.txt` file and the output file `converted_GLADE.h5` is the one that should be used with tilepy.

```bash
python ConvertGalaxyCatalog.py --input GLADE+.txt --output converted_GLADE.h5 --max-luminosity-distance 500
```

## Help
If any problem is found, please open an Issue in the project main page so it is documented. You can of course also directly create a Pull Request with the new features implemented.

Otherwise you can also contact us at astro.tilepy@gmail.com. 
