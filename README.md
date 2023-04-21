# Tilepy

## Installation

We clone the repo, create an environment to work, activate the enviroment and install the package.

```python
git clone git@drf-gitlab.cea.fr:multimessenger-IRFU/general-mwl-mm-transients-analyses/tilepy.git
conda env create -n tilepyenv -f environment.yml
conda activate tilepyenv
python -m pip install -e .
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

- examples: Example on how to use tilepy (work ongoing) 

## Creation of the reduced galaxy file

For using the 3D algorithm tilepy need to have access to a galaxy catalog. Currently, the only supported catalog is GLADE+. You could found the download link on this webpage : https://glade.elte.hu.
Then to prepare it for utilisation by tilepy you could use the `ConvertGalaxyCatalog.py` script to convert it to a file format compatible with tilepy. It is located in the `tilepy/tools` repository.

Example of use to keep only galaxies that are below 500 Mpc (recommended). It read the downloaded `GLADE+.txt` file and the output file `converted_GLADE.h5` is the one that should be used with tilepy.

```bash
python ConvertGalaxyCatalog.py --input GLADE+.txt --output converted_GLADE.h5 --max-luminosity-distance 500
```

## Help
If any problem is found, please open an Issue in the project main page so it is documented. 

Otherwise you can also contact me at seglararro@lapp.in2p3.fr, or you can create a Pull Request with the new features implemented.
