# Tilepy

## Installation

We clone the repo, create an environment to work, activate the enviroment and install the package.

```python
git clone  git@drf-gitlab.cea.fr:multimessenger-IRFU/cta/gw-follow-up-simulations.git 
cd gw-follow-up-simulations
conda env create -n tilepyenv -f environment.yml
conda activate tilepyenv
pip install .      
```

Requirements of the installation: 

- The current version of the package **only** runs with `python=3.6`. Careful as well with the versions of matplotlib and healpy, they should be the ones explicited in the requirements.yml, otherwise there will be conflicts between them when plotting skymaps.  
- Note that by creating the env from the environment.yml, the libraries and versions needed will be installed authomatically.
- Note that everytime we made changes to the package, you should update the installation of the package doing ```pip install .``` in the folder where the setup.py is. The changes will be only applied to the env where you are working. 

In the case you are working in CC-Lyon, the easiest solution is to do```ccenv conda ``` and then follow the instructions given above. 

## Description


Package including functions to perform GW follow-up scheduling and simulations in IACTS. The package can be found in the folder tilepy, which contains the following folders: 
 
- tools: Includes several scripts that have been used so far for different aims. 
- include: The main functions used by the two main scripts are in this folder. It includes the Pointing Tools specifically for CTA (the others are in GWHESSPointing tool which is imported by GWCTAPointingTools), the CTA observation scheduler, simulation tools and analysis tools (both using gammapy)
- dataset: Some files that have been used in the past and in ongoing simulation efforts. The important file in this folder finals2000.all. To have the entire dataset, 


## Help
If any problem is found, please open an Issue in the project main page so it is documented. 

Otherwise you can also contact me at seglararro@lapp.in2p3.fr, or you can create a Pull Request with the new features implemented.
