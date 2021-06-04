# Gamma GW Followup

## Installation


```python
git clone 
conda env create -n gammagwfollowup -f environment.yml
conda activate gammagwfollowup
pip install .      
```
Requirements of the installation: 

- The current version of the package **only** runs with `python=3.6`.
- The rest of the requirements are specified in the environment.yml, including the gammapy version of the codes. 

## Description


Package including functions to perform GW follow-up scheduling and simulations in IACTS. The package can be found in the folder gammagwfollowup, which contains the following folders: 
 
- tools: Includes several scripts that have been used so far for different aims. 
- include: The main functions used by the two main scripts are in this folder. It includes the Pointing Tools specifically for CTA (the others are in GWHESSPointing tool which is imported by GWCTAPointingTools), the CTA observation scheduler, simulation tools and analysis tools (both using gammapy)
- dataset: Some files that have been used in the past and in ongoing simulation efforts. The important file in this folder finals2000.all


## Help
If any problem is found, please open an Issue in the project main page so it is documented. 

Otherwise you can also contact me at seglararro@lapp.in2p3.fr, or you can create a Pull Request with the new features implemented 