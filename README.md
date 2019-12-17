This project includes the scheduling of GW-followup observation by CTA, simulation of the CTA observations and analysis of those same observation (both using  gammapy-0.9). The aim is to  give as an output the detectability of the GRB associated to the GW with CTA. 

The main two scripts are: 

- Pointing_PGWonFoV.py: Script to schedule GW follow-up observations with pointing telescopes by integrating the probability density covered per pointing in different regions of the GW C.R. 

- Pointing_PgalonFoV.py: Script to schedule GW follow-up observations with pointing telescopes by integrating the probability density from GW observation x galaxy distribution covered per pointing in different regions of the sky.

These two scripts take as input the number of the GWCosmos simulation. Then, other related information is selected by using the MergerID and RunNumber. These are:

- The associated GRB emission is selected. 
- The galaxy catalog
- The input time list
- Set of observation times computed to obtain 5-sigma detection post-trials.

In the same directory, there are other folders which are:
 
- tools: Includes several scripts that have been used so far for different aims. There are not related to the main code and contain paths which are not general. A short description of each of them will be included directly in the scripts.
- configs: The .ini files include the parameters used in each of the two scripts.
- include: The main functions used by the two main scripts are in this folder. It includes the Pointing Tools specifically for CTA (the others are in GWHESSPointing tool which is imported by GWCTAPointingTools), the CTA observation scheduler, simulation tools and analysis tools (both using gammapy)
- output: The results will be stored in this folder.
- testFiles: Not specially important. It contains files that were used in the past to make tests. Included since it could maybe be useful in the future.

Things to note: 

- In particular, the gammapy modules will obviously just work in a gammapy env., in particular v0.9.
- The scripts as they currently are use data from a dataset folder which should be placed relatively to this folder as  ../dataset/. 
- The script uses a GWHESSPointingTools module and PointingPlotting which should be also placed in the correct directory ('../../GW-followup/include')

In the dataset folder, one needs to include the following elements (to download):

- GW simulations: from GWCosmos by Patricelli et al. https://figshare.com/collections/GW_COSMoS_Gravitational_Wave_COmpact_binary_SysteM_Simulations/4243595
- BNS-GW.txt: lit of GWCosmos events with the different associated parameters.
- BNS-GW-new.txt same as previous but includes max and min distance
- GRB simulations (associated to GWCosmos): found in Redmine https://forge.in2p3.fr/projects/cta-gravitational-waves/files 
- CTA IRFs from https://www.cta-observatory.org/science/cta-performance/
- Galaxy catalogs: a galaxy catalog is associated to each of the run identifiers of the GW simulations (a sub-sample is included in the dataset folder)

Simulations have been performed of a subsample of events observed on-axis (theta<5deg) and for a priori selected times. These are included in the dataset folder.

- BNS-GW_onAxis5deg.txt
- BNS-GW_Time_onAxis5deg.txt
- Selected observation times: The observation intervals to obtain 5sigma post-trials detection obtained using cssens fro ctools.

