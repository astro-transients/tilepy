The directories in this folder include various examples of how tilepy can be used. These are organized in:

a `./config` folder containing configuration files:
    ExampleConfig.ini for regular 1-obs
    FollowupParameters_[ObsName].ini for the N-obs scheduler.

b. `./launcher` folder contains example scripts and jupyter notebooks to launch a scheduling for all cases. It also contains the jupyter notebook to reproduce the figures of tilepy paper

c. `./sciencecases` folder contains the outputs corresponding to the launcher scripts. The default values in the examples correspond to the science cases in the tilepy paper (GW 1-obs, GRB N-obs, optical N-obs, IPN GRB 1-obs)

d. `./visualization` folder contains useful jupyter notebooks for visualization of the scheduling and the cumulative probability achieved. It also contains jupyter notebooks to reproduce the figures of the tilepy paper using the input files in ./sciencases.

e. `./paperplots` folder contains the plots connected with the jupyter notebooks in the visualization.

Run the Jupyter notebook .ipynb or the .py scripts with Python > 3.9
The paths are self-contained and the skymap examples are self-downloaded. The only external dataset that the user should download is a galaxy catalog in dataset (see ../README.md for details).
