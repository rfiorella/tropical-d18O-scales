# README

## How to run the model
1. Update the path for input and run files. Go to `code/initialize.py` and update `run_path` to the location of the `run_name` directory.
2. Navigate to the DISTRIBUTE directory in Terminal.
3. Run `source .venv/bin/activate` to initialize virtual environment
4. Navigate the `code` directory within the run_name in Terminal.
5. Run the model via `Python3 RUN.py` 


## Relevant directories and files
### `code` directory
Contains three files for running the model: 
1. `RUN.py` is the model run script. You shouldn't have to change anything here. 
2. `initialize.py` is where users set the run parameters. Users will have to change the `run_path` to match the location of the `.../DISTRIBUTE/data` path on their machine. 
3. `attenuationMod_fxns.py` contains all the model code used in the run script. 
4. `accessory_fxns.py` contains some accessory functions separate from the main model.

### `inputs` directory
Contains two files that the model will access to run:
1. `*.nc` file containing the climatology data. Variables required and their names can be found in the `initialize.py` file. Some are not needed for running the model in isotope-only mode. 
2. `*.csv` file with coordinates of the upwind streamline data we want the model to save. Streamline data are memory-intensive, so we don't save the upwind streamline for each grid cell. But any lat, lon pair can be defined in this streamline .csv file – the model will find the closest grid cell to the defined point and save the results (this includes upwind tau plus climatology data).

### `results` directory
Contains the output of running the model. There are two files which correspond to the two input files:
1. `*.nc` file containing the climatology and the associated isotope data
2. `*.csv` file containing the upwind streamline data for each point listed in the input .csv file. 