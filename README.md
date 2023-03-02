# udm-heat
Convert UDM outputs for the HEAT-HARM workflow

## Description
The model performs a number of functions to convert UDM outputs into more derived datasets. Most notable are the methods to generate the new population at zone scale or 12km grid, and generate the number of new house types. More details on the methods employed can be found at the bottom of this page.

## Inputs variables
* calculate_new_population
  * Description: Calculate the population which has been housed
  * Type: boolean
  * Default: true
  * Output:
* demographic_breakdown
  * Description: Apply demographic ratios (passed via a file) to the population values 
  * Type: boolean
  * Default: true
  * Output:
* new_dwelling_totals
  * Description: Calculate the number of new dwellings of the 4 main types (detached, semi-detached, terraced and flat)
  * Type: boolean
  * Default: true
  * Output:
* dwelling_totals
  * Description: The total number of dwellings by type, existing and new
  * Type: boolean
  * Default: true
  * Output:
* rcm_population_outputs
  * Description: Convert the calculated population changes to raster outputs (12km, RCM grid)
  * Type: boolean
  * Default: true
  * Output:
* include_northern_ireland
  * Description: Include population values for Northern Ireland filling the gap from UDM which doesn't include NI data.
  * Type: boolean
  * Default: True
  * Output:

## Input files (data slots)
* layers
  * Description: This should be the files output by the UDM model. As a minimum it should include the following files: metadata.csv, out_cell_dph.asc, out_cell_build_type.asc
  * Location: /data/inputs/layers
* population
  * Description: A .csv of population data, broken down by zone and year.
  * Location: /data/inputs/population
* population_ratios
  * Description: A .csv file of demographic ratios for 4 age bands (0-64,65-74,75-84,84) broken down by zone code and SSP scenario. Column headings should take the form "SSP_YEAR_AGE-BAND".
  * Location: /data/inputs/population_ratios
* zones
  * Description: A spatial file containing the zones being used. These should contain codes. .shp or .gpkg expected.
  * Location: /data/inputs/zones
* base_house_types
  * Description: The baseline count of dwelling types on a 12km RCM grid
  * Location: /data/inputs/base_house_types


## Outputs
Outputs vary depending on passed parameters. Below summarises the potential set of outputs which can be generated with some details. A logfile is also generated. 
### Population
* (1) A geopackage of the 'new' population by local authority based on a UDM output
  * Name: population.gpkg 
  * Location: /data/outputs
* (2) Raster layer of the new population based on (1)
  * Name: population_total-12km.asc
  * Location: /data/outputs
* (3) A geopackage of the 'new' population based on (1) with additional demographic breakdowns
  * Name: population_demographics.gpkg
  * Location: /data/outputs
* (4) Raster layers per demographic category based on (3)
  * Name: population_demographics_#age band#.asc
  * Location: /data/outputs

### Dwellings
* (1) Raster layers of the counts of the 'new' dwellings based on a UDM output
  * Name: dwellings_#dwelling type#_new-12km.asc
  * Location: /data/outputs
* (2) Raster layers of the counts of the total dwellings based on a baseline and the 'new' dwellings based on a UDM output
  * Name: dwellings_#dwelling type#_total-12km.asc
  * Location: /data/outputs
   
### Other
* (1) A logfile for the processes undertaken
  * Name: udm-heat-#random code#.txt
  * Location: /data/outputs

## Running the model
### Docker image
For all image name is set to udm-heat, though this can be changed based on personal preference. The location of the folder with data is assumed to be called 'data', though if this is not the case on your personal system, this can be changed in the docker run command so long as the directory after the ':' remains unchanged('/data').
#### Build docker image
* docker build . -t udm-heat
#### Run for all outputs
* docker run --env calculate_new_population=True --env demographic_breakdown=True --env rasterise_population_outputs=True --env new_dwelling_totals=True --env dwelling_totals=True -v $PWD/data:/data -t udm-heat
#### Run for just population outputs
* docker run --env calculate_new_population=True -v $PWD/data:/data -t udm-heat
#### Run for population outputs and demographic breakdowns in raster format
* docker run --env calculate_new_population=True --env demographic_breakdown=True --env rasterise_population_outputs=True -v $PWD/data:/data -t udm-heat
#### Run for just dwelling outputs
* docker run --env new_dwelling_totals=True --env dwelling_totals=True -v $PWD/data:/data - udm-heat

## Processing methods
### Population calculations
#### Calculating the 'new' population
Using the pph (people per hectare) output from UDM the 'new' population can be found using the 1km gridded output. Where the pph file can't be found, the dph (density per hectare) is used, a value of the number of households per hectare and multiplied by 2.5 (the average number of people per household).

#### Calculating the total population
To derive the total population (baseline + 'new') the original SSP data for 2020 (data for SSP1 is used but all SSPs have the same baseline population in 2020) is found and used as the baseline values. Using the gridded values which are at 1km resolution, the values for the baseline population can be added straight to the new population from the UDM outputs.

#### Calculating the demographic breakdowns
The demographic breakdowns are calculated using pre-calculated ratios from the SSPs, on the same 1km grid as the population data. These are then applied to the total population values calculated as above.

#### Rescale to 12km RCM grid
The population data is calculated at a 1km resolution in raster form. To enable more easy analysis with data in the climate domain which can be in the climate RCM 12km grid, enabling this method will rasterise the 1km OS Grid population data to the climate 12km RCM grid.

#### Including population data for Northern Ireland
The outputs generated by UDM don't cover Northern Ireland, the SSP data is used to populate values for the zones which make up Northern Ireland. There values are then included in all population based outputs, but don't affect the dwelling results. This is enabled by default, but to remove, set the include_northern_ireland input variable to False.

### Dwelling calculations
It's important to note that this is not the number of buildings, but instead the number of dwellings (or addresses) of a building type. For example for flats a value of 100 means there are 100 flats in the cell, but the number of buildings these are contained within is not represented, so this could mean there are 10 buildings, or equally could be 1 building containing the 100 flats.

#### Calculating the number of new dwellings
Using the out_cell_build_type.asc from UDM which gives the type of building built in each cell and combining that with out_cell_dph.asc layer, a number per cell of the number of dwellings in that cell can be found. This can then be output as by house type as seperate layers giving a count of the number of houses of each type per cell.

#### Calculating the total number of dwellings
Using a baseline of 2017 for buildings the number of buildings by type (flat, detached and so on) are added to the 'new' buildings, as calculated above. The baseline data is prepared at the 12km scale on the RCM grid so can be added to the 'new' dwellings directly.
