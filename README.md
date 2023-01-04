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
* rasterise_population_outputs
  * Description: Convert the calculated population changes to raster outputs (12km, RCM grid)
  * Type: boolean
  * Default: true
  * Output:

## Input files
* layers
  * Description: This should be the files output by the UDM model
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

## Outputs
Outputs vary depending on functions run. Broadly speaking 


## Processing methods
### Population calculations
#### Calculating the 'new' population
Using the pph (people per hectare) output from UDM the 'new' population can be found using the 1km gridded output. Where the pph file can't be found, the dph (density per hectare) is used, a value of the number of households per hectare and multiplied by 2.5 (the average number of people per household).

#### Calculating the total population
To derive the total population (baseline + 'new') the original SSP data for 2020 (data for SSP1 is used but all SSPs have the same baseline population in 2020) is found and used as the baseline values. Using zonal statistics the new pph data is assigned to the appropriate local authority area and thus added to the baseline data from the SSPs which is also by local authority area. The population density is calculated.

#### Calculating the demographic breakdowns
The demographic breakdowns are calculated using pre-calculated ratios from the SSPs per LAD, per demographic category. These are then applied to the total population values calculated as above.

#### Raster outputs
The population data is calculated in vector form by local authority area; to output this as vector data, the data is rasterised at 1km scale initially giving population values per each 1km cell. This is then scaled to the 12km RCM grid by summing the values of all those cells that fall within each 12km grid cell.

### Dwelling calculations
It's important to note that this is not the number of buildings, but instead the number of dwellings (or addresses) of a building type. For example for flats a value of 100 means there are 100 flats in the cell, but the number of buildings these are contained within is not represented, so this could mean there are 10 buildings, or equally could be 1 building containing the 100 flats.

#### Calculating the number of new dwellings
Using the out_cell_build_type.asc from UDM which gives the type of building built in each cell and combining that with out_cell_dph.asc layer, a number per cell of the number of dwellings in that cell can be found. This can then be output as by house type as seperate layers giving a count of the number of houses of each type per cell.

#### Calculating the total number of dwellings
Using a baseline of 2017 for buildings the number of buildings by type (flat, detached and so on) are used already at the 12km scale.