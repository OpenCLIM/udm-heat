# udm-heat
Convert UDM outputs for the HEAT-HARM workflow

## Description
The model performs a number of functions to convert UDM outputs into more derived datasets. Most notable are the methods to generate the new population at zone scale or 12km grid, and generate the number of new house types.

## Inputs variables
* calculate_new_population
  * Description: Calculate the population which has been housed
  * Type: boolean
  * Default: false
* demographic_breakdown
  * Description: Apply demographic ratios (passed via a file) to the population values 
  * Type: boolean
  * Default: false
* new_dwelling_totals
  * Description: Calculate the number of new dwellings of the 4 main types (detached, semi-detached, terraced and flat)
  * Type: boolean
  * Default: false

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


