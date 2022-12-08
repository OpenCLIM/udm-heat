# import required libraries
import subprocess
from os.path import isfile, join, isdir
from os import listdir, remove, getenv, mkdir
from rasterstats import zonal_stats
import geopandas as gpd
import pandas as pd
import rasterio

def grid_file_to_12km_rcm(file):
    """
    Take a file and re-grid to the 12Km RCM grid
    """

    # get just the name of the file (remove the extension)
    file_name = file.split('.')[0]

    # re-grid to the 12km RCM grid, save as virtual raster in the temp directory
    subprocess.run(["gdalwarp", "-te", "0", "12000", "660000", "1212000", "-tr", "12000", "12000", "-r", "sum",
                    join(data_path, inputs_directory, input_data_directory, file),
                    join(data_path, temp_directory, '%s-%s.vrt' % (file_name, 'temp'))])

    # save the virtual raster as an .asc in the output directory
    subprocess.run(
        ["gdal_translate", "-of", "AAIGrid", join(data_path, temp_directory, "%s-%s.vrt" % (file_name, 'temp')),
         join(data_path, outputs_directory, "%s-12km-sum.asc" % file_name)])


def located_population(data_path, output_path, ssp_scenario, year, zone_id_column='id', total_population=False):
    """
    Uses zonal statistics to get the population in the newly developed cells per zone definition. Optional parameter to then also calculate the total population in the zone.

    :param data_path:
    :param output_path:
    :param ssp_scenario:
    :param year:
    :param total_population:
    :return:
    """
    # set parameters
    fz = False
    zone_id_column = 'DataZone'

    # load in zones
    input_files = [f for f in listdir(data_path) if isfile(join(data_path))]
    gdf = gpd.read_file(join(data_path, 'zones', input_files[0]))

    # get the name of the pph file
    if fz:
        file_name = 'out_cell_pph_ssp%s_fz_%s.asc' % (ssp_scenario, year)
    else:
        file_name = 'out_cell_pph_ssp%s_%s.asc' % (ssp_scenario, year)

    # run the zonal stats
    # get the number of new people per cell (doesn't include existing people)
    stats = zonal_stats(gdf, join(data_path, file_name), stats=['sum'])

    # get a list of zone IDs
    zone_ids = []
    for index, row in gdf.iterrows():
        #print(row[zone_id_column])
        zone_ids.append(row[zone_id_column])

    # assign zone ids to the stats
    zone_stats = {}
    pop_list = []
    index = []
    i = 0
    for stat in stats:
        zone_stats[zone_ids[i]] = stat
        pop_list.append(stat['sum'])
        index.append(i)
        i += 1

    # convert new pop to a series
    pop_series = pd.Series(pop_list, index=index)

    # create a column containing the population in the new cells of development for each zone
    gdf['POP.UDM.%s.%s' % (year, ssp_scenario)] = pop_series

    if total_population:
        # set the looping parameters for the SSPs
        decades = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100]
        ssps = [1, 2, 3, 4, 5]

        # remove data for the other SSPs
        for ssp_ in ssps:
            for decade in decades:
                if ssp_ != ssp_scenario:
                    col = 'POP.' + str(decade) + '.' + str(ssp_)
                    gdf = gdf.drop(columns=[col])

        # create a column of the total population in the zone by summing the new and the base/initial
        gdf['POP.UDM.TOTAL.%s.%s' % (year, ssp_scenario)] = gdf['POP.UDM.%s.%s' % (year, ssp_scenario)] + gdf[
            'POP.%s.%s' % (2020, ssp_scenario)]

    # remove other ssp columns
    for decade in decades:
        if year != decade:
            col = 'POP.' + str(decade) + '.' + str(ssp_scenario)
            gdf = gdf.drop(columns=[col])

    # save output
    if fz:
        gdf.to_file(join(output_path, "ffe-zones-clyde-ssps-udm-%s-fz-%s.gpkg" % (ssp_scenario, year)), layer='ssps', driver="GPKG")
    else:
        gdf.to_file(join(output_path, "ffe-zones-clyde-ssps-udm-%s-%s.gpkg" % (ssp_scenario, year)), layer='ssps', driver="GPKG")

    return


# set data path and directory names
data_path = '/data'
inputs_directory = 'inputs'
input_data_directory = 'layers'
temp_directory = 'temp'
outputs_directory = 'outputs'


# check if required folder structure in place
# if so and folders have files in, empty
## temp
if isdir(join(data_path, temp_directory)) is False:
    mkdir(join(data_path, temp_directory))
else:
    files = [f for f in listdir(join(data_path, temp_directory)) if isfile(join(data_path, temp_directory,f))]
    for file in files:
        remove(join(data_path, temp_directory,file))

## outputs
if isdir(join(data_path, outputs_directory)) is False:
    mkdir(join(data_path, outputs_directory))

# empty output dir
files = [f for f in listdir(join(data_path, outputs_directory)) if isfile(join(data_path, outputs_directory,f))]
for file in files:
    remove(join(data_path, outputs_directory,file))

#get list of input files to loop through
files = [f for f in listdir(join(data_path, inputs_directory, input_data_directory)) if isfile(join(data_path, inputs_directory, input_data_directory,f))]
print('files to loop through: %s' %files)

# loop through the files and process
for file in files:
    print('Processing: %s' %file)
    grid_file_to_12km_rcm(file)
