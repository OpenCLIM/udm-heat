# import required libraries
import subprocess
from os.path import isfile, join, isdir
from os import listdir, walk, remove, getenv, mkdir
from rasterstats import zonal_stats
import geopandas as gpd
import pandas as pd
import rasterio


def find_metadata_files():
    """
    Search all directories for any metadata files (metadat.csv)
    """

    suitable_extension_types = ['csv','']

    files_ = []

    for root, dirs, files in walk('/data/inputs'):
        print(root, files)
        for file in files:

            # check if extension acceptable
            extension = file.split('.')[-1]
            print('checking extension of:', file)
            if extension in suitable_extension_types:
                # check if metadata text in file name
                if 'metadata' in file:
                    # file is good and what we are looking for
                    files_.append(join(root, file))

    print(files_)
    return files_


def read_in_metadata():
    """
    Read in the metadata.csv file
    :return:
    """
    # search for the metadata.csv file
    metadata_files = find_metadata_files()

    # load in the metadata .csv file
    # assume the first one found is the correct one
    df = pd.read_csv(metadata_files[0], skipinitialspace=True)

    return df


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

def rasterise(file):
    """

    :return:
    """

    subprocess.call(['gdal_rasterize',
                     #'-burn', '1',  # fixed value to burn for all objects
                     '-a', 'population_total',
                     '-tr', '1000', '1000',  # target resolution <xres> <yres>
                     "-te", "0", "12000", "660000", "1212000",
                     '-co', 'COMPRESS=LZW', '-co', 'NUM_THREADS=ALL_CPUS',  # creation options
                     '-ot', 'UInt16',  # output data type
                     join(data_path, 'outputs', 'output.gpkg'),
                     join(data_path, 'outputs', 'rasterise.tif')])  # src_datasource, dst_filename
    return

def add_initial_population(gdf):
    """
    Gets the population in 2020 for the zones (LADs) of interest and adds them to the geodataframe
    :return:
    """
    # set the looping parameters for the SSPs
    ssp = 1 # this never changes, population in 2020 always the same

    #
    input_files = [f for f in listdir(join(data_path, 'inputs', 'population')) if isfile(join(data_path, 'inputs','population', f))]

    #population = pd.read_csv(next(inputs.glob('ssp/*.csv')),
    population = pd.read_csv(join(data_path, 'inputs', 'population', input_files[0]),
                             usecols=['ID', 'LAD19CD', 'LAD19NM', 'Age Class', 'Scenario', '2020'])


    # Filters the data by chosen SSP and locates the total population per LAD
    population = population.loc[population['Scenario'] == f'SSP{ssp}']
    population = population.loc[population['Age Class'] == 'Total']
    population = population.assign(ID=range(len(population))).set_index('ID')

    # Removes columns that are no needed (class and SSP)
    population = population.drop(['Age Class', 'Scenario'], axis=1)
    population.columns = ['code', 'Lad_Name', 'initial_population']

    # Identify which LADs are of interest and outputs only the population data for those LADs
    print(gdf.head())
    lads = gdf['code'].values.tolist()
    print(lads)


    clipped_pop = population[population["code"].isin(lads)]
    print(clipped_pop)
    clipped_pop.set_index('code')
    gdf.set_index('code')
    gdf = gdf.merge(clipped_pop, on='code', how='inner')

    gdf = gdf.drop(['Lad_Name'], axis=1)
    gdf = gdf.rename(columns={"population_total": "added_population"})

    # Convert population figures to 1000s
    gdf.loc[:, 'initial_population'] *= 1000

    # create a column of the total population in the zone by summing the new and the base/initial
    gdf['population_total'] = gdf['added_population'] + gdf['initial_population']


    # create a population density column
    gdf['population_density'] = gdf['population_total']/gdf['hectares']
    print(gdf)
    print(gdf.columns)

    return gdf

def located_population(file_name='out_cell_pph.asc', data_path='/data/inputs', output_path='/data/outputs', ssp_scenario=None, year=None, zone_id_column='id', total_population=False):
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
    zone_id_column = 'code'

    # load in zones
    input_files = [f for f in listdir(join(data_path,'zones')) if isfile(join(data_path,'zones', f))]
    if len(input_files) == 0:
        print('No input zones found')
    print('Zone files:', input_files)
    gdf = gpd.read_file(join(data_path, 'zones', input_files[0]))

    # run the zonal stats
    # get the number of new people per cell (doesn't include existing people)
    stats = zonal_stats(gdf, join(data_path, 'layers', file_name), stats=['sum'])

    print('Stats:', stats)

    # get a list of zone IDs
    zone_ids = []
    for index, row in gdf.iterrows():
        #print(row[zone_id_column])
        zone_ids.append(row[zone_id_column])
    print('Zone ids:', zone_ids)
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
    print('Pop data:', pop_series)

    # create a column containing the population in the new cells of development for each zone
    gdf['population_total'] = pop_series

    if total_population:
        gdf = total_population(gdf, ssp_scenario)

    # save output
    gdf.to_file(join(output_path, "output.gpkg"), layer='ssps', driver="GPKG")

    # this is where I add the base population
    # then use Katie's multipliers to adjust populations
    # then add population density column
    # create a population_per_cell value - 1km cells
    # above is then used when rasterising to 1km
    # re-scale to 12km using sum


    ## add population to existing LAD
    gdf = add_initial_population(gdf)


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


# get list of input files to loop through
files = [f for f in listdir(join(data_path, inputs_directory, input_data_directory)) if isfile(join(data_path, inputs_directory, input_data_directory,f))]
print('files to loop through: %s' %files)


# get the key parameters used for the UDM run
# this should count the year ('YEAR') and SSP ('SSP') as a minimum
parameters_dataframe = read_in_metadata()
# setting index on the parameter column
parameters_dataframe.set_index(['PARAMETER'], inplace=True)
# extract values for ssp and year
ssp = parameters_dataframe.loc['SSP']['VALUE']
year = parameters_dataframe.loc['YEAR']['VALUE']

# loop through the files and process
for file in files:
    print('Processing: %s' %file)
    grid_file_to_12km_rcm(file)

located_population()
rasterise('/data/outputs/output.gpkg')