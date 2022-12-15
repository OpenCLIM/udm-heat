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


    return gdf


def apply_demographic_rations(gdf, ssp='SSP1', year='2050'):
    """

    :param gdf:
    :return:
    """
    input_files = [f for f in listdir(join(data_path, 'inputs', 'population_ratios')) if isfile(join(data_path, 'inputs','population_ratios', f))]

    ratios = pd.read_csv(join(data_path, 'inputs', 'population_ratios', input_files[0]))#,
                             #usecols=['ID', 'LAD19CD', 'LAD19NM', 'Age Class', 'Scenario', '2020'])

    # get the list of zones of interest
    lads = gdf['code'].values.tolist()

    # filter the ratios to just the zones of interest
    ratios = ratios[ratios["LADcode"].isin(lads)]
    print(ratios.head())

    # filter the columns in the ratios to just the SSP of interest
    ratio_columns = []
    for col in ratios.columns:
        if ssp in col:
            if year in col:
                ratio_columns.append(col)

    # filter the ratios df to just the columns of interest
    print(ratio_columns)
    df_cols = ratio_columns
    df_cols.append('LADcode')
    ratios = ratios[df_cols]

    # merge gdf and ratio column
    ratios = ratios.rename(columns={"LADcode": "code"})
    ratios.set_index('code')
    gdf = gdf.merge(ratios, on='code', how='inner')

    # take the ratios and create new columns with ratios applied
    gdf['0-64'] = gdf['population_total'] * gdf['%s_%s_0-64'%(ssp, year)]
    gdf['65-74'] = gdf['population_total'] * gdf['%s_%s_65-74' % (ssp, year)]
    gdf['75-84'] = gdf['population_total'] * gdf['%s_%s_75-84' % (ssp, year)]
    gdf['85'] = gdf['population_total'] * gdf['%s_%s_85' % (ssp, year)]

    print(gdf.head())
    print(gdf.columns)

    return gdf


def create_house_type_layers():
    """
    Take the out_cell_build_type.asc from UDM and generate layers for each building type

    .asc from UDM has 4 values of interest. This just records the type of dwelling. Need to look at this and dwellings per hectare output (out_cell_dph.asc) to form a count of the number of dwellings in a cell
    - 1 = detached
    - 2 = semi-detached
    - 3 = terraced
    - 4 = flats

    :return:
    """
    # dwelling types as defined in UDM
    dwelling_types = {'1': 'detached', '2':'semi-detached', '3':'terraced', '4':'flat'}
    house_types = [1,2,3,4]

    input_files = [f for f in listdir(join(data_path, 'inputs', 'layers')) if
                   isfile(join(data_path, 'inputs', 'layers', f))]

    # get file name for raster containing building type
    for f in input_files:
        if 'build_type' in f:
            build_type_file = f
            break

    # get file name for raster containing dwelling density
    for f in input_files:
        if 'out_cell_dph' in f:
            dwelling_density_file = f
            break

    # read in base coverage raster
    building_types = rasterio.open(join(data_path,'inputs', 'layers', build_type_file))
    dwellings_per_hectare = rasterio.open(join(data_path,'inputs', 'layers', dwelling_density_file))

    # loop through and do a house type at a time
    for type in house_types:
        i = 0

        # copy raster so have raster to write to
        raster_outdev = building_types.read(1)

        dwellings = dwellings_per_hectare.read(1)

        while i < raster_outdev.shape[0]:
            j = 0
            while j < raster_outdev.shape[1]:
                cv = raster_outdev[i, j]
                if cv == type:
                    # get the density for the tile
                    density = dwellings[i, j]
                    # assign new value
                    raster_outdev[i, j] = density
                else:
                    raster_outdev[i,j] = 0
                j += 1
            i += 1

        new_dataset = rasterio.open(
            join('/data', 'outputs', 'dwellings_%s.asc'%dwelling_types[str(type)]), "w",
            # driver = "GTiff",
            height=raster_outdev.shape[0],
            width=raster_outdev.shape[1],
            count=1,
            nodata=-1,
            dtype=raster_outdev.dtype,
            crs=27700,
            transform=building_types.transform,
            compress='lzw'
        )

        # Z = raster_file.read(1)
        new_dataset.write(raster_outdev, 1)
        new_dataset.close()

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


## get passed variables
calc_new_population_total = getenv('calculate_new_population')
if calc_new_population_total is None:
    calc_new_population_total = False
new_population_demographic_breakdowns = getenv('demographic_breakdown')
if new_population_demographic_breakdowns is None:
    new_population_demographic_breakdowns = False
generate_new_dwelling_totals = getenv('new_dwelling_totals')
if generate_new_dwelling_totals is None:
    generate_new_dwelling_totals = False
dwellings_count_total = getenv('dwelling_totals')
if dwellings_count_total is None:
    dwellings_count_total = False


## start the processing
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
#for file in files:
#    print('Processing: %s' %file)
#    grid_file_to_12km_rcm(file)


# calculate the new population
if calc_new_population_total:
    gdf = located_population()

    if new_population_demographic_breakdowns:
        # create demographic breakdowns for the new populations
        apply_demographic_rations(gdf)

# calculate the total of new dwellings
if generate_new_dwelling_totals:
    create_house_type_layers()

    # get the total number of dwellings, old and new
    if dwellings_count_total:
        pass

#rasterise('/data/outputs/output.gpkg')

