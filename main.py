# import required libraries
import subprocess
from os.path import isfile, join, isdir
from os import listdir, walk, remove, getenv, mkdir
from rasterstats import zonal_stats
import geopandas as gpd
import pandas as pd
import rasterio
import logging
import random
import string
from pathlib import Path


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
    logger.info('Running read metadata method')
    # search for the metadata.csv file
    metadata_files = find_metadata_files()

    logger.info(f'Found files: {metadata_files}')
    # load in the metadata .csv file
    # assume the first one found is the correct one
    df = pd.read_csv(metadata_files[0], skipinitialspace=True)

    return df


def grid_file_to_12km_rcm(file, method='sum', output_name=None):
    """
    Take a file and re-grid to the 12Km RCM grid

    Inputs:
        - file: path to the raster layer to be processed
        - method: the method to use when converting the values for the new grid size

    Outputs:
        - new file: the path and name of the new file

    """
    logger.info('Running grid to 12km RCM grid')

    # get just the name of the file (remove the extension)
    file_name = file.split('.')[0]
    if output_name is None:
        output_name = f'{file_name}-12km-{method}'

    logger.info(f'Input: {join(data_path, inputs_directory, input_data_directory, file)}')
    logger.info('Output: %s' %(join(data_path, outputs_directory, "%s.asc" % output_name)))

    # re-grid to the 12km RCM grid, save as virtual raster in the temp directory
    subprocess.run(["gdalwarp",
                    "-te", "0", "12000", "660000", "1212000",
                    "-tr", "12000", "12000",
                    "-r", "sum",
                    '-ot', 'Float32',
                    '-co', 'COMPRESS=LZW',
                    '-co', 'NUM_THREADS=ALL_CPUS',
                    #join(data_path, inputs_directory, input_data_directory, file),
                    file,
                    join(data_path, temp_directory, f'{output_name}-temp.vrt')])

    # save the virtual raster as an .asc in the output directory
    subprocess.run(
        ["gdal_translate", "-of", "AAIGrid", join(data_path, temp_directory, f"{output_name}-temp.vrt"),
         join(data_path, outputs_directory, f"{output_name}.asc")])

    return join(data_path, outputs_directory, f"{output_name}.asc")

def rasterise(file, output_name='output.tif', attribute_name='value'):
    """
    Rasterise a vector layer to a tif

    Inputs:
        - file: the path and name of the vector file to be processed
        - output_name: the name to be given to the generated output
        - attribute_name: the attribute in the vector layer containing the values of interest

    Outputs:
        - the name of and path to the output file

    """

    subprocess.call(['gdal_rasterize',
                     #'-burn', '1',  # fixed value to burn for all objects
                     '-a', attribute_name,
                     '-tr', '1000', '1000',  # target resolution <xres> <yres>
                     "-te", "0", "12000", "660000", "1212000",
                     '-co', 'COMPRESS=LZW', '-co', 'NUM_THREADS=ALL_CPUS',  # creation options
                     '-ot', 'Float32',  # output data type
                     join(file),
                     join(data_path, 'outputs', f'{output_name}')])  # src_datasource, dst_filename

    return join(data_path, 'outputs', f'{output_name}')


def add_initial_population(gdf, ssp=1):
    """
    Gets the population in 2020 for the zones (LADs) of interest and adds them to the geodataframe

    Inputs: geodataframe
    Returns: geodatagrame
    """
    logger.info('Running add initial population method')

    # set the looping parameters for the SSPs

    # search for the population file
    input_files = [f for f in listdir(join(data_path, 'inputs', 'population')) if isfile(join(data_path, 'inputs','population', f))]

    # read in the population file into a dataframe
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
    lads = gdf['code'].values.tolist()
    print('Local authority codes to get population data for: ',lads)

    clipped_pop = population[population["code"].isin(lads)]
    clipped_pop.set_index('code')
    gdf.set_index('code')
    gdf = gdf.merge(clipped_pop, on='code', how='inner')

    gdf = gdf.drop(['Lad_Name'], axis=1)
    gdf = gdf.rename(columns={"population_total": "added_population"})

    # Convert population figures to 1000s
    gdf.loc[:, 'initial_population'] *= 1000

    # create a column of the total population in the zone by summing the new and the base/initial
    gdf['population_total'] = gdf['added_population'] + gdf['initial_population']

    # area in km
    gdf['area_h'] = gdf['hectares']/1000.0 # the hectares column is actually meters
    gdf['area_km'] = gdf['area_h'] / 100.0

    # create a population density column
    gdf['population_density'] = gdf['population_total']/gdf['area_km']

    # create population per 1km cell
    gdf['population_1km'] = gdf['population_density']

    print(f'Population stats {gdf.head()}')
    logger.info('Completed add initial population method')
    return gdf


def convert_dph_to_pph(file):
    """

    """
    print('Converting dph to pph')
    logger.info('Converting dph to pph')

    # constant
    people_per_household = 2.5

    #
    dph = rasterio.open(join(data_path, 'inputs', 'layers', 'out_cell_dph.asc'))

    # copy raster so have raster to write to
    pph = dph.read(1)

    i = 0
    while i < pph.shape[0]:
        j = 0
        while j < pph.shape[1]:
            d = pph[i, j]

            # assign new value
            pph[i, j] = d * people_per_household
            j += 1
        i += 1

    new_dataset = rasterio.open(
        join('/data', 'inputs', 'layers', 'out_cell_pph.asc'), "w",
        # driver = "GTiff",
        height=pph.shape[0],
        width=pph.shape[1],
        count=1,
        nodata=-1,
        dtype='float32',
        crs=27700,
        transform=dph.transform,
        compress='lzw'
    )

    new_dataset.write(pph, 1)
    new_dataset.close()
    print('Written new pph layer')
    logger.info('Written new pph layer')

    return


def located_population(file_name=None, data_path='/data/inputs', output_path='/data/outputs', ssp_scenario=None, year=None, zone_id_column='id', total_population=False, fill_northern_ireland=False):
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
    ssp_number = ssp_scenario[3:]

    logger.info('Running fetch population method')

    # looking for which input file should be used to find the new population
    # get a list of the input files to search through
    input_files = [f for f in listdir(join(data_path, 'layers')) if isfile(join(data_path, 'layers', f))]
    file_name = None

    # see if the out cell pph file is present
    for file in input_files:
        if 'cell_pph.asc' in file:
            file_name = 'out_cell_pph.asc'

    # if 1st choice is not present, look for out cell dph
    if file_name is None:
        for file in input_files:
            if 'cell_dph.asc' in file:
                file_name = 'out_cell_dph.asc'

    # if could not find either pph or dph files, return an error message and exit
    if file_name is None:
        msg = 'Error! Could not find either out_cell_dph or out_cell_pph.'
        print(msg)
        logger.info(msg)
        exit(2)
    else:
        # record which file was found and is going to be used
        logger.info(f'After scanning the input files, using {file_name} to get the allocated population')

    # if only the dph file was found, use to create the pph file and reset the file name to this
    if 'dph' in file_name:
        convert_dph_to_pph(file_name)
        file_name = 'out_cell_pph.asc'

    # load in zones
    input_files = [f for f in listdir(join(data_path,'zones')) if isfile(join(data_path,'zones', f))]
    if len(input_files) == 0:
        print('No input zones found')
    print('Zone files:', input_files)
    for file in input_files:
        ext = file.split('.')[-1]
        if ext == 'shp' or ext == 'gpkg':
            zone_file = file
    gdf = gpd.read_file(join(data_path, 'zones', zone_file))

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
    print(f'Population geodataframe: {gdf.head()}')

    if fill_northern_ireland:
        # UDM doesn't give results for NI, but can be filled using the SSP data

        # check if zones for NI are in the zone data
        # if zones, fetch values from population data
        # at this stage just get the new population (baseline - year)

        # check if row has a value of 0 and a code beginning with N
        # fetch values from population data for base year (2020) and current year

        # search for the population file
        input_files = [f for f in listdir(join(data_path, 'population')) if
                       isfile(join(data_path, 'population', f))]

        # read in the population file into a dataframe
        population = pd.read_csv(join(data_path, 'population', input_files[0]),
                                 usecols=['ID', 'LAD19CD', 'LAD19NM', 'Age Class', 'Scenario', '2020', year])

        # Filters the data by chosen SSP and locates the total population per LAD
        population = population.loc[population['Scenario'] == f'{ssp_scenario}']
        population = population.loc[population['Age Class'] == 'Total']
        # set index
        population = population.assign(ID=range(len(population))).set_index('ID')


        # Removes columns that are not needed (class and SSP)
        population = population.drop(['Age Class', 'Scenario'], axis=1)
        population.columns = ['code', 'Lad_Name', 'initial_population', 'year_population']

        # Identify which LADs are of interest and outputs only the population data for those LADs
        lads = gdf['code'].values.tolist()
        ni_lads = []
        for ld in lads:
            if 'N' in ld:
                ni_lads.append(ld)

        print('Local authority codes to get population data for: ', ni_lads)

        clipped_pop = population[population["code"].isin(ni_lads)]
        clipped_pop.set_index('code')
        gdf.set_index('code')

        # calculate the 'new' population
        clipped_pop['population_change'] = (clipped_pop['year_population']*1000) - (clipped_pop['initial_population']*1000)

        # loop through the lads and update the value in the main dataframe
        for lad in ni_lads:
            value_to_assign = clipped_pop.loc[clipped_pop['code'] == lad, 'population_change']
            gdf.loc[gdf.index[gdf['code']==lad].tolist()[0], 'population_total'] = value_to_assign.values[0]

        print('Calculated pop (after):', gdf.head())
        print('NI POP data:', clipped_pop)

    if total_population:
        ## add population to existing LAD
        gdf = add_initial_population(gdf)

    # save output
    gdf.to_file(join(output_path, "population.gpkg"), layer='ssps', driver="GPKG")
    logger.info('Written population gpkg to file - output.gpkg')

    # this is where I add the base population
    # then use Katie's multipliers to adjust populations
    # then add population density column
    # create a population_per_cell value - 1km cells
    # above is then used when rasterising to 1km
    # re-scale to 12km using sum


    logger.info('Completed population method(s)')
    return gdf


def apply_demographic_ratios(gdf, ssp='SSP1', year='2050', output_path='/data/outputs'):
    """

    Inputs:
        - geodataframe
        - ssp: the ssp for the population
        - year: the year of interest
        - output_path: the location to save the new geopackage of data

    Returns:
        - geodataframe
        - path to output
    """
    logger.info('Running apply demographic ratios method')
    print('Running apply demographic ratios method')

    # get list of input files
    input_files = [f for f in listdir(join(data_path, 'inputs', 'population_ratios')) if isfile(join(data_path, 'inputs','population_ratios', f))]

    # read in input file
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
    print('Ratio columns:', ratio_columns)
    logger.info(f'Using the following data for ratios: {ratio_columns}')
    df_cols = ratio_columns
    df_cols.append('LADcode')
    ratios = ratios[df_cols]

    # merge gdf and ratio column
    ratios = ratios.rename(columns={"LADcode": "code"})
    ratios.set_index('code')
    gdf = gdf.merge(ratios, on='code', how='inner')

    # take the ratios and create new columns with ratios applied
    gdf['0-64'] = gdf['population_total'] * gdf[f'{year}_0-64_{ssp}']
    gdf['65-74'] = gdf['population_total'] * gdf[f'{year}_65-74_{ssp}']
    gdf['75-84'] = gdf['population_total'] * gdf[f'{year}_75-84_{ssp}' ]
    gdf['85'] = gdf['population_total'] * gdf[f'{year}_85_{ssp}']

    # create a population density column and a per 1km column
    gdf['0-64_density'] = gdf['0-64'] / gdf['area_km']
    gdf['0-64_1km'] = gdf['0-64_density']# * 100
    gdf['65-74_density'] = gdf['65-74'] / gdf['area_km']
    gdf['65-74_1km'] = gdf['65-74_density']# * 100
    gdf['75-84_density'] = gdf['75-84'] / gdf['area_km']
    gdf['75-84_1km'] = gdf['75-84_density']# * 100
    gdf['85_density'] = gdf['85'] / gdf['area_km']
    gdf['85_1km'] = gdf['85_density']# * 100

    print('Saving output with ratio breakdowns')

    # save output
    gdf.to_file(join(output_path, "population_demographics.gpkg"), layer='ssps', driver="GPKG")
    logger.info('Written population gpkg to file - population_demographics.gpkg')

    logger.info('Completed apply demographic ratios method')

    return gdf, join(output_path, "population_demographics.gpkg")


def house_type_sum():
    """
    Sum the base house type count with the count of new houses. Reads in baseline data for buildings and adds the data from UDM outputs.
    """
    logger.info('Running house type sum method')
    # base data needs to be in RCM grid
    no_data_value = 99999999

    # dwelling types as defined in UDM
    dwelling_types = {'1': 'detached', '2': 'semi-detached', '3': 'terraced', '4': 'flat'}
    house_types = [1, 2, 3, 4]

    # loop through and do a house type at a time
    for type in house_types:
        # get dwelling type text
        type = dwelling_types[str(type)]
        logger.info(f'------ ------ Dwelling type: {type} ------ ------')

        # read in the new dwellings layer
        input_files = [f for f in listdir(join(data_path, 'outputs')) if
                       isfile(join(data_path, 'outputs', f))]
        print('Dwelling files:', input_files)
        dwelling_file = None
        for f in input_files:
            if f'{type}_new-12km.asc' in f:
                dwelling_file = f
                break
        if dwelling_file is None:
            msg = f'ERROR! Could not find the correct dwelling file for {type} in the outputs directory. Looking for file with "{type}_new-12km.asc" in the name.'
            logger.info(msg)
            print(msg)
            exit(2)

        # read in baseline building counts
        input_files = [f for f in listdir(join(data_path, 'inputs', 'base_house_types')) if
                       isfile(join(data_path, 'inputs', 'base_house_types', f))]
        print('Baseline files:', input_files)
        baseline_file = None
        for f in input_files:
            if f'gb-2017-{type}' in f:
                baseline_file = f
                break
        if baseline_file is None:
            msg = f'ERROR! Could not find the correct baseline file for {type} in the base_house_types directory'
            logger.info(msg)
            print(msg)
            exit(2)

        print('Doing house type:', type)
        print('Dwelling file:', dwelling_file)
        print('Baseline file:', baseline_file)
        logger.info(f'------ ------ Dwelling file: {dwelling_file}')
        logger.info(f'------ ------ Baseline file: {baseline_file}')

        # read in base coverage raster
        dwelling_values = rasterio.open(join(data_path, 'outputs', dwelling_file))
        baseline_values = rasterio.open(join(data_path, 'inputs', 'base_house_types', baseline_file))

        # copy raster so have raster to write to
        raster_outdev = baseline_values.read(1)

        baseline = baseline_values.read(1)
        dwellings = dwelling_values.read(1)

        # iterate over raster using baseline
        i = 0
        while i < raster_outdev.shape[0]:
            j = 0
            while j < raster_outdev.shape[1]:
                # get the baseline
                val_baseline = baseline[i, j]
                if val_baseline == no_data_value or val_baseline > no_data_value:
                    val_baseline = 0
                # get the new count
                val_new = dwellings[i,j]

                # assign new value
                raster_outdev[i, j] = val_baseline + val_new

                j += 1
            i += 1

        new_dataset = rasterio.open(
            join('/data', 'outputs', f'dwellings_{type}_total-12km.asc'), "w",
            # driver = "GTiff",
            height=raster_outdev.shape[0],
            width=raster_outdev.shape[1],
            count=1,
            nodata=-1,
            dtype=raster_outdev.dtype,
            crs=27700,
            transform=baseline_values.transform,
            compress='lzw'
        )

        # write new output
        new_dataset.write(raster_outdev, 1)
        new_dataset.close()

        dwelling_file = None
        baseline_file = None
        logger.info('------ ------ Written output: %s' %join('/data', 'outputs', 'dwellings_%s_total-12km.asc' % type))

    return


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
    logger.info('Running the create house type layers method')

    # dwelling types as defined in UDM
    dwelling_types = {'1': 'detached', '2':'semi-detached', '3':'terraced', '4':'flat'}
    house_types = [1,2,3,4]

    input_files = [f for f in listdir(join(data_path, 'inputs', 'layers')) if
                   isfile(join(data_path, 'inputs', 'layers', f))]

    logger.info(f'Searching for files: {input_files}')

    # get file name for raster containing building type
    build_type_file = None
    for f in input_files:
        if 'build_type' in f:
            build_type_file = f
            break
    if build_type_file is None:
        msg = 'ERROR! Could not find out_cell_build_type.asc file in layers directory'
        logger.info(msg)
        print(msg)
        exit(2)

    # get file name for raster containing dwelling density
    dwelling_density_file = None
    for f in input_files:
        if 'out_cell_dph' in f:
            dwelling_density_file = f
            break
    if dwelling_density_file is None:
        msg = 'ERROR! Could not find out_cell_dph.asc file in layers directory'
        logger.info(msg)
        print(msg)
        exit(2)

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
            join('/data', 'outputs', f'dwellings_{dwelling_types[str(type)]}_new.asc'), "w",
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

    # change new houses to 12km grid
    for type in house_types:
        grid_file_to_12km_rcm(file=join('/data', 'outputs', f'dwellings_{dwelling_types[str(type)]}_new.asc'), output_name=f'dwellings_{dwelling_types[str(type)]}_new-12km')

    logger.info('Completed create house type layers method')

    return

# set data path and directory names
data_path = '/data'
inputs_directory = 'inputs'
input_data_directory = 'layers'
temp_directory = 'temp'
outputs_directory = 'outputs'

# check if required folder structure in place
# if so and folders have files in, empty
# temp directory - create if does not exist
if isdir(join(data_path, temp_directory)) is False:
    mkdir(join(data_path, temp_directory))
else:
    files = [f for f in listdir(join(data_path, temp_directory)) if isfile(join(data_path, temp_directory,f))]
    for file in files:
        remove(join(data_path, temp_directory,file))

# input directory - create if does not exist
if isdir(join(data_path, inputs_directory)) is False:
    mkdir(join(data_path, inputs_directory))

# input data directory - create if does not exist
if isdir(join(data_path, inputs_directory, input_data_directory)) is False:
    mkdir(join(data_path, inputs_directory, input_data_directory))

# outputs directory - create if does not exist
if isdir(join(data_path, outputs_directory)) is False:
    mkdir(join(data_path, outputs_directory))

# empty output dir
files = [f for f in listdir(join(data_path, outputs_directory)) if isfile(join(data_path, outputs_directory,f))]
for file in files:
    remove(join(data_path, outputs_directory,file))

# stet up logger and log file
logger = logging.getLogger('udm-heat')
logger.setLevel(logging.INFO)
log_file_name = 'udm-heat-%s.log' %(''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)))
fh = logging.FileHandler( Path(join(data_path, outputs_directory)) / log_file_name)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info('Log file established!')
logger.info('------      ------')

## get passed variables
calc_new_population_total = getenv('calculate_new_population')
if calc_new_population_total is None or str(calc_new_population_total).lower() == 'false':
    calc_new_population_total = False
new_population_demographic_breakdowns = getenv('demographic_breakdown')
if new_population_demographic_breakdowns is None or str(new_population_demographic_breakdowns).lower() == 'false':
    new_population_demographic_breakdowns = False
generate_new_dwelling_totals = getenv('new_dwelling_totals')
if generate_new_dwelling_totals is None or str(generate_new_dwelling_totals).lower() == 'false':
    generate_new_dwelling_totals = False
dwellings_count_total = getenv('dwelling_totals')
if dwellings_count_total is None or str(generate_new_dwelling_totals).lower() == 'false':
    dwellings_count_total = False
rasterise_population_outputs = getenv('rasterise_population_outputs')
if rasterise_population_outputs is None or str(rasterise_population_outputs).lower() == 'false':
    rasterise_population_outputs = False
include_northern_ireland = getenv('include_northern_ireland')
if include_northern_ireland is None or str(include_northern_ireland).lower() == 'false':
    include_northern_ireland = False


logger.info('Fetched passed parameters')
logger.info(f'Calculate new population: {calc_new_population_total}')
logger.info(f'Calculate demographic breakdowns: {new_population_demographic_breakdowns}')
logger.info(f'Rasterise population outputs: {rasterise_population_outputs}')
logger.info(f'Calculate new dwelling totals: {generate_new_dwelling_totals}')
logger.info(f'Calculate total dwellings: {dwellings_count_total}')
logger.info(f'Include Northern Ireland in population outputs: {include_northern_ireland}')

## start the processing
# get list of input files to loop through
files = [f for f in listdir(join(data_path, inputs_directory, input_data_directory)) if isfile(join(data_path, inputs_directory, input_data_directory,f))]
print(f'Files to loop through: {files}')
logger.info('------      ------')
logger.info(f'Got files in input data directory ({files})')

# get the key parameters used for the UDM run
# this should count the year ('YEAR') and SSP ('SSP') as a minimum
parameters_dataframe = read_in_metadata()
# setting index on the parameter column
parameters_dataframe.set_index(['PARAMETER'], inplace=True)
# extract values for ssp and year
ssp = parameters_dataframe.loc['SSP']['VALUE']
year = parameters_dataframe.loc['YEAR']['VALUE']
logger.info('Read in metadata file and extracted key UDM parameter values')
logger.info(f'----- SSP:{ssp}')
logger.info(f'----- Year: {year}')


# calculate the new population
logger.info('------Population data------')
if calc_new_population_total:
    logger.info('Calculating the new population totals')
    gdf = located_population(year=year, ssp_scenario=ssp, total_population=True, fill_northern_ireland=include_northern_ireland )

    if rasterise_population_outputs:
        logger.info('Rasterising population output')
        # rasterise at 1km resolution, then convert to 12km RCM using sum method
        output = grid_file_to_12km_rcm(rasterise(file='/data/outputs/population.gpkg', attribute_name='population_1km'), output_name='population_total-12km')


    if new_population_demographic_breakdowns:
        logger.info('Creating new demographic profiles for new population')
        # create demographic breakdowns for the new populations
        gdf, output = apply_demographic_ratios(gdf, year=year, ssp=ssp)

        if rasterise_population_outputs:
            logger.info('Rasterising demographic population breakdown')
            # need to rasterise per demographic breakdown category
            grid_file_to_12km_rcm(rasterise(file=output, attribute_name='0-64_1km'), output_name='population_demographics_0-64')
            grid_file_to_12km_rcm(rasterise(file=output, attribute_name='65-74_1km'), output_name='population_demographics_65-74')
            grid_file_to_12km_rcm(rasterise(file=output, attribute_name='75-84_1km'), output_name='population_demographics_75-84')
            grid_file_to_12km_rcm(rasterise(file=output, attribute_name='85_1km'), output_name='population_demographics_85')

else:
    logger.info('Skipping population methods')

# calculate the total of new dwellings
logger.info('------Dwelling data------')
if generate_new_dwelling_totals:
    logger.info('Getting the totals for new dwellings')
    create_house_type_layers()

    # get the total number of dwellings, old and new
    if dwellings_count_total:
        logger.info('Calculating the total number of dwellings (old and new)')
        # load in the counts from the base data and add to the new
        house_type_sum()

else:
    logger.info('Skipping dwelling methods')

logger.info('Completed model')


