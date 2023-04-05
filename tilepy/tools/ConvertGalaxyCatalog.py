import pandas as pd
import numpy as np
import tqdm
import tables
import logging
import argparse
import time

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

# Define argument parser
parser = argparse.ArgumentParser(description='Tool to converte GLADE+ catalog to internal data format')
parser.add_argument("-i", "--input", nargs='?',
                    dest='input', action='store',
                    required=True, help="Input catalog file")
parser.add_argument("-o", "--output", nargs='?',
                    dest='output', action='store',
                    required=True, help="Output converted file")
parser.add_argument("-ca", "--compression-algorithm", nargs='?',
                    dest='compression_algorithm', action='store',
                    required=False, default=None, help="Algorithm used for compression")
parser.add_argument("-cl", "--compression-level", nargs='?',
                    dest='compression_level', action='store',
                    required=False, default=0, type=int, choices=range(0, 10),
                    help="Compression level")
parser.add_argument("-md", "--max-luminosity-distance", nargs='?',
                    dest='max_luminosity_distance', action='store',
                    required=False, default=500, type=float, help="Cuts to apply on the luminosity distance")
parser.add_argument("-sv", "--store-valid",
                    dest='store_valid', action='store_true',
                    default=False, help="Store only valid galaxies")
parser.add_argument("-t", "--text-format",
                    dest='text_format', action='store_true',
                    default=False, help="Store as a text format")


# Function to determine the best name for the galaxie
def best_name_galaxie(x):
    if x.name_SDSS_DR16Q != 'null':
        return 'SDSS ' + x.name_SDSS_DR16Q
    elif x.name_GWGC != 'null':
        return x.name_GWGC
    elif x.name_HyperLEDA != 'null':
        return x.name_HyperLEDA
    elif x.no_PGC != 'null':
        return 'PGC ' + x.no_PGC
    elif x.name_2MASS != 'null':
        return '2MASS ' + x.no_PGC
    elif x.name_WISExSCOS != 'null':
        return 'WISE ' + x.name_WISExSCOS
    else:
        return 'GLADE+ ' + str(x.no_GLADE)


# Function to convert flags from the catalog to integer
def flag_columns_converter(flag):
    if flag == 'null':
        return -1
    else:
        return int(flag)


# Define columns parameters for catalogs
name_columns_catalog = ['no_GLADE', 'no_PGC', 'name_GWGC', 'name_HyperLEDA', 'name_2MASS', 'name_WISExSCOS',
                        'name_SDSS_DR16Q', 'object_type', 'RA', 'Dec', 'B_mag', 'B_mag_err', 'B_mag_flag',
                        'B_mag_abs', 'J_mag', 'J_mag_err', 'H_mag', 'H_mag_err', 'K_mag', 'K_mag_err',
                        'W1_mag', 'W1_mag_err', 'W2_mag', 'W2_mag_err', 'W1_mag_flag', 'B_J_mag', 'B_J_mag_err',
                        'z_helio', 'z_cmb', 'z_flag', 'v_err', 'z_err', 'd_L', 'd_L_err', 'dist_flag', 'mass',
                        'mass_err', 'mass_flag', 'merger_rate', 'merger_rate_err']
dtype_columns_catalog = {'no_GLADE': int, 'RA': float, 'Dec': float,
                         'B_mag': float, 'B_mag_err': float, 'B_mag_abs': float, 'J_mag': float,
                         'J_mag_err': float, 'H_mag': float, 'H_mag_err': float, 'K_mag': float, 'K_mag_err': float,
                         'W1_mag': float, 'W1_mag_err': float, 'W2_mag': float, 'W2_mag_err': float,
                         'B_J_mag': float, 'B_J_mag_err': float, 'z_helio': float, 'z_cmb': float,
                         'v_err': float, 'z_err': float, 'd_L': float, 'd_L_err': float, 'mass': float,
                         'mass_err': float, 'merger_rate': float, 'merger_rate_err': float}
converter_columns_catalog = {'no_PGC': lambda x: str(x),
                             'name_GWGC': lambda x: str(x),
                             'name_HyperLEDA': lambda x: str(x),
                             'name_2MASS': lambda x: str(x),
                             'name_WISExSCOS': lambda x: str(x),
                             'name_SDSS_DR16Q': lambda x: str(x),
                             'object_type': lambda x: str(x),
                             'B_mag_flag': flag_columns_converter,
                             'W1_mag_flag': flag_columns_converter,
                             'z_flag': flag_columns_converter,
                             'dist_flag': flag_columns_converter,
                             'mass_flag': flag_columns_converter}


# Define table structure hdf5
class TableCatalog(tables.IsDescription):
    no_GLADE = tables.UInt32Col()
    valid_data = tables.BoolCol()
    name = tables.StringCol(30)
    RA = tables.Float64Col()
    Dec = tables.Float64Col()
    d_L = tables.Float64Col()
    d_L_err = tables.Float64Col()
    dist_flag = tables.Int8Col()
    mass = tables.Float64Col()
    mass_err = tables.Float64Col()
    mass_flag = tables.Int8Col()
    merger_rate = tables.Float64Col()
    merger_rate_err = tables.Float64Col()
    no_PGC = tables.StringCol(7)
    name_GWGC = tables.StringCol(28)
    name_HyperLEDA = tables.StringCol(29)
    name_2MASS = tables.StringCol(16)
    name_WISExSCOS = tables.StringCol(19)
    name_SDSS_DR16Q = tables.StringCol(18)
    object_type = tables.StringCol(1)
    B_mag = tables.Float64Col()
    B_mag_err = tables.Float64Col()
    B_mag_flag = tables.Int8Col()
    B_mag_abs = tables.Float64Col()
    J_mag = tables.Float64Col()
    J_mag_err = tables.Float64Col()
    H_mag = tables.Float64Col()
    H_mag_err = tables.Float64Col()
    K_mag = tables.Float64Col()
    K_mag_err = tables.Float64Col()
    W1_mag = tables.Float64Col()
    W1_mag_err = tables.Float64Col()
    W2_mag = tables.Float64Col()
    W2_mag_err = tables.Float64Col()
    W1_mag_flag = tables.Int8Col()
    B_J_mag = tables.Float64Col()
    B_J_mag_err = tables.Float64Col()
    z_helio = tables.Float64Col()
    z_cmb = tables.Float64Col()
    z_flag = tables.Int8Col()
    v_err = tables.Float64Col()
    z_err = tables.Float64Col()


# Load argument
args = parser.parse_args()

# Load catalog
logging.info('Start loading catalog file')
tstart = time.time()
catalog = pd.read_csv(args.input, sep=' ', header=None,
                      names=name_columns_catalog, dtype=dtype_columns_catalog, converters=converter_columns_catalog)
logging.info('Catalog files loaded with success in {0}s'.format(time.time() - tstart))

# Compute filter
logging.info('Start computing catalog filter')
tstart = time.time()

logging.info('Determine data validity')
# Remove galaxies with invalid distance or too far
catalog['valid_data'] = True
catalog['valid_data'] &= ~np.isnan(catalog['d_L'])
catalog['valid_data'] &= catalog['d_L'] < args.max_luminosity_distance

# Remove non valid data if requested
if args.store_valid:
    logging.info('Remove non valid data')
    catalog = catalog[catalog['valid_data']]

# Determine galaxies best names
logging.info('Determine best name for galaxies')
# catalog['name'] = catalog.apply(best_name_galaxie, axis=1)
catalog['name'] = catalog['no_GLADE'].apply(lambda x: 'GLADE+ ' + str(x))

logging.info('Catalog filter computed in {0}s'.format(time.time() - tstart))

# Create output files
logging.info('Start writing file')
tstart = time.time()

if args.text_format:
    catalog = catalog[catalog['valid_data']]
    catalog[['no_GLADE', 'RA', 'Dec', 'd_L', 'mass']].to_csv(args.output, index=False)
else:
    # if args.compression_algorithm is None:
    #    catalog.to_hdf(args.output, 'catalog', mode='w', index=False)
    # else:
    #    catalog.to_hdf(args.output, 'catalog', mode='w', index=False,
    #                   complib=args.compression_algorithm, complevel=args.compression_level)
    if args.compression_algorithm is None:
        filter_catalog = tables.Filters()
    else:
        filter_catalog = tables.Filters(complib=args.compression_algorithm, complevel=args.compression_level)

    h5file = tables.open_file(args.output, mode='w', filters=filter_catalog)
    table = h5file.create_table(h5file.root, 'catalog', TableCatalog, expectedrows=len(catalog))
    #table.append(catalog)
    entry = table.row
    for i in tqdm.tqdm(catalog.index):
        for key in catalog.columns:
            entry[key] = catalog.loc[i, key]
        entry.append()
    h5file.close()


logging.info('File written in {0}s'.format(time.time() - tstart))
