import argparse
import logging
import os.path
import time

import numpy as np
import pandas as pd
import tables

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

# Define argument parser
parser = argparse.ArgumentParser(
    description="Tool to convert GLADE+ catalog into internal data format"
)
parser.add_argument(
    "-i",
    "--input",
    nargs="?",
    dest="input",
    action="store",
    required=True,
    help="Input catalog file",
)
parser.add_argument(
    "-o",
    "--output",
    nargs="?",
    dest="output",
    action="store",
    required=True,
    help="Output converted file",
)
parser.add_argument(
    "-md",
    "--max-luminosity-distance",
    nargs="?",
    dest="max_luminosity_distance",
    action="store",
    required=False,
    default=500,
    type=float,
    help="Cuts to apply on the luminosity distance",
)
parser.add_argument(
    "-sa",
    "--store-all",
    dest="store_all",
    action="store_true",
    default=False,
    help="Store all galaxies included invalid ones",
)
parser.add_argument(
    "-t",
    "--text-format",
    dest="text_format",
    action="store_true",
    default=False,
    help="Store as a text format",
)
parser.add_argument(
    "-ca",
    "--compression-algorithm",
    nargs="?",
    dest="compression_algorithm",
    action="store",
    required=False,
    default=None,
    help="Algorithm used for compression",
)
parser.add_argument(
    "-cl",
    "--compression-level",
    nargs="?",
    dest="compression_level",
    action="store",
    required=False,
    default=0,
    type=int,
    choices=range(0, 10),
    help="Compression level",
)
parser.add_argument(
    "-sh",
    "--shuffle",
    dest="shuffle",
    action="store_true",
    default=False,
    help="Allow shuffle for file compression",
)
parser.add_argument(
    "-bs",
    "--bit-shuffle",
    dest="bit_shuffle",
    action="store_true",
    default=False,
    help="Allow bit shuffle for file compression",
)


# Function to determine the best name for the galaxie
def best_name_galaxie(x):
    if x.name_SDSS_DR16Q != "null":
        return "SDSS " + x.name_SDSS_DR16Q
    elif x.name_GWGC != "null":
        return x.name_GWGC
    elif x.name_HyperLEDA != "null":
        return x.name_HyperLEDA
    elif x.no_PGC != "null":
        return "PGC " + x.no_PGC
    elif x.name_2MASS != "null":
        return "2MASS " + x.no_PGC
    elif x.name_WISExSCOS != "null":
        return "WISE " + x.name_WISExSCOS
    else:
        return "GLADE+ " + str(x.no_GLADE)


# Function to convert flags from the catalog to integer
def flag_columns_converter(flag):
    if flag == "null":
        return -1
    else:
        return int(flag)


# Define columns parameters for catalogs
name_columns_catalog = [
    "no_GLADE",
    "no_PGC",
    "name_GWGC",
    "name_HyperLEDA",
    "name_2MASS",
    "name_WISExSCOS",
    "name_SDSS_DR16Q",
    "object_type",
    "RA",
    "Dec",
    "B_mag",
    "B_mag_err",
    "B_mag_flag",
    "B_mag_abs",
    "J_mag",
    "J_mag_err",
    "H_mag",
    "H_mag_err",
    "K_mag",
    "K_mag_err",
    "W1_mag",
    "W1_mag_err",
    "W2_mag",
    "W2_mag_err",
    "W1_mag_flag",
    "B_J_mag",
    "B_J_mag_err",
    "z_helio",
    "z_cmb",
    "z_flag",
    "v_err",
    "z_err",
    "d_L",
    "d_L_err",
    "dist_flag",
    "mass",
    "mass_err",
    "mass_flag",
    "merger_rate",
    "merger_rate_err",
]
dtype_columns_catalog = {
    "no_GLADE": int,
    "RA": float,
    "Dec": float,
    "B_mag": float,
    "B_mag_err": float,
    "B_mag_abs": float,
    "J_mag": float,
    "J_mag_err": float,
    "H_mag": float,
    "H_mag_err": float,
    "K_mag": float,
    "K_mag_err": float,
    "W1_mag": float,
    "W1_mag_err": float,
    "W2_mag": float,
    "W2_mag_err": float,
    "B_J_mag": float,
    "B_J_mag_err": float,
    "z_helio": float,
    "z_cmb": float,
    "v_err": float,
    "z_err": float,
    "d_L": float,
    "d_L_err": float,
    "mass": float,
    "mass_err": float,
    "merger_rate": float,
    "merger_rate_err": float,
}
converter_columns_catalog = {
    "no_PGC": lambda x: str(x),
    "name_GWGC": lambda x: str(x),
    "name_HyperLEDA": lambda x: str(x),
    "name_2MASS": lambda x: str(x),
    "name_WISExSCOS": lambda x: str(x),
    "name_SDSS_DR16Q": lambda x: str(x),
    "object_type": lambda x: str(x),
    "B_mag_flag": flag_columns_converter,
    "W1_mag_flag": flag_columns_converter,
    "z_flag": flag_columns_converter,
    "dist_flag": flag_columns_converter,
    "mass_flag": flag_columns_converter,
}

# Define column name and format conversion for data format from cosmo hub
cosmo_hub_converter_name_columns_catalog = {
    "glade_no": "no_GLADE",
    "pgc_no": "no_PGC",
    "gwgc_name": "name_GWGC",
    "hyperleda_name": "name_HyperLEDA",
    "2mass_name": "name_2MASS",
    "wisexscos": "name_WISExSCOS",
    "sdss_dr16q": "name_SDSS_DR16Q",
    "object_type_flag": "object_type",
    "ra": "RA",
    "dec": "Dec",
    "b": "B_mag",
    "b_err": "B_mag_err",
    "b_flag": "B_mag_flag",
    "b_abs": "B_mag_abs",
    "j": "J_mag",
    "j_err": "J_mag_err",
    "h": "H_mag",
    "h_err": "H_mag_err",
    "k": "K_mag",
    "k_err": "K_mag_err",
    "w1": "W1_mag",
    "w1_err": "W1_mag_err",
    "w2": "W2_mag",
    "w2_err": "W2_mag_err",
    "w1_flag": "W1_mag_flag",
    "b_j": "B_J_mag",
    "b_j_err": "B_J_mag_err",
    "z_helio": "z_helio",
    "z_cmb": "z_cmb",
    "z_flag": "z_flag",
    "v_err": "v_err",
    "z_err": "z_err",
    "d_l": "d_L",
    "d_l_err": "d_L_err",
    "dist_flag": "dist_flag",
    "m_star": "mass",
    "m_star_err": "mass_err",
    "m_star_flag": "mass_flag",
    "merger_rate": "merger_rate",
    "merger_rate_err": "merger_rate_err",
}
cosmo_hub_converter_columns_catalog = {
    "no_PGC": lambda x: "null" if np.isnan(x) else str(int(x)),
    "name_GWGC": lambda x: str(x),
    "name_HyperLEDA": lambda x: str(x),
    "name_2MASS": lambda x: str(x),
    "name_WISExSCOS": lambda x: str(x),
    "name_SDSS_DR16Q": lambda x: str(x),
    "object_type": lambda x: str(x),
    "B_mag_flag": lambda x: -1 if np.isnan(x) else int(x),
    "W1_mag_flag": lambda x: -1 if np.isnan(x) else int(x),
    "z_flag": lambda x: -1 if np.isnan(x) else int(x),
    "dist_flag": lambda x: -1 if np.isnan(x) else int(x),
    "mass_flag": lambda x: -1 if np.isnan(x) else int(x),
}

# Load argument
args = parser.parse_args()

# Load catalog
logging.info("Start loading catalog file")
tstart = time.time()
_, extension = os.path.splitext(args.input)
if extension == ".parquet":
    catalog = pd.read_parquet(args.input, engine="fastparquet")
    catalog = catalog.rename(cosmo_hub_converter_name_columns_catalog, axis=1)
    for col in cosmo_hub_converter_columns_catalog.keys():
        catalog[col] = catalog[col].apply(cosmo_hub_converter_columns_catalog[col])
else:
    catalog = pd.read_csv(
        args.input,
        sep=" ",
        header=None,
        names=name_columns_catalog,
        dtype=dtype_columns_catalog,
        converters=converter_columns_catalog,
    )
logging.info("Catalog files loaded with success in {0}s".format(time.time() - tstart))

# Compute filter
logging.info("Start computing catalog filter")
tstart = time.time()

logging.info("Determine data validity")
# Remove galaxies with invalid distance or too far
catalog["valid_data"] = True
catalog["valid_data"] &= ~np.isnan(catalog["d_L"])
catalog["valid_data"] &= catalog["d_L"] < args.max_luminosity_distance

# Remove non valid data if requested
if not args.store_all:
    logging.info("Remove non valid data")
    catalog = catalog[catalog["valid_data"]]

logging.info("Catalog filter computed in {0}s".format(time.time() - tstart))

# Create output files
logging.info("Start writing file")
tstart = time.time()

if args.text_format:
    catalog = catalog[catalog["valid_data"]]
    catalog[["no_GLADE", "RA", "Dec", "d_L", "mass"]].to_csv(args.output, index=False)
else:

    # Define table structure hdf5
    class TableCatalog(tables.IsDescription):
        no_GLADE = tables.UInt32Col()
        valid_data = tables.BoolCol()
        RA = tables.Float64Col()
        Dec = tables.Float64Col()
        d_L = tables.Float32Col()
        d_L_err = tables.Float32Col()
        dist_flag = tables.Int8Col()
        mass = tables.Float32Col()
        mass_err = tables.Float32Col()
        mass_flag = tables.Int8Col()
        merger_rate = tables.Float32Col()
        merger_rate_err = tables.Float32Col()
        no_PGC = tables.StringCol(np.max(catalog["no_PGC"].apply(lambda x: len(x))))
        name_GWGC = tables.StringCol(
            np.max(catalog["name_GWGC"].apply(lambda x: len(x)))
        )
        name_HyperLEDA = tables.StringCol(
            np.max(catalog["name_HyperLEDA"].apply(lambda x: len(x)))
        )
        name_2MASS = tables.StringCol(
            np.max(catalog["name_2MASS"].apply(lambda x: len(x)))
        )
        name_WISExSCOS = tables.StringCol(
            np.max(catalog["name_WISExSCOS"].apply(lambda x: len(x)))
        )
        name_SDSS_DR16Q = tables.StringCol(
            np.max(catalog["name_SDSS_DR16Q"].apply(lambda x: len(x)))
        )
        object_type = tables.StringCol(1)
        B_mag = tables.Float32Col()
        B_mag_err = tables.Float32Col()
        B_mag_flag = tables.Int8Col()
        B_mag_abs = tables.Float32Col()
        J_mag = tables.Float32Col()
        J_mag_err = tables.Float32Col()
        H_mag = tables.Float32Col()
        H_mag_err = tables.Float32Col()
        K_mag = tables.Float32Col()
        K_mag_err = tables.Float32Col()
        W1_mag = tables.Float32Col()
        W1_mag_err = tables.Float32Col()
        W2_mag = tables.Float32Col()
        W2_mag_err = tables.Float32Col()
        W1_mag_flag = tables.Int8Col()
        B_J_mag = tables.Float32Col()
        B_J_mag_err = tables.Float32Col()
        z_helio = tables.Float32Col()
        z_cmb = tables.Float32Col()
        z_flag = tables.Int8Col()
        v_err = tables.Float32Col()
        z_err = tables.Float32Col()

    # Create filters for data compression
    if args.compression_algorithm is None:
        filter_catalog = tables.Filters()
    else:
        filter_catalog = tables.Filters(
            complib=args.compression_algorithm,
            complevel=args.compression_level,
            shuffle=args.shuffle,
            bitshuffle=args.bit_shuffle,
        )

    # Create files and node
    h5file = tables.open_file(args.output, mode="w", filters=filter_catalog)
    table = h5file.create_table(
        h5file.root, "catalog", TableCatalog, expectedrows=len(catalog)
    )

    # Determine dtype conversion (pandas to file data model)
    dtype_data_model_file = h5file.root.catalog.dtype
    dtype_conversion = {}
    for k in dtype_data_model_file.names:
        dtype_conversion[k] = dtype_data_model_file[k]

    # Write data to file
    table.append(
        catalog[h5file.root.catalog.colnames].to_records(
            index=False, column_dtypes=dtype_conversion
        )
    )

    # Create index and close file
    h5file.root.catalog.cols.valid_data.create_index()
    h5file.root.catalog.cols.no_GLADE.create_index()
    h5file.close()

logging.info("File written in {0}s".format(time.time() - tstart))
