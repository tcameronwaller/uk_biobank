
"""

This module contains functions for general organization of the raw extraction
data from the UK Biobank accessions.

Combine information from multiple UK Biobank accessions.
Integrate identifier matchings.
Exclude individual cases (persons) who subsequently withdrew consent.
It is also necessary to handle the somewhat awkward format of array fields.

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


# Standard

import sys
#print(sys.path)
import os
import math
import statistics
import pickle
import copy
import random
import itertools

# Relevant

import numpy
import pandas
import scipy.stats

# Custom
import promiscuity.utility as utility
#import promiscuity.plot as plot

###############################################################################
# Functionality

# TODO: TCW 18 October 2021
# TODO: remove most of the old code that's now in the "assembly" module...



# TODO: read in "table_kinship_pairs"

##########
# Initialization


def initialize_directories(
    restore=None,
    path_dock=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        restore (bool): whether to remove previous versions of data
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_dock
    paths["merge_import"] = os.path.join(path_dock, "merge_import")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["merge_import"])
    # Initialize directories.
    utility.create_directories(
        path=paths["merge_import"]
    )
    # Return information.
    return paths


##########
# Read

# TODO: TCW 18 October 2021
# TODO: Currently need to update the name of import table file every accession date


# TODO: TCW 18 October 2021
# TODO: append a prefix to all imported column names...


def read_organize_uk_biobank_import_table(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    The UK Biobank "eid" designates unique persons in the cohort.
    The UK Biobank "IID" matches persons to their genotype information.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Specify directories and files.
    path_table_import = os.path.join(
        path_dock, "access", "ukbiobank_import",
        "waller_import_20211018.derived.csv.gz"
    )
    # Read information from file.
    table_import = pandas.read_csv(
        path_table_import,
        sep=",", # "," or "\t"
        header=0,
        dtype="string",
        na_values=["NA", "<NA>"],
        keep_default_na=True,
        compression="gzip",
    )
    # Organize table.
    table_import = table_import.loc[
        :, table_import.columns.isin([
            "eid",
            "icd_bipolar", "bipolar.diag2", "Smithbipolar", "SmithMood",
            "MHQ.bipolar.Definition", "bipolar", "bipolar.cc",
            "antipsychotics", "antidepressants", "sleepmeds", "opioids",
            "buprenorphine", "methadone", "lamotrigine", "lithium", "valproic_acid",
        ])
    ]

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("read_organize_uk_biobank_import_table()")
        print(table_import)
        utility.print_terminal_partition(level=2)
        pass
    # Return information.
    return table_import


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    The UK Biobank "eid" designates unique persons in the cohort.
    The UK Biobank "IID" matches persons to their genotype information.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Read and organize table of variables across UK Biobank persons for import.
    table_import = read_organize_uk_biobank_import_table(
        path_dock=path_dock,
        report=report,
    )
    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_dock, "assembly", "table_phenotypes.pickle"
    )
    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: read_source()")
        utility.print_terminal_partition(level=3)
        print(table_import)
        utility.print_terminal_partition(level=2)
        pass
    # Compile and return information.
    return {
        "table_import": table_import,
        "table_phenotypes": table_phenotypes,
    }


##########
# Merge


def merge_table_variables_identifiers(
    table_main=None,
    table_import=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        table_main (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_import (object): Pandas data frame of imported phenotype variables
            across UK Biobank cohort for import
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_main)
        utility.print_terminal_partition(level=2)
        print(table_import)
    # Remove rows with null values of merge identifier.
    table_main.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_import.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    # Organize data.

    table_main.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_main["eid"].astype("string")
    table_main.set_index(
        "eid",
        append=False,
        drop=True,
        inplace=True
    )

    table_import.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_import["eid"].astype("string")
    table_import.set_index(
        "eid",
        append=False,
        drop=True,
        inplace=True
    )

    # TODO: TCW 18 October 2021
    # TODO: keep only records from the "main" table...



    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table_identifier_pairs.merge(
        table_ukb_41826,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_pairs", "_41826"),
    )
    table_merge = table_merge.merge(
        table_ukb_43878,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_41826", "_43878"),
    )
    table_merge = table_merge.merge(
        table_import,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_merge", "_import"),
    )
    # Remove excess columns.

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_merge)
    # Return information.
    return table_merge


###############################################################################
# Procedure


# TODO: new parameter for whether to read in and merge the import table...


def execute_procedure(
    path_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Report version.
    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    # Read source information from file.
    # Exclusion identifiers are "eid".
    source = read_source(
        path_dock=path_dock,
        report=True,
    )

    # Merge tables.
    table_merge = merge_table_variables_identifiers(
        table_main=source["table_phenotypes"],
        table_import=source["table_import"],
        report=True,
    )

    pass
