
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
    paths["assembly"] = os.path.join(path_dock, "assembly")
    paths["raw"] = os.path.join(
        path_dock, "assembly", "raw"
    )
    paths["inspection"] = os.path.join(
        path_dock, "assembly", "inspection"
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["assembly"])
    # Initialize directories.
    utility.create_directories(
        path=paths["assembly"]
    )
    utility.create_directories(
        path=paths["raw"]
    )
    utility.create_directories(
        path=paths["inspection"]
    )
    # Return information.
    return paths


##########
# Read

# TODO: Currently need to update the name of import table file every accession date

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
        "waller_import_20210625.derived.csv.gz"
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


def read_ukbiobank_table_column_names(
    path_file=None,
    delimiter=None,
    start=None,
    stop=None,
):
    """
    Reads and extracts column names from UK Biobank data.

    The UK Biobank ukbconv tool's "txt" option exports to a text file with tab
    delimiters between columns; however, there are some quirks about the
    format.

    arguments:
        path_file (str): path to file of data export from UK Biobank ukbconv
            tool
        delimiter (str): delimiter between columns
        start (int): index of line at which to start reading
        stop (int): index of line at which to stop reading

    raises:

    returns:
        (list<str>): list of column names

    """

    # Read lines for samples with genes' signals.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=start,
        stop=stop,
    )
    line = lines[0]
    # Split line's content by delimiter.
    names = line.split(delimiter)
    # Remove excess quotation marks.
    names_simple = list()
    for name in names:
        name = name.replace("\'", "")
        name = name.replace("\"", "")
        names_simple.append(name)
    # Return information.
    return names_simple


def read_collect_ukbiobank_accession_column_names(
    path_table_ukb_41826=None,
    path_table_ukb_43878=None,
):
    """
    Reads and extracts column names from UK Biobank data.

    The UK Biobank ukbconv tool's "txt" option exports to a text file with tab
    delimiters between columns; however, there are some quirks about the
    format.

    arguments:
        path_table_ukb_41826 (str): path to file of data export from UK Biobank
            ukbconv tool
        path_table_ukb_43878 (str): path to file of data export from UK Biobank
            ukbconv tool

    raises:

    returns:
        (list<str>): list of column names

    """

    columns_accession = read_ukbiobank_table_column_names(
        path_file=path_table_ukb_41826,
        delimiter=",", # "," or "\t"
        start=0,
        stop=1,
    )
    columns_new = read_ukbiobank_table_column_names(
        path_file=path_table_ukb_43878,
        delimiter=",", # "," or "\t"
        start=0,
        stop=1,
    )
    columns_accession.extend(columns_new)
    columns_accession_unique = utility.collect_unique_elements(
        elements=columns_accession,
    )
    # Return information.
    return columns_accession_unique


def extract_organize_variables_types(
    table_ukbiobank_variables=None,
    columns=None,
    extra_pairs=None,
):
    """
    Organize information about data types of UK Biobank phenotype variables.

    arguments:
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        columns (list<str>): unique names of columns in UK Biobank tables
        extra_pairs (dict<str>): extra key value pairs to include

    raises:

    returns:
        (dict<str>): pairs of variables and data types

    """

    # Copy data.
    table_variables = table_ukbiobank_variables.copy(deep=True)
    # Organize information.
    table_variables = table_variables.loc[
        :, table_variables.columns.isin(["field", "type", "coverage"])
    ]
    records = utility.convert_dataframe_to_records(data=table_variables)
    # Iterate across records for rows.
    # Collect variables' names and types.
    variables_types = dict()
    variables_types.update(extra_pairs)
    for record in records:
        field = str(record["field"]).strip()
        type = str(record["type"]).strip()
        # Search UK Biobank table columns for matches.
        for column in columns:
            column_field = column.split("-")[0].strip()
            if (field == column_field):
                variables_types[column] = type
                pass
            pass
        pass
    # Return information.
    return variables_types


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
    path_table_ukbiobank_variables = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "table_ukbiobank_phenotype_variables.tsv"
    )
    path_exclusion_identifiers = os.path.join(
        path_dock, "access", "ukbiobank_phenotypes",
        "list_exclusion_identifiers.txt"
    )
    path_table_identifier_pairs = os.path.join(
        path_dock, "access", "ukbiobank_phenotypes",
        "table_identifier_pairs.csv"
    )
    path_table_kinship_pairs = os.path.join(
        path_dock, "access", "ukbiobank_phenotypes",
        "table_kinship_pairs.dat"
    )
    path_table_ukb_41826 = os.path.join(
        path_dock, "access", "ukbiobank_phenotypes",
        "ukb41826.raw.csv"
    )
    path_table_ukb_43878 = os.path.join(
        path_dock, "access", "ukbiobank_phenotypes",
        "ukb43878.raw.csv"
    )
    # Read all column names from UK Biobank tables.
    columns_accession = read_collect_ukbiobank_accession_column_names(
        path_table_ukb_41826=path_table_ukb_41826,
        path_table_ukb_43878=path_table_ukb_43878,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("unique column names: " + str(len(columns_accession)))
        print(columns_accession)
        utility.print_terminal_partition(level=2)
        pass
    # Determine variable types.
    table_ukbiobank_variables = pandas.read_csv(
        path_table_ukbiobank_variables,
        sep="\t",
        header=0,
    )
    table_ukbiobank_variables["field"].astype("string")
    variables_types = extract_organize_variables_types(
        table_ukbiobank_variables=table_ukbiobank_variables,
        columns=columns_accession,
        extra_pairs={
            "IID": "string",
            "eid": "string",
        },
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("31-0.0: " + str(variables_types["31-0.0"]))
        print("21022-0.0: " + str(variables_types["21022-0.0"]))
    # Read information from file.
    exclusion_identifiers = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_exclusion_identifiers
    )
    table_identifier_pairs = pandas.read_csv(
        path_table_identifier_pairs,
        sep=",",
        header=0,
        dtype="string",
    )
    table_kinship_pairs = pandas.read_csv(
        path_table_kinship_pairs,
        sep="\s+",
        header=0,
        dtype={
            "ID1": "string",
            "ID2": "string",
            "HetHet": "float32",
            "IBS0": "float32",
            "Kinship": "float32",
        },
    )
    table_ukb_41826 = pandas.read_csv(
        path_table_ukb_41826,
        sep=",", # "," or "\t"
        header=0,
        dtype=variables_types,
        na_values=["<NA>"],
        keep_default_na=True,
    )
    table_ukb_43878 = pandas.read_csv(
        path_table_ukb_43878,
        sep=",", # "," or "\t"
        header=0,
        dtype=variables_types,
        na_values=["<NA>"],
        keep_default_na=True,
    )
    # Compile and return information.
    return {
        "table_import": table_import,
        "table_ukbiobank_variables": table_ukbiobank_variables,
        "columns_accession": columns_accession,
        "exclusion_identifiers": exclusion_identifiers,
        "table_identifier_pairs": table_identifier_pairs,
        "table_kinship_pairs": table_kinship_pairs,
        "table_ukb_41826": table_ukb_41826,
        "table_ukb_43878": table_ukb_43878,
    }


##########
# Remove irrelevant field columns


def parse_field_instance_columns_to_keep(
    instances_raw=None,
):
    """
    Parse the field's instances to keep.

    arguments:
        instances_raw (str): raw string of field's instances to keep

    raises:

    returns:
        (list<str>): field's instances to keep

    """

    instances_simple = instances_raw.replace("\'", "")
    instances_simple = instances_simple.replace("\"", "")
    instances = list()
    for instance_raw in instances_simple.split(";"):
        instance = str(instance_raw).strip()
        if len(instance) > 0:
            instances.append(str(instance))
            pass
        pass
    return instances


def determine_keep_column_field_instance(
    column=None,
    table_ukbiobank_variables=None,
    report=None,
):
    """
    Determines whether to keep a column for UK Biobank field instance.

    arguments:
        column (str): column name from UK Biobank accessions
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        report (bool): whether to print reports

    raises:

    returns:
        (bool): whether to keep the column

    """

    # Determine whether column matches format of an original field instance.
    if ("-" in column):
        # Default value.
        keep = False
        # Organize information.
        column_field = column.split("-")[0].strip()
        column_instance = column.split("-")[1].strip()
        # int(column_field)
        instances_array = table_ukbiobank_variables.at[
            int(column_field), "instances_array",
        ]
        instances_keep_raw = table_ukbiobank_variables.at[
            int(column_field), "instances_keep",
        ]
        # Report.
        if (
            report and
            (len(str(column_field)) > 0) and
            (str(column_field) in ["31", "50", "22009"])
        ):
            utility.print_terminal_partition(level=4)
            print(column_field)
        # Determine whether to keep column for field's instance.
        if not pandas.isna(instances_array):
            if str(instances_array).strip().lower() == "yes":
                keep = True
        if not pandas.isna(instances_keep_raw):
            # Organize field instances to keep.
            instances_keep = parse_field_instance_columns_to_keep(
                instances_raw=instances_keep_raw,
            )
            # Report.
            if (
                report and
                (len(str(column_field)) > 0) and
                (str(column_field) in ["31", "50", "22009"])
            ):
                print(instances_keep_raw)
                print(instances_keep)
            if str(column_instance) in instances_keep:
                keep = True
                pass
            pass
        pass
    else:
        # Column is not an original field instance.
        # Keep the column.
        keep = True
    # Return information.
    return keep


def determine_ukbiobank_field_instance_columns_keep(
    columns_accession=None,
    table_ukbiobank_variables=None,
    extra_names=None,
    report=None,
):
    """
    Organizes column names for variable fields and instances.

    arguments:
        columns_accession (list<str>): column names from UK Biobank accessions
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        extra_names (list<str>): extra names to include
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): column names for variable fields and instances

    """

    # Copy data.
    table_ukbiobank_variables = table_ukbiobank_variables.copy(deep=True)
    # Organize information.
    table_ukbiobank_variables = table_ukbiobank_variables.loc[
        :, table_ukbiobank_variables.columns.isin([
            "field", "instances_array", "instances_keep"
        ])
    ]
    table_ukbiobank_variables["field"].astype("string")
    table_ukbiobank_variables.set_index(
        "field",
        drop=True,
        inplace=True,
    )
    # Iterate on actual column names from the accession table.
    # Determine whether to keep column.
    columns_keep = list()
    columns_keep.extend(extra_names)
    for column in columns_accession:
        column = str(column)
        # Determine whether to keep column.
        keep = determine_keep_column_field_instance(
            column=column,
            table_ukbiobank_variables=table_ukbiobank_variables,
            report=False,
        )
        if keep:
            columns_keep.append(column)
            pass
        pass
    # Return information.
    return columns_keep


def remove_table_irrelevant_field_instance_columns(
    table_ukbiobank_variables=None,
    columns_accession=None,
    table_ukb_41826=None,
    table_ukb_43878=None,
    report=None,
):
    """
    Removes irrelevant columns and rows from data.

    arguments:
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        columns_accession (list<str>): column names from UK Biobank accessions
        table_ukb_41826 (object): Pandas data frame of variables from UK
            Biobank phenotype data accession 41826
        table_ukb_43878 (object): Pandas data frame of variables from UK
            Biobank phenotype data accession 43878
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data frames after removal of
            irrelevant columns and rows

    """

    # Copy data.
    table_variables = table_ukbiobank_variables.copy(deep=True)
    table_ukb_41826 = table_ukb_41826.copy(deep=True)
    table_ukb_43878 = table_ukb_43878.copy(deep=True)
    # Determine which columns to keep for UK Biobank fields and instances.
    columns_keep = determine_ukbiobank_field_instance_columns_keep(
        columns_accession=columns_accession,
        table_ukbiobank_variables=table_variables,
        extra_names=["IID", "eid"],
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("columns to keep: " + str(len(columns_keep)))
        #print(columns_keep[0:50])
        print(columns_keep)
        utility.print_terminal_partition(level=3)
        print("...before pruning...")
        print("table_ukb_41826 shape: " + str(table_ukb_41826.shape))
        print("table_ukb_43878 shape: " + str(table_ukb_43878.shape))
        utility.print_terminal_partition(level=4)
    # Remove all irrelevant columns.
    table_ukb_41826 = table_ukb_41826.loc[
        :, table_ukb_41826.columns.isin(columns_keep)
    ]
    table_ukb_43878 = table_ukb_43878.loc[
        :, table_ukb_43878.columns.isin(columns_keep)
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after pruning...")
        print("table_ukb_41826 shape: " + str(table_ukb_41826.shape))
        print("table_ukb_43878 shape: " + str(table_ukb_43878.shape))
        utility.print_terminal_partition(level=4)
    # Collect information.
    bucket = dict()
    bucket["table_ukb_41826"] = table_ukb_41826
    bucket["table_ukb_43878"] = table_ukb_43878
    # Return information.
    return bucket


##########
# Organization


def simplify_field_instances_array_row(
    row=None,
    field=None,
    delimiter=None,
    report=None,
):
    """
    Simplify field instances for array values.

    arguments:
        row (object): Pandas data frame row
        field (str): identifier of UK Biobank field
        delimiter (str): delimiter for string representation of array values
        report (bool): whether to print reports

    raises:

    returns:
        (str): text array with semicolon delimiters

    """

    # Copy Pandas series for row.
    row = row.copy(deep=True)
    # Convert Pandas series row to dictionary.
    record = row.to_dict()
    # Collect all non-missing values from the field's columns for instances.
    values = list()
    fields = list()
    for key in record.keys():
        # Select original field-instance columns.
        if ("-" in key):
            key_field = key.split("-")[0].strip()
            if (str(field) == str(key_field)):
                fields.append(key)
                value = record[key]
                if (not pandas.isna(value)):
                    values.append(value)
                    pass
                pass
            pass
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("matching fields to " + str(field) + " :")
        print(fields)
    # Collect unique values.
    values_unique = utility.collect_unique_elements(
        elements=values,
    )
    # Combine values with text delimiter.
    if len(values_unique) > 0:
        text_array = delimiter.join(values_unique)
    else:
        text_array = ""
    # Return information.
    return text_array


def simplify_field_instances_array_columns(
    table_ukbiobank_variables=None,
    table_ukb_raw=None,
    delimiter=None,
    report=None,
):
    """
    Simplify field instances for array values.

    arguments:
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        table_ukb_raw (object): Pandas data frame of information from UK
            Biobank
        delimiter (str): delimiter for string representation of array values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information from UK Biobank

    """

    # Copy data.
    table_ukb = table_ukb_raw.copy(deep=True)
    table_variables = table_ukbiobank_variables.copy(deep=True)
    # Organize information.
    table_variables["field"].astype("string")
    table_variables = table_variables.loc[
        :, table_variables.columns.isin(["field", "type", "instances_array"])
    ]
    table_variables = table_variables.loc[
        ~pandas.isna(table_variables["instances_array"]), :
    ]
    fields_array_instances = table_variables["field"].to_list()
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("fields to simplify: " + str(len(fields_array_instances)))
        print(fields_array_instances)
    # Iterate on UK Biobank fields with array instances.
    for field in fields_array_instances:
        # Create new column with text array of all non-missing values from the
        # field's original columns.
        column_new = str(str(field) + "_array")
        table_ukb[column_new] = table_ukb.apply(
            lambda row:
                simplify_field_instances_array_row(
                    row=row,
                    field=field,
                    delimiter=delimiter,
                    report=False,
                ),
            axis="columns", # apply across rows
        )
        # Collect original columns for the field's array instances.
        columns_drop = list()
        for column in table_ukb.columns.to_list():
            # Select original field-instance columns.
            if ("-" in column):
                column_field = column.split("-")[0].strip()
                if (str(field) == str(column_field)):
                    # Column is an original instance of the field.
                    # Only original instance columns have the "-" delimiter.
                    columns_drop.append(column)
                    pass
                pass
            pass
        # Report.
        if report:
            utility.print_terminal_partition(level=3)
            print("dropping columns for field: " + str(field))
            print(columns_drop)
        # Drop columns.
        table_ukb.drop(
            labels=columns_drop,
            axis="columns",
            inplace=True
        )
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after simplification of array fields...")
        print("table_ukb shape: " + str(table_ukb.shape))
        utility.print_terminal_partition(level=4)
    # Return information.
    return table_ukb


def merge_table_variables_identifiers(
    table_identifier_pairs=None,
    table_ukb_41826=None,
    table_ukb_43878=None,
    table_import=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        table_identifier_pairs (object): Pandas data frame of associations
            between "IID" and "eid"
        table_ukb_41826 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 41826
        table_ukb_43878 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 43878
        table_import (object): Pandas data frame of phenotype variables across
            UK Biobank cohort for import
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_identifier_pairs)
        utility.print_terminal_partition(level=2)
        print(table_ukb_41826)
        utility.print_terminal_partition(level=2)
        print(table_ukb_43878)
        utility.print_terminal_partition(level=2)
        print(table_import)
    # Remove rows with null values of merge identifier.
    table_identifier_pairs.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_ukb_41826.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_ukb_43878.dropna(
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
    table_identifier_pairs.astype("string")
    table_identifier_pairs.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    table_ukb_41826["eid"].astype("string")
    table_ukb_41826.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    table_ukb_43878["eid"].astype("string")
    table_ukb_43878.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    table_import["eid"].astype("string")
    table_import.set_index(
        "eid",
        drop=True,
        inplace=True,
    )

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


def exclude_persons_ukbiobank_consent(
    exclusion_identifiers=None,
    table=None,
    report=None,
):
    """
    Removes entries with specific identifiers from data.

    arguments:
        exclusion_identifiers (list<str>): identifiers of persons who withdrew
            consent from UK Biobank
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Initial table dimensions: " + str(table.shape))
        utility.print_terminal_partition(level=4)
        print("Exclusion of persons: " + str(len(exclusion_identifiers)))
    # Copy data.
    table = table.copy(deep=True)
    # Filter data entries.
    table_exclusion = table.loc[
        ~table.index.isin(exclusion_identifiers), :
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Final data dimensions: " + str(table_exclusion.shape))
    # Return information.
    return table_exclusion


def drop_null_records_all_variables(
    table=None,
    columns_any=None,
    report=None,
):
    """
    Drops rows with null or missing values across all columns.

    arguments:
        table (object): Pandas data frame of information about UK Biobank
            phenotype variables
        columns_any (str): column names for which to drop rows with any missing
            values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about UK Biobank phenotype
            variables

    """

    # Copy data.
    table = table.copy(deep=True)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...before dropping null records...")
        print("table shape: " + str(table.shape))
        print(table)
        utility.print_terminal_partition(level=4)
    # Remove rows with null values across all columns.
    table.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    table.dropna(
        axis="index",
        how="any",
        subset=columns_any,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after dropping null records...")
        print("table shape: " + str(table.shape))
        print(table)
        utility.print_terminal_partition(level=4)
    # Return information.
    return table


##########
# Write


def write_product_raw(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory
            raises:

    returns:

    """

    # Specify directories and files.
    path_table_ukb_41826_text = os.path.join(
        path_parent, "table_ukb_41826.tsv"
    )
    # Write information to file.
    information["table_ukb_41826"].to_csv(
        path_or_buf=path_table_ukb_41826_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_inspection(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory
            raises:

    returns:

    """

    # Specify directories and files.
    path_table_text = os.path.join(
        path_parent, "table_phenotypes.tsv"
    )
    # Write information to file.
    information["table_phenotypes"].to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_assembly(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory
            raises:

    returns:

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_parent, "table_phenotypes.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        path_parent, "table_phenotypes.tsv"
    )
    path_table_kinship_pairs = os.path.join(
        path_parent, "table_kinship_pairs.pickle"
    )
    path_table_kinship_pairs_text = os.path.join(
        path_parent, "table_kinship_pairs.tsv"
    )
    # Write information to file.
    information["table_phenotypes"].to_pickle(
        path_table_phenotypes
    )
    information["table_phenotypes"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["table_kinship_pairs"].to_pickle(
        path_table_kinship_pairs
    )
    information["table_kinship_pairs"].to_csv(
        path_or_buf=path_table_kinship_pairs_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product(
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    write_product_raw(
        information=information["raw"],
        path_parent=paths["raw"],
    )
    write_product_inspection(
        information=information["inspection"],
        path_parent=paths["inspection"],
    )
    write_product_assembly(
        information=information["assembly"],
        path_parent=paths["assembly"],
    )
    pass

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

    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")

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

    # Remove data columns for irrelevant variable instances.
    prune = remove_table_irrelevant_field_instance_columns(
        table_ukbiobank_variables=source["table_ukbiobank_variables"],
        columns_accession=source["columns_accession"],
        table_ukb_41826=source["table_ukb_41826"],
        table_ukb_43878=source["table_ukb_43878"],
        report=True,
    )
    # Simplify UK Biobank fields with multiple instances.
    # Reduce these multiple field instance columns to arrays.
    table_ukb_41826_simple = simplify_field_instances_array_columns(
        table_ukbiobank_variables=source["table_ukbiobank_variables"],
        table_ukb_raw=prune["table_ukb_41826"],
        delimiter=";",
        report=True,
    )

    # Merge tables.
    table_merge = merge_table_variables_identifiers(
        table_identifier_pairs=source["table_identifier_pairs"],
        table_ukb_41826=table_ukb_41826_simple,
        table_ukb_43878=prune["table_ukb_43878"],
        table_import=source["table_import"],
        report=True,
    )

    # Exclude persons who withdrew consent from the UK Biobank.
    table_exclusion = exclude_persons_ukbiobank_consent(
        exclusion_identifiers=source["exclusion_identifiers"],
        table=table_merge,
        report=True,
    )
    # Drop any records (persons) with null values across all variables.
    table_valid = drop_null_records_all_variables(
        table=table_exclusion,
        columns_any=["31-0.0"], # "IID"
        report=True,
    )

    utility.print_terminal_partition(level=2)
    print("table_kinship_pairs")
    print(source["table_kinship_pairs"])

    # Write out raw tables for inspection.
    # Collect information.
    information = dict()
    information["raw"] = dict()
    information["inspection"] = dict()
    information["assembly"] = dict()
    information["raw"]["table_ukb_41826"] = (
        source["table_ukb_41826"].iloc[0:10000, :]
    )
    information["inspection"]["table_phenotypes"] = (
        table_valid.iloc[0:10000, :]
    )
    information["assembly"]["table_phenotypes"] = table_valid
    information["assembly"]["table_kinship_pairs"] = (
        source["table_kinship_pairs"]
    )
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )
    pass
