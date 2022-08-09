"""
Organize information and define variables from data of the U.K. Biobank.

Importation paths within this package "uk_biobank" require this package to be a
subpackage within a higher-level package that manages execution.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of UK_Biobank
    (https://github.com/tcameronwaller/uk_biobank/).

    UK_Biobank supports analyses on data from the U.K. Biobank.
    Copyright (C) 2022 Thomas Cameron Waller

    UK_Biobank is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    UK_Biobank is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with UK_Biobank. If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import sys
#print(sys.path)
import os
#import shutil
#import csv
import math
import statistics
import pickle
import copy
import random
import itertools
import textwrap
import time

# Relevant

import numpy
import scipy.stats
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"
import networkx

# Custom
import promiscuity.utility as utility
import promiscuity.decomposition as decomp
import uk_biobank.stratification as ukb_strat # problem when executing uk_biobank not as sub-directory...

###############################################################################
# Functionality


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
    paths["organization"] = os.path.join(path_dock, "organization")
    paths["cohorts_models"] = os.path.join(
        path_dock, "organization", "cohorts_models"
    )

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["organization"])
    # Initialize directories.
    utility.create_directories(
        path=paths["organization"]
    )
    #utility.create_directories(
    #    path=paths["cohorts_models"]
    #)
    # Return information.
    return paths


##########
# Read


def read_source(
    source=None,
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        source (str): name of directory for source table, either "importation"
            or "assembly"
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    if (str(source) == "importation"):
        path_table_phenotypes = os.path.join(
            path_dock, "importation", "table_phenotypes.pickle"
        )
    elif (str(source) == "assembly"):
        path_table_phenotypes = os.path.join(
            path_dock, "assembly", "table_phenotypes.pickle"
        )

    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    if False:
        path_table_ukb_samples = os.path.join(
            path_dock, "access", "ukb46237_imp_chr21_v3_s487320.sample"
        )
        table_ukb_samples = pandas.read_csv(
            path_table_ukb_samples,
            sep="\s+",
            header=0,
            dtype="string",
        )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        #"table_ukb_samples": table_ukb_samples,
    }


def read_source_fields_codes_interpretations(
    path_dock=None,
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

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_field_54 = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "table_ukbiobank_field_54_code_interpretation.tsv"
    )
    path_table_field_55 = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "table_ukbiobank_field_55_code_interpretation.tsv"
    )
    path_table_field_22000 = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "table_ukbiobank_field_22000_code_interpretation.tsv"
    )

    path_table_drug_class = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "table_ukbiobank_drug_class_translation.tsv"
    )

    # Organize code interpretations.

    # Data-field 54.
    table_field_54 = pandas.read_csv(
        path_table_field_54,
        sep="\t",
        header=0,
        dtype={
            "code": "string",
            "interpretation": "string",
            "coding": "string",
            "meaning": "string",
            "region_name": "string",
            "region": numpy.int32,
        },
    )
    table_field_54.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_field_54.set_index(
        "code",
        append=False,
        drop=True,
        inplace=True
    )
    field_54_codes_interpretations = table_field_54.to_dict(
        orient="index",
    )

    # Data-field 55.
    table_field_55 = pandas.read_csv(
        path_table_field_55,
        sep="\t",
        header=0,
        dtype={
            "code": "string",
            "interpretation": "string",
            "coding": "string",
            "meaning": "string",
            "season": numpy.int32,
            "season_strict": numpy.int32,
            "day_length": numpy.int32,
            "note": "string",
        },
    )
    table_field_55.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_field_55.set_index(
        "code",
        append=False,
        drop=True,
        inplace=True
    )
    field_55_codes_interpretations = table_field_55.to_dict(
        orient="index",
    )

    # Data-field 22000.
    table_field_22000 = pandas.read_csv(
        path_table_field_22000,
        sep="\t",
        header=0,
        dtype={
            "code": "string",
            "array": "string",
            "interpretation": "string",
            "coding": "string",
            "meaning": "string",
            "node_id": "string",
            "parent_id": "string",
        },
    )
    table_field_22000.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_field_22000.set_index(
        "code",
        append=False,
        drop=True,
        inplace=True
    )
    field_22000_codes_interpretations = table_field_22000.to_dict(
        orient="index",
    )

    # Drug classes.
    table_drug_class = pandas.read_csv(
        path_table_drug_class,
        sep="\t",
        header=0,
        dtype={
            "original": "string",
            "novel": "string",
        },
    )
    table_drug_class.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_drug_class.set_index(
        "original",
        append=False,
        drop=True,
        inplace=True
    )
    drug_class_translations = table_drug_class.to_dict(
        orient="index",
    )

    # Compile and return information.
    return {
        "field_54": field_54_codes_interpretations,
        "field_55": field_55_codes_interpretations,
        "field_22000": field_22000_codes_interpretations,
        "drug_class": drug_class_translations,
    }


def read_source_medication_codes_classes(
    path_dock=None,
):
    """
    Reads and organizes source information from file.

    This function reads and organizes information from tables that are
    derivatives of the sources below.
    1. "Supplementary Data 1" from Yeda Wu et al, Nature
    Communications, 2019 (PubMed: 31015401).
    2. "ATC Match List (Additional_file_2)" from Philip Duncan Appleby et al,
    Research Square, 2019 (https://www.researchsquare.com/article/rs-9729/v1).

    These tables match medication codes in UK Biobank data-field "20003" to
    medication classes in the Anatomical Therapeutic Chemical (ATC)
    classification system
    (https://www.who.int/tools/atc-ddd-toolkit/atc-classification).

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_wu = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "ukbiobank_field_20003_medications",
        "31015401_wu_2019_supplementary_data_1_edit.tsv"
    )
    path_table_appleby = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "ukbiobank_field_20003_medications",
        "appleby_2019_supplementary_data_2_edit.tsv"
    )
    # Read information from file.
    table_wu = pandas.read_csv(
        path_table_wu,
        sep="\t",
        header=0,
        dtype={
            "code_ukbiobank": "string",
            "medication_name": "string",
            "medication_category": "string",
            "code_atc": "string",
            "group_atc": "string",
        },
    )
    table_appleby = pandas.read_csv(
        path_table_appleby,
        sep="\t",
        header=0,
        dtype={
            "code_ukbiobank": "string",
            "medication_name": "string",
            "code_atc": "string",
            "group_atc": "string",
        },
    )
    # Compile and return information.
    return {
        "table_wu": table_wu,
        "table_appleby": table_appleby,
    }


##########
# Genotype


def match_column_field(
    column=None,
    field=None,
):
    """
    Determine whether a column represents an original instance of a field.

    arguments:
        column (str): name of column
        field (str): identifier of field for which to collect columns

    raises:

    returns:
        (bool): whether column name matches field

    """

    # Determine whether column matches format of an original field instance.
    if ("-" in str(column)):
        # Column is an original instance of the field.
        # Only original instance columns have the "-" delimiter.
        column_field = str(column).split("-")[0].strip()
        if (str(column_field) == str(field)):
            match = True
        else:
            match = False
    else:
        match = False
    # Return information.
    return match


def collect_sort_table_field_instance_columns(
    field=None,
    table=None,
    report=None,
):
    """
    Collects names of columns that represent instances of a specific field.

    arguments:
        field (str): identifier of field for which to collect columns
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): names of columns

    """

    # Extract table columns.
    columns = copy.deepcopy(table.columns.to_list())
    # Collect columns that represent instances of the field.
    columns_field = list(filter(
        lambda column: match_column_field(column=column, field=field),
        columns
    ))
    # Sort columns.
    # Be careful to sort by the integer value of the instance.
    #columns_sort = sorted(columns_field, reverse=False)
    columns_instance = list()
    for column in columns_field:
        column_instance = str(column).split("-")[1].strip()
        columns_instance.append(column_instance)
    columns_integer = list()
    for column in columns_instance:
        column_integer = str(column).split(".")[1].strip()
        columns_integer.append(int(column_integer))
    table_columns = pandas.DataFrame(
        data={
            "column_name": columns_field,
            "column_instance": columns_instance,
            "column_integer": columns_integer,
        }
    )
    # Sort table rows.
    table_columns.sort_values(
        by=["column_integer"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    columns_sort = table_columns["column_name"].to_list()
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Sorted column names")
        print(table_columns)
    # Return information.
    return columns_sort


def translate_column_field_instance_names(
    columns=None,
    prefix=None,
    base=None,
):
    """
    Translates names of columns.

    arguments:
        columns (list<str>): names of columns for field instances
        prefix (str): prefix for translations of column names
        base (int): initial integer for index

    raises:

    returns:
        (list<str>): names of columns

    """

    index = base
    translations = dict()
    for column in columns:
        translations[str(column)] = str(prefix + str(index))
        index += 1
    return translations


def convert_table_prefix_columns_variables_types_float(
    prefix=None,
    table=None,
):
    """
    Converts data variable types.

    The UK Biobank encodes several nominal variables with integers. Missing
    values for these variables necessitates some attention to type conversion
    from string to float for analysis.

    arguments:
        prefix (str): prefix for translations of column names
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    table = table.copy(deep=True)
    # Determine which columns to convert.
    columns = table.columns.to_list()
    columns_match = list(filter(
        lambda column: (str(prefix) in str(column)),
        columns
    ))
    # Convert data variable types.
    table_type = utility.convert_table_columns_variables_types_float(
        columns=columns_match,
        table=table,
    )
    # Return information.
    return table_type


def interpret_ancestry_white_british(
    field_22006=None,
):
    """
    Intepret UK Biobank's coding for field 22006.

    Data-Field 22006 relates to self report of racial and ethnic identity in
    Data-Field 21000.

    Data-Field "22006": "Genetic ethnic grouping"
    UK Biobank data coding "1002" for variable field "22006".
    1: "Caucasian"

    Accommodate inexact float values.

    Note:
    "
    Genetic ethnic group.
    Indicates samples who self-identified as 'White British' according to Field
    21000 and have very similar genetic ancestry based on a principal components
    analysis of the genotypes.
    "

    arguments:
        field_22006 (float): UK Biobank field 22006, whether person is in
            categorical ancestral group, "white british"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_22006)) and
        (0.5 <= field_22006 and field_22006 < 1.5)
    ):
        # 1: "Caucasian"
        value = 1
    else:
        # null
        value = float("nan")
    # Return.
    return value


def determine_ancestry_white_british(
    field_22006=None,
):
    """
    Determine whether persons experienced bilateral oophorectomy.

    arguments:
        field_22006 (float): UK Biobank field 22006, whether person is in
            categorical ancestral group, "white british"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret oophorectomy.
    white_british = interpret_ancestry_white_british(
        field_22006=field_22006,
    )
    # Comparison.
    if (
        (not pandas.isna(white_british)) and
        (white_british == 1)
    ):
        value = 1
    else:
        value = 0
    # Return information.
    return value


def organize_genotype_principal_component_variables(
    table=None,
    report=None,
):
    """
    Organizes information about principal components on genotypes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort


    """

    # Copy data.
    table = table.copy(deep=True)
    # Extract columns for genotype principal components.
    columns = collect_sort_table_field_instance_columns(
        field="22009",
        table=table,
        report=report,
    )
    # Translate column names.
    translations = translate_column_field_instance_names(
        columns=columns,
        prefix="genotype_pc_",
        base=1,
    )
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Convert variable types.
    table = convert_table_prefix_columns_variables_types_float(
        prefix="genotype_pc_",
        table=table,
    )
    columns_type = [
        "22006-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Organize indicator variable for genetic ancestry in the category of
    # "white british".
    table["white_british"] = table.apply(
        lambda row:
            determine_ancestry_white_british(
                field_22006=row["22006-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("translations of genotype PC column names...")
        for old in translations.keys():
            print("   " + old + ": " + translations[old])
        utility.print_terminal_partition(level=3)
        # Column names and values.
        table_report = table.loc[
            :, table.columns.str.startswith("genotype_pc_")
        ]
        utility.print_terminal_partition(level=2)
        print("Translation of columns for genotype principal components: ")
        print(table_report)
        utility.print_terminal_partition(level=3)
        # Variable types.
        utility.print_terminal_partition(level=2)
        print("After type conversion")
        print(table_report.dtypes)
        utility.print_terminal_partition(level=3)
        # Categorical ancestry and ethnicity.
        utility.print_terminal_partition(level=2)
        print("Categorical ancestry and ethnicity")
        table_white_british = table.loc[
            (table["white_british"] == 1), :
        ]
        print("total persons: " + str(table.shape[0]))
        print("white british persons: " + str(table_white_british.shape[0]))
    # Return information.
    return table


##########
# General assessment



# review: TCW, 25 January 2022
# review: TCW, 27 January 2022; data-field "54" coding and site names
def interpret_assessment_site(
    field_54=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field 54.

    Data-Field "54": "UK Biobank assessment centre"
    UK Biobank data-coding "10" for data-field "54".
    11012: "Barts"
    11021: "Birmingham"
    11011: "Bristol"
    11008: "Bury"
    11003: "Cardiff"
    11024: "Cheadle (revisit)"
    11020: "Croydon"
    11005: "Edinburgh"
    11004: "Glasgow"
    11018: "Hounslow"
    11010: "Leeds"
    11016: "Liverpool"
    11001: "Manchester"
    11017: "Middlesborough"
    11009: "Newcastle"
    11013: "Nottingham"
    11002: "Oxford"
    11007: "Reading"
    11014: "Sheffield"
    10003: "Stockport (pilot)"
    11006: "Stoke"
    11022: "Swansea"
    11023: "Wrexham"
    11025: "Cheadle (imaging)"
    11026: "Reading (imaging)"
    11027: "Newcastle (imaging)"
    11028: "Bristol (imaging)"

    Note:
    "
    The UK Biobank assessment centre at which participant consented.
    "

    Accommodate inexact float values.

    arguments:
        field_54 (float): UK Biobank field 54, assessment center
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_54)) and
        (10003 <= field_54 and field_54 < 11029)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_54 = copy.deepcopy(field_54)
        field_54_string = str(int(field_54))
        if (field_54_string in codes_interpretations.keys()):
            interpretation = (
                codes_interpretations[field_54_string]["interpretation"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


# review: TCW, 25 January 2022
# review: TCW, 27 January 2022; data-field "54" coding and site regions
def interpret_assessment_site_region(
    field_54=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field 54.

    Data-Field "54": "UK Biobank assessment centre"
    UK Biobank data-coding "10" for data-field "54".

    Note:
    "
    The UK Biobank assessment centre at which participant consented.
    "

    Accommodate inexact float values.

    arguments:
        field_54 (float): UK Biobank field 54, assessment center
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_54)) and
        (10003 <= field_54 and field_54 < 11029)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_54 = copy.deepcopy(field_54)
        field_54_string = str(int(field_54))
        if (field_54_string in codes_interpretations.keys()):
            interpretation = int(
                codes_interpretations[field_54_string]["region"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


def interpret_assessment_month_order(
    field_55=None,
):
    """
    Intepret UK Biobank's coding for data field 55.
    Data-Field "55": "Month of attending assessment centre"
    UK Biobank data coding "8" for variable field "55".
    1: "January"
    2: "February"
    3: "March"
    4: "April"
    5: "May"
    6: "June"
    7: "July"
    8: "August"
    9: "September"
    10: "October"
    11: "November"
    12: "December"
    Note:
    "
    Calendar month that participant attended a UK Biobank assessment centre.
    Automatically acquired at Reception stage.
    "
    Accommodate inexact float values.
    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
    raises:
    returns:
        (float): interpretation value
    """

    # Interpret field code.
    if (
        (not pandas.isna(field_55)) and
        (0.5 <= field_55 and field_55 < 12.5)
    ):
        # The variable has a valid value for month.
        value = float(field_55)
    else:
        # Missing value.
        value = float("nan")
    # Return.
    return value


# review: TCW, 16 March 2022
def interpret_assessment_month_category(
    field_55=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data field 55. Translate names of months to
    lower-case.

    Data-Field "55": "Month of attending assessment centre"
    UK Biobank data coding "8" for variable field "55".
    1: "January"
    2: "February"
    3: "March"
    4: "April"
    5: "May"
    6: "June"
    7: "July"
    8: "August"
    9: "September"
    10: "October"
    11: "November"
    12: "December"

    Note:
    "
    Calendar month that participant attended a UK Biobank assessment centre.
    Automatically acquired at Reception stage.
    "

    Accommodate inexact float values.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_55)) and
        (0.5 <= field_55 and field_55 < 12.5)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_55 = copy.deepcopy(field_55)
        field_55_string = str(int(field_55))
        if (field_55_string in codes_interpretations.keys()):
            interpretation = str(
                codes_interpretations[field_55_string]["interpretation"]
            ).strip().lower()
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


# review: TCW, 25 January 2022
def interpret_assessment_month_season(
    field_55=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field "55".

    0: Winter
    1: Spring or Autumn
    2: Summer

    Data-Field "55": "Month of attending assessment centre"
    UK Biobank data-coding "8" for data-field "55".
    1: "January"
    2: "February"
    3: "March"
    4: "April"
    5: "May"
    6: "June"
    7: "July"
    8: "August"
    9: "September"
    10: "October"
    11: "November"
    12: "December"

    Note:
    "
    Calendar month that participant attended a UK Biobank assessment centre.
    Automatically acquired at Reception stage.
    "

    Accommodate inexact float values.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_55)) and
        (0.5 <= field_55 and field_55 < 12.5)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_55 = copy.deepcopy(field_55)
        field_55_string = str(int(field_55))
        if (field_55_string in codes_interpretations.keys()):
            interpretation = int(
                codes_interpretations[field_55_string]["season_strict"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


# review: TCW, 30 January 2022
def interpret_assessment_month_day_length(
    field_55=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field "55".

    Data-Field "55": "Month of attending assessment centre"
    UK Biobank data-coding "8" for data-field "55".
    1: "January"
    2: "February"
    3: "March"
    4: "April"
    5: "May"
    6: "June"
    7: "July"
    8: "August"
    9: "September"
    10: "October"
    11: "November"
    12: "December"

    Note:
    "
    Calendar month that participant attended a UK Biobank assessment centre.
    Automatically acquired at Reception stage.
    "

    Accommodate inexact float values.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_55)) and
        (0.5 <= field_55 and field_55 < 12.5)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_55 = copy.deepcopy(field_55)
        field_55_string = str(int(field_55))
        if (field_55_string in codes_interpretations.keys()):
            interpretation = int(
                codes_interpretations[field_55_string]["day_length"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


# review: TCW, 16 March 2022
def interpret_assessment_month_season_summer_autumn(
    field_55=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field "55".

    0: Winter or Spring
    1: Summer or Autumn

    Data-Field "55": "Month of attending assessment centre"
    UK Biobank data-coding "8" for data-field "55".
    1: "January"
    2: "February"
    3: "March"
    4: "April"
    5: "May"
    6: "June"
    7: "July"
    8: "August"
    9: "September"
    10: "October"
    11: "November"
    12: "December"

    Note:
    "
    Calendar month that participant attended a UK Biobank assessment centre.
    Automatically acquired at Reception stage.
    "

    Accommodate inexact float values.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_55)) and
        (0.5 <= field_55 and field_55 < 12.5)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_55 = copy.deepcopy(field_55)
        field_55_string = str(int(field_55))
        if (field_55_string in codes_interpretations.keys()):
            interpretation = int(
                codes_interpretations[field_55_string]["season_summer_autumn"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


# review: TCW, 09 February 2022
def interpret_assessment_month_season_stratification(
    field_55=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field "55".

    Data-Field "55": "Month of attending assessment centre"
    UK Biobank data-coding "8" for data-field "55".
    1: "January"
    2: "February"
    3: "March"
    4: "April"
    5: "May"
    6: "June"
    7: "July"
    8: "August"
    9: "September"
    10: "October"
    11: "November"
    12: "December"

    Note:
    "
    Calendar month that participant attended a UK Biobank assessment centre.
    Automatically acquired at Reception stage.
    "

    Accommodate inexact float values.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_55)) and
        (0.5 <= field_55 and field_55 < 12.5)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_55 = copy.deepcopy(field_55)
        field_55_string = str(int(field_55))
        if (field_55_string in codes_interpretations.keys()):
            interpretation = int(
                codes_interpretations[field_55_string]["season_stratification"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing value.
        interpretation = str("nan")
    # Return.
    return interpretation


# review: TCW, 25 January 2022
def interpret_sex_self_report(
    field_31=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 31.

    Data-Field "31": "Sex"
    UK Biobank data-coding "9" for data-field "31".
    0: "Female"
    1: "Male"

    Accommodate inexact float values.

    Note:
    "Acquired from central registry at recruitment, but in some cases updated by
    the participant. Hence this field may contain a mixture of the sex the NHS
    had recorded for the participant and self-reported sex."

    arguments:
        field_31 (float): UK Biobank field 31, self-report sex

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_31)) and
        (-0.5 <= field_31 and field_31 < 1.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (-0.5 <= field_31 and field_31 < 0.5):
            # 0: "Female"
            interpretation = 0
        elif (0.5 <= field_31 and field_31 < 1.5):
            # 1: "Male"
            interpretation = 1
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


# review: TCW, 25 January 2022
def interpret_sex_genetic(
    field_22001=None,
):
    """
    Intepret UK Biobank's data-coding for data-field "22001".

    Data-Field "22001": "Genetic sex"
    UK Biobank data-coding "9" for data-field "22001".
    0: "Female"
    1: "Male"

    Accommodate inexact float values.

    Note:
    "Note that in over 300 cases the genetic sex differs from the self-reported
    value in Field 31."

    arguments:
        field_22001 (float): UK Biobank field 22001, genetic sex

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_22001)) and
        (-0.5 <= field_22001 and field_22001 < 1.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (-0.5 <= field_22001 and field_22001 < 0.5):
            # 0: "Female"
            interpretation = 0
        elif (0.5 <= field_22001 and field_22001 < 1.5):
            # 1: "Male"
            interpretation = 1
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


# review: TCW, 25 January 2022
def interpret_sex_chromosome_aneuploidy(
    field_22019=None,
):
    """
    Intepret UK Biobank's data-coding for data-field "22019".

    Data-Field "22019": "Sex chromosome aneuploidy"
    UK Biobank data-coding "1" for data-field "22019".
    1: "Yes"

    Accommodate inexact float values.

    Note:
    "Sex chromosome aneuploidy marker. This indicates samples which were
    identified as putatively carrying sex chromosome configurations that are not
    either XX or XY. These were identified by looking at average log2Ratios for
    Y and X chromosomes."

    arguments:
        field_22019 (float): UK Biobank data-field 22019, sex-chromosome
            aneuploidy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_22019)) and
        (0.5 <= field_22019 and field_22019 < 1.5)
    ):
        # 1: "Yes"
        interpretation = 1
    else:
        # Absence of affirmative indicator suggests opposite is true.
        interpretation = 0
    # Return.
    return interpretation


def interpret_age(
    field_raw=None,
):
    """
    Intepret UK Biobank's data-coding for a data-field for an age.

    Accommodate inexact float values.

    Note:

    arguments:
        field_raw (float): UK Biobank data-field for an age

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_raw))
    ):
        # The variable has a valid value.
        # Interpret the value.
        interpretation = float(field_raw)
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def interpret_body_mass_index(
    field_21001=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 21001.

    Accommodate inexact float values.

    arguments:
        field_21001 (float): UK Biobank field 21001, body mass index

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_21001))
    ):
        # The variable has a valid value.
        # Interpret the value.
        interpretation = float(field_21001)
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


# review: TCW, 08 June 2022
def interpret_self_report_ancestry_ethnicity(
    field_21000=None,
):
    """
    Intepret UK Biobank's data-coding for data-field "21000".

    Data-Field "21000": "Ethnic background"
    UK Biobank data-coding "1001" for data-field "21000".
    ...
    6: "Other ethnic group"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_21000 (float): UK Biobank data-field 21000, person's self-report
            of ancestry and or ethnicity

    raises:

    returns:
        (str): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_21000)) and
        (-3.5 <= field_21000 and field_21000 < 5000)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (0.5 <= field_21000 and field_21000 < 5000):
            # Value is one of many interpretable codes.
            value = str(int(copy.deepcopy(field_21000)))
        elif (-1.5 <= field_21000 and field_21000 < -0.5):
            # -1: "Do not know"
            value = "NA"
        elif (-3.5 <= field_21000 and field_21000 < -2.5):
            # -3: "Prefer not to answer"
            value = "NA"
        else:
            # Uninterpretable value
            value = "NA"
    else:
        # Missing or uninterpretable value
        value = "NA"
    # Return.
    return value


def interpret_genotype_batch(
    field_22000=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field 22000.

    Data-Field "22000": "Genotype measurement batch"
    UK Biobank data-coding "22000" for data-field "22000".

    11 genotype batches correspond to "BiLEVE" genotype array.
    95 genotype batches correspond to "Axiom" genotype array.

    Accommodate inexact float values.

    arguments:
        field_22000 (float): UK Biobank field 22000, genotype batch
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_22000)) and
        (-12 <= field_22000 and field_22000 < 2001)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_22000 = copy.deepcopy(field_22000)
        field_22000_string = str(int(field_22000))
        if (field_22000_string in codes_interpretations.keys()):
            interpretation = (
                codes_interpretations[field_22000_string]["interpretation"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing or uninterpretable value.
        interpretation = str("nan")
    # Return.
    return interpretation


def interpret_genotype_array(
    field_22000=None,
    codes_interpretations=None,
):
    """
    Intepret UK Biobank's coding for data-field 22000.

    Data-Field "22000": "Genotype measurement batch"
    UK Biobank data-coding "22000" for data-field "22000".

    11 genotype batches correspond to "BiLEVE" genotype array.
    95 genotype batches correspond to "Axiom" genotype array.

    Accommodate inexact float values.

    arguments:
        field_22000 (float): UK Biobank field 22000, genotype batch
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_22000)) and
        (-12 <= field_22000 and field_22000 < 2001)
    ):
        # The variable has a valid value.
        # Determine code interpretation.
        field_22000 = copy.deepcopy(field_22000)
        field_22000_string = str(int(field_22000))
        if (field_22000_string in codes_interpretations.keys()):
            interpretation = (
                codes_interpretations[field_22000_string]["array"]
            )
        else:
            # Uninterpretable value.
            interpretation = str("nan")
    else:
        # Missing or uninterpretable value.
        interpretation = str("nan")
    # Return.
    return interpretation


def determine_assessment_site_category(
    field_54=None,
    codes_interpretations=None,
):
    """
    Determine assessment site.

    arguments:
        field_54 (float): UK Biobank field 54, assessment center
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_site(
        field_54=field_54,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_assessment_site_region(
    field_54=None,
    codes_interpretations=None,
):
    """
    Determine assessment site.

    arguments:
        field_54 (float): UK Biobank field 54, assessment center
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_site_region(
        field_54=field_54,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_assessment_month_order(
    field_55=None,
):
    """
    Determine assessment month.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_month_order(
        field_55=field_55,
    )
    # Return information.
    return value


def determine_assessment_month_category(
    field_55=None,
    codes_interpretations=None,
):
    """
    Determine assessment month.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_month_category(
        field_55=field_55,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_assessment_month_season(
    field_55=None,
    codes_interpretations=None,
):
    """
    Determine assessment month.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_month_season(
        field_55=field_55,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_assessment_month_season_summer_autumn(
    field_55=None,
    codes_interpretations=None,
):
    """
    Determine assessment month.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_month_season_summer_autumn(
        field_55=field_55,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_assessment_month_season_stratification(
    field_55=None,
    codes_interpretations=None,
):
    """
    Determine assessment month.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_month_season_stratification(
        field_55=field_55,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_assessment_month_day_length(
    field_55=None,
    codes_interpretations=None,
):
    """
    Determine the ordinal approximate day length from the assessment month.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_assessment_month_day_length(
        field_55=field_55,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_genotype_batch(
    field_22000=None,
    codes_interpretations=None,
):
    """
    Determine genotype batch.

    arguments:
        field_22000 (float): UK Biobank field 22000, genotype batch
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_genotype_batch(
        field_22000=field_22000,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_genotype_array(
    field_22000=None,
    codes_interpretations=None,
):
    """
    Determine genotype array.

    arguments:
        field_22000 (float): UK Biobank field 22000, genotype batch
        codes_interpretations (dict<dict<string>>): interpretations for each
            code of UK Biobank field

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_genotype_array(
        field_22000=field_22000,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_genotype_array_axiom_logical_binary(
    genotype_array=None,
):
    """
    Determine whether genotype array was "Axiom" or "BiLEVE".

    1: "Axiom"
    0: "BiLEVE"

    arguments:
        genotype_array (str): name of genotype array, either "axiom" or "bileve"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(genotype_array)) and
        ((genotype_array == "axiom") or (genotype_array == "bileve"))
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (genotype_array == "axiom"):
            interpretation = 1
        elif (genotype_array == "bileve"):
            interpretation = 0
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return information.
    return interpretation


def determine_genotype_availability(
    genotype_array=None,
):
    """
    Determine whether genotype is available for the phenotype record.

    Code:
    0: genotype is unavailable
    1: genotype is available

    arguments:
        genotype_array (str): name of genotype array, either "axiom" or "bileve"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(genotype_array)) and
        ((genotype_array == "axiom") or (genotype_array == "bileve"))
    ):
        # The variable has a valid value.
        # Interpret the value.
        # 1: genotype is available
        interpretation = 1
    else:
        # Missing or uninterpretable value
        # Interpret the value.
        # 0: genotype is unavailable
        interpretation = 0
    # Return information.
    return interpretation


# review: TCW, 25 January 2022
def determine_consensus_biological_sex_y(
    field_31=None,
    field_22001=None,
    field_22019=None,
):
    """
    Determine consensus biological sex (not social gender).

    Prioritize interpretation of the genetic sex variable, and use self report
    to fill in any missing values. Exclude records for persons whose genotypes
    suggested aneuploidy in the sex chromosomes.

    Use a logical binary encoding for presence of Y chromosome.
    0 : female (false, XX)
    1 : male (true, XY)

    Code is same as UK Biobank for data-fields "31" and "22001".
    0: "Female"
    1: "Male"

    arguments:
        field_31 (float): UK Biobank data-field 31, person's self-reported sex
        field_22001 (float): UK Biobank data-field 22001, genetic sex
        field_22019 (float): UK Biobank data-field 22019, sex-chromosome
            aneuploidy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret self report of biological sex.
    sex_self_report = interpret_sex_self_report(
        field_31=field_31,
    )
    # Interpret genetic biological sex.
    sex_genetic = interpret_sex_genetic(
        field_22001=field_22001,
    )
    aneuploidy_sex = interpret_sex_chromosome_aneuploidy(
        field_22019=field_22019,
    )
    # Comparison.
    # Determine whether there was aneuploidy in the sex chromosomes.
    if (
        (not pandas.isna(aneuploidy_sex)) and
        (aneuploidy_sex == 0)
    ):
        # Prioritize genetic sex.
        if (not pandas.isna(sex_genetic)):
            # Genetic sex variable has a valid value.
            # Prioritize this variable.
            value = sex_genetic
        elif (not pandas.isna(sex_self_report)):
            # Person has missing value for genetic sex.
            # Self-report sex variable has a valid value.
            # Resort to self-report sex in absence of genetic sex.
            value = sex_self_report
        else:
            # Missing information.
            value = float("nan")
    else:
        # Aneuploidy in the sex chromosomes complicates definition of biological
        # sex.
        value = float("nan")
    # Return information.
    return value


def determine_biological_sex_x(
    sex_y=None,
):
    """
    Determine biological sex, using an integer encoding for the count of X
    chromosomes.
    2 : female (XX)
    1 : male (XY)

    arguments:
        sex_y (float): sex as logical presence of Y chromosome

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(sex_y)) and
        ((sex_y == 0) or (sex_y == 1))
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (sex_y == 0):
            # sex_y: 0, "female"
            # sex_x: 2, "female"
            interpretation = 2
        elif (sex_y == 1):
            # sex_y: 1, "male"
            # sex_x: 1, "male"
            interpretation = 1
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return information.
    return interpretation


def determine_biological_sex_text(
    sex_y=None,
):
    """
    Translate binary representation of sex to textual representation.

    arguments:
        sex_y (float): sex as logical presence of Y chromosome

    raises:

    returns:
        (str): text representation of person's sex

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(sex_y)) and
        ((sex_y == 0) or (sex_y == 1))
    ):
        # The variable has a valid value.
        if (sex_y == 0):
            # sex_y: 0, "female"
            sex_text = "female"
        elif (sex_y == 1):
            # sex_y: 1, "male"
            sex_text = "male"
        else:
            # Ambiguous, uninterpretable, or missing information.
            sex_text = ""
    else:
        # Ambiguous, uninterpretable, or missing information.
        sex_text = ""
    # Return information.
    return sex_text


def determine_age(
    field_raw=None,
):
    """
    Determine age with consideration of possible range.

    arguments:
        field_raw (float): UK Biobank data-field for an age

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret age.
    age = interpret_age(
        field_raw=field_raw,
    )
    # Comparison.
    if (
        (not pandas.isna(age)) and
        (0.0 <= age and age < 150.0)
    ):
        # Age variable has a valid value.
        value = age
    else:
        # Lack of information or unreasonable information.
        value = float("nan")
    # Return information.
    return value


def determine_body_mass_index(
    field_21001=None,
):
    """
    Determine body mass index (BMI) with consideration of possible range.

    Maximum realistic body mass index (BMI) is less than 185, which would be the
    BMI of a person with height 6 feet 1 inches and body weight 1,400 pounds.

    https://www.nhlbi.nih.gov/health/educational/lose_wt/BMI/bmicalc.htm

    arguments:
        field_21001 (float): UK Biobank field 21001, body mass index

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret age.
    body = interpret_body_mass_index(
        field_21001=field_21001,
    )
    # Comparison.
    if (
        (not pandas.isna(body)) and
        (5.0 <= body and body < 190.0)
    ):
        # Age variable has a valid value.
        value = body
    else:
        # Lack of information or unreasonable information.
        value = float("nan")
    # Return information.
    return value


# review: TCW; 08 June 2022
def determine_self_report_ancestry_ethnicity(
    field_21000=None,
):
    """
    Determine whether person self reported ancestry and or ethnicity as
    "white" or other.

    Code.
    0: "other"
    1: "white"

    arguments:
        field_21000 (float): UK Biobank data-field 21000, person's self-report
            of ancestry and or ethnicity

    raises:

    returns:
        (float): interpretation value

    """

    # Define codes.
    codes_white = ["1", "1001", "1002", "1003"]
    codes_other = [
        "2", "2001", "2002", "2003", "2004",
        "3", "3001", "3002", "3003", "3004",
        "4", "4001", "4002", "4003", "5", "6"
    ]

    # Interpret self report of biological sex.
    code = interpret_self_report_ancestry_ethnicity(
        field_21000=field_21000,
    )
    # Interpret code.
    if (
        (not pandas.isna(code)) and
        (code in codes_white)
    ):
        # Person self reported ancestry and or ethnicity as "white".
        value = 1
    elif (
        (not pandas.isna(code)) and
        (code in codes_other)
    ):
        # Person self reported ancestry and or ethnicity in one of several
        # categories "other".
        value = 0
    else:
        # Ambiguous, uninterpretable, or missing information.
        value = float("nan")
    # Return information.
    return value


# review: TCW, 26 January 2022
def create_categorical_variable_indicators(
    table=None,
    index=None,
    column=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Creates binary indicator (dummy) variables for values of a single
    categorical variable.

    Pandas function drops the original column for the categorical variable.
    Copy, split, and merge tables to preserve the original column with its
    indicator variables.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        index (str): name of table's index column by which to merge after create
            of columns for indicator variables
        column (str): name of column with categorical, character, string
            variable for which to define binary dummies and reduce by principal
            components
        prefix (str): prefix for names of new dummy and principal component
            columns in table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information about creation of categorical indicator variables

    """

    # Copy information.
    table = table.copy(deep=True)
    table_indicators = table.copy(deep=True)
    # Organize values of the categorical variable.
    column_values = str(prefix + "_values")
    table_indicators[column_values] = table_indicators.apply(
        lambda row:
            str(row[column]).replace(" ", "_"),
        axis="columns", # apply function to each row
    )
    # Count instances of all unique values of the categorical variable.
    counts_values = table_indicators[column_values].value_counts(
        normalize=False,
        sort=True,
        ascending=False, # whether to sort in order of ascending counts
        dropna=False,
    )
    # Organize tables.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )
    table_indicators.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_indicators = table_indicators.loc[
        :, table_indicators.columns.isin([index, column_values])
    ]
    # Drop any rows with missing values in the original column.
    table_indicators.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
    )
    # Create binary dummies for variable categories.
    table_indicators = pandas.get_dummies(
        table_indicators,
        prefix=prefix,
        prefix_sep=separator,
        dummy_na=False, # whether to create indicators for missing values
        columns=[column_values],
        drop_first=False, # whether to create "k - 1" dummies, adequate
        dtype=numpy.uint8,
    )
    table_indicators.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )
    # Collect names of new columns for indicator variables.
    columns_indicators = copy.deepcopy(table_indicators.columns.to_list())
    # Merge indicators to main table.
    table = table.merge(
        table_indicators,
        how="outer",
        left_on=index,
        right_on=index,
        suffixes=("_original", "_indicator"),
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: create_categorical_variable_indicators()")
        utility.print_terminal_partition(level=3)
        # Description.
        unique_values = table[column].unique()
        count_unique_values = unique_values.size
        count_indicators = len(columns_indicators)
        utility.print_terminal_partition(level=4)
        print("counts of instances of all original unique values...")
        print(counts_values)
        utility.print_terminal_partition(level=4)
        print("count of original unique values: " + str(count_unique_values))
        print(unique_values)
        print("count of dummy indicator variables: " + str(count_indicators))
        # Organize table.
        columns_report = copy.deepcopy(columns_indicators)
        columns_report.insert(0, column)
        columns_report.insert(0, "IID")
        #columns_report.insert(0, "eid")
        table_report = table.copy(deep=True)
        table_report = table_report.loc[
            :, table_report.columns.isin(columns_report)
        ]
        table_report = table_report[[*columns_report]]
        utility.print_terminal_partition(level=4)
        print("Here is table with dummy indicator variables...")
        print(table_report)
        # Aggregate by sum.
        table_report.drop(
            ["IID", column],
            axis=1,
            inplace=True,
        )
        series_aggregate = table_report.aggregate(
            lambda column_current: numpy.nansum(column_current.to_numpy()),
            axis="index", # Apply function to each column in table.
        )
        table_aggregate = series_aggregate.to_frame(name="sum")
        utility.print_terminal_partition(level=4)
        print("Here are sums for each category...")
        #print(series_aggregate)
        print(table_aggregate)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["columns_indicators"] = columns_indicators
    # Return information.
    return pail


# review: TCW, 26 January 2022
def reduce_categorical_variable_indicators(
    table=None,
    index=None,
    column=None,
    columns_indicators=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Reduces binary indicator (dummy) variables by Principal Component Analysis.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        index (str): name of table's index column by which to merge after create
            of columns for indicator variables
        column (str): name of column with categorical variable for which to
            define binary dummies and reduce by principal components
        columns_indicators (list<str>): names of columns with binary indicator
            (dummy) variables to reduce by principal components
        prefix (str): prefix for names of new principal component columns in
            table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information about creation of categorical indicator variables

    """

    # Copy information.
    table = table.copy(deep=True)
    table_indicators = table.copy(deep=True)
    # Organize tables.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )
    table_indicators.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    columns_keep = copy.deepcopy(columns_indicators)
    columns_keep.insert(0, index)
    table_indicators = table_indicators.loc[
        :, table_indicators.columns.isin(columns_keep)
    ]
    table_indicators.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )

    # Reduce dimensionality of indicator variables.
    pail_reduction = (
        decomp.organize_principal_components_by_singular_value_decomposition(
            table=table_indicators,
            index_name=index,
            prefix=prefix,
            separator=separator,
            report=False, # report is extensive
        )
    )
    # Organize table.
    table_component_scores = (
        pail_reduction["table_component_scores"].copy(deep=True)
    )
    table_component_scores.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_component_scores.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )
    columns_components = copy.deepcopy(table_component_scores.columns.to_list())
    table = table.merge(
        table_component_scores,
        how="outer",
        left_on=index,
        right_on=index,
        suffixes=("_original", "_component"),
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: reduce_categorical_variable_indicators()")
        utility.print_terminal_partition(level=3)
        # Description.
        unique_values = table[column].unique()
        count_unique_values = unique_values.size
        count_indicators = len(columns_indicators)
        count_components = len(columns_components)
        print("count of original unique values: " + str(count_unique_values))
        print(unique_values)
        print("count of dummy indicator variables: " + str(count_indicators))
        print("count of component variables: " + str(count_components))
        # Compare methods for Pincipal Component Analysis.
        decomp.compare_principal_components_methods(
            table=table_indicators,
            index_name=index,
            prefix=prefix,
            separator=separator,
            report=True,
        )
        # Organize table.
        columns_report = copy.deepcopy(columns_components)
        columns_report.insert(0, column)
        columns_report.insert(0, "IID")
        #columns_report.insert(0, "eid")
        table_report = table.copy(deep=True)
        table_report = table_report.loc[
            :, table_report.columns.isin(columns_report)
        ]
        table_report = table_report[[*columns_report]]
        utility.print_terminal_partition(level=4)
        print("Here is table with component variables...")
        print(table_report)
        utility.print_terminal_partition(level=4)
        print("Here are the new columns for principal components...")
        print(columns_components)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["columns_components"] = columns_components
    # Return information.
    return pail


def create_reduce_categorical_variable_indicators(
    table=None,
    index=None,
    column=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Creates binary indicator (dummy) variables for values of a single
    categorical variable and reduces these by principal components.

    Pandas function drops the original column for the categorical variable.
    Copy, split, and merge tables to preserve the original column with its
    indicator variables.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        index (str): name of table's index column by which to merge after create
            of columns for indicator variables
        column (str): name of column with categorical, character, string
            variable for which to define binary dummies and reduce by principal
            components
        prefix (str): prefix for names of new dummy indicator and principal
            component columns in table
        separator (str): separator for names of new columns, preferrably
            underscore "_" and not hyphen "-"
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Create dummy indicator variables for each category of original variable.
    pail_indicator = create_categorical_variable_indicators(
        table=table,
        index=index,
        column=column,
        prefix=str(prefix + "_" + "indicator"),
        separator=separator,
        report=report,
    )
    # Reduce dimensionality of indicator variables.
    pail_reduction = reduce_categorical_variable_indicators(
        table=pail_indicator["table"],
        index=index,
        column=column,
        columns_indicators=pail_indicator["columns_indicators"],
        prefix=str(prefix + "_" + "component"),
        separator=separator,
        report=report,
    )
    # Return information.
    return pail_reduction["table"]


def report_ordinal_stratifications_by_sex_continuous_variables(
    sex_text=None,
    variables=None,
    table=None,
):
    """
    Stratify persons to ordinal bins by their values of continuous variables.

    arguments:
        sex_text (str): name of column for textual representation of sex
        variables (list<str>): names of columns for continuous variables by
            which to define stratification variables
        table (object): Pandas data frame of features (columns) across
            observations (rows)

    raises:

    returns:

    """

    #print(pandas.qcut(
    #    table_female[variable], q=3, labels=[0, 1, 2,], retbins=True,
    #))
    #print(pandas.qcut(
    #    table_female[variable], q=3, labels=[0, 1, 2,],
    #).value_counts())

    # Copy information.
    table = table.copy(deep=True)
    # Report.
    utility.print_terminal_partition(level=3)
    print("report_ordinal_stratifications_by_sex_continuous_variable")

    # Iterate on variables and cohorts.
    for variable in variables:
        utility.print_terminal_partition(level=5)
        print("variable: " + str(variable))
        utility.print_terminal_partition(level=5)

        # Collect records of information about each cohort.
        records = list()

        record = dict()
        record["name"] = "female low"
        record["table"] = table.loc[
            (
                (table[sex_text] == "female") &
                (table[str(variable + "_grade_female")] == 0)
            ), :
        ]
        records.append(record)

        record = dict()
        record["name"] = "female middle"
        record["table"] = table.loc[
            (
                (table[sex_text] == "female") &
                (table[str(variable + "_grade_female")] == 1)
            ), :
        ]
        records.append(record)

        record = dict()
        record["name"] = "female high"
        record["table"] = table.loc[
            (
                (table[sex_text] == "female") &
                (table[str(variable + "_grade_female")] == 2)
            ), :
        ]
        records.append(record)

        record = dict()
        record["name"] = "male low"
        record["table"] = table.loc[
            (
                (table[sex_text] == "male") &
                (table[str(variable + "_grade_male")] == 0)
            ), :
        ]
        records.append(record)

        record = dict()
        record["name"] = "male middle"
        record["table"] = table.loc[
            (
                (table[sex_text] == "male") &
                (table[str(variable + "_grade_male")] == 1)
            ), :
        ]
        records.append(record)

        record = dict()
        record["name"] = "male high"
        record["table"] = table.loc[
            (
                (table[sex_text] == "male") &
                (table[str(variable + "_grade_male")] == 2)
            ), :
        ]
        records.append(record)

        # Iterate on tables.
        for record in records:
            array = copy.deepcopy(record["table"][variable].dropna().to_numpy())
            count_string = str(int(array.size))
            mean_string = str(round(numpy.nanmean(array), 3))
            minimum_string = str(round(numpy.nanmin(array), 3))
            maximum_string = str(round(numpy.nanmax(array), 3))
            print("..........")
            print("cohort: " + record["name"])
            print("count: " + count_string)
            print("mean: " + mean_string)
            print("minimum: " + minimum_string)
            print("maximum: " + maximum_string)
            pass
        pass
    pass


def define_ordinal_stratifications_by_sex_continuous_variables(
    index=None,
    sex_text=None,
    variables=None,
    table=None,
    report=None,
):
    """
    Stratify persons to ordinal bins by their values of continuous variables.

    arguments:
        index (str): name of column for identifier index
        sex_text (str): name of column for textual representation of sex
        variables (list<str>): names of columns for continuous variables by
            which to define stratification variables
        table (object): Pandas data frame of features (columns) across
            observations (rows)
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about persons

    """

    # Copy information.
    table = table.copy(deep=True)
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Females.
    table_female = table.copy(deep=True)
    table_female = table_female.loc[
        (table_female[sex_text] == "female"), :
    ]
    # Males.
    table_male = table.copy(deep=True)
    table_male = table_male.loc[
        (table_male[sex_text] == "male"), :
    ]
    # Collect records.
    table_collection = pandas.DataFrame()
    # Stratify separately by sex.
    for variable in variables:
        variable_female = str(variable + "_grade_female")
        variable_male = str(variable + "_grade_male")
        table_female[variable_female] = pandas.qcut(
            table_female[variable],
            q=3,
            labels=[0, 1, 2,],
        )
        table_male[variable_male] = pandas.qcut(
            table_male[variable],
            q=3,
            labels=[0, 1, 2,],
        )
    # Combine records for female and male persons.
    table_collection = table_collection.append(
        table_female,
        ignore_index=True,
    )
    table_collection = table_collection.append(
        table_male,
        ignore_index=True,
    )
    # Organize table.
    table_collection.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_collection.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )
    # Report.
    if report:
        report_ordinal_stratifications_by_sex_continuous_variables(
            sex_text=sex_text,
            variables=variables,
            table=table_collection,
        )
    # Return information.
    return table_collection


def report_genotype_arrays_batches(
    table=None,
):
    """
    Report counts of genotypes in each genotype array and genotype batch.

    arguments:
        table (object): Pandas data frame of features (columns) across
            observations (rows)

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Report.
    utility.print_terminal_partition(level=3)
    print("report_genotype_arrays_batches")
    # Summarize genotype arrays.
    utility.print_terminal_partition(level=4)
    print("genotype array: Axiom")
    utility.print_terminal_partition(level=4)
    table_axiom = table.loc[
        (table["genotype_array_axiom"] == 1), :
    ]
    count_axiom = copy.deepcopy(table_axiom.shape[0])
    print("records in genotype array: " + str(count_axiom))
    utility.print_terminal_partition(level=4)
    print("genotype array: BiLEVE")
    utility.print_terminal_partition(level=4)
    table_bileve = table.loc[
        (table["genotype_array_axiom"] == 0), :
    ]
    count_bileve = copy.deepcopy(table_bileve.shape[0])
    print("records in genotype array: " + str(count_bileve))

    # Extract unique batch identifiers.
    if False:
        batches = list(set(table["genotype_batch"].to_list()))
        batches = list(filter(
            lambda batch: (batch != "nan"),
            batches
        ))
        batches_sort = sorted(batches, reverse=False)
        count_batches = len(batches_sort)
        utility.print_terminal_partition(level=4)
        print("count of genotype batches : " + str(count_batches))
        utility.print_terminal_partition(level=4)
        # Iterate on batches.
        for batch in batches_sort:
            if (not pandas.isna(batch)):
                utility.print_terminal_partition(level=5)
                print("genotype batch: " + str(batch))
                utility.print_terminal_partition(level=5)
                table_batch = table.loc[
                    (table["genotype_batch"] == batch), :
                ]
                count_batch = copy.deepcopy(table_batch.shape[0])
                print("records in genotype batch: " + str(count_batch))
                pass
            pass
        pass
    pass


def report_assessment_season_region(
    column_name=None,
    table=None,
):
    """
    Reports counts of persons who had assessment in each ordinal season.

    "assessment_season"
    0: winter
    1: spring, autumn
    2: summer

    "assessment_region"
    0: Southern England
    1: Central England
    2: Northern England

    arguments:
        column_name (str): name of column for assessment season or region
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # assessment_season

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_female_0 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[column_name] == 0)
        ), :
    ]
    table_female_1 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[column_name] == 1)
        ), :
    ]
    table_female_2 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[column_name] == 2)
        ), :
    ]

    # Males.
    table_male = table.loc[
        (
            (table["sex_text"] == "male")
        ), :
    ]
    table_male_0 = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[column_name] == 0)
        ), :
    ]
    table_male_1 = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[column_name] == 1)
        ), :
    ]
    table_male_2 = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[column_name] == 2)
        ), :
    ]

    # Counts.
    count_female = table_female.shape[0]
    count_male = table_male.shape[0]
    count_female_0 = table_female_0.shape[0]
    count_female_1 = table_female_1.shape[0]
    count_female_2 = table_female_2.shape[0]
    count_male_0 = table_male_0.shape[0]
    count_male_1 = table_male_1.shape[0]
    count_male_2 = table_male_2.shape[0]
    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_assessment_season_region()")
    utility.print_terminal_partition(level=3)
    print("assessment categoriy variable: " + str(column_name))
    utility.print_terminal_partition(level=4)
    print("total female persons: " + str(count_female))
    print("female persons with value '0': " + str(count_female_0))
    print("female persons with value '1': " + str(count_female_1))
    print("female persons with value '2': " + str(count_female_2))
    utility.print_terminal_partition(level=4)
    print("total male persons: " + str(count_male))
    print("male persons with value '0': " + str(count_male_0))
    print("male persons with value '1': " + str(count_male_1))
    print("male persons with value '2': " + str(count_male_2))
    pass


def organize_assessment_basis_variables(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about general attributes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Read reference table of interpretations of variable codes.
    source = read_source_fields_codes_interpretations(
        path_dock=path_dock,
    )

    # Copy information.
    table = table.copy(deep=True)

    # TODO: TCW, 30 January 2022
    # - site
    # - latitude
    # - season
    # - day_length

    # Determine assessment site.
    table["site"] = table.apply(
        lambda row:
            determine_assessment_site_category(
                field_54=row["54-0.0"],
                codes_interpretations=source["field_54"],
            ),
        axis="columns", # apply function to each row
    )
    table["region"] = table.apply(
        lambda row:
            determine_assessment_site_region(
                field_54=row["54-0.0"],
                codes_interpretations=source["field_54"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine month of assessment.
    table["month_order"] = table.apply(
        lambda row:
            determine_assessment_month_order(
                field_55=row["55-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    table["month_name"] = table.apply(
        lambda row:
            determine_assessment_month_category(
                field_55=row["55-0.0"],
                codes_interpretations=source["field_55"],
            ),
        axis="columns", # apply function to each row
    )
    table["season"] = table.apply(
        lambda row:
            determine_assessment_month_season(
                field_55=row["55-0.0"],
                codes_interpretations=source["field_55"],
            ),
        axis="columns", # apply function to each row
    )
    table["day_length"] = table.apply(
        lambda row:
            determine_assessment_month_day_length(
                field_55=row["55-0.0"],
                codes_interpretations=source["field_55"],
            ),
        axis="columns", # apply function to each row
    )
    table["season_summer_autumn"] = table.apply(
        lambda row:
            determine_assessment_month_season_summer_autumn(
                field_55=row["55-0.0"],
                codes_interpretations=source["field_55"],
            ),
        axis="columns", # apply function to each row
    )
    table["season_stratification"] = table.apply(
        lambda row:
            determine_assessment_month_season_stratification(
                field_55=row["55-0.0"],
                codes_interpretations=source["field_55"],
            ),
        axis="columns", # apply function to each row
    )

    # Create binary indicators of categories and reduce their dimensionality.
    table = create_reduce_categorical_variable_indicators(
        table=table,
        index="eid",
        column="site",
        prefix="site",
        separator="_",
        report=True,
    )
    table = create_reduce_categorical_variable_indicators(
        table=table,
        index="eid",
        column="month_name",
        prefix="month",
        separator="_",
        report=True,
    )
    # Determine sex consensus between self-report and genotypic sex.
    # Reserve the variable names "sex" or "SEX" for recognition in PLINK2.
    # Use logical binary representation of presence of Y chromosome.
    table["sex_y"] = table.apply(
        lambda row:
            determine_consensus_biological_sex_y(
                field_31=row["31-0.0"],
                field_22001=row["22001-0.0"],
                field_22019=row["22019-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Use integer representation of count of X chromosomes.
    table["sex_x"] = table.apply(
        lambda row:
            determine_biological_sex_x(
                sex_y=row["sex_y"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine text representation of person's sex.
    table["sex_text"] = table.apply(
        lambda row:
            determine_biological_sex_text(
                sex_y=row["sex_y"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine age.
    table["age"] = table.apply(
        lambda row:
            determine_age(
                field_raw=row["21022-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine Body Mass Index (BMI).
    table["body"] = table.apply(
        lambda row:
            determine_body_mass_index(
                field_21001=row["21001-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine self report of ancestry and or ethnicity.
    table["ancestry_self_white"] = table.apply(
        lambda row:
            determine_self_report_ancestry_ethnicity(
                field_21000=row["21000-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    # Introduce stratification variables by values of continuous variables, in
    # particular age within female and male persons.
    table = define_ordinal_stratifications_by_sex_continuous_variables(
        index="eid",
        sex_text="sex_text",
        variables=["age", "body"],
        table=table,
        report=report,
    )
    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=["body"],
        table=table,
    )

    # Genotype arrays and batches.
    table["genotype_batch"] = table.apply(
        lambda row:
            determine_genotype_batch(
                field_22000=row["22000-0.0"],
                codes_interpretations=source["field_22000"],
            ),
        axis="columns", # apply function to each row
    )
    table["genotype_array"] = table.apply(
        lambda row:
            determine_genotype_array(
                field_22000=row["22000-0.0"],
                codes_interpretations=source["field_22000"],
            ),
        axis="columns", # apply function to each row
    )
    table["genotype_array_axiom"] = table.apply(
        lambda row:
            determine_genotype_array_axiom_logical_binary(
                genotype_array=row["genotype_array"],
            ),
        axis="columns", # apply function to each row
    )
    table["genotype_availability"] = table.apply(
        lambda row:
            determine_genotype_availability(
                genotype_array=row["genotype_array"],
            ),
        axis="columns", # apply function to each row
    )

    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "31-0.0", "50-0.0", "54-0.0", "55-0.0", "22000-0.0",
            "21001-0.0", "21002-0.0", "21022-0.0", "22001-0.0",
            "23104-0.0",
        ],
        axis="columns",
        inplace=True
    )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "site", "region",
        "month_name", "month_order", "season", "day_length",
        "sex_y", "sex_x", "sex_text",
        "age", "age_grade_female", "age_grade_male",
        "body", "body_grade_female", "body_grade_male", "body_log",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: organize_sex_age_body_variables()")
        utility.print_terminal_partition(level=3)
        print(table_report)
        utility.print_terminal_partition(level=3)
        # Genotype arrays and batches.
        report_genotype_arrays_batches(table=table)
        # Assessment season and region.
        report_assessment_season_region(
            column_name="season",
            table=table,
        )
        report_assessment_season_region(
            column_name="day_length",
            table=table,
        )
        report_assessment_season_region(
            column_name="region",
            table=table,
        )
        # Variable types.
        utility.print_terminal_partition(level=2)
        print("After type conversion")
        print(table_report.dtypes)
        utility.print_terminal_partition(level=3)

        table_female = table_report.loc[
            (table_report["sex_text"] == "female"), :
        ]
        table_female_young = table_female.loc[
            (table_female["age_grade_female"] == 0), :
        ]
        table_female_old = table_female.loc[
            (table_female["age_grade_female"] == 2), :
        ]
        age_mean_female_young = numpy.nanmean(
            table_female_young["age"].to_numpy()
        )
        age_mean_female_old = numpy.nanmean(table_female_old["age"].to_numpy())
        print("mean age in young females: " + str(age_mean_female_young))
        print("mean age in old females: " + str(age_mean_female_old))

        table_male = table_report.loc[
            (table_report["sex_text"] == "male"), :
        ]
        table_male_young = table_male.loc[
            (table_male["age_grade_male"] == 0), :
        ]
        table_male_old = table_male.loc[
            (table_male["age_grade_male"] == 2), :
        ]
        age_mean_male_young = numpy.nanmean(table_male_young["age"].to_numpy())
        age_mean_male_old = numpy.nanmean(table_male_old["age"].to_numpy())
        print("mean age in young males: " + str(age_mean_male_young))
        print("mean age in old males: " + str(age_mean_male_old))

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Sex hormones
# Review:
# 21 May 2021: TCW verified that the minimum values for imputation match the
# minimum values on UK Biobank Data Showcase for albumin, steroid globulin,
# oestradiol, testosterone, and vitamin D.
# 20 May 2021: TCW verified formulas for estimation of free testosterone,
# bioavailable testosterone, free oestradiol, and bioavailable oestradiol.
# 11 March 2021: TCW verified formulas for estimation of free oestradiol
# 11 March 2021: TCW verified formulas for estimation of free testosterone
# 11 March 2021: TCW verified UK Biobank fields and their codings.

# TODO: also implement the natural log definition of bioavailable testosterone...
# Chung, Pathology Informatics, 2017 (PubMed:28828199)


def define_biochemistry_fields_measurement_reportability_missingness():
    """
    Define data fields in U.K. Biobank for measurement, reportability, and
    missingness of hormones.

    arguments:

    raises:

    returns:
        (list<dict>): data fields for measurement, reportability, and
            missingness for each hormone

    """

    records = list()

    record = dict()
    record["name"] = "albumin"
    record["measurement"] = "30600-0.0"
    record["reportability"] = "30606-0.0"
    record["missingness"] = "30605-0.0"
    records.append(record)

    record = dict()
    record["name"] = "steroid_globulin"
    record["measurement"] = "30830-0.0"
    record["reportability"] = "30836-0.0"
    record["missingness"] = "30835-0.0"
    records.append(record)

    record = dict()
    record["name"] = "testosterone"
    record["measurement"] = "30850-0.0"
    record["reportability"] = "30856-0.0"
    record["missingness"] = "30855-0.0"
    records.append(record)

    record = dict()
    record["name"] = "oestradiol"
    record["measurement"] = "30800-0.0"
    record["reportability"] = "30806-0.0"
    record["missingness"] = "30805-0.0"
    records.append(record)

    record = dict()
    record["name"] = "vitamin_d"
    record["measurement"] = "30890-0.0"
    record["reportability"] = "30896-0.0"
    record["missingness"] = "30895-0.0"
    records.append(record)

    record = dict()
    record["name"] = "cholesterol"
    record["measurement"] = "30690-0.0"
    record["reportability"] = "30696-0.0"
    record["missingness"] = "30695-0.0"
    records.append(record)

    #record = dict()
    #record["name"] = "c_reactive_protein"
    #record["measurement"] = "30710-0.0"
    #record["reportability"] = "30716-0.0"
    #record["missingness"] = "30715-0.0"
    #records.append(record)

    #record = dict()
    #record["name"] = "rheumatoid_factor"
    #record["measurement"] = "30820-0.0"
    #record["reportability"] = "30826-0.0"
    #record["missingness"] = "30825-0.0"
    #records.append(record)

    # Return information.
    return records


# TODO: TCW; 07 July 2022
# TODO: "convert_biochemical_concentration_units_moles_per_liter()"
# TODO: use 66,500 g/mol for Albumin?

# review: TCW, 09 February 2022
def convert_hormone_concentration_units_moles_per_liter(
    table=None,
    factors_concentration=None,
):
    """
    Converts hormone concentrations to units of moles per liter (mol/L).

    Molar Mass, Molecular Weight
    Species     ...     Molar Mass     ...     Reference
    estradiol           272.4 g/mol            PubChem
    testosterone        288.4 g/mol            PubChem
    SHBG
    albumin             69,367 g/mol           UniProt
    albumin             66,500 g/mol           Drug Bank
    - use 66.5 kDa molar mass for albumin (anticipate post-translational
    - - cleavage)

    Metric Prefixes
    (https://www.nist.gov/pml/weights-and-measures/metric-si-prefixes)
    Prefix     ...     Abbreviation     ...     Factor
    deci               d                        1E-1
    centi              c                        1E-2
    milli              m                        1E-3
    micro              u                        1E-6
    nano               n                        1E-9
    pico               p                        1E-12

    UK Biobank data-field 30600, concentration in grams per liter (g/L) of
    albumin in blood

    UK Biobank data-field 30830, concentration in nanomoles per liter (nmol/L)
    of steroid hormone binding globulin (SHBG) in blood

    UK Biobank data-field 30800, concentration in picomoles per liter (pmol/L)
    of oestradiol in blood

    UK Biobank data-field 30850, concentration in nanomoles per liter (nmol/L)
    of total testosterone in blood

    UK Biobank data-field 30890, concentration in nanomoles per liter (nmol/L)
    of vitamin D in blood

    UK Biobank data-field 30690, concentration in millimoles per liter (mmol/L)
    of cholesterol in blood

    https://www.nist.gov/pml/weights-and-measures/metric-si-prefixes

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of visualization and analysis

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert concentrations to units of moles per liter (mol/L).
    table["albumin"] = table.apply(
        lambda row: float(
            (row["30600-0.0"] / 66472.2) * factors_concentration["albumin"]
        ),
        axis="columns", # apply function to each row
    ) # 1 mole = 66472.2 g = 66.5 kDa
    table["steroid_globulin"] = table.apply(
        lambda row: float(
            (row["30830-0.0"] / 1E9) * factors_concentration["steroid_globulin"]
        ),
        axis="columns", # apply function to each row
    ) # 1 mol = 1E9 nanomole
    table["oestradiol"] = table.apply(
        lambda row: float(
            (row["30800-0.0"] / 1E12) * factors_concentration["oestradiol"]
        ),
        axis="columns", # apply function to each row
    ) # 1 mol = 1E12 picomole
    table["testosterone"] = table.apply(
        lambda row: float(
            (row["30850-0.0"] / 1E9) * factors_concentration["testosterone"]
        ),
        axis="columns", # apply function to each row
    ) # 1 mol = 1E9 nanomole
    table["vitamin_d"] = table.apply(
        lambda row: float(
            (row["30890-0.0"] / 1E9) * factors_concentration["vitamin_d"]
        ),
        axis="columns", # apply function to each row
    ) # 1 mol = 1E9 nanomole
    table["cholesterol"] = table.apply(
        lambda row: float(
            (row["30690-0.0"] / 1E3) * factors_concentration["cholesterol"]
        ),
        axis="columns", # apply function to each row
    ) # 1 mol = 1E3 millimole
    # Return information.
    return table


# Imputation of missing biochemical measurements


# review: TCW, 09 February 2022
def interpret_biochemistry_missingness_beyond_detection_range(
    field_code_782=None,
):
    """
    Intepret UK Biobank's coding 782 for multiple data fields.

    Only interpret whether the variable indicates that the measurement was
    missing due to a measurement that was beyond the detection range (either
    less than or greater than).

    0: measurement was missing because of other problem
    1: measurement was not reportable and missing because value was less than or
    greater than the detection range

    Data-Field "30605": "Albumin missing reason"
    Data-Field "30835": "SHBG missing reason"
    Data-Field "30805": "Oestradiol missing reason"
    Data-Field "30855": "Testosterone missing reason"
    Data-Field "30895": "Vitamin D missing reason"
    Data-Field "30695": "Cholesterol missing reason"
    Data-Field "30715": "C-reactive protein missing reason"
    Data-Field "30825": "Rheumatoid factor missing reason"
    UK Biobank data-coding "782" for variable fields.
    1: "No data returned"
    2: "Original value above or below reportable limit"
    3: "Unrecoverable aliquot problem (dip)"
    4: "Unrecoverable aliquot problem (possible aliquot mis-classification)"
    5: "Aliquot 4 used"
    6: ?
    7: "Analyser deemed result not reportable for reason other than above or
       below reportable range"
    8: "Not reportable because dilution factor 0.9-0.99 or 1.01-1.1"
    9: "Not reportable because dilution factor <0.9 or >1.1"

    Accommodate inexact float values.

    arguments:
        field_code_4917 (float): UK Biobank fields in coding 4917, representing
            biochemistry reportability

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_code_782)) and
        (0.5 <= field_code_782 and field_code_782 < 9.5)
    ):
        # The variable has a valid value.
        if (0.5 <= field_code_782 and field_code_782 < 1.5):
            # 1: "No data returned"
            value = 0
        elif (1.5 <= field_code_782 and field_code_782 < 2.5):
            # 2: "Original value above or below reportable limit"
            value = 1
        elif (2.5 <= field_code_782 and field_code_782 < 3.5):
            # 3: "Unrecoverable aliquot problem (dip)"
            value = 0
        elif (3.5 <= field_code_782 and field_code_782 < 4.5):
            # 4: "Unrecoverable aliquot problem (possible aliquot
            # mis-classification)"
            value = 0
        elif (4.5 <= field_code_782 and field_code_782 < 5.5):
            # 5: "Aliquot 4 used"
            value = 0
        elif (5.5 <= field_code_782 and field_code_782 < 6.5):
            # Code "6" does not have an explanation.
            # 6: ""
            value = float("nan")
        elif (6.5 <= field_code_782 and field_code_782 < 7.5):
            # 7: "Analyser deemed result not reportable for reason other than
            # above or below reportable range"
            value = 0
        elif (7.5 <= field_code_782 and field_code_782 < 8.5):
            # 8: "Not reportable because dilution factor 0.9-0.99 or 1.01-1.1"
            value = 0
        elif (8.5 <= field_code_782 and field_code_782 < 9.5):
            # 9: "Not reportable because dilution factor <0.9 or >1.1"
            value = 0
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


# review: TCW, 09 February 2022
def interpret_biochemistry_reportability_detection_limit(
    field_code_4917=None,
):
    """
    Intepret UK Biobank's coding 4917 for multiple data fields.

    Interpret whether the variable indicates that the measurement was not
    reportable and missing due to a measurement that was less than or greater
    than the limit of detection.

    0: measurement was reportable and not missing
    1: measurement was not reportable and missing because value was less than
    detection limit
    2: measurement was not reportable and missing because value was greater than
    detection limit
    NA: not reportable and missing without definitive explanation (per 4917)

    Data-Field "30606": "Albumin reportability"
    Data-Field "30836": "SHBG reportability"
    Data-Field "30806": "Oestradiol reportability"
    Data-Field "30856": "Testosterone reportability"
    Data-Field "30896": "Vitamin D reportability"
    Data-Field "30696": "Cholesterol reportability"
    Data-Field "30716": "C-reactive protein reportability"
    Data-Field "30826": "Rheumatoid factor reportability"

    UK Biobank data-coding "4917" for variable fields.
    1: "Reportable at assay and after aliquot correction, if attempted"
    2: "Reportable at assay but not reportable after any corrections (too low)"
    3: "Reportable at assay but not reportable after any corrections (too high)"
    4: "Not reportable at assay (too low)"
    5: "Not reportable at assay (too high)"

    Note: "Absence of a value here means that a result has not been attempted."

    Accommodate inexact float values.

    arguments:
        field_code_4917 (float): UK Biobank fields in coding 4917, representing
            biochemistry reportability

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_code_4917)) and
        (0.5 <= field_code_4917 and field_code_4917 < 5.5)
    ):
        # The variable has a valid value.
        if (0.5 <= field_code_4917 and field_code_4917 < 1.5):
            # 1: "Reportable at assay and after aliquot correction, if
            # attempted"
            value = 0
        elif (1.5 <= field_code_4917 and field_code_4917 < 2.5):
            # 2: "Reportable at assay but not reportable after any corrections
            # (too low)"
            value = 1
        elif (2.5 <= field_code_4917 and field_code_4917 < 3.5):
            # 3: "Reportable at assay but not reportable after any corrections
            # (too high)"
            value = 2
        elif (3.5 <= field_code_4917 and field_code_4917 < 4.5):
            # 4: "Not reportable at assay (too low)"
            value = 1
        elif (4.5 <= field_code_4917 and field_code_4917 < 5.5):
            # 5: "Not reportable at assay (too high)"
            value = 2
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def determine_hormones_missingness_beyond_detection_range(
    hormone_fields=None,
    table=None,
):
    """
    Determine for each hormone whether missingness variable indicates that the
    measurement was not reportable because it was beyond the detection range.

    arguments:
        hormone_fields (list<dict>): data fields for measurement, reportability,
            and missingness for each hormone
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    table = table.copy(deep=True)
    # Determine missingness.
    for hormone in hormone_fields:
        missingness_range = str(
            str(hormone["name"]) + "_missingness_range"
        )
        table[missingness_range] = table.apply(
            lambda row:
                interpret_biochemistry_missingness_beyond_detection_range(
                    field_code_782=row[hormone["missingness"]],
                ),
            axis="columns", # apply function to each row
        )
        pass
    # Return information.
    return table


def determine_hormones_reportability_detection_limit(
    hormone_fields=None,
    table=None,
):
    """
    Determine for each hormone whether reportability variable indicates that the
    measurement was not reportable because it was less than detection limit.

    arguments:
        hormone_fields (list<dict>): data fields for measurement, reportability,
            and missingness for each hormone
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    table = table.copy(deep=True)
    # Determine reportability.
    for hormone in hormone_fields:
        reportability = str(
            str(hormone["name"]) + "_reportability_limit"
        )
        table[reportability] = table.apply(
            lambda row:
                interpret_biochemistry_reportability_detection_limit(
                    field_code_4917=row[hormone["reportability"]],
                ),
            axis="columns", # apply function to each row
        )
        pass
    # Return information.
    return table


def determine_hormone_binary_detection(
    hormone=None,
    measurement=None,
    missingness_range=None,
    reportability_limit=None,
):
    """
    Organizes binary representation of whether hormone was measurable or was
    missing due to the detection limit.

    1: hormone has a valid measurement
    0: hormone has a missing measurement due to the detection range (above or
    below)
    NA: hormone has a missing measurement without explicit annotation that the
    reason was detection limit

    Note:
    - hormones tend to have missing values either above or below the detection
    range (not both)

    arguments:
        hormone (str): name of a hormone
        measurement (float): value of hormone's actual measurement
        missingness_range (int): binary representation of whether measurement's
            value was missing because the value was out of detection range
        reportability_limit (int): integer categorical representation of whether
            measurement's value was missing because the value was less than or
            greater than the limit of detection

    raises:

    returns:
        (float): value for hormone's measurement

    """

    # Interpret field code.
    if (pandas.isna(measurement)):
        # Hormone's measurement has a missing value.
        if (
            (not pandas.isna(missingness_range)) and
            (missingness_range == 1) and
            (not pandas.isna(reportability_limit)) and
            ((reportability_limit == 1) or (reportability_limit == 2))
        ):
            # Record has valid explanations for the missing value of hormone's
            # measurement.
            # Measurement's missingness is due to detection range.
            value = 0
        else:
            # Without valid explanations of missing value for hormone's
            # measurement, do not assign definitive value.
            value = float("nan")
    else:
        # Hormone's measurement has a non-missing value.
        value = 1
    # Return information.
    return value


def determine_hormones_binary_detection(
    hormones=None,
    table=None,
    report=None,
):
    """
    Organizes binary representation of whether hormone was measurable or was
    missing due to the detection limit.

    arguments:
        hormones (list<str>): names of hormones for imputation
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Iterate on hormones.
    for hormone in hormones:
        # Determine column names.
        detection = str(str(hormone) + "_detection")
        missingness_range = str(str(hormone) + "_missingness_range")
        reportability_limit = str(str(hormone) + "_reportability_limit")
        # Impute missing measurements.
        table[detection] = table.apply(
            lambda row:
                determine_hormone_binary_detection(
                    hormone=hormone,
                    measurement=row[hormone],
                    missingness_range=row[missingness_range],
                    reportability_limit=row[reportability_limit],
                ),
            axis="columns", # apply function to each row
        )
        # Report.
        if report:
            # Column name translations.
            utility.print_terminal_partition(level=3)
            print("report: determine_hormones_binary_detection()")
            print("hormone: " + str(hormone))
        pass
    # Return information.
    return table


def impute_missing_hormone_detection_limit(
    hormone=None,
    sex_text=None,
    minima=None,
    maxima=None,
    measurement=None,
    missingness_range=None,
    reportability_limit=None,
):
    """
    Determine imputation value for a hormone.

    If record has a definitive explanation that the hormone's measurement was
    missing because its value was less than the assay's limit of detection, then
    impute to half the sex-specific minimal value for that hormone.
    (missingness_range = 1) and (reportability_limit = 1)

    If record has a definitive explanation that the hormone's measurement was
    missing because its value was greater than the assay's limit of detection,
    then impute to the sex-specific maximal value for that hormone.
    (missingness_range = 1) and (reportability_limit = 2)

    arguments:
        hormone (str): name of a hormone
        sex_text (str): textual representation of sex selection
        minima (dict<float>): sex-specific minimal values of hormone's
            measurements
        maxima (dict<float>): sex-specific maximal values of hormone's
            measurements
        measurement (float): value of hormone's actual measurement
        missingness_range (int): binary representation of whether measurement's
            value was missing because the value was out of detection range
        reportability_limit (int): integer categorical representation of whether
            measurement's value was missing because the value was less than or
            greater than the limit of detection

    raises:

    returns:
        (float): value for hormone's measurement

    """

    # Interpret field code.
    if (pandas.isna(measurement)):
        # Hormone's measurement has a missing value.
        if (
            (not pandas.isna(missingness_range)) and
            (missingness_range == 1) and
            (not pandas.isna(reportability_limit)) and
            ((reportability_limit == 1) or (reportability_limit == 2))
        ):
            # Record has valid explanations for the missing value of hormone's
            # measurement.
            if (reportability_limit == 1):
                # Hormone's measurement was less than the limit of detection.
                # Impute to half of the hormone's sex-specific minimal value of
                # measurement.
                minimum = minima[sex_text]
                value = float(minimum / 2)
            elif (reportability_limit == 2):
                # Hormone's measurement was greater than the limit of detection.
                # Impute to the hormone's sex-specific maximal value of
                # measurement.
                maximum = maxima[sex_text]
                value = maximum
            else:
                # Uninterpretable.
                value = measurement
        else:
            # Without valid explanations of missing value for hormone's
            # measurement, do not impute.
            value = measurement
    else:
        # Hormone's measurement has a not missing value.
        value = measurement
    # Return information.
    return value


def impute_missing_hormones_detection_limit(
    hormones=None,
    table=None,
    report=None,
):
    """
    Organizes imputation of missing measurements for hormones.

    arguments:
        hormones (list<str>): names of hormones for imputation
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Stratify table by sex to determine sex-specific minimal and maximal
    # values.
    table_female = table.loc[
        (table["sex_text"] == "female"), :
    ]
    table_male = table.loc[
        (table["sex_text"] == "male"), :
    ]
    # Iterate on hormones.
    for hormone in hormones:
        # Determine column names.
        imputation = str(str(hormone) + "_imputation")
        missingness_range = str(str(hormone) + "_missingness_range")
        reportability_limit = str(str(hormone) + "_reportability_limit")
        # Determine sex-specific minimal and maximal values.
        minima = dict()
        minima["female"] = numpy.nanmin(table_female[hormone].to_numpy())
        minima["male"] = numpy.nanmin(table_male[hormone].to_numpy())
        maxima = dict()
        maxima["female"] = numpy.nanmax(table_female[hormone].to_numpy())
        maxima["male"] = numpy.nanmax(table_male[hormone].to_numpy())
        # Impute missing measurements.
        table[imputation] = table.apply(
            lambda row:
                impute_missing_hormone_detection_limit(
                    hormone=hormone,
                    minima=minima,
                    maxima=maxima,
                    sex_text=row["sex_text"],
                    measurement=row[hormone],
                    missingness_range=row[missingness_range],
                    reportability_limit=row[reportability_limit],
                ),
            axis="columns", # apply function to each row
        )
        # Report.
        if report:
            # Column name translations.
            utility.print_terminal_partition(level=3)
            print("report: impute_missing_hormones_detection_limit()")
            print("hormone: " + str(hormone))
            print("minimum female: " + str(minima["female"]))
            print("minimum male: " + str(minima["male"]))
            print("maximum female: " + str(maxima["female"]))
            print("maximum male: " + str(maxima["male"]))
        pass
    # Return information.
    return table


def organize_hormones_missingness_detection_imputation(
    hormone_fields=None,
    table=None,
    report=None,
):
    """
    Organizes imputation of missing measurements for hormones.

    arguments:
        hormone_fields (list<dict>): data fields for measurement, reportability,
            and missingness for each hormone
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Interpret missingness and reportability explanation variables for
    # hormones.
    table = determine_hormones_missingness_beyond_detection_range(
        hormone_fields=hormone_fields,
        table=table,
    )
    table = determine_hormones_reportability_detection_limit(
        hormone_fields=hormone_fields,
        table=table,
    )
    # Impute missing values.
    hormones = list()
    for record in hormone_fields:
        hormones.append(record["name"])
    table = determine_hormones_binary_detection(
        hormones=hormones,
        table=table,
        report=report,
    )
    table = impute_missing_hormones_detection_limit(
        hormones=hormones,
        table=table,
        report=report,
    )
    # Return information.
    return table


# Cohort-specific ordinal representations


def determine_cohort_specific_hormone_ordinal(
    missingness_range=None,
    reportability_limit=None,
    measurement=None,
    median=None,
):
    """
    Determine for a hormone the cohort-specific ordinal representation with
    consideration of cohort-specific median and reportability due to detection
    limit.

    Ordinal code:
    0: below limit of detection
    1: measurement in range, below cohort-specific median
    2: measurement in range, equal to or above cohort-specific median
    3: above limit of detection

    arguments:
        missingness_range (int): binary representation of whether measurement's
            value was missing because the value was out of detection range
        reportability_limit (int): integer categorical representation of whether
            measurement's value was missing because the value was less than or
            greater than the limit of detection
        measurement (float): value of hormone's actual measurement
        median (float): cohort-specific median value of hormone's measurements
        missingness_range (int): binary representation of whether measurement's
            value was missing because the value was out of detection range

    raises:

    returns:
        (float): value for hormone's measurement

    """

    # Interpret field code.
    if (pandas.isna(measurement)):
        # Hormone's measurement has a missing value.
        if (
            (not pandas.isna(missingness_range)) and
            (missingness_range == 1) and
            (not pandas.isna(reportability_limit)) and
            ((reportability_limit == 1) or (reportability_limit == 2))
        ):
            # Record has valid explanations for the missing value of hormone's
            # measurement.
            if (reportability_limit == 1):
                # Hormone's measurement was less than the limit of detection.
                value = 0
            elif (reportability_limit == 2):
                # Hormone's measurement was greater than the limit of detection.
                value = 3
            else:
                # Uninterpretable.
                value = float("nan")
        else:
            # Without valid explanations of missing value for hormone's
            # measurement, do not assign ordinal value.
            value = float("nan")
    else:
        # Hormone's measurement has a not missing value.
        if (
            (not pandas.isna(median)) and
            (measurement < median)
        ):
            # Measurement is less than median.
            value = 1
        elif (
            (not pandas.isna(median)) and
            (measurement >= median)
        ):
            # Measurement is greater than or equal to median.
            value = 2
        else:
            # Uninterpretable.
            value = float("nan")
    # Return information.
    return value


def organize_cohort_specific_hormone_ordinal_representation(
    hormone=None,
    cohort=None,
    table_cohort=None,
    table=None,
):
    """
    Determine for a hormone the cohort-specific ordinal representation with
    consideration of cohort-specific median and reportability due to detection
    limit.

    It is important that the hormone's median value is cohort-specific to
    account for influences from sex and female menopause.

    arguments:
        hormone (str): name of a hormone
        cohort (str): name of a cohort
        table_cohort (object): Pandas data frame for a single cohort
        table (object): Pandas data frame collecting information for all cohorts
            and hormones

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information.
    table_cohort = table_cohort.copy(deep=True)
    table = table.copy(deep=True)
    # Define names of relevant columns.
    #hormone_name = hormone.replace("_", "-")
    column_hormone = str(hormone)
    column_missingness = str(str(hormone) + "_missingness_range")
    column_reportability = str(str(hormone) + "_reportability_limit")
    column_order = str(str(hormone) + "_" + str(cohort) + "_order")
    # Calculate cohort-specific median.
    array = copy.deepcopy(table_cohort[column_hormone].dropna().to_numpy())
    count = int(array.size)
    if (count > 10):
        median = numpy.nanmedian(array)
    else:
        median = float("nan")
    # Determine cohort-specific ordinal representation.
    table_cohort[column_order] = table_cohort.apply(
        lambda row:
            determine_cohort_specific_hormone_ordinal(
                missingness_range=row[column_missingness],
                reportability_limit=row[column_reportability],
                measurement=row[column_hormone],
                median=median,
            ),
        axis="columns", # apply function to each row
    )
    # Merge cohort table back to the collection table.
    table_cohort = table_cohort.loc[
        :, table_cohort.columns.isin([column_order])
    ]
    table = table.merge(
        table_cohort,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_original", "_cohort"),
    )
    # Return information.
    return table


def organize_cohorts_specific_hormones_ordinal_representations(
    table=None,
    report=None,
):
    """
    Organize the definition of cohort-specific ordinal representations of
    hormones with consideration of cohort-specific medians and reportabilities
    due to detection limits.

    It is important that each hormone's median value is cohort-specific to
    account for influences from sex and female menopause.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Define relevant cohorts.
    cohorts_records = ukb_strat.stratify_set_primary_sex_age_body_menopause(
        table=table
    )

    # Define relevant hormones.
    # These ordinal representations are only relevant to the raw hormones with
    # direct measurements.
    # Exclude estimates of bioavailable and free hormones.
    hormones = [
        "vitamin_d", "albumin", "steroid_globulin",
        "oestradiol", "testosterone",
    ]

    # Iterate across cohorts.
    for cohort_record in cohorts_records:
        cohort = cohort_record["cohort"]

        # Iterate across hormones.
        for hormone in hormones:
            table = organize_cohort_specific_hormone_ordinal_representation(
                hormone=hormone,
                cohort=cohort,
                table_cohort=cohort_record["table"],
                table=table,
            )
            pass
        pass

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_cohorts_specific_hormones_ordinal_representations()")
        utility.print_terminal_partition(level=3)
        # Organize table.
        columns_report = list()
        hormone = "oestradiol"
        #hormone_name = hormone.replace("_", "-")
        columns_report.append(hormone)
        for cohort_record in cohorts_records:
            columns_report.append(
                str(
                    str(hormone) + "_" +
                    str(cohort_record["cohort"]) + "_order"
                )
            )
        columns_report.insert(0, "IID")
        #columns_report.insert(0, "eid")
        table_report = table.copy(deep=True)
        table_report = table_report.loc[
            :, table_report.columns.isin(columns_report)
        ]
        table_report = table_report[[*columns_report]]
        utility.print_terminal_partition(level=3)
        print("Here is table with cohort-specific ordinal variables...")
        print(table_report)
    # Return information.
    return table


# Estimation of bioavailable and free hormones


def calculate_estimation_free_testosterone(
    testosterone=None,
    albumin=None,
    steroid_globulin=None,
    factors_concentration=None,
    associations=None,
):
    """
    Calculates an estimation of free testosterone (neither bound to albumin nor
    steroid hormone binding globulin).

    This function applies the formula of Sodergard, Journal of Steroid
    Biochemistry, 1982 (PubMed:7202083).

    The specific structure of this formula appears in the following articles:
    van den Beld, Clinical Endocrinology & Metabolism, 2000 (PubMed:10999822)
    de Ronde, Clinical Endocrinology & Metabolism, 2005 (PubMed:15509641)
    de Ronde, Clinical Chemistry, 2006 (PubMed:16793931)
    Chung, Pathology Informatics, 2017 (PubMed:28828199)

    The structure of this formula is slightly different, though presumably
    equivalent in the following articles:
    Vermeulen, Clinical Endocrinology & Metabolism, 1999 (PubMed:10523012)
    Emadi-Konjin, Clinical Biochemistry, 2003 (PubMed:14636872)

    Association constant for steroid hormone binding globulin (SHBG) to
    testosterone:
    5.97E8 - 1.4E9 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    Association constant for albumin to testosterone:
    1.3E4 - 4.06E4 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    arguments:
        testosterone (float): concentration in moles per liter (mol/L) of total
            testosterone in blood, with concentration factor
        albumin (float): concentration in moles per liter (mol/L) of albumin in
            blood, with concentration factor
        steroid_globulin (float): concentration in moles per liter (mol/L) of
            steroid hormone binding globulin in blood, with concentration factor
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of float storage and analysis
        associations (dict<float>): association constants in liters per mole for
            steroid hormone binding globulin (SHBG) and ablumin

    raises:

    returns:
        (float): estimate concentration in moles per liter of free testosterone,
            with concentration factor

    """

    # Raw values are stored with a concentration factor.
    # Convert concentrations to unit of moles per liter (mol/L).
    testosterone_unit = float(
        testosterone / factors_concentration["testosterone"]
    )
    albumin_unit = float(albumin / factors_concentration["albumin"])
    steroid_globulin_unit = float(
        steroid_globulin / factors_concentration["steroid_globulin"]
    )

    # Calculate simplification variables: a and b.
    a = (
        associations["albumin_test"] +
        associations["shbg_test"] +
        (
            (associations["albumin_test"] * associations["shbg_test"]) *
            (steroid_globulin_unit + albumin_unit - testosterone_unit)
        )
    ) # unit: L/mol
    b = (
        1 +
        (associations["shbg_test"] * steroid_globulin_unit) +
        (associations["albumin_test"] * albumin_unit) -
        (
            (associations["albumin_test"] + associations["shbg_test"]) *
            testosterone_unit
        )
    ) # unit: 1
    # Calculate free testosterone.
    testosterone_free_unit = (
        (-1 * b) + math.sqrt(math.pow(b, 2) + (4 * a * testosterone_unit))
    ) / (2 * a) # unit: mol/L
    # Convert units by factor.
    testosterone_free = float(
        testosterone_free_unit * factors_concentration["testosterone_free"]
    )
    # Return information.
    return testosterone_free


def calculate_estimation_bioavailable_testosterone(
    testosterone_free=None,
    albumin=None,
    factors_concentration=None,
    associations=None,
):
    """
    Calculates an estimation of bioavailable testosterone (not bound to steroid
    hormone binding globulin). Testosterone binds to SHBG more tightly
    (association constant) than it does to albumin. Testosterone bound to
    albumin is more "bioavailable". Hence the formula includes this testosterone
    bound to albumin in addition to the "free" testosterone.

    This function applies the formula of Sodergard, Journal of Steroid
    Biochemistry, 1982 (PubMed:7202083).

    The specific structure of this formula appears in the following articles:
    van den Beld, Clinical Endocrinology & Metabolism, 2000 (PubMed:10999822)
    de Ronde, Clinical Endocrinology & Metabolism, 2005 (PubMed:15509641)
    de Ronde, Clinical Chemistry, 2006 (PubMed:16793931)

    The structure of this formula is slightly different, though presumably
    equivalent in the following articles:
    Vermeulen, Clinical Endocrinology & Metabolism, 1999 (PubMed:10523012)
    Emadi-Konjin, Clinical Biochemistry, 2003 (PubMed:14636872)

    Association constant for steroid hormone binding globulin (SHBG) to
    testosterone:
    5.97E8 - 1.4E9 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    Association constant for albumin to testosterone:
    1.3E4 - 4.06E4 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    arguments:
        testosterone_free (float): estimate concentration in moles per liter
            (mol/L) of free testosterone, with concentration factor
        albumin (float): concentration in moles per liter (mol/L) of albumin in
            blood, with concentration factor
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of float storage and analysis
        associations (dict<float>): association constants in liters per mole for
            steroid hormone binding globulin (SHBG) and ablumin

    raises:

    returns:
        (float): estimate concentration in moles per liter of bioavailable
            testosterone, with concentration factor

    """

    # Raw values are stored with a concentration factor.
    # Convert concentrations to unit of moles per liter (mol/L).
    testosterone_free_unit = float(
        testosterone_free / factors_concentration["testosterone_free"]
    )
    albumin_unit = float(albumin / factors_concentration["albumin"])

    # Calculate simplification variables: numerator and denominator.
    numerator = (
        associations["albumin_test"] *
        albumin_unit *
        testosterone_free_unit
    ) # unit: mol/L
    denominator = (
        1 +
        (associations["albumin_test"] * testosterone_free_unit)
    ) # unit: 1
    # Calculate bioavailable testosterone.
    testosterone_bioavailable_unit = (
        (numerator / denominator) + testosterone_free_unit
    ) # unit: mol/L
    # Convert units by factor.
    testosterone_bioavailable = float(
        testosterone_bioavailable_unit *
        factors_concentration["testosterone_bioavailable"]
    )
    # Return information.
    return testosterone_bioavailable


def calculate_estimation_free_oestradiol(
    oestradiol=None,
    testosterone_free=None,
    albumin=None,
    steroid_globulin=None,
    factors_concentration=None,
    associations=None,
):
    """
    Calculates an estimation of free oestradiol (neither bound to albumin nor
    steroid hormone binding globulin).

    This function applies the formula of Sodergard, Journal of Steroid
    Biochemistry, 1982 (PubMed:7202083).

    The specific structure of this formula appears in the following articles:
    de Ronde, Clinical Endocrinology & Metabolism, 2005 (PubMed:15509641)

    Notice that there is a discrepancy with the formula in the following
    article:
    van den Beld, Clinical Endocrinology & Metabolism, 2000 (PubMed:10999822)
    This article seems to introduce an extra term for the concentration of SHBG
    in both the "b" and "c" formulas.
    This extra concentration term is likely erroneous as it would disrupt the
    units of the "b" and "c" formulas.
    It would also disrupt the final units (mol/L) of the total formula.
    Compare to the formula for free testosterone.

    Association constant for steroid hormone binding globulin (SHBG) to
    testosterone:
    5.97E8 - 1.4E9 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    Association constant for albumin to testosterone:
    1.3E4 - 4.06E4 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    Association constant for steroid hormone binding globulin (SHBG) to
    oestradiol:
    3.14E8 - 6.8E8 L/mol (PubMed:7202083, PubMed:7195404)

    Association constant for albumin to oestradiol:
    4.21E4 L/mol (PubMed:7202083)


    arguments:
        oestradiol (float): concentration in moles per liter (mol/L) of total
            oestradiol in blood, with concentration factor
        testosterone_free (float): estimate concentration in moles per liter
            (mol/L) of free testosterone, with concentration factor
        albumin (float): concentration in moles per liter (mol/L) of albumin in
            blood, with concentration factor
        steroid_globulin (float): concentration in moles per liter (mol/L) of
            steroid hormone binding globulin in blood, with concentration factor
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of float storage and analysis
        associations (dict<float>): association constants in liters per mole for
            steroid hormone binding globulin (SHBG) and ablumin

    raises:

    returns:
        (float): estimate concentration in moles per liter of free oestradiol,
            with concentration factor

    """

    # Raw values are stored with a concentration factor.
    # Convert concentrations to unit of moles per liter (mol/L).
    oestradiol_unit = float(
        oestradiol / factors_concentration["oestradiol"]
    )
    testosterone_free_unit = float(
        testosterone_free / factors_concentration["testosterone_free"]
    )
    albumin_unit = float(albumin / factors_concentration["albumin"])
    steroid_globulin_unit = float(
        steroid_globulin / factors_concentration["steroid_globulin"]
    )

    # Calculate simplification variables: a, b, and c
    a = (
        associations["shbg_oest"] *
        (1 + (associations["albumin_test"] * albumin_unit))
    ) # unit: L/mol
    b = (
        (oestradiol_unit * associations["shbg_oest"]) -
        (
            (1 + (associations["albumin_oest"] * albumin_unit)) *
            (1 + (associations["shbg_test"] * testosterone_free_unit))
        ) -
        (associations["shbg_oest"] * steroid_globulin_unit)
    ) # unit: 1
    c = (
        oestradiol_unit *
        (1 + (associations["shbg_test"] * testosterone_free_unit))
    ) # unit: mol/L
    # Calculate free oestradiol.
    oestradiol_free_unit = (
        (-1 * b) - math.sqrt(math.pow(b, 2) - (4 * a * c))
    ) / (2 * a)
    # Convert units by factor.
    oestradiol_free = float(
        oestradiol_free_unit * factors_concentration["oestradiol_free"]
    )
    # Return information.
    return oestradiol_free


def calculate_estimation_bioavailable_oestradiol(
    oestradiol=None,
    oestradiol_free=None,
    testosterone_free=None,
    steroid_globulin=None,
    factors_concentration=None,
    associations=None,
):
    """
    Calculates an estimation of free oestradiol (neither bound to albumin nor
    steroid hormone binding globulin).

    This function applies the formula of Sodergard, Journal of Steroid
    Biochemistry, 1982 (PubMed:7202083).

    The specific structure of this formula appears in the following articles:
    de Ronde, Clinical Endocrinology & Metabolism, 2005 (PubMed:15509641)

    Notice that there is a discrepancy with the formula in the following
    article:
    van den Beld, Clinical Endocrinology & Metabolism, 2000 (PubMed:10999822)
    This article seems to introduce an extra term for the concentration of SHBG
    in both the "b" and "c" formulas.
    This extra concentration term is likely erroneous as it would disrupt the
    units of the "b" and "c" formulas.
    It would also disrupt the final units (mol/L) of the total formula.
    Compare to the formula for free testosterone.

    Association constant for steroid hormone binding globulin (SHBG) to
    testosterone:
    5.97E8 - 1.4E9 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    Association constant for albumin to testosterone:
    1.3E4 - 4.06E4 L/mol (PubMed:7202083; PubMed:10523012; PubMed:14636872)

    Association constant for steroid hormone binding globulin (SHBG) to
    oestradiol:
    3.14E8 - 6.8E8 L/mol (PubMed:7202083, PubMed:7195404)

    Association constant for albumin to oestradiol:
    4.21E4 L/mol (PubMed:7202083)


    arguments:
        oestradiol (float): concentration in moles per liter (mol/L) of total
            oestradiol in blood, with concentration factor
        oestradiol_free (float): estimate concentration in moles per liter
            (mol/L) of free oestradiol, with concentration factor
        testosterone_free (float): estimate concentration in moles per liter
            (mol/L) of free testosterone, with concentration factor
        steroid_globulin (float): concentration in moles per liter (mol/L) of
            steroid hormone binding globulin in blood, with concentration factor
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of float storage and analysis
        associations (dict<float>): association constants in liters per mole for
            steroid hormone binding globulin (SHBG) and ablumin

    raises:

    returns:
        (float): estimate concentration in moles per liter of free oestradiol,
            with concentration factor

    """

    # Raw values are stored with a concentration factor.
    # Convert concentrations to unit of moles per liter (mol/L).
    oestradiol_unit = float(
        oestradiol / factors_concentration["oestradiol"]
    )
    oestradiol_free_unit = float(
        oestradiol_free / factors_concentration["oestradiol_free"]
    )
    testosterone_free_unit = float(
        testosterone_free / factors_concentration["testosterone_free"]
    )
    steroid_globulin_unit = float(
        steroid_globulin / factors_concentration["steroid_globulin"]
    )

    # Calculate simplification variables: numerator and denominator.
    numerator = (
        associations["shbg_oest"] *
        steroid_globulin_unit *
        oestradiol_free_unit
    ) # unit: mol/L
    denominator = (
        1 +
        (associations["shbg_oest"] * oestradiol_free_unit) +
        (associations["shbg_test"] * testosterone_free_unit)
    ) # unit: 1
    # Calculate bioavailable oestradiol.
    oestradiol_bioavailable_unit = (
        oestradiol_unit - (numerator / denominator)
    ) # unit: mol/L
    # Convert units by factor.
    oestradiol_bioavailable = float(
        oestradiol_bioavailable_unit *
        factors_concentration["oestradiol_bioavailable"]
    )
    # Return information.
    return oestradiol_bioavailable


def organize_calculation_estimate_bioavailable_free_hormones(
    factors_concentration=None,
    table=None,
    report=None,
):
    """
    Organizes calculation estimates of bioavailable and free concentrations
    of testosterone and oestradiol.

    arguments:
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of float storage and analysis
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Define association constants for calculation of free hormones.
    # units: L/mol
    associations = dict()
    associations["shbg_test"] = 5.97E8 #
    associations["shbg_oest"] = 3.14E8 #
    associations["albumin_test"] = 4.06E4 #
    associations["albumin_oest"] = 4.21E4 #

    # Measurements.

    # Calculate estimation of free and bioavailable testosterone.
    table["testosterone_free"] = table.apply(
        lambda row:
            calculate_estimation_free_testosterone(
                testosterone=row["testosterone"],
                albumin=row["albumin"],
                steroid_globulin=row["steroid_globulin"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    table["testosterone_bioavailable"] = table.apply(
        lambda row:
            calculate_estimation_bioavailable_testosterone(
                testosterone_free=row["testosterone_free"],
                albumin=row["albumin"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    # Calculate estimation of free, bioavailable oestradiol.
    table["oestradiol_free"] = table.apply(
        lambda row:
            calculate_estimation_free_oestradiol(
                oestradiol=row["oestradiol"],
                testosterone_free=row["testosterone_free"],
                albumin=row["albumin"],
                steroid_globulin=row["steroid_globulin"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    table["oestradiol_bioavailable"] = table.apply(
        lambda row:
            calculate_estimation_bioavailable_oestradiol(
                oestradiol=row["oestradiol"],
                oestradiol_free=row["oestradiol_free"],
                testosterone_free=row["testosterone_free"],
                steroid_globulin=row["steroid_globulin"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )

    # Imputations.

    # Calculate estimation of free and bioavailable testosterone.
    table["testosterone_free_imputation"] = table.apply(
        lambda row:
            calculate_estimation_free_testosterone(
                testosterone=row["testosterone_imputation"],
                albumin=row["albumin_imputation"],
                steroid_globulin=row["steroid_globulin_imputation"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    table["testosterone_bioavailable_imputation"] = table.apply(
        lambda row:
            calculate_estimation_bioavailable_testosterone(
                testosterone_free=row["testosterone_free_imputation"],
                albumin=row["albumin_imputation"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    # Calculate estimation of free, bioavailable oestradiol.
    table["oestradiol_free_imputation"] = table.apply(
        lambda row:
            calculate_estimation_free_oestradiol(
                oestradiol=row["oestradiol_imputation"],
                testosterone_free=row["testosterone_free_imputation"],
                albumin=row["albumin_imputation"],
                steroid_globulin=row["steroid_globulin_imputation"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    table["oestradiol_bioavailable_imputation"] = table.apply(
        lambda row:
            calculate_estimation_bioavailable_oestradiol(
                oestradiol=row["oestradiol_imputation"],
                oestradiol_free=row["oestradiol_free_imputation"],
                testosterone_free=row["testosterone_free_imputation"],
                steroid_globulin=row["steroid_globulin_imputation"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )

    # Return information.
    return table


# Hormonal medications


def organize_medication_code_atc_class_table(
    table_wu=None,
    table_appleby=None,
    report=None,
):
    """
    Organizes information about medications in UK Biobank data-field "20003" and
    their classes in the Anatomical Therapeutic Chemical (ATC) classification
    system (https://www.who.int/tools/atc-ddd-toolkit/atc-classification).

    arguments:
        table_wu (object): Pandas data frame of medication codes in UK Biobank
            data-field "20003" and their ATC classes from Wu et al, 2019
        table_wu (object): Pandas data frame of medication codes in UK Biobank
            data-field "20003" and their ATC classes from Appleby et al, 2019
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of medication codes in UK Biobank data-field
            "20003" and their ATC classes

    """

    # Copy information in table.
    table_wu = table_wu.copy(deep=True)
    table_appleby = table_appleby.copy(deep=True)
    # Select relevant columns.
    columns_relevant = [
        "code_ukbiobank",
        "medication_name",
        "medication_category",
        "code_atc",
        "group_atc",
    ]
    table_wu = table_wu.loc[
        :, table_wu.columns.isin(columns_relevant)
    ]
    table_appleby = table_appleby.loc[
        :, table_appleby.columns.isin(columns_relevant)
    ]
    # Organize table indices.
    table_wu.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_wu["code_ukbiobank"].astype("string")
    table_wu.dropna(
        axis="index",
        how="any",
        subset=["code_ukbiobank"],
        inplace=True,
    )
    table_wu.drop_duplicates(
        subset=["code_ukbiobank",],
        keep="first",
        inplace=True,
    )
    table_wu.set_index(
        "code_ukbiobank",
        drop=True,
        inplace=True,
    )
    table_appleby.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_appleby["code_ukbiobank"].astype("string")
    table_appleby.dropna(
        axis="index",
        how="any",
        subset=["code_ukbiobank"],
        inplace=True,
    )
    table_appleby.drop_duplicates(
        subset=["code_ukbiobank",],
        keep="first",
        inplace=True,
    )
    table_appleby.set_index(
        "code_ukbiobank",
        drop=True,
        inplace=True,
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = pandas.merge(
        table_wu, # left table
        table_appleby, # right table
        left_on="code_ukbiobank",
        right_on="code_ukbiobank",
        left_index=False,
        right_index=False,
        how="outer", # keep keys from union of both tables
        suffixes=("_wu", "_appleby"),
    )
    table_merge.reset_index(
        level=None,
        inplace=True,
        drop=False, # preserve index as a column
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_medication_code_atc_class_table()")
        utility.print_terminal_partition(level=3)
        print(table_merge)
        print(table_merge.columns.to_list())
    # Return information.
    return table_merge


def organize_atc_class_medication_codes(
    atc_classes=None,
    table_medication_codes_classes=None,
    report=None,
):
    """
    Organizes information about medications in UK Biobank data-field "20003" and
    their classes in the Anatomical Therapeutic Chemical (ATC) classification
    system (https://www.who.int/tools/atc-ddd-toolkit/atc-classification).

    arguments:
        atc_classes (list<str>): codes of relevant ATC classes
        table_medication_codes_classes (object): Pandas data frame of medication
            codes in UK Biobank data-field "20003" and their ATC classes
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<str>>): medication codes in relevant ATC classes

    """

    # Copy information in table.
    table_codes = table_medication_codes_classes.copy(deep=True)
    # Collect information.
    pail = dict()
    # Iterate on relevant ATC classes.
    for atc_class in atc_classes:
        # Filter tables by relevant ATC classes.
        table_class_wu = table_codes.loc[
            (
                ~pandas.isna(table_codes["group_atc_wu"]) &
                (table_codes["group_atc_wu"] == atc_class)
            ), :
        ]
        table_class_appleby = table_codes.loc[
            (
                ~pandas.isna(table_codes["group_atc_appleby"]) &
                (table_codes["group_atc_appleby"] == atc_class)
            ), :
        ]
        codes_class_wu = copy.deepcopy(
            table_class_wu["code_ukbiobank"].to_list()
        )
        codes_class_wu = list(map(str, codes_class_wu)) # string
        codes_class_wu = list(set(codes_class_wu)) # unique
        codes_class_appleby = copy.deepcopy(
            table_class_appleby["code_ukbiobank"].to_list()
        )
        codes_class_appleby = list(map(str, codes_class_appleby)) # string
        codes_class_appleby = list(set(codes_class_appleby)) # unique
        codes_class = copy.deepcopy(codes_class_wu)
        codes_class.extend(codes_class_appleby)
        codes_class = list(map(str, codes_class)) # string
        codes_class = list(set(codes_class)) # unique
        # Collect information.
        pail[atc_class] = dict()
        pail[atc_class]["codes_class_wu"] = codes_class_wu
        pail[atc_class]["codes_class_appleby"] = codes_class_appleby
        pail[atc_class]["codes_class"] = codes_class
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_atc_class_medication_codes()")
        utility.print_terminal_partition(level=3)
        for atc_class in pail.keys():
            count_codes_class_wu = len(pail[atc_class]["codes_class_wu"])
            count_codes_class_appleby = len(
                pail[atc_class]["codes_class_appleby"]
            )
            count_codes_class = len(pail[atc_class]["codes_class"])
            # Report.
            utility.print_terminal_partition(level=5)
            print("ATC class: " + str(atc_class))
            print("medications from Wu, 2019: " + str(count_codes_class_wu))
            print(
                "medications from Appleby, 2019: " +
                str(count_codes_class_appleby)
            )
            print("medications from consensus: " + str(count_codes_class))
            pass
        pass
    # Return information.
    return pail


def interpret_field_array_codes_match_any_relevant_codes(
    field_array_values_text=None,
    delimiter=None,
    relevant_query_codes=None,
):
    """
    Determines whether any actual codes from a UK Biobank data-field match
    relevant query codes.

    Match uses "any" logic.

    Implementation does not introduce missing values.

    arguments:
        field_array_values_text (str): textual representation of a list of codes
            from a UK Biobank array data-field
        delimiter (str): string delimiter between values in the textual
            representation of a list of codes
        relevant_query_codes (list<str>): relevant codes for which to query any
            matches

    raises:

    returns:
        (float): binary logical representation of whether any of record's codes
            match relevant query codes

    """

    # Determine whether any actual medication codes match the comparisons.
    values_actual = utility.parse_text_list_values(
        collection=field_array_values_text,
        delimiter=delimiter,
    )
    match = utility.determine_any_actual_values_match_comparisons(
        values_actual=values_actual,
        values_comparison=relevant_query_codes,
    )
    # Determine logical binary representation of whether there was any match.
    if (match):
        interpretation = 1
    else:
        interpretation = 0
        pass
    # Return information.
    return interpretation


def report_medication_class_groups_by_sex(
    medication_class_name=None,
    medication_use_indicator=None,
    table=None,
):
    """
    Reports counts of persons with stratification by sex who used medications in
    a specific class group.

    arguments:
        medication_class_name (str): name of medication class
        medication_use_indicator (str): name of column in table for binary
            logical indicator variable
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Sex total.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_male = table.loc[
        (
            (table["sex_text"] == "male")
        ), :
    ]
    # Sex medication use.
    table_female_medication_true = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[medication_use_indicator] == 1)
        ), :
    ]
    table_female_medication_false = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[medication_use_indicator] == 0)
        ), :
    ]
    table_male_medication_true = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[medication_use_indicator] == 1)
        ), :
    ]
    table_male_medication_false = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[medication_use_indicator] == 0)
        ), :
    ]
    # Counts.
    count_female = table_female.shape[0]
    count_male = table_male.shape[0]
    count_female_true = table_female_medication_true.shape[0]
    count_female_false = table_female_medication_false.shape[0]
    count_male_true = table_male_medication_true.shape[0]
    count_male_false = table_male_medication_false.shape[0]
    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_medication_class_groups_by_sex()")
    utility.print_terminal_partition(level=3)
    print("Name of medication class: " + str(medication_class_name))
    utility.print_terminal_partition(level=4)
    print("counts for: FEMALES")
    utility.print_terminal_partition(level=5)
    print("count of total: " + str(count_female))
    print("count who use medication class: " + str(count_female_true))
    print("count who do not use medication class: " + str(count_female_false))
    utility.print_terminal_partition(level=4)
    print("counts for: MALES")
    utility.print_terminal_partition(level=5)
    print("count of total: " + str(count_male))
    print("count who use medication class: " + str(count_male_true))
    print("count who do not use medication class: " + str(count_male_false))
    pass


def organize_hormonal_medications(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes interpretation of medications that alter hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Read reference table of interpretations of variable codes.
    source = read_source_medication_codes_classes(
        path_dock=path_dock,
    )
    table_medication_codes_classes = organize_medication_code_atc_class_table(
        table_wu=source["table_wu"],
        table_appleby=source["table_appleby"],
        report=report,
    )
    medication_classes = organize_atc_class_medication_codes(
        atc_classes=["G03", "A11CC",],
        table_medication_codes_classes=table_medication_codes_classes,
        report=report,
    )

    # Copy information in table.
    table = table.copy(deep=True)
    # Determine whether persons used medications that alter either Sex Hormones
    # or Vitamin D.
    table["medication_sex_hormone"] = table.apply(
        lambda row:
            interpret_field_array_codes_match_any_relevant_codes(
                field_array_values_text=row["20003_array"],
                delimiter=";",
                relevant_query_codes=medication_classes["G03"]["codes_class"],
            ),
        axis="columns", # apply function to each row
    )
    table["medication_vitamin_d"] = table.apply(
        lambda row:
            interpret_field_array_codes_match_any_relevant_codes(
                field_array_values_text=row["20003_array"],
                delimiter=";",
                relevant_query_codes=medication_classes["A11CC"]["codes_class"],
            ),
        axis="columns", # apply function to each row
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_hormonal_medications()")
        utility.print_terminal_partition(level=3)
        report_medication_class_groups_by_sex(
            medication_class_name="ATC class G03",
            medication_use_indicator="medication_sex_hormone",
            table=table,
        )
        report_medication_class_groups_by_sex(
            medication_class_name="ATC class A11CC",
            medication_use_indicator="medication_vitamin_d",
            table=table,
        )
        pass
    # Return information.
    return table


# Female Oral Contraception or Hormonal Therapy
# Review:
# 14 January 2022: TCW verified the documentations and interpretations of
# data-fields for recent use of Oral Contraception or Hormonal Therapy.


def interpret_recent_female_hormone_therapies(
    age=None,
    recent_range=None,
    ever_used=None,
    age_started=None,
    age_stopped=None,
):
    """
    Inteprets whether a female person used recently either oral contraception or
    hormone-replacement therapy.

    Relevant UK Biobank data fields and their codings:
    oral contraception
    2784 (100349)
    2794 (100291)
    2804 (100595)
    hormone-replacement therapy
    2814 (100349)
    3536 (100291)
    3546 (100598)

    Accommodate inexact float values.

    arguments:
        age (int): age of person in years
        recent_range (int): years within current age to consider recent
        ever_used (float): UK Biobank field 2784 or 2814, ever taken either oral
            contraception or hormone-replacement therapy
        age_started (float): UK Biobank field 2794 or 3536, age started either
            oral contraception or hormone-replacement therapy
        age_stopped (float): UK Biobank field 2804 or 3546, age stopped either oral
            contraception or hormone-replacement therapy

    raises:

    returns:
        (float): interpretation value (1: true, 0: false, NAN: missing, null)

    """

    if (
        (not pandas.isna(ever_used)) and
        (-3.5 <= ever_used and ever_used < 1.5)
    ):
        # Variable "ever_used" has a valid value.
        # Determine whether person has ever used therapy.
        if (0.5 <= ever_used and ever_used < 1.5):
            # 1: "Yes"
            # Person has used therapy.
            # Determine whether person used therapy recently.
            if (
                (not pandas.isna(age_stopped)) and
                (-11.5 <= age_stopped and age_stopped <= age)
            ):
                # Variable "age_stopped" has a valid value.
                if (-11.5 <= age_stopped and age_stopped < -10.5):
                    # -11: "Still taking the pill"
                    # -11: "Still taking HRT"
                    # Person used oral contraception currently.
                    value = 1
                elif (
                    (age_stopped >= (age - recent_range)) and
                    (age_stopped <= age)
                ):
                    # Person used therapy within recent range of current age.
                    value = 1
                elif ((age_stopped < (age - recent_range))):
                    # Person ceased use of therapy longer than recently.
                    value = 0
                elif (-1.5 <= age_stopped and age_stopped < -0.5):
                    # -1: "Do not know"
                    # Unable to determine how recently person used therapy.
                    value = float("nan")
                elif (-3.5 <= age_stopped and age_stopped < -2.5):
                    # -3: "Prefer not to answer"
                    # Unable to determine how recently person used therapy.
                    value = float("nan")
                else:
                    # uninterpretable
                    # Unable to determine how recently person used therapy.
                    value = float("nan")
            else:
                # Unable to determine how recently person used therapy.
                # null
                value = float("nan")
        elif (-0.5 <= ever_used and ever_used < 0.5):
            # 0: "No"
            # Person has not used therapy.
            value = 0
        elif (-1.5 <= ever_used and ever_used < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= ever_used and ever_used < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def interpret_recent_oral_contraception(
    age=None,
    recent_range=None,
    field_2784=None,
    field_2794=None,
    field_2804=None,
):
    """
    Inteprets recent use of oral contraception.

    Recent use of oral contraception will be a covariate in analyses.
    To preserve samples in analyses, do not perpetuate null values.
    Interpret null values as false.
    Hence the greatest confidence will be in true values, and these will drive
    much of the signal in analyses.

    Data-Field "2784": "Ever taken oral contraceptive pill"
    UK Biobank data coding "100349" for variable field "2784".
    1: "Yes"
    0: "No"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Data-Field "2794": "Age started oral contraceptive pill"
    UK Biobank data coding "100291" for variable field "2794".
    [5 to person's age]: "age in years"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Data-Field "2804": "Age when last used oral contraceptive pill"
    UK Biobank data coding "100595" for variable field "2804".
    [5 to person's age]: "age in years"
    -1: "Do not know"
    -3: "Prefer not to answer"
    -11: "Still taking the pill"

    Accommodate inexact float values.

    arguments:
        age (int): age of person in years
        recent_range (int): years within current age to consider recent
        field_2784 (float): UK Biobank field 2784, ever taken oral contraception
        field_2794 (float): UK Biobank field 2794, age started oral
            contraception
        field_2804 (float): UK Biobank field 2804, age stopped oral
            contraception

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    value = interpret_recent_female_hormone_therapies(
        age=age,
        recent_range=recent_range,
        ever_used=field_2784,
        age_started=field_2794,
        age_stopped=field_2804,
    )
    # Return.
    return value


def interpret_recent_hormone_replacement_therapy(
    age=None,
    recent_range=None,
    field_2814=None,
    field_3536=None,
    field_3546=None,
):
    """
    Inteprets recent use of hormone replacement therapy.

    Recent use of hormone replacement therapy will be a covariate in analyses.
    To preserve samples in analyses, do not perpetuate null values.
    Interpret null values as false.
    Hence the greatest confidence will be in true values, and these will drive
    much of the signal in analyses.

    Data-Field "2814": "Ever used hormone-replacement therapy (HRT)"
    UK Biobank data coding "100349" for variable field "2814".
    1: "Yes"
    0: "No"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Data-Field "3536": "Age started hormone-replacement therapy (HRT)"
    UK Biobank data coding "100291" for variable field "3536".
    [16 to person's age]: "age in years"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Data-Field "3546": "Age last used hormone-replacement therapy (HRT)"
    UK Biobank data coding "100598" for variable field "3546".
    [20 to person's age]: "age in years"
    -1: "Do not know"
    -3: "Prefer not to answer"
    -11: "Still taking HRT"

    Accommodate inexact float values.

    arguments:
        age (int): age of person in years
        recent_range (int): years within current age to consider recent
        field_2814 (float): UK Biobank field 2814, ever taken hormone
            replacement therapy (HRT)
        field_3536 (float): UK Biobank field 3536, age started hormone
            replacement therapy (HRT)
        field_3546 (float): UK Biobank field 3546, age stopped hormone
            replacement therapy (HRT)

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    value = interpret_recent_female_hormone_therapies(
        age=age,
        recent_range=recent_range,
        ever_used=field_2814,
        age_started=field_3536,
        age_stopped=field_3546,
    )
    # Return.
    return value


def determine_female_oral_contraception(
    sex_text=None,
    age=None,
    recent_range=None,
    field_2784=None,
    field_2794=None,
    field_2804=None,
    field_6153_values_text=None,
):
    """
    Determine whether person used Oral Contraception recently or currently.

    Data-Field "6153": "Medication for cholesterol, blood pressure, diabetes, or
    take exogenous hormones"
    UK Biobank data coding "100626" for data-field "6153".
    1: "Cholesterol lowering medication"
    2: "Blood pressure medication"
    3: "Insulin"
    4: "Hormone replacement therapy"
    5: "Oral contraceptive pill or minipill"
    -7: "None of the above"
    -1: "Do not know"
    -3: "Prefer not to answer"

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        recent_range (int): years within current age to consider recent
        field_2784 (float): UK Biobank field 2784, ever taken oral contraception
        field_2794 (float): UK Biobank field 2794, age started oral
            contraception
        field_2804 (float): UK Biobank field 2804, age stopped oral
            contraception
        field_6153_values_text (str): textual representation of a list of
            medication codes from UK Biobank data-field "6153"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret oral contraception.
    therapy = interpret_recent_oral_contraception(
        age=age,
        recent_range=recent_range, # recent range in years
        field_2784=field_2784,
        field_2794=field_2794,
        field_2804=field_2804,
    )
    # Interpret array codes.
    match = interpret_field_array_codes_match_any_relevant_codes(
        field_array_values_text=field_6153_values_text,
        delimiter=";",
        relevant_query_codes=["5",], # data-field "6153" coding "100626"
    )
    # Comparison.
    if (sex_text == "female"):
        value = utility.determine_logical_or_combination_binary_missing(
            first=therapy,
            second=match,
            single_false_sufficient=True,
        )
    else:
        # This specific variable is undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_hormone_replacement_therapy(
    sex_text=None,
    age=None,
    recent_range=None,
    field_2814=None,
    field_3536=None,
    field_3546=None,
    field_6153_values_text=None,
):
    """
    Determine whether person used Hormone Therapy recently or currently.

    Data-Field "6153": "Medication for cholesterol, blood pressure, diabetes, or
    take exogenous hormones"
    UK Biobank data coding "100626" for data-field "6153".
    1: "Cholesterol lowering medication"
    2: "Blood pressure medication"
    3: "Insulin"
    4: "Hormone replacement therapy"
    5: "Oral contraceptive pill or minipill"
    -7: "None of the above"
    -1: "Do not know"
    -3: "Prefer not to answer"

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        recent_range (int): years within current age to consider recent
        field_2814 (float): UK Biobank field 2814, ever taken hormone
            replacement therapy (HRT)
        field_3536 (float): UK Biobank field 3536, age started hormone
            replacement therapy (HRT)
        field_3546 (float): UK Biobank field 3546, age stopped hormone
            replacement therapy (HRT)
        field_6153_values_text (str): textual representation of a list of
            medication codes from UK Biobank data-field "6153"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret hormone replacement therapy.
    therapy = interpret_recent_hormone_replacement_therapy(
        age=age,
        recent_range=recent_range,
        field_2814=field_2814,
        field_3536=field_3536,
        field_3546=field_3546,
    )
    # Interpret array codes.
    match = interpret_field_array_codes_match_any_relevant_codes(
        field_array_values_text=field_6153_values_text,
        delimiter=";",
        relevant_query_codes=["4",], # data-field "6153" coding "100626"
    )
    # Comparison.
    if (sex_text == "female"):
        value = utility.determine_logical_or_combination_binary_missing(
            first=therapy,
            second=match,
            single_false_sufficient=True,
        )
    else:
        # This specific variable is undefined for males.
        value = float("nan")
    # Return information.
    return value


def organize_female_hormonal_therapy(
    table=None,
    report=None,
):
    """
    Organizes interpretation of female-specific use of Oral Contraception or
    Hormonal Therapy.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "2784-0.0", "2794-0.0", "2804-0.0",
        "2814-0.0", "3536-0.0", "3546-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )


    # Determine whether person was using oral contraception recently.
    table["oral_contraception"] = table.apply(
        lambda row:
            determine_female_oral_contraception(
                sex_text=row["sex_text"],
                age=row["age"],
                recent_range=1, # years within current age to consider recent
                field_2784=row["2784-0.0"],
                field_2794=row["2794-0.0"],
                field_2804=row["2804-0.0"],
                field_6153_values_text=row["6153_array"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person was using hormone replacement therapy recently.
    table["hormone_therapy"] = table.apply(
        lambda row:
            determine_female_hormone_replacement_therapy(
                sex_text=row["sex_text"],
                age=row["age"],
                recent_range=1, # years within current age to consider recent
                field_2814=row["2814-0.0"],
                field_3536=row["3536-0.0"],
                field_3546=row["3546-0.0"],
                field_6153_values_text=row["6153_array"],
            ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_female_hormonal_therapy()")
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


# Combination of Hormonal Alteration


def determine_sex_combination_alteration_sex_hormones(
    sex_text=None,
    female_oral_contraception=None,
    female_hormone_therapy=None,
    medication_sex_hormone=None,
):
    """
    Organizes interpretation of medications that alter hormones.

    arguments:
        sex_text (str): textual representation of sex, 'female' or 'male'
        female_oral_contraception (float): binary logical representation of
            whether female person used oral contraception recently
        female_hormone_therapy (float): binary logical representation of whether
            female person used hormone replacement therapy recently
        medication_sex_hormone (float): binary logical representation of whether
            person used any medications that influence sex hormones

    raises:

    returns:
        (float): binary logical representation of whether person used any
            medications that influence sex hormones, including oral
            contraception and other hormonal therapy in females

    """

    if (sex_text == "female"):
        combination_1 = utility.determine_logical_or_combination_binary_missing(
            first=female_oral_contraception,
            second=female_hormone_therapy,
            single_false_sufficient=True,
        )
        combination_2 = utility.determine_logical_or_combination_binary_missing(
            first=combination_1,
            second=medication_sex_hormone,
            single_false_sufficient=True,
        )
        interpretation = combination_2
    elif (sex_text == "male"):
        interpretation = medication_sex_hormone
    else:
        # Interpretation is null without definition of sex.
        interpretation = float("nan")
        pass
    # Return information.
    return interpretation


def report_combination_alteration_sex_hormones_by_sex(
    female_oral_contraception=None,
    female_hormone_therapy=None,
    medication_sex_hormone=None,
    alteration_sex_hormone=None,
    table=None,
):
    """
    Reports counts of persons with stratification by sex who used medications in
    a specific class group.

    arguments:
        female_oral_contraception (str): name of column for binary logical
            representation of whether female person used oral contraception
            recently
        female_hormone_therapy (str): name of column for binary logical
            representation of whether female person used hormone replacement
            therapy recently
        medication_sex_hormone (str): name of column for binary logical
            representation of whether person used any medications that influence
            sex hormones
        alteration_sex_hormone (str): name of column for binary logical
            representation of whether person used any medications that influence
            sex hormones, including oral contraception and other hormonal
            therapy in females
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Sex total.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_male = table.loc[
        (
            (table["sex_text"] == "male")
        ), :
    ]
    # Females.
    table_female_contraception = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[female_oral_contraception] == 1)
        ), :
    ]
    table_female_therapy = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[female_hormone_therapy] == 1)
        ), :
    ]
    table_female_medication = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[medication_sex_hormone] == 1)
        ), :
    ]
    table_female_alteration = table.loc[
        (
            (table["sex_text"] == "female") &
            (table[alteration_sex_hormone] == 1)
        ), :
    ]
    # Males.
    table_male_medication = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[medication_sex_hormone] == 1)
        ), :
    ]
    table_male_alteration = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[alteration_sex_hormone] == 1)
        ), :
    ]
    # Counts.
    count_female = table_female.shape[0]
    count_male = table_male.shape[0]
    count_female_contraception = table_female_contraception.shape[0]
    count_female_therapy = table_female_therapy.shape[0]
    count_female_medication = table_female_medication.shape[0]
    count_female_alteration = table_female_alteration.shape[0]
    count_male_medication = table_male_medication.shape[0]
    count_male_alteration = table_male_alteration.shape[0]
    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_combination_alteration_sex_hormones_by_sex()")
    utility.print_terminal_partition(level=3)
    utility.print_terminal_partition(level=4)
    print("counts for: FEMALES")
    utility.print_terminal_partition(level=5)
    print("count of total: " + str(count_female))
    print("count oral contraception: " + str(count_female_contraception))
    print("count hormone therapy: " + str(count_female_therapy))
    print("count medication: " + str(count_female_medication))
    print("count alteration: " + str(count_female_alteration))
    utility.print_terminal_partition(level=4)
    print("counts for: MALES")
    utility.print_terminal_partition(level=5)
    print("count of total: " + str(count_male))
    print("count medication: " + str(count_male_medication))
    print("count alteration: " + str(count_male_alteration))
    pass


def organize_combination_alteration_sex_hormones(
    table=None,
    female_oral_contraception=None,
    female_hormone_therapy=None,
    medication_sex_hormone=None,
    report=None,
):
    """
    Organizes interpretation of medications that alter hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        female_oral_contraception (str): name of column for binary logical
            representation of whether female person used oral contraception
            recently
        female_hormone_therapy (str): name of column for binary logical
            representation of whether female person used hormone replacement
            therapy recently
        medication_sex_hormone (str): name of column for binary logical
            representation of whether person used any medications that influence
            sex hormones
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Determine whether persons used any medications or other therapeutic
    # exposures that alter Sex Hormones.
    table["alteration_sex_hormone"] = table.apply(
        lambda row:
            determine_sex_combination_alteration_sex_hormones(
                sex_text=row["sex_text"],
                female_oral_contraception=row[female_oral_contraception],
                female_hormone_therapy=row[female_hormone_therapy],
                medication_sex_hormone=row[medication_sex_hormone],
            ),
        axis="columns", # apply function to each row
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_hormonal_medications()")
        utility.print_terminal_partition(level=3)
        report_combination_alteration_sex_hormones_by_sex(
            female_oral_contraception=female_oral_contraception,
            female_hormone_therapy=female_hormone_therapy,
            medication_sex_hormone=medication_sex_hormone,
            alteration_sex_hormone="alteration_sex_hormone",
            table=table,
        )
        pass
    # Return information.
    return table


# Report



# review: TCW on 20 January January 2022
def report_hormone_deficiency_by_sex_menopause_age(
    column_hormone=None,
    threshold_deficiency=None,
    table=None,
):
    """
    Reports counts and percentages of persons who were deficient in a hormone
    with stratification by sex, female menopause status, and age.

    arguments:
        column_hormone (str): name of column for hormone measurement
        threshold_deficiency (float): low threshold, concentrations below which
            qualify as deficiency
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_female_menstruation = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    # Counts.

    count_female = table_female.shape[0]

    count_irregularity = table_female_irregularity.shape[0]
    count_regularity = table_female_regularity.shape[0]

    count_valid_duration = table_female_valid_duration.shape[0]
    count_valid_current = table_female_valid_current.shape[0]
    count_valid_duration_or_current = (
        table_female_valid_duration_or_current.shape[0]
    )

    count_regular_range = table_female_regular_range.shape[0]
    count_duration_below = table_female_duration_below_range.shape[0]
    count_duration_above = table_female_duration_above_range.shape[0]
    count_current_above = table_female_current_days_above_range.shape[0]

    count_discrepancy_duration_current = (
        table_discrepancy_duration_current.shape[0]
    )
    count_discrepancy_regular_current = (
        table_discrepancy_regular_current.shape[0]
    )

    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_female_menstruation_regularity_duration_range()")
    utility.print_terminal_partition(level=3)
    print("... all counts specific to female persons ...")
    utility.print_terminal_partition(level=4)
    print("total: " + str(count_female))

    print("'irregular' menstrual cycle: " + str(count_irregularity))
    print("'regular' menstrual cycle: " + str(count_regularity))

    print(
        "non-missing regular duration of menstrual cycle: " +
        str(count_valid_duration)
    )
    print(
        "non-missing days of current menstrual cycle: " +
        str(count_valid_current)
    )
    print(
        "non-missing regular duration or days of current menstrual cycle: " +
        str(count_valid_duration_or_current)
    )

    print(
        "regular menstrual cycle of duration within threshold range: " +
        str(count_regular_range)
    )
    print(
        "regular menstrual cycle duration shorter than " +
        str(threshold_duration_low) +
        " days: " +
        str(count_duration_below)
    )
    print(
        "regular menstrual cycle duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_duration_above)
    )
    print(
        "current menstrual cycle at duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_current_above)
    )
    print(
        "regular menstrual cycle duration shorter than " +
        str(threshold_duration_high) +
        " days, but current menstrual cycle at duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_discrepancy_duration_current)
    )
    print(
        "regular menstrual cycle of duration within " +
        "threshold range, but current menstrual cycle at duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_discrepancy_regular_current)
    )
    pass




# Main driver


def organize_sex_hormone_variables(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy information in table.
    table = table.copy(deep=True)

    ##########
    # Organize raw hormone variables.

    # Define reportability and missingness variables for each hormone.
    hormone_fields = (
        define_biochemistry_fields_measurement_reportability_missingness()
    )
    # Convert variable types to float.
    columns_hormones = list()
    for record in hormone_fields:
        columns_hormones.append(record["measurement"])
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_hormones,
        table=table,
    )
    # Convert concentrations to units of moles per liter (mol/L) with adjustment
    # by specific factors for appropriate scale in analyses (and floats).
    factors_concentration = dict()
    factors_concentration["albumin"] = 1E6 # 1 umol / L
    factors_concentration["steroid_globulin"] = 1E9 # 1 nmol / L
    factors_concentration["oestradiol"] = 1E12 # 1 pmol / L
    factors_concentration["oestradiol_free"] = 1E12 # 1 pmol / L
    factors_concentration["oestradiol_bioavailable"] = 1E12 # 1 pmol / L
    factors_concentration["testosterone"] = 1E12 # 1 pmol / L
    factors_concentration["testosterone_free"] = 1E12 # 1 pmol / L
    factors_concentration["testosterone_bioavailable"] = 1E12 # 1 pmol / L
    factors_concentration["vitamin_d"] = 1E9 # 1 nmol / L
    factors_concentration["cholesterol"] = 1E3 # 1 mmol / L
    table = convert_hormone_concentration_units_moles_per_liter(
        table=table,
        factors_concentration=factors_concentration,
    )

    ##########
    # Determine binary representations of whether hormone measurements were
    # missing due to limit of detection.
    # Impute hormone measurements that were missing due to limit of detection.
    table = organize_hormones_missingness_detection_imputation(
        hormone_fields=hormone_fields,
        table=table,
        report=report,
    )

    ##########
    # Define cohort-specific ordinal representations of hormones with
    # consideration of cohort-specific medians and reportabilities due to
    # detection limits.
    if False:
        table = organize_cohorts_specific_hormones_ordinal_representations(
            table=table,
            report=report,
        )

    ##########
    # Calculate estimates of bioavailable and free hormones.
    table = organize_calculation_estimate_bioavailable_free_hormones(
        factors_concentration=factors_concentration,
        table=table,
        report=report,
    )

    ##########
    # Determine use of medications that alter hormones.
    # These medications are relevant to both females and males.
    table = organize_hormonal_medications(
        table=table,
        path_dock=path_dock,
        report=report,
    )

    ##########
    # Determine female-specific use of Oral Contraception or Hormonal Therapy.
    table = organize_female_hormonal_therapy(
        table=table,
        report=report,
    )

    ##########
    # Determine combination of exposures that alter sex hormones in females and
    # males.
    table = organize_combination_alteration_sex_hormones(
        table=table,
        female_oral_contraception="oral_contraception",
        female_hormone_therapy="hormone_therapy",
        medication_sex_hormone="medication_sex_hormone",
        report=report,
    )


    ##########
    # Transform variables' values to normalize distributions.
    columns_hormones_normalization = [
        "oestradiol", "oestradiol_imputation",
        "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
        "oestradiol_free", "oestradiol_free_imputation",
        "testosterone", "testosterone_imputation",
        "testosterone_bioavailable", "testosterone_bioavailable_imputation",
        "testosterone_free", "testosterone_free_imputation",
        "vitamin_d", "vitamin_d_imputation",
        "steroid_globulin", "steroid_globulin_imputation",
        "albumin", "albumin_imputation",
        "cholesterol", "cholesterol_imputation",
    ]
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=columns_hormones_normalization,
        table=table,
    )
    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "30600-0.0", "30830-0.0", "30850-0.0", "30800-0.0", "30890-0.0",
            "2784-0.0", "2794-0.0", "2804-0.0", # Oral Contraception
            "2814-0.0", "3536-0.0", "3546-0.0", # Hormonal Therapy
        ],
        axis="columns",
        inplace=True
    )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "vitamin_d",
        "vitamin_d_detection",
        "vitamin_d_imputation",
        "albumin",
        "albumin_detection",
        "albumin_imputation",
        "steroid_globulin",
        "steroid_globulin_detection",
        "steroid_globulin_imputation",
        "oestradiol",
        "oestradiol_detection",
        "oestradiol_imputation",
        "oestradiol_bioavailable",
        "oestradiol_bioavailable_imputation",
        "oestradiol_free",
        "oestradiol_free_imputation",
        "testosterone",
        "testosterone_detection",
        "testosterone_imputation",
        "testosterone_bioavailable",
        "testosterone_bioavailable_imputation",
        "testosterone_free",
        "testosterone_free_imputation",
        "cholesterol",
        "cholesterol_imputation",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: organize_sex_hormone_variables()")
        utility.print_terminal_partition(level=3)
        print("Translation of columns for hormones: ")
        print(table_report)
        utility.print_terminal_partition(level=3)
        # Variable types.
        utility.print_terminal_partition(level=2)
        print("After type conversion")
        print(table_report.dtypes)
        utility.print_terminal_partition(level=3)

        table_correlation = table.copy(deep=True)

        utility.print_terminal_partition(level=2)
        print("correlations in MALES!")
        table_male = table_correlation.loc[
            (table_correlation["sex_text"] == "male"), :
        ]
        utility.print_terminal_partition(level=3)
        utility.calculate_table_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_free",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_bioavailable",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="testosterone_free",
            column_two="testosterone_bioavailable",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_free",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_bioavailable",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="oestradiol_free",
            column_two="oestradiol_bioavailable",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="vitamin_d",
            column_two="oestradiol",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="vitamin_d",
            column_two="testosterone",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="cholesterol",
            column_two="vitamin_d",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="cholesterol",
            column_two="oestradiol",
            table=table_male,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="cholesterol",
            column_two="testosterone",
            table=table_male,
            report=True,
        )


        utility.print_terminal_partition(level=2)
        print("correlations in FEMALES!")
        table_female = table_correlation.loc[
            (table_correlation["sex_text"] == "female"), :
        ]
        utility.print_terminal_partition(level=3)
        utility.calculate_table_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_free",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_bioavailable",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="testosterone_free",
            column_two="testosterone_bioavailable",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_free",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_bioavailable",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="oestradiol_free",
            column_two="oestradiol_bioavailable",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="vitamin_d",
            column_two="oestradiol",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="vitamin_d",
            column_two="testosterone",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="cholesterol",
            column_two="vitamin_d",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="cholesterol",
            column_two="oestradiol",
            table=table_female,
            report=True,
        )
        utility.calculate_table_column_pair_correlations(
            column_one="cholesterol",
            column_two="testosterone",
            table=table_female,
            report=True,
        )

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Menstruation, pregnancy, menopause, contraception, hormone therapy
# Review:

# TODO: I need to review new definitions of menopause variables...

# 31 March 2021: TCW verified UK Biobank fields and their codings.
# 11 March 2021: TCW verified UK Biobank fields and their codings.

# TODO: review all "interpretation" functions below...


# review: TCW on 17 January 2022
def interpret_menstruation_days(
    field_3700=None,
):
    """
    Intepret UK Biobank's coding for data-field 3700.

    Data-Field "3700": "Time since last menstrual period"
    UK Biobank data-coding "100291" for data-field "3700".
    [0 - 365]: "days"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Note: "How many days since your last menstrual period?"
    Note: "Please count from the first day of your last menstrual period."

    Accommodate inexact float values.

    arguments:
        field_3700 (float): UK Biobank field 3700, count of days since previous
            menstruation (menstrual period)

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3700)) and
        (-3.5 <= field_3700 and field_3700 < 400)
    ):
        # The variable has a valid value.
        if (0 <= field_3700 and field_3700 < 400):
            # The variable has a valid value for days.
            value = float(field_3700)
        elif (-1.5 <= field_3700 and field_3700 < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= field_3700 and field_3700 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_menstruation_cycle_duration(
    field_3710=None,
):
    """
    Intepret UK Biobank's coding for data-field 3710.

    Data-Field "3710": "Length of menstrual cycle"
    UK Biobank data-coding "100582" for data-field "3710".
    [7 - 365]: "days"
    -6: "Irregular cycle"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Note: "How many days is your usual menstrual cycle? (The number of days
    between each menstrual period)"

    Accommodate inexact float values.

    arguments:
        field_3710 (float): UK Biobank field 3710, count of days in usual
            menstrual cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3710)) and
        (-6.5 <= field_3710 and field_3710 < 400)
    ):
        # The variable has a valid value.
        if (0 <= field_3710 and field_3710 < 400):
            # The variable has a valid value for days.
            value = float(field_3710)
        elif (-1.5 <= field_3710 and field_3710 < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= field_3710 and field_3710 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        elif (-6.5 <= field_3710 and field_3710 < -5.5):
            # -6: "Irregular cycle"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_menstruation_cycle_irregularity(
    field_3710=None,
):
    """
    Intepret UK Biobank's coding for data-field 3710.

    Only interpret whether the variable indicates that the person had an
    irregular menstrual cycle.

    Data-Field "3710": "Length of menstrual cycle"
    UK Biobank data-coding "100582" for data-field "3710".
    [7 - 365]: "days"
    -6: "Irregular cycle"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Note: "How many days is your usual menstrual cycle? (The number of days
    between each menstrual period)"

    Accommodate inexact float values.

    arguments:
        field_3710 (float): UK Biobank field 3710, count of days in usual
            menstrual cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3710)) and
        (-6.5 <= field_3710 and field_3710 < 400)
    ):
        # The variable has a valid value.
        if (0 <= field_3710 and field_3710 < 400):
            # The variable has a valid value for days.
            value = 0
        elif (-1.5 <= field_3710 and field_3710 < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= field_3710 and field_3710 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        elif (-6.5 <= field_3710 and field_3710 < -5.5):
            # -6: "Irregular cycle"
            value = 1
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_menstruation_period_current(
    field_3720=None,
):
    """
    Intepret UK Biobank's coding for data-field 3720.

    Data-Field "3720": "Menstruating today"
    UK Biobank data-coding "100349" for data-field "3720".
    0: "No"
    1: "Yes"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_3720 (float): UK Biobank field 3720, whether person was
            experiencing menstruation period on day of assessment

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3720)) and
        (-3.5 <= field_3720 and field_3720 < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_3720 and field_3720 < 0.5):
            # 0: "No"
            value = 0
        elif (0.5 <= field_3720 and field_3720 < 1.5):
            # 1: "Yes"
            value = 1
        elif (-1.5 <= field_3720 and field_3720 < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= field_3720 and field_3720 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_menopause_hysterectomy(
    field_2724=None,
):
    """
    Intepret UK Biobank's coding for data-field 2724.

    Only interpret whether the variable indicates that the person had a
    hysterectomy. The only valid value is whether person indicated hysterectomy.

    Data-Field "2724": "Had menopause"
    UK Biobank data-coding "100579" for data-field "2724".
    0: "No"
    1: "Yes"
    2: "Not sure - had a hysterectomy"
    3: "Not sure - other reason"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_2724 (float): UK Biobank field 2724, whether person self reported
            menopause

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_2724)) and
        (-3.5 <= field_2724 and field_2724 < 3.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_2724 and field_2724 < 0.5):
            # 0: "no"
            value = float("nan")
        elif (0.5 <= field_2724 and field_2724 < 1.5):
            # 1: "yes"
            value = float("nan")
        elif (1.5 <= field_2724 and field_2724 < 2.5):
            # 2: "not sure - had a hysterectomy"
            value = 1
        elif (2.5 <= field_2724 and field_2724 < 3.5):
            # 3: "not sure - other reason"
            value = float("nan")
        elif (-3.5 <= field_2724 and field_2724 < -2.5):
            # -3: "prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_menopause_self_report(
    field_2724=None,
):
    """
    Intepret UK Biobank's coding for data-field 2724.

    Only interpret whether the variable indicates that the person self reported
    natural menopause, not involving hysterectomy.

    Data-Field "2724": "Had menopause"
    UK Biobank data-coding "100579" for data-field "2724".
    0: "No"
    1: "Yes"
    2: "Not sure - had a hysterectomy"
    3: "Not sure - other reason"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_2724 (float): UK Biobank field 2724, whether person self reported
            menopause

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_2724)) and
        (-3.5 <= field_2724 and field_2724 < 3.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_2724 and field_2724 < 0.5):
            # 0: "No"
            value = 0
        elif (0.5 <= field_2724 and field_2724 < 1.5):
            # 1: "Yes"
            value = 1
        elif (1.5 <= field_2724 and field_2724 < 2.5):
            # 2: "Not sure - had a hysterectomy"
            value = float("nan")
        elif (2.5 <= field_2724 and field_2724 < 3.5):
            # 3: "Not sure - other reason"
            value = float("nan")
        elif (-3.5 <= field_2724 and field_2724 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_hysterectomy(
    field_3591=None,
):
    """
    Intepret UK Biobank's coding for data-field 3591.

    Data-Field "3591": "Ever had hysterectomy (womb removed)"
    UK Biobank data-coding "100599" for data-field "3591".
    0: "No"
    1: "Yes"
    -5: "Not sure"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_3591 (float): UK Biobank field 3591, whether person has had a
            hysterectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3591)) and
        (-5.5 <= field_3591 and field_3591 < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_3591 and field_3591 < 0.5):
            # 0: "No"
            value = 0
        elif (0.5 <= field_3591 and field_3591 < 1.5):
            # 1: "Yes"
            value = 1
        elif (-5.5 <= field_3591 and field_3591 < -4.5):
            # -5: "Not sure"
            value = float("nan")
        elif (-3.5 <= field_3591 and field_3591 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
def interpret_oophorectomy(
    field_2834=None,
):
    """
    Intepret UK Biobank's coding for data-field 2834.

    Data-Field "2834": "Bilateral oophorectomy (both ovaries
    removed)"
    UK Biobank data-coding "100599" for data-field "2834".
    0: "No"
    1: "Yes"
    -5: "Not sure"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_2834 (float): UK Biobank field 2834, whether person has had an
            oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_2834)) and
        (-5.5 <= field_2834 and field_2834 < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_2834 and field_2834 < 0.5):
            # 0: "No"
            value = 0
        elif (0.5 <= field_2834 and field_2834 < 1.5):
            # 1: "Yes"
            value = 1
        elif (-5.5 <= field_2834 and field_2834 < -4.5):
            # -5: "Not sure"
            value = float("nan")
        elif (-3.5 <= field_2834 and field_2834 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 17 January 2022
# review TCW on 15 February 2022
def interpret_pregnancy(
    field_3140=None,
):
    """
    Intepret UK Biobank's coding for data-field 3140.

    This method interprets self-report "unsure" pregnancy as potential
    pregnancy. This method assumes that persons who were likely or trying
    to become pregnant might be more likely to answer "unsure".

    Data-Field "3140": "Pregnant"
    UK Biobank data-coding "100267" for data-field "3140".
    0: "No"
    1: "Yes"
    2: "Unsure"

    Accommodate inexact float values.

    arguments:
        field_3140 (float): UK Biobank field 3140, whether person was pregnant

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3140)) and
        (-0.5 <= field_3140 and field_3140 < 2.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_3140 and field_3140 < 0.5):
            # 0: "No"
            value = 0
        elif (0.5 <= field_3140 and field_3140 < 1.5):
            # 1: "Yes"
            value = 1
        elif (1.5 <= field_3140 and field_3140 < 2.5):
            # 2: "Unsure"
            # Interpret "unsure" pregnancy as potential pregnancy.
            #value = float("nan")
            value = 1
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 23 February 2022
def interpret_live_births(
    field_2734=None,
):
    """
    Intepret UK Biobank's coding for data-field "2734".

    Data-Field "2734": "Number of live births"

    UK Biobank data-coding "100584" for data-field "2734".
    0-25: "Count of events"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_2734 (float): UK Biobank field 2734, count of live births

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_2734)) and
        (-3.5 <= field_2734 and field_2734 < 25.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_2734 and field_2734 < 25.5):
            # 0-25: "count of events"
            value = float(field_2734)
        elif (-3.5 <= field_2734 and field_2734 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 23 February 2022
def interpret_pregnancy_loss(
    field_2774=None,
):
    """
    Intepret UK Biobank's coding for data-field "2774".

    Data-Field "2774": "Ever had stillbirth, spontaneous miscarriage or
    termination"

    UK Biobank data-coding "100349" for data-field "2774".
    0: "No"
    1: "Yes"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_2774 (float): UK Biobank field 2774, whether person experienced
            any type of loss of a pregnancy that did not end in live birth

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_2774)) and
        (-3.5 <= field_2774 and field_2774 < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= field_2774 and field_2774 < 0.5):
            # 0: "No"
            value = 0
        elif (0.5 <= field_2774 and field_2774 < 1.5):
            # 1: "Yes"
            value = 1
        elif (-1.5 <= field_2774 and field_2774 < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= field_2774 and field_2774 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW on 23 February 2022
def interpret_pregnancy_terminations_miscarriages_stillbirths(
    code_100291=None,
):
    """
    Intepret UK Biobank's data-coding "100291" for data-fields "3849", "3839",
    and "3829".

    Data-Field "3849": "Number of pregnancy terminations"
    Data-Field "3839": "Number of spontaneous miscarriages"
    Data-Field "3829": "Number of stillbirths"

    UK Biobank data-coding "100291" for data-fields "3849", "3839", "3829".
    0-35: "Count of events"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        code_100291 (float): UK Biobank data-field "3849", "3839", or "3829", in
            data-coding "100291"

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(code_100291)) and
        (-3.5 <= code_100291 and code_100291 < 35.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= code_100291 and code_100291 < 35.5):
            # 0-35: "Count of events"
            value = float(code_100291)
        elif (-1.5 <= code_100291 and code_100291 < -0.5):
            # -1: "Do not know"
            value = float("nan")
        elif (-3.5 <= code_100291 and code_100291 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # missing or uninterpretable
        value = float("nan")
    # Return.
    return value


# review: TCW, 31 March 2022
def interpret_age_at_live_birth(
    field_code_100586=None,
):
    """
    Intepret UK Biobank's data-coding 100586 for multiple data-fields.

    Data-Field "2754": "Age at first live birth"
    Data-Field "2764": "Age at last live birth"
    Data-Field "3872": "Age of primiparous women at birth of child"
    UK Biobank data-coding "100586" for multiple data-fields.
    [8 - 65]: "Age in years"
    -4: "Do not remember"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_code_100586 (float): UK Biobank fields in coding 100586,
            representing age at live birth

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_code_100586)) and
        (-4.5 <= field_code_100586 and field_code_100586 < 65.5)
    ):
        # The variable has a valid value.
        if (7.5 <= field_code_100586 and field_code_100586 < 65.5):
            # [8 - 65]: "Age in years"
            value = float(field_code_100586)
        elif (-3.5 <= field_code_100586 and field_code_100586 < -2.5):
            # -3: "Prefer not to answer"
            value = float("nan")
        elif (-4.5 <= field_code_100586 and field_code_100586 < -3.5):
            # -4: "Do not remember"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


# review: TCW on 18 January 2022
def determine_female_menstruation_days(
    sex_text=None,
    imputation_days=None,
    field_3700=None,
    field_3720=None,
):
    """
    Determine count of days since previous menstruation (menstrual period).

    arguments:
        sex_text (str): textual representation of sex selection
        imputation_days (float): count of days for imputation in persons who
            were experiencing menstruation period currently
        field_3700 (float): UK Biobank field 3700, days since previous
            menstruation (menstrual period)
        field_3720 (float): UK Biobank field 3720, whether person was
            experiencing menstruation period on day of assessment

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret days since previous menstruation.
    menstruation_days = interpret_menstruation_days(
        field_3700=field_3700,
    )
    # Interpret whether person was experiencing menstruation period currently.
    menstruation_current = interpret_menstruation_period_current(
        field_3720=field_3720,
    )
    # Comparison.
    # Only define days since previous menstruation for female persons.
    if (
        (sex_text == "female")
    ):
        # Female.
        # Evaluate validity of relevant variables.
        if (not pandas.isna(menstruation_days)):
            # Person has valid value for days in menstruation cycle.
            # Prioritize this variable.
            value = menstruation_days
        elif (
            (not pandas.isna(menstruation_current)) and
            (menstruation_current == 1)
        ):
            # Person has missing value for days in menstruation cycle.
            # Person was experiencing menstruation period currently.
            # Impute value of days in menstruation cycle.
            value = imputation_days
        else:
            # Person does not qualify for any categories.
            value = float("nan")
    else:
        # Male.
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_menstruation_cycle_duration(
    sex_text=None,
    field_3710=None,
):
    """
    Determine count of days in female person's usual menstrual cycle.

    arguments:
        sex_text (str): textual representation of sex selection
        field_3710 (float): UK Biobank field 3710, count of days in usual
            menstrual cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret days since previous menstruation.
    menstruation_cycle_duration = interpret_menstruation_cycle_duration(
        field_3710=field_3710,
    )
    # Comparison.
    # Only define days since previous menstruation for female persons.
    if (
        (sex_text == "female")
    ):
        # Female.
        value = menstruation_cycle_duration
    else:
        # Male.
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_menstruation_cycle_irregularity(
    sex_text=None,
    field_3710=None,
):
    """
    Determine count of days in female person's usual menstrual cycle.

    arguments:
        sex_text (str): textual representation of sex selection
        field_3710 (float): UK Biobank field 3710, count of days in usual
            menstrual cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret days since previous menstruation.
    menstruation_cycle_irregularity = interpret_menstruation_cycle_irregularity(
        field_3710=field_3710,
    )
    # Comparison.
    # Only define days since previous menstruation for female persons.
    if (
        (sex_text == "female")
    ):
        # Female.
        value = menstruation_cycle_irregularity
    else:
        # Male.
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_hysterectomy(
    sex_text=None,
    field_2724=None,
    field_3591=None,
):
    """
    Determine whether female persons experienced hysterectomy.

    arguments:
        sex_text (str): textual representation of sex
        field_2724 (float): UK Biobank field 2724, whether person has
            experienced menopause
        field_3591 (float): UK Biobank field 3591, whether person has had a
            hysterectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret menopause.
    menopause_hysterectomy = interpret_menopause_hysterectomy(
        field_2724=field_2724,
    )
    # Interpret hysterectomy.
    hysterectomy = interpret_hysterectomy(
        field_3591=field_3591,
    )
    # Comparison.
    if (sex_text == "female"):
        value = utility.determine_logical_or_combination_binary_missing(
            first=menopause_hysterectomy,
            second=hysterectomy,
            single_false_sufficient=True, # both variables have similar meanings
        )
    else:
        # Hysterectomy is undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_oophorectomy(
    sex_text=None,
    field_2834=None,
):
    """
    Determine whether female persons experienced bilateral oophorectomy.

    arguments:
        sex_text (str): textual representation of sex
        field_2834 (float): UK Biobank field 2834, whether person has had an
            oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret oophorectomy.
    oophorectomy = interpret_oophorectomy(
        field_2834=field_2834,
    )
    # Comparison.
    if (sex_text == "female"):
        value = oophorectomy
    else:
        # Oophorectomy is undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_hysterectomy_or_oophorectomy(
    sex_text=None,
    hysterectomy=None,
    oophorectomy=None,
):
    """
    Determine whether female persons experienced hysterectomy.

    arguments:
        sex_text (str): textual representation of sex
        hysterectomy (float): binary logical representation of whether person
            has had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            has had a bilateral oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Comparison.
    if (sex_text == "female"):
        value = utility.determine_logical_or_combination_binary_missing(
            first=hysterectomy,
            second=oophorectomy,
            single_false_sufficient=False, # variables have different meanings
        )
    else:
        # Hysterectomy and oophorectomy are undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_menopause_self_report(
    sex_text=None,
    field_2724=None,
):
    """
    Determine whether female person self reported natural menopause.
    Natural menopause does not include hysterectomy or oophorectomy.

    arguments:
        sex_text (str): textual representation of sex selection
        field_2724 (float): UK Biobank field 2724, whether person self reported
            menopause

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret self report of natural menopause.
    menopause_self_report = interpret_menopause_self_report(
        field_2724=field_2724,
    )
    # Comparison.
    # Only define days since previous menstruation for female persons.
    if (
        (sex_text == "female")
    ):
        # Female.
        value = menopause_self_report
    else:
        # Male.
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 18 January 2022
def determine_female_pregnancy(
    sex_text=None,
    field_3140=None,
):
    """
    Determine whether female persons were pregnant.

    arguments:
        sex_text (str): textual representation of sex selection
        field_3140 (float): UK Biobank field 3140, whether person was pregnant

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret pregnancy.
    pregnancy = interpret_pregnancy(
        field_3140=field_3140,
    )
    # Comparison.
    if (sex_text == "female"):
        value = pregnancy
    else:
        # Pregnancy undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 20 January 2022
def determine_female_menstruation_cycle_regular_duration_range(
    sex_text=None,
    menstruation_irregularity=None,
    menstruation_duration=None,
    menstruation_days=None,
    threshold_duration_low=None,
    threshold_duration_high=None,
):
    """
    Determine whether female person experienced regular menstrual cycles of
    duration within low and high thresholds.

    Hence this definition combines information from multiple variables to
    indicate regularity in menstrual cycle and duration of the menstrual cycle
    within thresholds.

    Notice that variable "menstruation_irregularity" has value "0" if there is a
    valid value for "menstruation_duration".

    Notice that due to the repeating nature of the menstrual cycle, variable
    "menstruation_days" can have values below the lower threshold even when a
    regular cycle has a full duration within both thresholds; however, if
    "menstruation_days" has a value above the higher threshold then that
    information suffices to determine that the menstrual cycle does not have a
    regular duration within the threshold range.

    code:
    1: female person experienced regular menstrual cycles of duration within
    low and high thresholds
    (low <= full duration and full duration < high)
    (current days < high)
    0: female person experienced irregular menstrual cycles or regular cycles of
    duration beyond low or high thresholds
    (low > duration or duration >= high)
    (current days >= high)
    missing: missing information prevented interpretation

    arguments:
        sex_text (str): textual representation of sex selection
        menstruation_irregularity (float): binary logical representation of
            whether person had irregularity in menstrual cycle
        menstruation_duration (float): count of days in person's regular or
            usual menstrual cycle
        menstruation_days (float): count of days since person's previous
            menstrual period or start of menstrual cycle
        threshold_duration_low (float): low threshold in days on the
            duration of the menstrual cycle
        threshold_duration_high (float): high threshold in days on the
            duration of the menstrual cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Determine menopause status.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine whether female person experienced regular menstrual cycles
        # of duration within low and high thresholds.
        if (
            (
                (pandas.isna(menstruation_irregularity)) or
                (menstruation_irregularity == 0)
            ) and
            (
                (not pandas.isna(menstruation_duration)) and
                (menstruation_duration >= threshold_duration_low) and
                (menstruation_duration < threshold_duration_high)
            ) and
            (
                (pandas.isna(menstruation_days)) or
                (menstruation_days < threshold_duration_high)
            )
        ):
            # Female person experienced regular menstrual cycles of duration
            # within low and high thresholds.
            value = 1
        elif (
            (
                (not pandas.isna(menstruation_irregularity)) and
                (menstruation_irregularity == 1)
            ) or
            (
                (not pandas.isna(menstruation_duration)) and
                (
                    (menstruation_duration < threshold_duration_low) or
                    (menstruation_duration >= threshold_duration_high)
                )
            ) or
            (
                (not pandas.isna(menstruation_days)) and
                (menstruation_days >= threshold_duration_high)
            )
        ):
            # Female person experienced irregular menstrual cycles or regular
            # menstrual cycles of duration beyond low or high thresholds.
            value = 0
        else:
            # Female person did not match definitions due to missing or
            # uninterpretable information.
            value = float("nan")
    else:
        # Menopause is undefined for male persons.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 20 January 2022
def determine_female_menopause_ordinal(
    sex_text=None,
    age=None,
    menopause_self_report=None,
    hysterectomy=None,
    oophorectomy=None,
    menstruation_regular_range=None,
    threshold_age_low=None,
    threshold_age_high=None,
):
    """
    Determine whether female persons qualify for pre-menopause, peri-menopause,
    or post-menopause.

    The pre-, peri-, and post-menopause definitions are mutually exclusive.
    These definitions are permissive of missing values in all variables except
    for age.

    Pay close attention to "and", "or" logic in comparisons.

    Notice that age, self report, and oophorectomy criteria are sufficient for
    the definition of post-menopause due to the use of "or" logic.

    This definition uses an ordinal code.
    0: pre-menopause
    1: peri-menopause
    2: post-menopause

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        menopause_self_report (float): binary logical representation of whether
            person self reported having experienced menopause
        hysterectomy (float): binary logical representation of whether person
            had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            had a bilateral oophorectomy
        menstruation_regular_range (float): binary logical representation of
            whether person experienced regular menstrual cycles of duration
            within low and high thresholds
        threshold_age_low (int): threshold age in years, below which to consider
            female persons premenopausal
        threshold_age_high (int): threshold age in years, above which to
            consider female persons postmenopausal

    raises:

    returns:
        (float): interpretation value

    """

    # Determine menopause status.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine pre-menopause.
        if (
            (
                (not pandas.isna(age)) and
                (age < threshold_age_low)
            ) and
            (
                (pandas.isna(menopause_self_report)) or
                (menopause_self_report == 0)
            ) and
            (
                (pandas.isna(oophorectomy)) or
                (oophorectomy == 0)
            ) and
            (
                (pandas.isna(menstruation_regular_range)) or
                (menstruation_regular_range == 1)
            )
        ):
            # Person qualifies for premenopause.
            value = 0
        # Determine peri-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (threshold_age_low <= age and age < threshold_age_high)
            ) and
            (
                (pandas.isna(menopause_self_report)) or
                (menopause_self_report == 0)
            ) and
            (
                (pandas.isna(oophorectomy)) or
                (oophorectomy == 0)
            )
        ):
            # Person qualifies for peri-menopause.
            value = 1
        # Determine post-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (age >= threshold_age_high)
            ) or
            (
                (not pandas.isna(menopause_self_report)) and
                (menopause_self_report == 1)
            ) or
            (
                (not pandas.isna(oophorectomy)) and
                (oophorectomy == 1)
            )
        ):
            # Person qualifies for post-menopause.
            value = 2
        else:
            # Person does not qualify for any groups.
            value = float("nan")
    else:
        # Menopause is undefined for male persons.
        value = float("nan")
    # Return information.
    return value


#determine_female_menopause_ordinal_age_only
#determine_female_menopause_ordinal_age_self_only
#determine_female_menopause_ordinal_age_self_oophorectomy_only


def determine_female_menopause_ordinal_age_only(
    sex_text=None,
    age=None,
    menopause_self_report=None,
    hysterectomy=None,
    oophorectomy=None,
    menstruation_regular_range=None,
    threshold_age_low=None,
    threshold_age_high=None,
):
    """
    Determine whether female persons qualify for pre-menopause, peri-menopause,
    or post-menopause.

    The pre-, peri-, and post-menopause definitions are mutually exclusive.
    These definitions are permissive of missing values in all variables except
    for age.

    Pay close attention to "and", "or" logic in comparisons.

    Notice that age, self report, and oophorectomy criteria are sufficient for
    the definition of post-menopause due to the use of "or" logic.

    This definition uses an ordinal code.
    0: pre-menopause
    1: peri-menopause
    2: post-menopause

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        menopause_self_report (float): binary logical representation of whether
            person self reported having experienced menopause
        hysterectomy (float): binary logical representation of whether person
            had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            had a bilateral oophorectomy
        menstruation_regular_range (float): binary logical representation of
            whether person experienced regular menstrual cycles of duration
            within low and high thresholds
        threshold_age_low (int): threshold age in years, below which to consider
            female persons premenopausal
        threshold_age_high (int): threshold age in years, above which to
            consider female persons postmenopausal

    raises:

    returns:
        (float): interpretation value

    """

    # Determine menopause status.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine pre-menopause.
        if (
            (
                (not pandas.isna(age)) and
                (age < threshold_age_low)
            )
        ):
            # Person qualifies for premenopause.
            value = 0
        # Determine peri-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (threshold_age_low <= age and age < threshold_age_high)
            )
        ):
            # Person qualifies for peri-menopause.
            value = 1
        # Determine post-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (age >= threshold_age_high)
            )
        ):
            # Person qualifies for post-menopause.
            value = 2
        else:
            # Person does not qualify for any groups.
            value = float("nan")
    else:
        # Menopause is undefined for male persons.
        value = float("nan")
    # Return information.
    return value


def determine_female_menopause_ordinal_age_self_only(
    sex_text=None,
    age=None,
    menopause_self_report=None,
    hysterectomy=None,
    oophorectomy=None,
    menstruation_regular_range=None,
    threshold_age_low=None,
    threshold_age_high=None,
):
    """
    Determine whether female persons qualify for pre-menopause, peri-menopause,
    or post-menopause.

    The pre-, peri-, and post-menopause definitions are mutually exclusive.
    These definitions are permissive of missing values in all variables except
    for age.

    Pay close attention to "and", "or" logic in comparisons.

    Notice that age, self report, and oophorectomy criteria are sufficient for
    the definition of post-menopause due to the use of "or" logic.

    This definition uses an ordinal code.
    0: pre-menopause
    1: peri-menopause
    2: post-menopause

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        menopause_self_report (float): binary logical representation of whether
            person self reported having experienced menopause
        hysterectomy (float): binary logical representation of whether person
            had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            had a bilateral oophorectomy
        menstruation_regular_range (float): binary logical representation of
            whether person experienced regular menstrual cycles of duration
            within low and high thresholds
        threshold_age_low (int): threshold age in years, below which to consider
            female persons premenopausal
        threshold_age_high (int): threshold age in years, above which to
            consider female persons postmenopausal

    raises:

    returns:
        (float): interpretation value

    """

    # Determine menopause status.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine pre-menopause.
        if (
            (
                (not pandas.isna(age)) and
                (age < threshold_age_low)
            ) and
            (
                (pandas.isna(menopause_self_report)) or
                (menopause_self_report == 0)
            )
        ):
            # Person qualifies for premenopause.
            value = 0
        # Determine peri-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (threshold_age_low <= age and age < threshold_age_high)
            ) and
            (
                (pandas.isna(menopause_self_report)) or
                (menopause_self_report == 0)
            )
        ):
            # Person qualifies for peri-menopause.
            value = 1
        # Determine post-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (age >= threshold_age_high)
            ) or
            (
                (not pandas.isna(menopause_self_report)) and
                (menopause_self_report == 1)
            )
        ):
            # Person qualifies for post-menopause.
            value = 2
        else:
            # Person does not qualify for any groups.
            value = float("nan")
    else:
        # Menopause is undefined for male persons.
        value = float("nan")
    # Return information.
    return value


def determine_female_menopause_ordinal_age_self_oophorectomy_only(
    sex_text=None,
    age=None,
    menopause_self_report=None,
    hysterectomy=None,
    oophorectomy=None,
    menstruation_regular_range=None,
    threshold_age_low=None,
    threshold_age_high=None,
):
    """
    Determine whether female persons qualify for pre-menopause, peri-menopause,
    or post-menopause.

    The pre-, peri-, and post-menopause definitions are mutually exclusive.
    These definitions are permissive of missing values in all variables except
    for age.

    Pay close attention to "and", "or" logic in comparisons.

    Notice that age, self report, and oophorectomy criteria are sufficient for
    the definition of post-menopause due to the use of "or" logic.

    This definition uses an ordinal code.
    0: pre-menopause
    1: peri-menopause
    2: post-menopause

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        menopause_self_report (float): binary logical representation of whether
            person self reported having experienced menopause
        hysterectomy (float): binary logical representation of whether person
            had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            had a bilateral oophorectomy
        menstruation_regular_range (float): binary logical representation of
            whether person experienced regular menstrual cycles of duration
            within low and high thresholds
        threshold_age_low (int): threshold age in years, below which to consider
            female persons premenopausal
        threshold_age_high (int): threshold age in years, above which to
            consider female persons postmenopausal

    raises:

    returns:
        (float): interpretation value

    """

    # Determine menopause status.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine pre-menopause.
        if (
            (
                (not pandas.isna(age)) and
                (age < threshold_age_low)
            ) and
            (
                (pandas.isna(menopause_self_report)) or
                (menopause_self_report == 0)
            ) and
            (
                (pandas.isna(oophorectomy)) or
                (oophorectomy == 0)
            )
        ):
            # Person qualifies for premenopause.
            value = 0
        # Determine peri-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (threshold_age_low <= age and age < threshold_age_high)
            ) and
            (
                (pandas.isna(menopause_self_report)) or
                (menopause_self_report == 0)
            ) and
            (
                (pandas.isna(oophorectomy)) or
                (oophorectomy == 0)
            )
        ):
            # Person qualifies for peri-menopause.
            value = 1
        # Determine post-menopause.
        elif (
            (
                (not pandas.isna(age)) and
                (age >= threshold_age_high)
            ) or
            (
                (not pandas.isna(menopause_self_report)) and
                (menopause_self_report == 1)
            ) or
            (
                (not pandas.isna(oophorectomy)) and
                (oophorectomy == 1)
            )
        ):
            # Person qualifies for post-menopause.
            value = 2
        else:
            # Person does not qualify for any groups.
            value = float("nan")
    else:
        # Menopause is undefined for male persons.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 20 January 2022
def determine_female_menstruation_phase(
    sex_text=None,
    menstruation_regular_range=None,
    menstruation_days=None,
    threshold_ovulation=None,
):
    """
    Determine a female person's categorical menstruation phase.

    Only consider for female persons who had a regular menstrual cycle of
    duration within threshold range (menstruation_regular_range = 1).

    0: follicular phase of menstruation
    1: luteal phase of menstruation

    arguments:
        sex_text (str): textual representation of sex selection
        menstruation_regular_range (float): binary logical representation of
            whether person experienced regular menstrual cycles of duration
            within low and high thresholds
        menstruation_days (float): count of days since previous menstruation
            (menstrual period, start of menstrual cycle)
        threshold_ovulation (int): threshold in days for the approximate
            "middle" of menstruation cycle at ovulation

    raises:

    returns:
        (float): interpretation value

    """

    # Determine categorical phase of menstruation cycle.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(menstruation_regular_range)) and
            (not pandas.isna(menstruation_days))
        ):
            # Determine categorical phase of menstruation cycle.
            if (menstruation_regular_range == 1):
                if (menstruation_days < threshold_ovulation):
                    # Person qualifies for follicular phase of menstrual cycle.
                    value = 0
                elif (menstruation_days >= threshold_ovulation):
                    # Person qualitifies for luteal phase of menstrual cycle.
                    value = 1
                else:
                    # Persons does not qualify for any categories.
                    value = float("nan")
            else:
                # Person does not qualify for any categories.
                value = float("nan")
        else:
            # Necessary variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 20 January 2022
def determine_female_menstruation_phase_early_late(
    sex_text=None,
    menstruation_regular_range=None,
    menstruation_days=None,
    threshold_middle_follicular=None,
    threshold_ovulation=None,
    threshold_middle_luteal=None,
):
    """
    Determine a female person's categorical menstruation phase.

    Only consider for female persons who had a regular menstrual cycle of
    duration within threshold range (menstruation_regular_range = 1).

    0: early follicular phase of menstruation
    - - (day 0 to day threshold_middle_follicular)
    1: late follicular phase of menstruation
    - - (day threshold_middle_follicular to day threshold_ovulation)
    2: early luteal phase of menstruation
    - - (day threshold_ovulation to day threshold_middle_luteal)
    3: late luteal phase of menstruation
    - - (day threshold_middle_luteal to maximum day menstruation_days)

    arguments:
        sex_text (str): textual representation of sex selection
        menstruation_regular_range (float): binary logical representation of
            whether person experienced regular menstrual cycles of duration
            within low and high thresholds
        menstruation_days (float): count of days since previous menstruation
            (menstrual period, start of menstrual cycle)
        threshold_middle_follicular (int): threshold in days for middle of
            follicular phase of menstruation cycle
        threshold_ovulation (int): threshold in days for the approximate
            "middle" of menstruation cycle at ovulation
        threshold_middle_luteal (int): threshold in days for middle of luteal
            phase of menstruation cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Determine categorical phase of menstruation cycle.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(menstruation_regular_range)) and
            (not pandas.isna(menstruation_days))
        ):
            # Determine categorical phase of menstrual cycle.
            if (menstruation_regular_range == 1):
                if (
                    (menstruation_days < threshold_middle_follicular)
                ):
                    # Person qualifies for early follicular phase.
                    value = 0
                elif (
                    (menstruation_days >= threshold_middle_follicular) and
                    (menstruation_days < threshold_ovulation)
                ):
                    # Person qualitifies for late follicular phase.
                    value = 1
                elif (
                    (menstruation_days >= threshold_ovulation) and
                    (menstruation_days < threshold_middle_luteal)
                ):
                    # Person qualitifies for early luteal phase.
                    value = 2
                elif (
                    (menstruation_days >= threshold_middle_luteal)
                ):
                    # Person qualitifies for late luteal phase.
                    value = 3
                else:
                    # Persons does not qualify for any categories.
                    value = float("nan")
            else:
                # Person does not qualify for any categories.
                value = float("nan")
        else:
            # Necessary variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 20 January 2022
def determine_female_menstruation_phase_cycle(
    sex_text=None,
    menstruation_phase_early_late=None,
):
    """
    Determine a female person's menstruation phase from a cyclical perspective.

    0: early follicular or late luteal phase of menstrual cycle
    - "shoulder"
    1: late follicular or early luteal phase of menstrual cycle
    - "ovulation"

    arguments:
        sex_text (str): textual representation of sex selection
        menstruation_phase_early_late (float): ordinal representation of phases
            of the menstrual cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Determine categorical phase of menstruation cycle.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(menstruation_phase_early_late))
        ):
            # Determine categorical phase of menstrual cycle.
            if (
                (menstruation_phase_early_late == 0) or
                (menstruation_phase_early_late == 3)
            ):
                # Person qualifies for early follicular phase or late luteal
                # phase.
                # "shoulder"
                value = 0
            elif (
                (menstruation_phase_early_late == 1) or
                (menstruation_phase_early_late == 2)
            ):
                # Person qualitifies for late follicular or early luteal
                # phase.
                # "ovulation"
                value = 1
            else:
                # Persons does not qualify for any categories.
                value = float("nan")
        else:
            # Necessary variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 23 February 2022
def determine_female_births_live_still_count(
    sex_text=None,
    field_2734=None,
    field_2774=None,
    field_3829=None,
):
    """
    Determine total count of pregnancies that progressed to full term or nearly
    full term and ended in either a still birth or a live birth.

    arguments:
        sex_text (str): textual representation of sex selection
        field_2734 (float): UK Biobank field 2734, count of live births
        field_2774 (float): UK Biobank field 2774, whether person experienced
            any type of loss of a pregnancy that did not end in live birth
        field_3829 (float): UK Biobank field 3829, count of still births

    raises:

    returns:
        (float): interpretation value

    """

    # Interpretation.
    pregnancy_loss = interpret_pregnancy_loss(
        field_2774=field_2774,
    )
    births_still = interpret_pregnancy_terminations_miscarriages_stillbirths(
        code_100291=field_3829, # count still births
    )
    births_live = interpret_live_births(
        field_2734=field_2734,
    )
    # Comparison.
    if (
        (sex_text == "female")
    ):
        # Determine count of still births.
        if (
            (not pandas.isna(pregnancy_loss)) and
            (pregnancy_loss == 1)
        ):
            # Female person experienced a pregnancy loss.
            # Determine the count of pregnancies that ended in still birth.
            if (not pandas.isna(births_still)):
                count_still = births_still
            else:
                # There is inadequate explanation of the pregnancy loss.
                # Do not make any assumptions.
                # Lost pregnancy might have ended early in termination or
                # abortion or might have ended late in still birth.
                count_still = float("nan")
        elif (
            (not pandas.isna(pregnancy_loss)) and
            (pregnancy_loss == 0)
        ):
            # Female person did not experience any pregnancy loss.
            # Hence female person did not experience any still births.
            count_still = 0
        else:
            count_still = float("nan")
        # Determine count of live births.
        if (not pandas.isna(births_live)):
            count_live = births_live
        else:
            count_live = float("nan")
        # Determine count of total still and live births.
        if (
            (not pandas.isna(count_still)) and
            (not pandas.isna(count_live))
        ):
            value = (count_still + count_live)
        elif (not pandas.isna(count_still)):
            value = count_still
        elif (not pandas.isna(count_live)):
            value = count_live
        else:
            value = float("nan")
    else:
        # Pregnancy undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 30 March 2022
def determine_female_parity_births_any(
    sex_text=None,
    births=None,
):
    """
    Determine whether female person experienced any pregnancies that progressed
    to full term or nearly full term and ended in either a still birth or a live
    birth.

    Definition of this variable is the same as parity.

    code
    0: female person experienced zero still or live births
    1: female person experienced one or more still or live births

    arguments:
        sex_text (str): textual representation of sex selection
        births (float): total count of still or live births

    raises:

    returns:
        (float): interpretation value

    """

    # Comparison.
    if (
        (sex_text == "female")
    ):
        # Determine whether person experienced any births.
        if (
            (not pandas.isna(births)) and
            (births > 0)
        ):
            # Female person experienced one or more births.
            value = 1
        elif (
            (not pandas.isna(births)) and
            (births == 0)
        ):
            # Female person experienced zero births.
            value = 0
        else:
            value = float("nan")
    else:
        # Pregnancy undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 31 March 2022
def determine_female_age_at_most_recent_live_birth(
    sex_text=None,
    field_2754=None,
    field_2764=None,
    field_3872=None,
):
    """
    Determine the age in years at which female person experienced most recent
    live birth.

    UK Biobank does not include information about age at still birth.

    arguments:
        sex_text (str): textual representation of sex selection
        field_2754 (float): UK Biobank field 2754, age at first live birth
        field_2764 (float): UK Biobank field 2764, age at last live birth
        field_3872 (float): UK Biobank field 3872, age at only live birth

    raises:

    returns:
        (float): interpretation value

    """

    # Interpretation.
    age_first = interpret_age_at_live_birth(
        field_code_100586=field_2754,
    )
    age_last = interpret_age_at_live_birth(
        field_code_100586=field_2764,
    )
    age_only = interpret_age_at_live_birth(
        field_code_100586=field_3872,
    )

    # Combine ages.
    if (
        (not pandas.isna(age_first)) or
        (not pandas.isna(age_last)) or
        (not pandas.isna(age_only))
    ):
        #array = numpy.array([age_first, age_last, age_only])
        #value = numpy.nanmax(array)
        ages = list(filter(
            lambda age: (not math.isnan(age)),
            [age_first, age_last, age_only]
        ))
        value = max(ages)
    else:
        value = float("nan")
    # Return information.
    return value


# review: TCW on 31 March 2022
def determine_female_years_since_most_recent_live_birth(
    sex_text=None,
    age=None,
    age_last_live_birth=None,
):
    """
    Determine duration in years since female person experienced their most
    recent live birth.

    This function assigns missing values for female persons for whom information
    was missing about their age at most recent live birth.

    UK Biobank does not include information about age at still birth.

    arguments:
        sex_text (str): textual representation of sex selection
        age (float): age of person in years
        age_last_live_birth (float): age of female person at most recent live
            birth

    raises:

    returns:
        (float): interpretation value

    """

    # Comparison.
    if (
        (sex_text == "female")
    ):
        # Determine whether there is sufficient information to determine
        # duration in years since most recent live birth.
        if (
            (not pandas.isna(age)) and
            (not pandas.isna(age_last_live_birth)) and
            (age_last_live_birth <= age)
        ):
            value = (age - age_last_live_birth)
        else:
            # There is insufficient information.
            value = float("nan")
    else:
        # Pregnancy and birth undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 31 March 2022
def determine_female_birth_live_recent(
    sex_text=None,
    births_any=None,
    years_since_live_birth=None,
    threshold=None,
):
    """
    Determine whether female person experienced a live birth within specific
    threshold duration (in years) of current age (in years) at assessment.

    If a female person never experienced any still or live births, then that
    person also did not experience any recent live births.

    arguments:
        sex_text (str): textual representation of sex selection
        births_any (float): whether person had experienced any births, still or
            live
        years_since_live_birth (float): duration in years since most recent live
            birth
        threshold (float): threshold duration in years since a recent live birth

    raises:

    returns:
        (float): interpretation value

    """

    # Comparison.
    if (
        (sex_text == "female")
    ):
        # Determine whether person had experienced any births (live or still).
        if (
            (not pandas.isna(births_any)) and
            (births_any == 1)
        ):
            # Female person had experienced one or more births (live or still).
            # Determine whether there is sufficient information to determine
            # whether any live births were recent.
            if (not pandas.isna(years_since_live_birth)):
                if (years_since_live_birth < threshold):
                    value = 1
                else:
                    value = 0
            else:
                # There is insufficient information.
                value = float("nan")
        elif (
            (not pandas.isna(births_any)) and
            (births_any == 0)
        ):
            # Female person had not experienced any births (live or still).
            # Hence female person had not experienced a recent birth.
            value = 0
        else:
            value = float("nan")
    else:
        # Pregnancy and birth undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 23 February 2022
def determine_female_pregnancies_early_count(
    sex_text=None,
    field_2774=None,
    field_3849=None,
    field_3839=None,
):
    """
    Determine total count of pregnancies that ended early in either abortive
    termination or spontaneous miscarriage.

    arguments:
        sex_text (str): textual representation of sex selection
        field_2774 (float): UK Biobank field 2774, whether person experienced
            any type of loss of a pregnancy that did not end in live birth
        field_3849 (float): UK Biobank field 3849, count of pregnancy
            abortive terminations
        field_3839 (float): UK Biobank field 3839, count of pregnancy
            spontaneous miscarriages

    raises:

    returns:
        (float): interpretation value

    """

    # Interpretation.
    pregnancy_loss = interpret_pregnancy_loss(
        field_2774=field_2774,
    )
    terminations = interpret_pregnancy_terminations_miscarriages_stillbirths(
        code_100291=field_3849, # count terminations
    )
    miscarriages = interpret_pregnancy_terminations_miscarriages_stillbirths(
        code_100291=field_3839, # count miscarriages
    )

    # Comparison.
    if (
        (sex_text == "female")
    ):
        if (
            (not pandas.isna(pregnancy_loss)) and
            (pregnancy_loss == 1)
        ):
            # Female person experienced a pregnancy loss.
            # Determine the count of pregnancies that ended early.
            if (
                (not pandas.isna(terminations)) and
                (not pandas.isna(miscarriages))
            ):
                value = (terminations + miscarriages)
            elif (not pandas.isna(terminations)):
                value = terminations
            elif (not pandas.isna(miscarriages)):
                value = miscarriages
            else:
                # There is inadequate explanation of the pregnancy loss.
                # Do not make any assumptions.
                # Lost pregnancy might have ended early in termination or
                # abortion or might have ended late in still birth.
                value = float("nan")
        elif (
            (not pandas.isna(pregnancy_loss)) and
            (pregnancy_loss == 0)
        ):
            # Female person did not experience any pregnancy loss.
            # Hence female person did not experience any early pregnancy loss.
            value = 0
        else:
            value = float("nan")
    else:
        # Pregnancy undefined for males.
        value = float("nan")
    # Return information.
    return value


# review: TCW on 23 February 2022
def determine_female_pregnancies_total_count(
    sex_text=None,
    pregnancies_early=None,
    births=None,
):
    """
    Determine total count of pregnancies that female person experienced.

    arguments:
        sex_text (str): textual representation of sex selection
        pregnancies_early (float): total count of pregnancies that ended early
            in either abortive termination or spontaneous miscarriage
        births (float): total count of live and still births

    raises:

    returns:
        (float): interpretation value

    """

    # Comparison.
    if (sex_text == "female"):
        if (
            (not pandas.isna(pregnancies_early)) and
            (not pandas.isna(births))
        ):
            value = (pregnancies_early + births)
        elif (not pandas.isna(pregnancies_early)):
            value = pregnancies_early
        elif (not pandas.isna(births)):
            value = births
        else:
            # There is insufficient information.
            value = float("nan")
            pass
    else:
        # Pregnancy undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_age_event_simple(
    sex_text=None,
    field_raw_age=None,
):
    """
    Determines the age of a female person at a major life event.

    arguments:
        sex_text (str): textual representation of sex, 'female' or 'male'
        field_raw_age (float): raw data field for age at event

    raises:

    returns:
        (float): age of female person at event

    """

    if (sex_text == "female"):
        # Determine age.
        interpretation = determine_age(
            field_raw=field_raw_age,
        )
    elif (sex_text == "male"):
        # Interpretation is specific to female sex.
        interpretation = float("nan")
    else:
        # Interpretation is null without definition of sex.
        interpretation = float("nan")
        pass
    # Return information.
    return interpretation


def determine_female_age_menopause_never_oophorectomy(
    sex_text=None,
    oophorectomy=None,
    field_raw_age=None,
):
    """
    Determines the age of a female person at their self report of menopause,
    with consideration only of female persons who never experienced bilateral
    oophorectomy.

    arguments:
        sex_text (str): textual representation of sex, 'female' or 'male'
        oophorectomy (float): whether female person experienced oophorectomy
        field_raw_age (float): raw data field for age at event

    raises:

    returns:
        (float): age of female person at event

    """

    if (sex_text == "female"):
        if (oophorectomy == 0):
            # Female person who never experienced oophorectomy.
            # Determine age.
            interpretation = determine_age(
                field_raw=field_raw_age,
            )
        else:
            # Interpretation is null in female persons who experienced
            # oophorectomy.
            interpretation = float("nan")
    elif (sex_text == "male"):
        # Interpretation is specific to female sex.
        interpretation = float("nan")
    else:
        # Interpretation is null without definition of sex.
        interpretation = float("nan")
        pass
    # Return information.
    return interpretation


# review: TCW on 19 May 2022
def organize_female_ages_major_events_relevant_to_menstruation(
    table=None,
    report=None,
):
    """
    Organizes ages of female persons at major life events that are relevant to
    menstruation.

    TCW; 19 May 2022
    It was necessary to exclude "age_menopause" as the necessary data-field
    from UK Biobank was not available.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Determine simple ages at events.
    #table["age_menarche"] = table.apply(
    #    lambda row:
    #        determine_female_age_event_simple(
    #            sex_text=row["sex_text"],
    #            field_raw_age=row["2714-0.0"], # age at menarche
    #        ),
    #    axis="columns", # apply function to each row
    #)
    table["age_menopause_self_report"] = table.apply(
        lambda row:
            determine_female_age_event_simple(
                sex_text=row["sex_text"],
                field_raw_age=row["3581-0.0"], # age at menopause
            ),
        axis="columns", # apply function to each row
    )
    table["age_menopause_never_oophorectomy"] = table.apply(
        lambda row:
            determine_female_age_menopause_never_oophorectomy(
                sex_text=row["sex_text"],
                oophorectomy=row["oophorectomy"],
                field_raw_age=row["3581-0.0"], # age at menopause
            ),
        axis="columns", # apply function to each row
    )
    table["age_oophorectomy"] = table.apply(
        lambda row:
            determine_female_age_event_simple(
                sex_text=row["sex_text"],
                field_raw_age=row["3882-0.0"], # age at bilateral oophorectomy
            ),
        axis="columns", # apply function to each row
    )
    table["age_hysterectomy"] = table.apply(
        lambda row:
            determine_female_age_event_simple(
                sex_text=row["sex_text"],
                field_raw_age=row["2824-0.0"], # age at hysterectomy
            ),
        axis="columns", # apply function to each row
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_female_ages_major_events_relevant_to_menstruation()")
        utility.print_terminal_partition(level=3)
        # ...
        pass
    # Return information.
    return table


# review: TCW on 30 March 2022
def report_female_pregnancies_births(
    table=None,
):
    """
    Reports counts of female persons who experienced events of pregnancy and
    birth.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_pregnancies_0 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancies"] == 0)
        ), :
    ]
    table_pregnancies_1 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancies"] == 1)
        ), :
    ]
    table_pregnancies_2 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancies"] == 2)
        ), :
    ]
    table_pregnancies_3_or_more = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancies"] >= 3)
        ), :
    ]


    table_births_0 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["births"] == 0)
        ), :
    ]
    table_births_1 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["births"] == 1)
        ), :
    ]
    table_births_2 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["births"] == 2)
        ), :
    ]
    table_births_3_or_more = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["births"] >= 3)
        ), :
    ]

    table_births_any_0 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["births_any"] == 0)
        ), :
    ]
    table_births_any_1 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["births_any"] == 1)
        ), :
    ]

    table_birth_recent_0 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["birth_live_recent"] == 0)
        ), :
    ]
    table_birth_recent_1 = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["birth_live_recent"] == 1)
        ), :
    ]

    # Counts.
    count_female = table_female.shape[0]

    count_pregnancies_0 = table_pregnancies_0.shape[0]
    count_pregnancies_1 = table_pregnancies_1.shape[0]
    count_pregnancies_2 = table_pregnancies_2.shape[0]
    count_pregnancies_3_or_more = table_pregnancies_3_or_more.shape[0]

    count_births_0 = table_births_0.shape[0]
    count_births_1 = table_births_1.shape[0]
    count_births_2 = table_births_2.shape[0]
    count_births_3_or_more = table_births_3_or_more.shape[0]

    count_births_any_0 = table_births_any_0.shape[0]
    count_births_any_1 = table_births_any_1.shape[0]

    count_birth_recent_0 = table_birth_recent_0.shape[0]
    count_birth_recent_1 = table_birth_recent_1.shape[0]

    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_female_pregnancies_births()")
    utility.print_terminal_partition(level=3)
    print("... all counts specific to female persons ...")
    utility.print_terminal_partition(level=4)
    print("total: " + str(count_female))
    utility.print_terminal_partition(level=4)
    print("Pregnancies...")
    utility.print_terminal_partition(level=4)
    print("0: " + str(count_pregnancies_0))
    print("1: " + str(count_pregnancies_1))
    print("2: " + str(count_pregnancies_2))
    print("3 or more: " + str(count_pregnancies_3_or_more))
    utility.print_terminal_partition(level=4)
    print("Births...")
    utility.print_terminal_partition(level=4)
    print("0: " + str(count_births_0))
    print("1: " + str(count_births_1))
    print("2: " + str(count_births_2))
    print("3 or more: " + str(count_births_3_or_more))
    utility.print_terminal_partition(level=4)
    print("Births Any (Parity)...")
    utility.print_terminal_partition(level=4)
    print("0: " + str(count_births_any_0))
    print("1: " + str(count_births_any_1))
    utility.print_terminal_partition(level=4)
    print("Recent Live Birth...")
    utility.print_terminal_partition(level=4)
    print("0: " + str(count_birth_recent_0))
    print("1: " + str(count_birth_recent_1))
    pass


# review: TCW on 20 January January 2022
def report_female_menstruation_regularity_duration_range(
    threshold_duration_low=None,
    threshold_duration_high=None,
    table=None,
):
    """
    Reports counts of female persons who have regular or irregular menstrual
    cycles of durations below or above thresholds.

    arguments:
        threshold_duration_low (float): low threshold in days on the
            duration of the menstrual cycle
        threshold_duration_high (float): high threshold in days on the
            duration of the menstrual cycle
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_female_irregularity = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_irregularity"] == 1)
        ), :
    ]
    table_female_regularity = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_irregularity"] == 0)
        ), :
    ]

    table_female_valid_duration = table.loc[
        (
            (table["sex_text"] == "female") &
            (~pandas.isna(table["menstruation_duration"]))
        ), :
    ]
    table_female_valid_current = table.loc[
        (
            (table["sex_text"] == "female") &
            (~pandas.isna(table["menstruation_days"]))
        ), :
    ]
    table_female_valid_duration_or_current = table.loc[
        (
            (table["sex_text"] == "female") &
            ((~pandas.isna(table["menstruation_duration"])) |
            (~pandas.isna(table["menstruation_days"])))
        ), :
    ]

    table_female_regular_range = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    table_female_duration_below_range = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_duration"] < threshold_duration_low)
        ), :
    ]
    table_female_duration_above_range = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_duration"] >= threshold_duration_high)
        ), :
    ]
    table_female_current_days_above_range = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_days"] >= threshold_duration_high)
        ), :
    ]

    table_discrepancy_duration_current = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_duration"] < threshold_duration_high) &
            (table["menstruation_days"] >= threshold_duration_high)
        ), :
    ]
    table_discrepancy_regular_current = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menstruation_regular_range"] == 1) &
            (table["menstruation_days"] >= threshold_duration_high)
        ), :
    ]
    # Counts.

    count_female = table_female.shape[0]

    count_irregularity = table_female_irregularity.shape[0]
    count_regularity = table_female_regularity.shape[0]

    count_valid_duration = table_female_valid_duration.shape[0]
    count_valid_current = table_female_valid_current.shape[0]
    count_valid_duration_or_current = (
        table_female_valid_duration_or_current.shape[0]
    )

    count_regular_range = table_female_regular_range.shape[0]
    count_duration_below = table_female_duration_below_range.shape[0]
    count_duration_above = table_female_duration_above_range.shape[0]
    count_current_above = table_female_current_days_above_range.shape[0]

    count_discrepancy_duration_current = (
        table_discrepancy_duration_current.shape[0]
    )
    count_discrepancy_regular_current = (
        table_discrepancy_regular_current.shape[0]
    )

    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_female_menstruation_regularity_duration_range()")
    utility.print_terminal_partition(level=3)
    print("... all counts specific to female persons ...")
    utility.print_terminal_partition(level=4)
    print("total: " + str(count_female))

    print("'irregular' menstrual cycle: " + str(count_irregularity))
    print("'regular' menstrual cycle: " + str(count_regularity))

    print(
        "non-missing regular duration of menstrual cycle: " +
        str(count_valid_duration)
    )
    print(
        "non-missing days of current menstrual cycle: " +
        str(count_valid_current)
    )
    print(
        "non-missing regular duration or days of current menstrual cycle: " +
        str(count_valid_duration_or_current)
    )

    print(
        "regular menstrual cycle of duration within threshold range: " +
        str(count_regular_range)
    )
    print(
        "regular menstrual cycle duration shorter than " +
        str(threshold_duration_low) +
        " days: " +
        str(count_duration_below)
    )
    print(
        "regular menstrual cycle duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_duration_above)
    )
    print(
        "current menstrual cycle at duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_current_above)
    )
    print(
        "regular menstrual cycle duration shorter than " +
        str(threshold_duration_high) +
        " days, but current menstrual cycle at duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_discrepancy_duration_current)
    )
    print(
        "regular menstrual cycle of duration within " +
        "threshold range, but current menstrual cycle at duration of " +
        str(threshold_duration_high) +
        " days or longer: " +
        str(count_discrepancy_regular_current)
    )
    pass


# review: TCW on 18 January 2022
def report_female_self_report_menopause_hysterectomy_oophorectomy(
    table=None,
):
    """
    Reports counts of female persons who self-reported having experienced
    menopause or having received an oophorectomy.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_female_menopause_self_report = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["menopause_self_report"] == 1)
        ), :
    ]
    table_female_hysterectomy = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["hysterectomy"] == 1)
        ), :
    ]
    table_female_oophorectomy = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["oophorectomy"] == 1)
        ), :
    ]
    table_female_hyster_oophor = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["hysterectomy_or_oophorectomy"] == 1)
        ), :
    ]
    # Counts.
    count_female = table_female.shape[0]
    count_menopause_self = table_female_menopause_self_report.shape[0]
    count_hysterectomy = table_female_hysterectomy.shape[0]
    count_oophorectomy = table_female_oophorectomy.shape[0]
    count_hyster_oophor = table_female_hyster_oophor.shape[0]
    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_female_self_report_menopause_hysterectomy_oophorectomy()")
    utility.print_terminal_partition(level=3)
    print("... all counts specific to female persons ...")
    utility.print_terminal_partition(level=4)
    print("total: " + str(count_female))
    print(
        "self-reported menopause: " +
        str(count_menopause_self)
    )
    print(
        "had hysterectomy: " +
        str(count_hysterectomy)
    )
    print(
        "had oophorectomy: " +
        str(count_oophorectomy)
    )
    print(
        "had either hysterectomy or oophorectomy: " +
        str(count_hyster_oophor)
    )
    pass


# review: TCW on 18 January 2022
def report_female_menopause_three_groups(
    female_menopause=None,
    table=None,
):
    """
    Reports counts of female persons in groups by menopause status.

    arguments:
        female_menopause (str): name of column for definition of female
            menopause status
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_female_premenopause = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[female_menopause] == 0)
        ), :
    ]
    table_female_perimenopause = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[female_menopause] == 1)
        ), :
    ]
    table_female_postmenopause = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[female_menopause] == 2)
        ), :
    ]
    table_female_regular = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    # Counts.
    count_female = table_female.shape[0]
    count_premenopause = table_female_premenopause.shape[0]
    count_perimenopause = table_female_perimenopause.shape[0]
    count_postmenopause = table_female_postmenopause.shape[0]
    count_regular = table_female_regular.shape[0]
    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_female_menopause_three_groups()")
    utility.print_terminal_partition(level=3)
    print("menopause variable: " + str(female_menopause))
    print("... all counts specific to female persons ...")
    print("... all counts except for total females exclude pregnancy...")
    utility.print_terminal_partition(level=4)
    print("total: " + str(count_female))
    print("pre-menopause: " + str(count_premenopause))
    print("peri-menopause: " + str(count_perimenopause))
    print("post-menopause: " + str(count_postmenopause))
    print("regular menstruation: " + str(count_regular))
    pass


# review: TCW on 24 January 2022
def report_combination_female_regular_menstruation_by_menopause(
    menopause=None,
    menstruation_regular_range=None,
    table=None,
):
    """
    Reports counts of female persons in groups defined by combination of
    menopause status and regular menstruation.

    arguments:
        menopause (str): name of column for definition of female menopause
            status
        menstruation_regular_range (float): name of column for binary logical
            representation of whether person experienced regular menstrual
            cycles of duration within low and high thresholds
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    table_regular_premenopause = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[menstruation_regular_range] == 1) &
            (table[menopause] == 0)
        ), :
    ]
    table_regular_perimenopause = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[menstruation_regular_range] == 1) &
            (table[menopause] == 1)
        ), :
    ]
    table_regular_postmenopause = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[menstruation_regular_range] == 1) &
            (table[menopause] == 2)
        ), :
    ]
    # Counts.
    count_female = table_female.shape[0]
    count_regular_premenopause = table_regular_premenopause.shape[0]
    count_regular_perimenopause = table_regular_perimenopause.shape[0]
    count_regular_postmenopause = table_regular_postmenopause.shape[0]
    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    print("report_combination_female_regular_menstruation_by_menopause()")
    utility.print_terminal_partition(level=3)
    print("menopause variable: " + str(menopause))
    print("... all counts specific to female persons ...")
    print("... all counts except for total females exclude pregnancy...")
    utility.print_terminal_partition(level=4)
    print("total: " + str(count_female))
    print(
        "regular menstruation and pre-menopause: " +
        str(count_regular_premenopause)
    )
    print(
        "regular menstruation and peri-menopause: " +
        str(count_regular_perimenopause)
    )
    print(
        "regular menstruation and post-menopause: " +
        str(count_regular_postmenopause)
    )
    pass


def report_alteration_sex_hormones_by_female_menopause(
    menopause=None,
    menstruation_regular_range=None,
    alteration_sex_hormone=None,
    table=None,
):
    """
    Reports counts of persons with stratification by sex who used medications in
    a specific class group.

    arguments:
        menopause (str): name of column for definition of female menopause
            status
        menstruation_regular_range (float): name of column for binary logical
            representation of whether person experienced regular menstrual
            cycles of duration within low and high thresholds
        alteration_sex_hormone (str): name of column for binary logical
            representation of whether person used any medications that influence
            sex hormones, including oral contraception and other hormonal
            therapy in females
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Determine whether the necessary columns exist in the table.
    columns = table.columns.tolist()
    if (
        (menopause in columns) and
        (menstruation_regular_range in columns) and
        (alteration_sex_hormone in columns)
    ):
        # Copy information in table.
        table = table.copy(deep=True)

        # Female total.
        table_female = table.loc[
            (
                (table["sex_text"] == "female")
            ), :
        ]
        # Female menopause.
        table_female_premenopause_alteration = table.loc[
            (
                (table["sex_text"] == "female") &
                (table["pregnancy"] == 0) &
                (table[menopause] == 0) &
                (table[alteration_sex_hormone] == 1)
            ), :
        ]
        table_female_perimenopause_alteration = table.loc[
            (
                (table["sex_text"] == "female") &
                (table["pregnancy"] == 0) &
                (table[menopause] == 1) &
                (table[alteration_sex_hormone] == 1)
            ), :
        ]
        table_female_postmenopause_alteration = table.loc[
            (
                (table["sex_text"] == "female") &
                (table["pregnancy"] == 0) &
                (table[menopause] == 2) &
                (table[alteration_sex_hormone] == 1)
            ), :
        ]
        table_female_menstruation_alteration = table.loc[
            (
                (table["sex_text"] == "female") &
                (table["pregnancy"] == 0) &
                (table[menstruation_regular_range] == 1) &
                (table[alteration_sex_hormone] == 1)
            ), :
        ]
        # Counts.
        count_female = table_female.shape[0]
        count_pre_alter = table_female_premenopause_alteration.shape[0]
        count_peri_alter = table_female_perimenopause_alteration.shape[0]
        count_post_alter = table_female_postmenopause_alteration.shape[0]
        count_regular_alter = table_female_menstruation_alteration.shape[0]
        # Report.
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("report_alteration_sex_hormones_by_female_menopause()")
        utility.print_terminal_partition(level=3)
        utility.print_terminal_partition(level=4)
        print("counts for: FEMALES")
        print("... all counts except for total females exclude pregnancy...")
        utility.print_terminal_partition(level=5)
        print("count of total: " + str(count_female))
        print(
            "count premenopause sex hormone alteration: " + str(count_pre_alter)
        )
        print(
            "count perimenopause sex hormone alteration: " +
            str(count_peri_alter)
        )
        print(
            "count postmenopause sex hormone alteration: " +
            str(count_post_alter)
        )
        print(
            "count regular menstruation sex hormone alteration: " +
            str(count_regular_alter)
        )
        utility.print_terminal_partition(level=4)
        pass
    pass


def organize_female_menstruation_pregnancy_menopause_variables(
    table=None,
    report=None,
):
    """
    Organizes information about female persons' menopause, pregnancy,
    menstruation, oral contraception, and hormone therapy.

    These variables are specific to females and have null values for males.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "2724-0.0", "3591-0.0", "2834-0.0",
        "3140-0.0",
        "3700-0.0", "3710-0.0", "3720-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )

    ##########
    # Menstrual cycle

    # Determine count of days since person's last menstrual period.
    table["menstruation_days"] = table.apply(
        lambda row:
            determine_female_menstruation_days(
                sex_text=row["sex_text"],
                imputation_days=3, # imputation days to assign for current
                field_3700=row["3700-0.0"],
                field_3720=row["3720-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine count of days in person's usual menstrual cycle.
    table["menstruation_duration"] = table.apply(
        lambda row:
            determine_female_menstruation_cycle_duration(
                sex_text=row["sex_text"],
                field_3710=row["3710-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person reported irregularity in menstrual cycle.
    table["menstruation_irregularity"] = table.apply(
        lambda row:
            determine_female_menstruation_cycle_irregularity(
                sex_text=row["sex_text"],
                field_3710=row["3710-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether duration of regular menstrual cycle is in a normal
    # range.
    table["menstruation_regular_range"] = table.apply(
        lambda row:
            determine_female_menstruation_cycle_regular_duration_range(
                sex_text=row["sex_text"],
                menstruation_irregularity=row["menstruation_irregularity"],
                menstruation_duration=row["menstruation_duration"],
                menstruation_days=row["menstruation_days"],
                threshold_duration_low=21, # inclusive of lower threshold
                threshold_duration_high=36, # exclusive of higher threshold
            ),
        axis="columns", # apply function to each row
    )

    ##########
    # Parity and count of pregnancies and births.

    # Determine count of live or still births.
    table["births"] = table.apply(
        lambda row:
            determine_female_births_live_still_count(
                sex_text=row["sex_text"],
                field_2734=row["2734-0.0"], # count live births
                field_2774=row["2774-0.0"], # pregnancy loss
                field_3829=row["3829-0.0"], # count still births
            ),
        axis="columns", # apply function to each row
    )

    # Determine count of live or still births.
    table["births_any"] = table.apply(
        lambda row:
            determine_female_parity_births_any(
                sex_text=row["sex_text"],
                births=row["births"], # count live or still births
            ),
        axis="columns", # apply function to each row
    )

    # Determine whether female person experienced recent birth.
    table["age_last_live_birth"] = table.apply(
        lambda row:
            determine_female_age_at_most_recent_live_birth(
                sex_text=row["sex_text"],
                field_2754=row["2754-0.0"], # age at first live birth
                field_2764=row["2764-0.0"], # age at last live birth
                field_3872=row["3872-0.0"], # age at only live birth
            ),
        axis="columns", # apply function to each row
    )
    table["years_since_live_birth"] = table.apply(
        lambda row:
            determine_female_years_since_most_recent_live_birth(
                sex_text=row["sex_text"],
                age=row["age"],
                age_last_live_birth=row["age_last_live_birth"],
            ),
        axis="columns", # apply function to each row
    )
    table["birth_live_recent"] = table.apply(
        lambda row:
            determine_female_birth_live_recent(
                sex_text=row["sex_text"],
                births_any=row["births_any"],
                years_since_live_birth=row["years_since_live_birth"],
                threshold=6, # 5 or fewer years qualifies as recent
            ),
        axis="columns", # apply function to each row
    )

    # Determine count of pregnancies that ended early from abortive termination
    # or miscarriage.
    table["pregnancies_early"] = table.apply(
        lambda row:
            determine_female_pregnancies_early_count(
                sex_text=row["sex_text"],
                field_2774=row["2774-0.0"], # pregnancy loss
                field_3849=row["3849-0.0"], # count terminations
                field_3839=row["3839-0.0"], # count miscarriages
            ),
        axis="columns", # apply function to each row
    )

    # Determine count of total pregnancies.
    # These pregnancies might have ended in abortion, miscarriage, still birth,
    # or live birth.
    table["pregnancies"] = table.apply(
        lambda row:
            determine_female_pregnancies_total_count(
                sex_text=row["sex_text"],
                pregnancies_early=row["pregnancies_early"],
                births=row["births"],
            ),
        axis="columns", # apply function to each row
    )

    ##########
    # Menopause

    # Determine whether female persons experienced hysterectomy.
    table["hysterectomy"] = table.apply(
        lambda row:
            determine_female_hysterectomy(
                sex_text=row["sex_text"],
                field_2724=row["2724-0.0"],
                field_3591=row["3591-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether female persons experienced bilateral oophorectomy.
    table["oophorectomy"] = table.apply(
        lambda row:
            determine_female_oophorectomy(
                sex_text=row["sex_text"],
                field_2834=row["2834-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether female persons have experienced either hysterectomy or
    # oophorectomy.
    table["hysterectomy_or_oophorectomy"] = table.apply(
        lambda row:
            determine_female_hysterectomy_or_oophorectomy(
                sex_text=row["sex_text"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine whether female persons self reported natural menopause, which
    # neither involves hysterectomy nor oophorectomy.
    table["menopause_self_report"] = table.apply(
        lambda row:
            determine_female_menopause_self_report(
                sex_text=row["sex_text"],
                field_2724=row["2724-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine the menopause status of female persons.
    # This definition is ordinal.
    table["menopause_ordinal"] = table.apply(
        lambda row:
            determine_female_menopause_ordinal(
                sex_text=row["sex_text"],
                age=row["age"],
                menopause_self_report=row["menopause_self_report"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
                menstruation_regular_range=row["menstruation_regular_range"],
                threshold_age_low=47, # inclusive of lower threshold
                threshold_age_high=56, # exclusive of higher threshold
            ),
        axis="columns", # apply function to each row
    )

    table["menopause_age"] = table.apply(
        lambda row:
            determine_female_menopause_ordinal_age_only(
                sex_text=row["sex_text"],
                age=row["age"],
                menopause_self_report=row["menopause_self_report"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
                menstruation_regular_range=row["menstruation_regular_range"],
                threshold_age_low=47, # threshold age in years
                threshold_age_high=56, # threshold age in years
            ),
        axis="columns", # apply function to each row
    )
    table["menopause_age_self"] = table.apply(
        lambda row:
            determine_female_menopause_ordinal_age_self_only(
                sex_text=row["sex_text"],
                age=row["age"],
                menopause_self_report=row["menopause_self_report"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
                menstruation_regular_range=row["menstruation_regular_range"],
                threshold_age_low=47, # threshold age in years
                threshold_age_high=56, # threshold age in years
            ),
        axis="columns", # apply function to each row
    )
    table["menopause_age_self_oophorectomy"] = table.apply(
        lambda row:
            determine_female_menopause_ordinal_age_self_oophorectomy_only(
                sex_text=row["sex_text"],
                age=row["age"],
                menopause_self_report=row["menopause_self_report"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
                menstruation_regular_range=row["menstruation_regular_range"],
                threshold_age_low=47, # threshold age in years
                threshold_age_high=56, # threshold age in years
            ),
        axis="columns", # apply function to each row
    )

    ##########
    # Menstrual cycle

    # Determine female persons' categorical menstruation phase.
    table["menstruation_phase"] = table.apply(
        lambda row:
            determine_female_menstruation_phase(
                sex_text=row["sex_text"],
                menstruation_regular_range=row["menstruation_regular_range"],
                menstruation_days=row["menstruation_days"],
                threshold_ovulation=15, # threshold in days
            ),
        axis="columns", # apply function to each row
    )

    # Determine female persons' ordinal menstruation phase.
    table["menstruation_phase_early_late"] = table.apply(
        lambda row:
            determine_female_menstruation_phase_early_late(
                sex_text=row["sex_text"],
                menstruation_regular_range=row["menstruation_regular_range"],
                menstruation_days=row["menstruation_days"],
                threshold_middle_follicular=5, # threshold in days
                threshold_ovulation=15,
                threshold_middle_luteal=19, # threshold in days
            ),
        axis="columns", # apply function to each row
    )

    # Determine female persons' ordinal menstruation phase from cyclical
    # perspective.
    table["menstruation_phase_cycle"] = table.apply(
        lambda row:
            determine_female_menstruation_phase_cycle(
                sex_text=row["sex_text"],
                menstruation_phase_early_late=(
                    row["menstruation_phase_early_late"]
                ),
            ),
        axis="columns", # apply function to each row
    )

    # Determine whether female persons were pregnant.
    table["pregnancy"] = table.apply(
        lambda row:
            determine_female_pregnancy(
                sex_text=row["sex_text"],
                field_3140=row["3140-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    ##########
    # Ages at major events relevant to menstruation

    table = organize_female_ages_major_events_relevant_to_menstruation(
        table=table,
        report=report,
    )

    # Determine combination categories by menopause and hormone-atering therapy.
    # Report on 28 April 2021.
    # Category 1: 24158
    # Category 2: 178007
    # Category 3: 8438
    # Category 4: 51786
    #table = (
    #    utility.determine_binary_categorical_products_of_two_binary_variables(
    #        table=table,
    #        first="menopause_binary",
    #        second="hormone_therapy",
    #        prefix="menopause_hormone_category",
    #        report=report,
    #))

    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "2724-0.0", "3591-0.0", "2834-0.0",
            "3140-0.0",
            "3700-0.0", "3710-0.0", "3720-0.0",
        ],
        axis="columns",
        inplace=True
    )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "sex_text",
        "age",
        "menstruation_days",
        "menstruation_duration", "menstruation_irregularity",
        "menstruation_regular_range",
        "hysterectomy", "oophorectomy", "hysterectomy_or_oophorectomy",
        "menopause_ordinal",
        "menstruation_phase",
        "menstruation_phase_early_late", "menstruation_phase_cycle",
        "pregnancy",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    table_report.sort_values(
        by=["pregnancy", "menopause_ordinal"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print(
            "report: " +
            "organize_female_menstruation_pregnancy_menopause_variables()"
        )
        print(table_report)
        report_female_pregnancies_births(
            table=table,
        )
        report_female_menstruation_regularity_duration_range(
            threshold_duration_low=21, # inclusive of lower threshold
            threshold_duration_high=36, # exclusive of higher threshold
            table=table,
        )
        report_female_self_report_menopause_hysterectomy_oophorectomy(
            table=table,
        )
        report_female_menopause_three_groups(
            female_menopause="menopause_ordinal",
            table=table,
        )
        report_female_menopause_three_groups(
            female_menopause="menopause_age",
            table=table,
        )
        report_female_menopause_three_groups(
            female_menopause="menopause_age_self",
            table=table,
        )
        report_female_menopause_three_groups(
            female_menopause="menopause_age_self_oophorectomy",
            table=table,
        )
        report_combination_female_regular_menstruation_by_menopause(
            menopause="menopause_ordinal",
            menstruation_regular_range="menstruation_regular_range",
            table=table,
        )
        report_alteration_sex_hormones_by_female_menopause(
            menopause="menopause_ordinal",
            menstruation_regular_range="menstruation_regular_range",
            alteration_sex_hormone="alteration_sex_hormone",
            table=table,
        )
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


#####
# SCRAP
####


# review: TCW on _____
def obsolete_determine_female_menopause_ordinal_before_2022_01_18(
    sex_text=None,
    age=None,
    threshold_age_pre=None,
    threshold_age_post=None,
    menstruation_duration=None,
    menstruation_irregularity=None,
    menstruation_days=None,
    threshold_menstruation_pre=None,
    threshold_menstruation_post=None,
    menopause_self_report=None,
    hysterectomy=None,
    oophorectomy=None,
):
    """
    Determine whether female persons qualify for premenopause, perimenopause,
    or postmenopause categories.
    This definition's encoding is ordinal.

    This definition is permissive of missing values in some variables.

    Pay close attention to "and", "or" logic in comparisons.

    priority:
    1. self report of natural menopause; oophorectomy
    2. menstruation cycle irregularity; menstruation cycle duration; days since
    previous menstruation
    3. age

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        threshold_age_pre (int): threshold age in years, below which to consider
            females premenopausal
        threshold_age_post (int): threshold age in years, above which to
            consider all females postmenopausal
        menstruation_duration (int): count of days in person's usual
            menstrual cycle
        menstruation_irregularity (int): binary representation of whether
            person reported irregularity in menstrual cycle
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
        threshold_menstruation_pre (int): threshold in days since last
            menstrual period, below which to consider females premenopausal
        threshold_menstruation_post (int): threshold in days since last
            menstrual period, above which to consider all females postmenopausal
        menopause_self_report (float): binary representation of whether person
            self reported having experienced menopause
        hysterectomy (float): binary logical representation of whether person
            has had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            has had a bilateral oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Determine categories for menopause.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine postmenopause.
        if (
            (
                (not pandas.isna(menopause_self_report)) and
                (menopause_self_report == 1)
            ) or
            (
                (not pandas.isna(oophorectomy)) and
                (oophorectomy == 1)
            ) or
            (
                (not pandas.isna(menstruation_duration)) and
                (menstruation_duration >= threshold_menstruation_post)
            ) or
            (
                (not pandas.isna(menstruation_days)) and
                (menstruation_days >= threshold_menstruation_post)
            ) or
            (
                (not pandas.isna(age)) and
                (age >= threshold_age_post)
            )
        ):
            # Person qualifies for postmenopause.
            # Any valid variable that satisfies conditions can justify the
            # postmenopause definition.
            value = 2
        # Determine perimenopause.
        # Notice that "menstruation_days" suffices to qualify for perimenopause
        # category if value is greater than threshold for premenopause category.
        # Notice that "menstruation_irregularity" is either "1" or missing.
        elif (
            (
                (
                    (pandas.isna(menopause_self_report)) or
                    (menopause_self_report == 0)
                ) and
                (
                    (pandas.isna(oophorectomy)) or
                    (oophorectomy == 0)
                )
            ) and
            (
                (
                    (not pandas.isna(menstruation_irregularity)) and
                    (menstruation_irregularity == 1)
                ) or
                (
                    (not pandas.isna(menstruation_duration)) and
                    (menstruation_duration >= threshold_menstruation_pre) and
                    (menstruation_duration < threshold_menstruation_post)
                ) or
                (
                    (not pandas.isna(menstruation_days)) and
                    (menstruation_days >= threshold_menstruation_pre) and
                    (menstruation_days < threshold_menstruation_post)
                ) or
                (
                    (not pandas.isna(age)) and
                    (age >= threshold_age_pre) and
                    (age < threshold_age_post)
                )
            )
        ):
            # Person qualitifies for perimenopause.
            value = 1
        # Determine premenopause.
        # Notice that "menstruation_days" does not suffice to qualify for
        # premenopause category because even persons with longer cycles can have
        # values less than threshold.
        elif (
            (
                (
                    (pandas.isna(menopause_self_report)) or
                    (menopause_self_report == 0)
                ) and
                (
                    (pandas.isna(oophorectomy)) or
                    (oophorectomy == 0)
                )
            ) and
            (
                (
                    (not pandas.isna(menstruation_duration)) and
                    (menstruation_duration < threshold_menstruation_pre)
                ) or
                (
                    (not pandas.isna(age)) and
                    (age < threshold_age_pre)
                )
            )
        ):
            # Person qualifies for premenopause.
            value = 0
        else:
            # Person does not qualify for any categories.
            value = float("nan")
    else:
        # Menopause undefined for males.
        value = float("nan")
    # Return information.
    return value


# TODO: obsolete? consider removal?
def obsolete_determine_female_menopause_binary(
    sex_text=None,
    age=None,
    threshold_age=None,
    menstruation_days=None,
    threshold_menstruation_days=None,
    field_2724=None,
    hysterectomy=None,
    oophorectomy=None,
):
    """
    Determine whether female persons qualify for premenopause or postmenopause
    categories.
    This definition's encoding is binary.

    This definition is permissive of missing values in some variables.

    Pay close attention to "and", "or" logic in comparisons.

    priority:
    1. (self report of menopause) or (oophorectomy)
    2. (days since previous menstruation)
    3. (age)

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        threshold_age (int): threshold age in years, above which to
            consider all females postmenopausal
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
        threshold_menstruation_days (int): threshold in days since last
            menstrual period, above which to consider females postmenopausal
        field_2724 (float): UK Biobank field 2724, whether person has
            experienced menopause
        hysterectomy (float): binary logical representation of whether person
            has had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            has had a bilateral oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret menopause.
    # Interpret whether person self-reported menopause that did not involve
    # hysterectomy.
    menopause_natural = interpret_menopause_self_report(
        field_2724=field_2724,
    )
    # Determine categories for menopause.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine postmenopause.
        if (
            (
                (not pandas.isna(menopause_natural)) and
                (menopause_natural == 1)
            ) or
            (
                (not pandas.isna(oophorectomy)) and
                (oophorectomy == 1)
            ) or
            (
                (not pandas.isna(menstruation_days)) and
                (menstruation_days >= threshold_menstruation_days)
            ) or
            (
                (not pandas.isna(age)) and
                (age >= threshold_age)
            )
        ):
            # Person qualifies for postmenopause.
            # Any valid variable that satisfies conditions can justify the
            # postmenopause definition.
            value = 1
        # Determine premenopause.
        elif (
            (
                (
                    (pandas.isna(menopause_natural)) or
                    (menopause_natural == 0)
                ) and
                (
                    (pandas.isna(oophorectomy)) or
                    (oophorectomy == 0)
                )
            ) and
            (
                (
                    (not pandas.isna(menstruation_days)) and
                    (menstruation_days < threshold_menstruation_days)
                ) or
                (
                    (not pandas.isna(age)) and
                    (age < threshold_age)
                )
            )
        ):
            # Person qualifies for premenopause.
            value = 0
        else:
            # Person does not qualify for any categories.
            value = float("nan")
    else:
        # Menopause undefined for males.
        value = float("nan")
    # Return information.
    return value


# TODO: obsolete? consider removal?
def obsolete_determine_female_menopause_binary_strict(
    sex_text=None,
    age=None,
    threshold_age=None,
    menstruation_days=None,
    threshold_menstruation_days=None,
    field_2724=None,
    hysterectomy=None,
    oophorectomy=None,
):
    """
    Determine whether female persons qualify for premenopause or postmenopause
    categories.
    This definition's encoding is binary.

    This definition is permissive of null or missing values in self-report
    menopause ("menopause_natural" from field 2724) and in days since previous
    menstruation ("menstruation_days" from field 3700). Females who qualify for
    perimenopause might be more likely to reply "unsure" in response to the
    self-report menopause variable, giving it a null or missing value. Also,
    females who qualify for postmenopause might be more likely to have a null or
    missing value in the variable for days since previous menstruation.

    In contrast, this definition is not permissive of null or missing values in
    variables for oophorectomy or age.

    Pay close attention to "and", "or" logic in comparisons.

    priority:
    1. (self report of menopause) or (oophorectomy)
    2. (days since previous menstruation)
    3. (age)

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        threshold_age (int): threshold age in years, above which to
            consider all females postmenopausal
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
        threshold_menstruation_days (int): threshold in days since last
            menstrual period, above which to consider females postmenopausal
        field_2724 (float): UK Biobank field 2724, whether person has
            experienced menopause
        hysterectomy (float): binary logical representation of whether person
            has had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            has had a bilateral oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret menopause.
    # Interpret whether person self-reported menopause that did not involve
    # hysterectomy.
    menopause_natural = interpret_menopause_self_report(
        field_2724=field_2724,
    )
    # Determine categories for menopause.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(oophorectomy)) and
            (not pandas.isna(age))
        ):
            # Determine postmenopause.
            if (
                (
                    (not pandas.isna(menopause_natural)) and
                    (menopause_natural == 1)
                ) or
                (oophorectomy == 1) or
                (
                    (not pandas.isna(menstruation_days)) and
                    (menstruation_days >= threshold_menstruation_days)
                ) or
                (age >= threshold_age)
            ):
                # Person qualifies for postmenopause.
                value = 1
            elif (
                (
                    (
                        (pandas.isna(menopause_natural)) or
                        (menopause_natural == 0)
                    ) and
                    (
                        (oophorectomy == 0)
                    )
                ) and
                (
                    (
                        (not pandas.isna(menstruation_days)) and
                        (menstruation_days < threshold_menstruation_days)
                    ) or
                    (
                        (age < threshold_age)
                    )
                )
            ):
                # Person qualifies for premenopause.
                # In order even to consider this premenopause definition,
                # variable "menopause_natural" must be either missing or false.
                value = 0
            else:
                # Persons does not qualify for any categories.
                value = float("nan")
        else:
            # Variables have missing values.
            value = float("nan")
    else:
        # Menopause undefined for males.
        value = float("nan")
    # Return information.
    return value


# TODO: obsolete? consider removal?
def obsolete_determine_female_menopause_ordinal_strict(
    sex_text=None,
    age=None,
    threshold_age_pre=None,
    threshold_age_post=None,
    menstruation_days=None,
    threshold_menstruation_days_pre=None,
    threshold_menstruation_days_post=None,
    field_2724=None,
    hysterectomy=None,
    oophorectomy=None,
):
    """
    Determine whether female persons qualify for premenopause, perimenopause,
    or postmenopause categories.
    This definition's encoding is ordinal.

    This definition is not permissive of missing values in any variables, with
    the exception of self-reported menopause (data field 2724). Females who
    qualify for perimenopause might be likely to reply "unsure" in response to
    this variable, giving it a missing or uncertain value.

    Pay close attention to "and", "or" logic in comparisons.

    priority:
    1. (self report of menopause) or (oophorectomy)
    2. (days since previous menstruation)
    3. (age)

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        threshold_age_pre (int): threshold age in years, below which to consider
            females premenopausal
        threshold_age_post (int): threshold age in years, above which to
            consider all females postmenopausal
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
        threshold_menstruation_days_pre (int): threshold in days since last
            menstrual period, below which to consider females premenopausal
        threshold_menstruation_days_post (int): threshold in days since last
            menstrual period, above which to consider all females postmenopausal
        field_2724 (float): UK Biobank field 2724, whether person has
            experienced menopause
        hysterectomy (float): binary logical representation of whether person
            has had an hysterectomy
        oophorectomy (float): binary logical representation of whether person
            has had a bilateral oophorectomy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret menopause.
    # Interpret whether person self-reported menopause that did not involve
    # hysterectomy.
    menopause_natural = interpret_menopause_self_report(
        field_2724=field_2724,
    )
    # Assign aliases for variables with long names.
    menstruation_days_pre = threshold_menstruation_days_pre
    menstruation_days_post = threshold_menstruation_days_post
    # Determine categories for menopause.
    # Pay close attention to "and", "or" logic in comparisons.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(oophorectomy)) and
            (not pandas.isna(age))
        ):
            # Determine postmenopause.
            if (
                (
                    (not pandas.isna(menopause_natural)) and
                    (menopause_natural == 1)
                ) or
                (oophorectomy == 1) or
                (
                    (not pandas.isna(menstruation_days)) and
                    (menstruation_days >= menstruation_days_post)
                ) or
                (age >= threshold_age_post)
            ):
                # Person qualifies for postmenopause.
                value = 2
            # Determine perimenopause
            elif (
                (
                    (
                        (pandas.isna(menopause_natural)) or
                        (menopause_natural == 0)
                    ) and
                    (
                        (oophorectomy == 0)
                    )
                ) and
                (
                    (
                        (not pandas.isna(menstruation_days)) and
                        (menstruation_days >= menstruation_days_pre) and
                        (menstruation_days < menstruation_days_post)
                    ) or
                    (
                        (age >= threshold_age_pre) and
                        (age < threshold_age_post)
                    )
                )
            ):
                # Person qualitifies for perimenopause.
                value = 1
            # Determine premenopause.
            elif (
                (
                    (
                        (pandas.isna(menopause_natural)) or
                        (menopause_natural == 0)
                    ) and
                    (
                        (oophorectomy == 0)
                    )
                ) and
                (
                    (
                        (not pandas.isna(menstruation_days)) and
                        (menstruation_days < menstruation_days_pre)
                    ) or
                    (age < threshold_age_pre)
                )
            ):
                # Person qualifies for premenopause.
                value = 0
            else:
                # Persons does not qualify for any categories.
                value = float("nan")
        else:
            # Variables have missing values.
            value = float("nan")
    else:
        # Menopause undefined for males.
        value = float("nan")
    # Return information.
    return value


# TODO: TCW 20 January 2022 --> I think I want to remove this function... obsolete...
# review: TCW on 18 January 2022
def obsolete_determine_female_menstruation_days_threshold(
    sex_text=None,
    menopause_ordinal=None,
    menstruation_days=None,
    threshold_premenopause=None,
    threshold_perimenopause=None,
):
    """
    Determines a female person's days of the menstrual cycle by applying a
    threshold filter that is specific to menopause status.

    arguments:
        sex_text (str): textual representation of sex selection
        menopause_ordinal (float): ordinal representation of female person's
            menopause status (0: pre-; 1: peri-; 2: post-)
        menstruation_days (float): count of days since previous menstruation
            (menstrual period, start of menstrual cycle)
        threshold_premenopause (int): threshold in days for menstrual cycle in
            premenopausal female persons, above which to set menstruation days
            to missing
        threshold_perimenopause (int): threshold in days for menstrual cycle in
            perimenopausal female persons, above which to set menstruation days
            to missing

    raises:

    returns:
        (float): interpretation value

    """

    # Determine categorical phase of menstruation cycle.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(menopause_ordinal)) and
            (not pandas.isna(menstruation_days))
        ):
            # Determine threshold values for premenopausal females.
            if (menopause_ordinal == 0):
                if (menstruation_days < threshold_premenopause):
                    # Days of menstrual cycle are below threshold.
                    value = menstruation_days
                elif (menstruation_days >= threshold_premenopause):
                    # Days of menstrual cycle are above threshold.
                    value = float("nan")
                else:
                    # Confusion in interpretation.
                    value = float("nan")
            # Determine threshold values for perimenopausal females.
            elif (menopause_ordinal == 1):
                if (menstruation_days < threshold_perimenopause):
                    # Days of menstrual cycle are below threshold.
                    value = menstruation_days
                elif (menstruation_days >= threshold_perimenopause):
                    # Days of menstrual cycle are above threshold.
                    value = float("nan")
                else:
                    # Confusion in interpretation.
                    value = float("nan")
            # Determine threshold values for postmenopausal females.
            elif (menopause_ordinal == 2):
                # Menstruation undefined for postmenopause females.
                value = float("nan")
            else:
                # Person does not qualify for any categories.
                value = float("nan")
        else:
            # Necessary variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value



##########
# Psychological variables


# TODO: TCW 20 October 2021
# TODO: 1. I need to "import and interpret" medication class variables from Brandon J. Coombes
# TODO: 2. I need to define categorical dummy variables for the medication classes
# TODO: 3. I need to test those categorical dummy variables in phenotypic regressions

#    table = create_reduce_categorical_variable_indicators(
#        table=table,
#        index="eid",
#        column="assessment_site",
#        prefix="site",
#        separator="_",
#        report=True,
#    )

#    pail_indicator = create_categorical_variable_indicators(
#        table=table,
#        index=index,
#        column=column,
#        prefix=str(prefix + "-" + "indicator"),
#        separator=separator,
#        report=report,
#    )


def interpret_neuroticism(
    field_20127=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 20127.

    Accommodate inexact float values.

    Note:

    arguments:
        field_20127 (float): UK Biobank field 20127, neuroticism

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_20127)) and
        (0 <= float(field_20127) and float(field_20127) < 13)
    ):
        # The variable has a valid value.
        # Interpret the value.
        interpretation = float(field_20127)
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def determine_neuroticism(
    field_20127=None,
):
    """
    Determine neuroticism score.

    arguments:
        field_20127 (float): UK Biobank field 20127, neuroticism

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret codes.
    # Set value.
    value = interpret_neuroticism(
        field_20127=field_20127,
    )
    # Return information.
    return value


def translate_import_binary_to_boolean(
    import_binary=None,
):
    """
    Translate binary import variable to Boolean.

    arguments:
        import_boolean (str or bool): import Boolean variable

    raises:

    returns:
        (str): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(import_binary)) and
        (str(import_binary) != "<NA>") and
        (-0.5 <= float(import_binary) and float(import_binary) < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= float(import_binary) and float(import_binary) < 0.5):
            # 0: False
            value = "FALSE"
        elif (0.5 <= float(import_binary) and float(import_binary) < 1.5):
            # 1: True
            value = "TRUE"
        else:
            # uninterpretable
            value = "<NA>"
    else:
        # null
        value = "<NA>"
    # Return.
    return value


def translate_import_boolean_to_binary(
    import_boolean=None,
):
    """
    Translate Boolean import variable to binary.

    arguments:
        import_boolean (str or bool): import Boolean variable

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(import_boolean)) and
        (str(import_boolean) != "<NA>") and
        ((str(import_boolean) == "TRUE") or (str(import_boolean) == "FALSE"))
    ):
        # The variable has a valid value.
        if (str(import_boolean) == "FALSE"):
            # 0: "control"
            value = 0
        elif (str(import_boolean) == "TRUE"):
            # 1: "case"
            value = 1
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def interpret_import_binary_disorder_case(
    import_binary_disorder=None,
):
    """
    Intepret import binary variable for case definition.

    1: case
    0: not case

    Accommodate inexact float values.

    arguments:
        import_binary_disorder (float): import binary variable case definition

    raises:

    returns:
        (bool): interpretation value

    """

    #import_bipolar_disorder = float(import_bipolar_disorder)

    # Interpret field code.
    if (
        (not pandas.isna(import_binary_disorder)) and
        (-0.5 <= import_binary_disorder and import_binary_disorder < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= import_binary_disorder and import_binary_disorder < 0.5):
            # 0: "control"
            value = 0
        elif (0.5 <= import_binary_disorder and import_binary_disorder < 1.5):
            # 1: "case"
            value = 1
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def import_disorder_case_definitions(
    columns_import=None,
    table=None,
):
    """
    Organizes imputation of missing measurements for hormones.

    arguments:
        columns_import (list<dict>): names of import columns for Boolean
            disorder case definitions
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Translate import Boolean columns to binary.
    for column_import in columns_import:
        # Organize variable names.
        name_binary = str(column_import + "_binary")
        name_read = str(column_import + "_read")
        # Convert Boolean variables to binary.
        table[name_binary] = table.apply(
            lambda row:
                translate_import_boolean_to_binary(
                    import_boolean=row[column_import],
                ),
            axis="columns", # apply function to each row
        )
        # Interpret case definition.
        table[name_read] = table.apply(
            lambda row:
                interpret_import_binary_disorder_case(
                    import_binary_disorder=row[name_binary],
                ),
            axis="columns", # apply function to each row
        )
        pass
    # Return information.
    return table


def determine_loose_strict_disorder_control(
    case_loose=None,
    case_strict=None,
):
    """
    Determine controls who are cases neither by loose nor strict definitions.

    1: control
    0: not control

    arguments:
        case_loose (float): loose case definition
        case_strict (float): strict case definition

    raises:

    returns:
        (bool): interpretation value

    """

    # Interpret field code.
    if (
        (
            (not pandas.isna(case_loose)) and
            ((case_loose == 0) or (case_loose == 1))
        ) and
        (
            (not pandas.isna(case_strict)) and
            ((case_strict == 0) or (case_strict == 1))
        )
    ):
        # Require both loose and strict definitions to be non-missing and have
        # interpretable values.
        # The variable has a valid value.
        if (
            (case_loose == 0) and
            (case_strict == 0)
        ):
            # Both definitions designate as "not case".
            # Control.
            value = 1
        elif (
            (case_loose == 1) or
            (case_strict == 1)
        ):
            # One or two definitions designate as "case".
            # Not control.
            value = 0
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def organize_psychotropic_drug_class_categories(
    drug_classes=None,
    table=None,
    report=None,
):
    """
    Organizes categories for classes of psychotropic drugs (psychoactive
    medications).

    arguments:
        drug_classes (dict<dict<string>>): translations for names of drug
            classes
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Impute missing values.
    hormones = list()
    for drug_class_import in drug_classes.keys():
        drug_class = drug_classes[drug_class_import]
        # Interpret drug class encoding and translate name.
        table[drug_class] = table.apply(
            lambda row:
                interpret_import_drug_class(
                    code=row[drug_class_import],
                ),
            axis="columns", # apply function to each row
        )
    # Return information.
    return table


def determine_disorder_control_case_combination(
    control=None,
    case=None,
):
    """
    Determines whether person qualifies as a control or case for a disorder.

    NA: missing value indicates neither control nor case
    0: control
    1: case

    arguments:
        control (float): binary logical representation of control status
        case (float): binary logical representation of case status

    raises:

    returns:
        (float): binary logical representation of Bipolar Disorder control or
            case status

    """

    # Integrate information from all criteria.
    if (
        (
            (not pandas.isna(control)) and
            (control == 1)
        ) and
        (
            (pandas.isna(case)) or
            (case == 0)
        )
    ):
        # Control.
        value = 0
    elif (
        (
            (not pandas.isna(case)) and
            (case == 1)
        ) and
        (
            (pandas.isna(control)) or
            (control == 0)
        )
    ):
        # Case.
        value = 1
    else:
        # Missing information.
        value = float("nan")
    # Return information.
    return value


def organize_psychology_variables(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about general attributes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy data.
    table = table.copy(deep=True)

    # Read reference table of interpretations of variable codes.
    source = read_source_fields_codes_interpretations(
        path_dock=path_dock,
    )

    # Neuroticism score.
    table["neuroticism"] = table.apply(
        lambda row:
            determine_neuroticism(
                field_20127=row["20127-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=["neuroticism"],
        table=table,
    )

    # Import case definitions from Brandon J. Coombes for Bipolar Disorder and
    # Major Depressive Disorder.
    # Convert format of variable so that it is consistent with the others.
    table["import_bipolar.cc_boolean"] = table.apply(
        lambda row:
            translate_import_binary_to_boolean(
                import_binary=row["import_bipolar.cc"],
            ),
        axis="columns", # apply function to each row
    )
    columns_import = [
        "import_broad_depression", "import_probable_mdd",
        "import_bipolar", "import_icd_bipolar",
        "import_bipolar.cc_boolean", "import_icd_bipolar.cc",
    ]
    table = import_disorder_case_definitions(
        columns_import=columns_import,
        table=table,
    )
    # Translate column names.
    translations = dict()
    translations["import_broad_depression_read"] = "depression_case_loose"
    translations["import_probable_mdd_read"] = "depression_case_strict"
    ###translations["import_bipolar_read"] = "bipolar_case_loose"
    ###translations["import_icd_bipolar_read"] = "bipolar_case_strict"
    translations["import_bipolar.cc_boolean_read"] = "bipolar_case_loose"
    translations["import_icd_bipolar.cc_read"] = "bipolar_case_strict"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Determine controls for Major Depressive Disorder and Bipolar Disorder.
    table["depression_control"] = table.apply(
        lambda row:
            determine_loose_strict_disorder_control(
                case_loose=row["depression_case_loose"],
                case_strict=row["depression_case_strict"],
            ),
        axis="columns", # apply function to each row
    )
    table["bipolar_control"] = table.apply(
        lambda row:
            determine_loose_strict_disorder_control(
                case_loose=row["bipolar_case_loose"],
                case_strict=row["bipolar_case_strict"],
            ),
        axis="columns", # apply function to each row
    )


    # Organize case-control variables.
    table["depression_control_case_loose"] = table.apply(
        lambda row:
            determine_disorder_control_case_combination(
                control=row["depression_control"],
                case=row["depression_case_loose"],
            ),
        axis="columns", # apply function to each row
    )
    table["depression_control_case_strict"] = table.apply(
        lambda row:
            determine_disorder_control_case_combination(
                control=row["depression_control"],
                case=row["depression_case_strict"],
            ),
        axis="columns", # apply function to each row
    )
    table["bipolar_control_case_loose"] = table.apply(
        lambda row:
            determine_disorder_control_case_combination(
                control=row["bipolar_control"],
                case=row["bipolar_case_loose"],
            ),
        axis="columns", # apply function to each row
    )
    table["bipolar_control_case_strict"] = table.apply(
        lambda row:
            determine_disorder_control_case_combination(
                control=row["bipolar_control"],
                case=row["bipolar_case_strict"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine categorical classes of psychotropic drugs (medications).
    # Define names of classes of psychotropic drugs as defined by Brandon J.
    # Coombes.
    if False:
        table = organize_psychotropic_drug_class_categories(
            drug_classes=source["drug_class"],
            table=table,
            report=report,
        )

    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            #"20127-0.0",
            "import_broad_depression", "import_probable_mdd",
            "import_bipolar", "import_icd_bipolar"
        ],
        axis="columns",
        inplace=True
    )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "sex_text", "age", "neuroticism", "neuroticism_log",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: organize_psychology_variables()")

        count_depression_control_loose = int(
            table.loc[
                (
                    (table["depression_control_case_loose"] == 0)
                ), :
            ].shape[0]
        )
        count_depression_case_loose = int(
            table.loc[
                (
                    (table["depression_control_case_loose"] == 1)
                ), :
            ].shape[0]
        )
        count_depression_control_strict = int(
            table.loc[
                (
                    (table["depression_control_case_strict"] == 0)
                ), :
            ].shape[0]
        )
        count_depression_case_strict = int(
            table.loc[
                (
                    (table["depression_control_case_strict"] == 1)
                ), :
            ].shape[0]
        )
        count_bipolar_control_loose = int(
            table.loc[
                (
                    (table["bipolar_control_case_loose"] == 0)
                ), :
            ].shape[0]
        )
        count_bipolar_case_loose = int(
            table.loc[
                (
                    (table["bipolar_control_case_loose"] == 1)
                ), :
            ].shape[0]
        )
        count_bipolar_control_strict = int(
            table.loc[
                (
                    (table["bipolar_control_case_strict"] == 0)
                ), :
            ].shape[0]
        )
        count_bipolar_case_strict = int(
            table.loc[
                (
                    (table["bipolar_control_case_strict"] == 1)
                ), :
            ].shape[0]
        )
        # Categorical ancestry and ethnicity.
        utility.print_terminal_partition(level=2)
        print("Major Depressive Disorder cases and controls")
        utility.print_terminal_partition(level=5)
        print("Loose definition")
        print("controls: " + str(count_depression_control_loose))
        print("cases: " + str(count_depression_case_loose))
        utility.print_terminal_partition(level=5)
        print("Strict definition")
        print("controls: " + str(count_depression_control_strict))
        print("cases: " + str(count_depression_case_strict))
        utility.print_terminal_partition(level=2)
        print("Bipolar Disorder cases and controls")
        utility.print_terminal_partition(level=5)
        print("Loose definition")
        print("controls: " + str(count_bipolar_control_loose))
        print("cases: " + str(count_bipolar_case_loose))
        utility.print_terminal_partition(level=5)
        print("Strict definition")
        print("controls: " + str(count_bipolar_control_strict))
        print("cases: " + str(count_bipolar_case_strict))
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Depression


# ICD9 and ICD10 codes...
# Bespoke?
# Patient Health Questionnaire 9-Questions (PHQ-9)
# Composite International Diagnostic Interview - Short Form (CIDI-SF)


def organize_depression_variables(
    table=None,
    report=None,
):
    """
    Organizes information about general attributes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy data.
    table = table.copy(deep=True)


    # TODO: interpret ICD9 / ICD10 codes for depression...




    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    if False:
        table_clean.drop(
            labels=[
                "20127-0.0",
            ],
            axis="columns",
            inplace=True
        )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "sex_text", "age", "neuroticism", "neuroticism_log",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: organize_neuroticism_variables()")
        utility.print_terminal_partition(level=3)
        print("translations of general attribute column names...")
        for old in translations.keys():
            print("   " + old + ": " + translations[old])
        utility.print_terminal_partition(level=3)
        print(table_report)
        utility.print_terminal_partition(level=3)
        # Variable types.
        utility.print_terminal_partition(level=2)
        print("After type conversion")
        print(table_report.dtypes)
        utility.print_terminal_partition(level=3)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Psychiatric disorder diagnoses




##########
# Alcohol consumption


# review: TCW, 16 March 2022
def interpret_alcohol_consumption_frequency(
    field_1558=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 1558.

    Notice that the directionality of the UK Biobank's data-coding for
    data-field 1558 is the inverse of logical directionality.

    Data-Field "1558": "Alcohol intake frequency"
    UK Biobank data-coding "100402" for data-field "1558".
    1: "Daily or almost daily"
    2: "Three or four times a week"
    3: "Once or twice a week"
    4: "One to three times a month"
    5: "Special occasions only"
    6: "Never"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_1558 (float): UK Biobank field 1558, frequency of alcohol
            consumption

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_1558)) and
        (-3.5 <= field_1558 and field_1558 < 6.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (5.5 <= field_1558 and field_1558 < 6.5):
            # 6: "Never"
            interpretation = 0
        elif (4.5 <= field_1558 and field_1558 < 5.5):
            # 5: "Special occasions only"
            interpretation = 1
        elif (3.5 <= field_1558 and field_1558 < 4.5):
            # 4: "One to three times a month"
            interpretation = 2
        elif (2.5 <= field_1558 and field_1558 < 3.5):
            # 3: "Once or twice a week"
            interpretation = 3
        elif (1.5 <= field_1558 and field_1558 < 2.5):
            # 2: "Three or four times a week"
            interpretation = 4
        elif (0.5 <= field_1558 and field_1558 < 1.5):
            # 1: "Daily or almost daily"
            interpretation = 5
        elif (-3.5 <= field_1558 and field_1558 < -2.5):
            # -3: "Prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def interpret_alcohol_consumption_status(
    field_20117=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 20117.

    Data-Field "20117": "Alcohol drinker status"
    UK Biobank data-coding "90" for data-field "20117".
    0: "Never"
    1: "Previous"
    2: "Current"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_20117 (float): UK Biobank field 20117, current alcohol consumption
            status

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_20117)) and
        (-3.5 <= field_20117 and field_20117 < 2.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (-0.5 <= field_20117 and field_20117 < 0.5):
            # 0: "Never"
            interpretation = 0
        elif (0.5 <= field_20117 and field_20117 < 1.5):
            # 1: "Previous"
            interpretation = 1
        elif (1.5 <= field_20117 and field_20117 < 2.5):
            # 2: "Current"
            interpretation = 2
        elif (-3.5 <= field_20117 and field_20117 < -2.5):
            # -3: "Prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def interpret_former_alcohol_consumption(
    field_3731=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 3731.

    Data-Field "3731": "Former alcohol drinker"
    UK Biobank data-coding "100352" for data-field "3731".
    1: "Yes"
    0: "No"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_3731 (float): UK Biobank field 3731, former alcohol consumption

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_3731)) and
        (-3.5 <= field_3731 and field_3731 < 1.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (-0.5 <= field_3731 and field_3731 < 0.5):
            # 0: "No"
            interpretation = 0
        elif (0.5 <= field_3731 and field_3731 < 1.5):
            # 1: "Yes"
            interpretation = 1
        elif (-3.5 <= field_3731 and field_3731 < -2.5):
            # -3: "Prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def interpret_alcohol_consumption_trend(
    field_1628=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 1628.

    Data-Field "1628": "Alcohol intake versus 10 years previously"
    UK Biobank data-coding "100417" for data-field "1628".
    1: "More nowadays"
    2: "About the same"
    3: "Less nowadays"
    -1: "Do not know"
    -3: "Prefer not to answer"

    Accommodate inexact float values.

    arguments:
        field_1628 (float): UK Biobank field 1628, trend in alcohol consumption
            trend in alcohol consumption current versus previous

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_1628)) and
        (-3.5 <= field_1628 and field_1628 < 3.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (0.5 <= field_1628 and field_1628 < 1.5):
            # 1: "More nowadays"
            # "More nowadays" means previous was less than current.
            interpretation = 0
        elif (1.5 <= field_1628 and field_1628 < 2.5):
            # 2: "About the same"
            interpretation = 1
        elif (2.5 <= field_1628 and field_1628 < 3.5):
            # 3: "Less nowadays"
            # "Less nowadays" means previous was more than current.
            interpretation = 2
        elif (-1.5 <= field_1628 and field_1628 < -0.5):
            # -1: "Do not know"
            interpretation = float("nan")
        elif (-3.5 <= field_1628 and field_1628 < -2.5):
            # -3: "Prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def determine_alcohol_consumption_frequency(
    field_1558=None,
):
    """
    Determine frequency of alcohol consumption.

    # 0: "Never"
    # 1: "Special occasions only"
    # 2: "One to three times a month"
    # 3: "Once or twice a week"
    # 4: "Three or four times a week"
    # 5: "Daily or almost daily"

    arguments:
        field_1558 (float): UK Biobank field 1558, frequency of alcohol
            consumption

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret data field coding.
    interpretation = interpret_alcohol_consumption_frequency(
        field_1558=field_1558,
    )
    # Apply any relevant logic.
    value = interpretation
    # Return information.
    return value


def determine_alcohol_consumption_status(
    field_20117=None,
):
    """
    Determine frequency of alcohol consumption.

    arguments:
        field_20117 (float): UK Biobank field 20117, current alcohol consumption
            status

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret data field coding.
    interpretation = interpret_alcohol_consumption_status(
        field_20117=field_20117,
    )
    # Apply any relevant logic.
    value = interpretation
    # Return information.
    return value


def determine_former_alcohol_consumption(
    field_3731=None,
):
    """
    Determine former alcohol consumption.

    arguments:
        field_3731 (float): UK Biobank field 3731, former alcohol consumption

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret data field coding.
    interpretation = interpret_former_alcohol_consumption(
        field_3731=field_3731,
    )
    # Apply any relevant logic.
    value = interpretation
    # Return information.
    return value


def determine_alcohol_consumption_trend(
    field_1628=None,
):
    """
    Determine frequency of alcohol consumption.

    arguments:
        field_1628 (float): UK Biobank field 1628, comparison of current alcohol
            consumption to previous alcohol consumption

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret data field coding.
    interpretation = interpret_alcohol_consumption_trend(
        field_1628=field_1628,
    )
    # Apply any relevant logic.
    value = interpretation
    # Return information.
    return value


def determine_alcohol_consumption_current(
    alcohol_frequency=None,
    alcohol_status=None,
    alcohol_trend=None,
):
    """
    Determine current alcohol consumption status.

    arguments:
        alcohol_frequency (float): frequency of alcohol consumption
        alcohol_status (float): alcohol consumption status
        alcohol_trend (float): ordinal representation of previous alcohol
            consumption relative to current

    raises:

    returns:
        (float): interpretation value

    """

    # Apply any relevant logic.
    if (
        (
            (not pandas.isna(alcohol_frequency)) and
            (alcohol_frequency in [1, 2, 3, 4, 5])
        ) or (
            (not pandas.isna(alcohol_status)) and
            (alcohol_status == 2)
        ) or (
            (not pandas.isna(alcohol_trend)) and
            (alcohol_trend in [0, 1, 2])
        )
    ):
        # Current alcohol consumption.
        value = 1
    elif (
        (
            (
                (not pandas.isna(alcohol_frequency)) and
                (alcohol_frequency == 0)
            ) or (
                (not pandas.isna(alcohol_status)) and
                (alcohol_status in [0, 1])
            )
        ) and (
            (pandas.isna(alcohol_trend))
        )
    ):
        # Previous or never alcohol consumption.
        value = 0
    else:
        # Missing or uninterpretable values.
        value = float("nan")
    # Return information.
    return value


def determine_alcohol_consumption_ever(
    alcohol_frequency=None,
    alcohol_status=None,
    alcohol_former=None,
    alcohol_trend=None,
):
    """
    Determine whether person ever consumed alcohol, either previously or
    currently.

    Accommodate inexact float values.

    Current alcohol consumption implies previous alcohol consumption (unless
    person is literally in the process of taking their first sip).

    arguments:
        alcohol_frequency (float): frequency of alcohol consumption
        alcohol_status (float): status of alcohol consumption
        alcohol_former (float): former alcohol consumption
        alcohol_trend (float): ordinal representation of previous alcohol
            consumption relative to current

    raises:

    returns:
        (float): binary logical representation of whether person ever consumed
            alcohol (current or previous)

    """

    # Apply any relevant logic.
    if (
        (
            (not pandas.isna(alcohol_frequency)) and
            (alcohol_frequency in [1, 2, 3, 4, 5])
        ) or
        (
            (not pandas.isna(alcohol_status)) and
            (alcohol_status in [1, 2])
        ) or (
            (not pandas.isna(alcohol_former)) and
            (alcohol_former == 1)
        ) or (
            (not pandas.isna(alcohol_trend)) and
            (alcohol_trend in [0, 1, 2])
        )
    ):
        # Previous (or current) alcohol consumption.
        value = 1
    elif (
        (
            (not pandas.isna(alcohol_frequency)) and
            (alcohol_frequency == 0)
        ) and (
            (not pandas.isna(alcohol_status)) and
            (alcohol_status == 0)
        ) and (
            (pandas.isna(alcohol_former)) or
            (alcohol_former == 0)
        ) and (
            (pandas.isna(alcohol_trend))
        )
    ):
        # Never alcohol consumption.
        value = 0
    else:
        # Missing or uninterpretable values.
        value = float("nan")
    # Return information.
    return value


def determine_alcohol_consumption_never(
    alcohol_ever=None,
):
    """
    Determine whether person never consumed alcohol, neither previously nor
    currently.

    Accommodate inexact float values.

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)

    raises:

    returns:
        (float): binary representation of whether person never consumed alcohol

    """

    # Apply any relevant logic.
    if (
        (not pandas.isna(alcohol_ever)) and
        (alcohol_ever == 1)
    ):
        # Previous (or current) alcohol consumption.
        value = 0
    elif (
        (not pandas.isna(alcohol_ever)) and
        (alcohol_ever == 0)
    ):
        # Never alcohol consumption.
        value = 1
    else:
        # Missing or uninterpretable values.
        value = float("nan")
    # Return information.
    return value


def organize_alcohol_consumption_frequency_variables(
    table=None,
    report=None,
):
    """
    Organizes information about previous and current alcohol consumption.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "1558-0.0", "3731-0.0", "1628-0.0", "20117-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Determine person's frequency of alcohol consumption.
    table["alcohol_frequency"] = table.apply(
        lambda row:
            determine_alcohol_consumption_frequency(
                field_1558=row["1558-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine person's alcohol consumption status.
    table["alcohol_status"] = table.apply(
        lambda row:
            determine_alcohol_consumption_status(
                field_20117=row["20117-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine person's former alcohol consumption status.
    table["alcohol_former"] = table.apply(
        lambda row:
            determine_former_alcohol_consumption(
                field_3731=row["3731-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine person's former alcohol consumption status.
    table["alcohol_trend"] = table.apply(
        lambda row:
            determine_alcohol_consumption_trend(
                field_1628=row["1628-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    # Determine whether person consumes alcohol currently.
    table["alcohol_current"] = table.apply(
        lambda row:
            determine_alcohol_consumption_current(
                alcohol_frequency=row["alcohol_frequency"],
                alcohol_status=row["alcohol_status"],
                alcohol_trend=row["alcohol_trend"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person consumes alcohol not currently but previously.
    #table["alcohol_previous"]
    # Determine whether person ever consumed any alcohol.
    table["alcohol_ever"] = table.apply(
        lambda row:
            determine_alcohol_consumption_ever(
                alcohol_frequency=row["alcohol_frequency"],
                alcohol_status=row["alcohol_status"],
                alcohol_former=row["alcohol_former"],
                alcohol_trend=row["alcohol_trend"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person never consumes any alcohol, either currently or
    # previously.
    table["alcohol_never"] = table.apply(
        lambda row:
            determine_alcohol_consumption_never(
                alcohol_ever=row["alcohol_ever"],
            ),
        axis="columns", # apply function to each row
    )

    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=["1558-0.0", "20117-0.0", "3731-0.0", "1628-0.0",],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "1558-0.0", "3731-0.0", "1628-0.0", "20117-0.0",
            "alcohol_frequency",
            "alcohol_status",
            "alcohol_former",
            "alcohol_trend",
            "alcohol_current",
            "alcohol_ever",
            "alcohol_never",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol consumption variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


def calculate_sum_drinks(
    beer_cider=None,
    wine_red=None,
    wine_white=None,
    port=None,
    liquor=None,
    other=None,
):
    """
    Calculates the sum of drinks.

    arguments:
        beer_cider (int): count of drinks by type
        wine_red (int): count of drinks by type
        wine_white (int): count of drinks by type
        port (int): count of drinks by type
        liquor (int): count of drinks by type
        other (int): count of drinks by type

    raises:

    returns:
        (int): sum of counts of drinks

    """

    # Consider all types of alcoholic beverage.
    types = [beer_cider, wine_red, wine_white, port, liquor, other]
    # Determine whether any relevant variables have missing values.
    valid = False
    for type in types:
        if (not pandas.isna(type) and (type > -0.5)):
            valid = True
            pass
        pass
    if valid:
        # Variables have valid values.
        # Calculate sum of drinks.
        drinks = 0
        for type in types:
            if ((not pandas.isna(type)) and (type >= 0)):
                drinks = drinks + type
                pass
            pass
    else:
        drinks = float("nan")
    return drinks


def calculate_total_alcohol_consumption_monthly(
    drinks_weekly=None,
    drinks_monthly=None,
    weeks_per_month=None,
):
    """
    Calculate monthly alcohol consumption in drinks.

    arguments:
        drinks_weekly (float): sum of weekly drinks from weekly variables
        drinks_monthly (float): sum of monthly drinks from monthly variables
        weeks_per_month (float): factor to use for weeks per month

    raises:

    returns:
        (float): person's monthly alcohol consumption in drinks

    """

    # Use as much valid information as is available.
    if (not math.isnan(drinks_weekly) and not math.isnan(drinks_monthly)):
        alcohol_drinks_monthly = (
            drinks_monthly + (weeks_per_month * drinks_weekly)
        )
    elif (not math.isnan(drinks_weekly)):
        alcohol_drinks_monthly = (weeks_per_month * drinks_weekly)
    elif (not math.isnan(drinks_monthly)):
        alcohol_drinks_monthly = drinks_monthly
    else:
        # There is no valid information about alcohol consumption quantity.
        alcohol_drinks_monthly = float("nan")
        pass
    return alcohol_drinks_monthly


def determine_total_alcohol_consumption_monthly(
    alcohol_current=None,
    drinks_weekly=None,
    drinks_monthly=None,
    weeks_per_month=None,
):
    """
    Calculate monthly alcohol consumption in drinks.

    Accommodate inexact float values.

    arguments:
        alcohol_current (float): binary representation of whether person
            consumed alcohol currently at time of recruitment, assessment
        drinks_weekly (float): sum of weekly drinks from weekly variables
        drinks_monthly (float): sum of monthly drinks from monthly variables
        weeks_per_month (float): factor to use for weeks per month

    raises:

    returns:
        (float): person's monthly alcohol consumption in drinks

    """

    # Calculate alcohol consumption quantity.
    alcohol_monthly = calculate_total_alcohol_consumption_monthly(
        drinks_weekly=drinks_weekly,
        drinks_monthly=drinks_monthly,
        weeks_per_month=weeks_per_month,
    )
    # Consider current alcohol consumption.
    if (
        (not math.isnan(alcohol_current)) and
        (alcohol_current == 0)
    ):
        # Person never consumes any alcohol currently.
        # Still prioritize quantity variables in case they are more accurate.
        if (not math.isnan(alcohol_monthly)):
            alcohol_drinks_monthly = alcohol_monthly
        else:
            alcohol_drinks_monthly = 0.0
    else:
        if (not math.isnan(alcohol_monthly)):
            alcohol_drinks_monthly = alcohol_monthly
        else:
            alcohol_drinks_monthly = float("nan")
    # Return information.
    return alcohol_drinks_monthly


# review: TCW; 07 June 2022
def determine_alcohol_drinks_monthly_ordinal(
    alcohol_ever=None,
    alcohol_current=None,
    alcohol_drinks_monthly=None,
    threshold_low=None,
    threshold_high=None,
):
    """
    Determine ordinal representation of alcohol consumption quantity for
    convenient use in cohort stratification.

    This definition uses an ordinal code.
    0: person had never consumed alcohol in life
    1: person consumed alcohol previously in life but not current at recruitment
    2: person consumed drinks of alcohol per month fewer than lower threshold
    3: person consumed drinks of alcohol per month between thresholds
    4: person consumed drinks of alcohol per month greater than higher threshold

    arguments:
        alcohol_ever (float): binary representation of whether person had ever
            consumed alcohol in lifetime
        alcohol_current (float): binary representation of whether person
            consumed alcohol currently at time of recruitment, assessment
        alcohol_drinks_monthly (float): person's monthly alcohol consumption in
            drinks
        threshold_low (int): threshold count of drinks per month
        threshold_high (int): threshold count of drinks per month

    raises:

    returns:
        (float): interpretation value

    """

    # Consider ever alcohol consumption.
    if (
        (not pandas.isna(alcohol_ever)) and
        (alcohol_ever == 0)
    ):
        # Person had never consumed any alcohol.
        value = 0
    # Consider previous but not current alcohol consumption.
    elif (
        (
            (not pandas.isna(alcohol_ever)) and
            (alcohol_ever == 1)
        ) and
        (
            (not pandas.isna(alcohol_current)) and
            (alcohol_current == 0)
        )
    ):
        # Person had consumed alcohol previously but not currently.
        value = 1
    # Consider current alcohol consumption quantity.
    elif (
        (
            (not pandas.isna(alcohol_current)) and
            (alcohol_current == 1)
        ) and
        (
            (not pandas.isna(alcohol_drinks_monthly)) and
            (alcohol_drinks_monthly < threshold_low)
        )
    ):
        # Person consumed alcohol currently at quantity less than lower
        # threshold.
        value = 2
    elif (
        (
            (not pandas.isna(alcohol_current)) and
            (alcohol_current == 1)
        ) and
        (
            (not pandas.isna(alcohol_drinks_monthly)) and
            (
                (alcohol_drinks_monthly >= threshold_low) and
                (alcohol_drinks_monthly < threshold_high)
            )
        )
    ):
        # Person consumed alcohol currently at quantity between the lower and
        # higher thresholds.
        value = 3
    elif (
        (
            (not pandas.isna(alcohol_current)) and
            (alcohol_current == 1)
        ) and
        (
            (not pandas.isna(alcohol_drinks_monthly)) and
            (alcohol_drinks_monthly >= threshold_high)
        )
    ):
        # Person consumed alcohol currently at quantity greater than higher
        # threshold.
        value = 4
    else:
        # Ambiguous, uninterpretable, or missing information.
        value = float("nan")
    # Return information.
    return value


def organize_alcohol_consumption_quantity_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcohol consumption in standard drinks monthly.

    UK Biobank data coding 100291.
    "do not know": -1
    "prefer not to answer": -3

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "1588-0.0", "1568-0.0", "1578-0.0", "1608-0.0", "1598-0.0", "5364-0.0",
        "4429-0.0", "4407-0.0", "4418-0.0", "4451-0.0", "4440-0.0", "4462-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Calculate sum of drinks weekly.
    table["alcohol_drinks_weekly"] = table.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["1588-0.0"],
                wine_red=row["1568-0.0"],
                wine_white=row["1578-0.0"],
                port=row["1608-0.0"],
                liquor=row["1598-0.0"],
                other=row["5364-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Calculate sum of drinks monthly.
    table["alcohol_drinks_monthly"] = table.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["4429-0.0"],
                wine_red=row["4407-0.0"],
                wine_white=row["4418-0.0"],
                port=row["4451-0.0"],
                liquor=row["4440-0.0"],
                other=row["4462-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine sum of total drinks monthly.
    table["alcohol_drinks_monthly_combination"] = table.apply(
        lambda row:
            determine_total_alcohol_consumption_monthly(
                alcohol_current=row["alcohol_current"],
                drinks_weekly=row["alcohol_drinks_weekly"],
                drinks_monthly=row["alcohol_drinks_monthly"],
                weeks_per_month=4.345, # 52.143 weeks per year (12 months)
            ),
        axis="columns", # apply function to each row
    )
    # Determine ordinal variable on alcohol consumption quantity for
    # convenience in cohort stratification.
    table["alcohol_drinks_monthly_ordinal"] = table.apply(
        lambda row:
            determine_alcohol_drinks_monthly_ordinal(
                alcohol_ever=row["alcohol_ever"],
                alcohol_current=row["alcohol_current"],
                alcohol_drinks_monthly=(
                    row["alcohol_drinks_monthly_combination"]
                ),
                threshold_low=31, #
                threshold_high=61, #
            ),
        axis="columns", # apply function to each row
    )

    ##########
    # Transform variables' values to normalize distributions.
    columns_normalization = [
        "alcohol_drinks_weekly",
        "alcohol_drinks_monthly",
        "alcohol_drinks_monthly_combination",
    ]
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=columns_normalization,
        table=table,
    )

    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "1588-0.0", "1568-0.0", "1578-0.0", "1608-0.0", "1598-0.0",
            "5364-0.0",
            "4429-0.0", "4407-0.0", "4418-0.0", "4451-0.0", "4440-0.0",
            "4462-0.0",
        ],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "1588-0.0", "1568-0.0", "1578-0.0", "1608-0.0", "1598-0.0",
            "5364-0.0",
            "4429-0.0", "4407-0.0", "4418-0.0", "4451-0.0", "4440-0.0",
            "4462-0.0",
            "alcohol_drinks_weekly",
            "alcohol_drinks_monthly",
            "alcohol_drinks_monthly_combination",
            "alcohol_current",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol consumption quantity variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


def organize_alcohol_consumption_variables(
    table=None,
    report=None,
):
    """
    Organizes information about previous and current alcohol consumption.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Organize information about alcohol consumption and frequency.
    pail_consumption = organize_alcohol_consumption_frequency_variables(
        table=table,
        report=report,
    )
    # Organize information about current alcohol consumption quantity.
    if True:
        pail_quantity = organize_alcohol_consumption_quantity_variables(
            table=pail_consumption["table"],
            report=report,
        )
    # Return information.
    return pail_quantity


##########
# Alcohol AUDIT questionnaire


def interpret_alcohol_audit_one(
    value=None,
):
    """
    Intepret UK Biobank's coding for AUDIT questionnaire question 1.

    AUDIT 1
    Note: "Question asked: 'How often do you have a drink containing alcohol?'"

    Data-Field "20414": "Frequency of drinking alcohol"
    UK Biobank data-coding "521" for data-field "20414".
    0: "never"
    1: "monthly or less"
    2: "two to four times a month"
    3: "two to three times a week"
    4: "four or more times a week"
    -818: "prefer not to answer"

    Accommodate inexact float values.

    arguments:
        value (float): raw value from UK Biobank's coding

    raises:

    returns:
        (float): interpretation value

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(value)) and
        (-818.5 <= value and value < 4.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= value and value < 0.5):
            # 0: "never"
            interpretation = 0
        elif (0.5 <= value and value < 1.5):
            # 1: "monthly or less"
            interpretation = 1
        elif (1.5 <= value and value < 2.5):
            # 2: "two to four times a month"
            interpretation = 2
        elif (2.5 <= value and value < 3.5):
            # 3: "two to three times a week"
            interpretation = 3
        elif (3.5 <= value and value < 4.5):
            # 4: "four or more times a week"
            interpretation = 4
        elif (-818.5 <= value and value < -817.5):
            # -818: "prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # "prefer not to answer" or null
        interpretation = float("nan")
    # Return information.
    return interpretation


def interpret_alcohol_audit_two_to_eight(
    value=None,
):
    """
    Intepret UK Biobank's coding for AUDIT questionnaire questions 2 to 8.

    AUDIT 2
    Note (20403): "Question asked: 'In the next two questions, a "drink" is
    defined as one unit of alcohol. How many drinks containing alcohol do you
    have on a typical day when you are drinking?'"

    UK Biobank data-coding "522" for data-field "20403".
    1: "one or two"
    2: "three or four"
    3: "five or six"
    4: "seven, eight, or nine"
    5: "ten or more"
    -818: "prefer not to answer"

    AUDIT 3
    Note (20416): "Question asked: 'In the next two questions, a "drink" is
    defined as one unit of alcohol. How often do you have six or more drinks on
    one occasion?'"

    AUDIT 4
    Note (20413): "Question asked: 'How often during the last year have you
    found that you were not able to stop drinking once you had started?'"

    AUDIT 5
    Note (20407): "Question asked: 'How often during the last year have you
    failed to do what was normally expected from you because of drinking?'"

    AUDIT 6
    Note (20412): "Question asked: 'How often during the last year have you
    needed a first drink in the morning to get yourself going after a heavy
    drinking session?'"

    AUDIT 7
    Note (20409): "Question asked: 'How often during the last year have you had
    a feeling of guilt or remorse after drinking?'"

    AUDIT 8
    Note (20408): "Question asked: 'How often during the last year have you been
    unable to remember what happened the night before because you had been
    drinking?'"

    UK Biobank data-coding "523" for data-fields "20416", "20413", "20407",
    "20412", "20409", and "20408".
    1: "never"
    2: "less than monthly"
    3: "monthly"
    4: "weekly"
    5: "daily or almost daily"
    -818: "prefer not to answer"

    Notice that UK Biobank encodes these variables from 1 to 5; however, AUDIT
    encodes these from 0 to 4. Adjust the values accordingly.

    Accommodate inexact float values.

    arguments:
        value (float): raw value from UK Biobank's coding

    raises:

    returns:
        (float): interpretation value

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(value)) and
        (-818.5 <= value and value < 5.5)
    ):
        # The variable has a valid value.
        if (0.5 <= value and value < 1.5):
            # 1: "one or two" (code "522")
            # 1: "never" (code "523")
            interpretation = 0
        elif (1.5 <= value and value < 2.5):
            # 2: "three or four" (code "522")
            # 2: "less than monthly" (code "523")
            interpretation = 1
        elif (2.5 <= value and value < 3.5):
            # 3: "five or six" (code "522")
            # 3: "monthly" (code "523")
            interpretation = 2
        elif (3.5 <= value and value < 4.5):
            # 4: "seven, eight, or nine" (code "522")
            # 4: "weekly" (code "523")
            interpretation = 3
        elif (4.5 <= value and value < 5.5):
            # 5: "ten or more" (code "522")
            # 5: "daily or almost daily" (code "523")
            interpretation = 4
        elif (-818.5 <= value and value < -817.5):
            # -818: "prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # "prefer not to answer" or null
        interpretation = float("nan")
    # Return information.
    return interpretation


def interpret_alcohol_audit_nine_ten(
    value=None,
):
    """
    Intepret UK Biobank's coding for AUDIT questionnaire questions 9 and 10.

    AUDIT 9
    Note (20411): "Question asked: 'Have you or someone else been injured as a
    result of your drinking?'"

    AUDIT 10
    Note (20405): "Question asked: 'Has a relative or friend or a doctor or
    another health worker been concerned about your drinking or suggested you
    cut down?'"

    UK Biobank data-coding "524" for data-fields "20411", and "20405".
    0: "no"
    1: "yes, but not in the last year"
    2: "yes, during the last year"
    -818: "prefer not to answer"

    Notice that UK Biobank encodes these variables as 0, 1, or 2; however,
    AUDIT encodes these as 0, 2, or 4. Adjust the values accordingly.

    Accommodate inexact float values.

    arguments:
        value (float): raw value from UK Biobank's coding

    raises:

    returns:
        (float): interpretation value

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(value)) and
        (-818.5 <= value and value < 2.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= value and value < 0.5):
            # 0: "no"
            interpretation = 0
        elif (0.5 <= value and value < 1.5):
            # 1: "yes, but not in the last year"
            interpretation = 2
        elif (1.5 <= value and value < 2.5):
            # 2: "yes, during the last year"
            interpretation = 4
        elif (-818.5 <= value and value < -817.5):
            # -818: "prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # "prefer not to answer" or null
        interpretation = float("nan")
    # Return information.
    return interpretation


def determine_alcohol_auditc_score(
    audit_1=None,
    audit_2=None,
    audit_3=None,
):
    """
    Determine a person's score from the AUDIT-C questionnaire.

    Accommodate inexact float values.

    arguments:
        audit_1 (float): AUDIT questionnaire question 1, UK Biobank field
            20414
        audit_2 (float): AUDIT questionnaire question 2, UK Biobank field
            20403
        audit_3 (float): AUDIT questionnaire question 3, UK
            Biobank field 20416

    raises:

    returns:
        (float): combination score for AUDIT-C questionnaire

    """

    # Interpret raw variables.
    audit_1_clean = interpret_alcohol_audit_one(
        value=audit_1,
    )
    audit_2_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_2,
    )
    audit_3_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_3,
    )
    # Integrate information from multiple variables.
    if (
        (not pandas.isna(audit_1_clean)) and
        (not pandas.isna(audit_2_clean)) and
        (not pandas.isna(audit_3_clean))
    ):
        auditc_score = (audit_1_clean + audit_2_clean + audit_3_clean)
    else:
        auditc_score = float("nan")
    # Return information.
    return auditc_score


def organize_alcohol_auditc_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcohol AUDIT-C questionnaire.

    The official AUDIT-C questionnaire scores each question on 0 to 4 points.
    https://cde.drugabuse.gov/instrument/f229c68a-67ce-9a58-e040-bb89ad432be4

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "20414-0.0", "20403-0.0", "20416-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Determine person's AUDIT-C score.
    table["alcohol_auditc"] = table.apply(
        lambda row:
            determine_alcohol_auditc_score(
                audit_1=row["20414-0.0"],
                audit_2=row["20403-0.0"],
                audit_3=row["20416-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=["20414-0.0", "20403-0.0", "20416-0.0",],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "20414-0.0", "20403-0.0", "20416-0.0",
            "alcohol_auditc",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol AUDIT-C variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


def determine_alcohol_auditp_score(
    audit_4=None,
    audit_5=None,
    audit_6=None,
    audit_7=None,
    audit_8=None,
    audit_9=None,
    audit_10=None,
):
    """
    Determine a person's score from the AUDIT-P questionnaire.

    Accommodate inexact float values.

    arguments:
        audit_4 (float): AUDIT questionnaire question 4, UK Biobank field
            20413
        audit_5 (float): AUDIT questionnaire question 5, UK Biobank field
            20407
        audit_6 (float): AUDIT questionnaire question 6, UK Biobank field
            20412
        audit_7 (float): AUDIT questionnaire question 7, UK Biobank field
            20409
        audit_8 (float): AUDIT questionnaire question 8, UK Biobank field
            20408
        audit_9 (float): AUDIT questionnaire question 9, UK Biobank field
            20411
        audit_10 (float): AUDIT questionnaire question 10, UK Biobank field
            20405

    raises:

    returns:
        (float): combination score for AUDIT questionnaire

    """

    # Interpret raw variables.
    audit_4_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_4,
    )
    audit_5_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_5,
    )
    audit_6_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_6,
    )
    audit_7_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_7,
    )
    audit_8_clean = interpret_alcohol_audit_two_to_eight(
        value=audit_8,
    )
    audit_9_clean = interpret_alcohol_audit_nine_ten(
        value=audit_9,
    )
    audit_10_clean = interpret_alcohol_audit_nine_ten(
        value=audit_10,
    )

    # Integrate information from multiple variables.
    if (
        (not pandas.isna(audit_4_clean)) and
        (not pandas.isna(audit_5_clean)) and
        (not pandas.isna(audit_6_clean)) and
        (not pandas.isna(audit_7_clean)) and
        (not pandas.isna(audit_8_clean)) and
        (not pandas.isna(audit_9_clean)) and
        (not pandas.isna(audit_10_clean))
    ):
        auditp_score = (
            audit_4_clean + audit_5_clean + audit_6_clean + audit_7_clean +
            audit_8_clean + audit_9_clean + audit_10_clean
        )
    else:
        auditp_score = float("nan")
    # Return information.
    return auditp_score


def organize_alcohol_auditp_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcohol AUDIT-P questionnaire.

    The official AUDIT questionnaire scores questions 1-8 on 0 to 4 points and
    questions 9-10 on 0, 2, or 4 points.
    https://auditscreen.org/about/scoring-audit/

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "20413-0.0", "20407-0.0", "20412-0.0", "20409-0.0", "20408-0.0",
        "20411-0.0", "20405-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Determine person's AUDIT-C score.
    table["alcohol_auditp"] = table.apply(
        lambda row:
            determine_alcohol_auditp_score(
                audit_4=row["20413-0.0"],
                audit_5=row["20407-0.0"],
                audit_6=row["20412-0.0"],
                audit_7=row["20409-0.0"],
                audit_8=row["20408-0.0"],
                audit_9=row["20411-0.0"],
                audit_10=row["20405-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=columns_type,
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "20413-0.0", "20407-0.0", "20412-0.0", "20409-0.0", "20408-0.0",
            "20411-0.0", "20405-0.0",
            "alcohol_auditp",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol AUDIT-P variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


def determine_alcohol_audit_score(
    alcohol_auditc=None,
    alcohol_auditp=None,
):
    """
    Determine a person's score from the AUDIT questionnaire.

    Accommodate inexact float values.

    arguments:
        alcohol_auditc ((float): combination score for AUDIT-C questionnaire
        alcohol_auditp ((float): combination score for AUDIT-P questionnaire

    raises:

    returns:
        (float): combination score for AUDIT questionnaire

    """

    # Integrate information from multiple variables.
    if (
        (not pandas.isna(alcohol_auditc)) and
        (not pandas.isna(alcohol_auditp))
    ):
        audit_score = (alcohol_auditc + alcohol_auditp)
    else:
        audit_score = float("nan")
    # Return information.
    return audit_score


def organize_alcohol_audit_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcohol AUDIT-P questionnaire.

    The official AUDIT questionnaire scores questions 1-8 on 0 to 4 points and
    questions 9-10 on 0, 2, or 4 points.
    https://auditscreen.org/about/scoring-audit/

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Determine person's AUDIT score.
    table["alcohol_audit"] = table.apply(
        lambda row:
            determine_alcohol_audit_score(
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
            ),
        axis="columns", # apply function to each row
    )
    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "alcohol_auditc", "alcohol_auditp", "alcohol_audit",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol AUDIT variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


def organize_alcohol_audit_questionnaire_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcohol AUDIT questionnaire.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Organize information about alcohol AUDIT-C.
    pail_auditc = organize_alcohol_auditc_variables(
        table=table,
        report=report,
    )
    # Organize information about alcohol AUDIT.
    pail_auditp = organize_alcohol_auditp_variables(
        table=pail_auditc["table"],
        report=report,
    )
    # Organize information about alcohol AUDIT.
    pail_audit = organize_alcohol_audit_variables(
        table=pail_auditp["table"],
        report=report,
    )
    # Return information.
    return pail_audit


##########
# Alcoholism diagnoses


def specify_icd_alcoholism_diagnosis_groups_codes():
    """
    Specifies the codes for diagnostic groups relevant to alcoholism.

    arguments:

    raises:

    returns:
        (dict<dict<str>>): collections of codes for diagnostic groups

    """

    # Determine whether the variable has a valid (non-missing) value.
    codes = dict()

    codes["group_a"] = dict()
    codes["group_a"]["icd_9"] = [
        "3039", "291", "2910", "2918",
    ]
    codes["group_a"]["icd_10"] = [
        "F102", "F103", "F104", "F105", "F106", "F107", "F108",
    ]

    codes["group_b"] = dict()
    codes["group_b"]["icd_9"] = []
    codes["group_b"]["icd_10"] = [
        "F10", "F100", "F109",
    ]

    codes["group_c"] = dict()
    codes["group_c"]["icd_9"] = [
        "3050",
    ]
    codes["group_c"]["icd_10"] = [
        "F101",
    ]

    codes["group_d"] = dict()
    codes["group_d"]["icd_9"] = [
        "3575", "4255", "5353", "5710", "5711", "5712", "5713",
    ]
    codes["group_d"]["icd_10"] = [
        "G621", "I426", "K292", "K70", "K700", "K701", "K702", "K703", "K704",
        "K709",
    ]

    # Return information.
    return codes


def specify_self_alcoholism_diagnosis_codes():
    """
    Specifies the codes for diagnostic groups relevant to alcoholism.

    arguments:

    raises:

    returns:
        (list<str>): codes for self diagnosis

    """

    # Specify codes.
    codes = [
        "1408", "1604",
    ]
    # Return information.
    return codes


# TODO: TCW 14 September 2021
# TODO: I could possibly rewrite this more efficiently using "any()" and "filter()"
# TODO: the matching driver function in 'promiscuity.utility' could return TRUE as soon as it finds the first match


def parse_interpret_match_diagnosis_codes(
    row=None,
    fields=None,
    codes_match=None,
):
    """
    Parses and interprets ICD diagnosis codes.

    Any match in either ICD9 or ICD10 codes suffices for a diagnostic match.

    arguments:
        row (object): Pandas series corresponding to a row of a Pandas data
            frame
        fields (list<str>): names of columns for UK Biobank fields for
            diagnostic codes
        codes_match (list<str>): codes corresponding to a diagnostic
            group

    raises:

    returns:
        (bool): whether values in the current row match the diagnostic group

    """

    matches = list()
    for field in fields:
        collection = row[field]
        values_actual = utility.parse_text_list_values(
            collection=collection,
            delimiter=";",
        )
        match = utility.determine_any_actual_values_match_comparisons(
            values_actual=values_actual,
            values_comparison=codes_match,
        )
        matches.append(match)
    return any(matches)


def determine_diagnosis_icd_alcoholism_group(
    row=None,
    fields_icd_9=None,
    fields_icd_10=None,
    codes_icd_group=None,
):
    """
    Organizes information about diagnoses relevant to alcoholism.

    Do not propagate missing values in ICD codes as the absence of a code does
    not necessarily indicate missing information. It could indicate a valid
    null.

    arguments:
        row (object): Pandas series corresponding to a row of a Pandas data
            frame
        fields_icd_9 (list<str>): names of columns for UK Biobank fields for
            ICD9 diagnostic codes
        fields_icd_10 (list<str>): names of columns for UK Biobank fields for
            ICD10 diagnostic codes
        codes_icd_group (dict<list<str>>): ICD9 and ICD10 codes corresponding
            to a diagnostic group

    raises:

    returns:
        (bool): whether values in the current row match the diagnostic group

    """

    # Determine whether any codes in any ICD9 fields match diagnostic group.
    icd_9_group = parse_interpret_match_diagnosis_codes(
        row=row,
        fields=fields_icd_9,
        codes_match=codes_icd_group["icd_9"],
    )
    # Determine whether any codes in any ICD10 fields match diagnostic group.
    icd_10_group = parse_interpret_match_diagnosis_codes(
        row=row,
        fields=fields_icd_10,
        codes_match=codes_icd_group["icd_10"],
    )
    # Interpret matches from either ICD9 or ICD10 codes.
    if (icd_9_group or icd_10_group):
        match = 1
    else:
        match = 0
    # Return information.
    return match


def determine_diagnosis_alcoholism_self(
    row=None,
    fields=None,
    codes_match=None,
):
    """
    Organizes information about diagnoses relevant to alcoholism.

    Do not propagate missing values in ICD codes as the absence of a code does
    not necessarily indicate missing information. It could indicate a valid
    null.

    arguments:
        row (object): Pandas series corresponding to a row of a Pandas data
            frame
        fields (list<str>): names of columns for UK Biobank fields for self
            report of alcohol dependence
        codes_match (dict<list<str>>): diagnostic codes for self report of
            alcohol dependence

    raises:

    returns:
        (bool): whether values in the current row match the diagnostic group

    """

    # Determine whether any codes in any ICD9 fields match diagnostic group.
    self_report = parse_interpret_match_diagnosis_codes(
        row=row,
        fields=fields,
        codes_match=codes_match,
    )
    # Interpret matches from either ICD9 or ICD10 codes.
    if (self_report):
        match = 1
    else:
        match = 0
    # Return information.
    return match


def organize_alcoholism_diagnosis_variables(
    table=None,
    report=None,
):
    """
    Organizes information about diagnoses relevant to alcoholism.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about alcoholism diagnoses

    """

    # Copy data.
    table = table.copy(deep=True)
    # Specify ICD9 and ICD10 diagnostic codes relevant to alcoholism.
    codes_icd = specify_icd_alcoholism_diagnosis_groups_codes()
    # Specify relevant codes from self diagnosis.
    codes_self = specify_self_alcoholism_diagnosis_codes()

    # Determine whether person has diagnoses in group A.
    table["alcohol_diagnosis_a"] = table.apply(
        lambda row:
            determine_diagnosis_icd_alcoholism_group(
                row=row.copy(deep=True),
                fields_icd_9=["41271_array", "41203_array", "41205_array",],
                fields_icd_10=["41270_array", "41202_array", "41204_array",],
                codes_icd_group=codes_icd["group_a"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person has diagnoses in group B.
    table["alcohol_diagnosis_b"] = table.apply(
        lambda row:
            determine_diagnosis_icd_alcoholism_group(
                row=row.copy(deep=True),
                fields_icd_9=["41271_array", "41203_array", "41205_array",],
                fields_icd_10=["41270_array", "41202_array", "41204_array",],
                codes_icd_group=codes_icd["group_b"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person has diagnoses in group C.
    table["alcohol_diagnosis_c"] = table.apply(
        lambda row:
            determine_diagnosis_icd_alcoholism_group(
                row=row.copy(deep=True),
                fields_icd_9=["41271_array", "41203_array", "41205_array",],
                fields_icd_10=["41270_array", "41202_array", "41204_array",],
                codes_icd_group=codes_icd["group_c"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person has diagnoses in group D.
    table["alcohol_diagnosis_d"] = table.apply(
        lambda row:
            determine_diagnosis_icd_alcoholism_group(
                row=row.copy(deep=True),
                fields_icd_9=["41271_array", "41203_array", "41205_array",],
                fields_icd_10=["41270_array", "41202_array", "41204_array",],
                codes_icd_group=codes_icd["group_d"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine whether person has self diagnoses.
    table["alcohol_diagnosis_self"] = table.apply(
        lambda row:
            determine_diagnosis_alcoholism_self(
                row=row.copy(deep=True),
                fields=["20002_array"],
                codes_match=codes_self,
            ),
        axis="columns", # apply function to each row
    )

    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "41271_array", "41203_array", "41205_array",
            "41270_array", "41202_array", "41204_array",
            "20002_array",
        ],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    #table_report = table_report.loc[
    #    table_report["alcohol_diagnosis_a"] == True, :
    #]
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "41271_array", "41203_array", "41205_array", "41270_array",
            "41202_array", "41204_array", "20002_array",
            "alcohol_diagnosis_a",
            "alcohol_diagnosis_b",
            "alcohol_diagnosis_c",
            "alcohol_diagnosis_d",
            "alcohol_diagnosis_self",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol diagnosis variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Alcoholism cases and controls


def determine_alcoholism_control(
    alcohol_ever=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditp_control=None,
    threshold_audit_control=None,
):
    """
    Determines whether person qualifies as a control for alcoholism.

    Alcoholism control criteria:
    - person has consumed alcohol ("ever" consumer)
    - person did not self report alcohol dependence
    - person did not have any ICD9 or ICD10 codes indicative of alcohol
       dependence
    - person had moderate scores in AUDIT-C, AUDIT-P, or AUDIT questionnaires

    Missing values exist in AUDIT-C, AUDIT-P, and AUDIT variables.

    Missing values do not exist in diagnosis variables.

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)
        alcohol_diagnosis_a (bool): whether person has any diagnoses in
            alcoholism group A
        alcohol_diagnosis_b (bool): whether person has any diagnoses in
            alcoholism group B
        alcohol_diagnosis_c (bool): whether person has any diagnoses in
            alcoholism group C
        alcohol_diagnosis_d (bool): whether person has any diagnoses in
            alcoholism group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): control threshold (less than) for
            AUDIT-C score
        threshold_auditp_control (float): control threshold (less than) for
            AUDIT-P score
        threshold_audit_control (float): control threshold (less than) for
            AUDIT score

    raises:

    returns:
        (float): binary logical representation of alcoholism control status

    """

    # Determine whether person ever consumed alcohol (currently or previously).
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as controls for alcoholism.
    if (not pandas.isna(alcohol_ever)):
        if (alcohol_ever == 1):
            # Person consumed alcohol currently or previously.
            match_consumption_ever = True
        else:
            # Person never consumed alcohol.
            match_consumption_ever = False
    else:
        # Missing information.
        match_consumption_ever = False

    # Determine whether person did not qualify for any diagnostic definitions of
    # alcoholism.
    if (
        (not pandas.isna(alcohol_diagnosis_a)) and
        (not pandas.isna(alcohol_diagnosis_b)) and
        (not pandas.isna(alcohol_diagnosis_c)) and
        (not pandas.isna(alcohol_diagnosis_d)) and
        (not pandas.isna(alcohol_diagnosis_self))
    ):
        if (
            (alcohol_diagnosis_a == 0) and
            (alcohol_diagnosis_b == 0) and
            (alcohol_diagnosis_c == 0) and
            (alcohol_diagnosis_d == 0) and
            (alcohol_diagnosis_self == 0)
        ):
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_diagnosis = True
        else:
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_diagnosis = False
    else:
        # Missing information.
        # Missing information in diagnostic variables does disqualify control.
        # But diagnostic variables do not have missing values.
        match_diagnosis = False

    # Determine whether person had only moderate scores in AUDIT questionnaires.
    if (
        (not pandas.isna(alcohol_auditc)) and
        (not pandas.isna(alcohol_auditp)) and
        (not pandas.isna(alcohol_audit))
    ):
        if (
            (alcohol_auditc < threshold_auditc_control) and
            (alcohol_auditp < threshold_auditp_control) and
            (alcohol_audit < threshold_audit_control)
        ):
            # Person only had moderate scores in AUDIT questionnaires.
            match_questionnaire = True
        else:
            # Person had a more extreme score in at least one AUDIT
            # questionnaire.
            match_questionnaire = False
    else:
        # Missing information.
        # Missing information in AUDIT questionnaire does not disqualify
        # control.
        match_questionnaire = True

    # Integrate information from all criteria.
    if (match_consumption_ever and match_diagnosis and match_questionnaire):
        # Person qualifies as a control for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_case_one(
    alcohol_ever=None,
    alcohol_diagnosis_a=None,
):
    """
    Organizes information about alcoholism cases and controls.

    Alcoholism case definition one:
    - current or previous alcohol consumption ("ever")
    - ICD diagnostic group A

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)
        alcohol_diagnosis_a (bool): whether person has any diagnoses in
            alcoholism group A

    raises:

    returns:
        (float): binary logical representation of alcoholism case status

    """

    # Determine whether person ever consumed alcohol (currently or previously).
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as cases for alcoholism.
    if (not pandas.isna(alcohol_ever)):
        if (alcohol_ever == 1):
            # Person consumed alcohol currently or previously.
            match_consumption_ever = True
        else:
            # Person never consumed alcohol.
            match_consumption_ever = False
    else:
        # Missing information.
        match_consumption_ever = False

    # Determine whether person qualifies as a case for alcoholism.
    if (
        (not pandas.isna(alcohol_diagnosis_a))
    ):
        if (
            (alcohol_diagnosis_a == 1)
        ):
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_alcoholism = True
        else:
            # Person did not qualify by any diagnostic definitions of
            # alcoholism.
            match_alcoholism = False
    else:
        # Missing information.
        match_alcoholism = False

    # Integrate information from all criteria.
    if (match_consumption_ever and match_alcoholism):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_case_two(
    alcohol_ever=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
):
    """
    Organizes information about alcoholism cases and controls.

    Alcoholism case definition two:
    - current or previous alcohol consumption ("ever")
    - ICD diagnostic group A, B, C, or D

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)
        alcohol_diagnosis_a (bool): whether person has any diagnoses in
            alcoholism group A
        alcohol_diagnosis_b (bool): whether person has any diagnoses in
            alcoholism group B
        alcohol_diagnosis_c (bool): whether person has any diagnoses in
            alcoholism group C
        alcohol_diagnosis_d (bool): whether person has any diagnoses in
            alcoholism group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism

    raises:

    returns:
        (float): binary logical representation of alcoholism case status

    """

    # Determine whether person ever consumed alcohol (currently or previously).
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as cases for alcoholism.
    if (not pandas.isna(alcohol_ever)):
        if (alcohol_ever == 1):
            # Person consumed alcohol currently or previously.
            match_consumption_ever = True
        else:
            # Person never consumed alcohol.
            match_consumption_ever = False
    else:
        # Missing information.
        match_consumption_ever = False

    # Determine whether person qualifies as a case for alcoholism.
    if (
        (not pandas.isna(alcohol_diagnosis_a)) and
        (not pandas.isna(alcohol_diagnosis_b)) and
        (not pandas.isna(alcohol_diagnosis_c)) and
        (not pandas.isna(alcohol_diagnosis_d))# and
        #(not pandas.isna(alcohol_diagnosis_self))
    ):
        if (
            (alcohol_diagnosis_a == 1) or
            (alcohol_diagnosis_b == 1) or
            (alcohol_diagnosis_c == 1) or
            (alcohol_diagnosis_d == 1)# or
            #(alcohol_diagnosis_self == 1)
        ):
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_alcoholism = True
        else:
            # Person did not qualify by any diagnostic definitions of
            # alcoholism.
            match_alcoholism = False
    else:
        # Missing information.
        match_alcoholism = False

    # Integrate information from all criteria.
    if (match_consumption_ever and match_alcoholism):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_case_three(
    alcohol_ever=None,
    alcohol_auditc=None,
    threshold_auditc_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    Alcoholism case definition three:
    - current or previous alcohol consumption ("ever")
    - AUDIT-C beyond case threshold

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        threshold_auditc_case (float): case threshold (greater than or equal)
            for AUDIT-C score

    raises:

    returns:
        (float): binary logical representation of alcoholism case status

    """

    # Determine whether person ever consumed alcohol (currently or previously).
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as cases for alcoholism.
    if (not pandas.isna(alcohol_ever)):
        if (alcohol_ever == 1):
            # Person consumed alcohol currently or previously.
            match_consumption_ever = True
        else:
            # Person never consumed alcohol.
            match_consumption_ever = False
    else:
        # Missing information.
        match_consumption_ever = False

    # Determine whether person qualifies as a case for alcoholism.
    if (
        (not pandas.isna(alcohol_auditc))
    ):
        if (
            (alcohol_auditc >= threshold_auditc_case)
        ):
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_alcoholism = True
        else:
            # Person did not qualify by any diagnostic definitions of
            # alcoholism.
            match_alcoholism = False
    else:
        # Missing information.
        match_alcoholism = False

    # Integrate information from all criteria.
    if (match_consumption_ever and match_alcoholism):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_case_four(
    alcohol_ever=None,
    alcohol_auditp=None,
    threshold_auditp_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    Alcoholism case definition four:
    - current or previous alcohol consumption ("ever")
    - AUDIT-P beyond case threshold

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        threshold_auditp_case (float): case threshold (greater than or equal)
            for AUDIT-P score

    raises:

    returns:
        (float): binary logical representation of alcoholism case status

    """

    # Determine whether person ever consumed alcohol (currently or previously).
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as cases for alcoholism.
    if (not pandas.isna(alcohol_ever)):
        if (alcohol_ever == 1):
            # Person consumed alcohol currently or previously.
            match_consumption_ever = True
        else:
            # Person never consumed alcohol.
            match_consumption_ever = False
    else:
        # Missing information.
        match_consumption_ever = False

    # Determine whether person qualifies as a case for alcoholism.
    if (
        (not pandas.isna(alcohol_auditp))
    ):
        if (
            (alcohol_auditp >= threshold_auditp_case)
        ):
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_alcoholism = True
        else:
            # Person did not qualify by any diagnostic definitions of
            # alcoholism.
            match_alcoholism = False
    else:
        # Missing information.
        match_alcoholism = False

    # Integrate information from all criteria.
    if (match_consumption_ever and match_alcoholism):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_case_five(
    alcohol_ever=None,
    alcohol_audit=None,
    threshold_audit_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    Alcoholism case definition five:
    - current or previous alcohol consumption ("ever")
    - AUDIT beyond case threshold

    arguments:
        alcohol_ever (float): binary logical representation of whether person
            ever consumed alcohol (current or previous)
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_audit_case (float): case threshold (greater than or equal)
            for AUDIT score

    raises:

    returns:
        (float): binary logical representation of alcoholism case status

    """

    # Determine whether person ever consumed alcohol (currently or previously).
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as cases for alcoholism.
    if (not pandas.isna(alcohol_ever)):
        if (alcohol_ever == 1):
            # Person consumed alcohol currently or previously.
            match_consumption_ever = True
        else:
            # Person never consumed alcohol.
            match_consumption_ever = False
    else:
        # Missing information.
        match_consumption_ever = False

    # Determine whether person qualifies as a case for alcoholism.
    if (
        (not pandas.isna(alcohol_audit))
    ):
        if (
            (alcohol_audit >= threshold_audit_case)
        ):
            # Person qualified for at least one diagnostic definition of
            # alcoholism.
            match_alcoholism = True
        else:
            # Person did not qualify by any diagnostic definitions of
            # alcoholism.
            match_alcoholism = False
    else:
        # Missing information.
        match_alcoholism = False

    # Integrate information from all criteria.
    if (match_consumption_ever and match_alcoholism):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_case_any(
    alcoholism_case_1=None,
    alcoholism_case_2=None,
    alcoholism_case_3=None,
    alcoholism_case_4=None,
    alcoholism_case_5=None,
):
    """
    Determines whether person qualifies as a case for alcoholism.

    arguments:
        alcoholism_case_1 (float): binary logical representation of alcoholism
            case status type 1
        alcoholism_case_2 (float): binary logical representation of alcoholism
            case status type 2
        alcoholism_case_3 (float): binary logical representation of alcoholism
            case status type 3
        alcoholism_case_4 (float): binary logical representation of alcoholism
            case status type 4
        alcoholism_case_5 (float): binary logical representation of alcoholism
            case status type 5

    raises:

    returns:
        (float): binary logical representation of any alcoholism case status

    """

    # Integrate information from all criteria.
    if (
        (alcoholism_case_1 == 1) or
        (alcoholism_case_2 == 1) or
        (alcoholism_case_3 == 1) or
        (alcoholism_case_4 == 1) or
        (alcoholism_case_5 == 1)
    ):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        value = 0
    # Return information.
    return value


def determine_alcoholism_control_case(
    alcoholism_control=None,
    alcoholism_case=None,
):
    """
    Determines whether person qualifies as a control or case for alcoholism.

    Ensure mutual exclusivity of alcoholism control and case statuses.

    NA: missing value indicates neither control nor case
    0: control
    1: case

    arguments:
        alcoholism_control (float): binary logical representation of alcoholism
            control status
        alcoholism_case (float): binary logical representation of alcoholism
            case status

    raises:

    returns:
        (float): binary logical representation of alcoholism control or case
            status

    """

    # Integrate information from all criteria.
    if (
        (
            (not pandas.isna(alcoholism_control)) and
            (alcoholism_control == 1)
        ) and
        (
            (pandas.isna(alcoholism_case)) or
            (alcoholism_case == 0)
        )
    ):
        # Person qualifies as a control for alcoholism.
        value = 0
    elif (
        (
            (not pandas.isna(alcoholism_case)) and
            (alcoholism_case == 1)
        ) and
        (
            (pandas.isna(alcoholism_control)) or
            (alcoholism_control == 0)
        )
    ):
        # Person qualifies as a case for alcoholism.
        value = 1
    else:
        # Missing information.
        value = float("nan")
    # Return information.
    return value


def organize_report_cohort_by_sex_alcoholism_split_hormone(
    sex_text=None,
    alcohol_consumption=None,
    alcoholism=None,
    alcoholism_split=None,
    hormone=None,
    variables_names_valid=None,
    variables_prefixes_valid=None,
    table=None,
):
    """
    Organizes information about previous and current alcohol consumption.

    arguments:
        sex_text (str): textual representation of sex selection
        alcohol_consumption (bool): whether to filter to persons who are past or
            present consumers of alcohol
        alcoholism (str): name of column defining alcoholism cases and controls
        alcoholism_split (str): how to divide cohort by alcoholism, either
            "all", "case", or "control"
        hormone (str): name of column for hormone
        variables_names_valid (list<str>): names of columns for variables in
            which rows must have valid values
        variables_prefixes_valid (list<str>): prefixes of columns for variables
            in which rows must have valid values
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=4)
    print("sex: " + str(sex_text))
    print("alcoholism: " + str(alcoholism))
    print("alcoholism split: " + str(alcoholism_split))
    print("hormone: " + str(hormone))
    # Select cohort's variables and records with valid values.
    table_valid = select_sex_alcoholism_cohort_variables_valid_records(
        sex_text=sex_text,
        alcohol_consumption=alcohol_consumption,
        alcoholism=alcoholism,
        alcoholism_split=alcoholism_split,
        variables_names_valid=variables_names_valid,
        variables_prefixes_valid=variables_prefixes_valid,
        table=table,
    )
    # Report.
    print("table shape: " + str(table_valid.shape))
    pass


def organize_report_cohorts_by_sex_alcoholism_split_hormone(
    table=None,
):
    """
    Organizes report of cohorts by combinations of sex, alcoholism definition,
    alcoholism split (all, case, control), and valid hormone values.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Iterate on cohort permutations.
    sexes = ["female", "male",]
    alcoholisms = [
        "alcoholism_1", "alcoholism_2", "alcoholism_3", "alcoholism_4",
        "alcoholism_5",
    ]
    alcoholism_splits = [
        "all", "case", "control",
    ]
    hormones = ["oestradiol", "testosterone",]
    for sex in sexes:
        utility.print_terminal_partition(level=2)
        for alcoholism in alcoholisms:
            for alcoholism_split in alcoholism_splits:
                utility.print_terminal_partition(level=3)
                organize_report_cohort_by_sex_alcoholism_split_hormone(
                    sex_text=sex,
                    alcohol_consumption=True,
                    alcoholism=alcoholism,
                    alcoholism_split=alcoholism_split,
                    hormone="None... not specific",
                    variables_names_valid=[
                        "eid", "IID",
                        "sex_y", "sex_x", "sex_text", "age", "body",
                        "alcohol_none",
                        alcoholism,
                    ],
                    variables_prefixes_valid=["genotype_pc_",],
                    table=table,
                )
                for hormone in hormones:
                    organize_report_cohort_by_sex_alcoholism_split_hormone(
                        sex_text=sex,
                        alcohol_consumption=True,
                        alcoholism=alcoholism,
                        alcoholism_split=alcoholism_split,
                        hormone=hormone,
                        variables_names_valid=[
                            "eid", "IID",
                            "sex_y", "sex_x", "sex_text", "age", "body",
                            "alcohol_none",
                            alcoholism,
                            hormone,
                        ],
                        variables_prefixes_valid=["genotype_pc_",],
                        table=table,
                    )
    pass


def organize_alcoholism_case_control_definitions(
    table=None,
    report=None,
):
    """
    Organizes information about alcoholism cases and controls.

    Both alcoholism cases and controls are persons who have consumed alcohol at
    some point ("ever" alcohol consumers as opposed to "never").

    Alcoholism controls do not qualify as cases by any of the main definitions.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about alcoholism cases and controls

    """

    # Copy information.
    table = table.copy(deep=True)
    # Specify diagnostic thresholds for AUDIT-C and AUDIT scores.
    # https://auditscreen.org/about/scoring-audit/
    # https://cde.drugabuse.gov/instrument/f229c68a-67ce-9a58-e040-bb89ad432be4
    # Use "less than" logic (record < threshold) for control thresholds.
    # Use "greater than or equal" logic (record >= threshold) for case
    # thresholds.
    # Be a bit more permissive with the AUDIT-C score to allow for persons who
    # consume a lot of alcohol without psychological problems.
    threshold_auditc_control = 5 # include AUDIT-C scores 0-4
    threshold_auditc_case = 10
    threshold_auditp_control = 4 # include AUDIT-P scores 0-3
    threshold_auditp_case = 5
    threshold_audit_control = 8 # include AUDIT scores 0-7
    threshold_audit_case = 15 # AUDIT score 15 or more "likely alcohol dep"

    # Determine qualification for control status for alcohol dependence.
    table["alcoholism_control"] = table.apply(
        lambda row:
            determine_alcoholism_control(
                alcohol_ever=row["alcohol_ever"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
                alcohol_audit=row["alcohol_audit"],
                threshold_auditc_control=threshold_auditc_control,
                threshold_auditp_control=threshold_auditp_control,
                threshold_audit_control=threshold_audit_control,
            ),
        axis="columns", # apply function to each row
    )

    # Determine qualification for case status for alcoholism type 1.
    table["alcoholism_case_1"] = table.apply(
        lambda row:
            determine_alcoholism_case_one(
                alcohol_ever=row["alcohol_ever"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine qualification for case status for alcoholism type 2.
    table["alcoholism_case_2"] = table.apply(
        lambda row:
            determine_alcoholism_case_two(
                alcohol_ever=row["alcohol_ever"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine qualification for case status for alcoholism type 3.
    table["alcoholism_case_3"] = table.apply(
        lambda row:
            determine_alcoholism_case_three(
                alcohol_ever=row["alcohol_ever"],
                alcohol_auditc=row["alcohol_auditc"],
                threshold_auditc_case=threshold_auditc_case,
            ),
        axis="columns", # apply function to each row
    )
    # Determine qualification for case status for alcoholism type 4.
    table["alcoholism_case_4"] = table.apply(
        lambda row:
            determine_alcoholism_case_four(
                alcohol_ever=row["alcohol_ever"],
                alcohol_auditp=row["alcohol_auditp"],
                threshold_auditp_case=threshold_auditp_case,
            ),
        axis="columns", # apply function to each row
    )
    # Determine qualification for case status for alcoholism type 5.
    table["alcoholism_case_5"] = table.apply(
        lambda row:
            determine_alcoholism_case_five(
                alcohol_ever=row["alcohol_ever"],
                alcohol_audit=row["alcohol_audit"],
                threshold_audit_case=threshold_audit_case,
            ),
        axis="columns", # apply function to each row
    )

    # Determine combination of all cases statuses for alcoholism.
    table["alcoholism_case_any"] = table.apply(
        lambda row:
            determine_alcoholism_case_any(
                alcoholism_case_1=row["alcoholism_case_1"],
                alcoholism_case_2=row["alcoholism_case_2"],
                alcoholism_case_3=row["alcoholism_case_3"],
                alcoholism_case_4=row["alcoholism_case_4"],
                alcoholism_case_5=row["alcoholism_case_5"],
            ),
        axis="columns", # apply function to each row
    )

    # Combine alcoholism case and control definitions.
    # name: "alcoholism_control_case_any"
    # NA: missing value indicates neither control nor case
    # 0: control
    # 1: case
    table["alcoholism_control_case_any"] = table.apply(
        lambda row:
            determine_alcoholism_control_case(
                alcoholism_control=row["alcoholism_control"],
                alcoholism_case=row["alcoholism_case_any"],
            ),
        axis="columns", # apply function to each row
    )
    table["alcoholism_control_case_1"] = table.apply(
        lambda row:
            determine_alcoholism_control_case(
                alcoholism_control=row["alcoholism_control"],
                alcoholism_case=row["alcoholism_case_1"],
            ),
        axis="columns", # apply function to each row
    )
    table["alcoholism_control_case_2"] = table.apply(
        lambda row:
            determine_alcoholism_control_case(
                alcoholism_control=row["alcoholism_control"],
                alcoholism_case=row["alcoholism_case_2"],
            ),
        axis="columns", # apply function to each row
    )


    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "alcohol_diagnosis_a", "alcohol_diagnosis_b",
            "alcohol_diagnosis_c", "alcohol_diagnosis_d",
            "alcohol_diagnosis_self",
            #"alcohol_auditc",
            #"alcohol_audit",
        ],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "alcohol_ever",
            "alcoholism_control",
            "alcoholism_case_any",
            "alcoholism_control_case_any",
            "alcoholism_control_case_1",
            "alcoholism_control_case_2",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcoholism cases and controls variables: ")
        print(table_report)
        #organize_report_cohorts_by_sex_alcoholism_split_hormone(
        #    table=table_clean,
        #)
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Smoking, use of tobacco smoke


# review: TCW; 08 August 2022
def interpret_tobacco_smoke_ever(
    field_20160=None,
):
    """
    Intepret UK Biobank's data-coding for data-field '20160'.

    Data-Field '20160': "Ever smoked"
    UK Biobank data-coding '7' for data-field '20160'.
    1: 'Yes'
    0: 'No'

    Accommodate inexact float values.

    arguments:
        field_20160 (float): UK Biobank data-field '20160', whether person ever
            used tobacco smoke

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_20160)) and
        (-0.5 <= field_20160 and field_20160 < 1.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (0.5 <= field_20160 and field_20160 < 1.5):
            # 1: "Yes"
            interpretation = 1
        elif (-0.5 <= field_20160 and field_20160 < 0.5):
            # 0: "No"
            interpretation = 0
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


# review: TCW; 08 August 2022
def interpret_tobacco_smoke_status(
    field_20116=None,
):
    """
    Intepret UK Biobank's data-coding for data-field '20116'.

    Data-Field '20116': "Smoking status"
    UK Biobank data-coding '90' for data-field '20116'.
    2: 'Current'
    1: 'Previous'
    0: 'Never'
    -3: 'Prefer not to answer'

    Accommodate inexact float values.

    arguments:
        field_20116 (float): UK Biobank data-field '20116', person's status of
            using tobacco smoke

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(field_20116)) and
        (-3.5 <= field_20116 and field_20116 < 2.5)
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (-0.5 <= field_20116 and field_20116 < 0.5):
            # 0: "Never"
            interpretation = 0
        elif (0.5 <= field_20116 and field_20116 < 1.5):
            # 1: "Previous"
            interpretation = 1
        elif (1.5 <= field_20116 and field_20116 < 2.5):
            # 2: "Current"
            interpretation = 2
        elif (-3.5 <= field_20116 and field_20116 < -2.5):
            # -3: "Prefer not to answer"
            interpretation = float("nan")
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def organize_tobacco_smoke_variables(
    table=None,
    report=None,
):
    """
    Organizes information about previous and current alcohol consumption.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert variable types.
    columns_type = [
        "20160-0.0", "20116-0.0", "1239-0.0", "1249-0.0",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Determine person's use of tobacco smoke.
    table["smoke_ever"] = table.apply(
        lambda row:
            interpret_tobacco_smoke_ever(
                field_20160=row["20160-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    table["smoke_status"] = table.apply(
        lambda row:
            interpret_tobacco_smoke_status(
                field_20116=row["20116-0.0"],
            ),
        axis="columns", # apply function to each row
    )

    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=["20160-0.0", "20116-0.0", "1239-0.0", "1249-0.0",],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "20160-0.0", "20116-0.0", "1239-0.0", "1249-0.0",
            "smoke_ever",
            "smoke_status",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of tobacco smoking variables: ")
        print(table_report)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail





##########
##########
##########
##########
##########
# Scrap...


def scrap_translate_alcohol_none_auditc(
    frequency=None,
):
    """
    Translates information from AUDIT-C questionnaire about whether person
    never consumes any alcohol.

    "frequency", field 20414: "Frequency of drinking alcohol"
    UK Biobank data coding 521 for variable field 20414.
    "four or more times a week": 4
    "two to three times a week": 3
    "two to four times a month": 2
    "monthly or less": 1
    "never": 0
    "prefer not to answer": -818

    Accommodate inexact float values.

    arguments:
        frequency (float): AUDIT-C question 1, UK Biobank field 20414

    raises:

    returns:
        (bool): whether person never consumes any alcohol

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(frequency)) and
        (-0.5 <= frequency and frequency < 4.5)
    ):
        # The variable has a valid value.
        # Determine whether the person ever consumes any alcohol.
        if (-0.5 <= frequency and frequency < 0.5):
            # Person never consumes any alcohol.
            alcohol_none_auditc = True
        else:
            # Persons consumes alcohol.
            alcohol_none_auditc = False
    else:
        alcohol_none_auditc = float("nan")
    # Return information.
    return alcohol_none_auditc


def scrap_determine_auditc_questionnaire_alcoholism_score(
    frequency=None,
    quantity=None,
    binge=None,
):
    """
    Determine aggregate alcoholism score from responses to the AUDIT-C
    questionnaire.

    "frequency", field 20414: "Frequency of drinking alcohol"
    UK Biobank data coding 521 for variable field 20414.
    "four or more times a week": 4
    "two to three times a week": 3
    "two to four times a month": 2
    "monthly or less": 1
    "never": 0
    "prefer not to answer": -818

    UK Biobank only asked questions "quantity" and "binge" if question
    "frequency" was not "never".

    "quantity", field 20403: "Amount of alcohol drunk on a typical drinking
    day"
    UK Biobank data coding 522 for variable field 20403.
    "ten or more": 5
    "seven, eight, or nine": 4
    "five or six": 3
    "three or four": 2
    "one or two": 1
    "prefer not to answer": -818

    "binge", field 20416: "Frequency of consuming six or more units of alcohol"
    UK Biobank data coding 523 for variable field 20416.
    "daily or almost daily": 5
    "weekly": 4
    "monthly": 3
    "less than monthly": 2
    "never": 1
    "prefer not to answer": -818

    Accommodate inexact float values.

    arguments:
        frequency (float): AUDIT-C question 1, UK Biobank field 20414
        quantity (float): AUDIT-C question 2, UK Biobank field 20403
        binge (float): AUDIT-C question 3, UK Biobank field 20416

    raises:

    returns:
        (float): person's monthly alcohol consumption in drinks

    """

    # Only consider cases with valid responses to all three questions.
    if (
        (
            (not pandas.isna(frequency)) and
            (-0.5 <= frequency and frequency < 4.5)
        ) and
        (
            (not pandas.isna(quantity)) and
            (0.5 <= quantity and quantity < 5.5)
        ) and
        (
            (not pandas.isna(binge)) and
            (0.5 <= binge and binge < 5.5)
        )
    ):
        score = (frequency + quantity + binge)
    else:
        score = float("nan")
    # Return information.
    return score


def scrap_organize_auditc_questionnaire_alcoholism_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcoholism from the AUDIT-C questionnaire.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)
    # Determine whether person never consumes any alcohol.
    table["alcohol_none_auditc"] = table.apply(
        lambda row:
            translate_alcohol_none_auditc(
                frequency=row["20414-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Determine alcoholism aggregate score.
    table["alcoholism"] = table.apply(
        lambda row:
            determine_auditc_questionnaire_alcoholism_score(
                frequency=row["20414-0.0"],
                quantity=row["20403-0.0"],
                binge=row["20416-0.0"],
            ),
        axis="columns", # apply function to each row
    )
    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "20414-0.0", "20403-0.0", "20416-0.0",
        ],
        axis="columns",
        inplace=True
    )
    # Organize data for report.
    table_report = table.copy(deep=True)
    table_report = table_report.loc[
        :, table_report.columns.isin([
            "eid", "IID",
            "20414-0.0", "20403-0.0", "20416-0.0",
            "alcohol_none_auditc",
            "alcoholism",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcohol consumption quantity variables: ")
        print(table_report)
    # Collect information.
    bucket = dict()
    bucket["table"] = table
    bucket["table_clean"] = table_clean
    bucket["table_report"] = table_report
    # Return information.
    return bucket


# End scrap section...
##########
##########
##########
##########
##########


##########
# Write


def write_product_organization(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_directory, "table_phenotypes.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        path_directory, "table_phenotypes.tsv"
    )
    # Write information to file.
    pail_write["table_phenotypes"].to_pickle(
        path_table_phenotypes
    )
    pail_write["table_phenotypes"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product(
    pail_write=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Organization procedure main information.
    write_product_organization(
        pail_write=pail_write["organization"],
        path_directory=paths["organization"],
    )
    pass


###############################################################################
# Drivers
# The purpose of these driver functions is to make the module's functionality
# accessible in modular parts.

##########
# Data export


# This function is mostly obsolete... but still might come in handy
def organize_hormone_female_export_table(
    table=None,
    select_columns=None,
    report=None,
):
    """
    Organizes information for export.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        select_columns (bool): whether to select specific columns from table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank

    """

    table = table.copy(deep=True)
    if (select_columns):
        columns_export = [
            #"eid",
            "IID",
            "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
            "pregnancy",
            "hysterectomy", "oophorectomy", "hysterectomy_or_oophorectomy",
            "menopause_binary", "menopause_ordinal",
            "menstruation_days", "menstruation_phase",
            "oral_contraception", "hormone_replacement",
            "hormone_alteration",
            "albumin", "albumin_log", "steroid_globulin", "steroid_globulin_log",
            "oestradiol", "oestradiol_log",
            "oestradiol_free", "oestradiol_free_log",
            "oestradiol_bioavailable", "oestradiol_bioavailable_log",
            "testosterone", "testosterone_log",
            "testosterone_free", "testosterone_free_log",
            "testosterone_bioavailable", "testosterone_bioavailable_log",

        ]
        table = table.loc[
            :, table.columns.isin(columns_export)
        ]
        table = table[[*columns_export]]
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: organize_hormone_export_table()")
        utility.print_terminal_partition(level=3)
        print(table)
        print("table shape (rows, columns): " + str(table.shape))
        print(table.columns.to_list())
    # Return information.
    return table


# Definitions


def execute_genotype_assessment_basis(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about persons' genotypes, sex, age, and body mass
    index across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Organize information about genotype principal components.
    table_genotype = organize_genotype_principal_component_variables(
        table=table,
        report=report,
    )
    # Organize information about general attributes.
    pail_basis = organize_assessment_basis_variables(
        table=table_genotype,
        path_dock=path_dock,
        report=report,
    )

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_genotype_assessment_basis()")
        utility.print_terminal_partition(level=3)
        print(pail_basis["table"])
    # Return information.
    return pail_basis


def execute_sex_hormones(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Organize information about sex hormones.
    pail_hormone = organize_sex_hormone_variables(
        table=table,
        path_dock=path_dock,
        report=report,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_sex_hormones()")
        utility.print_terminal_partition(level=3)
        print(pail_hormone["table_report"])
    # Return information.
    return pail_hormone


def execute_female_menstruation(
    table=None,
    report=None,
):
    """
    Organizes information about female persons' menstruation across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Organize information about female persons' pregnancy and menopause.
    pail_female = organize_female_menstruation_pregnancy_menopause_variables(
        table=table,
        report=report,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_female_menstruation()")
        utility.print_terminal_partition(level=3)
        print(pail_female["table_report"])
    # Return information.
    return pail_female


def execute_substance_use(
    table=None,
    report=None,
):
    """
    Organizes information about person's use of alcohol and tobacco smoke in
    UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank

    """

    # Organize information about alcohol consumption.
    pail_alcohol_consumption = organize_alcohol_consumption_variables(
        table=table,
        report=True,
    )
    #print(pail_alcohol_consumption["table_report"])

    # Organize Alchol Use Disorders Identification Test (AUDIT) questionnaire
    # variables, including separate scores for AUDIT-Consumption (AUDIT-C) and
    # AUDIT-Problem (AUDIT-P) portions of the questionnaire.
    pail_audit = organize_alcohol_audit_questionnaire_variables(
        table=pail_alcohol_consumption["table"],
        report=True,
    )
    #print(pail_audit["table_report"])

    # Organize International Classification of Disease (ICD) ICD9 and ICD10
    # codes for diagnoses relevant to alcoholism.
    # Organize codes for self diagnoses relevant to alcoholism.
    pail_diagnosis = organize_alcoholism_diagnosis_variables(
        table=pail_audit["table"],
        report=True,
    )
    #print(pail_diagnosis["table_report"])

    # Organize alcoholism cases and controls.
    # Report females and males who consume alcohol and are candidates for
    # either controls or cases of alcoholism.
    pail_alcoholism = organize_alcoholism_case_control_definitions(
        table=pail_diagnosis["table"],
        report=True,
    )
    #print(pail_alcoholism["table_clean"])

    # Organize information about use of tobacco smoke.
    pail_tobacco_smoke = organize_tobacco_smoke_variables(
        table=pail_alcoholism["table"],
        report=True,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: execute_substance_use()")
        utility.print_terminal_partition(level=3)
        print(pail_tobacco_smoke["table_report"])
    # Return information.
    return pail_tobacco_smoke


def execute_psychology_psychiatry(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about persons' mental health across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Organize information about neuroticism.
    pail_psychology = organize_psychology_variables(
        table=table,
        path_dock=path_dock,
        report=report,
    )
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_psychology_psychiatry()")
        utility.print_terminal_partition(level=3)
        print(pail_psychology["table_clean"])
    # Return information.
    return pail_psychology



###############################################################################
# Procedure


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
    print("version check: TCW, 18 August 2021")
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
        source="importation",
        path_dock=path_dock,
        report=True,
    )
    # Organize variables for persons' genotypes, sex, age, and body mass index
    # across the UK Biobank.
    pail_basis = execute_genotype_assessment_basis(
        table=source["table_phenotypes"],
        path_dock=path_dock,
        report=True,
    )

    # Organize variables for persons' sex hormones across the UK Biobank.
    pail_hormone = execute_sex_hormones(
        table=pail_basis["table"],
        path_dock=path_dock,
        report=True,
    )
    # Organize variables for female menstruation across the UK Biobank.
    pail_female = execute_female_menstruation(
        table=pail_hormone["table"],
        report=True,
    )
    # Organize variables for persons' alcohol consumption across the UK Biobank.
    pail_substance_use = execute_substance_use(
        table=pail_female["table"],
        report=True,
    )
    # Organize variables for persons' mental health across the UK Biobank.
    pail_psychology = execute_psychology_psychiatry(
        table=pail_substance_use["table"],
        path_dock=path_dock,
        report=True,
    )

    # Collect information.
    pail_write = dict()
    pail_write["organization"] = dict()
    pail_write["organization"]["table_phenotypes"] = pail_psychology["table"]
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )
    pass

if (__name__ == "__main__"):
    print("procedure needs path to a dock directory...")
    #execute_procedure()


#
