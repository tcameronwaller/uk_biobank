"""
Describe variables from data of the U.K. Biobank.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of UK_Biobank
    (https://github.com/tcameronwaller/uk_biobank/).

    UK_Biobank supports analyses on data from the U.K. Biobank.
    Copyright (C) 2021 Thomas Cameron Waller

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
import promiscuity.regression as pro_reg
import promiscuity.plot as plot
import uk_biobank.stratification as ukb_strat
import uk_biobank.organization as ukb_organization

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
    paths["description"] = os.path.join(path_dock, "description")
    paths["description_tables"] = os.path.join(paths["description"], "tables")
    paths["description_plots"] = os.path.join(paths["description"], "plots")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["description"])
    # Initialize directories.
    utility.create_directories(
        path=paths["description"]
    )
    utility.create_directories(
        path=paths["description_tables"]
    )
    utility.create_directories(
        path=paths["description_plots"]
    )
    # Return information.
    return paths


##########
# Read


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_dock, "organization",
        "table_phenotypes.pickle",
    )

    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        #"table_ukb_samples": table_ukb_samples,
    }


def read_source_regression_summary_tables(
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
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Collect tables.
    pail = dict()
    # Define path to parent directory.
    path_directory_parent = os.path.join(
        path_dock, "regression",
    )
    # Read all files within parent directory and organize tables.
    pail = (
        utility.read_all_pandas_tables_files_within_parent_directory(
            path_directory_parent=path_directory_parent,
            types_pandas_table_read={
                "cohort": "string",
                "dependence": "string",
                "dependence_type": "string",
                "model_adjustment": "string",
                "model_context": "string",
                "variable": "string",
                "parameter": "float32",
                "error": "float32",
                "interval_95": "float32",
                "range_95_below": "float32",
                "range_95_above": "float32",
            },
            report=report,
    ))
    # Return information.
    return pail


def read_source_correlation_summary_tables(
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
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Collect tables.
    pail = dict()
    # Define path to parent directory.
    path_directory_parent = os.path.join(
        path_dock, "parameters", "uk_biobank", "genetic_correlation",
    )
    # Read all files within parent directory and organize tables.
    pail = (
        utility.read_all_pandas_tables_files_within_parent_directory(
            path_directory_parent=path_directory_parent,
            types_pandas_table_read={
                "cohort": "string",
                "phenotype_primary": "string",
                "dependence_type": "string",
                "model_adjustment": "string",
                "model_context": "string",
                "variable": "string",
                "correlation": "float32",
                "standard_error": "float32",
                "interval_95": "float32",
                "confidence_95_low": "float32",
                "confidence_95_high": "float32",
            },
            report=report,
    ))
    # Return information.
    return pail





# TODO: TCW, 23 March 2022
# TODO: the genotype cohort tables don't include all phenotype variables...
# TODO: 1. read in the genotype cohort tables and then extract "eid" identifiers
# TODO: 2. read in the full "table_phenotypes" (or accept as argument)
# TODO: 3. subset "table_phenotypes" by "eid" identifiers from genotype cohort tables


def read_organize_cohorts(
    path_directory=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_directory (str): path to container directory for files of
            stratified cohort phenotype tables
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """


    # Collect names of child files within parent directory.
    files_names = utility.extract_directory_file_names(path=path_directory)
    # Filter names of child files by suffix.
    files_names_keep = list(filter(
        lambda file_name: (".pickle" in file_name),
        files_names
    ))

    # Collect cohort records.
    records = list()

    # Iterate on files within directory.
    for file_name in files_names_keep:
        # Determine name of cohort.
        name_cohort = str(file_name).replace(".pickle", "").replace("table_", "")
        # Specify directories and files.
        path_file = os.path.join(
            path_directory, file_name,
        )
        # Read information from file.
        table_cohort = pandas.read_pickle(
            path_file
        )

        # Collect information about cohort.
        record = dict()
        record["name"] = name_cohort
        record["cohort"] = name_cohort
        record["table"] = table_cohort
        records.append(record)
        pass

    # Return information.
    return records


##########
# Attribution Table
# Counts of records in each cohort with each categorical or discrete value of
# variables
##########


def define_variables_description_table_attribution():
    """
    Defines categorical or discrete variables and their values for description
    in attribution table.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Collect records of information about each categorical varialbe and value.
    records = list()

    # Variable: "white_british"

    record = dict()
    record["name"] = "white_british_0"
    record["variable"] = "white_british" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "white_british_1"
    record["variable"] = "white_british" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "genotype_availability"

    record = dict()
    record["name"] = "genotype_availability_0"
    record["variable"] = "genotype_availability" # cat. or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "genotype_availability_1"
    record["variable"] = "genotype_availability" # cat. or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "genotype_array"

    record = dict()
    record["name"] = "genotype_array_axiom"
    record["variable"] = "genotype_array" # cat. or discrete variable
    record["value"] = "axiom" # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "genotype_array_bileve"
    record["variable"] = "genotype_array" # cat. or discrete variable
    record["value"] = "bileve" # categorical or discrete value of variable
    records.append(record)

    # Variable: "sex_text"

    record = dict()
    record["name"] = "sex_text_female"
    record["variable"] = "sex_text" # categorical or discrete variable
    record["value"] = "female" # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "sex_text_male"
    record["variable"] = "sex_text" # categorical or discrete variable
    record["value"] = "male" # categorical or discrete value of variable
    records.append(record)

    # Variable: "alteration_sex_hormone"

    record = dict()
    record["name"] = "alteration_sex_hormone_0"
    record["variable"] = "alteration_sex_hormone" # cat. or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alteration_sex_hormone_1"
    record["variable"] = "alteration_sex_hormone" # cat. or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)


    # Variable: "pregnancy"

    record = dict()
    record["name"] = "pregnancy_0"
    record["variable"] = "pregnancy" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "pregnancy_1"
    record["variable"] = "pregnancy" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)


    # Variable: "births_any"

    record = dict()
    record["name"] = "births_any_0"
    record["variable"] = "births_any" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "births_any_1"
    record["variable"] = "births_any" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)


    # Variable: "birth_live_recent"

    record = dict()
    record["name"] = "birth_live_recent_0"
    record["variable"] = "birth_live_recent" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "birth_live_recent_1"
    record["variable"] = "birth_live_recent" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)


    # Variable: "hysterectomy"

    record = dict()
    record["name"] = "hysterectomy_0"
    record["variable"] = "hysterectomy" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "hysterectomy_1"
    record["variable"] = "hysterectomy" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "oophorectomy"

    record = dict()
    record["name"] = "oophorectomy_0"
    record["variable"] = "oophorectomy" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "oophorectomy_1"
    record["variable"] = "oophorectomy" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "hysterectomy_or_oophorectomy"

    record = dict()
    record["name"] = "hysterectomy_or_oophorectomy_0"
    record["variable"] = "hysterectomy_or_oophorectomy"
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "hysterectomy_or_oophorectomy_1"
    record["variable"] = "hysterectomy_or_oophorectomy"
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "menstruation_regular_range"

    record = dict()
    record["name"] = "menstruation_regular_range_0"
    record["variable"] = "menstruation_regular_range"
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "menstruation_regular_range_1"
    record["variable"] = "menstruation_regular_range"
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "menopause_ordinal"

    record = dict()
    record["name"] = "menopause_ordinal_0"
    record["variable"] = "menopause_ordinal" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "menopause_ordinal_1"
    record["variable"] = "menopause_ordinal" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "menopause_ordinal_2"
    record["variable"] = "menopause_ordinal" # categorical or discrete variable
    record["value"] = 2 # categorical or discrete value of variable
    records.append(record)

    # Variable: "alcohol_ever"

    record = dict()
    record["name"] = "alcohol_ever_0"
    record["variable"] = "alcohol_ever" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcohol_ever_1"
    record["variable"] = "alcohol_ever" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "alcohol_frequency"

    record = dict()
    record["name"] = "alcohol_frequency_0"
    record["variable"] = "alcohol_frequency" # categorical or discrete variable
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency_1"
    record["variable"] = "alcohol_frequency" # categorical or discrete variable
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency_2"
    record["variable"] = "alcohol_frequency" # categorical or discrete variable
    record["value"] = 2 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency_3"
    record["variable"] = "alcohol_frequency" # categorical or discrete variable
    record["value"] = 3 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency_4"
    record["variable"] = "alcohol_frequency" # categorical or discrete variable
    record["value"] = 4 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency_5"
    record["variable"] = "alcohol_frequency" # categorical or discrete variable
    record["value"] = 5 # categorical or discrete value of variable
    records.append(record)

    # Variable: "alcoholism_control_case_1"

    record = dict()
    record["name"] = "alcoholism_control_case_1_0"
    record["variable"] = "alcoholism_control_case_1"
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcoholism_control_case_1_1"
    record["variable"] = "alcoholism_control_case_1"
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "alcoholism_control_case_any"

    record = dict()
    record["name"] = "alcoholism_control_case_any_0"
    record["variable"] = "alcoholism_control_case_any"
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "alcoholism_control_case_any_1"
    record["variable"] = "alcoholism_control_case_any"
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "depression_control_case_loose"

    record = dict()
    record["name"] = "depression_control_case_loose_0"
    record["variable"] = "depression_control_case_loose"
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "depression_control_case_loose_1"
    record["variable"] = "depression_control_case_loose"
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Variable: "depression_control_case_strict"

    record = dict()
    record["name"] = "depression_control_case_strict_0"
    record["variable"] = "depression_control_case_strict"
    record["value"] = 0 # categorical or discrete value of variable
    records.append(record)

    record = dict()
    record["name"] = "depression_control_case_strict_1"
    record["variable"] = "depression_control_case_strict"
    record["value"] = 1 # categorical or discrete value of variable
    records.append(record)

    # Return information
    return records


def organize_attribution_record(
    name_cohort=None,
    name_variable_value=None,
    variable=None,
    value=None,
    table=None,
):
    """
    Organize a record (single row in table) to describe attribution of
    categorical or discrete variable values across cohorts.

    arguments:
        name_cohort (str): name of cohort
        name_variable_value (str): name of variable's value for report
        variable (str): name of table's column for variable
        value (object): categorical or discrete value of variable
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable_value"] = str(name_variable_value)
    record["variable"] = str(variable)
    record["value"] = str(value)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    table_variable_value = table.loc[
        (
            (~pandas.isna(table[variable])) &
            (table[variable] == value)
        ), :
    ]

    # Count records.
    count_total = table.shape[0]
    count_variable_value = table_variable_value.shape[0]

    # Calculate percentages.
    if (count_total > 0):
        percentage_variable_value = round(
            ((count_variable_value / count_total) * 100), 3
        )
    else:
        percentage_variable_value = float("nan")
        pass

    # Collect information for record.
    record["count_cohort_total_records"] = count_total
    record["count_variable_value"] = count_variable_value
    record["count_variable_value_report"] = str(
        str(count_variable_value) +
        " (" + str(percentage_variable_value) + "%)"
    )
    # Return information.
    return record


def organize_description_table_attribution(
    records_cohorts=None,
    report=None,
):
    """
    Organizes a description table for attribution of categorical or discrete
    variable values across cohorts.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Define variables.
    records_variables = define_variables_description_table_attribution()

    # Collect summary records for rows within description table.
    records_description = list()
    # Iterate on cohorts.
    for record_cohort in records_cohorts:
        # Iterate on variables.
        for record_variable in records_variables:
            # Organize information for description record.
            record_description = organize_attribution_record(
                name_cohort=record_cohort["name"],
                name_variable_value=record_variable["name"],
                variable=record_variable["variable"],
                value=record_variable["value"],
                table=record_cohort["table"],
            )
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table_description = pandas.DataFrame(data=records_description)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_attribution()")
        utility.print_terminal_partition(level=3)
        print(table_description)
        pass
    # Return information.
    return table_description


##########
# Missingness Table
# Measurement missingness percentages across cohorts
##########


def organize_missingness_record_basis(
    name_cohort=None,
    name_variable=None,
    column_measurement=None,
    column_detection=None,
    column_missingness_range=None,
    column_reportability_limit=None,
    table=None,
):
    """
    Determine for each variable whether missingness indicates that the
    measurement was not reportable because it was beyond the detection range.

    "[variable]_missingness_range" 1: missing due to measurement beyond
    detection range

    "[variable]_reportability_limit" 1: not reportable due to measurement less
    than limit of detection
    "[variable]_reportability_limit" 2: not reportable due to measurement
    greater than limit of detection

    arguments:
        name_cohort (str): name of cohort
        name_variable (str): name of variable for report
        column_measurement (str): name of table's column for variable's
            measurements
        column_detection (str): name of table's column for measurement's
            detectability
        column_missingness_range (str): name of table's column for measurement's
            reason for missingness beyond detection range
        column_reportability_limit (str): name of table's column for
            measurement's reportability less than limit of detection
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable"] = str(name_variable)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    #array_measurement_total = copy.deepcopy(table[name].to_numpy())
    #array_measurement_valid = copy.deepcopy(table[name].dropna().to_numpy())
    table_not_missing = table.dropna(
        axis="index",
        how="any",
        subset=[column_measurement],
        inplace=False,
    )
    #table_missing = table[table[column_variable].isna()]
    table_missing = table.loc[
        (
            (pandas.isna(table[column_measurement]))
        ), :
    ]
    # Count records.
    #count_measurement_total = int(array_measurement_total.size)
    #count_measurement_valid = int(array_measurement_valid.size)
    count_total = table.shape[0]
    count_not_missing = table_not_missing.shape[0]
    count_missing = table_missing.shape[0]
    count_missing_range = float("nan")
    count_reportability_low = float("nan")
    count_reportability_high = float("nan")
    count_detectable = float("nan")
    count_undetectable = float("nan")
    count_undetectable_low = float("nan")
    count_undetectable_high = float("nan")

    # Calculate percentages.
    if (count_total > 0):
        percentage_not_missing = round(
            ((count_not_missing / count_total) * 100), 3
        )
        percentage_missing = round(
            ((count_missing / count_total) * 100), 3
        )
    else:
        percentage_not_missing = float("nan")
        percentage_missing = float("nan")
        pass
    percentage_missing_range = float("nan")
    percentage_reportability_low = float("nan")
    percentage_reportability_high = float("nan")
    percentage_detectable = float("nan")
    percentage_undetectable = float("nan")
    percentage_undetectable_low = float("nan")
    percentage_undetectable_high = float("nan")

    # Collect information for record.
    record["count_cohort_total_records"] = count_total
    record["not_missing"] = str(
        str(count_not_missing) +
        " (" + str(percentage_not_missing) + "%)"
    )
    record["missing"] = str(
        str(count_missing) +
        " (" + str(percentage_missing) + "%)"
    )
    record["missing_range"] = str(
        str(count_missing_range) +
        " (" + str(percentage_missing_range) + "%)"
    )
    record["unreportable_low"] = str(
        str(count_reportability_low) +
        " (" + str(percentage_reportability_low) + "%)"
    )
    record["unreportable_high"] = str(
        str(count_reportability_high) +
        " (" + str(percentage_reportability_high) + "%)"
    )
    record["detectable"] = str(
        str(count_detectable) +
        " (" + str(percentage_detectable) + "%)"
    )
    record["undetectable"] = str(
        str(count_undetectable) +
        " (" + str(percentage_undetectable) + "%)"
    )
    record["undetectable_low"] = str(
        str(count_undetectable_low) +
        " (" + str(percentage_undetectable_low) + "%)"
    )
    record["undetectable_high"] = str(
        str(count_undetectable_high) +
        " (" + str(percentage_undetectable_high) + "%)"
    )
    # Return information.
    return record


def organize_missingness_record_description(
    name_cohort=None,
    name_variable=None,
    column_measurement=None,
    column_detection=None,
    column_missingness_range=None,
    column_reportability_limit=None,
    table=None,
):
    """
    Determine for each variable whether missingness indicates that the
    measurement was not reportable because it was beyond the detection range.

    "[variable]_missingness_range" 1: missing due to measurement beyond
    detection range

    "[variable]_reportability_limit" 1: not reportable due to measurement less
    than limit of detection
    "[variable]_reportability_limit" 2: not reportable due to measurement
    greater than limit of detection

    arguments:
        name_cohort (str): name of cohort
        name_variable (str): name of variable for report
        column_measurement (str): name of table's column for variable's
            measurements
        column_detection (str): name of table's column for measurement's
            detectability
        column_missingness_range (str): name of table's column for measurement's
            reason for missingness beyond detection range
        column_reportability_limit (str): name of table's column for
            measurement's reportability less than limit of detection
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable"] = str(name_variable)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    #array_measurement_total = copy.deepcopy(table[name].to_numpy())
    #array_measurement_valid = copy.deepcopy(table[name].dropna().to_numpy())
    table_not_missing = table.dropna(
        axis="index",
        how="any",
        subset=[column_measurement],
        inplace=False,
    )
    #table_missing = table[table[column_variable].isna()]
    table_missing = table.loc[
        (
            (pandas.isna(table[column_measurement]))
        ), :
    ]
    table_missing_range = table.loc[
        (
            #(pandas.isna(table[column_variable]) &
            (table[column_missingness_range] == 1)
        ), :
    ]
    table_reportability_low = table.loc[
        (
            #(pandas.isna(table[column_variable]) &
            (table[column_reportability_limit] == 1)
        ), :
    ]
    table_reportability_high = table.loc[
        (
            #(pandas.isna(table[column_variable]) &
            (table[column_reportability_limit] == 2)
        ), :
    ]
    table_detectable = table.loc[
        (table[column_detection] == 1), :
    ]
    table_undetectable = table.loc[
        (table[column_detection] == 0), :
    ]
    table_undetectable_low = table.loc[
        (
            (table[column_detection] == 0) &
            (table[column_reportability_limit] == 1)
        ), :
    ]
    table_undetectable_high = table.loc[
        (
            (table[column_detection] == 0) &
            (table[column_reportability_limit] == 2)
        ), :
    ]

    # Count records.
    #count_measurement_total = int(array_measurement_total.size)
    #count_measurement_valid = int(array_measurement_valid.size)
    count_total = table.shape[0]
    count_not_missing = table_not_missing.shape[0]
    count_missing = table_missing.shape[0]
    count_missing_range = table_missing_range.shape[0]
    count_reportability_low = table_reportability_low.shape[0]
    count_reportability_high = table_reportability_high.shape[0]
    count_detectable = table_detectable.shape[0]
    count_undetectable = table_undetectable.shape[0]
    count_undetectable_low = table_undetectable_low.shape[0]
    count_undetectable_high = table_undetectable_high.shape[0]

    # Calculate percentages.
    if (count_total > 0):
        percentage_not_missing = round(
            ((count_not_missing / count_total) * 100), 3
        )
        percentage_missing = round(
            ((count_missing / count_total) * 100), 3
        )
        percentage_missing_range = round(
            ((count_missing_range / count_total) * 100), 3
        )
        percentage_reportability_low = round(
            ((count_reportability_low / count_total) * 100), 3
        )
        percentage_reportability_high = round(
            ((count_reportability_high / count_total) * 100), 3
        )
        percentage_detectable = round(
            ((count_detectable / count_total) * 100), 3
        )
        percentage_undetectable = round(
            ((count_undetectable / count_total) * 100), 3
        )
        percentage_undetectable_low = round(
            ((count_undetectable_low / count_total) * 100), 3
        )
        percentage_undetectable_high = round(
            ((count_undetectable_high / count_total) * 100), 3
        )
    else:
        percentage_not_missing = float("nan")
        percentage_missing = float("nan")
        percentage_missing_range = float("nan")
        percentage_reportability_low = float("nan")
        percentage_reportability_high = float("nan")
        percentage_detectable = float("nan")
        percentage_undetectable = float("nan")
        percentage_undetectable_low = float("nan")
        percentage_undetectable_high = float("nan")
        pass

    # Collect information for record.
    record["count_cohort_total_records"] = count_total
    record["not_missing"] = str(
        str(count_not_missing) +
        " (" + str(percentage_not_missing) + "%)"
    )
    record["missing"] = str(
        str(count_missing) +
        " (" + str(percentage_missing) + "%)"
    )
    record["missing_range"] = str(
        str(count_missing_range) +
        " (" + str(percentage_missing_range) + "%)"
    )
    record["unreportable_low"] = str(
        str(count_reportability_low) +
        " (" + str(percentage_reportability_low) + "%)"
    )
    record["unreportable_high"] = str(
        str(count_reportability_high) +
        " (" + str(percentage_reportability_high) + "%)"
    )
    record["detectable"] = str(
        str(count_detectable) +
        " (" + str(percentage_detectable) + "%)"
    )
    record["undetectable"] = str(
        str(count_undetectable) +
        " (" + str(percentage_undetectable) + "%)"
    )
    record["undetectable_low"] = str(
        str(count_undetectable_low) +
        " (" + str(percentage_undetectable_low) + "%)"
    )
    record["undetectable_high"] = str(
        str(count_undetectable_high) +
        " (" + str(percentage_undetectable_high) + "%)"
    )
    # Return information.
    return record


def organize_missingness_record(
    name_cohort=None,
    name_variable=None,
    column_measurement=None,
    column_detection=None,
    column_missingness_range=None,
    column_reportability_limit=None,
    table=None,
):
    """
    Determine for each variable whether missingness indicates that the
    measurement was not reportable because it was beyond the detection range.

    "[variable]_missingness_range" 1: missing due to measurement beyond
    detection range

    "[variable]_reportability_limit" 1: not reportable due to measurement less
    than limit of detection
    "[variable]_reportability_limit" 2: not reportable due to measurement
    greater than limit of detection

    arguments:
        name_cohort (str): name of cohort
        name_variable (str): name of variable for report
        column_measurement (str): name of table's column for variable's
            measurements
        column_detection (str): name of table's column for measurement's
            detectability
        column_missingness_range (str): name of table's column for measurement's
            reason for missingness beyond detection range
        column_reportability_limit (str): name of table's column for
            measurement's reportability less than limit of detection
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    if (
        (column_detection in table.columns.to_list()) and
        (column_missingness_range in table.columns.to_list()) and
        (column_reportability_limit in table.columns.to_list())
    ):
        record = organize_missingness_record_description(
            name_cohort=name_cohort,
            name_variable=name_variable,
            column_measurement=column_measurement,
            column_detection=column_detection,
            column_missingness_range=column_missingness_range,
            column_reportability_limit=column_reportability_limit,
            table=table,
        )
    else:
        record = organize_missingness_record_basis(
            name_cohort=name_cohort,
            name_variable=name_variable,
            column_measurement=column_measurement,
            column_detection=column_detection,
            column_missingness_range=column_missingness_range,
            column_reportability_limit=column_reportability_limit,
            table=table,
        )

    # Return information.
    return record


def organize_description_table_missingness(
    records_cohorts=None,
    report=None,
):
    """
    Organizes a description table for missingness of measurement variables
    across cohorts.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Define variables.
    variables = list()
    variables.append("albumin")
    variables.append("steroid_globulin")
    variables.append("cholesterol")
    variables.append("vitamin_d")
    variables.append("oestradiol")
    variables.append("oestradiol_bioavailable")
    variables.append("oestradiol_free")
    variables.append("testosterone")
    variables.append("testosterone_bioavailable")
    variables.append("testosterone_free")

    # Collect summary records for rows within description table.
    records_description = list()
    # Iterate on cohorts.
    for record_cohort in records_cohorts:
        # Iterate on variables.
        for variable in variables:
            # Determine column names.
            detection = str(str(variable) + "_detection")
            missingness_range = str(str(variable) + "_missingness_range")
            reportability_limit = str(str(variable) + "_reportability_limit")
            # Organize information for description record.
            record_description = organize_missingness_record(
                name_cohort=record_cohort["name"],
                name_variable=variable,
                column_measurement=variable,
                column_detection=detection,
                column_missingness_range=missingness_range,
                column_reportability_limit=reportability_limit,
                table=record_cohort["table"],
            )
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table_description = pandas.DataFrame(data=records_description)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_missingness()")
        utility.print_terminal_partition(level=3)
        print(table_description)
        pass
    # Return information.
    return table_description


##########
# Threshold Table
# Measurement threshold percentages across cohorts
##########


# TODO: TCW, 08 March 2022
# TODO: try to make this fairly versatile... the idea is supply a "column_name" and "threshold_low" and then report
# for the main cohort tables from "stratification" procedure...




def organize_cohort_hormone_deficiency_record(
    name_cohort=None,
    name_variable=None,
    column_variable=None,
    threshold=None,
    table=None,
):
    """
    Determine for each hormone whether missingness variable indicates that the
    measurement was not reportable because it was beyond the detection range.

    # cohort:
    # hormone:
    # total measurement missingness: <count> (<percentage>%)
    # missing with reason "above or below reportable limit": <count> (<percentage>%)
    # missing with reportability "not reportable ... too low": <count> (<percentage>%)
    # missing with both annotations: <count> (<percentage>%)

    "[variable]_missingness_range" 1: missing due to measurement beyond
    detection range

    "[variable]_reportability_limit" 1: not reportable due to measurement less
    than limit of detection
    "[variable]_reportability_limit" 2: not reportable due to measurement
    greater than limit of detection

    arguments:
        name_cohort (str): name of cohort
        name_variable (str): name of variable for report
        column_variable (str): name of table's column for variable's
            measurements
        threshold (float): threshold below which constitutes deficiency
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable"] = str(name_variable)
    record["column"] = str(column_variable)
    record["threshold"] = str(threshold)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    table_below_threshold = table.loc[
        (
            (~pandas.isna(table[column_variable])) &
            (table[column_variable] < threshold)
        ), :
    ]
    table_above_threshold = table.loc[
        (
            (~pandas.isna(table[column_variable])) &
            (table[column_variable] >= threshold)
        ), :
    ]

    # Count records.
    #count_measurement_total = int(array_measurement_total.size)
    #count_measurement_valid = int(array_measurement_valid.size)
    count_total = table.shape[0]
    count_below_threshold = table_below_threshold.shape[0]
    count_above_threshold = table_above_threshold.shape[0]

    # Calculate percentages.
    if (count_total > 0):
        percentage_below_threshold = round(
            ((count_below_threshold / count_total) * 100), 3
        )
        percentage_above_threshold = round(
            ((count_above_threshold / count_total) * 100), 3
        )
    else:
        percentage_below_threshold = float("nan")
        percentage_above_threshold = float("nan")
        pass

    # Collect information for record.
    record["count_cohort_samples"] = count_total
    record["below_threshold"] = str(
        str(count_below_threshold) +
        " (" + str(percentage_below_threshold) + "%)"
    )
    record["above_threshold"] = str(
        str(count_above_threshold) +
        " (" + str(percentage_above_threshold) + "%)"
    )
    # Return information.
    return record


def organize_cohorts_hormone_deficiency(
    name_variable=None,
    column_variable=None,
    threshold=None,
    table=None,
    report=None,
):
    """
    Organizes a summary table about missingness of hormone measurements in
    cohorts.

    arguments:
        name_variable (str): name variable of interest
        column_variable (str): name of column in table for hormone of interest
        threshold (float): threshold below which constitutes deficiency
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Copy information.
    table = table.copy(deep=True)

    # Define cohorts for description.

    # TODO: switch to cohorts that stratify by season...


    #records_cohorts = ukb_strat.stratify_phenotype_cohorts_regression(
    #    table=table,
    #)

    records_cohorts = (
        ukb_strat.stratify_phenotype_cohorts_season_sex_age_menopause(
            table=table,
    ))
    # Collect summary records and construct table.
    records = list()
    # Iterate on cohorts.
    for collection_cohort in records_cohorts:
        # Organize information in record.
        record = organize_cohort_hormone_deficiency_record(
            name_cohort=collection_cohort["name"],
            name_variable=name_variable,
            column_variable=column_variable,
            threshold=threshold,
            table=collection_cohort["table"],
        )
        # Collect records.
        records.append(record)
        pass
    # Organize table.
    table_deficiency = pandas.DataFrame(data=records)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_cohorts_hormone_deficiency()")
        utility.print_terminal_partition(level=3)
        print(table_deficiency)
        pass
    # Return information.
    return table_deficiency


##########
# Quantitation Table
# Cohort, model, phenotype description
##########


def organize_report_contingency_table_stratification_by_missingness(
    table=None,
    report=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)

    # Organize contingency table descriptions of missingness.

    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="sex_text",
        stratifications=["female", "male"],
        column_missingness="testosterone",
        table=table,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="sex_text",
        stratifications=["female", "male"],
        column_missingness="oestradiol",
        table=table,
        report=report,
    )
    # Filter to females who were not pregnant.
    table_female = table.loc[
        (table["sex_text"] == "female"), :
    ]
    table_female_not_pregnant = table_female.loc[
        (table_female["pregnancy"] == 0), :
    ]
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 1], # premenopause versus perimenopause
        column_missingness="testosterone",
        table=table_female_not_pregnant,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 1], # premenopause versus perimenopause
        column_missingness="oestradiol",
        table=table_female_not_pregnant,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 2], # premenopause versus postmenopause
        column_missingness="testosterone",
        table=table_female_not_pregnant,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 2], # premenopause versus postmenopause
        column_missingness="oestradiol",
        table=table_female_not_pregnant,
        report=report,
    )

    # TODO: also consider missingness in "younger" and "older" males
    table_male = table.loc[
        (table["sex_text"] == "male"), :
    ]
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="age_grade_male",
        stratifications=[0, 2], # young versus old
        column_missingness="testosterone",
        table=table_male,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="age_grade_male",
        stratifications=[0, 2], # young versus old
        column_missingness="oestradiol",
        table=table_male,
        report=report,
    )

    pass


# TODO: I need to test this...
def organize_text_report_column_pair_correlations(
    column_one=None,
    column_two=None,
    table=None,
):
    """
    Organizes information about previous and current alcohol consumption.

    arguments:
        column_one (str): name of first column
        column_two (str): name of second column
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (str): textual report

    """

    table = table.copy(deep=True)
    table.dropna(
        axis="index",
        how="any",
        subset=[column_one, column_two],
        inplace=True,
    )
    count_rows = table.shape[0]
    pearson_correlation, pearson_probability = scipy.stats.pearsonr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    spearman_correlation, spearman_probability = scipy.stats.spearmanr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    # Report.
    report_text = str()
    report_text += textwrap.dedent(
        """\

            --------------------------------------------------
            column one: {column_one}
            column two: {column_two}
            non-missing records: {count_rows}
            ----------
            Correlations
            ----------
            Pearson: {pearson_correlation} (p: {pearson_probability})
            Spearman: {spearman_correlation} (p: {spearman_probability})
            --------------------------------------------------
        """
    ).format(
        column_one=column_one,
        column_two=column_two,
        count_rows=count_rows,
        pearson_correlation=round(pearson_correlation, 3),
        pearson_probability=round(pearson_probability, 5),
        spearman_correlation=round(spearman_correlation, 3),
        spearman_probability=round(spearman_probability, 5),
    )
    # Return information.
    return report_text


# TODO: work on these text reports...
def organize_text_report_phenotype_correlations(
    table=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Organize textual report.
    report_text = str()
    report_text += organize_text_report_column_pair_correlations(
        column_one="testosterone",
        column_two="steroid_globulin",
        table=None,
    )

    # TODO: okay so now I need to append that new text report to the collection text string report...


    # Return information.
    return report_text


# TODO: TCW 30 September 2021
# TODO: this function has problems...
def organize_report_variables_summaries_record_hormone_cohort_ordinal(
    record=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        record (dict): information for summary table record on cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for cohort-specific ordinal representations of
    # hormones.
    hormones = [
        "albumin", "steroid_globulin", "vitamin_d",
        "oestradiol", "testosterone",
    ]
    cohorts = [
        "female", "female_age_low", "female_age_middle", "female_age_high",
        "male", "male_age_low", "male_age_middle", "male_age_high",
        "female_premenopause", "female_perimenopause", "female_postmenopause",
    ]
    # Iterate across cohorts.
    for cohort in cohorts:
        # Iterate across hormones.
        for hormone in hormones:
            column_ordinal = str(str(hormone) + "_" + str(cohort) + "_order")
            values = copy.deepcopy(table[column_ordinal].dropna().to_list())
            if (len(values) > 10):
                count_values_0 = len(list(filter(
                    lambda value: (value == 0),
                    values
                )))
                count_values_1 = len(list(filter(
                    lambda value: (value == 1),
                    values
                )))
                count_values_2 = len(list(filter(
                    lambda value: (value == 2),
                    values
                )))
                count_values_3 = len(list(filter(
                    lambda value: (value == 3),
                    values
                )))
            else:
                count_values_0 = 0
                count_values_1 = 0
                count_values_2 = 0
                count_values_3 = 0
            # Collect information for record.
            record[str(column_ordinal + "_0")] = str(count_values_0)
            record[str(column_ordinal + "_1")] = str(count_values_1)
            record[str(column_ordinal + "_2")] = str(count_values_2)
            record[str(column_ordinal + "_3")] = str(count_values_3)
            pass
        pass

    # Return information.
    return record


def organize_quantitation_record(
    name_cohort=None,
    variable=None,
    table=None,
):
    """
    Organize a record (single row in table) to describe for measures on
    quantitative variables across cohorts.

    arguments:
        name_cohort (str): name of cohort
        variable (str): name of table's column for variable
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable"] = str(variable)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    table_genotype = table.loc[
        (
            (~pandas.isna(table["genotype_availability"])) &
            (table["genotype_availability"] == 1)
        ), :
    ]

    # Count records.
    count_total = int(table.shape[0])
    count_total_genotype = int(table_genotype.shape[0])

    # Initialize missing values.
    count_variable = float("nan")
    count_genotype_variable = float("nan")
    mean = float("nan")
    standard_error = float("nan")
    interval_95 = float("nan")
    confidence_95_low = float("nan")
    confidence_95_high = float("nan")
    median = float("nan")
    standard_deviation = float("nan")
    minimum = float("nan")
    maximum = float("nan")
    # Determine whether table has the column.
    if (variable in table.columns.to_list()):
        array = copy.deepcopy(table[variable].dropna().to_numpy())
        array_genotype = copy.deepcopy(
            table_genotype[variable].dropna().to_numpy()
        )
        # Determine count of valid values.
        count_variable = int(array.size)
        count_genotype_variable = int(array_genotype.size)
        if (count_variable > 10):
            # Determine mean, median, standard deviation, and standard error of
            # values in array.
            mean = numpy.nanmean(array)
            standard_error = scipy.stats.sem(array)
            interval_95 = (1.96 * standard_error)
            confidence_95_low = (mean - interval_95)
            confidence_95_high = (mean + interval_95)
            median = numpy.nanmedian(array)
            standard_deviation = numpy.nanstd(array)
            minimum = numpy.nanmin(array)
            maximum = numpy.nanmax(array)
            pass
        pass
    # Collect information for record.
    record["count_cohort_total_records"] = count_total
    record["count_cohort_total_genotypes"] = count_total_genotype
    record["count_variable_non_missing"] = str(count_variable)
    record["count_variable_genotypes"] = str(count_genotype_variable)
    record["mean"] = str(round(mean, 7))
    record["standard_error"] = str(round(standard_error, 7))
    record["interval_95"] = str(round(interval_95, 7))
    record["confidence_95_low"] = str(round(confidence_95_low, 7))
    record["confidence_95_high"] = str(round(confidence_95_high, 7))
    record["range_confidence_95"] = str(
        str(round(confidence_95_low, 3)) + " ... " +
        str(round(confidence_95_high, 3))
    )
    record["median"] = str(round(median, 7))
    record["standard_deviation"] = str(round(standard_deviation, 7))
    record["minimum"] = str(round(minimum, 7))
    record["maximum"] = str(round(maximum, 7))

    # Return information.
    return record


def organize_description_table_quantitation(
    records_cohorts=None,
    report=None,
):
    """
    Organizes a description table for measures on quantitative variables. Most
    relevant for quantitative variables of Ratio Scale, but also informative for
    Interval and Ordinal Scales.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Define variables.
    variables = [
        "age", "body",
        "albumin", "albumin_imputation",
        "steroid_globulin", "steroid_globulin_imputation",
        "cholesterol", "cholesterol_imputation",
        "vitamin_d", "vitamin_d_imputation",
        "oestradiol", "oestradiol_imputation",
        "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
        "oestradiol_free", "oestradiol_free_imputation",
        "testosterone", "testosterone_imputation",
        "testosterone_bioavailable", "testosterone_bioavailable_imputation",
        "testosterone_free", "testosterone_free_imputation",
        "menstruation_days", "menstruation_duration",
        "pregnancies", "pregnancies_early", "births",
        "neuroticism",
        "alcohol_frequency",
        "alcohol_drinks_weekly", "alcohol_drinks_monthly",
        "alcohol_drinks_monthly_combination",
        "alcohol_auditc", "alcohol_auditp", "alcohol_audit",
    ]

    # Collect summary records for rows within description table.
    records_description = list()
    # Iterate on cohorts.
    for record_cohort in records_cohorts:
        # Iterate on variables.
        for variable in variables:
            # Organize information for description record.
            record_description = organize_quantitation_record(
                name_cohort=record_cohort["name"],
                variable=variable,
                table=record_cohort["table"],
            )
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table_description = pandas.DataFrame(data=records_description)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_quantitation()")
        utility.print_terminal_partition(level=3)
        print(table_description)
        pass
    # Return information.
    return table_description


##########
# Plot
##########


def plot_variable_values_histogram(
    label_title=None,
    array=None,
    bins=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        label_title (str): label title for plot
        array (object): NumPy array of values to bin and plot in histogram
        bins (int): count of bins for histogram

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()
    # Determine bin method.
    if bins is None:
        bin_method = "auto"
    else:
        bin_method = "count"
    # Create figure.
    figure = plot.plot_distribution_histogram(
        array=array,
        title="",
        bin_method=bin_method, # "auto" or "count"
        bin_count=bins,
        label_bins="values",
        label_counts="counts per bin",
        fonts=fonts,
        colors=colors,
        line=True,
        line_position=numpy.nanmean(array),
        label_title=label_title, # ""
        label_report=True,
    )
    # Return.
    return figure


def plot_variable_means_dot_trajectory(
    label_title=None,
    column_phenotype=None,
    column_trajectory=None,
    threshold_trajectory=None,
    title_abscissa=None,
    title_ordinate=None,
    table=None,
):
    """
    Plots charts from the analysis process.

    The trajectory variable must have discrete values.

    arguments:
        label_title (str): label title for plot
        column_phenotype (str): name of column in table for continuous phenotype
        column_trajectory (str): name of column in table for trajectory variable
        threshold_trajectory (int): maximal value of trajectory variable
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Organize information for plot.
    # Copy information.
    table = table.copy(deep=True)
    # Select relevant information.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    columns = [
        #"eid",
        "IID",
        column_trajectory, column_phenotype
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table.loc[
        (table[column_trajectory] < threshold_trajectory), :
    ]
    # Aggregate phenotype values by day.
    table.set_index(
        [column_trajectory],
        append=False,
        drop=True,
        inplace=True
    )
    groups = table.groupby(level=[column_trajectory], axis="index",)
    #groups.aggregate(
    #    mean=(column_phenotype, numpy.nanmean),
    #    deviation=(column_phenotype, numpy.nanstd)
    #)
    records = list()
    for name, group in groups:
        # Initialize missing values.
        days = float("nan")
        count = float("nan")
        mean = float("nan")
        standard_error = float("nan")
        confidence_95 = float("nan")
        confidence_95_low = float("nan")
        confidence_95_high = float("nan")

        table_group = group.reset_index(
            level=None,
            inplace=False,
            drop=False,
        )
        trajectory_units = table_group[column_trajectory].dropna().to_list()[0]
        array = copy.deepcopy(table_group[column_phenotype].dropna().to_numpy())
        # Determine count of valid values.
        count = int(array.size)
        if (count > 5):
            # Determine summary statistics for values in array.
            #mean = round(numpy.nanmean(array), 3)
            mean = numpy.nanmean(array)
            standard_error = scipy.stats.sem(array)
            confidence_95 = (1.96 * standard_error)
            confidence_95_low = (mean - confidence_95)
            confidence_95_high = (mean + confidence_95)
            pass
        # Collect information.
        record = dict()
        record["trajectory_units"] = trajectory_units
        record["count"] = count
        record["mean"] = mean
        record["error"] = standard_error
        record["confidence_95"] = confidence_95
        record["confidence_95_low"] = confidence_95_low
        record["confidence_95_high"] = confidence_95_high
        records.append(record)
        pass
    table_aggregate = pandas.DataFrame(data=records)

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Create figure.
    figure = plot.plot_scatter_points_discrete_abscissa_ordinate_error_bars(
        table=table_aggregate,
        abscissa="trajectory_units",
        ordinate="mean",
        ordinate_error_low="confidence_95",
        ordinate_error_high="confidence_95",
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        fonts=fonts,
        colors=colors,
        size=15,
        label_title=label_title,
    )
    # Return.
    return figure


def organize_phenotypes_parameters_for_plots_histogram():
    """
    Organizes parameters for histogram plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information about each cohort and model.
    records = list()

    #"age",
    #"body",
    #"month",
    #"albumin",
    #"albumin_imputation",
    #"steroid_globulin",
    #"steroid_globulin_imputation", "steroid_globulin_imputation_log",
    #"oestradiol",
    #"oestradiol_imputation", "oestradiol_imputation_log",
    #"oestradiol_bioavailable",
    #"oestradiol_bioavailable_imputation",
    #"oestradiol_bioavailable_imputation_log",
    #"oestradiol_free",
    #"oestradiol_free_imputation", "oestradiol_free_imputation_log",
    #"testosterone", "testosterone_imputation",
    #"testosterone_bioavailable", "testosterone_bioavailable_imputation",
    #"testosterone_free", "testosterone_free_imputation",
    #"vitamin_d", "vitamin_d_imputation",
    #"alcohol_frequency",
    #"alcohol_drinks_weekly",
    #"alcohol_drinks_monthly",
    #"alcohol_drinks_monthly_combination",

    #record = dict()
    #record["name"] = "alcohol_frequency"
    #record["variable"] = "alcohol_frequency"
    #record["threshold"] = 15 # None or upper threshold
    #record["bins"] = None # None or count of bins
    #records.append(record)


    #record = dict()
    #record["name"] = "alcohol_drinks_weekly"
    #record["variable"] = "alcohol_drinks_weekly"
    #record["threshold"] = 100 # None or upper threshold
    #record["bins"] = 70 # None or count of bins
    #records.append(record)

    #record = dict()
    #record["name"] = "alcohol_drinks_monthly"
    #record["variable"] = "alcohol_drinks_monthly"
    #record["threshold"] = 40 # None or upper threshold
    #record["bins"] = 40 # None or count of bins
    #records.append(record)

    #record = dict()
    #record["name"] = "alcohol_drinks_monthly_combination"
    #record["variable"] = "alcohol_drinks_monthly_combination"
    #record["threshold"] = 400 # None or upper threshold
    #record["bins"] = 70 # None or count of bins
    #records.append(record)

    #record = dict()
    #record["name"] = "alcohol_drinks_weekly_log"
    #record["variable"] = "alcohol_drinks_weekly_log"
    #record["threshold"] = None # None or upper threshold
    #record["bins"] = 70 # None or count of bins
    #records.append(record)

    #record = dict()
    #record["name"] = "alcohol_drinks_monthly_log"
    #record["variable"] = "alcohol_drinks_monthly_log"
    #record["threshold"] = None # None or upper threshold
    #record["bins"] = 40 # None or count of bins
    #records.append(record)

    #record = dict()
    #record["name"] = "alcohol_drinks_monthly_combination_log"
    #record["variable"] = "alcohol_drinks_monthly_combination_log"
    #record["threshold"] = None # None or upper threshold
    #record["bins"] = 70 # None or count of bins
    #records.append(record)

    record = dict()
    record["name"] = "albumin_imputation"
    record["variable"] = "albumin_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "albumin_imputation_log"
    record["variable"] = "albumin_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "steroid_globulin_imputation"
    record["variable"] = "steroid_globulin_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "steroid_globulin_imputation_log"
    record["variable"] = "steroid_globulin_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol"
    record["variable"] = "oestradiol"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_log"
    record["variable"] = "oestradiol_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_imputation"
    record["variable"] = "oestradiol_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_imputation_log"
    record["variable"] = "oestradiol_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_bioavailable_imputation"
    record["variable"] = "oestradiol_bioavailable_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_bioavailable_imputation_log"
    record["variable"] = "oestradiol_bioavailable_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_free_imputation"
    record["variable"] = "oestradiol_free_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_free_imputation_log"
    record["variable"] = "oestradiol_free_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone"
    record["variable"] = "testosterone"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_log"
    record["variable"] = "testosterone_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_imputation"
    record["variable"] = "testosterone_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_imputation_log"
    record["variable"] = "testosterone_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_bioavailable_imputation"
    record["variable"] = "testosterone_bioavailable_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_bioavailable_imputation_log"
    record["variable"] = "testosterone_bioavailable_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_free_imputation"
    record["variable"] = "testosterone_free_imputation"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    record = dict()
    record["name"] = "testosterone_free_imputation_log"
    record["variable"] = "testosterone_free_imputation_log"
    record["threshold"] = None # None or upper threshold
    record["bins"] = 70 # None or count of bins
    records.append(record)

    # Return information
    return records


def organize_plots_histogram(
    records_cohorts=None,
    report=None,
):
    """
    Organizes description histogram plots.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Organize phenotypes and parameters for plots.
    records_phenotypes = organize_phenotypes_parameters_for_plots_histogram()

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record_cohort in records_cohorts:
        #if (record_cohort["menstruation"]):
        #    #phenotypes.append("menstruation_days_threshold")
        #    pass
        # Iterate on phenotypes.
        for record_phenotype in records_phenotypes:
            # Organize information for plot.
            name = record_phenotype["name"]
            threshold = record_phenotype["threshold"]
            variable = record_phenotype["variable"]
            bins = record_phenotype["bins"]
            # Copy information.
            table_cohort = record_cohort["table"].copy(deep=True)
            # Apply threshold.
            if (threshold is not None):
                table_cohort = table_cohort.loc[
                    (table_cohort[variable] < threshold), :
                ]
            # Histogram.
            name_label = str(str(record_cohort["name"]) + "_" + str(name))
            name_plot = str(name_label + "_histogram")
            pail[name_plot] = plot_variable_values_histogram(
                label_title=name_label,
                array=table_cohort[variable].dropna().to_numpy(),
                bins=bins, # None or count of bins
            )
            pass
        pass
    # Return information.
    return pail


def organize_phenotypes_plots_dot_trajectory_month(
    table=None,
):
    """
    Organizes information and plots.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Copy information.
    table = table.copy(deep=True)
    # Prepare table to summarize phenotype variables across cohorts and models.
    # These cohorts and models are simple and do not include multiple covariates
    # for genetic analyses.
    # Collect records of information about each cohort and model.
    pails_cohorts = ukb_strat.stratify_phenotype_cohorts_set_sex_age_menopause(
        table=table,
    )

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for pail_cohort in pails_cohorts:
        # Define phenotypes.
        phenotypes = [
            "cholesterol", "cholesterol_imputation",
            "oestradiol", "oestradiol_imputation",
            "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
            "oestradiol_free", "oestradiol_free_imputation",
            "testosterone", "testosterone_imputation",
            "testosterone_bioavailable",
            "testosterone_bioavailable_imputation",
            "testosterone_free", "testosterone_free_imputation",
            "vitamin_d", "vitamin_d_imputation",
            "steroid_globulin", "steroid_globulin_imputation",
            "albumin", "albumin_imputation",
            "alcohol_frequency",
            "alcohol_drinks_monthly",
            "alcohol_audit",
            "alcohol_auditc",
            "alcohol_auditp",
        ]
        # Iterate on phenotypes.
        for phenotype in phenotypes:
            # Plot hormone trajectories across month of assessment.
            # Chart.
            name_label = str(str(pail_cohort["name"]) + "_" + str(phenotype))
            name_plot = str(name_label + "_dot_month")
            pail[name_plot] = plot_variable_means_dot_trajectory(
                label_title=name_label,
                column_phenotype=phenotype,
                column_trajectory="month_order",
                threshold_trajectory=13,
                title_abscissa="month of assessment and blood draw",
                title_ordinate="mean concentration (95% C.I.)",
                table=pail_cohort["table"],
            )
            pass
        pass
    # Return information.
    return pail


def organize_phenotypes_plots_dot_trajectory_menstruation(
    table=None,
):
    """
    Organizes information and plots.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Copy information.
    table = table.copy(deep=True)
    # Prepare table to summarize phenotype variables across cohorts and models.
    # These cohorts and models are simple and do not include multiple covariates
    # for genetic analyses.
    # Collect records of information about each cohort and model.
    pails_cohorts = ukb_strat.stratify_phenotype_cohorts_set_sex_age_menopause(
        table=table,
    )

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for pail_cohort in pails_cohorts:
        # Determine whether the current cohort and model table is relevant to
        # the plot of phenotypes across days of menstrual cycle.
        if (pail_cohort["menstruation"]):
            # Define phenotypes.
            phenotypes = [
                "cholesterol", "cholesterol_imputation",
                "oestradiol", "oestradiol_imputation",
                "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
                "oestradiol_free", "oestradiol_free_imputation",
                "testosterone", "testosterone_imputation",
                "testosterone_bioavailable",
                "testosterone_bioavailable_imputation",
                "testosterone_free", "testosterone_free_imputation",
                "vitamin_d", "vitamin_d_imputation",
                "steroid_globulin", "steroid_globulin_imputation",
                "albumin", "albumin_imputation",
                "alcohol_frequency",
            ]

            # Iterate on phenotypes.
            for phenotype in phenotypes:
                # Plot hormone trajectories across days of the menstrual cycle.
                # Chart.
                name_label = str(
                    str(pail_cohort["name"]) + "_" + str(phenotype)
                )
                name_plot = str(name_label + "_dot_menstruation")
                pail[name_plot] = plot_variable_means_dot_trajectory(
                    label_title=name_label,
                    column_phenotype=phenotype,
                    column_trajectory="menstruation_days",
                    threshold_trajectory=31,
                    title_abscissa="days of menstrual cycle",
                    title_ordinate="mean concentration (95% C.I.)",
                    table=pail_cohort["table"],
                )
                pass
            pass
        pass
    # Return information.
    return pail


##########
# Write
##########


def write_product_table(
    name=None,
    table=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        table (object): table of information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_parent, str(name + ".tsv")
    )
    # Write information to file.
    table.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_tables(
    pail=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        pail (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    for name in pail.keys():
        write_product_table(
            name=name,
            table=pail[name],
            path_parent=path_parent,
        )
    pass


###############################################################################
# Drivers
# The purpose of these driver functions is to make the module's functionality
# accessible in modular parts.


def execute_description_tables(
    set_cohorts=None,
    set_tables=None,
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive summary tables about phenotypes across cohorts
    from the UK Biobank.

    1. define stratified Cohort Tables
    1.1. read Cohort Tables from File in Stratification directory, such as to use same specific Cohort Tables from GWAS
    1.2. define new Cohort Tables for Phenotypes alone (not limited for available genotypes or kinship)
    1.3. define new Cohort Tables for Phenotypes and Genotypes (ancestry, available genotypes, kinship)
    2. define Description records for each Cohort Table
    2.1. "Missingness Table"
    2.2. "Attribution Table" Count, Percentage
    - - For example: how many female persons have value "3" for "alcohol_frequency" variable
    2.3. "Threshold Table"
    - - For example: percentage of Vitamin D measures above and below threshold for "deficiency"
    2.4. "Quantitation Table"
    - - For example: Min, Max, Median, Mean, Standard Deviation, Standard Error, 95% CI for Age

    arguments:
        set_cohorts (str): name of the set of cohorts for which to create
            description tables; "genotype", "phenotype", or "read"
        set_tables (list<str>): names of description tables to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Define new or read old stratified cohort tables.
    if (set_cohorts == "phenotype"):
        # Read source information from file.
        source = read_source(
            path_dock=paths["dock"],
            report=True,
        )
        records_cohorts = (
            ukb_strat.stratify_phenotype_cohorts_set_description_tables(
                table=source["table_phenotypes"],
        ))
    elif (set_cohorts == "read"):
        path_directory = os.path.join(
            paths["dock"], "stratification", "vitamin_d_linear",
        )
        records_cohorts = read_organize_cohorts(
            path_directory=path_directory,
            report=report,
        )
        pass

    # Create description tables across cohorts.

    # Attribution table.
    if ("attribution" in set_tables):
        table_attribution = organize_description_table_attribution(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        table_attribution = pandas.DataFrame()

    # Missingness table.
    if ("missingness" in set_tables):
        table_missingness = organize_description_table_missingness(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        table_missingness = pandas.DataFrame()

    # Threshold table.

    # Quantitation table.
    if ("quantitation" in set_tables):
        table_quantitation = organize_description_table_quantitation(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        table_quantitation = pandas.DataFrame()

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_description_tables()")
        utility.print_terminal_partition(level=3)

    # Collect information.
    pail = dict()
    pail["table_attribution"] = table_attribution
    pail["table_missingness"] = table_missingness
    #pail["table_threshold"] = table_threshold
    pail["table_quantitation"] = table_quantitation
    # Write product information to file.
    write_product_tables(
        pail=pail,
        path_parent=paths["description_tables"],
    )
    pass


def create_description_plots_from_cohorts(
    set_cohorts=None,
    set_plots=None,
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_cohorts (str): name of the set of cohorts for which to create
            description tables; "genotype", "phenotype", or "read"
        set_plots (list<str>): names of description plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Define new or read old stratified cohort tables.
    if (set_cohorts == "phenotype"):
        # Read source information from file.
        source = read_source(
            path_dock=paths["dock"],
            report=True,
        )
        # TODO: TCW, 08 April 2022
        # TODO: switch from alcohol "current" to alcohol "ever"
        # ukb_strat.stratify_phenotype_cohorts_alcohol_current_sex_age_menopause
        records_cohorts = (
            ukb_strat.stratify_phenotype_cohorts_set_sex_age_menopause(
                table=source["table_phenotypes"],
        ))
        pass

    # Create description plots across cohorts.

    # Histograms.
    if ("histogram" in set_plots):
        pail_histogram = organize_plots_histogram(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        pail_histogram = dict()
        pass

    # Box plots.
    # TODO: TCW 3 August 2021
    # TODO: box plots for hormones in different groups (assessment center?)

    # Dot plots of phenotypes ordinal variables.
    #pail["dot_trajectory_month"] = (
    #    organize_phenotypes_plots_dot_trajectory_month(
    #        table=table,
    #))
    #pail["dot_trajectory_menstruation"] = (
    #    organize_phenotypes_plots_dot_trajectory_menstruation(
    #        table=table,
    #))

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_description_tables()")
        utility.print_terminal_partition(level=3)

    # Collect information.
    pail_write = dict()
    pail_write["histogram"] = pail_histogram
    # Write product information to file.
    plot.write_product_plots_child_directories(
        pail_write=pail_write,
        path_parent=paths["description_plots"],
    )
    pass


def define_parameters_regression_summaries():
    """
    Defines parameters for organization of information from regression
    summaries.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about plots for regression
            summaries

    """

    # Collect records of information about regression summaries and their
    # Forest Plots.
    records = list()

    record = dict()
    record["name"] = "alcohol_use_disorder"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.6
    record["abscissa_maximum"] = 0.6
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.15
    record["abscissa_maximum"] = 0.19
    records.append(record)

    record = dict()
    record["name"] = "alcohol_auditc"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.30
    record["abscissa_maximum"] = 0.25
    records.append(record)

    record = dict()
    record["name"] = "alcohol_quantity"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.13
    record["abscissa_maximum"] = 0.13
    records.append(record)

    # Return information
    return records


def define_parameters_correlation_summaries():
    """
    Defines parameters for organization of information from regression
    summaries.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about plots for regression
            summaries

    """

    # Collect records of information about regression summaries and their
    # Forest Plots.
    records = list()

    record = dict()
    record["name"] = "alcohol_use_disorder"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.7
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "alcohol_quantity"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.3
    record["abscissa_maximum"] = 0.3
    records.append(record)

    record = dict()
    record["name"] = "depression"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.5
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "schizophrenia"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.5
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "bipolar_disorder"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.5
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "bipolar_disorder_type_1"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.5
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "bipolar_disorder_type_2"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.5
    record["abscissa_maximum"] = 0.5
    records.append(record)


    # Return information
    return records


def organize_regression_summary_tables_for_forest_plots(
    pail_regression_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Organize the format of regression summary tables for creation of Forest
    Plots.

    arguments:
        pail_regression_tables (dict<object>): collection of Pandas data-frame
            tables for regression summaries
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of Pandas data-frame tables for
            plots

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            pro_reg.organize_regression_summary_table_for_forest_plots(
                type=record_parameter["regression_type"],
                model_contexts=["marginal"],
                model_adjustments=["adjust", "unadjust",],
                column_variable="variable",
                variables=[
                    "oestradiol_imputation",
                    "oestradiol_bioavailable_imputation",
                    "oestradiol_free_imputation",
                    "testosterone_imputation",
                    "testosterone_bioavailable_imputation",
                    "testosterone_free_imputation",
                    "steroid_globulin_imputation",
                    "albumin_imputation",
                ],
                columns_translations={
                    "model_adjustment": "group",
                    "variable": "category",
                    "parameter": "value",
                    "interval_95": "interval_below",
                },
                labels_categories={
                    "oestradiol_imputation": "ESTR-T",
                    "oestradiol_bioavailable_imputation": "ESTR-B",
                    "oestradiol_free_imputation": "ESTR-F",
                    "testosterone_imputation": "TEST-T",
                    "testosterone_bioavailable_imputation": "TEST-B",
                    "testosterone_free_imputation": "TEST-F",
                    "steroid_globulin_imputation": "SHBG",
                    "albumin_imputation": "ALBU",
                },
                sorts_categories={
                    "oestradiol_imputation": 1,
                    "oestradiol_bioavailable_imputation": 2,
                    "oestradiol_free_imputation": 3,
                    "testosterone_imputation": 4,
                    "testosterone_bioavailable_imputation": 5,
                    "testosterone_free_imputation": 6,
                    "steroid_globulin_imputation": 7,
                    "albumin_imputation": 8,
                },
                column_stratification="cohort",
                table=pail_regression_tables[name_table],
                report=report,
        ))
        pass
    # Return information.
    return pail


# TODO: TCW; 18 April 2022
# TODO: 1. read in the combined and extracted correlation table from collection procedure
# TODO: 2. new function(s) to automate

# TODO: Implement new functions to read the genetic correlation tables from "collection" procedure.
# TODO: Aggregate the tables and parse 1. cohort 2. model_adjustment 3. phenotype_secondary.
# TODO: Also need to interpret "phenotype_primary". (simple translation)
# TODO: Introduce new columns to fit the organization function for regression summary tables.
# TODO: Stratify by phenotype_primary.
# TODO: Write out tables with missing values represented as "nan".
# TODO: I probably will still need to edit the table(s) manually.

# TODO: TCW; 19 April 2022
# TODO: variable 'oestradiol_priority' will depend on the cohort
# TODO: use 'oestradiol_detection' for cohorts below...
# TODO: 'female_postmenopause', 'male_age_low', 'male_age_middle', 'male_age_high'
# TODO: use 'schmitz_2021' for cohorts below...
# TODO: 'female', 'male'


def organize_correlation_summary_tables_for_forest_plots(
    pail_correlation_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Organize the format of regression summary tables for creation of Forest
    Plots.

    arguments:
        pail_correlation_tables (dict<object>): collection of Pandas data-frame
            tables for correlation summaries
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of Pandas data-frame tables for
            plots

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            pro_reg.organize_regression_summary_table_for_forest_plots(
                type=record_parameter["regression_type"],
                model_contexts=["joint"],
                model_adjustments=["adjust", "unadjust",],
                column_variable="variable",
                variables=[
                    "oestradiol_priority",
                    "oestradiol_bioavailable_imputation_log",
                    "oestradiol_free_imputation_log",
                    "testosterone_imputation_log",
                    "testosterone_bioavailable_imputation_log",
                    "testosterone_free_imputation_log",
                    "steroid_globulin_imputation_log",
                    "albumin_imputation",
                ],
                columns_translations={
                    "model_adjustment": "group",
                    "variable": "category",
                    "correlation": "value",
                    "interval_95": "interval_below",
                },
                labels_categories={
                    "oestradiol_priority": "ESTR-T",
                    "oestradiol_bioavailable_imputation_log": "ESTR-B",
                    "oestradiol_free_imputation_log": "ESTR-F",
                    "testosterone_imputation_log": "TEST-T",
                    "testosterone_bioavailable_imputation_log": "TEST-B",
                    "testosterone_free_imputation_log": "TEST-F",
                    "steroid_globulin_imputation_log": "SHBG",
                    "albumin_imputation": "ALBU",
                },
                sorts_categories={
                    "oestradiol_priority": 1,
                    "oestradiol_bioavailable_imputation_log": 2,
                    "oestradiol_free_imputation_log": 3,
                    "testosterone_imputation_log": 4,
                    "testosterone_bioavailable_imputation_log": 5,
                    "testosterone_free_imputation_log": 6,
                    "steroid_globulin_imputation_log": 7,
                    "albumin_imputation": 8,
                },
                column_stratification="cohort",
                table=pail_correlation_tables[name_table],
                report=report,
        ))
        pass
    # Return information.
    return pail


def create_regression_summary_forest_plots(
    pail_plot_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        pail_plot_tables (dict<dict<object>>): collection tree of Pandas
            data-frame tables for plots
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of figures

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            plot.drive_iterate_plot_forest_two_groups(
                pail_tables=pail_plot_tables[record_parameter["name"]],
                column_group="group",
                column_ordinate_label="category_label",
                column_ordinate_sort="category_sort",
                column_abscissa_value="value",
                column_abscissa_interval_below="interval_below",
                column_abscissa_interval_above="interval_above",
                group_one="adjust", # markers for group one are above center
                group_two="unadjust", # markers for group two are below center
                abscissa_minimum=record_parameter["abscissa_minimum"],
                abscissa_maximum=record_parameter["abscissa_maximum"],
                ordinate_title="Sex Hormone or Protein",
                abscissa_title="Marginal Model Coefficient (95% C.I.)",
                label_chart_prefix=str(record_parameter["name"]),
                label_size_ordinate_categories="one",
                size_marker=40,
                space_groups=0.075, # vertical space between groups' markers
        ))
        pass
    # Return information.
    return pail


def read_organize_regression_summaries_create_forest_plots(
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_plots (list<str>): names of description plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plots

    """

    # Read regression summary tables from file.
    pail_regression_tables = read_source_regression_summary_tables(
        path_dock=paths["dock"],
        report=True,
    )
    # Define parameters for regression summaries.
    records_parameters = define_parameters_regression_summaries()
    # Organize tables for charts.
    pail_plot_tables = organize_regression_summary_tables_for_forest_plots(
        pail_regression_tables=pail_regression_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Create Forest Plots.
    pail_plots = create_regression_summary_forest_plots(
        pail_plot_tables=pail_plot_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Return information.
    return pail_plots


def read_organize_correlation_summaries_create_forest_plots(
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_plots (list<str>): names of description plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plots

    """

    # Read regression summary tables from file.
    pail_correlation_tables = read_source_correlation_summary_tables(
        path_dock=paths["dock"],
        report=True,
    )
    # Define parameters for regression summaries.
    records_parameters = define_parameters_correlation_summaries()
    # Organize tables for charts.
    pail_plot_tables = organize_correlation_summary_tables_for_forest_plots(
        pail_correlation_tables=pail_correlation_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Create Forest Plots.
    pail_plots = create_regression_summary_forest_plots(
        pail_plot_tables=pail_plot_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Return information.
    return pail_plots


def create_description_plots_from_summaries(
    set_plots=None,
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_plots (list<str>): names of description plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Forest Plots for summary tables from phenotype regressions.
    if ("forest_regressions" in set_plots):
        # Read regression summary tables from file.
        # Organize information from regression summary tables.
        # Create Forest Plots.
        pail_forest_regressions = (
            read_organize_regression_summaries_create_forest_plots(
                paths=paths,
                report=True,
            )
        )
    else:
        pail_forest_regressions = dict()
        pass

    # Forest Plots for summary tables from genetic correlations.
    if ("forest_correlations" in set_plots):
        # Read regression summary tables from file.
        # Organize information from regression summary tables.
        # Create Forest Plots.
        pail_forest_correlations = (
            read_organize_correlation_summaries_create_forest_plots(
                paths=paths,
                report=True,
            )
        )
    else:
        pail_forest_correlations = dict()
        pass

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: create_description_plots_from_summaries()")
        utility.print_terminal_partition(level=3)

    # Collect information.
    pail_write = dict()
    pail_write["forest_regressions"] = pail_forest_regressions
    pail_write["forest_correlations"] = pail_forest_correlations
    # Write product information to file.
    plot.write_product_plots_child_child_directories(
        pail_write=pail_write,
        path_parent=paths["description_plots"],
    )
    pass


def execute_description_plots(
    set_cohorts=None,
    set_cohort_plots=None,
    set_summary_plots=None,
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_cohorts (str): name of the set of cohorts for which to create
            description tables; "genotype", "phenotype", or "read"
        set_cohort_plots (list<str>): names of description plots to create from
            phenotype tables for cohorts
        set_summary_plots (list<str>): names of description plots to create from
            summary tables for previous analyses
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Create new description plots from phenotype tables for cohorts.
    create_description_plots_from_cohorts(
        set_cohorts=set_cohorts,
        set_plots=set_cohort_plots,
        paths=paths,
        report=report,
    )
    # Create new description plots from summary tables for previous analyses.
    create_description_plots_from_summaries(
        set_plots=set_summary_plots,
        paths=paths,
        report=report,
    )

    pass


################################################################################
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

    if False:
        # Create description tables.
        execute_description_tables(
            set_cohorts="phenotype",
            set_tables=["attribution", "missingness", "threshold", "quantitation",],
            paths=paths,
            report=True,
        )

    if True:
        # Create description plots.
        execute_description_plots(
            set_cohorts="phenotype",
            set_cohort_plots=[], # ["histogram",]
            set_summary_plots=["forest_regressions", "forest_correlations",],
            paths=paths,
            report=True,
        )

    pass


# TODO: old? maybe still useful... probably OBSOLETE ... use calls for old plots
def execute_describe_cohorts_models_phenotypes(
    table=None,
    genotype_cohorts=None,
    set=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        genotype_cohorts (bool): whether to generate description for genotype
            cohorts (slow)
        set (str): name of set of cohorts and models to select and organize
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    ##########

    # TODO: store the contingency table reports in a text file?

    # Organize reports on contingency tables by missingness.
    #organize_report_contingency_table_stratification_by_missingness(
    #    table=table,
    #    report=report,
    #)

    # Organize report summary for missingness in biochemical measurements.
    table_missingness = organize_cohorts_phenotypes_hormones_missingness(
        table=table,
        report=report,
    )

    ##########
    #"column_name" and "threshold_low"

    # Organize report summary for threshold in biochemical measurements.
    table_deficiency = organize_cohorts_hormone_deficiency(
        name_variable="vitamin_d",
        column_variable="vitamin_d_imputation",
        threshold=25.0, # nmol / L
        table=table,
        report=report,
    )

    # Prepare table to summarize phenotype variables across cohorts and models.
    # These cohorts and models are simple and do not include multiple covariates
    # for genetic analyses.

    # Define cohorts for description.
    pails_cohorts = ukb_strat.stratify_phenotype_cohorts_regression(
        table=table,
    )

    # Collect summary records and construct table.
    records = list()
    for collection in pails_cohorts:
        records_cohort = organize_cohort_model_variables_summary_long_records(
            name=collection["name"],
            cohort_model=collection["cohort_model"],
            category=collection["category"],
            table=collection["table"],
        )
        records.extend(records_cohort)
    # Organize table.
    table_phenotypes = pandas.DataFrame(data=records)

    ##########

    if (genotype_cohorts):
        # Read source information from file.
        table_kinship_pairs = ukb_strat.read_source_table_kinship_pairs(
            path_dock=path_dock,
            report=report,
        )
        # Prepare table to summarize phenotype variables across cohorts and models
        # for genetic analyses.
        if (set == "sex_hormones"):
            pail_genotypes = (
                ukb_strat.stratify_genotype_cohorts_linear_set_sex_hormones(
                    table=table,
                    table_kinship_pairs=table_kinship_pairs,
                    report=report,
            ))
        elif (set == "bipolar_disorder_body"):
            pail_genotypes = (
                ukb_strat.stratify_cohorts_genotypes_set_bipolar_body(
                    table=table,
                    table_kinship_pairs=table_kinship_pairs,
                    report=report,
            ))
        else:
            print("set of cohorts and models unrecognizable...")
            pail_genotypes = list()
        # Collect summary records and construct table.
        records = list()
        for collection in pail_genotypes:
            record = organize_cohort_model_variables_summary_wide_record(
                name=collection["name"],
                cohort_model=collection["cohort_model"],
                category=collection["category"],
                table=collection["table"],
            )
            records.append(record)
        # Organize table.
        table_genotypes = pandas.DataFrame(data=records)
    else:
        table_genotypes = pandas.DataFrame()

    ##########

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_analyze_sex_cohorts_hormones()")
        utility.print_terminal_partition(level=3)

    # Collect information.
    pail = dict()
    pail["table_cohorts_measurements_missingness"] = table_missingness
    pail["table_cohorts_vitamin_d_deficiency"] = table_deficiency
    pail["table_summary_cohorts_models_phenotypes"] = table_phenotypes
    pail["table_summary_cohorts_models_genotypes"] = table_genotypes
    # Return information.
    return pail


#
