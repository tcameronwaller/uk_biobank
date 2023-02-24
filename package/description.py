drive_stratify_phenotype_cohorts_set_main"""
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
import promiscuity.regression as pro_reg
import promiscuity.plot as plot
import promiscuity.scale as pscale
import promiscuity.description as pdesc
import uk_biobank.stratification as ukb_strat
import uk_biobank.organization as ukb_organization

###############################################################################
# Functionality

#####
# TODO: TCW; 13 December 2022
# It might be reasonable to move some of the functionality herein to the "promiscuity.utility" module
# or to another module such as "promiscuity.description" (descriptive statistics, report summary tables, etc).
# For example, the main record function for the quantitation table would be generally useful.
#####

# TODO: TCW; 15 December 2022
# Restructure to move the "attribution record" and "quantitation record" functions to "promiscuity".
# Also move the driver functions that assemble the tables to "promiscuity".
# These driver functions can accept the list of cohorts and the parameters for reports on variables.



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

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["description"])
    # Initialize directories.
    utility.create_directories(path=paths["description"])
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
                "cohort_name": "string",
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

    Summary tables about Genetic Correlations between pairs of studies.

    Genetic Correlation summary tables from collection procedure need some
    manual editing.
    1. Designate secondary study '34255042_schmitz_2021_female' as cohort
        'female', model context 'joint', and model adjustment 'adjust'.
    2. Designate cohort-specific phenotypes for variable 'oestradiol_priority'.
        Cohorts 'female_menstruation_regular', 'female_premenopause', and
        'female_perimenopause' have sufficient non-missing measurements for
        linear GWAS on oestradiol measurements with imputation.
        All other cohorts ought to use the logistic GWAS on detectability of
        oestradiol.
    3. Fill in missing values for comparisons that did not complete
        successfully.
        For some reason, several comparisons did not complete for primary
        phenotype 'schizophrenia' and secondary phenotype bioavailable
        oestradiol.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

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
        path_dock, "parameters", "uk_biobank",
        "genetic_correlation_summary_tables",
    )
    # Read all files within parent directory and organize tables.
    pail = (
        utility.read_all_pandas_tables_files_within_parent_directory(
            path_directory_parent=path_directory_parent,
            types_pandas_table_read={
                "study_primary": "string",
                "study_secondary": "string",
                "cohort_name": "string",
                "phenotype_primary": "string",
                "dependence_type": "string",
                "model_adjustment": "string",
                "model_context": "string",
                "variable": "string",
                "correlation": "float32",
                "standard_error": "float32",
                "interval_95": "float32",
                "range_95_low": "float32",
                "range_95_high": "float32",
                "probability": "float32",
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
        cohort_name = str(file_name).replace(".pickle", "").replace("table_", "")
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
        record["name"] = cohort_name
        record["cohort_name"] = cohort_name
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

# Variables:
# "oestradiol_imputation_available"
# "testosterone_imputation_available"
# "steroid_globulin_imputation_available"
# "albumin_imputation_available"


def define_variables_table_attribution():
    """
    Defines values of nominal, categorical, or discrete variables for
    description in attribution table.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about values of variables

    """

    # Collect records of information about each categorical varialbe and value.
    records = list()

    # General variables.
    if True:

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

        record = dict()
        record["name"] = "sex_chromosome_aneuploidy_1"
        record["variable"] = "sex_chromosome_aneuploidy" # categorical or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "sex_discrepancy_identity_genetic_1"
        record["variable"] = "sex_discrepancy_identity_genetic" # categorical or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

        # Variable: "race_white"

        record = dict()
        record["name"] = "race_white"
        record["variable"] = "race_white" # cat. or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "race_non_white"
        record["variable"] = "race_white" # cat. or discrete variable
        record["value"] = 0 # categorical or discrete value of variable
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
    if False:

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

    # Variable: "ancestry_white_british"
    if False:

        record = dict()
        record["name"] = "ancestry_white_british_0"
        record["variable"] = "ancestry_white_british" # categorical or discrete variable
        record["value"] = 0 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "ancestry_white_british_1"
        record["variable"] = "ancestry_white_british" # categorical or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

    # Variable: "alteration_sex_hormone"
    if True:

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
    if True:

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
    if False:

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
    if False:

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
    if True:

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
    if True:

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
    if False:

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
    if False:

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
    if False:

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

    # Variable: "menstruation_phase"
    # 0: follicular phase of menstruation
    # 1: luteal phase of menstruation
    if False:

        record = dict()
        record["name"] = "menstruation_follicular"
        record["variable"] = "menstruation_phase" # categorical or discrete variable
        record["value"] = 0 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "menstruation_luteal"
        record["variable"] = "menstruation_phase" # categorical or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

    # Variable: "depression_control_case_loose"
    if False:

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

    # Variable: "depression_control_case_loose"
    if False:

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

    # Variable: "bipolar_control_case_strict"
    if False:

        record = dict()
        record["name"] = "bipolar_control_case_strict_0"
        record["variable"] = "bipolar_control_case_strict"
        record["value"] = 0 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "bipolar_control_case_strict_1"
        record["variable"] = "bipolar_control_case_strict"
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

    # Variable: "bipolar_control_case_loose"
    if False:

        record = dict()
        record["name"] = "bipolar_control_case_loose_0"
        record["variable"] = "bipolar_control_case_loose"
        record["value"] = 0 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "bipolar_control_case_loose_1"
        record["variable"] = "bipolar_control_case_loose"
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

    # Variable: "alcohol_ever"
    if False:

        record = dict()
        record["name"] = "alcohol_never"
        record["variable"] = "alcohol_ever" # categorical or discrete variable
        record["value"] = 0 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "alcohol_ever"
        record["variable"] = "alcohol_ever" # categorical or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

    # Variable: "alcohol_frequency"
    # 0: "Never"
    # 1: "Special occasions only"
    # 2: "One to three times a month"
    # 3: "Once or twice a week"
    # 4: "Three or four times a week"
    # 5: "Daily or almost daily"
    if False:

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

    # Variable: "alcohol_drinks_monthly_ordinal"
    # 0: person had never consumed alcohol in life
    # 1: person consumed alcohol previously in life but not current at recruitment
    # 2: person consumed 30 or fewer drinks of alcohol per month
    # 3: person consumed 31 to 60 drinks of alcohol per month
    # 4: person consumed 61 or more drinks of alcohol per month
    if False:

        record = dict()
        record["name"] = "alcohol_consumption_0"
        record["variable"] = "alcohol_drinks_monthly_ordinal" # cat. or discrete variable
        record["value"] = 0 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "alcohol_consumption_1"
        record["variable"] = "alcohol_drinks_monthly_ordinal" # cat. or discrete variable
        record["value"] = 1 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "alcohol_consumption_2"
        record["variable"] = "alcohol_drinks_monthly_ordinal" # cat. or discrete variable
        record["value"] = 2 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "alcohol_consumption_3"
        record["variable"] = "alcohol_drinks_monthly_ordinal" # cat. or discrete variable
        record["value"] = 3 # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "alcohol_consumption_4"
        record["variable"] = "alcohol_drinks_monthly_ordinal" # cat. or discrete variable
        record["value"] = 4 # categorical or discrete value of variable
        records.append(record)

    # Variable: "alcoholism_control_case_1"
    if False:

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

    # Return information
    return records


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

    # Define values of variables for attribution.
    records_attribution = define_variables_table_attribution()
    # Create a table from records of attribution.
    table_attribution = pdesc.drive_assemble_attribution_table(
        records_attribution=records_attribution,
        records_cohorts=records_cohorts,
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_attribution()")
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table_attribution



##########
# Missingness Table
# Measurement missingness percentages across cohorts
##########


def organize_missingness_record_basis(
    cohort_name=None,
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
        cohort_name (str): name of cohort
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
    record["cohort_name"] = str(cohort_name)
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
    cohort_name=None,
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
        cohort_name (str): name of cohort
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
    record["cohort_name"] = str(cohort_name)
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
    cohort_name=None,
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
        cohort_name (str): name of cohort
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
            cohort_name=cohort_name,
            name_variable=name_variable,
            column_measurement=column_measurement,
            column_detection=column_detection,
            column_missingness_range=column_missingness_range,
            column_reportability_limit=column_reportability_limit,
            table=table,
        )
    else:
        record = organize_missingness_record_basis(
            cohort_name=cohort_name,
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
                cohort_name=record_cohort["name"],
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
    cohort_name=None,
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
        cohort_name (str): name of cohort
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
    record["cohort_name"] = str(cohort_name)
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
            cohort_name=collection_cohort["name"],
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


def define_variables_table_quantitation():
    """
    Defines discrete or continuous variables on ordinal, interval, or ratio
    scales for description in quantitation table.

    arguments:

    raises:

    returns:
        (list<str>): names of variables

    """

    # Define variables.
    variables = [
        "age", "body",
        #"albumin",
        "albumin_imputation",
        #"steroid_globulin",
        "steroid_globulin_imputation",
        #"cholesterol", "cholesterol_imputation",
        #"vitamin_d",
        "vitamin_d_imputation",
        #"oestradiol",
        "oestradiol_imputation",
        #"oestradiol_bioavailable",
        "oestradiol_bioavailable_imputation",
        #"oestradiol_free",
        "oestradiol_free_imputation",
        #"testosterone",
        "testosterone_imputation",
        #"testosterone_bioavailable",
        "testosterone_bioavailable_imputation",
        #"testosterone_free",
        "testosterone_free_imputation",
        #"menstruation_days", "menstruation_duration",
        #"pregnancies", "pregnancies_early", "births",
        #"age_menarche", # necessary data-field is unavailable
        #"age_menopause_self_report",
        #"age_menopause_never_oophorectomy",
        #"age_oophorectomy", "age_hysterectomy",
        #"neuroticism",
        #"alcohol_frequency",
        #"alcohol_drinks_weekly", "alcohol_drinks_monthly",
        "alcohol_drinks_monthly_combination",
        "alcohol_auditc", "alcohol_auditp", "alcohol_audit",
    ]

    # Return information
    return variables


def organize_description_table_quantitation(
    records_cohorts=None,
    report=None,
):
    """
    Drives the assembly of a description table from records of quantitative
    descriptive statistics on variables of interest.

    These descriptive statistics are most appropriate for continuous variables
    on interval, or ratio scales, but they can also be informative for discrete
    variables on ordinal scales.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Define variables.
    variables = define_variables_table_quantitation()
    # Create a table from records of quantitation.
    table_quantitation = pdesc.drive_assemble_quantitation_table(
        variables=variables,
        variable_attribution="genotype_availability",
        value_attribution=1,
        records_cohorts=records_cohorts,
        report=report,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_quantitation()")
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table_quantitation


##########
# Other Tables
# Not currently in use
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



##########
# Management of tables on phenotype variables within stratification cohorts
##########


def prepare_phenotype_variables_in_stratification_cohorts(
    set_tables=None,
    table=None,
    report=None,
):
    """
    Prepare phenotype variables within stratification cohorts for relevant types
    of plots.

    arguments:
        set_tables (list<str>): names of description tables to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<list<dict>>): collection of records with information about cohorts

    """

    # Collect information.
    pail_cohorts = dict()

    # Prepare phenotype variables in stratification cohorts specific to each
    # type of plot.
    if ("attribution" in set_tables):
        # Initialize collection.
        records_cohorts = list()
        # Attribution Table.
        # Stratify records within separate tables for cohorts.
        records_cohorts = (
            ukb_strat.drive_stratify_phenotype_cohorts_set_main(
                table=table,
        ))
        # Filter relevant cohorts.
        if False:
            names_cohorts = [
                #"female_male",
                #"female",
                #"female_menstruation_regular",
                #"female_premenopause",
                #"female_perimenopause",
                #"female_postmenopause",
                #"male",
                #"male_age_low",
                #"male_age_middle",
                #"male_age_high",
                #"race_white_female_male",
                #"race_white_female",
                #"race_white_female_menstruation_regular",
                #"race_white_female_premenopause",
                #"race_white_female_perimenopause",
                #"race_white_female_postmenopause",
                #"race_white_male",
                #"race_white_male_age_low",
                #"race_white_male_age_middle",
                #"race_white_male_age_high",
                "bipolar_case_female_male",
                "bipolar_case_female",
                "bipolar_case_male",
                "bipolar_control_female_male",
                "bipolar_control_female",
                "bipolar_control_male",
                #"race_white_alcohol_current_female_male",
                #"race_white_alcohol_current_female",
                #"race_white_alcohol_current_female_menstruation_regular",
                #"race_white_alcohol_current_female_premenopause",
                #"race_white_alcohol_current_female_perimenopause",
                #"race_white_alcohol_current_female_postmenopause",
                #"race_white_alcohol_current_male",
                #"race_white_alcohol_current_male_age_low",
                #"race_white_alcohol_current_male_age_middle",
                #"race_white_alcohol_current_male_age_high",
            ]
            records_cohorts = utility.filter_records_by_name(
                names=names_cohorts,
                records=records_cohorts,
                report=True,
            )
        # Collect information.
        pail_cohorts["attribution"] = copy.deepcopy(records_cohorts)
    if ("missingness" in set_tables):
        # Initialize collection.
        records_cohorts = list()
        # Box plots for groups.
        records_cohorts = (
            ukb_strat.drive_stratify_phenotype_cohorts_set_main(
                table=table,
        ))
        # Collect information.
        pail_cohorts["missingness"] = copy.deepcopy(records_cohorts)
    if ("quantitation" in set_tables):
        # Initialize collection.
        records_cohorts = list()
        # Attribution Table.
        # Stratify records within separate tables for cohorts.
        records_cohorts = (
            ukb_strat.drive_stratify_phenotype_cohorts_set_main(
                table=table,
        ))
        # Filter relevant cohorts.
        if False:
            names_cohorts = [
                #"female_male",
                #"female",
                #"female_menstruation_regular",
                #"female_premenopause",
                #"female_perimenopause",
                #"female_postmenopause",
                #"male",
                #"male_age_low",
                #"male_age_middle",
                #"male_age_high",
                #"race_white_female_male",
                #"race_white_female",
                #"race_white_female_menstruation_regular",
                #"race_white_female_premenopause",
                #"race_white_female_perimenopause",
                #"race_white_female_postmenopause",
                #"race_white_male",
                #"race_white_male_age_low",
                #"race_white_male_age_middle",
                #"race_white_male_age_high",
                "bipolar_case_female_male",
                "bipolar_case_female",
                "bipolar_case_male",
                "bipolar_control_female_male",
                "bipolar_control_female",
                "bipolar_control_male",
                #"race_white_alcohol_current_female_male",
                #"race_white_alcohol_current_female",
                #"race_white_alcohol_current_female_menstruation_regular",
                #"race_white_alcohol_current_female_premenopause",
                #"race_white_alcohol_current_female_perimenopause",
                #"race_white_alcohol_current_female_postmenopause",
                #"race_white_alcohol_current_male",
                #"race_white_alcohol_current_male_age_low",
                #"race_white_alcohol_current_male_age_middle",
                #"race_white_alcohol_current_male_age_high",
            ]
            records_cohorts = utility.filter_records_by_name(
                names=names_cohorts,
                records=records_cohorts,
                report=True,
            )
        # Apply Distribution Scale Transformations to variables of interest in
        # each cohort.
        if False:
            records_cohorts = (
                pscale.drive_transformations_on_multiple_variables_in_cohorts(
                    variables=[
                        #"body",
                        #"alcohol_drinks_monthly_combination",
                        #"alcohol_frequency",
                        #"alcohol_auditc",
                        #"vitamin_d_imputation",
                        #"oestradiol_imputation",
                        #"oestradiol_bioavailable_imputation",
                        #"oestradiol_free_imputation",
                        #"testosterone_imputation",
                        #"testosterone_bioavailable_imputation",
                        #"testosterone_free_imputation",
                        #"steroid_globulin_imputation",
                        #"albumin_imputation",
                    ],
                    records_cohorts=records_cohorts,
                    report=True,
            ))
        # Collect information.
        pail_cohorts["quantitation"] = copy.deepcopy(records_cohorts)
        pass
    # Return information
    return pail_cohorts


def create_tables_for_phenotype_variables_in_cohorts(
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
        # Prepare information about phenotype variables within stratification
        # cohorts.
        pail_cohorts = prepare_phenotype_variables_in_stratification_cohorts(
            set_tables=set_tables,
            table=source["table_phenotypes"],
            report=report,
        )
        pass
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
            records_cohorts=pail_cohorts["attribution"],
            report=report,
        )
    else:
        table_attribution = pandas.DataFrame()

    # Missingness table.
    if ("missingness" in set_tables):
        table_missingness = organize_description_table_missingness(
            records_cohorts=pail_cohorts["missingness"],
            report=report,
        )
    else:
        table_missingness = pandas.DataFrame()

    # Threshold table.

    # Quantitation table.
    if ("quantitation" in set_tables):
        table_quantitation = organize_description_table_quantitation(
            records_cohorts=pail_cohorts["quantitation"],
            report=report,
        )
    else:
        table_quantitation = pandas.DataFrame()

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: create_tables_for_phenotype_variables_in_cohorts()")
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
        path_parent=paths["description"],
    )
    pass


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
    print("version check: TCW, 12 December 2022")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )

    if True:
        # Create description tables.
        create_tables_for_phenotype_variables_in_cohorts(
            set_cohorts="phenotype",
            set_tables=[
                "attribution",
            ], # "attribution", "missingness", "threshold", "quantitation",
            paths=paths,
            report=True,
        )

    pass


# TODO: old? maybe still useful... probably OBSOLETE ... check calls to description tables
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
