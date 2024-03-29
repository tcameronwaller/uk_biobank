"""
Organize regression analyses on data from the U.K. Biobank.

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
import partner.utility as utility
import partner.regression as pro_reg
import partner.decomposition as decomp
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
    paths["regression"] = os.path.join(path_dock, "regression")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["regression"])
    # Initialize directories.
    utility.create_directories(
        path=paths["regression"]
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


def read_source_cohort_model_reference(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        path_dock, "parameters", "uk_biobank", "regression_cohorts_models",
    )
    # Read all files within parent directory and organize tables.
    pail = utility.read_all_pandas_tables_files_within_parent_directory(
        path_directory_parent=path_directory_parent,
        types_pandas_table_read={
            "execution": "int",
            "cohort": "string",
            "cohort_sort": "int",
            "dependence": "string",
            "dependence_sort": "int",
            "dependence_type": "string",
            "independence": "string",
            "model": "string",
            "model_sort": "int",
            "model_adjustment": "string",
            "model_context": "string",
            "model_note": "string",
        },
        report=report,
    )
    # Return information.
    return pail


##########
# Need a coherent name for this section...

# TODO: TCW 5 January 2022
# TODO: I introduced column "execution" to the cohort-dependence-model table to turn regressions "on" or "off"
# TODO: need to implement code to read that "execution" and handle accordingly... probably just filter to "1"

# TODO: TCW 9 November 2021
# TODO: this structure seems quite inefficient to me...
# TODO: It might be more efficient to define these "dependent" variables and
# TODO: regression models of interest in a separate table?
def define_model_dependence_records_hormones():
    """
    Define records for regression model and dependent variables.

    arguments:

    raises:

    returns:
        (list<dict>): information from regressions

    """

    records = [

        {"dependence": "albumin", "model": "complex"},
        {"dependence": "albumin", "model": "unadjust"},
        {"dependence": "albumin_imputation", "model": "complex"},
        {"dependence": "albumin_imputation", "model": "unadjust"},

        {"dependence": "steroid_globulin_log", "model": "complex"},
        {"dependence": "steroid_globulin_log", "model": "unadjust"},
        {"dependence": "steroid_globulin_imputation_log", "model": "complex"},
        {"dependence": "steroid_globulin_imputation_log", "model": "unadjust"},

        {"dependence": "oestradiol_log", "model": "complex"},
        {"dependence": "oestradiol_log", "model": "unadjust"},
        {"dependence": "oestradiol_imputation", "model": "complex"},
        {"dependence": "oestradiol_imputation", "model": "unadjust"},

        {"dependence": "oestradiol_bioavailable_log", "model": "complex"},
        {"dependence": "oestradiol_bioavailable_log", "model": "unadjust"},
        {
            "dependence": "oestradiol_bioavailable_imputation",
            "model": "complex"
        },
        {
            "dependence": "oestradiol_bioavailable_imputation",
            "model": "unadjust"
        },

        {"dependence": "oestradiol_free_log", "model": "complex"},
        {"dependence": "oestradiol_free_log", "model": "unadjust"},
        {"dependence": "oestradiol_free_imputation", "model": "complex"},
        {"dependence": "oestradiol_free_imputation", "model": "unadjust"},

        {"dependence": "testosterone_log", "model": "complex"},
        {"dependence": "testosterone_log", "model": "unadjust"},
        {"dependence": "testosterone_imputation", "model": "complex"},
        {"dependence": "testosterone_imputation", "model": "unadjust"},

        {"dependence": "testosterone_bioavailable_log", "model": "complex"},
        {"dependence": "testosterone_bioavailable_log", "model": "unadjust"},
        {
            "dependence": "testosterone_bioavailable_imputation",
            "model": "complex"
        },
        {
            "dependence": "testosterone_bioavailable_imputation",
            "model": "unadjust"
        },

        {"dependence": "testosterone_free_log", "model": "complex"},
        {"dependence": "testosterone_free_log", "model": "unadjust"},
        {"dependence": "testosterone_free_imputation", "model": "complex"},
        {"dependence": "testosterone_free_imputation", "model": "unadjust"},

        {"dependence": "vitamin_d_log", "model": "complex"},
        {"dependence": "vitamin_d_log", "model": "unadjust"},
        {"dependence": "vitamin_d_imputation_log", "model": "complex"},
        {"dependence": "vitamin_d_imputation_log", "model": "unadjust"},

    ]
    return records


def drive_regressions_female_cohorts_models_hormones(
    table=None,
    table_cohort_model=None,
    report=None,
):
    """
    Organize regressions.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_cohort_model (object): Pandas data frame of cohorts, models,
            dependent variables, and independent variables for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Define relevant cohorts.
    cohorts_records = ukb_strat.stratify_set_primary_sex_age_body_menopause(
        table=table
    )
    cohorts_relevant = [
        #"female", "female_premenopause", "female_perimenopause",
        #"female_postmenopause",
        "male", "male_age_low", "male_age_middle", "male_age_high",
    ]
    cohorts_records = list(filter(
        lambda cohort_record: (cohort_record["cohort"] in cohorts_relevant),
        cohorts_records
    ))

    # Define relevant dependent variables.
    records_models = define_model_dependence_records_hormones()
    # Iterate across cohorts.
    for cohort_record in cohorts_records:
        cohort = cohort_record["cohort"]
        menstruation = cohort_record["menstruation"]
        table_cohort = cohort_record["table"]
        # Iterate across outcomes (dependent variables).
        for record_model in records_models:
            # Define cohort-specific ordinal representation.
            dependence_ordinal = str(
                str(record_model["dependence"]) + "_" + str(cohort) + "_order"
            )
            pail_regression = (
                pro_regression.drive_cohort_model_linear_regression(
                    table=table_cohort,
                    table_cohort_model=table_cohort_model,
                    cohort=cohort,
                    model=record_model["model"],
                    dependence=record_model["dependence"],
                    report=report,

            ))
            pass
        pass

    # Compile information.
    pail = dict()
    # Return information.
    return pail


def drive_linear_regressions_hormones_alcoholism(
    table=None,
    table_cohort_model=None,
    report=None,
):
    """
    Organize regressions.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_cohort_model (object): Pandas data frame of cohorts, models,
            dependent variables, and independent variables for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Define relevant cohorts.
    #cohorts_records = ukb_strat.stratify_set_primary_sex_age_body_menopause(
    #    table=table
    #)
    cohorts_records = ukb_strat.stratify_set_primary_sex_age_body_menopause(
        table=table,
    )

    cohorts_relevant = [
        "female",
        "female_premenopause",
        "female_perimenopause",
        "female_postmenopause",
        "male",
        "male_age_low",
        "male_age_middle",
        "male_age_high",
    ]
    cohorts_records = list(filter(
        lambda cohort_record: (cohort_record["cohort"] in cohorts_relevant),
        cohorts_records
    ))

    # Collect summary records for each regression.
    records = list()
    # Define relevant dependent variables.
    records_models = define_model_dependence_records_hormones()
    # Iterate across cohorts.
    for cohort_record in cohorts_records:
        cohort = cohort_record["cohort"]
        menstruation = cohort_record["menstruation"]
        table_cohort = cohort_record["table"]
        if report:
            utility.print_terminal_partition(level=2)
            print(cohort)
            print(table_cohort)
        # Iterate across outcomes (dependent variables).
        for record_model in records_models:
            dependence = record_model["dependence"]
            model = record_model["model"]
            # Organize record.
            record = dict()
            record["cohort"] = cohort
            record["dependence"] = dependence
            record["model"] = model
            record["name"] = str(
                record["cohort"] + "_" +
                record["dependence"] + "_" +
                record["model"]
            )
            # Define cohort-specific ordinal representation.
            #dependence_ordinal = str(
            #    str(record_model["dependence"]) + "_" + str(cohort) + "_order"
            #)
            pail_regression = (
                pro_regression.drive_cohort_model_linear_regression(
                    table=table_cohort,
                    table_cohort_model=table_cohort_model,
                    cohort=cohort,
                    model=model,
                    dependence=dependence,
                    report=report,

            ))
            record.update(pail_regression["summary"])
            # Collect records.
            records.append(record)
            pass
        pass

    # Organize table.
    table_regressions_raw = pandas.DataFrame(data=records)
    table_regressions = pro_regression.organize_table_regression_summaries(
        independence=["alcoholism_control_case_1"],
        table=table_regressions_raw,
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_regressions)
    # Compile information.
    pail = dict()
    pail["table"] = table_regressions
    # Return information.
    return pail


def organize_regressions_site_month(
    table=None,
    report=None,
):
    """
    Organize regressions.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Define relevant cohorts.
    cohorts_records = ukb_strat.stratify_set_primary_sex_menopause_age(
        table=table
    )
    cohorts_relevant = [
        "female-premenopause", "female-perimenopause", "female-postmenopause",
        "male", "male-younger", "male-older"
    ]
    cohorts_records = list(filter(
        lambda cohort_record: (cohort_record["cohort"] in cohorts_relevant),
        cohorts_records
    ))

    # TODO: TCW 19 August 2021
    # TODO: eventually... I will need cohort-specific models...

    # Define outcome dependent variables.
    hormones = [
        "vitamin_d", "albumin", "steroid_globulin",
        "oestradiol", "testosterone",
    ]
    # Define predictor independent variables.
    predictors_region = ["assessment_region",]
    predictors_season = ["assessment_season",]
    predictors_site_components = [
        "site-component_1", "site-component_2", "site-component_3",
        "site-component_4", "site-component_5", "site-component_6",
        "site-component_7", "site-component_8", "site-component_9",
        "site-component_10", "site-component_11", "site-component_12",
        "site-component_13", "site-component_14", "site-component_15",
        "site-component_16", "site-component_17", "site-component_18",
        "site-component_19", "site-component_20",
    ]
    predictors_month_components = [
        "month-component_1", "month-component_2", "month-component_3",
        "month-component_4", "month-component_5", "month-component_6",
        "month-component_7", "month-component_8", "month-component_9",
        "month-component_10",
    ]
    predictors_month_indicators = [
        "month-indicator_January", "month-indicator_February",
        "month-indicator_March", #"month-indicator_April",
        "month-indicator_May", "month-indicator_June",
        "month-indicator_July", "month-indicator_August",
        "month-indicator_September", "month-indicator_October",
        "month-indicator_November", "month-indicator_December",
    ]

    # Iterate across cohorts.
    for cohort_record in cohorts_records:
        cohort = cohort_record["cohort"]
        menstruation = cohort_record["menstruation"]
        table_cohort = cohort_record["table"]
        # Iterate across outcomes (dependent variables).
        for hormone in hormones:
            # Define cohort-specific ordinal representation.
            hormone_ordinal = str(str(hormone) + "_" + str(cohort) + "_order")

            # Specify outcome and predictors.

            #outcome = hormone
            outcome = hormone_ordinal

            #predictors = predictors_month_components
            #predictors = predictors_site_components
            predictors = predictors_region
            #predictors = predictors_season


            # Report.
            if report:
                utility.print_terminal_partition(level=3)
                print("report: ")
                print("organize_cohorts_models_phenotypes_regressions()")
                utility.print_terminal_partition(level=5)
                print("cohort: " + str(cohort))
                print("outcome: " + str(outcome))
                utility.print_terminal_partition(level=5)
            pail_regression = regression.regress_linear_ordinary_least_squares(
                dependence=outcome, # parameter
                independence=predictors, # parameter
                threshold_samples=100,
                table=table_cohort,
                report=report,
            )
            pass
        pass

    # Compile information.
    pail = dict()
    # Return information.
    return pail



##########
# Driver


def stratify_cohorts_call_run_regressions(
    table=None,
    table_cohorts_models=None,
    independences_summary=None,
    filter_execution=None,
    type=None,
    report=None,
):
    """
    Stratify cohorts for phenotype tables and call driver for regressions.

    Format of "table"...
    Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows, with an explicit index

    Format of "table_cohorts_models"...
    columns: "execution", "cohort", "cohort_sort", "dependence",
    "dependence_sort", "dependence_type", "model", "model_sort", "independence",
    "model_note"

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_cohorts_models (object): Pandas data frame that specifies cohorts,
            dependent variables, and independent variables (models) for
            regression
        independences_summary (list<str>): names of independent variables for
            which to include information in the summary table, or "None" to
            include information for all original independent variables
        filter_execution (bool): whether to filter records in cohort-model table
            by logical binary "execution" variable
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Stratify phenotypes in cohorts.
    records_cohorts = (
        ukb_strat.drive_stratify_phenotype_cohorts_set_description_tables(
            table=table,
    ))
    entries_cohorts = (
        utility.structure_from_records_to_entries(
            records=records_cohorts,
    ))
    # Call driver for regressions.
    pail = (
        pro_reg.drive_linear_logistic_regressions_cohorts_models(
            entries_cohorts=entries_cohorts,
            table_cohorts_models=table_cohorts_models,
            independences_summary=independences_summary,
            filter_execution=filter_execution,
            type=type,
            report=report,
    ))
    # Return information.
    return pail


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
    print("version check: TCW, 24 October 2022")
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
    source_reference = read_source_cohort_model_reference(
        path_dock=path_dock,
        report=True,
    )

    # Collect information.
    pail_write = dict()
    pail_write["tables"] = dict()

    # Drive regressions.

    if True:
        pail_regression = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_bipolar_disorder_sex_hormones"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic", # "linear" or "logistic"
            report=True,
        )
        pail_write["tables"]["table_bipolar_disorder_sex_hormones"] = (
            pail_regression["table"]
        )
        pass

    # Write product information to file.
    pro_reg.write_product(
        pail_write=pail_write,
        path_parent=paths["regression"],
    )
    pass





#
