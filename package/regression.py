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
import promiscuity.utility as utility
import promiscuity.regression as pro_regression
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

    # Collect tables.
    pail = dict()
    # Iterate on tables.
    file_names = [
        "table_linear_hormones_sex_age_menopause.tsv",
        "table_linear_vitamin-d_sex_age_body.tsv",
        "table_linear_hormones_alcoholism_any.tsv",
        "table_linear_hormones_alcoholism_1.tsv",
        "table_linear_hormones_bipolar_loose.tsv",
        "table_linear_hormones_bipolar_strict.tsv",
    ]
    for file_name in file_names:
        # Determine table name.
        name = file_name.replace(".tsv", "")
        # Specify directories and files.
        path_table = os.path.join(
            path_dock, "parameters", "uk_biobank", "regression_cohorts_models",
            file_name,
        )
        # Read information from file.
        table = pandas.read_csv(
            path_table,
            sep="\t",
            header=0,
            dtype={
                "cohort": "string",
                "model": "string",
                "dependence": "string",
                "independence": "string",
            },
        )
        table.reset_index(
            level=None,
            inplace=True,
            drop=True,
        )
        # Collect table.
        pail[name] = table.copy(deep=True)
        pass
    # Return information.
    return pail


##########
# Need a coherent name for this section...


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

# TODO: TCW 9 November 2021
# TODO: this function is ALMOST identical to "drive_linear_regressions_hormones_alcoholism"
# TODO: consolidate these functions by passing an argument for the main independent
# TODO: variable of interest

def drive_linear_regressions_cohorts_models_dependences(
    table=None,
    table_cohort_model=None,
    independences_summary=None,
    report=None,
):
    """
    Organize regressions.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_cohort_model (object): Pandas data frame of cohorts, models,
            dependent variables, and independent variables for regression
        independences_summary (list<str>): names of independent variables for
            which to include statistics in the summary table
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Define relevant cohorts.
    cohorts_records = ukb_strat.stratify_set_primary_sex_age_body_menopause(
        table=table,
    )
    cohorts_relevant = [
        "female_male",
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
        independence=independences_summary,
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
# Write


def write_product_table(
    name=None,
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_parent, str(name + ".tsv")
    )
    # Write information to file.
    information.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_tables(
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

    for name in information.keys():
        write_product_table(
            name=name,
            information=information[name],
            path_parent=path_parent,
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

    # Export information.
    write_product_tables(
        information=information["tables"],
        path_parent=paths["regression"],
    )
    pass


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
        path_dock=path_dock,
        report=True,
    )
    source_reference = read_source_cohort_model_reference(
        path_dock=path_dock,
    )

    # table_linear_vitamin-d_sex_age_body

    # Drive regressions.
    if True:
        pass
    if False:
        pail_bip_loose = drive_linear_regressions_cohorts_models_dependences(
            table=source["table_phenotypes"],
            table_cohort_model=(
                source_reference["table_linear_hormones_bipolar_loose"]
            ),
            independences_summary=["bipolar_control_case_loose"],
            report=True
        )
        pail_bip_strict = drive_linear_regressions_cohorts_models_dependences(
            table=source["table_phenotypes"],
            table_cohort_model=(
                source_reference["table_linear_hormones_bipolar_strict"]
            ),
            independences_summary=["bipolar_control_case_strict"],
            report=True
        )
    if False:
        pail_regression = drive_regressions_female_cohorts_models_hormones(
            table=source["table_phenotypes"],
            table_cohort_model=(
                source_reference["table_linear_hormones_sex_age_menopause"]
            ),
            report=True
        )
    if False:
        pail_regression = drive_linear_regressions_hormones_alcoholism(
            table=source["table_phenotypes"],
            table_cohort_model=(
                source_reference["table_linear_hormones_alcoholism_1"]
            ),
            independences_summary=["alcoholism_control_case_1"],
            report=True
        )

    # Collect information.
    information = dict()
    information["tables"] = dict()
    information["tables"]["table_regressions_bipolar_loose"] = (
        pail_bip_loose["table"]
    )
    information["tables"]["table_regressions_bipolar_strict"] = (
        pail_bip_strict["table"]
    )
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )
    pass





#
