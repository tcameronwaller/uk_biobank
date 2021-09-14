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

    # Specify directories and files.
    path_table_linear_hormones_sex_age_menopause = os.path.join(
        path_dock, "parameters", "uk_biobank", "regression_cohorts_models",
        "table_linear_hormones_sex_age_menopause.tsv"
    )
    path_table_linear_hormones_alcoholism = os.path.join(
        path_dock, "parameters", "uk_biobank", "regression_cohorts_models",
        "table_linear_hormones_alcoholism.tsv"
    )

    # Organize code interpretations.
    table_linear_hormones_sex_age_menopause = pandas.read_csv(
        path_table_linear_hormones_sex_age_menopause,
        sep="\t",
        header=0,
        dtype={
            "cohort": "string",
            "model": "string",
            "dependence": "string",
            "independence": "string",
        },
    )
    table_linear_hormones_sex_age_menopause.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )

    table_linear_hormones_alcoholism = pandas.read_csv(
        path_table_linear_hormones_alcoholism,
        sep="\t",
        header=0,
        dtype={
            "cohort": "string",
            "model": "string",
            "dependence": "string",
            "independence": "string",
        },
    )
    table_linear_hormones_alcoholism.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )


    # Compile and return information.
    return {
        "table_linear_hormones_sex_age_menopause": (
            table_linear_hormones_sex_age_menopause
        ),
        "table_linear_hormones_alcoholism": table_linear_hormones_alcoholism,
    }


##########
# Need a coherent name for this section...


def define_model_dependence_records_hormones():
    """
    Define records for regression model and dependent variables.

    arguments:

    raises:

    returns:
        (list<dict>): information from regressions

    """

    records = [
        {"dependence": "vitamin_d_log", "model": "complex"},
        #{"dependence": "vitamin_d_log", "model": "simple"},
        {"dependence": "vitamin_d_imputation_log", "model": "complex"},
        {"dependence": "albumin", "model": "complex"},
        #{"dependence": "albumin", "model": "alternate_one"},
        #{"dependence": "albumin", "model": "simple"},
        {"dependence": "albumin_imputation", "model": "complex"},
        {"dependence": "steroid_globulin_log", "model": "complex"},
        #{"dependence": "steroid_globulin_log", "model": "simple"},
        {"dependence": "steroid_globulin_imputation_log", "model": "complex"},
        {"dependence": "oestradiol_log", "model": "complex"},
        #{"dependence": "oestradiol_log", "model": "alternate_one"},
        #{"dependence": "oestradiol_log", "model": "alternate_two"},
        #{"dependence": "oestradiol_log", "model": "simple"},
        {"dependence": "oestradiol_imputation", "model": "complex"},
        {"dependence": "oestradiol_bioavailable_log", "model": "complex"},
        {
            "dependence": "oestradiol_bioavailable_imputation",
            "model": "complex"
        },
        {"dependence": "oestradiol_free_log", "model": "complex"},
        {"dependence": "oestradiol_free_imputation", "model": "complex"},
        {"dependence": "testosterone_log", "model": "complex"},
        #{"dependence": "testosterone_log", "model": "simple"},
        {"dependence": "testosterone_imputation", "model": "complex"},
        {"dependence": "testosterone_bioavailable_log", "model": "complex"},
        {
            "dependence": "testosterone_bioavailable_imputation",
            "model": "complex"
        },
        {"dependence": "testosterone_free_log", "model": "complex"},
        {"dependence": "testosterone_free_imputation", "model": "complex"},
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
    cohorts_records = ukb_strat.stratify_set_primary_sex_age_body_menopause(
        table=table
    )
    cohorts_relevant = [
        "female", "female_premenopause", "female_perimenopause",
        "female_postmenopause",
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

    # Drive regressions.
    if False:
        pail_female = drive_regressions_female_cohorts_models_hormones(
            table=source["table_phenotypes"],
            table_cohort_model=(
                source_reference["table_linear_hormones_sex_age_menopause"]
            ),
            report=True
        )
    if True:
        pail_alcohol = drive_linear_regressions_hormones_alcoholism(
            table=source["table_phenotypes"],
            table_cohort_model=(
                source_reference["table_linear_hormones_alcoholism"]
            ),
            report=True
        )



    pass





#
