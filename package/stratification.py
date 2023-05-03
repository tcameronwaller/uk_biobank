"""
Stratify cohorts for analyses on phenotypes and genotypes from the U.K. Biobank.

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
import promiscuity.scale as pscale

###############################################################################
# Functionality

# "female"
# variables:
# "eid", "IID", "assessment_region", "assessment_season", "white_british", "sex_y", "sex_x", "sex_text", "age", "body", "body_log", "pregnancy", "menopause_ordinal", "hormone_alteration",

# "female_premenopause"
# variables:
# "eid", "IID", "assessment_region", "assessment_season", "white_british", "sex_y", "sex_x", "sex_text", "age", "body", "body_log", "pregnancy", "menopause_binary", "menopause_ordinal", "menstruation_days", "menstruation_phase", "menstruation_phase_cycle", "hormone_alteration",

# "female_perimenopause"
# variables:
# "eid", "IID", "assessment_region", "assessment_season", "white_british", "sex_y", "sex_x", "sex_text", "age", "body", "body_log", "pregnancy", "menopause_binary", "menopause_ordinal", "menstruation_days", "menstruation_phase", "menstruation_phase_cycle", "hormone_alteration",

# "female_postmenopause"
# variables:
# "eid", "IID", "assessment_region", "assessment_season", "white_british", "sex_y", "sex_x", "sex_text", "age", "body", "body_log", "pregnancy", "menopause_binary", "menopause_ordinal", "hormone_alteration",

# "male"
# variables:
# "eid", "IID", "assessment_region", "assessment_season", "white_british", "sex_y", "sex_x", "sex_text", "age", "age_grade_male", "body", "body_log",

#"vitamin_d", "vitamin_d_log",
#"vitamin_d_imputation", "vitamin_d_imputation_log",
#"albumin", "albumin_log",
#"albumin_imputation", "albumin_imputation_log",
#"steroid_globulin", #"steroid_globulin_log",
#"steroid_globulin_imputation", #"steroid_globulin_imputation_log",
#"oestradiol", "oestradiol_log",
#"oestradiol_imputation", "oestradiol_imputation_log",
#"oestradiol_bioavailable", "oestradiol_bioavailable_log",
#"oestradiol_bioavailable_imputation",
#"oestradiol_bioavailable_imputation_log",
#"oestradiol_free", "oestradiol_free_log",
#"oestradiol_free_imputation", "oestradiol_free_imputation_log",
#"testosterone", #"testosterone_log",
#"testosterone_imputation", #"testosterone_imputation_log",
#"testosterone_bioavailable", #"testosterone_bioavailable_log",
#"testosterone_bioavailable_imputation",
#"testosterone_bioavailable_imputation_log",
#"testosterone_free", #"testosterone_free_log",
#"testosterone_free_imputation", #"testosterone_free_imputation_log",

#"albumin_detection",
#"steroid_globulin_detection",
#"oestradiol_detection",
#"testosterone_detection",

######
# Cleanup; TCW; 9 November 2022

# Moved to "genotype" module and modified.
# organize_phenotype_covariate_table_plink_format
# stratify_by_nonmissing_values_specific_variables
# filter_table_phenotype_persons_by_kinship

# Deleted:
# execute_organize_stratification_genotype_cohorts_models
# stratify_genotype_cohorts_models_set_reference_population
# stratify_genotype_cohorts_models_set
# stratify_genotype_cohort_model_instance
# select_records_by_ancestry_case_control_valid_variables_values
# select_records_by_ancestry_sex_specific_valid_variables_values
# select_records_by_female_specific_valid_variables_values
# select_records_by_male_specific_valid_variables_values

###########################################
############################################
############################################






###########################################
############################################
############################################

# TODO: TCW; 23 February 2023

# TODO: MAJOR RESTRUCTURE!!!
# TODO: Each and every cohort record needs the information below:


# I need some sort of "main" function to drive cohort stratifications
# This function needs to accept a parameter for whether or not run exclusions

# Where reasonable, I need some sort of stratification hierarchy...
# "cohort_genotypes" can definitely be a parameter passed to a high-level driver function
# "cohort_genotypes", "cohort_race", "cohort_exclusions" can all be applied either before or after others
#  - - These will be simple, one-step "filter" functions applied on each cohort record passed to them
#  - - use the ".update()" dictionary method to change the corresponding variable


##########
# Phenotype cohorts
# Cohort, model selection: sets for descriptions of phenotypes within cohorts
# Review: TCW; 3 May 2023

# review: TCW; 23 February 2023
def stratify_phenotype_cohorts_set_sex_age_menopause(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Full.

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "any"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "any"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort
    records.append(record)


    # Sex.

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "any"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female")
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "male"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "any"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "male")
        ), :
    ]
    records.append(record)

    # Menstruation

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "menstruation"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["menstruation_regular_range"] == 1)
        ), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "premenopause"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["menopause_ordinal"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "perimenopause"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["menopause_ordinal"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "postmenopause"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["menopause_ordinal"] == 2)
        ), :
    ]
    records.append(record)


    # Age

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "young"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["age_grade_female"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "middle"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["age_grade_female"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "female"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "old"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "female") &
            (table_cohort["age_grade_female"] == 2)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "male"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "young"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "male") &
            (table_cohort["age_grade_male"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "male"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "middle"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "male") &
            (table_cohort["age_grade_male"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["cohort_name"] = "space_holder"
    record["cohort_phenotypes"] = "yes"
    record["cohort_genotypes"] = "any"
    record["cohort_sex"] = "male"
    record["cohort_race"] = "any"
    record["cohort_ancestry"] = "any"
    record["cohort_life_stage"] = "old"
    record["cohort_exclusions"] = "none"
    table_cohort = table.copy(deep=True)
    record["table"] = table_cohort.loc[
        (
            (table_cohort["sex_text"] == "male") &
            (table_cohort["age_grade_male"] == 2)
        ), :
    ]
    records.append(record)

    # Return information
    return records


# review: TCW; 23 February 2023
def filter_stratification_cohort_records_by_single_variable_values(
    filter_variable=None,
    filter_values=None,
    cohort_record_variable=None,
    cohort_record_value=None,
    records_cohorts=None,
):
    """
    Filters tables within stratification cohort records by values of a single
    nominal, categorical, or discrete variable.

    arguments:
        filter_variable (str): name of variable (column in table) by which to
            filter records in the stratification cohort table
        filter_values (list): values of variable by which to filter
        cohort_record_variable (str): relevant variable in cohort record to
            change
        cohort_record_value (str): new value to assign to the relevant cohort
            record variable
        records_cohorts (list<dict>): records with information about cohorts

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Iterate on cohort cohorts.
    for record_cohort in records_cohorts:
        # Copy information in table.
        table = record_cohort["table"].copy(deep=True)
        record_cohort["table"] = table.loc[
            (table[filter_variable].isin(filter_values)), :
        ]
        record_cohort[cohort_record_variable] = copy.deepcopy(
            cohort_record_value
        )
        pass
    # Return information.
    return records_cohorts


# review: TCW; 23 February 2023
def filter_stratification_cohort_records_by_exclusions(
    cohort_record_variable=None,
    cohort_record_value=None,
    records_cohorts=None,
):
    """
    Filters tables within stratification cohort records by values of a single
    nominal, categorical, or discrete variable.

    arguments:
        cohort_record_variable (str): relevant variable in cohort record to
            change
        cohort_record_value (str): new value to assign to the relevant cohort
            record variable
        records_cohorts (list<dict>): records with information about cohorts

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Iterate on cohort cohorts.
    for record_cohort in records_cohorts:
        # Copy information in table.
        table = record_cohort["table"].copy(deep=True)
        record_cohort["table"] = table.loc[
            (
                (table["sex_chromosome_aneuploidy"] == 0) &
                (table["sex_discrepancy_identity_genetic"] == 0) &
                (pandas.isna(table["pregnancy"]) | (table["pregnancy"] == 0))
            ), :
        ]
        record_cohort[cohort_record_variable] = copy.deepcopy(
            cohort_record_value
        )
        pass
    # Return information.
    return records_cohorts

# TODO: TCW; 23 February 2023
# TODO: Need a new variable within the cohort records
# TODO: new variable will be "cohort_disorders": "bipolar_disorder", "alcoholism", "alcohol_ever", etc


def filter_previous_stratification_cohort_records(
    records_cohorts=None,
    cohort_name=None,
    cohort_phenotypes=None,
    cohort_genotypes=None,
    cohort_sex=None,
    cohort_race=None,
    cohort_ancestry=None,
    cohort_life_stage=None,
    cohort_exclusions=None,
):
    """
    Filters tables within stratification cohort records by values of a single
    nominal, categorical, or discrete variable.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        cohort_name (list<str>): "None" or values to keep
        cohort_phenotypes (list<str>): "None" or values to keep
        cohort_genotypes (list<str>): "None" or values to keep
        cohort_sex (list<str>): "None" or values to keep
        cohort_race (list<str>): "None" or values to keep
        cohort_ancestry (list<str>): "None" or values to keep
        cohort_life_stage (list<str>): "None" or values to keep
        cohort_exclusions (list<str>): "None" or values to keep

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Prepare entry pairs.
    entry_pairs = dict()
    if cohort_name is not None:
        entry_pairs["cohort_name"] = cohort_name
    if cohort_phenotypes is not None:
        entry_pairs["cohort_phenotypes"] = cohort_phenotypes
    if cohort_genotypes is not None:
        entry_pairs["cohort_genotypes"] = cohort_genotypes
    if cohort_sex is not None:
        entry_pairs["cohort_sex"] = cohort_sex
    if cohort_race is not None:
        entry_pairs["cohort_race"] = cohort_race
    if cohort_ancestry is not None:
        entry_pairs["cohort_ancestry"] = cohort_ancestry
    if cohort_life_stage is not None:
        entry_pairs["cohort_life_stage"] = cohort_life_stage
    if cohort_exclusions is not None:
        entry_pairs["cohort_exclusions"] = cohort_exclusions
    # Filter records.
    records_filter = utility.filter_records_by_multiple_keys(
        entry_pairs=entry_pairs,
        records=records_cohorts,
        report=True,
    )
    # Return information.
    return records_filter


# "ancestry_white_british" <-- will need eventually
# review: TCW; 23 February 2023
def drive_stratify_phenotype_cohorts_set_main(
    table=None,
    report=None,
):
    """
    Stratify phenotype records in cohorts specifically for description tables.

    Variables in cohort record.
    cohort_name: name of cohort   # potential convenience for descriptions
    cohort_phenotypes: "any", "yes", "no"   # availability of records
    cohort_genotypes: "any", "yes", "no"   # availability of records
    cohort_sex: "any", "female", "male"
    cohort_race: "any", "white", etc   # ancestry, race, or ethnicity
    cohort_ancestry: "any", "white", etc   # ancestry from formal genetics
    cohort_life_stage: "any", "young", "middle", "old", "menstruation",
      "premenopause", "perimenopause", "postmenopause", etc
    cohort_exclusions: "none", "sex_aneuploidy;sex_discrepancy;pregnancy", etc

    Variables (columns in table) relevant to each variable in cohort record.
    cohort_phenotypes: always "yes"
    cohort_genotypes: "genotype_availability"
    cohort_sex: "sex_text"
    cohort_race: "race_white"
    cohort_ancestry: "ancestry_white_british"
    cohort_life_stage: "menopause_ordinal", "menstruation_regular_range",
      "age_grade_male"
    cohort_exclusions: "sex_chromosome_aneuploidy",
      "sex_discrepancy_identity_genetic", "pregnancy",

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True) # not actually necessary

    # Collect records of information about each cohort.
    records = list()

    ##########
    # Base

    # Standard, base stratifications by sex and stage of life.
    records_novel = stratify_phenotype_cohorts_set_sex_age_menopause(
        table=table,
    )
    records.extend(records_novel)
    records_base = copy.deepcopy(records_novel)

    # Base; Exclusions;
    records_novel = (
        filter_stratification_cohort_records_by_exclusions(
            cohort_record_variable="cohort_exclusions",
            cohort_record_value="sex_aneuploidy;sex_discrepancy;pregnancy",
            records_cohorts=copy.deepcopy(records_base),
        )
    )
    records.extend(records_novel)
    records_base_exclusions = copy.deepcopy(records_novel)

    # Base; Genotypes;
    records_novel = (
        filter_stratification_cohort_records_by_single_variable_values(
            filter_variable="genotype_availability",
            filter_values=[1,],
            cohort_record_variable="cohort_genotypes",
            cohort_record_value="yes",
            records_cohorts=copy.deepcopy(records_base),
        )
    )
    records.extend(records_novel)
    records_base_genotypes = copy.deepcopy(records_novel)

    # Base; Exclusions; Genotypes;
    records_novel = (
        filter_stratification_cohort_records_by_single_variable_values(
            filter_variable="genotype_availability",
            filter_values=[1,],
            cohort_record_variable="cohort_genotypes",
            cohort_record_value="yes",
            records_cohorts=copy.deepcopy(records_base_exclusions),
        )
    )
    records.extend(records_novel)
    records_base_exclusions_geno = copy.deepcopy(records_novel)

    # Base; Exclusions; Genotypes; Ancestry-White
    records_novel = (
        filter_stratification_cohort_records_by_single_variable_values(
            filter_variable="ancestry_white_british",
            filter_values=[1,],
            cohort_record_variable="cohort_ancestry",
            cohort_record_value="white",
            records_cohorts=copy.deepcopy(records_base_exclusions_geno),
        )
    )
    records.extend(records_novel)
    records_base_excl_geno_antrywhite = copy.deepcopy(records_novel)



    ##########
    # Parent: Base
    if False:

        # Base; Race white.
        records_novel = (
            filter_stratification_cohort_records_by_single_variable_values(
                filter_variable="race_white",
                filter_values=[1,],
                cohort_record_variable="cohort_race",
                cohort_record_value="white",
                records_cohorts=copy.deepcopy(records_base),
            )
        )
        records.extend(records_novel)

        # Base; Race non-white.
        records_novel = (
            filter_stratification_cohort_records_by_single_variable_values(
                filter_variable="race_white",
                filter_values=[0,],
                cohort_record_variable="cohort_race",
                cohort_record_value="non_white",
                records_cohorts=copy.deepcopy(records_base),
            )
        )
        records.extend(records_novel)

    ##########
    # Parent: Base; Genotypes
    if False:

        # Base; genotypes; race white.
        records_novel = (
            filter_stratification_cohort_records_by_single_variable_values(
                filter_variable="race_white",
                filter_values=[1,],
                cohort_record_variable="cohort_race",
                cohort_record_value="white",
                records_cohorts=copy.deepcopy(records_base_genotypes),
            )
        )
        records.extend(records_novel)

        # Base; genotypes; race non-white.
        records_novel = (
            filter_stratification_cohort_records_by_single_variable_values(
                filter_variable="race_white",
                filter_values=[0,],
                cohort_record_variable="cohort_race",
                cohort_record_value="non_white",
                records_cohorts=copy.deepcopy(records_base_genotypes),
            )
        )
        records.extend(records_novel)

    ##########
    # Parent: Base; Genotypes; Exclusions
    if False:

        # Base; exclusions; race white.
        records_novel = (
            filter_stratification_cohort_records_by_single_variable_values(
                filter_variable="race_white",
                filter_values=[1,],
                cohort_record_variable="cohort_race",
                cohort_record_value="white",
                records_cohorts=copy.deepcopy(records_base_geno_exclusions),
            )
        )
        records.extend(records_novel)

        # Base; exclusions; race non-white.
        records_novel = (
            filter_stratification_cohort_records_by_single_variable_values(
                filter_variable="race_white",
                filter_values=[0,],
                cohort_record_variable="cohort_race",
                cohort_record_value="non_white",
                records_cohorts=copy.deepcopy(records_base_geno_exclusions),
            )
        )
        records.extend(records_novel)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "drive_stratify_phenotype_cohorts_set_main()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information
    return records


##########
# Phenotype cohorts
# Special set for plotting functions. Probably temporary until plotting cohort
# specification becomes more sophisticated as in the description tables.
# Review: TCW; 3 May 2023





##########
# Old stuff below here??? ###


def stratify_phenotype_cohorts_set_special_sex_age_menopause(
    column_special=None,
    values_special=None,
    prefix_name=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        column_special (str): name of column in table for variable in special
            selection
        values_special (list): values of variable for special selection
        prefix_name (str): prefix for cohort names
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    record = dict()
    record["name"] = str(str(prefix_name) + "female_male")
    record["cohort"] = str(str(prefix_name) + "female_male")
    record["cohort_model"] = "female_male"
    record["category"] = "sex_together"
    record["phenotype"] = "null"
    record["menstruation"] = False
    #record["table"] = table
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (pandas.isna(table["pregnancy"]) | (table["pregnancy"] == 0))
        ), :
    ]
    records.append(record)

    # Sex

    record = dict()
    record["name"] = str(str(prefix_name) + "female")
    record["cohort"] = str(str(prefix_name) + "female")
    record["cohort_model"] = "female"
    record["category"] = "sex"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male")
    record["cohort"] = str(str(prefix_name) + "male")
    record["cohort_model"] = "male"
    record["category"] = "sex"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "male")
        ), :
    ]
    records.append(record)

    # Menstruation

    record = dict()
    record["name"] = str(str(prefix_name) + "female_menstruation_regular")
    record["cohort"] = "female_menstruation_regular"
    record["cohort_model"] = "female_menstruation_regular"
    record["category"] = "menstruation"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = str(str(prefix_name) + "female_premenopause")
    record["cohort"] = "female_premenopause"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_perimenopause")
    record["cohort"] = "female_perimenopause"
    record["cohort_model"] = "female_perimenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_postmenopause")
    record["cohort"] = "female_postmenopause"
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2)
        ), :
    ]
    records.append(record)

    # Age

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_low")
    record["cohort"] = "female_age_low"
    record["cohort_model"] = "female_age_low"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_middle")
    record["cohort"] = "female_age_middle"
    record["cohort_model"] = "female_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_high")
    record["cohort"] = "female_age_high"
    record["cohort_model"] = "female_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 2)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_low")
    record["cohort"] = "male_age_low"
    record["cohort_model"] = "male_age_low"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_middle")
    record["cohort"] = "male_age_middle"
    record["cohort_model"] = "male_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_high")
    record["cohort"] = "male_age_high"
    record["cohort_model"] = "male_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_phenotype_cohorts_set_special_genotype_sex_age_menopause(
    column_special=None,
    values_special=None,
    prefix_name=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        column_special (str): name of column in table for variable in special
            selection
        values_special (list): values of variable for special selection
        prefix_name (str): prefix for cohort names
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    record = dict()
    record["name"] = str(str(prefix_name) + "female_male")
    record["cohort"] = str(str(prefix_name) + "female_male")
    record["cohort_model"] = "female_male"
    record["category"] = "sex_together"
    record["phenotype"] = "null"
    record["menstruation"] = False
    #record["table"] = table
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (pandas.isna(table["pregnancy"]) | (table["pregnancy"] == 0))
        ), :
    ]
    records.append(record)

    # Sex

    record = dict()
    record["name"] = str(str(prefix_name) + "female")
    record["cohort"] = str(str(prefix_name) + "female")
    record["cohort_model"] = "female"
    record["category"] = "sex"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male")
    record["cohort"] = str(str(prefix_name) + "male")
    record["cohort_model"] = "male"
    record["category"] = "sex"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "male")
        ), :
    ]
    records.append(record)

    # Menstruation

    record = dict()
    record["name"] = str(str(prefix_name) + "female_menstruation_regular")
    record["cohort"] = "female_menstruation_regular"
    record["cohort_model"] = "female_menstruation_regular"
    record["category"] = "menstruation"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = str(str(prefix_name) + "female_premenopause")
    record["cohort"] = "female_premenopause"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_perimenopause")
    record["cohort"] = "female_perimenopause"
    record["cohort_model"] = "female_perimenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_postmenopause")
    record["cohort"] = "female_postmenopause"
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2)
        ), :
    ]
    records.append(record)

    # Age

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_low")
    record["cohort"] = "female_age_low"
    record["cohort_model"] = "female_age_low"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_middle")
    record["cohort"] = "female_age_middle"
    record["cohort_model"] = "female_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_high")
    record["cohort"] = "female_age_high"
    record["cohort_model"] = "female_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 2)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_low")
    record["cohort"] = "male_age_low"
    record["cohort_model"] = "male_age_low"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_middle")
    record["cohort"] = "male_age_middle"
    record["cohort_model"] = "male_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_high")
    record["cohort"] = "male_age_high"
    record["cohort_model"] = "male_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["genotype_availability"] == 1) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_phenotype_cohorts_set_special_white_sex_age_menopause(
    column_special=None,
    values_special=None,
    prefix_name=None,
    table=None,
):
    """
    Stratifies cohorts by sex and stage of life in addition to self report of
    'white' ancestry and/or ethnicity in addition to a special criterion.

    arguments:
        column_special (str): name of column in table for variable in special
            selection
        values_special (list): values of variable for special selection
        prefix_name (str): prefix for cohort names
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    record = dict()
    record["name"] = str(str(prefix_name) + "female_male")
    record["cohort"] = str(str(prefix_name) + "female_male")
    record["cohort_model"] = "female_male"
    record["category"] = "sex_together"
    record["phenotype"] = "null"
    record["menstruation"] = False
    #record["table"] = table
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (pandas.isna(table["pregnancy"]) | (table["pregnancy"] == 0))
        ), :
    ]
    records.append(record)

    # Sex

    record = dict()
    record["name"] = str(str(prefix_name) + "female")
    record["cohort"] = str(str(prefix_name) + "female")
    record["cohort_model"] = "female"
    record["category"] = "sex"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male")
    record["cohort"] = str(str(prefix_name) + "male")
    record["cohort_model"] = "male"
    record["category"] = "sex"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "male")
        ), :
    ]
    records.append(record)

    # Menstruation

    record = dict()
    record["name"] = str(str(prefix_name) + "female_menstruation_regular")
    record["cohort"] = str(str(prefix_name) + "female_menstruation_regular")
    record["cohort_model"] = "female_menstruation_regular"
    record["category"] = "menstruation"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = str(str(prefix_name) + "female_premenopause")
    record["cohort"] =  str(str(prefix_name) + "female_premenopause")
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_perimenopause")
    record["cohort"] = str(str(prefix_name) + "female_perimenopause")
    record["cohort_model"] = "female_perimenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_postmenopause")
    record["cohort"] = str(str(prefix_name) + "female_postmenopause")
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2)
        ), :
    ]
    records.append(record)

    # Age

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_low")
    record["cohort"] = str(str(prefix_name) + "female_age_low")
    record["cohort_model"] = "female_age_low"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_middle")
    record["cohort"] = str(str(prefix_name) + "female_age_middle")
    record["cohort_model"] = "female_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_age_high")
    record["cohort"] = str(str(prefix_name) + "female_age_high")
    record["cohort_model"] = "female_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 2)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_low")
    record["cohort"] = str(str(prefix_name) + "male_age_low")
    record["cohort_model"] = "male_age_low"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_middle")
    record["cohort"] = str(str(prefix_name) + "male_age_middle")
    record["cohort_model"] = "male_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "male_age_high")
    record["cohort"] = str(str(prefix_name) + "male_age_high")
    record["cohort_model"] = "male_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["race_white"] == 1) &
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_phenotype_cohorts_set_special_female_menstruation(
    column_special=None,
    values_special=None,
    prefix_name=None,
    table=None,
):
    """
    Organizes phenotypes for stratification cohorts.

    This collection of stratification cohorts is specific to include female sex,
    to exclude current pregnancy, and to include regular menstruation.

    arguments:
        column_special (str): name of column in table for variable in special
            selection
        values_special (list): values of variable for special selection
        prefix_name (str): prefix for cohort names
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    record = dict()
    record["name"] = str(str(prefix_name) + "female_menstruation_regular")
    record["cohort"] = "female_menstruation_regular"
    record["cohort_model"] = "female_menstruation_regular"
    record["category"] = "menstruation"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_regular_range"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str(str(prefix_name) + "female_premenopause")
    record["cohort"] = "female_premenopause"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table[column_special].isin(values_special)) &
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0)
        ), :
    ]
    records.append(record)

    # Return information
    return records


##########
# Driver functions that combine stratifications of multiple relevant cohorts


def drive_stratify_phenotype_cohorts_set_female_menstruation(
    table=None,
):
    """
    Stratify phenotype records in cohorts.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)

    # Collect records of information about each cohort and model.
    records = list()

    records_novel = (
        stratify_phenotype_cohorts_set_special_female_menstruation(
            column_special="alteration_sex_hormone",
            values_special=[0,],
            prefix_name="alter_hormone_no_",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_female_menstruation(
            column_special="alteration_sex_hormone",
            values_special=[1,],
            prefix_name="alter_hormone_yes_",
            table=table,
        )
    )
    records.extend(records_novel)

    # Return information
    return records


# TODO: TCW; 22 February 2023
# TODO: for manuscript, need to know counts in any ancestry-race-ethnicity, "race_white", and "ancestry_white_british"

# TODO: TCW; 23 February 2023
# implement a separate function to handle all general exclusions:
# 1. "sex_chromosome_aneuploidy"
# 2. "sex_discrepancy_identity_genetic"
# 3. "pregnancy"
# I need to be able to turn the exclusions on and off...




# review: TCW; 24 October 2022
def drive_stratify_phenotype_cohorts_set_description_tables(
    table=None,
    report=None,
):
    """
    Stratify phenotype records in cohorts specifically for description tables.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Collect records of information about each cohort.
    records = list()

    # Standard stratifications by sex, age, and menopause.
    records_novel = (
        stratify_phenotype_cohorts_set_sex_age_menopause(
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="race_white",
            values_special=[1,],
            prefix_name="race_white_",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="race_white",
            values_special=[0,],
            prefix_name="race_other_",
            table=table,
        )
    )
    records.extend(records_novel)

    # Standard stratifications by sex, age, and menopause in records with
    # genotypes available.
    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="genotype_availability",
            values_special=[1,],
            prefix_name="genotype_",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="ancestry_white_british",
            values_special=[1,],
            prefix_name="ancestry_white_british_",
            table=table,
        )
    )
    records.extend(records_novel)

    if True:
        records_novel = (
            stratify_phenotype_cohorts_set_special_sex_age_menopause(
                column_special="bipolar_control_case_strict",
                values_special=[1,],
                prefix_name="bipolar_case_",
                table=table,
            )
        )
        records.extend(records_novel)

        records_novel = (
            stratify_phenotype_cohorts_set_special_sex_age_menopause(
                column_special="bipolar_control_case_strict",
                values_special=[0,],
                prefix_name="bipolar_control_",
                table=table,
            )
        )
        records.extend(records_novel)



    # Stratify cohorts by sex and stage of life in addition to self report of
    # "white" ancestry and/or ethnicity and 'ever' consumption of alcohol.

    records_novel = (
        stratify_phenotype_cohorts_set_special_white_sex_age_menopause(
            column_special="alcohol_ever",
            values_special=[1,],
            prefix_name="race_white_alcohol_ever_",
            table=table,
        )
    )
    records.extend(records_novel)

    # Stratify cohorts by sex and stage of life in addition to self report of
    # "white" ancestry and/or ethnicity and 'ever' consumption of alcohol.

    records_novel = (
        stratify_phenotype_cohorts_set_special_white_sex_age_menopause(
            column_special="alcohol_current",
            values_special=[1,],
            prefix_name="race_white_alcohol_current_",
            table=table,
        )
    )
    records.extend(records_novel)

    # Standard stratifications by sex, age, and menopause in records with
    # genotypes available and with consideration of special categories.

    records_novel = (
        stratify_phenotype_cohorts_set_special_genotype_sex_age_menopause(
            column_special="race_white",
            values_special=[1,],
            prefix_name="race_white_genotype_",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_genotype_sex_age_menopause(
            column_special="race_white",
            values_special=[0,],
            prefix_name="self_other_genotype_",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_genotype_sex_age_menopause(
            column_special="ancestry_white_british",
            values_special=[1,],
            prefix_name="ancestry_white_british_genotype_",
            table=table,
        )
    )
    records.extend(records_novel)


    if False:

        records_novel = (
            stratify_phenotype_cohorts_set_special_sex_age_menopause(
                column_special="alcohol_ever",
                values_special=[1,],
                prefix_name="alcohol_ever_",
                table=table,
            )
        )
        records.extend(records_novel)

        records_novel = (
            stratify_phenotype_cohorts_set_special_sex_age_menopause(
                column_special="alcoholism_control_case_1",
                values_special=[1,],
                prefix_name="alcoholism_case_",
                table=table,
            )
        )
        records.extend(records_novel)

        records_novel = (
            stratify_phenotype_cohorts_set_special_sex_age_menopause(
                column_special="alcoholism_control_case_1",
                values_special=[0,],
                prefix_name="alcoholism_control_",
                table=table,
            )
        )
        records.extend(records_novel)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "drive_stratify_phenotype_cohorts_set_description_tables()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information
    return records


#############
# Old stratification functions now mostly obsolete
##############

# TODO: TCW, 09 February 2022
# TODO: some of the functions below might be obsolete


# obsolete, or not in use
def stratify_phenotype_cohorts_set_sex_body(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Body Mass Index (BMI)

    record = dict()
    record["name"] = "female_body_low"
    record["cohort"] = "female_body_low"
    record["cohort_model"] = "female_body_low"
    record["category"] = "body"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["body_grade_female"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_body_middle"
    record["cohort"] = "female_body_middle"
    record["cohort_model"] = "female_body_middle"
    record["category"] = "body"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["body_grade_female"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_body_high"
    record["cohort"] = "female_body_high"
    record["cohort_model"] = "female_body_high"
    record["category"] = "body"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["body_grade_female"] == 2)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_body_low"
    record["cohort"] = "male_body_low"
    record["cohort_model"] = "male_body_low"
    record["category"] = "body"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["body_grade_male"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_body_middle"
    record["cohort"] = "male_body_middle"
    record["cohort_model"] = "male_body_middle"
    record["category"] = "body"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["body_grade_male"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_body_high"
    record["cohort"] = "male_body_high"
    record["cohort_model"] = "male_body_high"
    record["category"] = "body"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["body_grade_male"] == 2)
        ), :
    ]
    records.append(record)

    # Return information
    return records


# obsolete, or not in use
def stratify_phenotype_cohorts_set_regression(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)

    # Collect records of information about each cohort and model.
    records = list()
    records_novel = (
        stratify_phenotype_cohorts_set_sex_age_menopause(
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="season_summer_autumn",
            values_special=[1,],
            name_special="season_summer_autumn",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="season_summer_autumn",
            values_special=[0,],
            name_special="season_fall_winter",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_phenotype_cohorts_set_special_sex_age_menopause(
            column_special="alcohol_ever",
            values_special=[1,],
            name_special="alcohol_ever",
            table=table,
        )
    )
    records.extend(records_novel)

    # Return information
    return records


def stratify_set_secondary_menopause(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Menopause binary

    record = dict()
    record["name"] = "female_premenopause_binary"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_binary"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_binary"
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_binary"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary"] == 1)
        ), :
    ]
    records.append(record)

    # Menopause binary strict

    record = dict()
    record["name"] = "female_premenopause_binary_strict"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_binary_strict"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary_strict"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_binary_strict"
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_binary_strict"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary_strict"] == 1)
        ), :
    ]
    records.append(record)


    # Menopause ordinal strict

    record = dict()
    record["name"] = "female_premenopause_ordinal_strict"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal_strict"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal_strict"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_ordinal_strict"
    record["cohort_model"] = "female_perimenopause"
    record["category"] = "menopause_ordinal_strict"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal_strict"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_ordinal_strict"
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_ordinal_strict"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal_strict"] == 2)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_set_female_menstruation(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Menstruation

    record = dict()
    record["name"] = "female_premenopause_menstruation_follicular"
    record["cohort_model"] = "female_premenopause_menstruation_follicular"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["menstruation_phase"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_menstruation_luteal"
    record["cohort_model"] = "female_premenopause_menstruation_luteal"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["menstruation_phase"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_menstruation_cycle_shoulder"
    record["cohort_model"] = "female_premenopause_menstruation_cycle_shoulder"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["menstruation_phase_cycle"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_menstruation_cycle_ovulation"
    record["cohort_model"] = "female_premenopause_menstruation_cycle_ovulation"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["menstruation_phase_cycle"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_menstruation_follicular"
    record["cohort_model"] = "female_perimenopause_menstruation_follicular"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["menstruation_phase"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_menstruation_luteal"
    record["cohort_model"] = "female_perimenopause_menstruation_luteal"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["menstruation_phase"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_menstruation_cycle_shoulder"
    record["cohort_model"] = "female_perimenopause_menstruation_cycle_shoulder"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["menstruation_phase_cycle"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_menstruation_cycle_ovulation"
    record["cohort_model"] = "female_perimenopause_menstruation_cycle_ovulation"
    record["category"] = "menstruation_phase"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["menstruation_phase_cycle"] == 1)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_set_female_hormone_alteration(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Hormone-alteration therapies

    record = dict()
    record["name"] = "female_oral_contraception_yes"
    record["cohort_model"] = "female_oral_contraception_yes"
    record["category"] = "oral_contraception"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["oral_contraception"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_oral_contraception_no"
    record["cohort_model"] = "female_oral_contraception_no"
    record["category"] = "oral_contraception"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["oral_contraception"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_hormone_replacement_yes"
    record["cohort_model"] = "female_hormone_replacement_yes"
    record["category"] = "hormone_replacement"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_replacement"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_hormone_replacement_no"
    record["cohort_model"] = "female_hormone_replacement_no"
    record["category"] = "hormone_replacement"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_replacement"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_hormone_alteration_yes"
    record["cohort_model"] = "female_hormone_alteration_yes"
    record["category"] = "hormone_alteration"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_hormone_alteration_no"
    record["cohort_model"] = "female_hormone_alteration_no"
    record["category"] = "hormone_alteration"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    records.append(record)

    # Hormone-altering therapies in premenopausal females

    record = dict()
    record["name"] = "female_premenopause_oral_contraception_yes"
    record["cohort_model"] = "female_premenopause_oral_contraception_yes"
    record["category"] = "oral_contraception"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["oral_contraception"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_oral_contraception_no"
    record["cohort_model"] = "female_premenopause_oral_contraception_no"
    record["category"] = "oral_contraception"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["oral_contraception"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_hormone_replacement_yes"
    record["cohort_model"] = "female_premenopause_hormone_replacement_yes"
    record["category"] = "hormone_replacement"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["hormone_replacement"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_hormone_replacement_no"
    record["cohort_model"] = "female_premenopause_hormone_replacement_no"
    record["category"] = "hormone_replacement"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["hormone_replacement"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_hormone_alteration_yes"
    record["cohort_model"] = "female_premenopause_hormone_alteration_yes"
    record["category"] = "hormone_alteration"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_premenopause_hormone_alteration_no"
    record["cohort_model"] = "female_premenopause_hormone_alteration_no"
    record["category"] = "hormone_alteration"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    records.append(record)

    # Hormone-altering therapies in perimenopausal females

    record = dict()
    record["name"] = "female_perimenopause_oral_contraception_yes"
    record["cohort_model"] = "female_perimenopause_oral_contraception_yes"
    record["category"] = "oral_contraception"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["oral_contraception"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_oral_contraception_no"
    record["cohort_model"] = "female_perimenopause_oral_contraception_no"
    record["category"] = "oral_contraception"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["oral_contraception"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_hormone_replacement_yes"
    record["cohort_model"] = "female_perimenopause_hormone_replacement_yes"
    record["category"] = "hormone_replacement"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["hormone_replacement"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_hormone_replacement_no"
    record["cohort_model"] = "female_perimenopause_hormone_replacement_no"
    record["category"] = "hormone_replacement"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["hormone_replacement"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_hormone_alteration_yes"
    record["cohort_model"] = "female_perimenopause_hormone_alteration_yes"
    record["category"] = "hormone_alteration"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_hormone_alteration_no"
    record["cohort_model"] = "female_perimenopause_hormone_alteration_no"
    record["category"] = "hormone_alteration"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    records.append(record)

    # Hormone-altering therapies in postmenopausal females

    record = dict()
    record["name"] = "female_postmenopause_oral_contraception_yes"
    record["cohort_model"] = "female_postmenopause_oral_contraception_yes"
    record["category"] = "oral_contraception"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["oral_contraception"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_oral_contraception_no"
    record["cohort_model"] = "female_postmenopause_oral_contraception_no"
    record["category"] = "oral_contraception"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["oral_contraception"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_hormone_replacement_yes"
    record["cohort_model"] = "female_postmenopause_hormone_replacement_yes"
    record["category"] = "hormone_replacement"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["hormone_replacement"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_hormone_replacement_no"
    record["cohort_model"] = "female_postmenopause_hormone_replacement_no"
    record["category"] = "hormone_replacement"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["hormone_replacement"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_hormone_alteration_yes"
    record["cohort_model"] = "female_postmenopause_hormone_alteration_yes"
    record["category"] = "hormone_alteration"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_hormone_alteration_no"
    record["cohort_model"] = "female_postmenopause_hormone_alteration_no"
    record["category"] = "hormone_alteration"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_set_female_pregnancy(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Pregnancy

    record = dict()
    record["name"] = "female_pregnancy_definite"
    record["cohort_model"] = "female_pregnancy_definite"
    record["category"] = "pregnancy"
    record["menstruation"] = False
    record["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["3140-0.0"] == 1)), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_pregnancy_unsure"
    record["cohort_model"] = "female_pregnancy_unsure"
    record["category"] = "pregnancy"
    record["menstruation"] = False
    record["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["3140-0.0"] == 2)), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_pregnancy_yes"
    record["cohort_model"] = "female_pregnancy_yes"
    record["category"] = "pregnancy"
    record["menstruation"] = False
    record["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["pregnancy"] == 1)), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_pregnancy_no"
    record["cohort_model"] = "female_pregnancy_no"
    record["category"] = "pregnancy"
    record["menstruation"] = False
    record["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["pregnancy"] == 0)), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = "female_menopause_unsure"
    record["cohort_model"] = "female_menopause_unsure"
    record["category"] = "menopause"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["2724-0.0"] == 3)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_menopause_blank"
    record["cohort_model"] = "female_menopause_blank"
    record["category"] = "menopause"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["2724-0.0"] == -3)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_set_alcohol_sex_menopause_age(
    alcohol_variable=None,
    alcohol_value=None,
    cohort_suffix=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    This function allows selection of cohorts on the basis of a variety of
    alcohol variables.
    Examples include current or previous alcohol consumption and designation as
    a control or case for alcoholism.

    arguments:
        alcohol_variable (str): name of column for an alcohol variable of
            interest for cohort stratification
        alcohol_value (float): a single specific value of the alcohol variable
            for cohort stratification
        cohort_suffix (str): suffix for cohort name to indicate alcohol variable
            stratification
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    #####
    # Sex

    record = dict()
    record["name"] = str("female" + cohort_suffix)
    record["cohort"] = str("female" + cohort_suffix)
    record["cohort_model"] = "female"
    record["category"] = "sex_alcohol"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male" + cohort_suffix)
    record["cohort"] = str("male" + cohort_suffix)
    record["cohort_model"] = "male"
    record["category"] = "sex_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    #####
    # Menopause

    record = dict()
    record["name"] = str("female_premenopause" + cohort_suffix)
    record["cohort"] = str("female_premenopause" + cohort_suffix)
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal_alcohol"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("female_perimenopause" + cohort_suffix)
    record["cohort"] = str("female_perimenopause" + cohort_suffix)
    record["cohort_model"] = "female_perimenopause"
    record["category"] = "menopause_ordinal_alcohol"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("female_postmenopause" + cohort_suffix)
    record["cohort"] = str("female_postmenopause" + cohort_suffix)
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_ordinal_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    #####
    # Age

    record = dict()
    record["name"] = str("female_age_low" + cohort_suffix)
    record["cohort"] = str("female_age_low" + cohort_suffix)
    record["cohort_model"] = "female_age_low"
    record["category"] = "age_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 0) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("female_age_middle" + cohort_suffix)
    record["cohort"] = str("female_age_middle" + cohort_suffix)
    record["cohort_model"] = "female_age_middle"
    record["category"] = "age_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 1) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("female_age_high" + cohort_suffix)
    record["cohort"] = str("female_age_high" + cohort_suffix)
    record["cohort_model"] = "female_age_high"
    record["category"] = "age_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 2) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male_age_low" + cohort_suffix)
    record["cohort"] = str("male_age_low" + cohort_suffix)
    record["cohort_model"] = "male_age_low"
    record["category"] = "age_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 0) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male_age_middle" + cohort_suffix)
    record["cohort"] = str("male_age_middle" + cohort_suffix)
    record["cohort_model"] = "male_age_middle"
    record["category"] = "age_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 1) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male_age_high" + cohort_suffix)
    record["cohort"] = str("male_age_high" + cohort_suffix)
    record["cohort_model"] = "male_age_high"
    record["category"] = "age_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    #####
    # Body Mass Index (BMI)

    record = dict()
    record["name"] = str("female_body_low" + cohort_suffix)
    record["cohort"] = str("female_body_low" + cohort_suffix)
    record["cohort_model"] = "female_body_low"
    record["category"] = "body_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["body_grade_female"] == 0) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("female_body_middle" + cohort_suffix)
    record["cohort"] = str("female_body_middle" + cohort_suffix)
    record["cohort_model"] = "female_body_middle"
    record["category"] = "body_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["body_grade_female"] == 1) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("female_body_high" + cohort_suffix)
    record["cohort"] = str("female_body_high" + cohort_suffix)
    record["cohort_model"] = "female_body_high"
    record["category"] = "body_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["body_grade_female"] == 2) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male_body_low" + cohort_suffix)
    record["cohort"] = str("male_body_low" + cohort_suffix)
    record["cohort_model"] = "male_body_low"
    record["category"] = "body_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["body_grade_male"] == 0) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male_body_middle" + cohort_suffix)
    record["cohort"] = str("male_body_middle" + cohort_suffix)
    record["cohort_model"] = "male_body_middle"
    record["category"] = "body_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["body_grade_male"] == 1) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = str("male_body_high" + cohort_suffix)
    record["cohort"] = str("male_body_high" + cohort_suffix)
    record["cohort_model"] = "male_body_high"
    record["category"] = "body_alcohol"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["body_grade_male"] == 2) &
            (table[alcohol_variable] == alcohol_value)
        ), :
    ]
    records.append(record)

    # Return information
    return records


def stratify_cohorts_models_phenotypes_sets(
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()
    records_novel = (
        stratify_set_primary_sex_age_body_menopause(
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcohol_ever",
            alcohol_value=1,
            cohort_suffix="_alcohol_ever",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcohol_ever",
            alcohol_value=0,
            cohort_suffix="_alcohol_never",
            table=table,
        )
    )
    records.extend(records_novel)

    # name: "alcoholism_control_case_any"
    # NA: missing value indicates neither control nor case
    # 0: control
    # 1: case
    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcoholism_control_case_any",
            alcohol_value=0,
            cohort_suffix="_alcoholism_any_control",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcoholism_control_case_any",
            alcohol_value=1,
            cohort_suffix="_alcoholism_any_case",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcoholism_control_case_1",
            alcohol_value=0,
            cohort_suffix="_alcoholism_1_control",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcoholism_control_case_1",
            alcohol_value=1,
            cohort_suffix="_alcoholism_1_case",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcoholism_control_case_2",
            alcohol_value=0,
            cohort_suffix="_alcoholism_2_control",
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            alcohol_variable="alcoholism_control_case_2",
            alcohol_value=1,
            cohort_suffix="_alcoholism_2_case",
            table=table,
        )
    )
    records.extend(records_novel)

    ##########

    records_novel = (
        stratify_set_secondary_menopause(
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_female_menstruation(
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_female_hormone_alteration(
            table=table,
        )
    )
    records.extend(records_novel)

    records_novel = (
        stratify_set_female_pregnancy(
            table=table,
        )
    )
    records.extend(records_novel)

    # Return information
    return records



###############################################################################
# Procedure



#
