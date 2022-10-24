"""
Stratify cohorts for analyses on genotypes from the U.K. Biobank.

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

# TODO: TCW; 11 October 2022
# TODO: I think it would be helpful to have a new module named "genotype".
# TODO: Within this module, control the stratification of phenotype records with matching gentoype records...
# TODO: Try to keep the stratification and preparation of these genotype cohort tables modular and versatile.
# TODO: 1. stratify the tables with selection of non-missing values for variables of interest.
# TODO: - - drive this stratification for a group of cohort tables AND variables of interest... "stratification parameter" records?
# TODO: - - each record should include a list of variables to include in the table
# TODO: - - each record should also include information about selection of female or male cohorts (menopause, stage of life, pregnancy, etc)
# TODO: - - this should probably look like a separate dictionary with variable names (keys) and acceptable values within a list
# TODO: - - or this could include a specific and set collection of cohort variables (sex, pregnancy, menopause, life stage)
# TODO: - - each record should also include a flag for whether or not to select unrelated persons to represent the cohort



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
    paths["genotype"] = os.path.join(path_dock, "genotype")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["genotype"])
    # Initialize directories.
    utility.create_directories(
        path=paths["genotype"]
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


def read_source_table_kinship_pairs(
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
        (object): Pandas data frame of kinship coefficients across pairs of
            persons in UK Biobank cohort

    """

    # Specify directories and files.
    path_table_kinship_pairs = os.path.join(
        path_dock, "assembly", "table_kinship_pairs.pickle"
    )
    # Read information from file.
    table_kinship_pairs = pandas.read_pickle(
        path_table_kinship_pairs
    )

    # Report.
    if report:
        # Print report.
        utility.print_terminal_partition(level=2)
        print(
            "report: " +
            "read_source_table_kinship_pairs()"
        )
        describe_report_table_kinship_pairs(
            threshold_kinship=0.1,
            table_kinship_pairs=table_kinship_pairs,
        )
        pass
    # Return information.
    return table_kinship_pairs


def read_source_table_stratification_cohorts_models(
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
        (object): Pandas data frame of kinship coefficients across pairs of
            persons in UK Biobank cohort

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "stratification_genotype_cohorts_models",
        "table_stratification_genotype_cohorts_models.tsv",
    )
    # Read information from file.
    table = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype={
            "execution": "int", # logical binary of whether to execute
            "group": "string",
            "name": "string", # name of cohort table
            "identity_self_white": "string", # 'None' or strat values
            "ancestry_white_british": "string", # 'None' or strat values
            "bipolar_disorder": "string", # 'None' or stratification values
            "alcohol_ever": "string", # 'None' or stratification values
            "alcohol_current": "string", # 'None' or stratification values
            "alcohol_moderate": "string", # 'None' or stratification values
            "sex_text": "string", # 'None' or stratification values
            "female_age_grade": "string", # stratification variable values
            "female_menopause": "string", # stratification variable values
            "female_menstruation": "string", # stratification variable values
            "female_pregnancy": "string", # stratification variable values
            "male_age_grade": "string", # stratification variable values
            "kinship_unrelated": "int", # whether to select unrelated kinship
            "priority_female": "int", # whether to prioritize female kinship
            "dependence": "string", # names of dependent variables
            "dependence_scale_transforms": "int", # whether to transform
            "dependence_logical_binary": "int", # whether to translate to PLINK
            "variables_nonmissing": "string", # variables to select non-missing
            "variables_prefix_nonmissing": "string", # variables by prefix
            "independence": "string", # names of independent variables
            "independence_prefix": "string", # independent variables by prefix
            "independence_scale": "int", # whether to transform to z standard
        },
    )

    # Report.
    if report:
        # Print report.
        utility.print_terminal_partition(level=2)
        print(
            "report: " +
            "read_source_table_stratification_cohorts_models()"
        )
        print(table)
        pass
    # Return information.
    return table


##########
# Kinship


def describe_report_table_kinship_pairs(
    threshold_kinship=None,
    table_kinship_pairs=None,
):
    """
    Describes a table of kinship between pairs of persons.

    arguments:
        threshold_kinship (float): maximal value of kinship coefficient
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons

    raises:

    returns:

    """

    # Copy information in table.
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)
    # Summarize values of the Kinship Coefficient.
    array_kinship = copy.deepcopy(
        table_kinship_pairs["Kinship"].dropna().to_numpy()
    )
    median = numpy.nanmedian(array_kinship)
    minimum = numpy.nanmin(array_kinship)
    maximum = numpy.nanmax(array_kinship)
    # Summarize pairs for different filters on Kinship Coefficient.
    count_pairs_total = copy.deepcopy(table_kinship_pairs.shape[0])

    table_negative = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] < 0), :
    ]
    count_pairs_negative = table_negative.shape[0]

    table_positive = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] >= 0), :
    ]
    count_pairs_positive = table_positive.shape[0]

    table_below = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] < threshold_kinship), :
    ]
    count_pairs_below = table_below.shape[0]


    table_above = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] >= threshold_kinship), :
    ]
    count_pairs_above = table_above.shape[0]

    # Summarize values of the Kinship Coefficient that are positive.
    array_positive = copy.deepcopy(
        table_positive["Kinship"].dropna().to_numpy()
    )
    median_positive = numpy.nanmedian(array_positive)
    minimum_positive = numpy.nanmin(array_positive)
    maximum_positive = numpy.nanmax(array_positive)

    # Print report.
    utility.print_terminal_partition(level=2)
    print(
        "report: " +
        "read_source_table_kinship_pairs()"
    )
    utility.print_terminal_partition(level=5)
    print("... kinship coefficients across total pairs ...")
    print("median kinship coefficient: " + str(median))
    print("minimum kinship coefficient: " + str(minimum))
    print("maximum kinship coefficient: " + str(maximum))
    utility.print_terminal_partition(level=5)
    print("... kinship coefficients across pairs with kinship >= 0 ...")
    print("median kinship coefficient: " + str(median_positive))
    print("minimum kinship coefficient: " + str(minimum_positive))
    print("maximum kinship coefficient: " + str(maximum_positive))
    utility.print_terminal_partition(level=5)
    print("... kinship pairs ...")
    print("total pairs: " + str(count_pairs_total))
    print("pairs with kinship < 0: " + str(count_pairs_negative))
    print("pairs with kinship >= 0: " + str(count_pairs_positive))
    utility.print_terminal_partition(level=5)
    print("Kinship Coefficient threshold: " + str(threshold_kinship))
    print("pairs with kinship < threshold: " + str(count_pairs_below))
    print("pairs with kinship >= threshold: " + str(count_pairs_above))

    pass


##########
# Preparation


def extract_stratification_record_parameters(
    record=None,
):
    """
    Extracts information about parameters for stratification of a single cohort
    for analyses on genotypes.

    Some columns accommodate a single integer value for logical binary
    interpretation.
    Some columns accommodate multiple integer values indicative of categorical
    groups.
    Some columns accommodate multiple string values indicative of categorical
    groups or column names for specific variables.

    arguments:
        record (dict): information about stratification parameters

    raises:

    returns:
        (dict): parameters for stratification of a single cohort

    """

    # Nested function.
    def extract_variable_values_string(
        field=None,
    ):
        if (
            (not pandas.isna(field)) and
            (len(str(field)) > 0) and
            (str(field).strip() != "None")
        ):
            values = str(field).split(";")
        else:
            values = None
        return values

    def extract_variable_values_integer(
        field=None,
    ):
        values_string = extract_variable_values_string(field=field)
        if (values_string is not None):
            values_integer = list(map(
                lambda value_string: int(value_string),
                values_string
            ))
        else:
            values_integer = None
        return values_integer

    # Copy information.
    record = copy.deepcopy(record)
    # Extract and collect parameters from record.
    pail = dict()
    pail["group"] = str(record["group"])
    pail["name"] = str(record["name"])
    pail["identity_self_white"] = extract_variable_values_integer(
        field=record["identity_self_white"]
    )
    pail["ancestry_white_british"] = extract_variable_values_integer(
        field=record["ancestry_white_british"]
    )
    pail["bipolar_disorder"] = extract_variable_values_integer(
        field=record["bipolar_disorder"]
    )
    pail["alcohol_ever"] = extract_variable_values_integer(
        field=record["alcohol_ever"]
    )
    pail["alcohol_current"] = extract_variable_values_integer(
        field=record["alcohol_current"]
    )
    pail["alcohol_moderate"] = extract_variable_values_integer(
        field=record["alcohol_moderate"]
    )
    pail["sex_text"] = extract_variable_values_string(
        field=record["sex_text"]
    )
    pail["female_age_grade"] = extract_variable_values_integer(
        field=record["female_age_grade"]
    )
    pail["female_menopause"] = extract_variable_values_integer(
        field=record["female_menopause"]
    )
    pail["female_menstruation"] = extract_variable_values_integer(
        field=record["female_menstruation"]
    )
    pail["female_pregnancy"] = extract_variable_values_integer(
        field=record["female_pregnancy"]
    )
    pail["male_age_grade"] = extract_variable_values_integer(
        field=record["male_age_grade"]
    )
    pail["kinship_unrelated"] = int(record["kinship_unrelated"])
    pail["priority_female"] = int(record["priority_female"])
    pail["dependence"] = extract_variable_values_string(
        field=record["dependence"]
    )
    pail["dependence_scale_transforms"] = int(
        record["dependence_scale_transforms"]
    )
    pail["dependence_logical_binary"] = int(record["dependence_logical_binary"])
    pail["variables_nonmissing"] = extract_variable_values_string(
        field=record["variables_nonmissing"]
    )
    pail["variables_prefix_nonmissing"] = extract_variable_values_string(
        field=record["variables_prefix_nonmissing"]
    )
    pail["independence"] = extract_variable_values_string(
        field=record["independence"]
    )
    pail["independence_prefix"] = extract_variable_values_string(
        field=record["independence_prefix"]
    )
    pail["independence_scale"] = int(record["independence_scale"])

    # Return information.
    return pail


##########
# Selections on samples for both sexes



def control_stratify_by_sex_together_categories(
    identity_self_white=None,
    ancestry_white_british=None,
    bipolar_disorder=None,
    alcohol_ever=None,
    alcohol_current=None,
    alcohol_moderate=None,
    table=None,
    report=None,
):
    """
    Control filters, selections, stratifications on all samples (rows) within a
    table regardless of sex.

    arguments:
        identity_self_white (list<int>): categorical values for selection
        ancestry_white_british (list<int>): categorical values for selection
        bipolar_disorder (list<int>): categorical values for selection
        alcohol_ever (list<int>): categorical values for selection
        alcohol_current (list<int>): categorical values for selection
        alcohol_moderate (list<int>): categorical values for selection
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Select samples (rows) within a table.
    selections = dict()
    selections["ancestry_self_white"] = identity_self_white
    selections["white_british"] = ancestry_white_british
    selections["bipolar_control_case_strict"] = bipolar_disorder
    selections["alcohol_ever"] = alcohol_ever
    selections["alcohol_current"] = alcohol_current
    selections["alcohol_moderate"] = alcohol_moderate
    for selection in selections.keys():
        if (
            (not pandas.isna(selections[selection])) and
            (selections[selection] is not None)
        ):
            table = table.loc[
                (table[selection].isin(selections[selection])), :
            ]
            pass
        pass

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "control_stratify_by_sex_together_categories()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table




##########
# Control

# 2. select records by ancestry ("ancestry_white_british") if variable is not None.
# 2. select records by sex (variable column "sex_text")
# 3. select records by other variable values on both sexes
# - - Remember that the selection parameters are strings! Convert these strings to int before comparison.
# - - Example: Bipolar Disorder cases and controls
# - - Only filter by these variables at all if they do not have values of "None" in the parameter table.
# - - Important to avoid losing records due to any missing values in selection variables.
# 3. if selection includes both sexes, then split phenotype table records by sex
# 4. stratify cohorts by sex-specific variables and values (such as menopause)
# - - Only filter by these variables at all if they do not have values of "None" in the parameter table.
# - - Important to avoid losing records due to any missing values in selection variables.
# - - Example: apply each filter selection sequentially after checking the parameter variable in an "if then" block.
# 5. recombine records for female and male sex
# 6. select columns from union of "dependence", "variables_nonmissing", and "variables_prefix_nonmissing"
# 7. select records by non-missing values in union of "dependence", "variables_nonmissing", and "variables_prefix_nonmissing"
# 8. select records for unrelated kinship with any priority in random selection
# 9. apply any transformations (scale or logical binary) to variables "dependence"



def control_stratify_genotype_cohorts(
    group=None,
    name=None,
    identity_self_white=None,
    ancestry_white_british=None,
    bipolar_disorder=None,
    alcohol_ever=None,
    alcohol_current=None,
    alcohol_moderate=None,
    sex_text=None,
    female_age_grade=None,
    female_menopause=None,
    female_menstruation=None,
    female_pregnancy=None,
    male_age_grade=None,
    kinship_unrelated=None,
    priority_female=None,
    dependence=None,
    dependence_scale_transforms=None,
    dependence_logical_binary=None,
    variables_nonmissing=None,
    variables_prefix_nonmissing=None,
    independence=None,
    independence_prefix=None,
    independence_scale=None,
    table=None,
    report=None,
):
    """
    Control the preparation of a table and record for single genotype
    stratification cohort.

    arguments:
        group (str): name of group of stratification cohort tables
        name (str): name of stratification cohort table
        identity_self_white (list<int>): categorical values for selection
        ancestry_white_british (list<int>): categorical values for selection
        bipolar_disorder (list<int>): categorical values for selection
        alcohol_ever (list<int>): categorical values for selection
        alcohol_current (list<int>): categorical values for selection
        alcohol_moderate (list<int>): categorical values for selection
        sex_text (list<str>): categorical values for selection
        female_age_grade (list<int>): categorical values for selection
        female_menopause (list<int>): categorical values for selection
        female_pregnancy (list<int>): categorical values for selection
        male_age_grade (list<int>): categorical values for selection
        kinship_unrelated (int): logical binary whether to select persons with
            unrelated kinship
        priority_female (int): logical binary whether to give priority to
            persons with female sex in selection of unrelated kinship
        dependence (list<str>): names of table's columns for dependent variables
        dependence_scale_transforms (int): logical binary whether to apply
            distribution scale transformations to values of dependent variables
        dependence_logical_binary (int): logical binary whether to translate
            values of logical binary dependent variable to code for PLINK2 with
            zero mapping to one (0 --> 1; control) and one mapping to two
            (1 --> 2; case)
        variables_nonmissing (list<str>): names of table's columns for variables
            for which to select non-missing values
        variables_prefix_nonmissing (list<str>): prefixes of names of table's
            columns for variables for which to select non-missing values
        independence (list<str>): names of table's columns for independent
            variables
        independence_prefix (list<str>): prefixes of names of table's columns
            for independent variables
        independence_scale (int): logical binary whether to apply Standard Z
            Score transformation to values of independent variables
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["group"] = str(group)
    record["name"] = str(name)

    # Control filters, selections, stratifications on all samples (rows)
    # regardless of sex.
    table = control_stratify_by_sex_together_categories(
        identity_self_white=identity_self_white,
        ancestry_white_british=ancestry_white_british,
        bipolar_disorder=bipolar_disorder,
        alcohol_ever=alcohol_ever,
        alcohol_current=alcohol_current,
        alcohol_moderate=alcohol_moderate,
        table=table,
        report=report,
    )

    # Control filters, selections, stratifications on samples (rows) that are
    # specific to female or male sex.
    table = control_stratify_by_sex_separate_categories(
        sex_text=sex_text,
        female_age_grade=female_age_grade,
        female_menopause=female_menopause,
        female_menstruation=female_menstruation,
        female_pregnancy=pregnancy,
        male_age_grade=male_age_grade,
        table=table,
        report=report,
    )


    # Control filters, selections, stratifications on samples (rows) by the
    # non-missingness of values of variables (columns).


    # Control filters, selections, stratifications on samples (rows) by the
    # kinship coefficient (relatedness) between persons in the cohort.



    # Control filters, selections, stratifications on variables (columns) in the
    # final table.

    # Collect in record the cohort table after filters, selections, and
    # stratifications on samples (rows) and variables (columns).
    record["table"] = table

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "control_stratify_genotype_cohorts()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return record




##########
# Write


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
    print("version check: TCW; 11 October 2022")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    # Read source information from file.
    source = read_source(
        path_dock=path_dock,
        report=True,
    )
    table_kinship_pairs = read_source_table_kinship_pairs(
        path_dock=path_dock,
        report=True,
    )
    table_stratification = read_source_table_stratification_cohorts_models(
        path_dock=path_dock,
        report=True,
    )
    # Filter stratification table by "execution".
    if True:
        table_stratification = table_stratification.loc[
            (
                (table_stratification["execution"] == 1)
            ), :
        ]
        pass

    # Collect records of information about each stratification cohort.
    records_cohorts = list()
    # Iterate on records in stratification parameter table.
    # Each record (row) in stratification parameter table specifies parameters
    # for preparation of a table and record for a single genotype stratification
    # cohort.
    records_stratification = table_stratification.to_dict(
        orient="records",
    )
    for record_stratification in records_stratification:
        # Extract stratification parameter information from record.
        pail = extract_stratification_record_parameters(
            record=record_stratification,
        )
        # Drive regression on single cohort and model.
        record_cohort = control_stratify_genotype_cohorts(
            group=pail["group"],
            name=pail["name"],
            identity_self_white=pail["identity_self_white"],
            ancestry_white_british=pail["ancestry_white_british"],
            bipolar_disorder=pail["bipolar_disorder"],
            alcohol_ever=pail["alcohol_ever"],
            alcohol_current=pail["alcohol_current"],
            alcohol_moderate=pail["alcohol_moderate"],
            sex_text=pail["sex_text"],
            female_age_grade=pail["female_age_grade"],
            female_menopause=pail["female_menopause"],
            female_menstruation=pail["female_menstruation"],
            female_pregnancy=pail["female_pregnancy"],
            male_age_grade=pail["male_age_grade"],
            kinship_unrelated=pail["kinship_unrelated"],
            priority_female=pail["priority_female"],
            dependence=pail["dependence"],
            dependence_scale_transforms=pail["dependence_scale_transforms"],
            dependence_logical_binary=pail["dependence_logical_binary"],
            variables_nonmissing=pail["variables_nonmissing"],
            variables_prefix_nonmissing=pail["variables_prefix_nonmissing"],
            independence=pail["independence"],
            independence_prefix=pail["independence_prefix"],
            independence_scale=pail["independence_scale"],
            table=source["table_phenotypes"],
            report=False,
        )
        # Collect records.
        records_cohorts.append(record_cohort)
        pass

    # TODO: TCW; 11 October 2022
    # TODO: now write those 'records_cohorts' to file...
    # TODO: the record's "group" value specifies the name of the directory
    # TODO: the record's "name" value specifies the name of the "table_*.tsv" file


    #############old procedure###################

    if False:

        # Extract names of sets.
        stratification_sets = list(set(table_stratification["set"].to_list()))
        if (True):
            stratification_sets.insert(0, "reference_population")
        # Iterate on sets of stratification cohorts and models.
        for stratification_set in stratification_sets:
            # Stratify cohorts and models in set.
            pail_set = execute_organize_stratification_genotype_cohorts_models(
                table=source["table_phenotypes"],
                table_kinship_pairs=table_kinship_pairs,
                table_stratification=table_stratification,
                stratification_set=stratification_set,
                variance_scale=False, # whether to standardize variance scale
                format_plink=True, # whether to convert table formats for PLINK2
                report=True,
            )
            # Write product information to file.
            write_genotype_product(
                pail_write=pail_set,
                directory=stratification_set,
                paths=paths,
            )
        pass
    pass


#
