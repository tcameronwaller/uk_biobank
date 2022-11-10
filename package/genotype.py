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
            "priority_male": "int", # whether to prioritize male kinship
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
            (len(str(field).strip()) > 0) and
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
    pail["priority_male"] = int(record["priority_male"])
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
# Selections on samples for both sexes together


# TODO: TCW; 7 November 2022
# TODO: change variable name "ancestry_self_white" to "identity_self_white"
def stratify_by_sex_together_categories(
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
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

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
        #if (selections[selection] is None):
        #    print(str(selection) + " is None!")
        pass

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "stratify_by_sex_together_categories()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


def control_report_category_values_in_cohort(
    cohort=None,
    variable=None,
    records_cohorts=None,
):
    """
    Control report on the values of a categorical variable that remain within a
    specific cohort after stratification.

    arguments:
        cohort (str): name of a single cohort
        variable (str): name of a single variable with categorical values
        records_cohorts (list<dict>): records of information about cohorts

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=2)
    print("report: ")
    name_function = (
        "control_report_category_values_in_cohort()"
    )
    print(name_function)
    print("cohort: " + str(cohort))
    print("variable: " + str(variable))
    utility.print_terminal_partition(level=3)

    # Organize hashable entries for cohorts.
    entries_cohorts = utility.structure_from_records_to_entries(
        records=records_cohorts,
    )
    # Extract table for a specific cohort.
    if (cohort in entries_cohorts.keys()):
        # Copy information in table.
        table = entries_cohorts[cohort]["table"].copy(deep=True)
    else:
        table = pandas.DataFrame()
    # Determine whether the variable exists in the table's columns.
    columns = table.columns.tolist()
    if (variable in columns):
        # Report.
        utility.report_table_column_categories_rows(
            column=variable,
            table=table,
        )
    else:
        print("The variable " + str(variable) + " is not available.")
    pass


##########
# Selections on samples for both sexes separate


def stratify_by_female_categories(
    female_age_grade=None,
    female_menopause=None,
    female_menstruation=None,
    female_pregnancy=None,
    table=None,
    report=None,
):
    """
    Control filters, selections, stratifications on samples (rows) within a
    table for persons of female sex.

    arguments:
        female_age_grade (list<int>): categorical values for selection
        female_menopause (list<int>): categorical values for selection
        female_menstruation (list<int>): categorical values for selection
        female_pregnancy (list<int>): categorical values for selection
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Select records for females.
    table = table.loc[
        (table["sex_text"] == "female"), :
    ]
    # Select samples (rows) within a table.
    selections = dict()
    selections["age_grade_female"] = female_age_grade
    selections["menopause_ordinal"] = female_menopause
    selections["menstruation_regular_range"] = female_menstruation
    selections["pregnancy"] = female_pregnancy
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
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "stratify_by_female_categories()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


def stratify_by_male_categories(
    male_age_grade=None,
    table=None,
    report=None,
):
    """
    Control filters, selections, stratifications on samples (rows) within a
    table for persons of male sex.

    arguments:
        male_age_grade (list<int>): categorical values for selection
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Select records for females.
    table = table.loc[
        (table["sex_text"] == "male"), :
    ]
    # Select samples (rows) within a table.
    selections = dict()
    selections["age_grade_male"] = male_age_grade
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
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "stratify_by_male_categories()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


def control_stratify_by_sex_separate_categories(
    sex_text=None,
    female_age_grade=None,
    female_menopause=None,
    female_menstruation=None,
    female_pregnancy=None,
    male_age_grade=None,
    table=None,
    report=None,
):
    """
    Control filters, selections, stratifications on samples (rows) within a
    table separately by sex.

    arguments:
        sex_text (list<str>): categorical values for selection
        female_age_grade (list<int>): categorical values for selection
        female_menopause (list<int>): categorical values for selection
        female_menstruation (list<int>): categorical values for selection
        female_pregnancy (list<int>): categorical values for selection
        male_age_grade (list<int>): categorical values for selection
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Collect records.
    table_collection = pandas.DataFrame()
    # Select records for females.
    if ("female" in sex_text):
        table_female = stratify_by_female_categories(
            female_age_grade=female_age_grade,
            female_menopause=female_menopause,
            female_menstruation=female_menstruation,
            female_pregnancy=female_pregnancy,
            table=table,
            report=report,
        )
        table_collection = table_collection.append(
            table_female,
            ignore_index=True,
        )
        pass
    # Select records for males.
    if ("male" in sex_text):
        table_male = stratify_by_male_categories(
            male_age_grade=male_age_grade,
            table=table,
            report=report,
        )
        table_collection = table_collection.append(
            table_male,
            ignore_index=True,
        )
        pass
    # Organize table.
    table_collection.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_collection.set_index(
        "eid",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "control_stratify_by_sex_separate_categories()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table_collection


##########
# Selections on samples by nonmissing values in specific variables


def stratify_by_nonmissing_values_specific_variables(
    names=None,
    prefixes=None,
    table=None,
    drop_columns=None,
    report=None,
):
    """
    Control filters, selections, stratifications on samples (rows) within a
    table by valid, non-missing values of specific variables.

    arguments:
        names (list<str>): explicit names of variables (columns) for which to
            require nonmissing values in all samples (rows)
        prefixes (list<str>): prefixes of names of variables (columns) for which
            to require nonmissing values in all samples (rows)
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        drop_columns (bool): whether to drop variables (columns) other than
            those for which to require nonmissing values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract table columns.
    columns_all = copy.deepcopy(table.columns.to_list())
    # Collect columns to keep.
    columns_keep = list()
    columns_names = list(filter(
        lambda column: (str(column) in names),
        columns_all
    ))
    columns_keep.extend(columns_names)
    for prefix in prefixes:
        columns_prefix = list(filter(
            lambda column: (str(prefix) in str(column)),
            columns_all
        ))
        columns_keep.extend(columns_prefix)
    # Select columns.
    if drop_columns:
        table = table.loc[
            :, table.columns.isin(columns_keep)
        ]
    # Remove any rows with null missing values in specific columns.
    table.dropna(
        axis="index",
        how="any",
        subset=columns_keep,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "stratify_by_nonmissing_values_specific_variables()"
        )
        print(name_function)
        print("Columns for which to require nonmissing values: ")
        print(columns_keep)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


##########
# Cohort, model selection: kinship


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


def filter_kinship_pairs_by_threshold_relevance(
    persons_relevance=None,
    threshold_kinship=None,
    table_kinship_pairs=None,
    report=None,
):
    """
    Filters kinship pairs of related persons by a threshold on the Kinship
    Coefficient and by genotype identifiers for persons who are relevant to a
    specific cohort table.

    This function returns a filtered table of kinship coefficients between pairs
    of genotype identifiers. After filtering, this table only includes pairs of
    genotypes that are:
    1. relevant to a specific cohort for genetic analysis
    2. related with a kinship coefficient beyond (greater than or equal to) a
       specific threshold

    Greater values of the Kinship Coefficient indicate a stronger genetic
    relationship between a pair of persons.

    Ignore kinship pairs of related persons if their Kinship Coefficient is
    below the threshold.

    This function returns pairs of related persons who are:
    1. relevant to a specific analysis cohort
    2. have Kinship Coefficient greater than or equal to (>=) threshold

    arguments:
        persons_relevance (list<str>): identifiers of persons who, on basis of
            phenotypes and other filters, are relevant to a genetic analysis
        threshold_kinship (float): maximal value of kinship coefficient
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of kinship coefficients across pairs of
            persons

    """

    # Copy information.
    persons_relevance = copy.deepcopy(persons_relevance)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)
    # Filter kinship pairs to exclude any redundant pairs.
    table_kinship_pairs.drop_duplicates(
        subset=["ID1", "ID2",],
        keep="first",
        inplace=True,
    )
    # Filter kinship pairs to exclude persons with kinship below threshold.
    # These pairs do not have strong enough relation for further consideration.
    table_kinship_pairs_below = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] < threshold_kinship), :
    ]
    table_kinship_pairs_above = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] >= threshold_kinship), :
    ]
    # Filter kinship pairs to exclude persons who are not in the
    # analysis-specific cohort table.
    # Consider both persons from each kinship pair.
    table_kinship_pairs_above.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_kinship_pairs_above.set_index(
        "ID1",
        append=False,
        drop=True,
        inplace=True
    )
    table_kinship_pairs_relevance = table_kinship_pairs_above.loc[
        table_kinship_pairs_above.index.isin(persons_relevance), :
    ]
    table_kinship_pairs_relevance.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_kinship_pairs_relevance.set_index(
        "ID2",
        append=False,
        drop=True,
        inplace=True
    )
    table_kinship_pairs_relevance = table_kinship_pairs_relevance.loc[
        table_kinship_pairs_relevance.index.isin(persons_relevance), :
    ]
    table_kinship_pairs_relevance.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Report.
    if report:
        # Organization information.
        count_pairs_original = copy.deepcopy(table_kinship_pairs.shape[0])
        count_pairs_kinship_below = copy.deepcopy(
            table_kinship_pairs_below.shape[0]
        )
        count_pairs_kinship_above = copy.deepcopy(
            table_kinship_pairs_above.shape[0]
        )
        count_pairs_relevance = copy.deepcopy(
            table_kinship_pairs_relevance.shape[0]
        )
        persons_kin = copy.deepcopy(
            table_kinship_pairs_relevance["ID1"].to_list()
        )
        persons_kin_second = copy.deepcopy(
            table_kinship_pairs_relevance["ID2"].to_list()
        )
        persons_kin.extend(persons_kin_second)
        persons_kin = list(map(str, persons_kin)) # string
        persons_kin = list(set(persons_kin)) # unique
        count_persons_kin_relevance = len(persons_kin)
        # Report.
        utility.print_terminal_partition(level=3)
        print(
            "report: " +
            "filter_kinship_pairs_by_threshold_relevance()"
        )
        utility.print_terminal_partition(level=5)
        print("... kinship pairs ...")
        print("original pairs: " + str(count_pairs_original))
        print("Kinship Coefficient threshold: " + str(threshold_kinship))
        print("pairs kinship < threshold: " + str(count_pairs_kinship_below))
        print("pairs kinship >= threshold: " + str(count_pairs_kinship_above))
        print(
            "pairs relevant to analysis cohort: " + str(count_pairs_relevance)
        )
        print("... persons ...")
        print(
            "persons relevant and with kinship beyond threshold: " +
            str(count_persons_kin_relevance)
        )
        utility.print_terminal_partition(level=5)
    # Return information.
    return table_kinship_pairs_relevance


def select_kinship_network_component_representatives(
    persons_priority=None,
    table_kinship_pairs=None,
    report=None,
):
    """
    Selects persons who are priority representatives from distinct connected
    components (families) in kinship network.

    arguments:
        persons_priority (list<str>): identifiers of persons to select
            preferentially in kinship filter
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of persons who are priority representatives
            from connected components of kinship network

    """

    # Copy information.
    persons_priority = copy.deepcopy(persons_priority)
    persons_priority_set = set(persons_priority)
    # Define an unweighted, undirected network graph between pairs of related
    # persons.
    columns_pairs = ["ID1", "ID2"]
    table_kinship_pairs = table_kinship_pairs.loc[
        :, table_kinship_pairs.columns.isin(columns_pairs)
    ]
    table_kinship_pairs = table_kinship_pairs[[*columns_pairs]]
    network = networkx.from_pandas_edgelist(
        table_kinship_pairs,
        "ID1",
        "ID2",
        edge_attr=None,
        create_using=networkx.Graph(),
    )
    # Extract a single node (genotype) at random from each connected component
    # from the network graph of kinship pairs.
    # representative = random.choice(component.nodes())
    representatives = list()
    for component in networkx.connected_components(network):
        network_component = network.subgraph(component).copy()
        component_nodes = copy.deepcopy(list(network_component.nodes()))
        component_nodes_set = set(component_nodes)
        priority_nodes = list(
            persons_priority_set.intersection(component_nodes_set)
        )
        if (len(priority_nodes) > 0):
            representative = random.choice(priority_nodes)
        else:
            representative = random.choice(component_nodes)
            pass
        representatives.append(representative)
        pass
    # Convert genotype identifiers to strings.
    persons_represent = list(map(str, representatives))
    # Report.
    if report:
        # Organization information.
        count_persons_represent = len(persons_represent)
        # Report.
        utility.print_terminal_partition(level=3)
        print(
            "report: " +
            "select_kinship_network_component_representatives()"
        )
        utility.print_terminal_partition(level=5)
        print("... kinship network graph ...")
        print("network order: " + str(network.order()))
        print("network size: " + str(network.size()))
        print(
            "count connected components: " +
            str(networkx.number_connected_components(network))
        )
        print(
            "count representative nodes: " + str(count_persons_represent)
        )
        pass
    # Return information.
    return persons_represent


def filter_cohort_relevance_persons_by_priority_kinship(
    persons_relevance=None,
    persons_priority=None,
    threshold_kinship=None,
    table_kinship_pairs=None,
    report=None,
):
    """
    Filters persons who are relevant to a cohort for genetic analysis on basis
    of their kinship relatedness.

    arguments:
        persons_relevance (list<str>): identifiers of persons who, on basis of
            phenotypes and other filters, are relevant to a genetic analysis
        persons_priority (list<str>): identifiers of persons to select
            preferentially in kinship filter
        threshold_kinship (float): maximal value of kinship coefficient
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of persons who are relevant to a cohort for
            genetic analysis and are also "unrelated" on basis of threshold on
            kinship

    """

    # Copy information.
    persons_relevance = copy.deepcopy(persons_relevance)
    # Filter kinship pairs.
    # After filter, kinship pairs only include persons who are relevant to
    # the specific cohort for genetic analysis and also have relatedness with
    # other relevant persons at a kinship coefficient greater than or equal to
    # a threshold.
    table_kinship_pairs = filter_kinship_pairs_by_threshold_relevance(
        persons_relevance=persons_relevance,
        threshold_kinship=threshold_kinship,
        table_kinship_pairs=table_kinship_pairs,
        report=report, # report procedure is quite slow
    )
    persons_kin = copy.deepcopy(table_kinship_pairs["ID1"].to_list())
    persons_kin_second = copy.deepcopy(table_kinship_pairs["ID2"].to_list())
    persons_kin.extend(persons_kin_second)
    persons_kin = list(map(str, persons_kin)) # string
    persons_kin = list(set(persons_kin)) # unique
    # Select priority representative persons from connected components
    # (families) in kinship network.
    persons_represent = select_kinship_network_component_representatives(
        persons_priority=persons_priority,
        table_kinship_pairs=table_kinship_pairs,
        report=report,
    )

    # 1. persons_relevance: complete collection of persons relevant to cohort
    #    and genetic analysis
    # 2. persons_kin: persons relevant to cohort who also are related beyond
    #    threshold
    # 3. persons_represent: persons who are single representatives from each
    #    family in "persons_kin"

    # Exclude all related persons from the cohort, with the exception of a
    # single representative from each family (connected component of the kinship
    # pair network).

    # Combination:
    # - remove all "persons_kin" from "persons_relevance" with the exception of
    #   "persons_represent"
    # 1. remove "persons_represent" from "persons_kin"
    # 2. remove "persons_kin" from "persons_relevance"

    # Computation note:
    # At this scale, computing a set difference is much more efficient than
    # applying a filter across elements in two lists.
    #persons_keep = list(filter(
    #    lambda value: (value not in persons_exclusion),
    #    persons_relevance
    #))

    # Remove "persons_represent" from "persons_kin".
    persons_exclusion = list(set(persons_kin) - set(persons_represent))
    # Remove "persons_exclusion" from "persons_relevance".
    persons_keep = list(set(persons_relevance) - set(persons_exclusion))

    # Report.
    if report:
        # Organization information.
        count_persons_relevance = len(persons_relevance)
        count_persons_priority = len(persons_priority)
        count_persons_kin = len(persons_kin)
        count_persons_represent = len(persons_represent)
        count_persons_keep = len(persons_keep)

        count_loss_kinship = (count_persons_kin - count_persons_represent)
        percentage_loss_kinship = round(
            ((count_loss_kinship / count_persons_relevance) * 100), 2
        )
        count_loss_check = (
            count_persons_relevance - count_persons_keep
        )
        percentage_loss_check = round(
            ((count_loss_check / count_persons_relevance) * 100), 2
        )

        # Report.
        utility.print_terminal_partition(level=3)
        print(
            "report: " +
            "filter_cohort_relevance_persons_by_priority_kinship()"
        )
        utility.print_terminal_partition(level=5)
        print("... counts of persons ...")
        print("count_persons_relevance: " + str(count_persons_relevance))
        print("count_persons_priority: " + str(count_persons_priority))
        print("count_persons_kin: " + str(count_persons_kin))
        print("count_persons_represent: " + str(count_persons_represent))
        print("count_persons_keep: " + str(count_persons_keep))
        utility.print_terminal_partition(level=5)
        print(
            "cohort loss to kinship: " + str(count_loss_kinship) + " of " +
            str(count_persons_relevance) + " (" + str(percentage_loss_kinship) +
            "%)"
        )
        print(
            "check cohort loss to kinship: " + str(count_loss_check) + " of " +
            str(count_persons_relevance) + " (" + str(percentage_loss_check) +
            "%)"
        )

        pass
    # Return information.
    return persons_keep


def test_persons_relevance_maximal_kinship(
    persons_relevance=None,
    table_kinship_pairs=None,
    report=None,
):
    """
    Filters kinship pairs of related persons by genotype identifiers for persons
    who are relevant to a specific cohort table. Returns the maximal value of
    kinship coefficient between any pairs of these persons.

    arguments:
        persons_relevance (list<str>): identifiers of persons who, on basis of
            phenotypes and other filters, are relevant to a genetic analysis
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        report (bool): whether to print reports

    raises:

    returns:
        (float): maximal value of Kinship Coefficient between pairs of persons
            who are relevant to analysis cohort

    """

    # Copy information.
    persons_relevance = copy.deepcopy(persons_relevance)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)
    # Filter kinship pairs to exclude any redundant pairs.
    table_kinship_pairs.drop_duplicates(
        subset=["ID1", "ID2",],
        keep="first",
        inplace=True,
    )
    # Filter kinship pairs to exclude persons who are not in the
    # analysis-specific cohort table.
    # Consider both persons from each kinship pair.
    table_kinship_pairs.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_kinship_pairs.set_index(
        "ID1",
        append=False,
        drop=True,
        inplace=True
    )
    table_kinship_pairs_relevance = table_kinship_pairs.loc[
        table_kinship_pairs.index.isin(persons_relevance), :
    ]
    table_kinship_pairs_relevance.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_kinship_pairs_relevance.set_index(
        "ID2",
        append=False,
        drop=True,
        inplace=True
    )
    table_kinship_pairs_relevance = table_kinship_pairs_relevance.loc[
        table_kinship_pairs_relevance.index.isin(persons_relevance), :
    ]
    table_kinship_pairs_relevance.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Extract values of kinship coefficient.
    array_kinship = copy.deepcopy(
        table_kinship_pairs_relevance["Kinship"].dropna().to_numpy()
    )
    maximum_kinship = numpy.nanmax(array_kinship)
    # Report.
    if report:
        # Organization information.
        count_pairs_relevance = copy.deepcopy(
            table_kinship_pairs_relevance.shape[0]
        )
        persons_kin = copy.deepcopy(
            table_kinship_pairs_relevance["ID1"].to_list()
        )
        persons_kin_second = copy.deepcopy(
            table_kinship_pairs_relevance["ID2"].to_list()
        )
        persons_kin.extend(persons_kin_second)
        persons_kin = list(map(str, persons_kin)) # string
        persons_kin = list(set(persons_kin)) # unique
        count_persons_kin_relevance = len(persons_kin)
        # Report.
        utility.print_terminal_partition(level=3)
        print(
            "report: " +
            "test_persons_relevance_maximal_kinship()"
        )
        utility.print_terminal_partition(level=5)
        print("... kinship pairs ...")
        print(
            "pairs relevant to analysis cohort: " + str(count_pairs_relevance)
        )
        print("... persons ...")
        print(
            "relevant persons with kinship: " +
            str(count_persons_kin_relevance)
        )
        utility.print_terminal_partition(level=5)
        print("... maximal value of kinship coefficient in cohort ...")
        print("maximal kinship: " + str(maximum_kinship))
    # Return information.
    return maximum_kinship


def filter_table_phenotype_persons_by_kinship(
    name=None,
    priority_values=None,
    priority_variable=None,
    threshold_kinship=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Selects records of persons who are unrelated.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        priority_values (list<int>): which values of "priority_variable" to use
            to select priority persons for kinship filter
        priority_variable (list<str>): name of column for variable to consider
            to select priority persons for kinship filter
        threshold_kinship (float): maximal value of kinship coefficient
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Report.
    if report:
        # Report.
        utility.print_terminal_partition(level=2)
        print("... Begin cohort filter for kinship ...")
        utility.print_terminal_partition(level=5)
        print("cohort: " + str(name))
        utility.print_terminal_partition(level=2)

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract from table identifiers of persons or genotypes.
    persons_relevance = copy.deepcopy(table["IID"].to_list())
    # Determine priority persons for kinship filter.
    if (
        (priority_variable is not None) and (len(priority_values) > 0)
    ):
        # Select identifiers of priority persons.
        table_priority = table.loc[
            (table[priority_variable].isin(priority_values)), :
        ]
        persons_priority = copy.deepcopy(table_priority["IID"].to_list())
    else:
        persons_priority = []
        pass
    # Select unrelated persons to keep in the cohort table.
    persons_keep = filter_cohort_relevance_persons_by_priority_kinship(
        persons_relevance=persons_relevance,
        persons_priority=persons_priority,
        threshold_kinship=threshold_kinship,
        table_kinship_pairs=table_kinship_pairs,
        report=report,
    )
    # Select table records for relevant, unrelated persons in cohort.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table.set_index(
        "IID",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )
    table = table.loc[
        table.index.isin(persons_keep), :
    ]
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Report.
    if report:
        # Organization information.
        count_persons_keep = copy.deepcopy(table.shape[0])
        # Report.
        utility.print_terminal_partition(level=3)
        print(
            "report: " +
            "filter_table_phenotype_persons_by_kinship()"
        )
        utility.print_terminal_partition(level=5)
        print("count_persons_keep: " + str(count_persons_keep))
        # Test kinship filter.
        persons_test = copy.deepcopy(table["IID"].to_list())
        kinship_maximum = test_persons_relevance_maximal_kinship(
            persons_relevance=persons_test,
            table_kinship_pairs=table_kinship_pairs,
            report=report,
        )
        pass
    # Return information.
    return table


def report_kinship_filter_priority_selection(
    name=None,
    priority_values=None,
    priority_variable=None,
    table_full=None,
    table_simple=None,
    table_priority=None,
    report=None,
):
    """
    Reports counts of records for priority persons in cohort tables after
    Kinship Filter for genetic analyses.

    Counts are specific to the "IID" identifier to represent records that match to
    genotypes.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        priority_values (list<int>): which values of "priority_variable" to use
            to select priority persons in kinship filter
        priority_variable (list<str>): name of column for variable to consider
            to select priority persons in kinship filter
        threshold_kinship (float): maximal value of kinship coefficient
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        table (object): Pandas data frame of values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    """

    # Report.
    if report:
        # Report.
        utility.print_terminal_partition(level=2)
        print(
            "report: " +
            "report_kinship_filter_priority_selection()"
        )
        utility.print_terminal_partition(level=3)
        print(name)
        pass
    # Copy information in tables.
    table_full = table_full.copy(deep=True)
    table_simple = table_simple.copy(deep=True)
    table_priority = table_priority.copy(deep=True)
    # Organize tables.
    records = list()
    record = {
        "table": table_full,
        "description": "full phenotype table before Kinship Filter",
    }
    records.append(record)
    record = {
        "table": table_simple,
        "description": "table after Kinship Filter without priority",
    }
    records.append(record)
    record = {
        "table": table_priority,
        "description": "table after Kinship Filter with priority",
    }
    records.append(record)
    # Collect information for tables.
    for record in records:
        # Remove any table rows for records with null values in genotype
        # identifier ("IID") or priority variable.
        table_valid = record["table"].dropna(
            axis="index",
            how="any",
            subset=["IID", priority_variable],
            inplace=False,
        )
        # Select table rows for records representing priority persons.
        table_preference = table_valid.loc[
            (table_valid[priority_variable].isin(priority_values)), :
        ]
        table_not_preference = table_valid.loc[
            (~table_valid[priority_variable].isin(priority_values)), :
        ]
        # Extract identifiers of genotypes for priority persons.
        genotypes_preference = copy.deepcopy(
            table_preference["IID"].to_list()
        )
        genotypes_preference = list(map(str, genotypes_preference)) # string
        genotypes_preference = list(set(genotypes_preference)) # unique
        genotypes_not_preference = copy.deepcopy(
            table_not_preference["IID"].to_list()
        )
        genotypes_not_preference = list(map(str, genotypes_not_preference)) # string
        genotypes_not_preference = list(set(genotypes_not_preference)) # unique

        # Count genotypes.
        count_genotypes_priority = len(genotypes_preference)
        count_genotypes_not_priority = len(genotypes_not_preference)
        # Report.
        if report:
            # Report.
            utility.print_terminal_partition(level=4)
            print(record["description"])
            utility.print_terminal_partition(level=5)
            print(
                "Count of priority genotypes in table: " +
                str(count_genotypes_priority)
            )
            print(
                "Count of not priority genotypes in table: " +
                str(count_genotypes_not_priority)
            )
            pass
        pass
    pass


def control_stratify_by_unrelated_kinship(
    name=None,
    kinship_unrelated=None,
    priority_female=None,
    priority_male=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Control filters, selections, stratifications on samples (rows) within a
    table separately by sex.

    arguments:
        name (str): name of stratification cohort table
        kinship_unrelated (int): logical binary whether to select persons with
            unrelated kinship
        priority_female (int): logical binary whether to give priority to
            persons with female sex in selection of unrelated kinship
        priority_male (int): logical binary whether to give priority to
            persons with male sex in selection of unrelated kinship
        table_kinship_pairs (object): Pandas data frame table of kinship
            coefficients across pairs of persons
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Stratify by unrelated kinship.
    if (kinship_unrelated == 1):
        if (
            (priority_female == 1) and
            (priority_male == 0)
        ):
            priority_values = ["female", ]
            priority_variable = "sex_text"
        elif (
            (priority_male == 1) and
            (priority_female == 0)
        ):
            priority_values = ["male", ]
            priority_variable = "sex_text"
        else:
            priority_values = []
            priority_variable = None
            pass
        # Apply filters on kinship.
        table = filter_table_phenotype_persons_by_kinship(
            name=name,
            priority_values=priority_values,
            priority_variable=priority_variable,
            threshold_kinship=0.1, # pairs with kinship >= threshold for exclusion
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
        )
        pass
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.set_index(
        "eid",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "control_stratify_by_unrelated_kinship()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


def control_report_unrelated_kinship_priority(
    records_cohorts=None,
    table_full=None,
):
    """
    Control report on the prioritization of specific groups of people in the
    otherwise random selection of persons with unrelated kinship.

    arguments:
        records_cohorts (list<dict>): records of information about cohorts
        table_full (object): Pandas data frame table of variables (features)
            across columns and samples (records) across rows

    raises:

    returns:

    """

    # Organize hashable entries for cohorts.
    entries_cohorts = utility.structure_from_records_to_entries(
        records=records_cohorts,
    )
    # Extract tables for specific cohorts.
    name = "white_british_unrelated_female_male"
    if (name in entries_cohorts.keys()):
        table_simple = entries_cohorts[name]["table"]
    else:
        table_simple = pandas.DataFrame()
    name = "white_british_unrelated_female_male_kinship_priority_male"
    if (name in entries_cohorts.keys()):
        table_priority = entries_cohorts[name]["table"]
    else:
        table_priority = pandas.DataFrame()

    # Report.
    report_kinship_filter_priority_selection(
        name="... Comparison of priority to male persons in cohort ...",
        priority_values=["male",],
        priority_variable="sex_text",
        table_full=table_full,
        table_simple=table_simple,
        table_priority=table_priority,
        report=True,
    )
    pass


##########
# Selection of columns for relevant variables


def select_table_columns_relevant_variables(
    dependence=None,
    variables_nonmissing=None,
    variables_prefix_nonmissing=None,
    independence=None,
    independence_prefix=None,
    table=None,
    report=None,
):
    """
    Control selection on relevan variables (columns) within a table.

    arguments:
        dependence (list<str>): names of table's columns for dependent variables
        variables_nonmissing (list<str>): names of table's columns for variables
            for which to select non-missing values
        variables_prefix_nonmissing (list<str>): prefixes of names of table's
            columns for variables for which to select non-missing values
        independence (list<str>): names of table's columns for independent
            variables
        independence_prefix (list<str>): prefixes of names of table's columns
            for independent variables
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract table columns.
    columns_all = copy.deepcopy(table.columns.to_list())
    # Collect columns to keep.
    columns_keep = list()

    # Columns with specific names.
    columns_names = list(filter(
        lambda column: (str(column) in variables_nonmissing),
        columns_all
    ))
    columns_keep.extend(columns_names)
    columns_names = list(filter(
        lambda column: (str(column) in dependence),
        columns_all
    ))
    columns_keep.extend(columns_names)
    columns_names = list(filter(
        lambda column: (str(column) in independence),
        columns_all
    ))
    columns_keep.extend(columns_names)

    # Columns with prefixes.
    for prefix in variables_prefix_nonmissing:
        columns_prefix = list(filter(
            lambda column: (str(prefix) in str(column)),
            columns_all
        ))
        columns_keep.extend(columns_prefix)
    for prefix in independence_prefix:
        columns_prefix = list(filter(
            lambda column: (str(prefix) in str(column)),
            columns_all
        ))
        columns_keep.extend(columns_prefix)

    # Select columns for relevant variables.
    table = table.loc[
        :, table.columns.isin(columns_keep)
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "select_table_columns_relevant_variables()"
        )
        print(name_function)
        print("Columns to keep: ")
        print(columns_keep)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


##########
# Transformations on the scales of variables


def control_transformations_on_variables(
    dependence=None,
    dependence_logical_binary=None,
    dependence_scale_transforms=None,
    independence=None,
    independence_prefix=None,
    independence_scale=None,
    table=None,
    report=None,
):
    """
    Control selection on relevan variables (columns) within a table.

    arguments:
        dependence (list<str>): names of table's columns for dependent variables
        dependence_logical_binary (int): logical binary whether to translate
            values of logical binary dependent variable to code for PLINK2 with
            zero mapping to one (0 --> 1; control) and one mapping to two
            (1 --> 2; case)
        dependence_scale_transforms (int): logical binary whether to apply
            distribution scale transformations to values of dependent variables
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
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract table columns.
    columns_all = copy.deepcopy(table.columns.to_list())

    # Confirm names of columns for dependent and independent variables.
    columns_dependence = list(filter(
        lambda column: (str(column) in dependence),
        columns_all
    ))
    columns_independence = list(filter(
        lambda column: (str(column) in independence),
        columns_all
    ))
    for prefix in independence_prefix:
        columns_prefix = list(filter(
            lambda column: (str(prefix) in str(column)),
            columns_all
        ))
        columns_independence.extend(columns_prefix)

    # Confirm that all continuous-scale dependent and independent variables are
    # in a float numerical type.
    # Dependent variables or independent covariates that need conversion to
    # float numerical type include measurement scale binary discrete, ordinal
    # discrete, interval continuous, or ratio continuous.
    if (
        (not pandas.isna(dependence_logical_binary)) and
        (dependence_logical_binary == 0)
    ):
        table = utility.convert_table_columns_variables_types_float(
            columns=columns_dependence,
            table=table,
        )
        pass
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_independence,
        table=table,
    )

    # Apply transformations on the scales of dependent variables.
    if (
        (not pandas.isna(dependence_scale_transforms)) and
        (dependence_scale_transforms == 1)
    ):
        for variable in columns_dependence:
            table = pscale.apply_transformations_to_variable_distribution_scale(
                column=str(variable),
                logarithm_e=True,
                standard_z_score=True,
                rank_inverse=True,
                suffix_logarithm_e="_log",
                suffix_standard_z_score="_z",
                suffix_rank_inverse="_rank",
                table=table,
                report=False,
            )
        pass
    # Apply Standard Z-Score transformation on independent variables.
    if (
        (not pandas.isna(independence_scale)) and
        (independence_scale == 1)
    ):
        table = pscale.drive_transform_variables_distribution_scale_z_score(
            table=table,
            columns=columns_independence,
            report=False,
        )
        pass

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "control_transformations_on_variables()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


##########
# Adjustment of table format for interpretation in PLINK2


def translate_boolean_phenotype_plink2(
    boolean_value=None,
):
    """
    Translate information from simple binary representation to plink
    representation.

    Accommodate inexact float values and null values.

    https://www.cog-genomics.org/plink/1.9/input
    Section: "Phenotypes"
    Subsection: "Phenotype encoding"
        "Missing phenotypes are normally expected to be encoded as -9. You can
        change this to another integer with --missing-phenotype. (This is a
        slight change from PLINK 1.07: floating point values are now disallowed
        due to rounding issues, and nonnumeric values such as 'NA' are rejected
        since they're treated as missing phenotypes no matter what. Note that
        --output-missing-phenotype can be given a nonnumeric string.)

        Case/control phenotypes are expected to be encoded as 1=unaffected
        (control), 2=affected (case); 0 is accepted as an alternate missing
        value encoding. If you use the --1 flag, 0 is interpreted as unaffected
        status instead, while 1 maps to affected. This also forces phenotypes to
        be interpreted as case/control."

    arguments:
        boolean_value (float): boolean (False, True) representation of a
            phenotype

    raises:

    returns:
        (float): plink binary representation of a phenotype

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (not pandas.isna(boolean_value)):
        # The variable has a valid value.
        if (not boolean_value):
            # control, original value: False
            # control, Plink translation: 1
            value = 1
        elif (boolean_value):
            # case, original value: True
            # case, Plink translation: 2
            value = 2
    else:
        # null
        #value = -9
        value = "nan"
    # Return information.
    return value


# review: TCW; 9 November 2022
def translate_binary_phenotype_plink2(
    binary_value=None,
):
    """
    Translate a dependent variable from logical binary representation to format
    for interpretation in PLINK2.

    Original code:
    0: control
    1: case

    PLINK2 code:
    1: control (unaffected)
    2: case (affected)

    https://www.cog-genomics.org/plink/2.0/input#phenotypes
    Section: "Phenotypes"
    Subsection: "Phenotype encoding"
        "
        Missing case/control or quantitative phenotypes are expected to be
        encoded as 'NA'/'nan' (any capitalization) or -9. (Other strings which
        don't start with a number are now interpreted as categorical
        phenotype/covariate values.) You can change the numeric missing
        phenotype code to another integer with --input-missing-phenotype, or
        just disable -9 with --no-input-missing-phenotype.

        Case/control phenotypes are expected to be encoded as 1=unaffected
        (control), 2=affected (case); 0 is accepted as an alternate missing
        value encoding. If you use the --1 flag, 0 is interpreted as unaffected
        status instead, while 1 maps to affected. Note that this only affects
        interpretation of input files; output files still use 1=control/2=case
        encoding.
        (Unlike PLINK 1.x, this does not force all phenotypes to be interpreted
        as case/control.)
        "

    arguments:
        binary_value (float): logical binary (0, 1) representation of a
            case-control phenotype

    raises:

    returns:
        (float): PLINK2 representation of a case-control phenotype

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(binary_value)) and
        (-0.5 <= binary_value and binary_value < 1.5)
    ):
        # The variable has a valid value.
        if (binary_value == 0):
            # control, original value: 0
            # control, PLINK2 translation: 1
            value = 1
        elif (binary_value == 1):
            # case, original value: 1
            # case, PLINK2 translation: 2
            value = 2
        else:
            # missing value
            #value = -9
            value = "nan"
    else:
        # missing value
        #value = -9
        value = "nan"
    # Return information.
    return value


def translate_binary_biological_sex_plink2(
    sex_y=None,
):
    """
    Translate biological sex to format for PLINK2.
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
            # PLINK2: 2, "female"
            interpretation = 2
        elif (sex_y == 1):
            # sex_y: 1, "male"
            # PLINK2: 1, "male"
            interpretation = 1
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return information.
    return interpretation


# review: TCW; 9 November 2022
def translate_text_biological_sex_plink2(
    sex_text=None,
):
    """
    Translate biological sex to format for PLINK2.

    PLINK2 code:
    2: female (XX)
    1: male (XY)
    NA: missing or unknown

    https://www.cog-genomics.org/plink/2.0/formats#cov
    Section: ".cov (covariate table)"

    arguments:
        sex_text (str): textual representation of biological sex, female or male

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(sex_text)) and
        ((sex_text == "female") or (sex_text == "male"))
    ):
        # The variable has a valid value.
        # Interpret the value.
        if (sex_text == "female"):
            # sex_text: "female"
            # PLINK2: 2, "female"
            interpretation = 2
        elif (sex_text == "male"):
            # sex_text: "male"
            # PLINK2: 1, "male"
            interpretation = 1
        else:
            # Uninterpretable value
            interpretation = float("nan")
    else:
        # Missing or uninterpretable value
        interpretation = float("nan")
    # Return information.
    return interpretation


# review: TCW; 9 November 2022
def translate_format_phenotype_covariate_table_plink2(
    dependence_boolean=None,
    dependence_binary=None,
    sex_text=None,
    remove_missing_rows=None,
    table=None,
):
    """
    Translate the format of a table for phenotype and covariate variables for
    interpretation in PLINK2.

    1. Remove any rows with missing, empty values.
    PLINK2 cannot accommodate rows with empty cells.
    This function removes table's rows with any missing values in any columns.
    It might be advantageous to apply a separate function to remove missing
    values only in the minimal essential columns, such as those for dependent
    variables and independent covariates.

    2. Introduce family identifiers.
    Family (FID) and individual (IID) identifiers must match the ID_1 and ID_2
    columns in the sample table.

    3. If "SID", "PAT", and "MAT" are irrelevant, then omit these columns.
    PLINK2 throws an error if these columns exist with empty values.

    4. Sort column sequence.
    PLINK2 requires FID and IID columns to come first.

    https://www.cog-genomics.org/plink/2.0/formats#cov
    Section: ".cov (covariate table)"
        "
        .cov (covariate table)
        Produced by --write-covar, --make-[b]pgen/--make-bed, and --export when
        covariates have been loaded/specified. Valid input for --covar.

        A text file with a header line, and one line per sample with the
            following columns:

        Header Column set Contents
        FID maybefid, fid Family ID
        IID (required) Individual ID
        SID maybesid, sid Source ID
        PAT maybeparents, parents Paternal individual ID
        MAT maybeparents, parents Maternal individual ID
        SEX sex Sex (1 = male, 2 = female, 'NA' = unknown)
        PHENO1 pheno1 All-missing phenotype column, if none loaded
        <Pheno name>, ... pheno1, phenos Phenotype value(s) (only first if just
            'pheno1')
        <Covar name>, ... (required) Covariate values

        (Note that --covar can also be used with files lacking a header row.)
    "

    arguments:
        dependence_boolean (list<str>): names of columns with discrete boolean
            (False: control or negative; True: case or positive) dependent
            variables that need conversion to PLINK format (1: control; 2: case)
        dependence_binary (list<str>): names of columns with discrete logical
            binary (0: control or negative; 1: case or positive) dependent
            variables that need conversion to PLINK format (control: 1, case: 2)
        sex_text (str): name of column for textual representation of biological
            sex ("female" or "male") for clarity
        remove_missing_rows (bool): whether to remove table's rows with any
            missing values in any columns
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows in format for GWAS in PLINK2

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Translate binary phenotype variables.
    for boolean_variable in dependence_boolean:
        name_plink = str(boolean_variable + "_plink")
        table[name_plink] = table.apply(
            lambda row:
                translate_boolean_phenotype_plink2(
                    boolean_value=row[boolean_variable],
                ),
            axis="columns", # apply across rows
        )
        pass
    # Translate binary phenotype variables.
    for binary_variable in dependence_binary:
        name_plink = str(binary_variable + "_plink")
        table[name_plink] = table.apply(
            lambda row:
                translate_binary_phenotype_plink2(
                    binary_value=row[binary_variable],
                ),
            axis="columns", # apply across rows
        )
        pass

    # Sex.
    # 2 : female (XX)
    # 1 : male (XY)
    # "NA": missing or unknown
    table["SEX"] = table.apply(
        lambda row:
            translate_text_biological_sex_plink2(
                sex_text=row[sex_text],
            ),
        axis="columns", # apply function to each row
    )
    # Organize.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Remove table's rows with any missing values in any columns.
    if remove_missing_rows:
        table.dropna(
            axis="index",
            how="any",
            subset=None,
            inplace=True,
        )
        pass

    # Introduce family identifier.
    table["FID"] = table["IID"]
    # Introduce other PLINK2 columns.
    #table["SID"] = ""
    #table["PAT"] = ""
    #table["MAT"] = ""
    columns_special = [
        "eid", "SEX", "MAT", "PAT", "SID", "IID", "FID",
    ]
    # Sort column sequence.
    columns = copy.deepcopy(table.columns.to_list())
    columns_sequence = list(filter(
        lambda element: element not in columns_special,
        columns
    ))
    columns_sequence.insert(0, "eid") # column 7
    columns_sequence.insert(0, "SEX") # column 6
    #columns_sequence.insert(0, "MAT") # column 5
    #columns_sequence.insert(0, "PAT") # column 4
    #columns_sequence.insert(0, "SID") # column 3
    columns_sequence.insert(0, "IID") # column 2
    columns_sequence.insert(0, "FID") # column 1
    table = table.loc[
        :, table.columns.isin(columns_sequence)
    ]
    table = table[[*columns_sequence]]
    # Return information.
    return table


def control_translate_table_format_for_plink2(
    dependence=None,
    dependence_logical_binary=None,
    table=None,
    report=None,
):
    """
    Control format translation of phenotype table after cohort stratification
    for interpretation in PLINK2.

    arguments:
        dependence (list<str>): names of table's columns for dependent variables
        dependence_logical_binary (int): logical binary whether to translate
            values of logical binary dependent variable to code for PLINK2 with
            zero mapping to one (0 --> 1; control) and one mapping to two
            (1 --> 2; case)
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of variables (features) across columns
            and samples (records) across rows

    """

    # Extract table columns.
    columns_all = copy.deepcopy(table.columns.to_list())
    # Confirm names of columns for dependent variables.
    columns_dependence = list(filter(
        lambda column: (str(column) in dependence),
        columns_all
    ))
    # Translate table format for interpretation in PLINK2.
    if (
        (not pandas.isna(dependence_logical_binary)) and
        (dependence_logical_binary == 1)
    ):
        table = translate_format_phenotype_covariate_table_plink2(
            dependence_boolean=[],
            dependence_binary=columns_dependence,
            sex_text="sex_text",
            remove_missing_rows=False, # remove NAs for specific columns only
            table=table,
        )
    else:
        table = translate_format_phenotype_covariate_table_plink2(
            dependence_boolean=[],
            dependence_binary=[],
            sex_text="sex_text",
            remove_missing_rows=False, # remove NAs for specific columns only
            table=table,
        )
        pass

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "control_translate_table_format_for_plink2()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


##########
# Control


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
    priority_male=None,
    dependence=None,
    dependence_scale_transforms=None,
    dependence_logical_binary=None,
    variables_nonmissing=None,
    variables_prefix_nonmissing=None,
    independence=None,
    independence_prefix=None,
    independence_scale=None,
    table_kinship_pairs=None,
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
        female_menstruation (list<int>): categorical values for selection
        female_pregnancy (list<int>): categorical values for selection
        male_age_grade (list<int>): categorical values for selection
        kinship_unrelated (int): logical binary whether to select persons with
            unrelated kinship
        priority_female (int): logical binary whether to give priority to
            persons with female sex in selection of unrelated kinship
        priority_male (int): logical binary whether to give priority to
            persons with male sex in selection of unrelated kinship
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
        table_kinship_pairs (object): Pandas data frame table of kinship
            coefficients across pairs of persons
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["group"] = str(group).strip()
    record["name"] = str(name).strip()

    # Control filters, selections, stratifications on all samples (rows)
    # regardless of sex.
    table = stratify_by_sex_together_categories(
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
        female_pregnancy=female_pregnancy,
        male_age_grade=male_age_grade,
        table=table,
        report=report,
    )

    # Control filters, selections, stratifications on samples (rows) by the
    # non-missingness of values of variables (columns).
    table = stratify_by_nonmissing_values_specific_variables(
        names=variables_nonmissing,
        prefixes=variables_prefix_nonmissing,
        table=table,
        drop_columns=False,
        report=report,
    )

    # Control filters, selections, stratifications on samples (rows) by the
    # kinship coefficient (relatedness) between persons in the cohort.
    if (
        (str("female_male_bipolar_case_body") in str(name)) or
        (str("female_male_bipolar_control_body") in str(name))
    ):
        report_kinship = True
    else:
        report_kinship = report
    table = control_stratify_by_unrelated_kinship(
        name=name,
        kinship_unrelated=kinship_unrelated,
        priority_female=priority_female,
        priority_male=priority_male,
        table_kinship_pairs=table_kinship_pairs,
        table=table,
        report=report_kinship,
    )

    # Control filters, selections, stratifications on variables (columns) in the
    # final table.
    table = select_table_columns_relevant_variables(
        dependence=dependence,
        variables_nonmissing=variables_nonmissing,
        variables_prefix_nonmissing=variables_prefix_nonmissing,
        independence=independence,
        independence_prefix=independence_prefix,
        table=table,
        report=report,
    )

    # Apply any transformations to the scale or distribution of values in
    # dependent and independent variables.
    table = control_transformations_on_variables(
        dependence=dependence,
        dependence_logical_binary=dependence_logical_binary,
        dependence_scale_transforms=dependence_scale_transforms,
        independence=independence,
        independence_prefix=independence_prefix,
        independence_scale=independence_scale,
        table=table,
        report=report,
    )

    # Translate format of table for interpretation in PLINK2.
    table = control_translate_table_format_for_plink2(
        dependence=dependence,
        dependence_logical_binary=dependence_logical_binary,
        table=table,
        report=report,
    )

    # Report count of samples in tables for special cohorts.
    if (
        (str("female_male_bipolar_case_body") in str(name)) or
        (str("female_male_bipolar_control_body") in str(name))
    ):
        utility.print_terminal_partition(level=2)
        print("Sample Counts for Special Cohort")
        print("group: " + str(group))
        print("cohort: " + str(name))
        print("sample counts:")
        print(str(table.shape[0]))
        utility.print_terminal_partition(level=3)
        pass

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


def write_genotype_product_table(
    name=None,
    table=None,
    path_directory=None,
):
    """
    Writes product information in a Pandas data frame table to file.

    arguments:
        name (str): base name for file
        table (object): Pandas data frame table
        path_directory (str): path to directory in which to write file

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_directory, str(name + ".pickle")
    )
    path_table_text = os.path.join(
        path_directory, str(name + ".tsv")
    )
    # Write information to file.
    table.to_pickle(
        path_table
    )
    table.to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=False,
        na_rep="NA",
    )
    pass


def control_write_genotype_product(
    pail_write=None,
    directory=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        directory (str): name of child directory to create within parent
            directory
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Define path to directories.
    path_directory = os.path.join(path_parent, directory)
    # Initialize directories.
    utility.create_directories(
        path=path_directory,
    )
    # Write genotype cohort, model tables in PLINK format.
    write_genotype_product_table(
        name=pail_write["name"],
        table=pail_write["table"],
        path_directory=path_directory,
    )
    pass


###############################################################################
# Procedure

# TODO: TCW; 8 November 2022
# TODO: may I simply drop indicator variables for all "site" categories with fewer than 10 samples?

#"genotype_array_axiom"

#Assessment Site
#- Example: "site_indicator_manchester"
#- Example: "site_indicator_oxford"
#- Include all of these site indicator variables in the phenotype table for GWAS.
#- Also include the "site" variable for convenience.
#- Summarize counts for each site.
#- - For each value of "site" variable, report counts of records.
#- Select independent covariates accordingly.


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
    print("version check: TCW; 9 November 2022")
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
            priority_male=pail["priority_male"],
            dependence=pail["dependence"],
            dependence_scale_transforms=pail["dependence_scale_transforms"],
            dependence_logical_binary=pail["dependence_logical_binary"],
            variables_nonmissing=pail["variables_nonmissing"],
            variables_prefix_nonmissing=pail["variables_prefix_nonmissing"],
            independence=pail["independence"],
            independence_prefix=pail["independence_prefix"],
            independence_scale=pail["independence_scale"],
            table_kinship_pairs=table_kinship_pairs,
            table=source["table_phenotypes"],
            report=False,
        )
        # Collect records.
        records_cohorts.append(record_cohort)
        pass
    # Report on priority kinship selection.
    control_report_unrelated_kinship_priority(
        records_cohorts=records_cohorts,
        table_full=source["table_phenotypes"],
    )
    # Report on categorical values within a cohort of interest.
    control_report_category_values_in_cohort(
        cohort="white_british_unrelated_female_male_bipolar_control_body",
        variable="site",
        records_cohorts=records_cohorts,
    )
    control_report_category_values_in_cohort(
        cohort="white_british_unrelated_female_male_bipolar_case_body",
        variable="site",
        records_cohorts=records_cohorts,
    )

    # Write product information to file.
    for record_cohort in records_cohorts:
        control_write_genotype_product(
            pail_write=record_cohort,
            directory=record_cohort["group"],
            path_parent=paths["genotype"],
        )
        pass
    pass



#
