"""
Organize information and define variables from data of the U.K. Biobank.

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
    utility.create_directories(
        path=paths["cohorts_models"]
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
    # Organize code interpretations.
    table_field_54 = pandas.read_csv(
        path_table_field_54,
        sep="\t",
        header=0,
        dtype={
            "code": "string",
            "interpretation": "string",
            "coding": "string",
            "meaning": "string",
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
    # Compile and return information.
    return {
        "field_54": field_54_codes_interpretations,
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
    "Caucasian": 1

    Accommodate inexact float values.

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
        axis="columns", # apply across rows
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


def interpret_assessment_center(
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


def interpret_assessment_month(
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


def interpret_sex_consensus(
    field_31=None,
    field_22001=None,
):
    """
    Determine consensus sex (biological sex, not social gender).
    Prioritize interpretation of the genetic sex variable.

    Data-Field "31": "sex"
    UK Biobank data coding "9" for variable field "31".
    "female": 0
    "male": 1

    Data-Field "22001": "genetic sex"
    UK Biobank data coding "9" for variable field "22001".
    "female": 0
    "male": 1

    arguments:
        field_31 (float): UK Biobank field 31, person's self-reported sex
        field_22001 (float): UK Biobank field 22001, ...

    raises:

    returns:
        (float): interpretation value

    """

    if (
        (
            (not pandas.isna(field_22001)) and
            (-0.5 <= field_22001 and field_22001 < 1.5)
        )
    ):
        # Genetic sex variable has a valid value.
        # Prioritize interpretation of the genetic sex variable.
        if (-0.5 <= field_22001 and field_22001 < 0.5):
            # "female": 0
            value = 0
        elif (0.5 <= field_22001 and field_22001 < 1.5):
            # "male": 1
            value = 1
    elif (
        (
            (not pandas.isna(field_31)) and
            (-0.5 <= field_31 and field_31 < 1.5)
        )
    ):
        # Self-reported sex variable has a valid value.
        if (-0.5 <= field_31 and field_31 < 0.5):
            # "female": 0
            value = 0
        elif (0.5 <= field_31 and field_31 < 1.5):
            # "male": 1
            value = 1
    else:
        # Sex is missing in both variables.
        value = float("nan")
    # Return information.
    return value


def determine_assessment_center_category(
    field_54=None,
    codes_interpretations=None,
):
    """
    Determine whether female persons experienced bilateral oophorectomy.

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
    value = interpret_assessment_center(
        field_54=field_54,
        codes_interpretations=codes_interpretations,
    )
    # Return information.
    return value


def determine_sex_text(
    sex=None,
):
    """
    Translate binary representation of sex to textual representation.

    arguments:
        sex (float): binary representation of person's sex

    raises:

    returns:
        (str): text representation of person's sex

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(sex)) and
        (
            (sex == 0) or
            (sex == 1)
        )
    ):
        # The variable has a valid value.
        if (sex == 0):
            # "female"
            sex_text = "female"
        elif (sex == 1):
            # "male"
            sex_text = "male"
    else:
        # null
        sex_text = "nan"
    # Return information.
    return sex_text


def determine_assessment_month_order(
    field_55=None,
):
    """
    Determine whether female persons experienced bilateral oophorectomy.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret oophorectomy.
    month = interpret_assessment_month(
        field_55=field_55,
    )
    # Set value.
    value = month
    # Return information.
    return value


def determine_assessment_month_category(
    field_55=None,
):
    """
    Determine whether female persons experienced bilateral oophorectomy.

    arguments:
        field_55 (float): UK Biobank field 55, month of assessment

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret oophorectomy.
    month = interpret_assessment_month(
        field_55=field_55,
    )
    # Set value.
    if (
        (not pandas.isna(month))
    ):
        # The variable has a valid value for month.
        if (month == 1):
            value = "January"
        elif (month == 2):
            value = "February"
        elif (month == 3):
            value = "March"
        elif (month == 4):
            value = "April"
        elif (month == 5):
            value = "May"
        elif (month == 6):
            value = "June"
        elif (month == 7):
            value = "July"
        elif (month == 8):
            value = "August"
        elif (month == 9):
            value = "September"
        elif (month == 10):
            value = "October"
        elif (month == 11):
            value = "November"
        elif (month == 12):
            value = "December"
        else:
            # Uninterpretable value.
            value = "nan"
    else:
        # Missing value.
        value = "nan"
    # Return information.
    return value


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
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("define_ordinal_stratifications_by_sex_continuous_variables")
            print("variable: " + str(variable))
            print("female tertiles")
            print(pandas.qcut(
                table_female[variable],
                q=3,
                labels=[0, 1, 2,],
                retbins=True,
            ))
            print(pandas.qcut(
                table_female[variable],
                q=3,
                labels=[0, 1, 2,],
            ).value_counts())
            print("male tertiles")
            print(pandas.qcut(
                table_male[variable],
                q=3,
                labels=[0, 1, 2,],
                retbins=True,
            ))
            print(pandas.qcut(
                table_male[variable],
                q=3,
                labels=[0, 1, 2,],
            ).value_counts())
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
        print(table_collection)
    # Return information.
    return table_collection


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
        column (str): name of column with categorical variable for which to
            define binary dummies and reduce by principal components
        prefix (str): prefix for names of new dummy and principal component
            columns in table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    def match_column_dummy(
        name=None,
        prefix=None,
        separator=None,
    ):
        if (separator in str(name)):
            name_prefix = str(name).split(separator)[0].strip()
            if (str(prefix) == str(name_prefix)):
                match = True
            else:
                match = False
        else:
            match = False
        return match

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
    table_indicators = table_indicators.loc[
        :, table_indicators.columns.isin([index, column])
    ]
    # Create binary dummies for variable categories.
    table_indicators = pandas.get_dummies(
        table_indicators,
        prefix=prefix,
        prefix_sep=separator,
        dummy_na=True,
        columns=[column],
        drop_first=False, # whether to create "k - 1" dummies, adequate
        dtype=numpy.uint8,
    )
    table_indicators.set_index(
        index,
        append=False,
        drop=True,
        inplace=True
    )
    table = table.merge(
        table_indicators,
        how="outer",
        left_on=index,
        right_on=index,
        suffixes=("_original", "_indicator"),
    )

    # Extract names of columns for dummies.
    columns = copy.deepcopy(table.columns.to_list())
    columns_dummies = list(filter(
        lambda column_trial: match_column_dummy(
            name=column_trial,
            prefix=prefix,
            separator=separator,
        ),
        columns
    ))

    columns_example = copy.deepcopy(columns_dummies)
    columns_example.extend(["eid", "IID", column])

    table_demo = table.copy(deep=True)
    table_demo = table_demo.loc[
        :, table_demo.columns.isin(columns_example)
    ]
    print("here is the table after dummy stuff...")
    print(table_demo)

    # Report.
    print("here are the columns in the table...")
    print(table.columns.to_list())
    unique_values = table[column].unique()
    count_unique_values = unique_values.size
    count_dummies = len(columns_dummies)
    utility.print_terminal_partition(level=2)
    print("report: create_reduce_categorical_variable_indicators()")
    utility.print_terminal_partition(level=3)
    print("count of original unique values: " + str(count_unique_values))
    print(unique_values)
    print("count of dummy variables: " + str(count_dummies))

    # Return information.
    return table


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

    # Copy information.
    table = table.copy(deep=True)

    # Read reference table of interpretations of variable codes.
    source = read_source_fields_codes_interpretations(
        path_dock=path_dock,
    )

    # Determine assessment center.
    table["assessment_site"] = table.apply(
        lambda row:
            determine_assessment_center_category(
                field_54=row["54-0.0"],
                codes_interpretations=source["field_54"],
            ),
        axis="columns", # apply across rows
    )
    # Determine month of assessment.
    table["month_order"] = table.apply(
        lambda row:
            determine_assessment_month_order(
                field_55=row["55-0.0"],
            ),
        axis="columns", # apply across rows
    )
    table["month"] = table.apply(
        lambda row:
            determine_assessment_month_category(
                field_55=row["55-0.0"],
            ),
        axis="columns", # apply across rows
    )

    # Create binary indicators of categories and reduce their dimensionality.

    # TODO: TCW 10 August 2021
    # TODO: I need to introduce a function to create categorical dummy variables
    # TODO: AND to prepare PCA reductions of those dummy variables
    # TODO: AND report variance for each PC

    table = create_reduce_categorical_variable_indicators(
        table=table,
        index="eid",
        column="assessment_site",
        prefix="site",
        separator="_",
        report=True,
    )


    # TODO: TCW 10 August 2021
    # TODO: follow pattern of interpret --> determine

    # Convert variable types.
    columns_type = [
        "31-0.0", "22001-0.0", "21022-0.0", "21001-0.0"
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )

    # Translate column names.
    translations = dict()
    translations["21022-0.0"] = "age"
    translations["21001-0.0"] = "body_mass_index"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Determine sex consensus between self-report and genotypic sex.
    table["sex"] = table.apply(
        lambda row:
            interpret_sex_consensus(
                field_31=row["31-0.0"],
                field_22001=row["22001-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine text representation of person's sex.
    table["sex_text"] = table.apply(
        lambda row:
            determine_sex_text(
                sex=row["sex"],
            ),
        axis="columns", # apply across rows
    )

    # Determine age.


    # Determine body mass index (BMI).


    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=["body_mass_index"],
        table=table,
    )
    # Introduce stratification variables by values of continuous variables, in
    # particular age within female and male persons.
    table = define_ordinal_stratifications_by_sex_continuous_variables(
        index="eid",
        sex_text="sex_text",
        variables=["age"],
        table=table,
        report=report,
    )
    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "31-0.0",
            "55-0.0",
            "22001-0.0", #"21022-0.0",
            "21002-0.0", "50-0.0",
            #"21001-0.0",
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
        "sex", "sex_text",
        "age", "age_grade_female", "age_grade_male",
        "body_mass_index", "body_mass_index_log",
        "month",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: organize_sex_age_body_variables()")
        utility.print_terminal_partition(level=3)
        print("translations of general attribute column names...")
        for old in translations.keys():
            print("   " + old + ": " + translations[old])
        utility.print_terminal_partition(level=3)
        utility.print_terminal_partition(level=2)
        print("Translation of columns for general attributes: ")
        print(table_report)
        utility.print_terminal_partition(level=3)
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
# Development code







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



def define_hormone_fields_measurement_reportability_missingness():
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

    # Return information.
    return records


def convert_hormone_concentration_units_moles_per_liter(
    table=None,
    factors_concentration=None,
):
    """
    Converts hormone concentrations to units of moles per liter (mol/L).

    UK Biobank field 30600, concentration in grams per liter (g/L) of albumin in
    blood

    UK Biobank field 30830, concentration in nanomoles per liter (nmol/L) of
    steroid hormone binding globulin (SHBG) in blood

    UK Biobank field 30800, concentration in picomoles per liter (pmol/L) of
    oestradiol in blood

    UK Biobank field 30850, concentration in nanomoles per liter (nmol/L) of
    total testosterone in blood

    UK Biobank field 30890, concentration in nanomoles per liter (nmol/L) of
    vitamin D in blood

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
        axis="columns", # apply across rows
    ) # 1 mole = 66472.2 g
    table["steroid_globulin"] = table.apply(
        lambda row: float(
            (row["30830-0.0"] / 1E9) * factors_concentration["steroid_globulin"]
        ),
        axis="columns", # apply across rows
    ) # 1 mol = 1E9 nmol
    table["oestradiol"] = table.apply(
        lambda row: float(
            (row["30800-0.0"] / 1E12) * factors_concentration["oestradiol"]
        ),
        axis="columns", # apply across rows
    ) # 1 mol = 1E12 pmol
    table["testosterone"] = table.apply(
        lambda row: float(
            (row["30850-0.0"] / 1E9) * factors_concentration["testosterone"]
        ),
        axis="columns", # apply across rows
    ) # 1 mol = 1E9 nmol
    table["vitamin_d"] = table.apply(
        lambda row: float(
            (row["30890-0.0"] / 1E9) * factors_concentration["vitamin_d"]
        ),
        axis="columns", # apply across rows
    ) # 1 mol = 1E9 nmol
    # Return information.
    return table

# TODO: TCW 30 July 2021
# TODO: I'm working on replacing this function with a function that returns a text report

def organize_report_column_pair_correlations(
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

    """

    table = table.copy(deep=True)
    table.dropna(
        axis="index",
        how="any",
        subset=[column_one, column_two],
        inplace=True,
    )
    pearson_correlation, pearson_probability = scipy.stats.pearsonr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    spearman_correlation, spearman_probability = scipy.stats.spearmanr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    # Report.
    utility.print_terminal_partition(level=2)
    print("Correlations between pair of columns")
    print("valid value pairs for correlation: " + str(table.shape[0]))
    print("column one: " + str(column_one))
    print("column_two: " + str(column_two))
    print("Pearson correlation: " + str(pearson_correlation))
    print("Pearson probability: " + str(pearson_probability))
    print("Spearman correlation: " + str(spearman_correlation))
    print("Spearman probability: " + str(spearman_probability))

    pass


# Imputation of missing hormone measurements


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
    UK Biobank data coding "782" for variable fields.
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

    Data-Field "30606": "Albumin reportability"
    Data-Field "30836": "SHBG reportability"
    Data-Field "30806": "Oestradiol reportability"
    Data-Field "30856": "Testosterone reportability"
    Data-Field "30896": "Vitamin D reportability"
    UK Biobank data coding "4917" for variable fields.
    1: "Reportable at assay and after aliquot correction, if attempted"
    2: "Reportable at assay but not reportable after any corrections (too low)"
    3: "Reportable at assay but not reportable after any corrections (too high)"
    4: "Not reportable at assay (too low)"
    5: "Not reportable at assay (too high)"

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
            axis="columns", # apply across rows
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
        reportability_low = str(
            str(hormone["name"]) + "_reportability_limit"
        )
        table[reportability_low] = table.apply(
            lambda row:
                interpret_biochemistry_reportability_detection_limit(
                    field_code_4917=row[hormone["reportability"]],
                ),
            axis="columns", # apply across rows
        )
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
            axis="columns", # apply across rows
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


def organize_hormones_missingness_imputation(
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
    table = determine_hormones_reportability_detection_limit(
        hormone_fields=hormone_fields,
        table=table,
    )
    table = determine_hormones_missingness_beyond_detection_range(
        hormone_fields=hormone_fields,
        table=table,
    )
    # Impute missing values.
    hormones = list()
    for record in hormone_fields:
        hormones.append(record["name"])
    table = impute_missing_hormones_detection_limit(
        hormones=hormones,
        table=table,
        report=report,
    )
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
        axis="columns", # apply across rows
    )
    table["testosterone_bioavailable"] = table.apply(
        lambda row:
            calculate_estimation_bioavailable_testosterone(
                testosterone_free=row["testosterone_free"],
                albumin=row["albumin"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
    )
    table["testosterone_bioavailable_imputation"] = table.apply(
        lambda row:
            calculate_estimation_bioavailable_testosterone(
                testosterone_free=row["testosterone_free_imputation"],
                albumin=row["albumin_imputation"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
    )

    # Return information.
    return table


# Main driver


def organize_sex_hormone_variables(
    table=None,
    report=None,
):
    """
    Organizes information about sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
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
        define_hormone_fields_measurement_reportability_missingness()
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
    table = convert_hormone_concentration_units_moles_per_liter(
        table=table,
        factors_concentration=factors_concentration,
    )

    ##########
    # Impute hormone measurements that were missing due to limit of detection.
    table = organize_hormones_missingness_imputation(
        hormone_fields=hormone_fields,
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
    # Transform variables' values to normalize distributions.
    columns_hormones_normalization = [
        "albumin", "albumin_imputation",
        "steroid_globulin", "steroid_globulin_imputation",
        "oestradiol", "oestradiol_imputation",
        "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
        "oestradiol_free", "oestradiol_free_imputation",
        "testosterone", "testosterone_imputation",
        "testosterone_bioavailable", "testosterone_bioavailable_imputation",
        "testosterone_free", "testosterone_free_imputation",
        "vitamin_d", "vitamin_d_imputation",
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
        ],
        axis="columns",
        inplace=True
    )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "albumin",
        "albumin_imputation",
        "steroid_globulin",
        "steroid_globulin_imputation",
        "oestradiol",
        "oestradiol_imputation",
        "oestradiol_bioavailable",
        "oestradiol_bioavailable_imputation",
        "oestradiol_free",
        "oestradiol_free_imputation",
        "testosterone",
        "testosterone_imputation",
        "testosterone_bioavailable",
        "testosterone_bioavailable_imputation",
        "testosterone_free",
        "testosterone_free_imputation",
        "vitamin_d",
        "vitamin_d_imputation",
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

        utility.print_terminal_partition(level=2)
        print("correlations in MALES!")
        table_male = table_clean.loc[
            (table_clean["sex_text"] == "male"), :
        ]
        utility.print_terminal_partition(level=3)
        organize_report_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_free",
            table=table_male,
        )
        organize_report_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_bioavailable",
            table=table_male,
        )
        organize_report_column_pair_correlations(
            column_one="testosterone_free",
            column_two="testosterone_bioavailable",
            table=table_male,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_free",
            table=table_male,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_bioavailable",
            table=table_male,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol_free",
            column_two="oestradiol_bioavailable",
            table=table_male,
        )

        utility.print_terminal_partition(level=2)
        print("correlations in FEMALES!")
        table_female = table.loc[
            (table["sex_text"] == "female"), :
        ]
        utility.print_terminal_partition(level=3)
        organize_report_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_free",
            table=table_female,
        )
        organize_report_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_bioavailable",
            table=table_female,
        )
        organize_report_column_pair_correlations(
            column_one="testosterone_free",
            column_two="testosterone_bioavailable",
            table=table_female,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_free",
            table=table_female,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_bioavailable",
            table=table_female,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol_free",
            column_two="oestradiol_bioavailable",
            table=table_female,
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


def interpret_menstruation_days(
    field_3700=None,
):
    """
    Intepret UK Biobank's coding for field 3700.

    Data-Field "3700": "Time since last menstrual period"
    UK Biobank data coding "100291" for variable field "3700".
    0 - 365: "days"
    -1: "Do not know"
    -3: "Prefer not to answer"

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
        (-3.5 <= field_3700 and field_3700 < 1000)
    ):
        # The variable has a valid value.
        if (0 <= field_3700 and field_3700 < 1000):
            # The variable has a valid value for days.
            value = float(field_3700)
        elif (-1.5 <= field_3700 and field_3700 < -0.5):
            # -1: "do not know"
            value = float("nan")
        elif (-3.5 <= field_3700 and field_3700 < -2.5):
            # -3: "prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def interpret_menstruation_cycle_duration(
    field_3710=None,
):
    """
    Intepret UK Biobank's coding for field 3710.

    Data-Field "3710": "Length of menstrual cycle"
    UK Biobank data coding "100582" for variable field "3710".
    7 - 365: "days"
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
        # null
        value = float("nan")
    # Return.
    return value


def interpret_menstruation_cycle_irregularity(
    field_3710=None,
):
    """
    Intepret UK Biobank's coding for field 3710.

    Only interpret whether the variable indicates that the person had an
    irregular menstrual cycle.

    Data-Field "3710": "Length of menstrual cycle"
    UK Biobank data coding "100582" for variable field "3710".
    7 - 365: "days"
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
            value = float("nan")
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
        # null
        value = float("nan")
    # Return.
    return value


def interpret_menstruation_period_current(
    field_3720=None,
):
    """
    Intepret UK Biobank's coding for field 3720.

    Data-Field "3720": "Menstruating today"
    UK Biobank data coding "100349" for variable field "3720".
    1: "Yes"
    0: "No"
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
        if (0.5 <= field_3720 and field_3720 < 1.5):
            # 1: "Yes"
            value = 1
        elif (-0.5 <= field_3720 and field_3720 < 0.5):
            # 0: "No"
            value = 0
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
        # null
        value = float("nan")
    # Return.
    return value


def interpret_menopause_hysterectomy(
    field_2724=None,
):
    """
    Intepret UK Biobank's coding for field 2724.

    Only interpret whether the variable indicates that the person had a
    hysterectomy.

    The only valid value is whether person indicated hysterectomy.

    Data-Field "2724": "Had menopause"
    UK Biobank data coding "100579" for variable field "2724".
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
        # null
        value = float("nan")
    # Return.
    return value


def interpret_menopause_natural_self_report(
    field_2724=None,
):
    """
    Intepret UK Biobank's coding for field 2724.

    Only interpret whether the variable indicates that the person self reported
    natural menopause, not involving hysterectomy.

    Data-Field "2724": "Had menopause"
    UK Biobank data coding "100579" for variable field "2724".
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
        # null
        value = float("nan")
    # Return.
    return value


def interpret_hysterectomy(
    field_3591=None,
):
    """
    Intepret UK Biobank's coding for field 3591.

    Data-Field "3591": "Ever had hysterectomy (womb removed)"
    UK Biobank data coding "100599" for variable field "3591".
    "no": 0
    "yes": 1
    "not sure": -5
    "prefer not to answer": -3

    Accommodate inexact float values.

    Interpret missing or null values as False.

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
            # 0: "no"
            value = 0
        elif (0.5 <= field_3591 and field_3591 < 1.5):
            # 1: "yes"
            value = 1
        elif (-5.5 <= field_3591 and field_3591 < -4.5):
            # -5: "not sure"
            value = float("nan")
        elif (-3.5 <= field_3591 and field_3591 < -2.5):
            # -3: "prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def interpret_oophorectomy(
    field_2834=None,
):
    """
    Intepret UK Biobank's coding for field 2834.

    Data-Field "2834": "Bilateral oophorectomy (both ovaries
    removed)"
    UK Biobank data coding "100599" for variable field "2834".
    "no": 0
    "yes": 1
    "not sure": -5
    "prefer not to answer": -3

    Accommodate inexact float values.

    Interpret missing or null values as False.

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
            # 0: "no"
            value = 0
        elif (0.5 <= field_2834 and field_2834 < 1.5):
            # 1: "yes"
            value = 1
        elif (-5.5 <= field_2834 and field_2834 < -4.5):
            # -5: "not sure"
            value = float("nan")
        elif (-3.5 <= field_2834 and field_2834 < -2.5):
            # -3: "prefer not to answer"
            value = float("nan")
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def interpret_pregnancy(
    field_3140=None,
):
    """
    Intepret UK Biobank's coding for field 3140.

    Interpret "unsure" pregnancy as potential pregnancy.

    Data-Field "3140": "Pregnant"
    UK Biobank data coding "100267" for variable field "3140".
    "no": 0
    "yes": 1
    "unsure": 2

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
            # 0: "no"
            value = 0
        elif (0.5 <= field_3140 and field_3140 < 1.5):
            # 1: "yes"
            value = 1
        elif (1.5 <= field_3140 and field_3140 < 2.5):
            # 2: "unsure"
            # Interpret "unsure" pregnancy as potential pregnancy.
            value = 1
        else:
            # uninterpretable
            value = float("nan")
    else:
        # null
        value = float("nan")
    # Return.
    return value


def interpret_recent_hormone_alteration_therapies(
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
            # 1: "yes"
            # Person has used therapy.
            # Determine whether person used therapy recently.
            if (
                (not pandas.isna(age_stopped)) and
                (-11.5 <= age_stopped and age_stopped <= age)
            ):
                # Variable "age_stopped" has a valid value.
                if (-11.5 <= age_stopped and age_stopped < -10.5):
                    # -11: "still taking the pill"
                    # -11: "still taking HRT"
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
                    # -1: "do not know"
                    value = float("nan")
                elif (-3.5 <= age_stopped and age_stopped < -2.5):
                    # -3: "prefer not to answer"
                    value = float("nan")
                else:
                    # uninterpretable
                    value = float("nan")
            else:
                # null
                value = float("nan")
        elif (-0.5 <= ever_used and ever_used < 0.5):
            # 0: "no"
            # Person has not used therapy.
            value = 0
        elif (-1.5 <= ever_used and ever_used < -0.5):
            # -1: "do not know"
            value = float("nan")
        elif (-3.5 <= ever_used and ever_used < -2.5):
            # -3: "prefer not to answer"
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

    Data-Field "2784": "ever taken oral contraceptive pill"
    UK Biobank data coding "100349" for variable field "2784".
    "yes": 1
    "no": 0
    "do not know": -1
    "prefer not to answer": -3

    Data-Field "2794": "age started oral contraceptive pill"
    UK Biobank data coding "100291" for variable field "2794".
    "age in years": 5 - person's age
    "do not know": -1
    "prefer not to answer": -3

    Data-Field "2804": "age when last used oral contraceptive pill"
    UK Biobank data coding "100595" for variable field "2804".
    "age in years": 5 - person's age
    "do not know": -1
    "prefer not to answer": -3
    "still taking the pill": -11

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
    value = interpret_recent_hormone_alteration_therapies(
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

    Data-Field "2814": "ever used hormone-replacement therapy (HRT)"
    UK Biobank data coding "100349" for variable field "2814".
    "yes": 1
    "no": 0
    "do not know": -1
    "prefer not to answer": -3

    Data-Field "3536": "age started hormone-replacement therapy (HRT)"
    UK Biobank data coding "100291" for variable field "3536".
    "age in years": 16 - person's age
    "do not know": -1
    "prefer not to answer": -3

    Data-Field "3546": "age last used hormone-replacement therapy (HRT)"
    UK Biobank data coding "100598" for variable field "3546".
    "age in years": 20 - person's age
    "do not know": -1
    "prefer not to answer": -3
    "still taking HRT": -11

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
    value = interpret_recent_hormone_alteration_therapies(
        age=age,
        recent_range=recent_range,
        ever_used=field_2814,
        age_started=field_3536,
        age_stopped=field_3546,
    )
    # Return.
    return value


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
            single_false_sufficient=True, # both variables have similar meaning
        )
    else:
        # Hysterectomy is undefined for males.
        value = float("nan")
    # Return information.
    return value


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
            single_false_sufficient=False,
        )
    else:
        # Hysterectomy and oophorectomy are undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_menopause_natural_self_report(
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
    menopause_natural_self_report = interpret_menopause_natural_self_report(
        field_2724=field_2724,
    )
    # Comparison.
    # Only define days since previous menstruation for female persons.
    if (
        (sex_text == "female")
    ):
        # Female.
        value = menopause_natural_self_report
    else:
        # Male.
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_menopause_binary(
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
    menopause_natural = interpret_menopause_natural_self_report(
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


def determine_female_menopause_binary_strict(
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
    menopause_natural = interpret_menopause_natural_self_report(
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


def determine_female_menopause_ordinal(
    sex_text=None,
    age=None,
    threshold_age_pre=None,
    threshold_age_post=None,
    menstruation_duration=None,
    menstruation_irregularity=None,
    menstruation_days=None,
    threshold_menstruation_pre=None,
    threshold_menstruation_post=None,
    menopause_natural_self_report=None,
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
        menopause_natural_self_report (float): binary representation of whether
            person self reported natural menopause
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
                (not pandas.isna(menopause_natural_self_report)) and
                (menopause_natural_self_report == 1)
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
                    (pandas.isna(menopause_natural_self_report)) or
                    (menopause_natural_self_report == 0)
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
                    (pandas.isna(menopause_natural_self_report)) or
                    (menopause_natural_self_report == 0)
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


def determine_female_menopause_ordinal_strict(
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
    menopause_natural = interpret_menopause_natural_self_report(
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


def determine_female_menstruation_days_threshold(
    sex_text=None,
    menopause_ordinal=None,
    menstruation_days=None,
    threshold_premenopause=None,
    threshold_perimenopause=None,
):
    """
    Determine a female person's categorical menstruation phase.

    0: follicular phase of menstruation
    1: luteal phase of menstruation

    arguments:
        sex_text (str): textual representation of sex selection
        menopause_ordinal (float): ordinal representation of whether female
            person has experienced menopause
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
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
            # Variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_menstruation_phase(
    sex_text=None,
    menopause_ordinal=None,
    menstruation_days=None,
    threshold_premenopause=None,
    threshold_perimenopause=None,
):
    """
    Determine a female person's categorical menstruation phase.

    0: follicular phase of menstruation
    1: luteal phase of menstruation

    arguments:
        sex_text (str): textual representation of sex selection
        menopause_ordinal (float): ordinal representation of whether female
            person has experienced menopause
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
        threshold_premenopause (int): threshold in days for premenopause female
            persons, below which to consider follicular phase of menstruation
            cycle
        threshold_perimenopause (int): threshold in days for perimenopause
            female persons, below which to consider follicular phase of
            menstruation cycle

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
            # Determine categorical phase of menstruation cycle for premenopause
            # females.
            if (menopause_ordinal == 0):
                if (menstruation_days < threshold_premenopause):
                    # Person qualifies for follicular phase of menstruation
                    # cycle.
                    value = 0
                elif (menstruation_days >= threshold_premenopause):
                    # Person qualitifies for luteal phase of menstruation cycle.
                    value = 1
                else:
                    # Persons does not qualify for any categories.
                    value = float("nan")
            # Determine categorical phase of menstruation cycle for
            # perimenopause females.
            elif (menopause_ordinal == 1):
                if (menstruation_days < threshold_perimenopause):
                    # Person qualifies for follicular phase of menstruation
                    # cycle.
                    value = 0
                elif (menstruation_days >= threshold_perimenopause):
                    # Person qualitifies for luteal phase of menstruation cycle.
                    value = 1
                else:
                    # Persons does not qualify for any categories.
                    value = float("nan")
            # Determine categorical phase of menstruation cycle for
            # postmenopause females.
            elif (menopause_ordinal == 2):
                # Menstruation undefined for postmenopause females.
                value = float("nan")
            else:
                # Person does not qualify for any categories.
                value = float("nan")
        else:
            # Variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_menstruation_phase_early_late(
    sex_text=None,
    menopause_ordinal=None,
    menstruation_days=None,
    threshold_middle_follicular=None,
    threshold_ovulation=None,
    threshold_middle_luteal=None,
):
    """
    Determine a female person's categorical menstruation phase.

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
        menopause_ordinal (float): ordinal representation of whether female
            person has experienced menopause
        menstruation_days (int): count of days since previous menstruation
            (menstrual period)
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
            (not pandas.isna(menopause_ordinal)) and
            (not pandas.isna(menstruation_days))
        ):
            # Determine categorical phase of menstruation cycle for premenopause
            # and perimenopause females.
            if ((menopause_ordinal == 0) or (menopause_ordinal == 1)):
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
            # Determine categorical phase of menstruation cycle for
            # postmenopause females.
            elif (menopause_ordinal == 2):
                # Menstruation undefined for postmenopause females.
                value = float("nan")
            else:
                # Person does not qualify for any categories.
                value = float("nan")
        else:
            # Variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_menstruation_phase_cycle(
    sex_text=None,
    menopause_ordinal=None,
    menstruation_phase_early_late=None,
):
    """
    Determine a female person's menstruation phase from a cyclical perspective.

    0: late luteal or early follicular phase of menstruation cycle
    - "shoulder"
    1: late follicular or early luteal phase of menstruation cycle
    - "ovulation"

    arguments:
        sex_text (str): textual representation of sex selection
        menopause_ordinal (float): ordinal representation of whether female
            person has experienced menopause
        menstruation_phase_early_late (int): ordinal representation of phases of
            the menstruation cycle

    raises:

    returns:
        (float): interpretation value

    """

    # Determine categorical phase of menstruation cycle.
    if (sex_text == "female"):
        # Determine whether any variables have missing values.
        if (
            (not pandas.isna(menopause_ordinal)) and
            (not pandas.isna(menstruation_phase_early_late))
        ):
            # Determine categorical phase of menstruation cycle for premenopause
            # and perimenopause females.
            if ((menopause_ordinal == 0) or (menopause_ordinal == 1)):
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
            # Determine categorical phase of menstruation cycle for
            # postmenopause females.
            elif (menopause_ordinal == 2):
                # Menstruation undefined for postmenopause females.
                value = float("nan")
            else:
                # Person does not qualify for any categories.
                value = float("nan")
        else:
            # Variables have missing values.
            value = float("nan")
    else:
        # Menstruation undefined for males.
        value = float("nan")
    # Return information.
    return value


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


def determine_female_oral_contraception(
    sex_text=None,
    age=None,
    recent_range=None,
    field_2784=None,
    field_2794=None,
    field_2804=None,
):
    """
    Determine whether person used oral contraception recently or currently.

    arguments:
        sex_text (str): textual representation of sex selection
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

    # Interpret oral contraception.
    contraception = interpret_recent_oral_contraception(
        age=age,
        recent_range=recent_range, # recent range in years
        field_2784=field_2784,
        field_2794=field_2794,
        field_2804=field_2804,
    )
    # Comparison.
    if (sex_text == "female"):
        value = contraception
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
):
    """
    Determine count of days since previous menstruation (menstrual period).

    This function uses a specific definition of pregnancy that considers
    menopause and does not include uncertain cases.

    arguments:
        sex_text (str): textual representation of sex selection
        age (int): age of person in years
        recent_range (int): years within current age to consider recent
        field_2814 (float): UK Biobank field 2784, ever taken hormone
            replacement therapy
        field_3536 (float): UK Biobank field 2794, age started hormone
            replacement therapy (HRT)
        field_3546 (float): UK Biobank field 2804, age stopped hormone
            replacement therapy (HRT)

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
    # Comparison.
    if (sex_text == "female"):
        value = therapy
    else:
        # This specific variable is undefined for males.
        value = float("nan")
    # Return information.
    return value


def determine_female_any_hormone_alteration_medication(
    sex_text=None,
    oral_contraception=None,
    hormone_therapy=None,
):
    """
    Determine count of days since previous menstruation (menstrual period).

    This function uses a specific definition of pregnancy that considers
    menopause and does not include uncertain cases.

    arguments:
        sex_text (str): textual representation of sex selection
        oral_contraception (float): binary logical representation of whether
            person used oral contraception recently
        hormone_therapy (float): binary logical representation of whether
            person used hormone replacement therapy recently

    raises:

    returns:
        (float): binary logical representation of whether person used any
            hormone-alteration medications recently

    """

    # Comparison.
    if (sex_text == "female"):
        value = utility.determine_logical_or_combination_binary_missing(
            first=oral_contraception,
            second=hormone_therapy,
            single_false_sufficient=False,
        )
    else:
        # This specific variable is undefined for males.
        value = float("nan")
    # Return information.
    return value


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
        "2784-0.0", "2794-0.0", "2804-0.0",
        "2814-0.0", "3536-0.0", "3546-0.0",
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
        axis="columns", # apply across rows
    )

    # Determine count of days in person's usual menstrual cycle.
    table["menstruation_cycle_duration"] = table.apply(
        lambda row:
            determine_female_menstruation_cycle_duration(
                sex_text=row["sex_text"],
                field_3710=row["3710-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person reports irregularity in menstrual cycle.
    table["menstruation_cycle_irregularity"] = table.apply(
        lambda row:
            determine_female_menstruation_cycle_irregularity(
                sex_text=row["sex_text"],
                field_3710=row["3710-0.0"],
            ),
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
    )
    # Determine whether female persons experienced bilateral oophorectomy.
    table["oophorectomy"] = table.apply(
        lambda row:
            determine_female_oophorectomy(
                sex_text=row["sex_text"],
                field_2834=row["2834-0.0"],
            ),
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
    )

    # Determine whether female persons self reported natural menopause, which
    # neither involves hysterectomy nor oophorectomy.
    table["menopause_natural_self_report"] = table.apply(
        lambda row:
            determine_female_menopause_natural_self_report(
                sex_text=row["sex_text"],
                field_2724=row["2724-0.0"],
            ),
        axis="columns", # apply across rows
    )

    # Determine whether female persons have experienced menopause.
    # This definition is binary.
    table["menopause_binary"] = table.apply(
        lambda row:
            determine_female_menopause_binary(
                sex_text=row["sex_text"],
                age=row["age"],
                threshold_age=51, # threshold age in years
                menstruation_days=row["menstruation_days"],
                threshold_menstruation_days=90, # threshold in days
                field_2724=row["2724-0.0"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
            ),
        axis="columns", # apply across rows
    )
    table["menopause_binary_strict"] = table.apply(
        lambda row:
            determine_female_menopause_binary_strict(
                sex_text=row["sex_text"],
                age=row["age"],
                threshold_age=51, # threshold age in years
                menstruation_days=row["menstruation_days"],
                threshold_menstruation_days=90, # threshold in days
                field_2724=row["2724-0.0"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
            ),
        axis="columns", # apply across rows
    )

    # Determine whether female persons have experienced menopause.
    # This definition is ordinal.
    table["menopause_ordinal"] = table.apply(
        lambda row:
            determine_female_menopause_ordinal(
                sex_text=row["sex_text"],
                age=row["age"],
                threshold_age_pre=47, # threshold age in years
                threshold_age_post=55, # threshold age in years
                menstruation_duration=row["menstruation_cycle_duration"],
                menstruation_irregularity=(
                    row["menstruation_cycle_irregularity"]
                ),
                menstruation_days=row["menstruation_days"],
                threshold_menstruation_pre=60, # threshold in days
                threshold_menstruation_post=365, # threshold in days
                menopause_natural_self_report=(
                    row["menopause_natural_self_report"]
                ),
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
            ),
        axis="columns", # apply across rows
    )
    table["menopause_ordinal_strict"] = table.apply(
        lambda row:
            determine_female_menopause_ordinal_strict(
                sex_text=row["sex_text"],
                age=row["age"],
                threshold_age_pre=47, # threshold age in years
                threshold_age_post=55, # threshold age in years
                menstruation_days=row["menstruation_days"],
                threshold_menstruation_days_pre=60, # threshold in days
                threshold_menstruation_days_post=365, # threshold in days
                field_2724=row["2724-0.0"],
                hysterectomy=row["hysterectomy"],
                oophorectomy=row["oophorectomy"],
            ),
        axis="columns", # apply across rows
    )

    ##########
    # Menstrual cycle

    # Set thresholds on days of menstrual cycle that are specific to
    # premenopause or perimenopause status.
    # Notice that these thresholds differ from those for definition of
    # menopause status.
    # Rather, these thresholds clean up the menopause cohorts after definition.
    table["menstruation_days_threshold"] = table.apply(
        lambda row:
            determine_female_menstruation_days_threshold(
                sex_text=row["sex_text"],
                menopause_ordinal=row["menopause_ordinal"],
                menstruation_days=row["menstruation_days"],
                threshold_premenopause=31, # threshold in days
                threshold_perimenopause=41, # threshold in days
            ),
        axis="columns", # apply across rows
    )

    # Determine female persons' categorical menstruation phase.
    table["menstruation_phase"] = table.apply(
        lambda row:
            determine_female_menstruation_phase(
                sex_text=row["sex_text"],
                menopause_ordinal=row["menopause_ordinal"],
                menstruation_days=row["menstruation_days_threshold"],
                threshold_premenopause=15, # threshold in days
                threshold_perimenopause=15, # threshold in days
            ),
        axis="columns", # apply across rows
    )

    # Determine female persons' ordinal menstruation phase.
    table["menstruation_phase_early_late"] = table.apply(
        lambda row:
            determine_female_menstruation_phase_early_late(
                sex_text=row["sex_text"],
                menopause_ordinal=row["menopause_ordinal"],
                menstruation_days=row["menstruation_days_threshold"],
                threshold_middle_follicular=5, # threshold in days
                threshold_ovulation=15,
                threshold_middle_luteal=20, # threshold in days
            ),
        axis="columns", # apply across rows
    )

    # Determine female persons' ordinal menstruation phase from cyclical
    # perspective.
    table["menstruation_phase_cycle"] = table.apply(
        lambda row:
            determine_female_menstruation_phase_cycle(
                sex_text=row["sex_text"],
                menopause_ordinal=row["menopause_ordinal"],
                menstruation_phase_early_late=(
                    row["menstruation_phase_early_late"]
                ),
            ),
        axis="columns", # apply across rows
    )

    # Determine whether female persons were pregnant.
    table["pregnancy"] = table.apply(
        lambda row:
            determine_female_pregnancy(
                sex_text=row["sex_text"],
                field_3140=row["3140-0.0"],
            ),
        axis="columns", # apply across rows
    )

    ##########
    # Medications and therapies that alter hormone levels

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
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person was using hormone replacement therapy recently.
    table["hormone_replacement"] = table.apply(
        lambda row:
            determine_female_hormone_replacement_therapy(
                sex_text=row["sex_text"],
                age=row["age"],
                recent_range=1, # years within current age to consider recent
                field_2814=row["2814-0.0"],
                field_3536=row["3536-0.0"],
                field_3546=row["3546-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person was using any hormone-altering medications
    # recently.
    table["hormone_alteration"] = table.apply(
        lambda row:
            determine_female_any_hormone_alteration_medication(
                sex_text=row["sex_text"],
                oral_contraception=row["oral_contraception"],
                hormone_therapy=row["hormone_replacement"],
             ),
        axis="columns", # apply across rows
    )
    # Determine combination categories by menopause and hormone-atering therapy.
    # Report on 28 April 2021.
    # Category 1: 24158
    # Category 2: 178007
    # Category 3: 8438
    # Category 4: 51786
    table = (
        utility.determine_binary_categorical_products_of_two_binary_variables(
            table=table,
            first="menopause_binary",
            second="hormone_alteration",
            prefix="menopause_hormone_category",
            report=report,
    ))

    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "2724-0.0", "3591-0.0", "2834-0.0",
            "3140-0.0",
            "3700-0.0", "3710-0.0", "3720-0.0",
            "2784-0.0", "2794-0.0", "2804-0.0",
            "2814-0.0", "3536-0.0", "3546-0.0",
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
        "hysterectomy", "oophorectomy", "hysterectomy_or_oophorectomy",
        "menopause_binary", "menopause_binary_strict",
        "menopause_ordinal", "menopause_ordinal_strict",
        "menstruation_days", "menstruation_days_threshold",
        "menstruation_cycle_duration", "menstruation_cycle_irregularity",
        "menstruation_phase",
        "menstruation_phase_early_late", "menstruation_phase_cycle",
        "pregnancy",
        "oral_contraception", "hormone_replacement", "hormone_alteration",
        "menopause_hormone_category_1", "menopause_hormone_category_2",
        "menopause_hormone_category_3", "menopause_hormone_category_4",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    table_report.sort_values(
        by=["pregnancy", "menopause_binary"],
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
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Psychological variables


def interpret_import_bipolar_disorder(
    import_bipolar_disorder=None,
):
    """
    Intepret import variable for bipolar disorder cases and controls.

    Accommodate inexact float values.

    arguments:
        import_bipolar_disorder (float): import variable "bipolar.cc" from
            variable definitions by Brandon J. Coombes

    raises:

    returns:
        (bool): interpretation value

    """

    # Interpret field code.
    if (
        (not pandas.isna(import_bipolar_disorder)) and
        (-0.5 <= import_bipolar_disorder and import_bipolar_disorder < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= import_bipolar_disorder and import_bipolar_disorder < 0.5):
            # 0: "control"
            value = 0
        elif (0.5 <= import_bipolar_disorder and import_bipolar_disorder < 1.5):
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


def organize_psychology_variables(
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
    # Translate column names.
    translations = dict()
    translations["20127-0.0"] = "neuroticism"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Convert variable types.
    columns_type = [
        "neuroticism",
        "bipolar.cc",
    ]
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=["neuroticism"],
        table=table,
    )
    # Determine whether persons qualify as cases or controls for bipolar
    # disorder.
    table["bipolar_disorder"] = table.apply(
        lambda row:
            interpret_import_bipolar_disorder(
                import_bipolar_disorder=row["bipolar.cc"],
            ),
        axis="columns", # apply across rows
    )

    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            #"20127-0.0",
            "icd_bipolar", "bipolar.diag2", "Smithbipolar", "SmithMood",
            "MHQ.bipolar.Definition", "bipolar", "bipolar.cc",
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
        # Categorical ancestry and ethnicity.
        utility.print_terminal_partition(level=2)
        print("Bipolar Disorder cases and controls")
        table_cases = table.loc[
            (table["bipolar_disorder"] == 1), :
        ]
        table_controls = table.loc[
            (table["bipolar_disorder"] == 0), :
        ]
        print("bipolar disorder cases: " + str(table_cases.shape[0]))
        print("bipolar disorder controls: " + str(table_controls.shape[0]))

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


def interpret_alcohol_consumption_frequency(
    field_1558=None,
):
    """
    Intepret UK Biobank's data-coding for data-field 1558.

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
        (float): binary representation of whether person ever consumed alcohol

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
        alcohol_ever (float): whether person ever consumed alcohol

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
        axis="columns", # apply across rows
    )
    # Determine person's alcohol consumption status.
    table["alcohol_status"] = table.apply(
        lambda row:
            determine_alcohol_consumption_status(
                field_20117=row["20117-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine person's former alcohol consumption status.
    table["alcohol_former"] = table.apply(
        lambda row:
            determine_former_alcohol_consumption(
                field_3731=row["3731-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine person's former alcohol consumption status.
    table["alcohol_trend"] = table.apply(
        lambda row:
            determine_alcohol_consumption_trend(
                field_1628=row["1628-0.0"],
            ),
        axis="columns", # apply across rows
    )

    # Determine whether person consumes alcohol currently.
    table["alcohol_current"] = table.apply(
        lambda row:
            determine_alcohol_consumption_current(
                alcohol_frequency=row["alcohol_frequency"],
                alcohol_status=row["alcohol_status"],
                alcohol_trend=row["alcohol_trend"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person consumes alcohol not currently but previously.
    #table["alcohol_previous"]
    # Determine whether person never consumes any alcohol.
    table["alcohol_ever"] = table.apply(
        lambda row:
            determine_alcohol_consumption_ever(
                alcohol_frequency=row["alcohol_frequency"],
                alcohol_status=row["alcohol_status"],
                alcohol_former=row["alcohol_former"],
                alcohol_trend=row["alcohol_trend"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person never consumes any alcohol, either currently or
    # previously.
    table["alcohol_never"] = table.apply(
        lambda row:
            determine_alcohol_consumption_never(
                alcohol_ever=row["alcohol_ever"],
            ),
        axis="columns", # apply across rows
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
            consumes alcohol currently
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
        (-0.5 <= alcohol_current and alcohol_current < 0.5)
    ):
        # Person never consumes any alcohol currently.
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
    table["drinks_weekly"] = table.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["1588-0.0"],
                wine_red=row["1568-0.0"],
                wine_white=row["1578-0.0"],
                port=row["1608-0.0"],
                liquor=row["1598-0.0"],
                other=row["5364-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Calculate sum of drinks monthly.
    table["drinks_monthly"] = table.apply(
        lambda row:
            calculate_sum_drinks(
                beer_cider=row["4429-0.0"],
                wine_red=row["4407-0.0"],
                wine_white=row["4418-0.0"],
                port=row["4451-0.0"],
                liquor=row["4440-0.0"],
                other=row["4462-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine sum of total drinks monthly.
    table["alcohol_drinks_monthly"] = table.apply(
        lambda row:
            determine_total_alcohol_consumption_monthly(
                alcohol_current=row["alcohol_current"],
                drinks_weekly=row["drinks_weekly"],
                drinks_monthly=row["drinks_monthly"],
                weeks_per_month=4.345, # 52.143 weeks per year (12 months)
            ),
        axis="columns", # apply across rows
    )
    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "1588-0.0", "1568-0.0", "1578-0.0", "1608-0.0", "1598-0.0",
            "5364-0.0",
            "4429-0.0", "4407-0.0", "4418-0.0", "4451-0.0", "4440-0.0",
            "4462-0.0",
            "drinks_weekly",
            "drinks_monthly",
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
            "drinks_weekly",
            "drinks_monthly",
            "alcohol_drinks_monthly",
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
    if False:
        pail_quantity = organize_alcohol_consumption_quantity_variables(
            table=pail_consumption["table"],
            report=report,
        )
    # Return information.
    return pail_consumption


##########
# Alcohol AUDIT questionnaire


# TODO: TCW 3 August 2021
# TODO: this code block of AUDIT-C and AUDIT-P looks pretty good... but good idea to double check.

def interpret_alcohol_audit_one(
    value=None,
):
    """
    Intepret UK Biobank's coding for AUDIT questionnaire question 1.

    "audit_1", field "20414": "Frequency of drinking alcohol"
    UK Biobank data coding "521" for variable field "20414".
    "never": 0
    "monthly or less": 1
    "two to four times a month": 2
    "two to three times a week": 3
    "four or more times a week": 4
    "prefer not to answer": -818

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
        (-0.5 <= value and value < 4.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= value and value < 0.5):
            # "never"
            value_clean = 0
        elif (0.5 <= value and value < 1.5):
            # "monthly or less"
            value_clean = 1
        elif (1.5 <= value and value < 2.5):
            # "two to four times a month"
            value_clean = 2
        elif (2.5 <= value and value < 3.5):
            # "two to three times a week"
            value_clean = 3
        elif (3.5 <= value and value < 4.5):
            # "four or more times a week"
            value_clean = 4
    else:
        # "prefer not to answer" or null
        value_clean = float("nan")
    # Return information.
    return value_clean


def interpret_alcohol_audit_two_to_eight(
    value=None,
):
    """
    Intepret UK Biobank's coding for AUDIT questionnaire questions 2 to 8.

    UK Biobank data coding "522" for variable field "20403".
    "one or two": 1
    "three or four": 2
    "five or six": 3
    "seven, eight, or nine": 4
    "ten or more": 5
    "prefer not to answer": -818

    UK Biobank data coding "523" for variable fields "20416", "20413", "20407",
    "20412", "20409", and "20408".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

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
        (0.5 <= value and value < 5.5)
    ):
        # The variable has a valid value.
        if (0.5 <= value and value < 1.5):
            # code "522": "one or two"
            # code "523": "never"
            value_clean = 0
        elif (1.5 <= value and value < 2.5):
            # code "522": "three or four"
            # code "523": "less than monthly"
            value_clean = 1
        elif (2.5 <= value and value < 3.5):
            # code "522": "five or six"
            # code "523": "monthly"
            value_clean = 2
        elif (3.5 <= value and value < 4.5):
            # code "522": "seven, eight, or nine"
            # code "523": "weekly"
            value_clean = 3
        elif (4.5 <= value and value < 5.5):
            # code "522": "ten or more"
            # code "523": "daily or almost daily"
            value_clean = 4
    else:
        # "prefer not to answer" or null
        value_clean = float("nan")
    # Return information.
    return value_clean


def interpret_alcohol_audit_nine_ten(
    value=None,
):
    """
    Intepret UK Biobank's coding for AUDIT questionnaire questions 9 and 10.

    UK Biobank data coding "524" for variable fields "20411", and "20405".
    "no": 0
    "yes, but not in the last year": 1
    "yes, during the last year": 2
    "prefer not to answer": -818

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
        (-0.5 <= value and value < 2.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= value and value < 0.5):
            # code "524": "no"
            value_clean = 0
        elif (0.5 <= value and value < 1.5):
            # code "524": "yes, but not in the last year"
            value_clean = 2
        elif (1.5 <= value and value < 2.5):
            # code "524": "yes, during the last year"
            value_clean = 4
    else:
        # "prefer not to answer" or null
        value_clean = float("nan")
    # Return information.
    return value_clean


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
        (not math.isnan(audit_1_clean)) and
        (not math.isnan(audit_2_clean)) and
        (not math.isnan(audit_3_clean))
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

    "audit_1", field "20414": "Frequency of drinking alcohol"
    UK Biobank data coding "521" for variable field "20414".
    "never": 0
    "monthly or less": 1
    "two to four times a month": 2
    "two to three times a week": 3
    "four or more times a week": 4
    "prefer not to answer": -818

    "audit_2", field "20403": "Amount of alcohol drunk on a typical drinking
    day"
    UK Biobank data coding "522" for variable field "20403".
    "one or two": 1
    "three or four": 2
    "five or six": 3
    "seven, eight, or nine": 4
    "ten or more": 5
    "prefer not to answer": -818

    "audit_3", field "20416": "Frequency of consuming six or more units of
    alcohol"
    UK Biobank data coding "523" for variable field "20416".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

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
        axis="columns", # apply across rows
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
        (not math.isnan(audit_4_clean)) and
        (not math.isnan(audit_5_clean)) and
        (not math.isnan(audit_6_clean)) and
        (not math.isnan(audit_7_clean)) and
        (not math.isnan(audit_8_clean)) and
        (not math.isnan(audit_9_clean)) and
        (not math.isnan(audit_10_clean))
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

    "audit_4", field "20413": "Frequency of inability to cease drinking in last
    year"
    UK Biobank data coding "523" for variable field "20413".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

    "audit_5", field "20407": "Frequency of failure to fulfil normal
    expectations due to drinking alcohol in last year"
    UK Biobank data coding "523" for variable field "20407".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

    "audit_6", field "20412": "Frequency of needing morning drink of alcohol
    after heavy drinking session in last year"
    UK Biobank data coding "523" for variable field "20412".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

    "audit_7", field "20409": "Frequency of feeling guilt or remorse after
    drinking alcohol in last year"
    UK Biobank data coding "523" for variable field "20409".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

    "audit_8", field "20408": "Frequency of memory loss due to drinking alcohol
    in last year"
    UK Biobank data coding "523" for variable field "20408".
    "never": 1
    "less than monthly": 2
    "monthly": 3
    "weekly": 4
    "daily or almost daily": 5
    "prefer not to answer": -818

    "audit_9", field "20411": "Ever been injured or injured someone else
    through drinking alcohol"
    UK Biobank data coding "524" for variable field "20411".
    "no": 0
    "yes, but not in the last year": 1
    "yes, during the last year": 2
    "prefer not to answer": -818

    "audit_10", field "20405": "Ever had known person concerned about, or
    recommend reduction of, alcohol consumption"
    UK Biobank data coding "524" for variable field "20405".
    "no": 0
    "yes, but not in the last year": 1
    "yes, during the last year": 2
    "prefer not to answer": -818

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
        axis="columns", # apply across rows
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
        (not math.isnan(alcohol_auditc)) and
        (not math.isnan(alcohol_auditp))
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
        axis="columns", # apply across rows
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
        table=pail_auditc["table_clean"],
        report=report,
    )
    # Organize information about alcohol AUDIT.
    pail_audit = organize_alcohol_audit_variables(
        table=pail_auditp["table_clean"],
        report=report,
    )
    # Collect information.
    pail = dict()
    pail["auditc"] = pail_auditc
    pail["auditp"] = pail_auditp
    pail["audit"] = pail_audit
    # Return information.
    return pail


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


def parse_field_array_codes(
    collection=None,
    delimiter=None,
):
    """
    Parse the field's array codes.

    arguments:
        collection (str): raw string of field's array codes
        delimiter (str): delimiter between codes

    raises:

    returns:
        (list<str>): field's array codes

    """

    if (
        (len(str(collection)) > 0) and (delimiter in str(collection))
    ):
        codes_raw = str(collection).split(delimiter)
        codes = list()
        for code_raw in codes_raw:
            code = str(code_raw).strip()
            if (len(code) > 0):
                codes.append(str(code))
            pass
    else:
        codes = list()
    return codes


def parse_interpret_match_diagnosis_codes(
    row=None,
    fields=None,
    codes_match=None,
):
    """
    Parses and interprets ICD diagnosis codes.

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

    match = False
    for field in fields:
        collection = row[field]
        codes = parse_field_array_codes(
            collection=collection,
            delimiter=";",
        )
        if (len(codes) > 0):
            for code in codes:
                if (code in codes_match):
                    match = True
            pass
        pass
    return match


def determine_diagnosis_icd_alcoholism_group(
    row=None,
    fields_icd_9=None,
    fields_icd_10=None,
    codes_icd_group=None,
):
    """
    Organizes information about diagnoses relevant to alcoholism.

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
        match = True
    else:
        match = False
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
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
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
        axis="columns", # apply across rows
    )
    # Determine whether person has self diagnoses.
    table["alcohol_diagnosis_self"] = table.apply(
        lambda row:
            parse_interpret_match_diagnosis_codes(
                row=row.copy(deep=True),
                fields=["20002_array"],
                codes_match=codes_self,
            ),
        axis="columns", # apply across rows
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

# TODO: change alcohol_none to True/False/None like the diagnosis flags
def determine_control_alcoholism(
    alcohol_none=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditc_case=None,
    threshold_auditp_control=None,
    threshold_auditp_case=None,
    threshold_audit_control=None,
    threshold_audit_case=None,
):
    """
    Determines whether person qualifies as a control for alcoholism.

    arguments:
        alcohol_none (float): binary representation of whether person has never
            consumed alcohol, formerly or currently
        alcohol_diagnosis_a (bool): whether person has diagnosis in alcoholism
            group A
        alcohol_diagnosis_b (bool): whether person has diagnosis in alcoholism
            group B
        alcohol_diagnosis_c (bool): whether person has diagnosis in alcoholism
            group C
        alcohol_diagnosis_d (bool): whether person has diagnosis in alcoholism
            group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-C score
        threshold_auditc_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-C score
        threshold_auditp_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-P score
        threshold_auditp_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-P score
        threshold_audit_control (float): diagnostic control threshold (less
            than or equal) for AUDIT score
        threshold_audit_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT score

    raises:

    returns:
        (bool): whether person qualifies as a control for alcoholism

    """

    # Determine whether person has consumed alcohol, previously or currently.
    # Only persons who have consumed alcohol, previously or currently, ought to
    # qualify as controls for alcoholism.
    if (math.isnan(alcohol_none)):
        match_consumption = False
    else:
        if (alcohol_none < 0.5):
            # Person has consumed alcohol, previously or currently.
            match_consumption = True
        else:
            # Person has never consumed alcohol, previously or currently.
            match_consumption = False

    # Determine whether person's AUDIT scores are within thresholds.
    if (math.isnan(alcohol_audit)):
        comparison_audit = False
    else:
        if (alcohol_audit <= threshold_audit_control):
            comparison_audit = True
        else:
            comparison_audit = False
    # Interpret comparisons from all AUDIT socre combinations.
    if (comparison_audit):
        match_audit = True
    else:
        match_audit = False

    # Determine whether person has diagnoses in any relevant groups.
    # Diagnoses from ICD9 and ICD10 cannot be null or missing under current
    # interpretation.
    if (
        (not alcohol_diagnosis_a) and
        (not alcohol_diagnosis_b) and
        (not alcohol_diagnosis_c) and
        (not alcohol_diagnosis_d) and
        (not alcohol_diagnosis_self)
    ):
        match_diagnosis = True
    else:
        match_diagnosis = False
    # Integrate information from both criteria.
    if (match_consumption and match_audit and match_diagnosis):
        # Person qualifies as a control for alcoholism.
        # Person's AUDIT scores are below diagnostic thresholds.
        # Person does not have any ICD9 or ICD10 diagnostic codes
        # indicative of alcoholism.
        match = True
    else:
        match = False

    # Return information.
    return match


def determine_case_control_alcoholism_one(
    alcohol_none=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditc_case=None,
    threshold_auditp_control=None,
    threshold_auditp_case=None,
    threshold_audit_control=None,
    threshold_audit_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    arguments:
        alcohol_none (float): binary representation of whether person has never
            consumed alcohol, formerly or currently
        alcohol_diagnosis_a (bool): whether person has diagnosis in alcoholism
            group A
        alcohol_diagnosis_b (bool): whether person has diagnosis in alcoholism
            group B
        alcohol_diagnosis_c (bool): whether person has diagnosis in alcoholism
            group C
        alcohol_diagnosis_d (bool): whether person has diagnosis in alcoholism
            group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-C score
        threshold_auditc_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-C score
        threshold_auditp_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-P score
        threshold_auditp_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-P score
        threshold_audit_control (float): diagnostic control threshold (less
            than or equal) for AUDIT score
        threshold_audit_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT score

    raises:

    returns:
        (float): binary representation whether person qualifies as a case or
            control by specific definition of alcoholism

    """

    # Determine whether person qualifies as a case of alcoholism.
    #case = bool(alcohol_diagnosis_a)
    if (
        (alcohol_diagnosis_a)
    ):
        case = True
    else:
        case = False

    # Determine whether person qualifies as a control of alcoholism.
    # Person's control status can have null or missing values.
    control = determine_control_alcoholism(
        alcohol_none=alcohol_none,
        alcohol_diagnosis_a=alcohol_diagnosis_a,
        alcohol_diagnosis_b=alcohol_diagnosis_b,
        alcohol_diagnosis_c=alcohol_diagnosis_c,
        alcohol_diagnosis_d=alcohol_diagnosis_d,
        alcohol_diagnosis_self=alcohol_diagnosis_self,
        alcohol_auditc=alcohol_auditc,
        alcohol_auditp=alcohol_auditp,
        alcohol_audit=alcohol_audit,
        threshold_auditc_control=threshold_auditc_control,
        threshold_auditc_case=threshold_auditc_case,
        threshold_auditp_control=threshold_auditp_control,
        threshold_auditp_case=threshold_auditp_case,
        threshold_audit_control=threshold_audit_control,
        threshold_audit_case=threshold_audit_case,
    )

    # Interpret case and control and assign value.
    # Interpretation variables "case" and "control" do not have null values.
    # Assign missing value to persons who qualify neither as case nor control.
    if case:
        value = True
    elif control:
        value = False
    else:
        value = None
    # Return information.
    return value


def determine_case_control_alcoholism_two(
    alcohol_none=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditc_case=None,
    threshold_auditp_control=None,
    threshold_auditp_case=None,
    threshold_audit_control=None,
    threshold_audit_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    arguments:
        alcohol_none (float): binary representation of whether person has never
            consumed alcohol, formerly or currently
        alcohol_diagnosis_a (bool): whether person has diagnosis in alcoholism
            group A
        alcohol_diagnosis_b (bool): whether person has diagnosis in alcoholism
            group B
        alcohol_diagnosis_c (bool): whether person has diagnosis in alcoholism
            group C
        alcohol_diagnosis_d (bool): whether person has diagnosis in alcoholism
            group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-C score
        threshold_auditc_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-C score
        threshold_auditp_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-P score
        threshold_auditp_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-P score
        threshold_audit_control (float): diagnostic control threshold (less
            than or equal) for AUDIT score
        threshold_audit_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT score

    raises:

    returns:
        (float): binary representation whether person qualifies as a case or
            control by specific definition of alcoholism

    """

    # Determine whether person qualifies as a case of alcoholism.
    if (
        (alcohol_diagnosis_a) or
        (alcohol_diagnosis_b) or
        (alcohol_diagnosis_c) or
        (alcohol_diagnosis_d)
    ):
        case = True
    else:
        case = False

    # Determine whether person qualifies as a control of alcoholism.
    # Person's control status can have null or missing values.
    control = determine_control_alcoholism(
        alcohol_none=alcohol_none,
        alcohol_diagnosis_a=alcohol_diagnosis_a,
        alcohol_diagnosis_b=alcohol_diagnosis_b,
        alcohol_diagnosis_c=alcohol_diagnosis_c,
        alcohol_diagnosis_d=alcohol_diagnosis_d,
        alcohol_diagnosis_self=alcohol_diagnosis_self,
        alcohol_auditc=alcohol_auditc,
        alcohol_auditp=alcohol_auditp,
        alcohol_audit=alcohol_audit,
        threshold_auditc_control=threshold_auditc_control,
        threshold_auditc_case=threshold_auditc_case,
        threshold_auditp_control=threshold_auditp_control,
        threshold_auditp_case=threshold_auditp_case,
        threshold_audit_control=threshold_audit_control,
        threshold_audit_case=threshold_audit_case,
    )

    # Interpret case and control and assign value.
    # Interpretation variables "case" and "control" do not have null values.
    # Assign missing value to persons who qualify neither as case nor control.
    if case:
        value = True
    elif control:
        value = False
    else:
        value = None
    # Return information.
    return value


def determine_case_control_alcoholism_three(
    alcohol_none=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditc_case=None,
    threshold_auditp_control=None,
    threshold_auditp_case=None,
    threshold_audit_control=None,
    threshold_audit_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    arguments:
        alcohol_none (float): binary representation of whether person has never
            consumed alcohol, formerly or currently
        alcohol_diagnosis_a (bool): whether person has diagnosis in alcoholism
            group A
        alcohol_diagnosis_b (bool): whether person has diagnosis in alcoholism
            group B
        alcohol_diagnosis_c (bool): whether person has diagnosis in alcoholism
            group C
        alcohol_diagnosis_d (bool): whether person has diagnosis in alcoholism
            group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-C score
        threshold_auditc_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-C score
        threshold_auditp_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-P score
        threshold_auditp_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-P score
        threshold_audit_control (float): diagnostic control threshold (less
            than or equal) for AUDIT score
        threshold_audit_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT score

    raises:

    returns:
        (float): binary representation whether person qualifies as a case or
            control by specific definition of alcoholism

    """

    # Determine whether person qualifies as a case of alcoholism.
    if (math.isnan(alcohol_auditc)):
        case = False
    else:
        if (alcohol_auditc >= threshold_auditc_case):
            case = True
        else:
            case = False

    # Determine whether person qualifies as a control of alcoholism.
    # Person's control status can have null or missing values.
    control = determine_control_alcoholism(
        alcohol_none=alcohol_none,
        alcohol_diagnosis_a=alcohol_diagnosis_a,
        alcohol_diagnosis_b=alcohol_diagnosis_b,
        alcohol_diagnosis_c=alcohol_diagnosis_c,
        alcohol_diagnosis_d=alcohol_diagnosis_d,
        alcohol_diagnosis_self=alcohol_diagnosis_self,
        alcohol_auditc=alcohol_auditc,
        alcohol_auditp=alcohol_auditp,
        alcohol_audit=alcohol_audit,
        threshold_auditc_control=threshold_auditc_control,
        threshold_auditc_case=threshold_auditc_case,
        threshold_auditp_control=threshold_auditp_control,
        threshold_auditp_case=threshold_auditp_case,
        threshold_audit_control=threshold_audit_control,
        threshold_audit_case=threshold_audit_case,
    )

    # Interpret case and control and assign value.
    # Interpretation variables "case" and "control" do not have null values.
    # Assign missing value to persons who qualify neither as case nor control.
    if case:
        value = True
    elif control:
        value = False
    else:
        value = None
    # Return information.
    return value


def determine_case_control_alcoholism_four(
    alcohol_none=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditc_case=None,
    threshold_auditp_control=None,
    threshold_auditp_case=None,
    threshold_audit_control=None,
    threshold_audit_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    arguments:
        alcohol_none (float): binary representation of whether person has never
            consumed alcohol, formerly or currently
        alcohol_diagnosis_a (bool): whether person has diagnosis in alcoholism
            group A
        alcohol_diagnosis_b (bool): whether person has diagnosis in alcoholism
            group B
        alcohol_diagnosis_c (bool): whether person has diagnosis in alcoholism
            group C
        alcohol_diagnosis_d (bool): whether person has diagnosis in alcoholism
            group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-C score
        threshold_auditc_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-C score
        threshold_auditp_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-P score
        threshold_auditp_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-P score
        threshold_audit_control (float): diagnostic control threshold (less
            than or equal) for AUDIT score
        threshold_audit_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT score

    raises:

    returns:
        (float): binary representation whether person qualifies as a case or
            control by specific definition of alcoholism

    """

    # Determine whether person qualifies as a case of alcoholism.
    if (math.isnan(alcohol_auditp)):
        case = False
    else:
        if (alcohol_auditp >= threshold_auditp_case):
            case = True
        else:
            case = False
    # Determine whether person qualifies as a control of alcoholism.
    # Person's control status can have null or missing values.
    control = determine_control_alcoholism(
        alcohol_none=alcohol_none,
        alcohol_diagnosis_a=alcohol_diagnosis_a,
        alcohol_diagnosis_b=alcohol_diagnosis_b,
        alcohol_diagnosis_c=alcohol_diagnosis_c,
        alcohol_diagnosis_d=alcohol_diagnosis_d,
        alcohol_diagnosis_self=alcohol_diagnosis_self,
        alcohol_auditc=alcohol_auditc,
        alcohol_auditp=alcohol_auditp,
        alcohol_audit=alcohol_audit,
        threshold_auditc_control=threshold_auditc_control,
        threshold_auditc_case=threshold_auditc_case,
        threshold_auditp_control=threshold_auditp_control,
        threshold_auditp_case=threshold_auditp_case,
        threshold_audit_control=threshold_audit_control,
        threshold_audit_case=threshold_audit_case,
    )

    # Interpret case and control and assign value.
    # Interpretation variables "case" and "control" do not have null values.
    # Assign missing value to persons who qualify neither as case nor control.
    if case:
        value = True
    elif control:
        value = False
    else:
        value = None
    # Return information.
    return value


def determine_case_control_alcoholism_five(
    alcohol_none=None,
    alcohol_diagnosis_a=None,
    alcohol_diagnosis_b=None,
    alcohol_diagnosis_c=None,
    alcohol_diagnosis_d=None,
    alcohol_diagnosis_self=None,
    alcohol_auditc=None,
    alcohol_auditp=None,
    alcohol_audit=None,
    threshold_auditc_control=None,
    threshold_auditc_case=None,
    threshold_auditp_control=None,
    threshold_auditp_case=None,
    threshold_audit_control=None,
    threshold_audit_case=None,
):
    """
    Organizes information about alcoholism cases and controls.

    arguments:
        alcohol_none (float): binary representation of whether person has never
            consumed alcohol, formerly or currently
        alcohol_diagnosis_a (bool): whether person has diagnosis in alcoholism
            group A
        alcohol_diagnosis_b (bool): whether person has diagnosis in alcoholism
            group B
        alcohol_diagnosis_c (bool): whether person has diagnosis in alcoholism
            group C
        alcohol_diagnosis_d (bool): whether person has diagnosis in alcoholism
            group D
        alcohol_diagnosis_self (bool): whether person has self diagnosis of
            alcoholism
        alcohol_auditc (float): combination score for AUDIT-C questionnaire
        alcohol_auditp (float): combination score for AUDIT-P questionnaire
        alcohol_audit (float): combination score for AUDIT questionnaire
        threshold_auditc_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-C score
        threshold_auditc_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-C score
        threshold_auditp_control (float): diagnostic control threshold (less
            than or equal) for AUDIT-P score
        threshold_auditp_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT-P score
        threshold_audit_control (float): diagnostic control threshold (less
            than or equal) for AUDIT score
        threshold_audit_case (float): diagnostic case threshold (greater than
            or equal) for AUDIT score

    raises:

    returns:
        (float): binary representation whether person qualifies as a case or
            control by specific definition of alcoholism

    """

    # Determine whether person qualifies as a case of alcoholism.
    if (math.isnan(alcohol_audit)):
        case = False
    else:
        if (alcohol_audit >= threshold_audit_case):
            case = True
        else:
            case = False
    # Determine whether person qualifies as a control of alcoholism.
    # Person's control status can have null or missing values.
    control = determine_control_alcoholism(
        alcohol_none=alcohol_none,
        alcohol_diagnosis_a=alcohol_diagnosis_a,
        alcohol_diagnosis_b=alcohol_diagnosis_b,
        alcohol_diagnosis_c=alcohol_diagnosis_c,
        alcohol_diagnosis_d=alcohol_diagnosis_d,
        alcohol_diagnosis_self=alcohol_diagnosis_self,
        alcohol_auditc=alcohol_auditc,
        alcohol_auditp=alcohol_auditp,
        alcohol_audit=alcohol_audit,
        threshold_auditc_control=threshold_auditc_control,
        threshold_auditc_case=threshold_auditc_case,
        threshold_auditp_control=threshold_auditp_control,
        threshold_auditp_case=threshold_auditp_case,
        threshold_audit_control=threshold_audit_control,
        threshold_audit_case=threshold_audit_case,
    )

    # Interpret case and control and assign value.
    # Interpretation variables "case" and "control" do not have null values.
    # Assign missing value to persons who qualify neither as case nor control.
    if case:
        value = True
    elif control:
        value = False
    else:
        value = None
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
                        "sex", "sex_text", "age", "body_mass_index",
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
                            "sex", "sex_text", "age", "body_mass_index",
                            "alcohol_none",
                            alcoholism,
                            hormone,
                        ],
                        variables_prefixes_valid=["genotype_pc_",],
                        table=table,
                    )
    pass


def organize_alcoholism_cases_controls_variables(
    table=None,
    report=None,
):
    """
    Organizes information about alcoholism cases and controls.

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
    # Use less than or equal for control thresholds.
    # Use greater than or equal for case thresholds.
    threshold_auditc_control = 4
    threshold_auditc_case = 10
    threshold_auditp_control = 3
    threshold_auditp_case = 5
    threshold_audit_control = 7
    threshold_audit_case = 15

    # Determine whether person is a case or control for alcoholism type 1.
    # case: ICD9 or ICD10 codes in diagnostic group A
    # control:
    # - only persons who have consumed alcohol, previously or currently
    # - no ICD9 or ICD10 codes in diagnostic groups A, B, C, or D
    # - no self diagnoses of alcoholism
    # - below AUDIT control thresholds
    table["alcoholism_1"] = table.apply(
        lambda row:
            determine_case_control_alcoholism_one(
                alcohol_none=row["alcohol_none"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
                alcohol_audit=row["alcohol_audit"],
                threshold_auditc_control=threshold_auditc_control,
                threshold_auditc_case=threshold_auditc_case,
                threshold_auditp_control=threshold_auditp_control,
                threshold_auditp_case=threshold_auditp_case,
                threshold_audit_control=threshold_audit_control,
                threshold_audit_case=threshold_audit_case,
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person is a case or control for alcoholism type 2.
    # case: ICD9 or ICD10 codes in diagnostic group A, B, C, or D
    # control:
    # - only persons who have consumed alcohol, previously or currently
    # - no ICD9 or ICD10 codes in diagnostic groups A, B, C, or D
    # - no self diagnoses of alcoholism
    # - below AUDIT control thresholds
    table["alcoholism_2"] = table.apply(
        lambda row:
            determine_case_control_alcoholism_two(
                alcohol_none=row["alcohol_none"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
                alcohol_audit=row["alcohol_audit"],
                threshold_auditc_control=threshold_auditc_control,
                threshold_auditc_case=threshold_auditc_case,
                threshold_auditp_control=threshold_auditp_control,
                threshold_auditp_case=threshold_auditp_case,
                threshold_audit_control=threshold_audit_control,
                threshold_audit_case=threshold_audit_case,
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person is a case or control for alcoholism type 3.
    # case:
    # AUDIT-C score > threshold_auditc
    # control:
    # - only persons who have consumed alcohol, previously or currently
    # - no ICD9 or ICD10 codes in diagnostic groups A, B, C, or D
    # - no self diagnoses of alcoholism
    # - below AUDIT control thresholds
    table["alcoholism_3"] = table.apply(
        lambda row:
            determine_case_control_alcoholism_three(
                alcohol_none=row["alcohol_none"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
                alcohol_audit=row["alcohol_audit"],
                threshold_auditc_control=threshold_auditc_control,
                threshold_auditc_case=threshold_auditc_case,
                threshold_auditp_control=threshold_auditp_control,
                threshold_auditp_case=threshold_auditp_case,
                threshold_audit_control=threshold_audit_control,
                threshold_audit_case=threshold_audit_case,
            ),
        axis="columns", # apply across rows
    )

    # Determine whether person is a case or control for alcoholism type 4.
    # case:
    # AUDIT-P score > threshold_auditp_case
    # control:
    # - only persons who have consumed alcohol, previously or currently
    # - no ICD9 or ICD10 codes in diagnostic groups A, B, C, or D
    # - no self diagnoses of alcoholism
    # - below AUDIT control thresholds
    table["alcoholism_4"] = table.apply(
        lambda row:
            determine_case_control_alcoholism_four(
                alcohol_none=row["alcohol_none"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
                alcohol_audit=row["alcohol_audit"],
                threshold_auditc_control=threshold_auditc_control,
                threshold_auditc_case=threshold_auditc_case,
                threshold_auditp_control=threshold_auditp_control,
                threshold_auditp_case=threshold_auditp_case,
                threshold_audit_control=threshold_audit_control,
                threshold_audit_case=threshold_audit_case,
            ),
        axis="columns", # apply across rows
    )

    # Determine whether person is a case or control for alcoholism type 5.
    # case:
    # AUDIT score > threshold_audit
    # control:
    # - only persons who have consumed alcohol, previously or currently
    # - no ICD9 or ICD10 codes in diagnostic groups A, B, C, or D
    # - no self diagnoses of alcoholism
    # - below AUDIT control thresholds
    table["alcoholism_5"] = table.apply(
        lambda row:
            determine_case_control_alcoholism_five(
                alcohol_none=row["alcohol_none"],
                alcohol_diagnosis_a=row["alcohol_diagnosis_a"],
                alcohol_diagnosis_b=row["alcohol_diagnosis_b"],
                alcohol_diagnosis_c=row["alcohol_diagnosis_c"],
                alcohol_diagnosis_d=row["alcohol_diagnosis_d"],
                alcohol_diagnosis_self=row["alcohol_diagnosis_self"],
                alcohol_auditc=row["alcohol_auditc"],
                alcohol_auditp=row["alcohol_auditp"],
                alcohol_audit=row["alcohol_audit"],
                threshold_auditc_control=threshold_auditc_control,
                threshold_auditc_case=threshold_auditc_case,
                threshold_auditp_control=threshold_auditp_control,
                threshold_auditp_case=threshold_auditp_case,
                threshold_audit_control=threshold_audit_control,
                threshold_audit_case=threshold_audit_case,
            ),
        axis="columns", # apply across rows
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
            "alcohol_none",
            "alcohol_diagnosis_a",
            "alcohol_diagnosis_b",
            "alcohol_diagnosis_c",
            "alcohol_diagnosis_d",
            "alcohol_diagnosis_self",
            "alcohol_auditc",
            "alcohol_audit",
            "alcoholism_1", "alcoholism_2", "alcoholism_3", "alcoholism_4",
            "alcoholism_5",
        ])
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Summary of alcoholism cases and controls variables: ")
        print(table_report)
        organize_report_cohorts_by_sex_alcoholism_split_hormone(
            table=table_clean,
        )
        pass

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
# ... in progress...

# Scrap now...


def translate_alcohol_none_auditc(
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


def determine_auditc_questionnaire_alcoholism_score(
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


def organize_auditc_questionnaire_alcoholism_variables(
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
        axis="columns", # apply across rows
    )
    # Determine alcoholism aggregate score.
    table["alcoholism"] = table.apply(
        lambda row:
            determine_auditc_questionnaire_alcoholism_score(
                frequency=row["20414-0.0"],
                quantity=row["20403-0.0"],
                binge=row["20416-0.0"],
            ),
        axis="columns", # apply across rows
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

    # Organization procedure main information.
    write_product_organization(
        information=information["organization"],
        path_parent=paths["organization"],
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
            "sex", "sex_text", "age", "body_mass_index", "body_mass_index_log",
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
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Organize information about sex hormones.
    pail_hormone = organize_sex_hormone_variables(
        table=table,
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


def execute_psychology_psychiatry(
    table=None,
    report=None,
):
    """
    Organizes information about persons' mental health across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Organize information about neuroticism.
    pail_psychology = organize_psychology_variables(
        table=table,
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


def execute_alcohol(
    table=None,
    report=None,
):
    """
    Organizes information about persons' alcohol consumption and dependence
    across UK Biobank.

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
        report=False,
    )
    #print(pail_alcohol_consumption["quantity"]["table_clean"])

    if False:
        # Organize Alchol Use Disorders Identification Test (AUDIT) questionnaire
        # variables, including separate scores for AUDIT-Consumption (AUDIT-C) and
        # AUDIT-Problem (AUDIT-P) portions of the questionnaire.
        pail_audit = organize_alcohol_audit_questionnaire_variables(
            table=pail_alcohol_consumption["quantity"]["table_clean"],
            report=False,
        )
        #print(pail_audit["audit"]["table_clean"])

        # Organize International Classification of Disease (ICD) ICD9 and ICD10
        # codes for diagnoses relevant to alcoholism.
        # Organize codes for self diagnoses relevant to alcoholism.
        pail_diagnosis = organize_alcoholism_diagnosis_variables(
            table=pail_audit["audit"]["table_clean"],
            report=False,
        )
        #print(pail_diagnosis["table_clean"])

        # Organize alcoholism cases and controls.
        # Report females and males who consume alcohol and are candidates for
        # either controls or cases of alcoholism.
        pail_alcoholism = organize_alcoholism_cases_controls_variables(
            table=pail_diagnosis["table_clean"],
            report=True,
        )
        #print(pail_alcoholism["table_clean"])

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: execute_alcohol()")
        utility.print_terminal_partition(level=3)
        print(pail_alcohol_consumption["table_report"])
    # Return information.
    return pail_alcohol_consumption




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
    print("version check: 1")
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
    # Organize variables for persons' genotypes, sex, age, and body mass index
    # across the UK Biobank.
    pail_basis = execute_genotype_assessment_basis(
        table=source["table_phenotypes"],
        path_dock=path_dock,
        report=True,
    )
    # Organize variables for female menstruation across the UK Biobank.
    pail_female = execute_female_menstruation(
        table=pail_basis["table"], # pail_hormone["table_clean"]
        report=True,
    )
    # Organize variables for persons' sex hormones across the UK Biobank.
    pail_hormone = execute_sex_hormones(
        table=pail_female["table"], # pail_basis["table_clean"]
        report=True,
    )

    # Collect information.
    information = dict()
    information["organization"] = dict()
    #information["organization"]["table_phenotypes"] = pail_basis["table"]
    information["organization"]["table_phenotypes"] = pail_hormone["table"]
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )
    pass
