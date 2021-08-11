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

###############################################################################
# Functionality

##########
# Read


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
    # Return information.
    return table_kinship_pairs


##########
# Cohort, model selection: kinship


def filter_kinship_pairs_by_threshold_relevance(
    name=None,
    threshold_kinship=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Selects records by sex and by sex-specific criteria and variables.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        threshold_kinship (float): maximal value of kinship coefficient
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        table (object): Pandas data frame of values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of kinship coefficients
            across pairs of persons

    """

    # Copy information.
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)
    table = table.copy(deep=True)
    # Filter kinship pairs to exclude any redundant pairs.
    table_kinship_pairs.drop_duplicates(
        subset=["ID1", "ID2",],
        keep="first",
        inplace=True,
    )
    # Filter kinship pairs to exclude persons with kinship below threshold.
    # These pairs do not have strong enough relation for further consideration.
    count_pairs_original = copy.deepcopy(table_kinship_pairs.shape[0])
    table_kinship_pairs = table_kinship_pairs.loc[
        (table_kinship_pairs["Kinship"] < threshold_kinship), :
    ]
    count_pairs_kinship = copy.deepcopy(table_kinship_pairs.shape[0])
    # Filter kinship pairs to exclude persons who are not in the
    # analysis-specific cohort table.
    # Consider both persons from each kinship pair.
    genotypes_relevant = copy.deepcopy(table["IID"].to_list())
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
    table_kinship_pairs = table_kinship_pairs.loc[
        table_kinship_pairs.index.isin(genotypes_relevant), :
    ]
    table_kinship_pairs.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_kinship_pairs.set_index(
        "ID2",
        append=False,
        drop=True,
        inplace=True
    )
    table_kinship_pairs = table_kinship_pairs.loc[
        table_kinship_pairs.index.isin(genotypes_relevant), :
    ]
    count_pairs_relevant = copy.deepcopy(table_kinship_pairs.shape[0])
    table_kinship_pairs.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Report.
    if report:
        # Determine count of persons relevant to the analysis-specific cohort who
        # belong to kinship pairs beyond threshold.
        genotypes_kinship = copy.deepcopy(table_kinship_pairs["ID1"].to_list())
        genotypes_kinship_second = copy.deepcopy(
            table_kinship_pairs["ID2"].to_list()
        )
        genotypes_kinship.extend(genotypes_kinship_second)
        genotypes_kinship_unique = list(set(genotypes_kinship))
        genotypes_common = utility.filter_common_elements(
            list_minor=genotypes_kinship_unique,
            list_major=genotypes_relevant,
        )
        count_genotypes_relevant = len(genotypes_relevant)
        count_genotypes_kinship = len(genotypes_common)
        percentage = round(
            ((count_genotypes_kinship / count_genotypes_relevant) * 100), 2
        )
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print(
            "report: " +
            "filter_kinship_pairs_by_threshold_relevance()"
        )
        utility.print_terminal_partition(level=5)
        print("cohort: " + str(name))
        utility.print_terminal_partition(level=5)
        print("... kinship pairs ...")
        print("original: " + str(count_pairs_original))
        print("kinship threshold: " + str(count_pairs_kinship))
        print("relevant to analysis cohort: " + str(count_pairs_relevant))
        print("..........")
        print("... genotypes (persons) ...")
        print(
            "total relevant genotypes in cohort: " +
            str(count_genotypes_relevant)
        )
        print(
            "relevant genotypes with kinship beyond threshold: " +
            str(count_genotypes_kinship) + " (" + str(percentage) + "%)"
        )
        utility.print_terminal_partition(level=5)
    # Return information.
    return table_kinship_pairs


def filter_persons_ukbiobank_by_kinship(
    name=None,
    threshold_kinship=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Selects records by sex and by sex-specific criteria and variables.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
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

    # Copy information.
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)
    table = table.copy(deep=True)
    # Filter kinship pairs to exclude persons with kinship below threshold and
    # those who are not in the analysis-specific cohort table.
    table_kinship_pairs = filter_kinship_pairs_by_threshold_relevance(
        name=name,
        threshold_kinship=threshold_kinship,
        table_kinship_pairs=table_kinship_pairs,
        table=table,
        report=report, # report procedure is quite slow
    )
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
    representative_nodes = list()
    for component in networkx.connected_components(network):
        network_component = network.subgraph(component).copy()
        representative = random.choice(list(network_component.nodes()))
        representative_nodes.append(representative)
    # Convert genotype identifiers to strings.
    representative_nodes = list(map(str, representative_nodes))
    count_kinship_representatives = len(representative_nodes)
    # Exclude all related persons from the cohort, with the exception of a
    # single representative from each family (connected component of the kinship
    # pair network).
    genotypes_kinship = copy.deepcopy(table_kinship_pairs["ID1"].to_list())
    genotypes_kinship_second = copy.deepcopy(
        table_kinship_pairs["ID2"].to_list()
    )
    genotypes_kinship.extend(genotypes_kinship_second)
    genotypes_kinship_unique = list(set(genotypes_kinship))
    genotypes_kinship_unique = list(map(str, genotypes_kinship_unique))
    count_kinship_unique = len(genotypes_kinship_unique)
    genotypes_kinship_exclusion = list(filter(
        lambda value: (value not in representative_nodes),
        genotypes_kinship_unique
    ))
    count_persons_original = copy.deepcopy(table.shape[0])
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table.set_index(
        "IID",
        append=False,
        drop=True,
        inplace=True
    )
    table = table.loc[
        ~table.index.isin(genotypes_kinship_exclusion), :
    ]
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    count_persons_novel = copy.deepcopy(table.shape[0])

    # Report.
    if report:
        # Calculate percentage of persons lost to kinship threshold.
        count_loss_kinship = (
            count_kinship_unique - count_kinship_representatives
        )
        percentage_loss_kinship = round(
            ((count_loss_kinship / count_persons_original) * 100), 2
        )
        count_loss_check = (
            count_persons_original - count_persons_novel
        )
        percentage_loss_check = round(
            ((count_loss_check / count_persons_original) * 100), 2
        )
        # Report.
        utility.print_terminal_partition(level=2)
        print(
            "report: " +
            "filter_persons_ukbiobank_by_kinship()"
        )
        utility.print_terminal_partition(level=5)
        print("cohort: " + str(name))
        utility.print_terminal_partition(level=5)
        print("kinship network graph")
        print("network order: " + str(network.order()))
        print("network size: " + str(network.size()))
        print(
            "count connected components: " +
            str(networkx.number_connected_components(network))
        )
        print(
            "count representative nodes: " + str(count_kinship_representatives)
        )
        print(
            "cohort loss to kinship: " + str(count_loss_kinship) + " of " +
            str(count_persons_original) + " (" + str(percentage_loss_kinship) +
            "%)"
        )
        print(
            "check cohort loss to kinship: " + str(count_loss_check) + " of " +
            str(count_persons_original) + " (" + str(percentage_loss_check) +
            "%)"
        )
    # Return information.
    return table


##########
# Cohort, model selection: general


def select_valid_records_all_specific_variables(
    names=None,
    prefixes=None,
    table=None,
    drop_columns=None,
    report=None,
):
    """
    Selects variable columns and record rows with valid values across all
    variables.

    arguments:
        names (list<str>): explicit names of columns to keep
        prefixes (list<str>): prefixes of names of columns to keep
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        drop_columns (bool): whether to drop other columns
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

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
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Columns to keep: ")
        print(columns_keep)
    # Select columns.
    if drop_columns:
        table = table.loc[
            :, table.columns.isin(columns_keep)
        ]
    # Remove any record rows with null values.
    table.dropna(
        axis="index",
        how="any",
        subset=columns_keep,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("After dropping rows with null values in specific columns.")
        print("Selection of table columns and rows: ")
        print(table)
    # Return information.
    return table


def select_records_by_female_specific_valid_variables_values(
    pregnancy=None,
    menopause_binary=None,
    menopause_ordinal=None,
    variables=None,
    prefixes=None,
    table=None,
):
    """
    Selects records for females with sex-specific criteria.

    When filtering either by "menopause_binary" or "menopause_ordinal", be sure
    to require that these variables have valid, non-null values.

    arguments:
        pregnancy (list<int>): which values of pregnancy definition to include
            for females
        menopause_binary (list<int>): which values of binary menopause
            definition to include for females
        menopause_ordinal (list<int>): which values of ordinal menopause
            definition to include for females
        variables (list<str>): names of columns for variables in which
            rows must have valid values to keep records
        prefixes (list<str>): prefixes of columns for variables in which
            rows must have valid values to keep records
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Select records with valid (non-null) values of relevant variables.
    # Exclude missing values first to avoid interpretation of "None" as False.
    table = select_valid_records_all_specific_variables(
        names=variables,
        prefixes=prefixes,
        table=table,
        drop_columns=True,
        report=False,
    )
    # Select records for females.
    table = table.loc[
        (table["sex_text"] == "female"), :
    ]
    # Determine whether to filter by pregnancy.
    if (
        (0 not in pregnancy) or
        (1 not in pregnancy)
    ):
        # Select records.
        table = table.loc[
            (table["pregnancy"].isin(pregnancy)), :
        ]
    # Determine whether to filter by binary definition of menopause.
    if (
        (0 not in menopause_binary) or
        (1 not in menopause_binary)
    ):
        # Select records.
        table = table.loc[
            (table["menopause_binary"].isin(menopause_binary)), :
        ]
    # Determine whether to filter by ordinal definition of menopause.
    if (
        (0 not in menopause_ordinal) or
        (1 not in menopause_ordinal) or
        (2 not in menopause_ordinal)
    ):
        # Select records.
        table = table.loc[
            (table["menopause_ordinal"].isin(menopause_ordinal)), :
        ]
    # Return information.
    return table


def select_records_by_male_specific_valid_variables_values(
    age_grade_male=None,
    variables=None,
    prefixes=None,
    table=None,
):
    """
    Selects records for males with sex-specific criteria.

    arguments:
        age_grade_male (list<int>): which values of ordinal age variable to
            include for males
        variables (list<str>): names of columns for variables in which
            rows must have valid values to keep records
        prefixes (list<str>): prefixes of columns for variables in which
            rows must have valid values to keep records
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Select records with valid (non-null) values of relevant variables.
    # Exclude missing values first to avoid interpretation of "None" as False.
    table = select_valid_records_all_specific_variables(
        names=variables,
        prefixes=prefixes,
        table=table,
        drop_columns=True,
        report=False,
    )
    # Select records for males.
    table = table.loc[
        (table["sex_text"] == "male"), :
    ]
    # Determine whether to filter by ordinal definition of age.
    if (
        (0 not in age_grade_male) or
        (1 not in age_grade_male) or
        (2 not in age_grade_male)
    ):
        # Select records.
        table = table.loc[
            (table["age_grade_male"].isin(age_grade_male)), :
        ]
    # Return information.
    return table


def select_records_by_ancestry_sex_specific_valid_variables_values(
    name=None,
    white_british=None,
    female=None,
    female_pregnancy=None,
    female_menopause_binary=None,
    female_menopause_ordinal=None,
    female_variables=None,
    female_prefixes=None,
    male=None,
    age_grade_male=None,
    male_variables=None,
    male_prefixes=None,
    table_kinship_pairs=None,
    table=None,
):
    """
    Selects records by sex and by sex-specific criteria and variables.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        white_british (list<int>): which values of white british categorical
            ancestry variable to include
        female (bool): whether to include records for females
        female_pregnancy (list<int>): which values of pregnancy definition to
            include for females
        female_menopause_binary (list<int>): which values of binary menopause
            definition to include for females
        female_menopause_ordinal (list<int>): which values of ordinal menopause
            definition to include for females
        female_variables (list<str>): names of columns for variables in which
            rows must have valid values to keep records for females
        female_prefixes (list<str>): prefixes of columns for variables in which
            rows must have valid values to keep records for females
        male (bool): whether to include records for males
        age_grade_male (list<int>): which values of ordinal age variable to
            include for males
        male_variables (list<str>): names of columns for variables in which rows
            must have valid values to keep records for males
        male_prefixes (list<str>): prefixes of columns for variables in which
            rows must have valid values to keep records for males
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        table (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect records.
    table_collection = pandas.DataFrame()
    # Select records for females.
    if female:
        table_female = select_records_by_female_specific_valid_variables_values(
            pregnancy=female_pregnancy,
            menopause_binary=female_menopause_binary,
            menopause_ordinal=female_menopause_ordinal,
            variables=female_variables,
            prefixes=female_prefixes,
            table=table,
        )
        table_collection = table_collection.append(
            table_female,
            ignore_index=True,
        )
        pass
    # Select records for males.
    if male:
        table_male = select_records_by_male_specific_valid_variables_values(
            age_grade_male=age_grade_male,
            variables=male_variables,
            prefixes=male_prefixes,
            table=table,
        )
        table_collection = table_collection.append(
            table_male,
            ignore_index=True,
        )
        pass
    # Determine whether to filter by white british categorical ancestry.
    if (
        (0 not in white_british) or
        (1 not in white_british)
    ):
        # Select records.
        table_collection = table_collection.loc[
            (table_collection["white_british"].isin(white_british)), :
        ]
    # Filter by kinship relatedness.
    table_unrelated = filter_persons_ukbiobank_by_kinship(
        name=name,
        threshold_kinship=0.1, # pairs with kinship >= threshold for exclusion
        table_kinship_pairs=table_kinship_pairs,
        table=table_collection,
        report=False,
    )
    # Organize table.
    table_unrelated.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_unrelated.set_index(
        "eid",
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return table_unrelated


def translate_boolean_phenotype_plink(
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
        value = float("nan")
    # Return information.
    return value


def translate_binary_phenotype_plink(
    binary_value=None,
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
        binary_value (float): binary (0, 1) representation of a phenotype

    raises:

    returns:
        (float): plink binary representation of a phenotype

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not math.isnan(binary_value)) and
        (-0.5 <= binary_value and binary_value < 1.5)
    ):
        # The variable has a valid value.
        if (-0.5 <= binary_value and binary_value < 0.5):
            # control, original value: 0
            # control, Plink translation: 1
            value = 1
        elif (0.5 <= binary_value and binary_value < 1.5):
            # case, original value: 1
            # case, Plink translation: 2
            value = 2
    else:
        # null
        value = float("nan")
    # Return information.
    return value


def organize_phenotype_covariate_table_plink_format(
    boolean_phenotypes=None,
    binary_phenotypes=None,
    continuous_variables=None,
    remove_null_records=None,
    table=None,
):
    """
    Organize table for phenotypes and covariates in format for PLINK.

    1. Remove any rows with missing, empty values.
    PLINK cannot accommodate rows with empty cells.

    2. Introduce family identifiers.
    Family (FID) and individual (IID) identifiers must match the ID_1 and ID_2
    columns in the sample table.

    3. Sort column sequence.
    PLINK requires FID and IID columns to come first.

    arguments:
        boolean_phenotypes (list<str>): names of columns with discrete boolean
            case-control phenotypes (control: False, case: True) that need
            conversion to PLINK format (control: 1, case: 2)
        binary_phenotypes (list<str>): names of columns with discrete binary
            case-control phenotypes (control: 0, case: 1) that need conversion
            to PLINK format (control: 1, case: 2)
        continuous_variables (list<str>): names of columns that PLINK ought
            to interpret as continuous variables (ordinal discrete, interval
            continuous, or ratio continuous)
        remove_null_records (bool): whether to remove records with any null
            values
        table (object): Pandas data frame of information about phenotype and
            covariate variables for GWAS

    raises:

    returns:
        (object): Pandas data frame of information about phenotype and
            covariate variables in format for GWAS in PLINK

    """

    # Copy data.
    table = table.copy(deep=True)
    # Translate binary phenotype variables.
    for boolean_phenotype in boolean_phenotypes:
        table[boolean_phenotype] = table.apply(
            lambda row:
                translate_boolean_phenotype_plink(
                    boolean_value=row[boolean_phenotype],
                ),
            axis="columns", # apply across rows
        )
        pass
    # Translate binary phenotype variables.
    for binary_phenotype in binary_phenotypes:
        table[binary_phenotype] = table.apply(
            lambda row:
                translate_binary_phenotype_plink(
                    binary_value=row[binary_phenotype],
                ),
            axis="columns", # apply across rows
        )
        pass
    # Translate continuous variables, whether they are phenotypes or covariates.
    table_type = convert_table_columns_variables_types_float(
        columns=continuous_variables,
        table=table,
    )
    # Organize.
    table.reset_index(
        level=None,
        inplace=True
    )
    # Remove table rows with any empty cells or missing values.
    if remove_null_records:
        table.dropna(
            axis="index",
            how="any",
            subset=None,
            inplace=True,
        )
    # Introduce family identifier.
    table["FID"] = table["IID"]
    # Sort column sequence.
    columns = table.columns.to_list()
    columns_sequence = list(filter(
        lambda element: element not in ["eid", "IID", "FID"],
        columns
    ))
    columns_sequence.insert(0, "eid") # third column
    columns_sequence.insert(0, "IID") # second column
    columns_sequence.insert(0, "FID") # first column
    table_columns = table.loc[
        :, table.columns.isin(columns_sequence)
    ]
    table_sequence = table_columns[[*columns_sequence]]
    # Return information.
    return table_sequence


def select_records_by_ancestry_case_control_valid_variables_values(
    name=None,
    white_british=None,
    case_control=None,
    case_control_values=None,
    variables=None,
    prefixes=None,
    table_kinship_pairs=None,
    table=None,
):
    """
    Selects records for males with sex-specific criteria.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        white_british (list<int>): which values of white british categorical
            ancestry variable to include
        case_control (str): name of column for relevant case control definition
        case_control_values (list<int>): which values of case control variable
            to include
        variables (list<str>): names of columns for variables in which
            rows must have valid values to keep records
        prefixes (list<str>): prefixes of columns for variables in which
            rows must have valid values to keep records
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    # Select records with valid (non-null) values of relevant variables.
    # Exclude missing values first to avoid interpretation of "None" as False.
    table = select_valid_records_all_specific_variables(
        names=variables,
        prefixes=prefixes,
        table=table,
        drop_columns=True,
        report=False,
    )
    # Determine whether to filter by white british categorical ancestry.
    if (
        (0 not in white_british) or
        (1 not in white_british)
    ):
        # Select records.
        table = table.loc[
            (table["white_british"].isin(white_british)), :
        ]
    # Determine whether to filter by definition of cases and controls.
    if (
        (0 not in case_control_values) or
        (1 not in case_control_values)
    ):
        # Select records.
        table = table.loc[
            (table[case_control].isin(case_control_values)), :
        ]
    # Filter by kinship relatedness.
    table_unrelated = filter_persons_ukbiobank_by_kinship(
        name=name,
        threshold_kinship=0.1, # pairs with kinship >= threshold for exclusion
        table_kinship_pairs=table_kinship_pairs,
        table=table,
        report=False,
    )
    # Organize table.
    table_unrelated.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table_unrelated.set_index(
        "eid",
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return table_unrelated


##########
# Cohort, model selection: sets for genetic analyses on phenotypes and genotypes


def select_organize_cohorts_variables_by_sex_hormone(
    hormone=None,
    table_kinship_pairs=None,
    table=None,
):
    """
    Organizes tables for specific cohorts, phenotypes, and model covariates
    for genetic analysis.

    Notice that different variables are relevant for females and males.
    Combination of tables for females and males introduces missing or null
    values for these sex-specific variables.
    Remove records with null values separately for females and males but not
    after combination.

    arguments:
        hormone (str): name of column for hormone variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: non-pregnant females and males together
    if False:
        record = dict()
        record["name"] = str("table_female_male_" + hormone)
        record["cohort_model"] = "female_male"
        record["category"] = "female_male"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[0, 1,],
                female_menopause_ordinal=[0, 1, 2,],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=True,
                age_grade_male=[0, 1, 2,],
                male_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "age_grade_male",
                    "body_mass_index_log",
                    hormone,
                ],
                male_prefixes=["genotype_pc_",],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: all non-pregnant females together
    if True:
        record = dict()
        record["name"] = str("table_female_" + hormone)
        record["cohort_model"] = "female"
        record["category"] = "female"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[0, 1,],
                female_menopause_ordinal=[0, 1, 2,],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy", "menopause_binary", "menopause_ordinal",
                    "hormone_alteration",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=False,
                age_grade_male=[0, 1, 2,],
                male_variables=[],
                male_prefixes=[],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: premenopausal females by binary menopause definition
    if False:
        record = dict()
        record["name"] = str("table_female_premenopause_binary_" + hormone)
        record["cohort_model"] = "female_premenopause_binary"
        record["category"] = "female_menopause_binary"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[0,],
                female_menopause_ordinal=[0, 1, 2],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy", "menopause_binary", "menopause_ordinal",
                    "menstruation_days",
                    "menstruation_phase", "menstruation_phase_cycle",
                    "hormone_alteration",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=False,
                age_grade_male=[0, 1, 2,],
                male_variables=[],
                male_prefixes=[],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: postmenopausal females by binary menopause definition
    if False:
        record = dict()
        record["name"] = str("table_female_postmenopause_binary_" + hormone)
        record["cohort_model"] = "female_postmenopause_binary"
        record["category"] = "female_menopause_binary"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[1,],
                female_menopause_ordinal=[0, 1, 2,],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy", "menopause_binary", "menopause_ordinal",
                    "hormone_alteration",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=False,
                age_grade_male=[0, 1, 2,],
                male_variables=[],
                male_prefixes=[],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: premenopausal females by ordinal menopause definition
    if True:
        record = dict()
        record["name"] = str("table_female_premenopause_ordinal_" + hormone)
        record["cohort_model"] = "female_premenopause_ordinal"
        record["category"] = "female_menopause_ordinal"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[0, 1,],
                female_menopause_ordinal=[0,],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy", "menopause_binary", "menopause_ordinal",
                    "menstruation_days",
                    "menstruation_phase", "menstruation_phase_cycle",
                    "hormone_alteration",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=False,
                age_grade_male=[0, 1, 2,],
                male_variables=[],
                male_prefixes=[],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: perimenopausal females by ordinal menopause definition
    if True:
        record = dict()
        record["name"] = str("table_female_perimenopause_ordinal_" + hormone)
        record["cohort_model"] = "female_perimenopause_ordinal"
        record["category"] = "female_menopause_ordinal"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[0, 1,],
                female_menopause_ordinal=[1,],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy", "menopause_binary", "menopause_ordinal",
                    "menstruation_days",
                    "menstruation_phase", "menstruation_phase_cycle",
                    "hormone_alteration",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=False,
                age_grade_male=[0, 1, 2,],
                male_variables=[],
                male_prefixes=[],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: postmenopausal females by ordinal menopause definition
    if True:
        record = dict()
        record["name"] = str("table_female_postmenopause_ordinal_" + hormone)
        record["cohort_model"] = "female_postmenopause_ordinal"
        record["category"] = "female_menopause_ordinal"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=True,
                female_pregnancy=[0,],
                female_menopause_binary=[0, 1,],
                female_menopause_ordinal=[2,],
                female_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "body_mass_index_log",
                    "pregnancy", "menopause_binary", "menopause_ordinal",
                    "hormone_alteration",
                    hormone,
                ],
                female_prefixes=["genotype_pc_",],
                male=False,
                age_grade_male=[0, 1, 2,],
                male_variables=[],
                male_prefixes=[],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: males
    if True:
        record = dict()
        record["name"] = str("table_male_" + hormone)
        record["cohort_model"] = "male"
        record["category"] = "male"
        record["phenotype"] = hormone
        record["table"] = (
            select_records_by_ancestry_sex_specific_valid_variables_values(
                name=record["name"],
                white_british=[1,],
                female=False,
                female_pregnancy=[0,],
                female_menopause_binary=[0, 1,],
                female_menopause_ordinal=[0, 1, 2,],
                female_variables=[],
                female_prefixes=[],
                male=True,
                age_grade_male=[0, 1, 2,],
                male_variables=[
                    "eid", "IID",
                    "white_british",
                    "sex", "sex_text", "age", "age_grade_male",
                    "body_mass_index_log",
                    hormone,
                ],
                male_prefixes=["genotype_pc_",],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
        ))
        records.append(record)

    # Cohort: young males

    record = dict()
    record["name"] = str("table_male_young_" + hormone)
    record["cohort_model"] = "male_young"
    record["category"] = "male"
    record["phenotype"] = hormone
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            white_british=[1,],
            female=False,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[],
            female_prefixes=[],
            male=True,
            age_grade_male=[0,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex", "sex_text", "age", "age_grade_male",
                "body_mass_index_log",
                hormone,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
    ))
    records.append(record)

    # Cohort: old males

    record = dict()
    record["name"] = str("table_male_old_" + hormone)
    record["cohort_model"] = "male_old"
    record["category"] = "male"
    record["phenotype"] = hormone
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            white_british=[1,],
            female=False,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[],
            female_prefixes=[],
            male=True,
            age_grade_male=[2,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex", "sex_text", "age", "age_grade_male",
                "body_mass_index_log",
                hormone,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
    ))
    records.append(record)

    # Return information.
    return records


def stratify_cohorts_genotypes_set_sex_hormones(
    table=None,
    table_kinship_pairs=None,
    report=None,
):
    """
    Organizes tables for specific cohorts and model covariates for genetic
    analyses.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Compile information.
    records = list()
    # Select and organize variables across cohorts.
    hormones = [
        "albumin", "albumin_log", "albumin_imputation_log",
        "steroid_globulin", "steroid_globulin_log",
        "steroid_globulin_imputation_log",
        "oestradiol", "oestradiol_log", "oestradiol_imputation_log",
        "oestradiol_free", "oestradiol_free_log",
        "oestradiol_bioavailable", "oestradiol_bioavailable_log",
        "testosterone", "testosterone_log", "testosterone_imputation_log",
        "testosterone_free", "testosterone_free_log",
        "testosterone_bioavailable", "testosterone_bioavailable_log",
        "vitamin_d", "vitamin_d_log", "vitamin_d_imputation_log",
    ]
    for hormone in hormones:
        records_hormone = select_organize_cohorts_variables_by_sex_hormone(
            hormone=hormone,
            table_kinship_pairs=table_kinship_pairs,
            table=table,
        )
        records.extend(records_hormone)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "select_organize_cohorts_models_genotypes_analyses_set_sex_hormones"
        )
        print(
            "report: " + function_name
        )
        for record in records:
            utility.print_terminal_partition(level=5)
            print(record["name"])
            print(
                "Count records: " + str(record["table"].shape[0])
            )
    # Return information.
    return records


def select_organize_cohorts_variables_by_body(
    phenotype=None,
    table_kinship_pairs=None,
    table=None,
):
    """
    Organizes tables for specific cohorts and model covariates for genetic
    analyses.

    arguments:
        phenotype (str): name of column for phenotype variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: all ancestries, bipolar disorder controls
    record = dict()
    record["name"] = str(
        "table_all_bipolar_disorder_control_" + phenotype
    )
    record["cohort_model"] = "all"
    record["category"] = "general"
    record["phenotype"] = phenotype
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            white_british=[0, 1,],
            case_control="bipolar_disorder",
            case_control_values=[0,],
            variables=[
                "eid", "IID",
                #"white_british",
                "sex", "sex_text", "age",
                "body_mass_index", "body_mass_index_log",
                "bipolar_disorder",
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
    ))
    records.append(record)

    # Cohort: "White British" ancestry, bipolar disorder controls
    record = dict()
    record["name"] = str(
        "table_white_bipolar_disorder_control_" + phenotype
    )
    record["cohort_model"] = "white_british"
    record["category"] = "general"
    record["phenotype"] = phenotype
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            white_british=[1,],
            case_control="bipolar_disorder",
            case_control_values=[0,],
            variables=[
                "eid", "IID",
                "white_british",
                "sex", "sex_text", "age",
                "body_mass_index", "body_mass_index_log",
                "bipolar_disorder",
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
    ))
    records.append(record)

    # Cohort: all ancestries, bipolar disorder cases
    record = dict()
    record["name"] = str(
        "table_all_bipolar_disorder_case_" + phenotype
    )
    record["cohort_model"] = "all"
    record["category"] = "general"
    record["phenotype"] = phenotype
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            white_british=[0, 1,],
            case_control="bipolar_disorder",
            case_control_values=[1,],
            variables=[
                "eid", "IID",
                #"white_british",
                "sex", "sex_text", "age",
                "body_mass_index", "body_mass_index_log",
                "bipolar_disorder",
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
    ))
    records.append(record)

    # Cohort: "White British" ancestry, bipolar disorder cases
    record = dict()
    record["name"] = str(
        "table_white_bipolar_disorder_case_" + phenotype
    )
    record["cohort_model"] = "white_british"
    record["category"] = "general"
    record["phenotype"] = phenotype
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            white_british=[1,],
            case_control="bipolar_disorder",
            case_control_values=[1,],
            variables=[
                "eid", "IID",
                "white_british",
                "sex", "sex_text", "age",
                "body_mass_index", "body_mass_index_log",
                "bipolar_disorder",
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
    ))
    records.append(record)

    # Return information.
    return records


def stratify_cohorts_genotypes_set_bipolar_body(
    table=None,
    table_kinship_pairs=None,
    report=None,
):
    """
    Organizes tables for specific cohorts and model covariates for genetic
    analyses.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Compile information.
    records = list()
    # Select and organize variables across cohorts.
    phenotypes = [
        "body_mass_index", "body_mass_index_log",
    ]
    for phenotype in phenotypes:
        records_phenotype = select_organize_cohorts_variables_by_body(
            phenotype=phenotype,
            table_kinship_pairs=table_kinship_pairs,
            table=table,
        )
        records.extend(records_phenotype)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "select_organize_cohorts_models_genotypes_analyses_set_bipolar_body"
        )
        print(
            "report: " + function_name
        )
        for record in records:
            utility.print_terminal_partition(level=5)
            print(record["name"])
            print(
                "Count records: " + str(record["table"].shape[0])
            )
    # Return information.
    return records


##########
# Cohort, model selection: sets for descriptions of phenotypes within cohorts


def stratify_set_primary_sex_menopause_age(
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

    # Sex

    record = dict()
    record["name"] = "female"
    record["cohort_model"] = "female"
    record["category"] = "sex"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (table["sex_text"] == "female"), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male"
    record["cohort_model"] = "male"
    record["category"] = "sex"
    record["menstruation"] = False
    record["table"] = table.loc[
        (table["sex_text"] == "male"), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = "female_premenopause_ordinal"
    record["cohort_model"] = "female_premenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_ordinal"
    record["cohort_model"] = "female_perimenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_ordinal"
    record["cohort_model"] = "female_postmenopause"
    record["category"] = "menopause_ordinal"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2)
        ), :
    ]
    records.append(record)

    # Age

    record = dict()
    record["name"] = "female_younger"
    record["cohort_model"] = "female_younger"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_older"
    record["cohort_model"] = "female_older"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 2)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_younger"
    record["cohort_model"] = "male_younger"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_older"
    record["cohort_model"] = "male_older"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2)
        ), :
    ]
    records.append(record)

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
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    Persons who consumed alcohol currently.

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

    # Sex

    record = dict()
    record["name"] = "female_alcohol_current"
    record["cohort_model"] = "female_alcohol_current"
    record["category"] = "alcohol_current"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_alcohol_current"
    record["cohort_model"] = "male_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = "female_premenopause_ordinal_alcohol_current"
    record["cohort_model"] = "female_premenopause_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_perimenopause_ordinal_alcohol_current"
    record["cohort_model"] = "female_perimenopause_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = True
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_postmenopause_ordinal_alcohol_current"
    record["cohort_model"] = "female_postmenopause_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    # Age

    record = dict()
    record["name"] = "female_younger_alcohol_current"
    record["cohort_model"] = "female_younger_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 0) &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_older_alcohol_current"
    record["cohort_model"] = "female_older_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 2) &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_younger_alcohol_current"
    record["cohort_model"] = "male_younger_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 0) &
            (table["alcohol_current"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_older_alcohol_current"
    record["cohort_model"] = "male_older_alcohol_current"
    record["category"] = "alcohol_current"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2) &
            (table["alcohol_current"] == 1)
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
        stratify_set_primary_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)

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

    records_novel = (
        stratify_set_alcohol_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)

    # Return information
    return records


###############################################################################
# Drivers
# The purpose of these driver functions is to make the module's functionality
# accessible in modular parts.

# Stratification and formatting for genetic analyses


def execute_cohorts_models_genetic_analysis(
    table=None,
    set=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about cohorts and models for genetic analysis.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        set (str): name of set of cohorts and models to select and organize
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collecton of Pandas data frames of phenotype variables across
            UK Biobank cohort

    """

    # Read source information from file.
    table_kinship_pairs = read_source_table_kinship_pairs(
        path_dock=path_dock,
        report=report,
    )
    # Determine which set of cohorts and models to select and organize.
    if (set == "sex_hormones"):
        # Collect summary records and construct table.
        pail_records = (
            stratify_cohorts_genotypes_set_sex_hormones(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for collection in pail_records:
            pail[collection["name"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_phenotypes=[],
                    binary_phenotypes=[],
                    continuous_variables=[collection["phenotype"]],
                    remove_null_records=False,
                    table=collection["table"],
            ))
            pass
        pass
    elif (set == "bipolar_disorder_body"):
        pail_records = (
            stratify_cohorts_genotypes_set_bipolar_body(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for collection in pail_records:
            pail[collection["name"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_phenotypes=[],
                    binary_phenotypes=[],
                    continuous_variables=[collection["phenotype"]],
                    remove_null_records=False,
                    table=collection["table"],
            ))
            pass
        pass
    else:
        print("set of cohorts and models unrecognizable...")
        pail = dict()

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        function_name = str(
            "execute_cohorts_models_genetic_analysis()"
        )
        print("report: " + function_name)
        utility.print_terminal_partition(level=3)
        for name in pail.keys():
            print(name)
    # Return information.
    return pail





#
