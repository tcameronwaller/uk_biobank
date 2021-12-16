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
    paths["stratification"] = os.path.join(path_dock, "stratification")
    paths["reference_population"] = os.path.join(
        path_dock, "stratification", "reference_population"
    )
    paths["hormones_linear"] = os.path.join(
        path_dock, "stratification", "hormones_linear"
    )
    paths["hormones_logistic"] = os.path.join(
        path_dock, "stratification", "hormones_logistic"
    )
    paths["bipolar_body_linear"] = os.path.join(
        path_dock, "stratification", "bipolar_body_linear"
    )
    paths["bipolar_body_logistic"] = os.path.join(
        path_dock, "stratification", "bipolar_body_logistic"
    )

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["stratification"])
    # Initialize directories.
    utility.create_directories(
        path=paths["stratification"]
    )
    utility.create_directories(
        path=paths["reference_population"]
    )
    utility.create_directories(
        path=paths["hormones_linear"]
    )
    utility.create_directories(
        path=paths["hormones_logistic"]
    )
    utility.create_directories(
        path=paths["bipolar_body_linear"]
    )
    utility.create_directories(
        path=paths["bipolar_body_logistic"]
    )

    # Return information.
    return paths


##########
# Read


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


##########
# Cohort, model selection: kinship


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
        print("... Begin cohort filter for kinship ...")
        utility.print_terminal_partition(level=5)
        print("cohort: " + str(name))
        utility.print_terminal_partition(level=2)

    # Copy information.
    table = table.copy(deep=True)
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
        drop=False,
    )
    table.set_index(
        "IID",
        append=False,
        drop=True,
        inplace=True
    )
    table = table.loc[
        table.index.isin(persons_keep), :
    ]
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
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
    # Copy information.
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
        # Extract identifiers of genotypes for priority persons.
        genotypes = copy.deepcopy(
            table_preference["IID"].to_list()
        )
        genotypes = list(map(str, genotypes)) # string
        genotypes = list(set(genotypes)) # unique
        # Count genotypes.
        count_genotypes_priority = len(genotypes)
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
            pass
        pass
    pass



##########
# Format adjustment for PLINK2


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
        #value = -9
        value = "nan"
    # Return information.
    return value


def translate_binary_phenotype_plink(
    binary_value=None,
):
    """
    Translate information from simple binary representation to plink
    representation.

    Accommodate inexact float values and null values.

    https://www.cog-genomics.org/plink/2.0/input
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
        binary_value (float): binary (0, 1) representation of a phenotype

    raises:

    returns:
        (float): plink binary representation of a phenotype

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(binary_value)) and
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
        #value = -9
        value = "nan"
    # Return information.
    return value


def translate_biological_sex_plink(
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


def organize_phenotype_covariate_table_plink_format(
    boolean_variables=None,
    binary_variables=None,
    continuous_variables=None,
    sex_y=None,
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
        PAT maybe parents, parents Paternal individual ID
        MAT maybe parents, parents Maternal individual ID
        SEX sex Sex (1 = male, 2 = female, 'NA' = unknown)
        PHENO1 pheno1 All-missing phenotype column, if none loaded
        <Pheno name>, ... pheno1, phenos Phenotype value(s) (only first if just
            'pheno1')
        <Covar name>, ... (required) Covariate values

        (Note that --covar can also be used with files lacking a header row.)
    "

    arguments:
        boolean_variables (list<str>): names of columns with discrete boolean
            (False: control or negative; True: case or positive) dependent
            variables that need conversion to PLINK format (1: control; 2: case)
        binary_variables (list<str>): names of columns with discrete logical
            binary (0: control or negative; 1: case or positive) dependent
            variables that need conversion to PLINK format (control: 1, case: 2)
        continuous_variables (list<str>): names of columns that PLINK ought
            to interpret as continuous variables (ordinal discrete, interval
            continuous, or ratio continuous)
        sex_y (str): name of column for logical binary representation of
            biological sex in terms of presence of Y chromosome (0: female, XX;
            1: male, XY)
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
    for boolean_variable in boolean_variables:
        name_plink = str(boolean_variable + "_plink")
        table[name_plink] = table.apply(
            lambda row:
                translate_boolean_phenotype_plink(
                    boolean_value=row[boolean_variable],
                ),
            axis="columns", # apply across rows
        )
        pass
    # Translate binary phenotype variables.
    for binary_variable in binary_variables:
        name_plink = str(binary_variable + "_plink")
        table[name_plink] = table.apply(
            lambda row:
                translate_binary_phenotype_plink(
                    binary_value=row[binary_variable],
                ),
            axis="columns", # apply across rows
        )
        pass
    # Translate continuous variables, whether they are phenotypes or covariates.
    table = utility.convert_table_columns_variables_types_float(
        columns=continuous_variables,
        table=table,
    )
    # Sex.
    # 2 : female (XX)
    # 1 : male (XY)
    # "NA": missing
    table["SEX"] = table.apply(
        lambda row:
            translate_biological_sex_plink(
                sex_y=row[sex_y],
            ),
        axis="columns", # apply function to each row
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
        pass

    # TODO: TCW 14 December 2021
    # 1. introduce "SEX" column and place in proper order
    # 2. specify "dependent variable" column in proper order
    # 3. then come covariates

    # Introduce family identifier.
    table["FID"] = table["IID"]
    # Introduce other PLINK2 columns.
    table["SID"] = ""
    table["PAT"] = ""
    table["MAT"] = ""
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
    columns_sequence.insert(0, "MAT") # column 5
    columns_sequence.insert(0, "PAT") # column 4
    columns_sequence.insert(0, "SID") # column 3
    columns_sequence.insert(0, "IID") # column 2
    columns_sequence.insert(0, "FID") # column 1
    table_columns = table.loc[
        :, table.columns.isin(columns_sequence)
    ]
    table_sequence = table_columns[[*columns_sequence]]
    # Return information.
    return table_sequence


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
    priority_values=None,
    priority_variable=None,
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
    report=None,
):
    """
    Selects records by sex and by sex-specific criteria and variables.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        priority_values (list<int>): which values of "priority_variable" to use
            to select priority persons in kinship filter
        priority_variable (list<str>): name of column for variable to consider
            to select priority persons in kinship filter
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
        report (bool): whether to print reports

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
    table_unrelated = filter_table_phenotype_persons_by_kinship(
        name=name,
        priority_values=priority_values,
        priority_variable=priority_variable,
        threshold_kinship=0.1, # pairs with kinship >= threshold for exclusion
        table_kinship_pairs=table_kinship_pairs,
        table=table_collection,
        report=report,
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


def select_records_by_ancestry_case_control_valid_variables_values(
    name=None,
    priority_values=None,
    priority_variable=None,
    white_british=None,
    case_control=None,
    case_control_values=None,
    variables=None,
    prefixes=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Selects records for females and males without sex-specific criteria.

    arguments:
        name (str): unique name for the relevant cohort, model, and phenotype
        priority_values (list<int>): which values of "priority_variable" to use
            to select priority persons in kinship filter
        priority_variable (list<str>): name of column for variable to consider
            to select priority persons in kinship filter
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
        report (bool): whether to print reports

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
        report=report,
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
    table_unrelated = filter_table_phenotype_persons_by_kinship(
        name=name,
        priority_values=priority_values,
        priority_variable=priority_variable,
        threshold_kinship=0.1, # pairs with kinship >= threshold for exclusion
        table_kinship_pairs=table_kinship_pairs,
        table=table,
        report=report,
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


def stratify_genotype_cohorts_ordinal_hormones(
    dependence=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Organizes tables for specific cohorts, phenotypes, and model covariates
    for genetic analysis.

    These definitions are for the cohort-specific ordinal representations of
    hormones-proteins. This representation is an attempt to represent limit of
    detection missingness in the hormones-proteins by means of a cohort-specific
    median.

    arguments:
        dependence (str): name of table's column for dependent (outcome)
            variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)

    # Organize concise name for dependent variable.
    dependence_name = str(dependence + "_order")

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: non-pregnant females and males together
    record = dict()
    record["category"] = "female_male"
    record["cohort"] = "female_male"
    record["cohort_model"] = "female_male"
    dependence_specific = str(
            str(dependence) + "_" + str(record["cohort"]) + "_order"
    )
    record["dependence"] = dependence_specific
    record["name"] = str(
        record["cohort_model"] + "_" + dependence_name
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body_log",
                "pregnancy",
                dependence_specific,
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body_log",
                dependence_specific,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    record["table"][dependence_name] = record["table"][dependence_specific]
    records.append(record)

    # Cohort: all non-pregnant females together
    record = dict()
    record["category"] = "female"
    record["cohort"] = "female"
    record["cohort_model"] = "female"
    dependence_specific = str(
            str(dependence) + "_" + str(record["cohort"]) + "_order"
    )
    record["dependence"] = dependence_specific
    record["name"] = str(
        record["cohort_model"] + "_" + dependence_name
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy", "menopause_ordinal", "hormone_alteration",
                dependence_specific,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            age_grade_male=[0, 1, 2,],
            male_variables=[],
            male_prefixes=[],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    record["table"][dependence_name] = record["table"][dependence_specific]
    records.append(record)

    # Cohort: males
    record = dict()
    record["category"] = "male"
    record["cohort"] = "male"
    record["cohort_model"] = "male"
    dependence_specific = str(
            str(dependence) + "_" + str(record["cohort"]) + "_order"
    )
    record["dependence"] = dependence_specific
    record["name"] = str(
        record["cohort_model"] + "_" + dependence_name
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
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
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "age_grade_male",
                "body", "body_log",
                dependence_specific,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    record["table"][dependence_name] = record["table"][dependence_specific]
    records.append(record)

    # Return information.
    return records


def stratify_genotype_cohorts_hormones_by_sex_together(
    dependence=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
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
        dependence (str): name of table's column for dependent (outcome)
            variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): collection of information about phenotype variables in
            separate cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: non-pregnant females and males together
    # Kinship priority: none
    record = dict()
    record["category"] = "female_male"
    record["cohort"] = "female_male"
    record["cohort_model"] = "female_male"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_none = record["table"].copy(deep=True)
    records.append(record)

    # Cohort: non-pregnant females and males together
    # Kinship priority: female
    record = dict()
    record["category"] = "female_male"
    record["cohort"] = "female_male"
    record["cohort_model"] = "female_male_priority_female"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=["female",],
            priority_variable="sex_text",
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_female = record["table"].copy(deep=True)
    records.append(record)

    # Cohort: non-pregnant females and males together
    # Kinship priority: male
    record = dict()
    record["category"] = "female_male"
    record["cohort"] = "female_male"
    record["cohort_model"] = "female_male_priority_male"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=["male",],
            priority_variable="sex_text",
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_male = record["table"].copy(deep=True)
    records.append(record)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "stratify_genotype_cohorts_hormones_by_sex_together"
        )
        print(
            "report: " + function_name
        )
        utility.print_terminal_partition(level=3)
        print("dependence: " + str(dependence))
        report_kinship_filter_priority_selection(
            name="... Comparison of priority to female persons in population ...",
            priority_values=["female",],
            priority_variable="sex_text",
            table_full=table,
            table_simple=table_priority_none,
            table_priority=table_priority_female,
            report=report,
        )
        report_kinship_filter_priority_selection(
            name="... Comparison of priority to male persons in population ...",
            priority_values=["male",],
            priority_variable="sex_text",
            table_full=table,
            table_simple=table_priority_none,
            table_priority=table_priority_male,
            report=report,
        )
    # Return information.
    return records


def stratify_genotype_cohorts_hormones_by_female_menopause(
    dependence=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
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
        dependence (str): name of table's column for dependent (outcome)
            variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): collection of information about phenotype variables in
            separate cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: all non-pregnant females together
    record = dict()
    record["category"] = "female"
    record["cohort"] = "female"
    record["cohort_model"] = "female"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy", "menopause_ordinal", "hormone_alteration",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            age_grade_male=[0, 1, 2,],
            male_variables=[],
            male_prefixes=[],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: premenopausal females by ordinal menopause definition
    record = dict()
    record["category"] = "female_menopause_ordinal"
    record["cohort"] = "female_premenopause"
    record["cohort_model"] = "female_premenopause"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "menstruation_days",
                "menstruation_phase", "menstruation_phase_cycle",
                "hormone_alteration",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            age_grade_male=[0, 1, 2,],
            male_variables=[],
            male_prefixes=[],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: perimenopausal females by ordinal menopause definition
    record = dict()
    record["category"] = "female_menopause_ordinal"
    record["cohort"] = "female_perimenopause"
    record["cohort_model"] = "female_perimenopause"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[1,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "menstruation_days",
                "menstruation_phase", "menstruation_phase_cycle",
                "hormone_alteration",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            age_grade_male=[0, 1, 2,],
            male_variables=[],
            male_prefixes=[],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: postmenopausal females by ordinal menopause definition
    record = dict()
    record["category"] = "female_menopause_ordinal"
    record["cohort"] = "female_postmenopause"
    record["cohort_model"] = "female_postmenopause"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[2,],
            female_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "body", "body_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "hormone_alteration",
                dependence,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            age_grade_male=[0, 1, 2,],
            male_variables=[],
            male_prefixes=[],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Return information.
    return records


def stratify_genotype_cohorts_hormones_by_male_age(
    dependence=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
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
        dependence (str): name of table's column for dependent (outcome)
            variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): collection of information about phenotype variables in
            separate cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: males
    record = dict()
    record["category"] = "male"
    record["cohort"] = "male"
    record["cohort_model"] = "male"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
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
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "age_grade_male",
                "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: young males
    record = dict()
    record["category"] = "male_age"
    record["cohort"] = "male_age_low"
    record["cohort_model"] = "male_age_low"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
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
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "age_grade_male",
                "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    record = dict()
    record["category"] = "male_age"
    record["cohort"] = "male_age_middle"
    record["cohort_model"] = "male_age_middle"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=False,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[],
            female_prefixes=[],
            male=True,
            age_grade_male=[1,],
            male_variables=[
                "eid", "IID",
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "age_grade_male",
                "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: old males

    record = dict()
    record["category"] = "male_age"
    record["cohort"] = "male_age_high"
    record["cohort_model"] = "male_age_high"
    record["dependence"] = dependence
    record["name"] = str(record["cohort_model"] + "_" + dependence)
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
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
                "assessment_region",
                "assessment_season",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age", "age_grade_male",
                "body", "body_log",
                dependence,
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Return information.
    return records


def stratify_genotype_cohorts_hormones_by_sex_menopause_age(
    dependence=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
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
        dependence (str): name of table's column for dependent (outcome)
            variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): collection of information about phenotype variables in
            separate cohorts

    """

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohorts with females and males together and Kinship selection priority for
    # either females or males.
    records_novel = (
        stratify_genotype_cohorts_hormones_by_sex_together(
            dependence=dependence,
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.extend(records_novel)

    if False:
        # Cohorts with females in different stages before or after menopause.
        records_novel = (
            stratify_genotype_cohorts_hormones_by_female_menopause(
                dependence=dependence,
                table_kinship_pairs=table_kinship_pairs,
                table=table,
                report=report,
        ))
        records.extend(records_novel)

        # Cohorts with males in different tertiles by age.
        records_novel = (
            stratify_genotype_cohorts_hormones_by_male_age(
                dependence=dependence,
                table_kinship_pairs=table_kinship_pairs,
                table=table,
                report=report,
        ))
        records.extend(records_novel)

    # Return information.
    return records


def stratify_genotype_cohorts_set_reference_population(
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

    # Copy information in table.
    table = table.copy(deep=True)
    table_kinship_pairs = table_kinship_pairs.copy(deep=True)

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    # Cohort: "White British" ancestry, unrelated, females and males together
    record = dict()
    record["category"] = "white_unrelated_female_male"
    record["cohort"] = "white_unrelated_female_male"
    record["cohort_model"] = "white_unrelated_female_male"
    record["dependence"] = "age"
    record["name"] = str(record["cohort_model"])
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,1,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_none = record["table"].copy(deep=True)
    records.append(record)

    # Cohort: "White British" ancestry, unrelated, females and males together
    record = dict()
    record["category"] = "white_unrelated_female_male"
    record["cohort"] = "white_unrelated_female_male"
    record["cohort_model"] = "white_unrelated_female_male_priority_female"
    record["dependence"] = "age"
    record["name"] = str(record["cohort_model"])
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=["female",],
            priority_variable="sex_text",
            white_british=[1,],
            female=True,
            female_pregnancy=[0,1,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_female = record["table"].copy(deep=True)
    records.append(record)

    # Cohort: "White British" ancestry, unrelated, females and males together
    record = dict()
    record["category"] = "white_unrelated_female_male"
    record["cohort"] = "white_unrelated_female_male"
    record["cohort_model"] = "white_unrelated_female_male_priority_male"
    record["dependence"] = "age"
    record["name"] = str(record["cohort_model"])
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=["male",],
            priority_variable="sex_text",
            white_british=[1,],
            female=True,
            female_pregnancy=[0,1,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_male = record["table"].copy(deep=True)
    records.append(record)

    # Cohort: "White British" ancestry, unrelated, females
    record = dict()
    record["category"] = "white_unrelated_female"
    record["cohort"] = "white_unrelated_female"
    record["cohort_model"] = "white_unrelated_female"
    record["dependence"] = "age"
    record["name"] = str(record["cohort_model"])
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=True,
            female_pregnancy=[0,1,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: "White British" ancestry, unrelated, males
    record = dict()
    record["category"] = "white_unrelated_male"
    record["cohort"] = "white_unrelated_male"
    record["cohort_model"] = "white_unrelated_male"
    record["dependence"] = "age"
    record["name"] = str(record["cohort_model"])
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_sex_specific_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            female=False,
            female_pregnancy=[0,1,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            age_grade_male=[0, 1, 2,],
            male_variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
            ],
            male_prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "stratify_genotype_cohorts_set_reference_population"
        )
        print(
            "report: " + function_name
        )
        utility.print_terminal_partition(level=3)
        # Describe effect of preferential selection of females and males in the
        # Kinship Filter for genetic analyses.
        report_kinship_filter_priority_selection(
            name="... Comparison of priority to female persons in population ...",
            priority_values=["female",],
            priority_variable="sex_text",
            table_full=table,
            table_simple=table_priority_none,
            table_priority=table_priority_female,
            report=report,
        )
        report_kinship_filter_priority_selection(
            name="... Comparison of priority to male persons in population ...",
            priority_values=["male",],
            priority_variable="sex_text",
            table_full=table,
            table_simple=table_priority_none,
            table_priority=table_priority_male,
            report=report,
        )
        for record in records:
            utility.print_terminal_partition(level=5)
            print(record["name"])
            print(
                "Count records: " + str(record["table"].shape[0])
            )
    # Return information.
    return records


def stratify_genotype_cohorts_linear_set_hormones_proteins(
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
    dependences = [
        "vitamin_d", #"vitamin_d_log",
        "vitamin_d_imputation", #"vitamin_d_imputation_log",
        "albumin", #"albumin_log",
        "albumin_imputation", #"albumin_imputation_log",
        "steroid_globulin", #"steroid_globulin_log",
        "steroid_globulin_imputation", #"steroid_globulin_imputation_log",
        "oestradiol", #"oestradiol_log",
        "oestradiol_imputation", #"oestradiol_imputation_log",
        "oestradiol_bioavailable", #"oestradiol_bioavailable_log",
        "oestradiol_bioavailable_imputation",
        #"oestradiol_bioavailable_imputation_log",
        "oestradiol_free", #"oestradiol_free_log",
        "oestradiol_free_imputation", #"oestradiol_free_imputation_log",
        "testosterone", #"testosterone_log",
        "testosterone_imputation", #"testosterone_imputation_log",
        "testosterone_bioavailable", #"testosterone_bioavailable_log",
        "testosterone_bioavailable_imputation",
        #"testosterone_bioavailable_imputation_log",
        "testosterone_free", #"testosterone_free_log",
        "testosterone_free_imputation", #"testosterone_free_imputation_log",
    ]
    for dependence in dependences:
        records_hormone = (
            stratify_genotype_cohorts_hormones_by_sex_menopause_age(
                dependence=dependence,
                table_kinship_pairs=table_kinship_pairs,
                table=table,
                report=report,
        ))
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


def stratify_genotype_cohorts_logistic_set_hormones_proteins(
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
        "vitamin_d_detection",
        "albumin_detection",
        "steroid_globulin_detection",
        "oestradiol_detection",
        "testosterone_detection",
    ]
    for phenotype in phenotypes:
        records_phenotype = (
            stratify_genotype_cohorts_hormones_by_sex_menopause_age(
                dependence=phenotype,
                table_kinship_pairs=table_kinship_pairs,
                table=table,
                report=report,
        ))
        records.extend(records_phenotype)
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


def stratify_linear_set_bipolar_body_by_bipolar(
    case_control=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Organizes tables for specific cohorts and model covariates for genetic
    analyses.

    arguments:
        case_control (str): name of table's column for case control definition
            for cohort stratification or dependent variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Organize concise name for phenotype response variable.
    if "loose" in case_control:
        case_control_abbreviation = "bipolar_loose"
    elif "strict" in case_control:
        case_control_abbreviation = "bipolar_strict"
    else:
        case_control_abbreviation = "unknown"

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    ##########
    # Bipolar Disorder definition: loose

    if False:
        # Cohort: all ancestries, loose bipolar disorder controls
        record = dict()
        record["category"] = "general"
        record["cohort"] = "all_bipolar_control"
        record["cohort_model"] = "all_bipolar_control"
        record["dependence"] = "body"
        record["name"] = str(
            "all_" + case_control_abbreviation + "_control"
        )
        record["name_table"] = str("table_" + record["name"])
        record["table"] = (
            select_records_by_ancestry_case_control_valid_variables_values(
                name=record["name"],
                priority_values=[],
                priority_variable=None,
                white_british=[0, 1,],
                case_control=case_control,
                case_control_values=[0,],
                variables=[
                    "eid", "IID",
                    #"white_british",
                    "sex_y", "sex_x", "sex_text", "age",
                    "body", "body_log",
                    case_control,
                ],
                prefixes=["genotype_pc_",],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
                report=report,
        ))
        records.append(record)
        pass

    # Cohort: "White British" ancestry, loose Bipolar Disorder controls
    record = dict()
    record["category"] = "general"
    record["cohort"] = "white_bipolar_control"
    record["cohort_model"] = "white_bipolar_control"
    record["dependence"] = "body"
    record["name"] = str(
        "white_" + case_control_abbreviation + "_control_body"
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            case_control=case_control,
            case_control_values=[0,],
            variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
                "body", "body_log",
                case_control,
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Cohort: "White British" ancestry, Bipolar Disorder cases
    record = dict()
    record["category"] = "general"
    record["cohort"] = "white_bipolar_case"
    record["cohort_model"] = "white_bipolar_case"
    record["dependence"] = "body"
    record["name"] = str(
        "white_" + case_control_abbreviation + "_case_body"
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            case_control=case_control,
            case_control_values=[1,],
            variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
                "body", "body_log",
                case_control,
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    records.append(record)

    # Return information.
    return records


def stratify_genotype_cohorts_linear_set_bipolar_body(
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
        "bipolar_control_case_loose",
        "bipolar_control_case_strict",
    ]
    for phenotype in phenotypes:
        records_phenotype = stratify_linear_set_bipolar_body_by_bipolar(
            case_control=phenotype,
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
        )
        records.extend(records_phenotype)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "stratify_genotype_cohorts_linear_set_bipolar_body()"
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


def stratify_logistic_set_bipolar_body_by_bipolar(
    case_control=None,
    table_kinship_pairs=None,
    table=None,
    report=None,
):
    """
    Organizes tables for specific cohorts and model covariates for genetic
    analyses.

    arguments:
        case_control (str): name of table's column for case control definition
            for cohort stratification or dependent variable
        table_kinship_pairs (object): Pandas data frame of kinship coefficients
            across pairs of persons in UK Biobank cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Collect records of information about each cohort, model, and phenotype.
    # Select and organize variables across cohorts.
    records = list()

    if False:
        # Cohort: all ancestries
        record = dict()
        record["category"] = "general"
        record["cohort"] = "all"
        record["cohort_model"] = "all"
        record["dependence"] = case_control
        record["name"] = str(
            "all_" + case_control
        )
        record["name_table"] = str("table_" + record["name"])
        record["table"] = (
            select_records_by_ancestry_case_control_valid_variables_values(
                name=record["name"],
                priority_values=[],
                priority_variable=None,
                white_british=[0, 1,],
                case_control=case_control,
                case_control_values=[0, 1,],
                variables=[
                    "eid", "IID",
                    #"white_british",
                    "sex_y", "sex_x", "sex_text", "age",
                    "body", "body_log",
                    case_control,
                ],
                prefixes=["genotype_pc_",],
                table_kinship_pairs=table_kinship_pairs,
                table=table,
                report=report,
        ))
        records.append(record)
        pass

    # Cohort: "White British" ancestry
    # Note: Without priority for Bipolar Disorder cases in Kinship Filter
    record = dict()
    record["category"] = "general"
    record["cohort"] = "white"
    record["cohort_model"] = "white"
    record["dependence"] = case_control
    record["name"] = str(
        "white_" + case_control
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            priority_values=[],
            priority_variable=None,
            white_british=[1,],
            case_control=case_control,
            case_control_values=[0, 1,],
            variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
                "body", "body_log",
                case_control,
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_none = record["table"].copy(deep=True)
    records.append(record)

    # Cohort: "White British" ancestry
    # Note: With priority for Bipolar Disorder cases in Kinship Filter
    record = dict()
    record["category"] = "general"
    record["cohort"] = "white"
    record["cohort_model"] = "white"
    record["dependence"] = case_control
    record["name"] = str(
        "white_" + case_control + "_priority_case"
    )
    record["name_table"] = str("table_" + record["name"])
    record["table"] = (
        select_records_by_ancestry_case_control_valid_variables_values(
            name=record["name"],
            priority_values=[1,],
            priority_variable=case_control,
            white_british=[1,],
            case_control=case_control,
            case_control_values=[0, 1,],
            variables=[
                "eid", "IID",
                "white_british",
                "sex_y", "sex_x", "sex_text", "age",
                "body", "body_log",
                case_control,
            ],
            prefixes=["genotype_pc_",],
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
    ))
    table_priority_case = record["table"].copy(deep=True)
    records.append(record)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "stratify_logistic_set_bipolar_body_by_bipolar"
        )
        print(
            "report: " + function_name
        )
        utility.print_terminal_partition(level=3)
        # Describe effect of preferential selection of females and males in the
        # Kinship Filter for genetic analyses.
        report_kinship_filter_priority_selection(
            name="... Comparison of case priority for Bipolar Disorder ...",
            priority_values=[1,],
            priority_variable=case_control,
            table_full=table,
            table_simple=table_priority_none,
            table_priority=table_priority_case,
            report=report,
        )
        for record in records:
            utility.print_terminal_partition(level=5)
            print(record["name"])
            print(
                "Count records: " + str(record["table"].shape[0])
            )
    # Return information.
    return records


def stratify_genotype_cohorts_logistic_set_bipolar_body(
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
        "bipolar_control_case_loose",
        "bipolar_control_case_strict",
    ]
    for phenotype in phenotypes:
        records_phenotype = stratify_logistic_set_bipolar_body_by_bipolar(
            case_control=phenotype,
            table_kinship_pairs=table_kinship_pairs,
            table=table,
            report=report,
        )
        records.extend(records_phenotype)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        function_name = str(
            "stratify_genotype_cohorts_logistic_set_bipolar_body()"
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
# Drivers for genotype cohorts


def execute_stratify_genotype_cohorts_plink_format_set(
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
    # Collect stratification information and tables.
    if (set == "reference_population"):
        records = (
            stratify_genotype_cohorts_set_reference_population(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for record in records:
            pail[record["name_table"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_variables=[],
                    binary_variables=[],
                    continuous_variables=[],
                    sex_y="sex_y",
                    remove_null_records=False,
                    table=record["table"],
            ))
            pass
    elif (set == "hormones_linear"):
        records = (
            stratify_genotype_cohorts_linear_set_hormones_proteins(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for record in records:
            pail[record["name_table"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_variables=[],
                    binary_variables=[],
                    continuous_variables=[], # record["dependence"]
                    sex_y="sex_y",
                    remove_null_records=False,
                    table=record["table"],
            ))
            pass
    elif (set == "hormones_logistic"):
        records = (
            stratify_genotype_cohorts_logistic_set_hormones_proteins(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for record in records:
            pail[record["name_table"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_variables=[],
                    binary_variables=[record["dependence"]],
                    continuous_variables=[],
                    sex_y="sex_y",
                    remove_null_records=False,
                    table=record["table"],
            ))
            pass
    elif (set == "bipolar_body_linear"):
        records = (
            stratify_genotype_cohorts_linear_set_bipolar_body(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for record in records:
            pail[record["name_table"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_variables=[],
                    binary_variables=[],
                    continuous_variables=[],
                    sex_y="sex_y",
                    remove_null_records=False,
                    table=record["table"],
            ))
            pass
    elif (set == "bipolar_body_logistic"):
        records = (
            stratify_genotype_cohorts_logistic_set_bipolar_body(
                table=table,
                table_kinship_pairs=table_kinship_pairs,
                report=report,
        ))
        # Translate variable encodings and table format for analysis in PLINK2.
        pail = dict()
        for record in records:
            pail[record["name_table"]] = (
                organize_phenotype_covariate_table_plink_format(
                    boolean_variables=[],
                    binary_variables=[record["dependence"]],
                    continuous_variables=[],
                    sex_y="sex_y",
                    remove_null_records=False,
                    table=record["table"],
            ))
            pass
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        function_name = str(
            "execute_stratify_genotype_cohorts_plink_format_set()"
        )
        print("report: " + function_name)
        utility.print_terminal_partition(level=3)
        for name in pail.keys():
            print(name)
    # Return information.
    return pail



##########
# Phenotype cohorts
# Cohort, model selection: sets for descriptions of phenotypes within cohorts


def stratify_set_primary_sex_age_body_menopause(
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
    record["cohort"] = "female"
    record["cohort_model"] = "female"
    record["category"] = "sex"
    record["phenotype"] = "null"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male"
    record["cohort"] = "male"
    record["cohort_model"] = "male"
    record["category"] = "sex"
    record["menstruation"] = False
    record["table"] = table.loc[
        (table["sex_text"] == "male"), :
    ]
    records.append(record)

    # Menopause

    record = dict()
    record["name"] = "female_premenopause"
    record["cohort"] = "female_premenopause"
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
    record["name"] = "female_perimenopause"
    record["cohort"] = "female_perimenopause"
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
    record["name"] = "female_postmenopause"
    record["cohort"] = "female_postmenopause"
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
    record["name"] = "female_age_low"
    record["cohort"] = "female_age_low"
    record["cohort_model"] = "female_age_low"
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
    record["name"] = "female_age_middle"
    record["cohort"] = "female_age_middle"
    record["cohort_model"] = "female_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["age_grade_female"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_age_high"
    record["cohort"] = "female_age_high"
    record["cohort_model"] = "female_age_high"
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
    record["name"] = "male_age_low"
    record["cohort"] = "male_age_low"
    record["cohort_model"] = "male_age_low"
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
    record["name"] = "male_age_middle"
    record["cohort"] = "male_age_middle"
    record["cohort_model"] = "male_age_middle"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_age_high"
    record["cohort"] = "male_age_high"
    record["cohort_model"] = "male_age_high"
    record["category"] = "age"
    record["menstruation"] = False
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["age_grade_male"] == 2)
        ), :
    ]
    records.append(record)

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



##########
# Write


def write_genotype_product_cohort_model_table(
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
        path_parent, str(name + ".pickle")
    )
    path_table_text = os.path.join(
        path_parent, str(name + ".tsv")
    )
    # Write information to file.
    information.to_pickle(
        path_table
    )
    information.to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=False,
        na_rep="NA",
    )
    pass


def write_genotype_product_cohorts_models(
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
        write_genotype_product_cohort_model_table(
            name=name,
            information=information[name],
            path_parent=path_parent,
        )
    pass


def write_genotype_product(
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

    # Cohort tables in PLINK format.
    write_genotype_product_cohorts_models(
        information=information["reference_population"],
        path_parent=paths["reference_population"],
    )
    write_genotype_product_cohorts_models(
        information=information["hormones_linear"],
        path_parent=paths["hormones_linear"],
    )
    write_genotype_product_cohorts_models(
        information=information["hormones_logistic"],
        path_parent=paths["hormones_logistic"],
    )
    write_genotype_product_cohorts_models(
        information=information["bipolar_body_linear"],
        path_parent=paths["bipolar_body_linear"],
    )
    write_genotype_product_cohorts_models(
        information=information["bipolar_body_logistic"],
        path_parent=paths["bipolar_body_logistic"],
    )
    pass




###############################################################################
