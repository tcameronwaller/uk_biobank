"""
Prepare plots from data of the U.K. Biobank.

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
import promiscuity.plot as plot
import promiscuity.scale as scale
import promiscuity.regression as pro_reg
import uk_biobank.stratification as ukb_strat
import uk_biobank.organization as ukb_organization

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
    paths["plot"] = os.path.join(path_dock, "plot")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["plot"])
    # Initialize directories.
    utility.create_directories(path=paths["plot"])
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
                "cohort": "string",
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
        name_cohort = str(file_name).replace(".pickle", "").replace("table_", "")
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
        record["name"] = name_cohort
        record["cohort"] = name_cohort
        record["table"] = table_cohort
        records.append(record)
        pass

    # Return information.
    return records



##########
# Controls for plotting functions
##########


# Histogram


def create_plot_histogram(
    label_title=None,
    array=None,
    bins=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        label_title (str): label title for plot
        array (object): NumPy array of values to bin and plot in histogram
        bins (int): count of bins for histogram

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()
    # Determine bin method.
    if bins is None:
        bin_method = "auto"
    else:
        bin_method = "count"
    # Create figure.
    figure = plot.plot_distribution_histogram(
        array=array,
        title="",
        bin_method=bin_method, # "auto" or "count"
        bin_count=bins,
        bar_width=0.5, # 0.35, 0.50
        label_bins="values",
        label_counts="counts per bin",
        fonts=fonts,
        colors=colors,
        line=True,
        line_position=numpy.nanmean(array),
        label_title=label_title, # ""
        label_report=False, # report count, median, and mean on plot
    )
    # Return.
    return figure


def define_phenotypes_parameters_for_plots_histogram():
    """
    Organizes parameters for histogram plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information about each cohort and model.
    records = list()

    #"age",
    #"body",
    #"month",
    #"alcohol_drinks_monthly_combination",
    #"alcohol_frequency",
    #"alcohol_auditc",

    #"albumin",
    #"albumin_imputation",
    #"steroid_globulin",
    #"steroid_globulin_imputation", "steroid_globulin_imputation_log",
    #"oestradiol",
    #"oestradiol_imputation", "oestradiol_imputation_log",
    #"oestradiol_bioavailable",
    #"oestradiol_bioavailable_imputation",
    #"oestradiol_bioavailable_imputation_log",
    #"oestradiol_free",
    #"oestradiol_free_imputation", "oestradiol_free_imputation_log",
    #"testosterone", "testosterone_imputation",
    #"testosterone_bioavailable", "testosterone_bioavailable_imputation",
    #"testosterone_free", "testosterone_free_imputation",
    #"vitamin_d", "vitamin_d_imputation",
    #"alcohol_frequency",
    #"alcohol_drinks_weekly",
    #"alcohol_drinks_monthly",
    #"alcohol_drinks_monthly_combination",

    if False:

        # Sex-steroid hormones.

        record = dict()
        record["name"] = "oestradiol_imputation"
        record["variable"] = "oestradiol_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_imputation_log"
        record["variable"] = "oestradiol_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_imputation_z"
        record["variable"] = "oestradiol_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_imputation_rank"
        record["variable"] = "oestradiol_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_bioavailable_imputation"
        record["variable"] = "oestradiol_bioavailable_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_bioavailable_imputation_log"
        record["variable"] = "oestradiol_bioavailable_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_bioavailable_imputation_z"
        record["variable"] = "oestradiol_bioavailable_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_bioavailable_imputation_rank"
        record["variable"] = "oestradiol_bioavailable_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_free_imputation"
        record["variable"] = "oestradiol_free_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_free_imputation_log"
        record["variable"] = "oestradiol_free_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_free_imputation_z"
        record["variable"] = "oestradiol_free_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_free_imputation_rank"
        record["variable"] = "oestradiol_free_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_imputation"
        record["variable"] = "testosterone_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_imputation_log"
        record["variable"] = "testosterone_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_imputation_z"
        record["variable"] = "testosterone_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_imputation_rank"
        record["variable"] = "testosterone_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_bioavailable_imputation"
        record["variable"] = "testosterone_bioavailable_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_bioavailable_imputation_log"
        record["variable"] = "testosterone_bioavailable_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_bioavailable_imputation_z"
        record["variable"] = "testosterone_bioavailable_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_bioavailable_imputation_rank"
        record["variable"] = "testosterone_bioavailable_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_free_imputation"
        record["variable"] = "testosterone_free_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_free_imputation_log"
        record["variable"] = "testosterone_free_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_free_imputation_z"
        record["variable"] = "testosterone_free_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_free_imputation_rank"
        record["variable"] = "testosterone_free_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        pass

    if True:

        # Binding Proteins for Sex-steroid hormones.

        record = dict()
        record["name"] = "albumin_imputation"
        record["variable"] = "albumin_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "albumin_imputation_log"
        record["variable"] = "albumin_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "albumin_imputation_z"
        record["variable"] = "albumin_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "albumin_imputation_rank"
        record["variable"] = "albumin_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "steroid_globulin_imputation"
        record["variable"] = "steroid_globulin_imputation"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "steroid_globulin_imputation_log"
        record["variable"] = "steroid_globulin_imputation_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "steroid_globulin_imputation_z"
        record["variable"] = "steroid_globulin_imputation_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "steroid_globulin_imputation_rank"
        record["variable"] = "steroid_globulin_imputation_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        pass


    if False:

        # Measures of alcohol consumption.

        record = dict()
        record["name"] = "alcohol_drinks_monthly_combination"
        record["variable"] = "alcohol_drinks_monthly_combination"
        record["threshold"] = 400 # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_monthly_combination_log"
        record["variable"] = "alcohol_drinks_monthly_combination_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_monthly_combination_z"
        record["variable"] = "alcohol_drinks_monthly_combination_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_monthly_combination_rank"
        record["variable"] = "alcohol_drinks_monthly_combination_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_frequency"
        record["variable"] = "alcohol_frequency"
        record["threshold"] = 15 # None or upper threshold
        record["bins"] = 10 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_frequency_log"
        record["variable"] = "alcohol_frequency_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 10 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_frequency_z"
        record["variable"] = "alcohol_frequency_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 10 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_frequency_rank"
        record["variable"] = "alcohol_frequency_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 10 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_auditc"
        record["variable"] = "alcohol_auditc"
        record["threshold"] = 15 # None or upper threshold
        record["bins"] = 15 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_auditc_log"
        record["variable"] = "alcohol_auditc_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 15 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_auditc_z"
        record["variable"] = "alcohol_auditc_z"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 15 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_auditc_rank"
        record["variable"] = "alcohol_auditc_rank"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 15 # None or count of bins
        records.append(record)

        pass

    if False:

        record = dict()
        record["name"] = "age"
        record["variable"] = "age"
        record["threshold"] = None # None or upper threshold
        record["bins"] = None # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "body"
        record["variable"] = "body"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_weekly"
        record["variable"] = "alcohol_drinks_weekly"
        record["threshold"] = 100 # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_monthly"
        record["variable"] = "alcohol_drinks_monthly"
        record["threshold"] = 40 # None or upper threshold
        record["bins"] = 40 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_weekly_log"
        record["variable"] = "alcohol_drinks_weekly_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_monthly_log"
        record["variable"] = "alcohol_drinks_monthly_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 40 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "alcohol_drinks_monthly_combination_log"
        record["variable"] = "alcohol_drinks_monthly_combination_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol"
        record["variable"] = "oestradiol"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "oestradiol_log"
        record["variable"] = "oestradiol_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone"
        record["variable"] = "testosterone"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        record = dict()
        record["name"] = "testosterone_log"
        record["variable"] = "testosterone_log"
        record["threshold"] = None # None or upper threshold
        record["bins"] = 70 # None or count of bins
        records.append(record)

        pass

    # Return information
    return records


def drive_create_plots_histogram(
    records_cohorts=None,
    report=None,
):
    """
    Organizes description histogram plots.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Organize phenotypes and parameters for plots.
    records_phenotypes = define_phenotypes_parameters_for_plots_histogram()

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record_cohort in records_cohorts:
        #if (record_cohort["menstruation"]):
        #    #phenotypes.append("menstruation_days_threshold")
        #    pass
        # Iterate on phenotypes.
        for record_phenotype in records_phenotypes:
            # Organize information for plot.
            name = record_phenotype["name"]
            threshold = record_phenotype["threshold"]
            variable = record_phenotype["variable"]
            bins = record_phenotype["bins"]
            # Copy information.
            table_cohort = record_cohort["table"].copy(deep=True)
            # Apply threshold.
            if (threshold is not None):
                table_cohort = table_cohort.loc[
                    (table_cohort[variable] < threshold), :
                ]
            # Histogram plot.
            name_label = str(str(record_cohort["name"]) + "_" + str(name))
            name_plot = str(name_label + "_histogram")
            pail[name_plot] = create_plot_histogram(
                label_title="", # name_label,
                array=table_cohort[variable].dropna().to_numpy(),
                bins=bins, # None or count of bins
            )
            pass
        pass
    # Return information.
    return pail


# Box Plots


def define_parameters_box_plots_female_male():
    """
    Organizes parameters for box plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information for each plot.
    records = list()

    # Plots of Female and Male cohorts together.

    record = dict()
    record["name_figure"] = "estradiol_total_female_male"
    record["variable"] = "oestradiol_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Total (pmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "estradiol_bioavailable_female_male"
    record["variable"] = "oestradiol_bioavailable_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Bioavailable (pmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "estradiol_free_female_male"
    record["variable"] = "oestradiol_free_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Free (pmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_total_female_male"
    record["variable"] = "testosterone_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Total (pmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_bioavailable_female_male"
    record["variable"] = "testosterone_bioavailable_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Bioavailable (pmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_free_female_male"
    record["variable"] = "testosterone_free_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Free (pmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "steroid_globulin_female_male"
    record["variable"] = "steroid_globulin_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "SHBG (nmol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "albumin_female_male"
    record["variable"] = "albumin_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Albumin (umol/L)"
    record["aspect"] = "portrait"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    # Return information
    return records


def define_parameters_box_plots_female():
    """
    Organizes parameters for box plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information for each plot.
    records = list()

    # Plots of Female and Male cohorts together.

    record = dict()
    record["name_figure"] = "estradiol_total_female"
    record["variable"] = "oestradiol_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Total (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
        #"ancestry_white_female_menstruation_regular",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
        #"Female-Menstruation",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "estradiol_bioavailable_female"
    record["variable"] = "oestradiol_bioavailable_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Bioavailable (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "estradiol_free_female"
    record["variable"] = "oestradiol_free_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Free (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_total_female"
    record["variable"] = "testosterone_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Total (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_bioavailable_female"
    record["variable"] = "testosterone_bioavailable_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Bioavailable (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_free_female"
    record["variable"] = "testosterone_free_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Free (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "steroid_globulin_female"
    record["variable"] = "steroid_globulin_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "SHBG (nmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "albumin_female"
    record["variable"] = "albumin_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Albumin (umol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_female",
        "ancestry_white_female_premenopause",
        "ancestry_white_female_perimenopause",
        "ancestry_white_female_postmenopause",
    ]
    record["titles_abscissa_groups"] = [
        "Female-All",
        "Female-Pre...",
        "Female-Peri...",
        "Female-Post...",
    ]
    record["colors_groups"] = [
        (0.039, 0.196, 0.588, 1.0), # "blue_navy"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
        (0.196, 0.588, 1.0, 1.0), # "blue_sky"
    ]
    records.append(record)

    # Return information
    return records


def define_parameters_box_plots_male():
    """
    Organizes parameters for box plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information for each plot.
    records = list()

    # Plots of Female and Male cohorts together.

    record = dict()
    record["name_figure"] = "estradiol_total_male"
    record["variable"] = "oestradiol_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Total (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "estradiol_bioavailable_male"
    record["variable"] = "oestradiol_bioavailable_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Bioavailable (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "estradiol_free_male"
    record["variable"] = "oestradiol_free_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Estradiol, Free (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_total_male"
    record["variable"] = "testosterone_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Total (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_bioavailable_male"
    record["variable"] = "testosterone_bioavailable_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Bioavailable (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "testosterone_free_male"
    record["variable"] = "testosterone_free_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Testosterone, Free (pmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "steroid_globulin_male"
    record["variable"] = "steroid_globulin_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "SHBG (nmol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    record = dict()
    record["name_figure"] = "albumin_male"
    record["variable"] = "albumin_imputation"
    record["title_ordinate"] = ""
    record["title_abscissa"] = "" # "Albumin (umol/L)"
    record["aspect"] = "portrait_half_width"
    record["orientation_box"] = "vertical"
    record["cohorts_groups"] = [
        "ancestry_white_male",
        "ancestry_white_male_age_low",
        "ancestry_white_male_age_middle",
        "ancestry_white_male_age_high",
    ]
    record["titles_abscissa_groups"] = [
        "Male-All",
        "Male-Age-Low",
        "Male-Age-Middle",
        "Male-Age-High",
    ]
    record["colors_groups"] = [
        (0.510, 0.039, 0.510, 1.0), # "purple"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
        (0.784, 0.275, 0.784, 1.0), # "magenta"
    ]
    records.append(record)

    # Return information
    return records


def define_parameters_for_plots_groups_box():
    """
    Organizes parameters for box plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information for each plot.
    records = list()

    # Plots of Female and Male cohorts together.
    records_both = define_parameters_box_plots_female_male()
    records.extend(records_both)

    # Plots of Female cohorts separate.
    records_female = define_parameters_box_plots_female()
    records.extend(records_female)

    # Plots of Male cohorts separate.
    records_male = define_parameters_box_plots_male()
    records.extend(records_male)

    # Return information
    return records


def create_plot_box(
    name_figure=None,
    variable=None,
    title_ordinate=None,
    title_abscissa=None,
    cohorts_groups=None,
    titles_abscissa_groups=None,
    aspect=None,
    orientation_box=None,
    axis_linear_minimum=None,
    colors_groups=None,
    records_cohorts=None,
):
    """
    Drive creation of box plots for comparison of a quantitative variable across
    categorical groups.

    arguments:
        name_figure (str): name of figure
        variable (str): name of quantitative continuous variable for values on
            plot
        title_ordinate (str): title for ordinate or vertical axis
        title_abscissa (str): title for abscissa or horizontal axis
        cohorts_groups (list<str>): names of stratification cohorts from which
            to collect values of variable for groups on plot
        titles_abscissa_groups (list<str>): titles for groups on abscissa or
            horizontal axis in same sort order as 'cohorts_groups'
        aspect (str): orientation and aspect ratio of figure: either 'portrait',
            'portrait_half_width', 'landscape', or 'landscape_half_height'
        orientation_box (str): whether the orientation of boxes is 'horizontal'
            or 'vertical'
        axis_linear_minimum (float): minimal value for range of linear axis
        colors_groups (list<tuple>): color parameters for boxes of groups in
            same sort order as 'cohorts_groups'
        records_cohorts (list<dict>): records with information about cohorts

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Organize records and entries for cohorts.
    entries_cohorts = (
        utility.structure_from_records_to_entries(
            records=records_cohorts,
    ))
    # Collect values of variable for each group.
    values_groups = list()
    # Iterate on cohorts for groups in figure.
    for cohort_group in cohorts_groups:
        # Access cohort table for current group.
        # Copy information.
        table_cohort = (entries_cohorts[cohort_group]["table"]).copy(deep=True)
        #print(cohort_group)
        #print(table_cohort)
        # Extract values of variable for cohort.
        values = copy.deepcopy(table_cohort[variable].dropna().to_numpy())
        # Collect values for variable for group cohort.
        values_groups.append(values)
        pass
    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()
    # Create figure.
    figure = plot.plot_boxes_groups(
        values_groups=values_groups,
        title_ordinate=title_ordinate,
        title_abscissa=title_abscissa,
        titles_abscissa_groups=titles_abscissa_groups,
        colors_groups=colors_groups,
        label_top_center="",
        label_top_left="", # name_figure
        label_top_right="",
        aspect=aspect,
        orientation_box=orientation_box,
        axis_linear_minimum=axis_linear_minimum,
        fonts=fonts,
        colors=colors,
    )
    # Return information.
    return figure


def drive_create_plots_box(
    records_cohorts=None,
    report=None,
):
    """
    Drive creation of box plots for comparison of a quantitative variable across
    categorical groups.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Organize parameters for plotting figures.
    records_parameters = define_parameters_for_plots_groups_box()
    # Collect information for figures.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record_parameter in records_parameters:
        # Create box plot for a specific cohort, specific groups within that
        # cohort, and other specific parameters.
        pail[record_parameter["name_figure"]] = create_plot_box(
            name_figure=record_parameter["name_figure"],
            variable=record_parameter["variable"],
            title_ordinate=record_parameter["title_ordinate"],
            title_abscissa=record_parameter["title_abscissa"],
            cohorts_groups=record_parameter["cohorts_groups"],
            titles_abscissa_groups=record_parameter["titles_abscissa_groups"],
            aspect=record_parameter["aspect"],
            orientation_box=record_parameter["orientation_box"],
            axis_linear_minimum=0.0,
            colors_groups=record_parameter["colors_groups"],
            records_cohorts=records_cohorts,
        )
        pass

    # Return information.
    return pail


def plot_variable_means_dot_trajectory(
    label_title=None,
    column_phenotype=None,
    column_trajectory=None,
    threshold_trajectory=None,
    title_abscissa=None,
    title_ordinate=None,
    table=None,
):
    """
    Plots charts from the analysis process.

    The trajectory variable must have discrete values.

    arguments:
        label_title (str): label title for plot
        column_phenotype (str): name of column in table for continuous phenotype
        column_trajectory (str): name of column in table for trajectory variable
        threshold_trajectory (int): maximal value of trajectory variable
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    # Organize information for plot.
    # Copy information.
    table = table.copy(deep=True)
    # Select relevant information.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    columns = [
        #"eid",
        "IID",
        column_trajectory, column_phenotype
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table.loc[
        (table[column_trajectory] < threshold_trajectory), :
    ]
    # Aggregate phenotype values by day.
    table.set_index(
        [column_trajectory],
        append=False,
        drop=True,
        inplace=True
    )
    groups = table.groupby(level=[column_trajectory], axis="index",)
    #groups.aggregate(
    #    mean=(column_phenotype, numpy.nanmean),
    #    deviation=(column_phenotype, numpy.nanstd)
    #)
    records = list()
    for name, group in groups:
        # Initialize missing values.
        days = float("nan")
        count = float("nan")
        mean = float("nan")
        standard_error = float("nan")
        confidence_95 = float("nan")
        confidence_95_low = float("nan")
        confidence_95_high = float("nan")

        table_group = group.reset_index(
            level=None,
            inplace=False,
            drop=False,
        )
        trajectory_units = table_group[column_trajectory].dropna().to_list()[0]
        array = copy.deepcopy(table_group[column_phenotype].dropna().to_numpy())
        # Determine count of valid values.
        count = int(array.size)
        if (count > 5):
            # Determine summary statistics for values in array.
            #mean = round(numpy.nanmean(array), 3)
            mean = numpy.nanmean(array)
            standard_error = scipy.stats.sem(array)
            confidence_95 = (1.96 * standard_error)
            confidence_95_low = (mean - confidence_95)
            confidence_95_high = (mean + confidence_95)
            pass
        # Collect information.
        record = dict()
        record["trajectory_units"] = trajectory_units
        record["count"] = count
        record["mean"] = mean
        record["error"] = standard_error
        record["confidence_95"] = confidence_95
        record["confidence_95_low"] = confidence_95_low
        record["confidence_95_high"] = confidence_95_high
        records.append(record)
        pass
    table_aggregate = pandas.DataFrame(data=records)

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Create figure.
    figure = plot.plot_scatter_points_discrete_abscissa_ordinate_error_bars(
        table=table_aggregate,
        abscissa="trajectory_units",
        ordinate="mean",
        ordinate_error_low="confidence_95",
        ordinate_error_high="confidence_95",
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        fonts=fonts,
        colors=colors,
        size_marker=25,
        label_title=label_title,
    )
    # Return.
    return figure



def organize_phenotypes_plots_dot_trajectory_month(
    table=None,
):
    """
    Organizes information and plots.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Copy information.
    table = table.copy(deep=True)
    # Prepare table to summarize phenotype variables across cohorts and models.
    # These cohorts and models are simple and do not include multiple covariates
    # for genetic analyses.
    # Collect records of information about each cohort and model.
    pails_cohorts = ukb_strat.stratify_phenotype_cohorts_set_sex_age_menopause(
        table=table,
    )

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for pail_cohort in pails_cohorts:
        # Define phenotypes.
        phenotypes = [
            "cholesterol", "cholesterol_imputation",
            "oestradiol", "oestradiol_imputation",
            "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
            "oestradiol_free", "oestradiol_free_imputation",
            "testosterone", "testosterone_imputation",
            "testosterone_bioavailable",
            "testosterone_bioavailable_imputation",
            "testosterone_free", "testosterone_free_imputation",
            "vitamin_d", "vitamin_d_imputation",
            "steroid_globulin", "steroid_globulin_imputation",
            "albumin", "albumin_imputation",
            "alcohol_frequency",
            "alcohol_drinks_monthly",
            "alcohol_audit",
            "alcohol_auditc",
            "alcohol_auditp",
        ]
        # Iterate on phenotypes.
        for phenotype in phenotypes:
            # Plot hormone trajectories across month of assessment.
            # Chart.
            name_label = str(str(pail_cohort["name"]) + "_" + str(phenotype))
            name_plot = str(name_label + "_dot_month")
            pail[name_plot] = plot_variable_means_dot_trajectory(
                label_title=name_label,
                column_phenotype=phenotype,
                column_trajectory="month_order",
                threshold_trajectory=13,
                title_abscissa="month of assessment and blood draw",
                title_ordinate="mean concentration (95% C.I.)",
                table=pail_cohort["table"],
            )
            pass
        pass
    # Return information.
    return pail


def organize_phenotypes_parameters_for_plots_dot_trajectory_menstruation():
    """
    Organizes parameters for histogram plots.

    arguments:

    raises:

    returns:
        (dict<object>): collection of parameters for plots

    """

    # Collect records of information about each cohort and model.
    records = list()

    record = dict()
    record["name"] = "albumin"
    record["variable"] = "albumin"
    records.append(record)

    record = dict()
    record["name"] = "albumin_imputation"
    record["variable"] = "albumin_imputation"
    records.append(record)

    record = dict()
    record["name"] = "steroid_globulin"
    record["variable"] = "steroid_globulin"
    records.append(record)

    record = dict()
    record["name"] = "steroid_globulin_imputation"
    record["variable"] = "steroid_globulin_imputation"
    records.append(record)

    record = dict()
    record["name"] = "vitamin_d"
    record["variable"] = "vitamin_d"
    records.append(record)

    record = dict()
    record["name"] = "vitamin_d_imputation"
    record["variable"] = "vitamin_d_imputation"
    records.append(record)

    record = dict()
    record["name"] = "oestradiol"
    record["variable"] = "oestradiol"
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_imputation"
    record["variable"] = "oestradiol_imputation"
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_bioavailable_imputation"
    record["variable"] = "oestradiol_bioavailable_imputation"
    records.append(record)

    record = dict()
    record["name"] = "oestradiol_free_imputation"
    record["variable"] = "oestradiol_free_imputation"
    records.append(record)

    record = dict()
    record["name"] = "testosterone"
    record["variable"] = "testosterone"
    records.append(record)

    record = dict()
    record["name"] = "testosterone_imputation"
    record["variable"] = "testosterone_imputation"
    records.append(record)

    record = dict()
    record["name"] = "testosterone_bioavailable_imputation"
    record["variable"] = "testosterone_bioavailable_imputation"
    records.append(record)

    record = dict()
    record["name"] = "testosterone_free_imputation"
    record["variable"] = "testosterone_free_imputation"
    records.append(record)

    # Return information
    return records


def create_plots_dot_trajectory_menstruation(
    records_cohorts=None,
    report=None,
):
    """
    Organizes description histogram plots.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plot objects

    """

    # Organize phenotypes and parameters for plots.
    records_phenotypes = (
        organize_phenotypes_parameters_for_plots_dot_trajectory_menstruation()
    )

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record_cohort in records_cohorts:
        # Iterate on phenotypes.
        for record_phenotype in records_phenotypes:
            # Organize information for plot.
            name = record_phenotype["name"]
            variable = record_phenotype["variable"]
            # Copy information.
            table_cohort = record_cohort["table"].copy(deep=True)
            # Dot trajectory plot.
            name_label = str(str(record_cohort["name"]) + "_" + str(name))
            name_plot = str(name_label + "_dot_menstruation")
            pail[name_plot] = plot_variable_means_dot_trajectory(
                label_title="", # name_label
                column_phenotype=variable,
                column_trajectory="menstruation_days",
                threshold_trajectory=31,
                title_abscissa="", # "days of menstrual cycle",
                title_ordinate="", # "mean concentration (95% C.I.)"
                table=table_cohort,
            )
            pass
        pass
    # Return information.
    return pail


##########
# Management of plots on phenotype variables within stratification cohorts
##########


def prepare_phenotype_variables_in_stratification_cohorts(
    set_plots=None,
    table=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_cohorts (str): name of the set of cohorts for which to create
            description tables; "genotype", "phenotype", or "read"
        set_plots (list<str>): names of plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Prepare phenotype variables in stratification cohorts specific to each
    # type of plot.
    if ("histogram" in set_plots):
        # Histograms.
        # Stratify records within separate tables for cohorts.
        records_cohorts = (
            ukb_strat.drive_stratify_phenotype_cohorts_set_description_tables(
                table=table,
        ))
        # Filter relevant cohorts.
        names_cohorts = [
            "self_white_female_male",
            "self_white_female",
            "self_white_female_menstruation_regular",
            "self_white_female_premenopause",
            "self_white_female_perimenopause",
            "self_white_female_postmenopause",
            "self_white_male",
            "self_white_male_age_low",
            "self_white_male_age_middle",
            "self_white_male_age_high",
            #"self_white_alcohol_current_female_male",
            #"self_white_alcohol_current_female",
            #"self_white_alcohol_current_female_menstruation_regular",
            #"self_white_alcohol_current_female_premenopause",
            #"self_white_alcohol_current_female_perimenopause",
            #"self_white_alcohol_current_female_postmenopause",
            #"self_white_alcohol_current_male",
            #"self_white_alcohol_current_male_age_low",
            #"self_white_alcohol_current_male_age_middle",
            #"self_white_alcohol_current_male_age_high",
        ]
        records_cohorts = utility.filter_records_by_name(
            names=names_cohorts,
            records=records_cohorts,
            report=True,
        )
        # Apply Distribution Scale Transformations to variables of interest in
        # each cohort.
        records_cohorts = (
            scale.drive_transformations_on_multiple_variables_in_cohorts(
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
                    "steroid_globulin_imputation",
                    "albumin_imputation",
                ],
                records_cohorts=records_cohorts,
                report=True,
        ))
    elif ("box" in set_plots):
        # Box plots for groups.
        records_cohorts = (
            ukb_strat.drive_stratify_phenotype_cohorts_set_description_tables(
                table=table,
        ))
    elif ("dot_trajectory_menstruation" in set_plots):
        # Dot trajectory plots of menstrual cycle.
        records_cohorts = (
            ukb_strat.drive_stratify_phenotype_cohorts_set_female_menstruation(
                table=table,
        ))
    else:
        records_cohorts = []
        pass
    # Return information
    return records_cohorts


def create_plots_for_phenotype_variables_in_cohorts(
    set_cohorts=None,
    set_plots=None,
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_cohorts (str): name of the set of cohorts for which to create
            description tables; "genotype", "phenotype", or "read"
        set_plots (list<str>): names of plots to create
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
        records_cohorts = prepare_phenotype_variables_in_stratification_cohorts(
            set_plots=set_plots,
            table=source["table_phenotypes"],
            report=report,
        )
        pass

    # Create description plots across cohorts.
    # Histograms.
    if ("histogram" in set_plots):
        pail_histogram = drive_create_plots_histogram(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        pail_histogram = dict()
        pass
    # Box plots for groups.
    if ("box" in set_plots):
        pail_box = drive_create_plots_box(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        pail_box = dict()
        pass
    # Dot trajectory plots of menstrual cycle.
    if ("dot_trajectory_menstruation" in set_plots):
        pail_dot_menstruation = create_plots_dot_trajectory_menstruation(
            records_cohorts=records_cohorts,
            report=report,
        )
    else:
        pail_dot_menstruation = dict()
        pass

    # Dot plots of phenotypes ordinal variables.
    #pail["dot_trajectory_month"] = (
    #    organize_phenotypes_plots_dot_trajectory_month(
    #        table=table,
    #))

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_description_tables()")
        utility.print_terminal_partition(level=3)

    # Collect information.
    pail_write = dict()
    pail_write["histogram"] = pail_histogram
    pail_write["box"] = pail_box
    pail_write["dot_trajectory_menstruation"] = pail_dot_menstruation
    # Write product information to file.
    plot.write_product_plots_child_directories(
        pail_write=pail_write,
        path_parent=paths["plot"],
    )
    pass



##########
# Management of plots on sumaries of genetic correlation analyses
##########


# TODO: TCW; 6 October 2022
# TODO: Will need to update for newer analyses...
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
                "cohort": "string",
                "phenotype_primary": "string",
                "dependence_type": "string",
                "model_adjustment": "string",
                "model_context": "string",
                "variable": "string",
                "correlation": "float32",
                "standard_error": "float32",
                "interval_95": "float32",
                "confidence_95_low": "float32",
                "confidence_95_high": "float32",
                "probability": "float32",
            },
            report=report,
    ))
    # Return information.
    return pail


def define_parameters_correlation_summaries():
    """
    Defines parameters for organization of information from regression
    summaries.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about plots for regression
            summaries

    """

    # Collect records of information about regression summaries and their
    # Forest Plots.
    records = list()

    record = dict()
    record["name"] = "alcohol_dependence"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.6
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "alcohol_dependence_genotype"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.5
    record["abscissa_maximum"] = 0.5
    records.append(record)

    record = dict()
    record["name"] = "alcohol_quantity"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.24
    record["abscissa_maximum"] = 0.24
    records.append(record)

    record = dict()
    record["name"] = "alcohol_quantity_no_ukb"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.24
    record["abscissa_maximum"] = 0.24
    records.append(record)

    if False:

        record = dict()
        record["name"] = "depression"
        record["regression_type"] = "logistic"
        record["abscissa_minimum"] = -0.5
        record["abscissa_maximum"] = 0.5
        records.append(record)

        record = dict()
        record["name"] = "schizophrenia"
        record["regression_type"] = "logistic"
        record["abscissa_minimum"] = -0.5
        record["abscissa_maximum"] = 0.5
        records.append(record)

        record = dict()
        record["name"] = "bipolar"
        record["regression_type"] = "logistic"
        record["abscissa_minimum"] = -0.5
        record["abscissa_maximum"] = 0.5
        records.append(record)

        record = dict()
        record["name"] = "bipolar_type_1"
        record["regression_type"] = "logistic"
        record["abscissa_minimum"] = -0.5
        record["abscissa_maximum"] = 0.5
        records.append(record)

        record = dict()
        record["name"] = "bipolar_type_2"
        record["regression_type"] = "logistic"
        record["abscissa_minimum"] = -0.5
        record["abscissa_maximum"] = 0.5
        records.append(record)

        pass


    # Return information
    return records


def organize_correlation_summary_tables_for_forest_plots(
    pail_correlation_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Organize the format of regression summary tables for creation of Forest
    Plots.

    arguments:
        pail_correlation_tables (dict<object>): collection of Pandas data-frame
            tables for correlation summaries
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of Pandas data-frame tables for
            plots

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            pro_reg.organize_regression_summary_table_for_forest_plots(
                type=record_parameter["regression_type"],
                model_contexts=["joint"],
                model_adjustments=["adjust", "unadjust",],
                column_variable="variable",
                variables=[
                    "oestradiol_detection",
                    "oestradiol_bioavailable_imputation_log",
                    "oestradiol_free_imputation_log",
                    "testosterone_imputation_log",
                    "testosterone_bioavailable_imputation_log",
                    "testosterone_free_imputation_log",
                    "steroid_globulin_imputation_log",
                    "albumin_imputation",
                ],
                columns_translations={
                    "model_adjustment": "group",
                    "variable": "category",
                    "correlation": "value",
                    #"interval_99": "interval_above",
                    "interval_99": "interval_below",
                },
                labels_categories={
                    "oestradiol_detection": "EST-T",
                    "oestradiol_bioavailable_imputation_log": "EST-B",
                    "oestradiol_free_imputation_log": "EST-F",
                    "testosterone_imputation_log": "TST-T",
                    "testosterone_bioavailable_imputation_log": "TST-B",
                    "testosterone_free_imputation_log": "TST-F",
                    "steroid_globulin_imputation_log": "SHBG",
                    "albumin_imputation": "ALBU",
                },
                sorts_categories={
                    "oestradiol_detection": 1,
                    "oestradiol_bioavailable_imputation_log": 2,
                    "oestradiol_free_imputation_log": 3,
                    "testosterone_imputation_log": 4,
                    "testosterone_bioavailable_imputation_log": 5,
                    "testosterone_free_imputation_log": 6,
                    "steroid_globulin_imputation_log": 7,
                    "albumin_imputation": 8,
                },
                column_stratification="cohort",
                table=pail_correlation_tables[name_table],
                report=report,
        ))
        pass
    # Return information.
    return pail


def create_correlation_summary_forest_plots(
    pail_plot_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        pail_plot_tables (dict<dict<object>>): collection tree of Pandas
            data-frame tables for plots
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of figures

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            plot.drive_iterate_plot_forest_two_groups(
                pail_tables=pail_plot_tables[record_parameter["name"]],
                column_group="group",
                column_ordinate_label="category_label",
                column_ordinate_sort="category_sort",
                column_abscissa_value="value",
                column_abscissa_interval_below="interval_below",
                column_abscissa_interval_above="interval_above",
                group_one="adjust", # group one markers: circles above center
                group_two="unadjust", # group two markers: triangles below
                abscissa_minimum=record_parameter["abscissa_minimum"],
                abscissa_maximum=record_parameter["abscissa_maximum"],
                ordinate_title="", # "Sex Hormone or Protein"
                abscissa_title="", # "Marginal Model Coefficient (95% C.I.)"
                print_label_on_chart=False,
                label_chart_prefix=str(record_parameter["name"]),
                label_size_ordinate_categories="one",
                label_size_abscissa_values="one",
                size_marker_one=70,
                size_marker_two=35,
                color_marker_one="blue_navy", # "blue_sky", "purple", "magenta"
                color_marker_two="orange", # "green"
                space_groups=0.15, # vertical space between groups' markers
        ))
        pass
    # Return information.
    return pail





def read_organize_correlation_summaries_create_forest_plots(
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_plots (list<str>): names of plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plots

    """

    # Read regression summary tables from file.
    pail_correlation_tables = read_source_correlation_summary_tables(
        path_dock=paths["dock"],
        report=True,
    )
    # Define parameters for regression summaries.
    records_parameters = define_parameters_correlation_summaries()
    # Organize tables for charts.
    pail_plot_tables = organize_correlation_summary_tables_for_forest_plots(
        pail_correlation_tables=pail_correlation_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Create Forest Plots.
    pail_plots = create_correlation_summary_forest_plots(
        pail_plot_tables=pail_plot_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Return information.
    return pail_plots


def create_plots_for_analysis_summaries(
    set_plots=None,
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_plots (list<str>): names of plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Forest Plots for summary tables from phenotype regressions.
    if ("forest_regressions" in set_plots):
        # Read regression summary tables from file.
        # Organize information from regression summary tables.
        # Create Forest Plots.
        pail_forest_regressions = (
            read_organize_regression_summaries_create_forest_plots(
                paths=paths,
                report=True,
            )
        )
    else:
        pail_forest_regressions = dict()
        pass

    # Forest Plots for summary tables from genetic correlations.
    if ("forest_correlations" in set_plots):
        # Read regression summary tables from file.
        # Organize information from regression summary tables.
        # Create Forest Plots.
        pail_forest_correlations = (
            read_organize_correlation_summaries_create_forest_plots(
                paths=paths,
                report=True,
            )
        )
    else:
        pail_forest_correlations = dict()
        pass

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: create_plots_for_analysis_summaries()")
        utility.print_terminal_partition(level=3)

    # Collect information.
    pail_write = dict()
    pail_write["forest_regressions"] = pail_forest_regressions
    pail_write["forest_correlations"] = pail_forest_correlations
    # Write product information to file.
    plot.write_product_plots_child_child_directories(
        pail_write=pail_write,
        path_parent=paths["plot"],
    )
    pass



##########
# Write
##########



###############################################################################
# Drivers
# The purpose of these driver functions is to make the module's functionality
# accessible in modular parts.

# TODO: TCW; 01 June 2022
# TODO: Eventually, adapt this function so that it can create plots for cohors
# TODO: from multiple sources, including genotype GWAS cohorts read from file.


def define_parameters_regression_summaries():
    """
    Defines parameters for organization of information from regression
    summaries.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about plots for regression
            summaries

    """

    # Collect records of information about regression summaries and their
    # Forest Plots.
    records = list()

    record = dict()
    record["name"] = "alcohol_dependence"
    record["regression_type"] = "logistic"
    record["abscissa_minimum"] = -0.6
    record["abscissa_maximum"] = 0.6
    records.append(record)

    record = dict()
    record["name"] = "alcohol_frequency"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.2
    record["abscissa_maximum"] = 0.2
    records.append(record)

    record = dict()
    record["name"] = "alcohol_auditc"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.30
    record["abscissa_maximum"] = 0.24
    records.append(record)

    record = dict()
    record["name"] = "alcohol_quantity"
    record["regression_type"] = "linear"
    record["abscissa_minimum"] = -0.149
    record["abscissa_maximum"] = 0.149
    records.append(record)

    # Return information
    return records


def organize_regression_summary_tables_for_forest_plots(
    pail_regression_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Organize the format of regression summary tables for creation of Forest
    Plots.

    arguments:
        pail_regression_tables (dict<object>): collection of Pandas data-frame
            tables for regression summaries
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of Pandas data-frame tables for
            plots

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            pro_reg.organize_regression_summary_table_for_forest_plots(
                type=record_parameter["regression_type"],
                model_contexts=["marginal"],
                model_adjustments=["adjust", "unadjust",],
                column_variable="variable",
                variables=[
                    "oestradiol_imputation",
                    "oestradiol_bioavailable_imputation",
                    "oestradiol_free_imputation",
                    "testosterone_imputation",
                    "testosterone_bioavailable_imputation",
                    "testosterone_free_imputation",
                    "steroid_globulin_imputation",
                    "albumin_imputation",
                ],
                columns_translations={
                    "model_adjustment": "group",
                    "variable": "category",
                    "parameter": "value",
                    "interval_95": "interval_below",
                },
                labels_categories={
                    "oestradiol_imputation": "EST-T",
                    "oestradiol_bioavailable_imputation": "EST-B",
                    "oestradiol_free_imputation": "EST-F",
                    "testosterone_imputation": "TST-T",
                    "testosterone_bioavailable_imputation": "TST-B",
                    "testosterone_free_imputation": "TST-F",
                    "steroid_globulin_imputation": "SHBG",
                    "albumin_imputation": "ALBU",
                },
                sorts_categories={
                    "oestradiol_imputation": 1,
                    "oestradiol_bioavailable_imputation": 2,
                    "oestradiol_free_imputation": 3,
                    "testosterone_imputation": 4,
                    "testosterone_bioavailable_imputation": 5,
                    "testosterone_free_imputation": 6,
                    "steroid_globulin_imputation": 7,
                    "albumin_imputation": 8,
                },
                column_stratification="cohort",
                table=pail_regression_tables[name_table],
                report=report,
        ))
        pass
    # Return information.
    return pail


def create_regression_summary_forest_plots(
    pail_plot_tables=None,
    records_parameters=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        pail_plot_tables (dict<dict<object>>): collection tree of Pandas
            data-frame tables for plots
        records_parameters (list<dict>): records with information about plots
            for regression summaries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<dict<object>>): collection tree of figures

    """

    # Collect tables for plots.
    pail = dict()
    # Iterate on phenotypes.
    for record_parameter in records_parameters:
        # Organize tables.
        name_table = str("table_" + record_parameter["name"])
        pail[record_parameter["name"]] = (
            plot.drive_iterate_plot_forest_two_groups(
                pail_tables=pail_plot_tables[record_parameter["name"]],
                column_group="group",
                column_ordinate_label="category_label",
                column_ordinate_sort="category_sort",
                column_abscissa_value="value",
                column_abscissa_interval_below="interval_below",
                column_abscissa_interval_above="interval_above",
                group_one="adjust", # group one markers: circles above center
                group_two="unadjust", # group two markers: triangles below
                abscissa_minimum=record_parameter["abscissa_minimum"],
                abscissa_maximum=record_parameter["abscissa_maximum"],
                ordinate_title="", # "Sex Hormone or Protein"
                abscissa_title="", # "Marginal Model Coefficient (95% C.I.)"
                print_label_on_chart=False,
                label_chart_prefix=str(record_parameter["name"]),
                label_size_ordinate_categories="one",
                label_size_abscissa_values="one",
                size_marker_one=70,
                size_marker_two=35,
                color_marker_one="purple",
                color_marker_two="orange", # "green"
                space_groups=0.15, # vertical space between groups' markers
        ))
        pass
    # Return information.
    return pail


def read_organize_regression_summaries_create_forest_plots(
    paths=None,
    report=None,
):
    """
    Drive creation of descriptive charts about phenotypes across cohorts from
    the UK Biobank.

    arguments:
        set_plots (list<str>): names of plots to create
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of plots

    """

    # Read regression summary tables from file.
    pail_regression_tables = read_source_regression_summary_tables(
        path_dock=paths["dock"],
        report=True,
    )
    # Define parameters for regression summaries.
    records_parameters = define_parameters_regression_summaries()
    # Organize tables for charts.
    pail_plot_tables = organize_regression_summary_tables_for_forest_plots(
        pail_regression_tables=pail_regression_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Create Forest Plots.
    pail_plots = create_regression_summary_forest_plots(
        pail_plot_tables=pail_plot_tables,
        records_parameters=records_parameters,
        report=report,
    )
    # Return information.
    return pail_plots


################################################################################
# Procedure

# TODO: TCW; 01 June 2022
# TODO: activate description plots for cohorts
# TODO: trajectory plots of estradiol, testosterone (total, bioavailable, free)
# TODO: across female-premenopause and female-menstruation-regular cohorts
# TODO: consider excluding persons who used hormone-altering contraception or therapy


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
    print("version check: TCW, 6 October 2022")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )

    # Create new description plots from phenotype tables for cohorts.
    create_plots_for_phenotype_variables_in_cohorts(
        set_cohorts="phenotype", # could change this...
        set_plots=[], # ["box", "histogram", "dot_trajectory_menstruation",]
        paths=paths,
        report=True,
    )
    # Create new description plots from summary tables for previous analyses.
    create_plots_for_analysis_summaries(
        set_plots=["forest_correlations",], # ["forest_regressions", "forest_correlations",],
        paths=paths,
        report=True,
    )

    pass



#
