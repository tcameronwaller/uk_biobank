"""
Describe variables from data of the U.K. Biobank.

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
import stratification

###############################################################################
# Functionality


##########
# Cohort, model, phenotype description


def organize_report_contingency_table_stratification_by_missingness(
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
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)

    # Organize contingency table descriptions of missingness.

    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="sex_text",
        stratifications=["female", "male"],
        column_missingness="testosterone",
        table=table,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="sex_text",
        stratifications=["female", "male"],
        column_missingness="oestradiol",
        table=table,
        report=report,
    )
    # Filter to females who were not pregnant.
    table_female = table.loc[
        (table["sex_text"] == "female"), :
    ]
    table_female_not_pregnant = table_female.loc[
        (table_female["pregnancy"] == 0), :
    ]
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 1], # premenopause versus perimenopause
        column_missingness="testosterone",
        table=table_female_not_pregnant,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 1], # premenopause versus perimenopause
        column_missingness="oestradiol",
        table=table_female_not_pregnant,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 2], # premenopause versus postmenopause
        column_missingness="testosterone",
        table=table_female_not_pregnant,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="menopause_ordinal",
        stratifications=[0, 2], # premenopause versus postmenopause
        column_missingness="oestradiol",
        table=table_female_not_pregnant,
        report=report,
    )

    # TODO: also consider missingness in "younger" and "older" males
    table_male = table.loc[
        (table["sex_text"] == "male"), :
    ]
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="age_grade_male",
        stratifications=[0, 2], # young versus old
        column_missingness="testosterone",
        table=table_male,
        report=report,
    )
    utility.report_contingency_table_stratification_by_missingness(
        column_stratification="age_grade_male",
        stratifications=[0, 2], # young versus old
        column_missingness="oestradiol",
        table=table_male,
        report=report,
    )

    pass


# TODO: I need to test this...
def organize_text_report_column_pair_correlations(
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
        (str): textual report

    """

    table = table.copy(deep=True)
    table.dropna(
        axis="index",
        how="any",
        subset=[column_one, column_two],
        inplace=True,
    )
    count_rows = table.shape[0]
    pearson_correlation, pearson_probability = scipy.stats.pearsonr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    spearman_correlation, spearman_probability = scipy.stats.spearmanr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    # Report.
    report_text = str()
    report_text += textwrap.dedent(
        """\

            --------------------------------------------------
            column one: {column_one}
            column two: {column_two}
            non-missing records: {count_rows}
            ----------
            Correlations
            ----------
            Pearson: {pearson_correlation} (p: {pearson_probability})
            Spearman: {spearman_correlation} (p: {spearman_probability})
            --------------------------------------------------
        """
    ).format(
        column_one=column_one,
        column_two=column_two,
        count_rows=count_rows,
        pearson_correlation=round(pearson_correlation, 3),
        pearson_probability=round(pearson_probability, 5),
        spearman_correlation=round(spearman_correlation, 3),
        spearman_probability=round(spearman_probability, 5),
    )
    # Return information.
    return report_text


# TODO: work on these text reports...
def organize_text_report_phenotype_correlations(
    table=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    # Copy information.
    table = table.copy(deep=True)
    # Organize textual report.
    report_text = str()
    report_text += organize_text_report_column_pair_correlations(
        column_one="testosterone",
        column_two="steroid_globulin",
        table=None,
    )

    # TODO: okay so now I need to append that new text report to the collection text string report...


    # Return information.
    return report_text


def organize_report_cohort_model_variables_summaries_record(
    name=None,
    cohort_model=None,
    category=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        name (str): name for cohort, model, and phenotype
        cohort_model (str): name for cohort and model
        category (str): name of category for cohort and model to facilitate
            table sorts and comparisons
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for cohort.
    record = dict()
    record["name"] = str(name)
    record["cohort_model"] = str(cohort_model)
    record["category"] = str(category)
    record["cohort_count"] = int(table.shape[0])
    # Collect information for general columns.
    columns = [
        "age",
        "body_mass_index", "body_mass_index_log",
        "menstruation_days", "menstruation_days_threshold",
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

    # Iterate on relevant columns.
    # Collect information for record.
    for column in columns:
        # Initialize missing values.
        count = float("nan")
        mean = float("nan")
        standard_error = float("nan")
        confidence_95_low = float("nan")
        confidence_95_high = float("nan")
        median = float("nan")
        standard_deviation = float("nan")
        minimum = float("nan")
        maximum = float("nan")
        # Determine whether table has the column.
        if (column in table.columns.to_list()):
            array = copy.deepcopy(table[column].dropna().to_numpy())
            # Determine count of valid values.
            count = int(array.size)
            if (count > 10):
                # Determine mean, median, standard deviation, and standard error of
                # values in array.
                mean = numpy.nanmean(array)
                standard_error = scipy.stats.sem(array)
                confidence_95_low = (mean - (1.96 * standard_error))
                confidence_95_high = (mean + (1.96 * standard_error))
                median = numpy.nanmedian(array)
                standard_deviation = numpy.nanstd(array)
                minimum = numpy.nanmin(array)
                maximum = numpy.nanmax(array)
                pass
            pass
        # Collect information for record.
        record[str(column + "_count")] = str(count)
        record[str(column + "_mean")] = str(round(mean, 3))
        record[str(column + "_stderr")] = str(round(standard_error, 3))
        record[str(column + "_confidence_95")] = str(
            str(round(confidence_95_low, 3)) + " ... " +
            str(round(confidence_95_high, 3))
        )
        record[str(column + "_median")] = str(round(median, 3))
        record[str(column + "_stdev")] = str(round(standard_deviation, 3))
        record[str(column + "_min")] = str(round(minimum, 3))
        record[str(column + "_max")] = str(round(maximum, 3))
        pass
    # Return information.
    return record



# TODO: TCW 30 July 2021
# TODO: work on this... report separate percentages for reportability LOW and reportability HIGH

def organize_cohort_hormone_missingness_record(
    name_cohort=None,
    name_hormone=None,
    column_hormone=None,
    column_reportability_limit=None,
    column_missingness_range=None,
    table=None,
):
    """
    Determine for each hormone whether missingness variable indicates that the
    measurement was not reportable because it was beyond the detection range.

    # cohort:
    # hormone:
    # total measurement missingness: <count> (<percentage>%)
    # missing with reason "above or below reportable limit": <count> (<percentage>%)
    # missing with reportability "not reportable ... too low": <count> (<percentage>%)
    # missing with both annotations: <count> (<percentage>%)

    # "<hormone>_reportability_limit" 1: not reportable due to measurement less than limit of detection
    # "<hormone>_missingness_range" 1: missing due to measurement beyond detection range

    arguments:
        name_cohort (str): name of cohort
        name_hormone (str): name of hormone
        column_hormone (str): name of table's column for hormone's measurements
        column_reportability_limit (str): name of table's column for hormone's
            reportability less than limit of detection
        column_missingness_range (str): name of table's column for hormone's
            reason for missingness beyond detection range
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["hormone"] = str(name_hormone)
    # Copy information.
    table = table.copy(deep=True)

    # Select relevant rows of the table.
    #array_measurement_total = copy.deepcopy(table[name].to_numpy())
    #array_measurement_valid = copy.deepcopy(table[name].dropna().to_numpy())
    table_measurement_valid = table.dropna(
        axis="index",
        how="any",
        subset=[column_hormone],
        inplace=False,
    )
    table_measurement_missing = table[table[column_hormone].isna()]
    table_missingness_range = table_measurement_missing.loc[
        (
            (table_measurement_missing[column_missingness_range] == 1)
        ), :
    ]
    table_reportability_limit = table_measurement_missing.loc[
        (
            (table_measurement_missing[column_reportability_limit] == 1)
        ), :
    ]
    table_both = table_measurement_missing.loc[
        (
            (table_measurement_missing[column_missingness_range] == 1) &
            (table_measurement_missing[column_reportability_limit] == 1)
        ), :
    ]

    # Count records.
    #count_measurement_total = int(array_measurement_total.size)
    #count_measurement_valid = int(array_measurement_valid.size)
    count_cohort_total = table.shape[0]
    count_measurement_valid = table_measurement_valid.shape[0]
    count_measurement_missing = table_measurement_missing.shape[0]
    count_missingness_range = table_missingness_range.shape[0]
    count_reportability_limit = table_reportability_limit.shape[0]
    count_both = table_both.shape[0]

    # Calculate percentages.
    percentage_measurement_valid = round(
        ((count_measurement_valid / count_cohort_total) * 100), 3
    )
    percentage_measurement_missing = round(
        ((count_measurement_missing / count_cohort_total) * 100), 3
    )
    percentage_missingness_range = round(
        ((count_missingness_range / count_measurement_missing) * 100), 3
    )
    percentage_reportability_limit = round(
        ((count_reportability_limit / count_measurement_missing) * 100), 3
    )
    percentage_both = round(
        ((count_both / count_measurement_missing) * 100), 3
    )

    # Collect information for record.
    record["count_cohort_samples"] = count_cohort_total
    record["measurement_valid"] = str(
        str(count_measurement_valid) +
        " (" + str(percentage_measurement_valid) + "%)"
    )
    record["measurement_missingness"] = str(
        str(count_measurement_missing) +
        " (" + str(percentage_measurement_missing) + "%)"
    )
    record["missingness_range"] = str(
        str(count_missingness_range) +
        " (" + str(percentage_missingness_range) + "%)"
    )
    record["reportability_low"] = str(
        str(count_reportability_limit) +
        " (" + str(percentage_reportability_limit) + "%)"
    )
    record["both_missingness_reportability"] = str(
        str(count_both) +
        " (" + str(percentage_both) + "%)"
    )
    # Return information.
    return record


def organize_cohorts_phenotypes_hormones_missingness(
    table=None,
):
    """
    Organizes a summary table about missingness of hormone measurements in
    cohorts.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Copy information.
    table = table.copy(deep=True)
    # Define hormones.
    hormones = list()
    hormones.append("albumin")
    hormones.append("steroid_globulin")
    hormones.append("testosterone")
    hormones.append("oestradiol")
    hormones.append("vitamin_d")
    # Define cohorts for description.
    pail_phenotypes = select_organize_cohorts_models_phenotypes_sets(
        table=table,
    )
    # Collect summary records and construct table.
    records = list()
    # Iterate on cohorts.
    for collection_cohort in pail_phenotypes:
        # Iterate on hormones.
        for hormone in hormones:
            reportability_limit = str(
                str(hormone) + "_reportability_limit"
            )
            missingness_range = str(
                str(hormone) + "_missingness_range"
            )
            # Organize information in record.
            record = organize_cohort_hormone_missingness_record(
                name_cohort=collection_cohort["name"],
                name_hormone=hormone,
                column_hormone=hormone,
                column_reportability_limit=reportability_limit,
                column_missingness_range=missingness_range,
                table=collection_cohort["table"],
            )
            # Collect records.
            records.append(record)
            pass
        pass
    # Organize table.
    table_missingness = pandas.DataFrame(data=records)
    # Return information.
    return table_missingness





##########
# Plot

# TODO: it would be NICE to to organize these plots within sub-directories


def plot_variable_values_histogram(
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
        label_bins="values",
        label_counts="counts per bin",
        fonts=fonts,
        colors=colors,
        line=True,
        line_position=numpy.nanmean(array),
        label_title=label_title, # ""
        label_report=True,
    )
    # Return.
    return figure


# TODO: change name to "ordinal_group" or something like that to differentiate from "categorical_group"

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
        size=15,
        label_title=label_title,
    )
    # Return.
    return figure


def organize_phenotypes_plots_histogram(
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
    records = list()
    records_novel = (
        stratificiation.stratify_set_primary_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record in records:
        # Define phenotypes.
        phenotypes = [
            "age",
            "body_mass_index",
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
        if (record["menstruation"]):
            phenotypes.append("menstruation_days_threshold")
        # Iterate on phenotypes.
        for phenotype in phenotypes:
            # Histogram.
            name_label = str(str(record["name"]) + "_" + str(phenotype))
            name_plot = str(name_label + "_histogram")
            pail[name_plot] = plot_variable_values_histogram(
                label_title=name_label,
                array=record["table"][phenotype].dropna().to_numpy(),
                bins=None,
            )
            pass
        pass
    # Return information.
    return pail


def organize_phenotypes_plots_dot_trajectory(
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
    records = list()
    records_novel = (
        stratification.stratify_set_primary_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record in records:
        # Define phenotypes.
        phenotypes = [
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
        # Iterate on phenotypes.
        for phenotype in phenotypes:
            # Determine whether the current cohort and model table is relevant to
            # the plot of phenotypes across days of menstrual cycle.
            if (record["menstruation"]):
                # Chart.
                name_label = str(str(record["name"]) + "_" + str(phenotype))
                name_plot = str(name_label + "_dot_menstruation")
                pail[name_plot] = plot_variable_means_dot_trajectory(
                    label_title=name_label,
                    column_phenotype=phenotype,
                    column_day="menstruation_days_threshold",
                    threshold_days=45,
                    title_abscissa="days of menstrual cycle",
                    title_ordinate="mean concentration",
                    table=record["table"],
                )
                pass
            pass
        pass
    # Return information.
    return pail


###############################################################################
# Drivers
# The purpose of these driver functions is to make the module's functionality
# accessible in modular parts.


# Descriptions

# TODO: divide this into 2 functions and then decide whether to call each
# TODO: in particular, the descriptions for genotype cohorts are slow (due to kinship)
#

# TODO: also write out "print" reports to text files for Chi-square tests and correlations


def execute_describe_cohorts_models_phenotypes(
    table=None,
    genotype_cohorts=None,
    set=None,
    path_dock=None,
    report=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        genotype_cohorts (bool): whether to generate description for genotype
            cohorts (slow)
        set (str): name of set of cohorts and models to select and organize
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports


    raises:

    returns:
        (dict): collection of information about phenotype variables across
            UK Biobank cohort

    """

    ##########

    # TODO: store the contingency table reports in a text file?

    # Organize reports on contingency tables by missingness.
    organize_report_contingency_table_stratification_by_missingness(
        table=table,
        report=report,
    )

    table_missingness = organize_cohorts_phenotypes_hormones_missingness(
        table=table,
    )

    ##########

    # Prepare table to summarize phenotype variables across cohorts and models.
    # These cohorts and models are simple and do not include multiple covariates
    # for genetic analyses.
    pail_phenotypes = (
        stratification.stratify_cohorts_models_phenotypes_sets(
        table=table,
    ))
    # Collect summary records and construct table.
    records = list()
    for collection in pail_phenotypes:
        record = organize_report_cohort_model_variables_summaries_record(
            name=collection["name"],
            cohort_model=collection["cohort_model"],
            category=collection["category"],
            table=collection["table"],
        )
        records.append(record)
    # Organize table.
    table_phenotypes = pandas.DataFrame(data=records)

    ##########

    if (genotype_cohorts):
        # Read source information from file.
        table_kinship_pairs = stratification.read_source_table_kinship_pairs(
            path_dock=path_dock,
            report=False,
        )
        # Prepare table to summarize phenotype variables across cohorts and models
        # for genetic analyses.
        if (set == "sex_hormones"):
            pail_genotypes = (
                stratification.stratify_cohorts_genotypes_set_sex_hormones(
                    table=table,
                    table_kinship_pairs=table_kinship_pairs,
                    report=False,
            ))
        elif (set == "bipolar_disorder_body"):
            pail_genotypes = (
                stratification.stratify_cohorts_genotypes_set_bipolar_body(
                    table=table,
                    table_kinship_pairs=table_kinship_pairs,
                    report=report,
            ))
        else:
            print("set of cohorts and models unrecognizable...")
            pail_genotypes = list()
        # Collect summary records and construct table.
        records = list()
        for collection in pail_genotypes:
            record = organize_report_cohort_model_variables_summaries_record(
                name=collection["name"],
                cohort_model=collection["cohort_model"],
                category=collection["category"],
                table=collection["table"],
            )
            records.append(record)
        # Organize table.
        table_genotypes = pandas.DataFrame(data=records)
    else:
        table_genotypes = pandas.DataFrame()

    ##########

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_analyze_sex_cohorts_hormones()")
        utility.print_terminal_partition(level=3)
    # Collect information.
    pail = dict()
    pail["table_summary_cohorts_models_phenotypes"] = table_phenotypes
    pail["table_summary_cohorts_models_genotypes"] = table_genotypes
    pail["table_cohorts_hormones_missingness"] = table_missingness
    # Return information.
    return pail


# Plots

# TODO: TCW 2 August 2021
# TODO: include a new plot for hormones across months of assessment (should be simple)

def execute_plot_cohorts_models_phenotypes(
    table=None,
    report=None,
):
    """
    Organizes information and plots.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of collections of plot objects

    """

    # TODO: match the cohorts-models to the phenotypes and chart type...

    # Collect plots.
    #pail.update(pail_collection)
    pail = dict()
    # Histograms.
    pail["histogram"] = organize_phenotypes_plots_histogram(
        table=table,
    )
    # Box plots.
    # TODO: TCW 3 August 2021
    # TODO: box plots for hormones in different groups (assessment center?)


    # Dot plots of phenotypes ordinal variables.
    pail["dot_trajectory"] = organize_phenotypes_plots_dot_trajectory(
        table=table,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: execute_plot_cohorts_models_phenotypes()")
        utility.print_terminal_partition(level=3)

    # Return information.
    return pail





#
