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
import uk_biobank.stratification as ukb_strat

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


# TODO: TCW 30 September 2021
# TODO: this function has problems...
def organize_report_variables_summaries_record_hormone_cohort_ordinal(
    record=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        record (dict): information for summary table record on cohort
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for cohort-specific ordinal representations of
    # hormones.
    hormones = [
        "albumin", "steroid_globulin", "vitamin_d",
        "oestradiol", "testosterone",
    ]
    cohorts = [
        "female", "female_age_low", "female_age_middle", "female_age_high",
        "male", "male_age_low", "male_age_middle", "male_age_high",
        "female_premenopause", "female_perimenopause", "female_postmenopause",
    ]
    # Iterate across cohorts.
    for cohort in cohorts:
        # Iterate across hormones.
        for hormone in hormones:
            column_ordinal = str(str(hormone) + "_" + str(cohort) + "_order")
            values = copy.deepcopy(table[column_ordinal].dropna().to_list())
            if (len(values) > 10):
                count_values_0 = len(list(filter(
                    lambda value: (value == 0),
                    values
                )))
                count_values_1 = len(list(filter(
                    lambda value: (value == 1),
                    values
                )))
                count_values_2 = len(list(filter(
                    lambda value: (value == 2),
                    values
                )))
                count_values_3 = len(list(filter(
                    lambda value: (value == 3),
                    values
                )))
            else:
                count_values_0 = 0
                count_values_1 = 0
                count_values_2 = 0
                count_values_3 = 0
            # Collect information for record.
            record[str(column_ordinal + "_0")] = str(count_values_0)
            record[str(column_ordinal + "_1")] = str(count_values_1)
            record[str(column_ordinal + "_2")] = str(count_values_2)
            record[str(column_ordinal + "_3")] = str(count_values_3)
            pass
        pass

    # Return information.
    return record



def organize_cohort_model_variables_summary_wide_record(
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
        "body", "body_log",
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
        interval_95 = float("nan")
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
                interval_95 = (1.96 * standard_error)
                confidence_95_low = (mean - interval_95)
                confidence_95_high = (mean + interval_95)
                median = numpy.nanmedian(array)
                standard_deviation = numpy.nanstd(array)
                minimum = numpy.nanmin(array)
                maximum = numpy.nanmax(array)
                pass
            pass
        # Collect information for record.
        record[str(column + "_count")] = str(count)
        record[str(column + "_mean")] = str(round(mean, 7))
        record[str(column + "_stderr")] = str(round(standard_error, 7))
        record[str(column + "_interval_95")] = str(round(interval_95, 7))
        record[str(column + "_confidence_95_low")] = str(round(
            confidence_95_low, 7
        ))
        record[str(column + "_confidence_95_high")] = str(round(
            confidence_95_high, 7
        ))
        record[str(column + "_confidence_95")] = str(
            str(round(confidence_95_low, 3)) + " ... " +
            str(round(confidence_95_high, 3))
        )
        record[str(column + "_median")] = str(round(median, 7))
        record[str(column + "_stdev")] = str(round(standard_deviation, 7))
        record[str(column + "_min")] = str(round(minimum, 7))
        record[str(column + "_max")] = str(round(maximum, 7))
        pass
    # Collect information for cohort-specific ordinal representations of
    # hormones.
    if False:
        record = (
            organize_report_variables_summaries_record_hormone_cohort_ordinal(
                record=record,
                table=table,
        ))
    # Return information.
    return record


# TODO: TCW 30 September 2021
# TODO: create dictionary...
# TODO: keys are NAMES of phenotypes, values are SORT ORDER of those phenotypes
# TODO: include these sort order values in the table


def organize_cohort_model_variables_summary_long_records(
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
        (list<dict>): information for summary table record on cohort

    """

    # Collect information for general columns.
    phenotypes = [
        "age", "body",
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

    records = list()

    # Iterate on relevant columns.
    # Collect information for record.
    for phenotype in phenotypes:
        # Collect information for cohort.
        record = dict()
        record["name"] = str(name)
        record["cohort_model"] = str(cohort_model)
        record["category"] = str(category)
        record["cohort_count"] = int(table.shape[0])
        record["phenotype"] = phenotype

        # Initialize missing values.
        count = float("nan")
        mean = float("nan")
        standard_error = float("nan")
        interval_95 = float("nan")
        confidence_95_low = float("nan")
        confidence_95_high = float("nan")
        median = float("nan")
        standard_deviation = float("nan")
        minimum = float("nan")
        maximum = float("nan")
        # Determine whether table has the column.
        if (phenotype in table.columns.to_list()):
            array = copy.deepcopy(table[phenotype].dropna().to_numpy())
            # Determine count of valid values.
            count = int(array.size)
            if (count > 10):
                # Determine mean, median, standard deviation, and standard error of
                # values in array.
                mean = numpy.nanmean(array)
                standard_error = scipy.stats.sem(array)
                interval_95 = (1.96 * standard_error)
                confidence_95_low = (mean - interval_95)
                confidence_95_high = (mean + interval_95)
                median = numpy.nanmedian(array)
                standard_deviation = numpy.nanstd(array)
                minimum = numpy.nanmin(array)
                maximum = numpy.nanmax(array)
                pass
            pass
        # Collect information for record.
        record[str("count")] = str(count)
        record[str("mean")] = str(round(mean, 7))
        record[str("standard_error")] = str(round(standard_error, 7))
        record[str("interval_95")] = str(round(interval_95, 7))
        record[str("confidence_95_low")] = str(round(confidence_95_low, 7))
        record[str("confidence_95_high")] = str(round(confidence_95_high, 7))
        record[str("range_confidence_95")] = str(
            str(round(confidence_95_low, 3)) + " ... " +
            str(round(confidence_95_high, 3))
        )
        record[str("median")] = str(round(median, 7))
        record[str("standard_deviation")] = str(round(standard_deviation, 7))
        record[str("minimum")] = str(round(minimum, 7))
        record[str("maximum")] = str(round(maximum, 7))

        records.append(record)
        pass
    # Return information.
    return records


# TODO: TCW 30 July 2021
# TODO: work on this... report separate percentages for reportability LOW and reportability HIGH


# TODO: TCW 15 February 2022
# TODO: report ALL percentage relative to the TOTAL size of the cohort!!!

def organize_cohort_hormone_missingness_record(
    name_cohort=None,
    name_variable=None,
    column_variable=None,
    column_detection=None,
    column_missingness_range=None,
    column_reportability_limit=None,
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

    "[variable]_missingness_range" 1: missing due to measurement beyond
    detection range

    "[variable]_reportability_limit" 1: not reportable due to measurement less
    than limit of detection
    "[variable]_reportability_limit" 2: not reportable due to measurement
    greater than limit of detection

    arguments:
        name_cohort (str): name of cohort
        name_variable (str): name of variable for report
        column_variable (str): name of table's column for variable's
            measurements
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
    record["variable"] = str(name_variable)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    #array_measurement_total = copy.deepcopy(table[name].to_numpy())
    #array_measurement_valid = copy.deepcopy(table[name].dropna().to_numpy())
    table_not_missing = table.dropna(
        axis="index",
        how="any",
        subset=[column_variable],
        inplace=False,
    )
    #table_missing = table[table[column_variable].isna()]
    table_missing = table.loc[
        (
            (pandas.isna(table[column_variable])
        ), :
    ]
    table_missing_range = table.loc[
        (
            #(pandas.isna(table[column_variable]) &
            (table[column_missingness_range] == 1)
        ), :
    ]
    table_reportability_low = table.loc[
        (
            #(pandas.isna(table[column_variable]) &
            (table[column_reportability_limit] == 1)
        ), :
    ]
    table_reportability_high = table.loc[
        (
            #(pandas.isna(table[column_variable]) &
            (table[column_reportability_limit] == 2)
        ), :
    ]
    table_detectable = table.loc[
        (table[column_detection] == 1), :
    ]
    table_undetectable = table.loc[
        (table[column_detection] == 0), :
    ]
    table_undetectable_low = table.loc[
        (
            (table[column_detection] == 0) &
            (table[column_reportability_limit] == 1)
        ), :
    ]
    table_undetectable_high = table.loc[
        (
            (table[column_detection] == 0) &
            (table[column_reportability_limit] == 2)
        ), :
    ]

    # Count records.
    #count_measurement_total = int(array_measurement_total.size)
    #count_measurement_valid = int(array_measurement_valid.size)
    count_total = table.shape[0]
    count_not_missing = table_not_missing.shape[0]
    count_missing = table_missing.shape[0]
    count_missing_range = table_missing_range.shape[0]
    count_reportability_low = table_reportability_low.shape[0]
    count_reportability_high = table_reportability_high.shape[0]
    count_detectable = table_detectable.shape[0]
    count_undetectable = table_undetectable.shape[0]
    count_undetectable_low = table_undetectable_low.shape[0]
    count_undetectable_high = table_undetectable_high.shape[0]

    # Calculate percentages.
    if (count_total > 0):
        percentage_not_missing = round(
            ((count_not_missing / count_total) * 100), 3
        )
        percentage_missing = round(
            ((count_missing / count_total) * 100), 3
        )
        percentage_missing_range = round(
            ((count_missing_range / count_total) * 100), 3
        )
        percentage_reportability_low = round(
            ((count_reportability_low / count_total) * 100), 3
        )
        percentage_reportability_high = round(
            ((count_reportability_high / count_total) * 100), 3
        )
        percentage_detectable = round(
            ((count_detectable / count_total) * 100), 3
        )
        percentage_undetectable = round(
            ((count_undetectable / count_total) * 100), 3
        )
        percentage_undetectable_low = round(
            ((count_undetectable_low / count_total) * 100), 3
        )
        percentage_undetectable_high = round(
            ((count_undetectable_high / count_total) * 100), 3
        )
    else:
        percentage_not_missing = float("nan")
        percentage_missing = float("nan")
        percentage_missing_range = float("nan")
        percentage_reportability_low = float("nan")
        percentage_reportability_high = float("nan")
        percentage_detectable = float("nan")
        percentage_undetectable = float("nan")
        percentage_undetectable_low = float("nan")
        percentage_undetectable_high = float("nan")
        pass

    # Collect information for record.
    record["count_cohort_samples"] = count_total
    record["not_missing"] = str(
        str(count_not_missing) +
        " (" + str(percentage_not_missing) + "%)"
    )
    record["missing"] = str(
        str(count_missing) +
        " (" + str(percentage_missing) + "%)"
    )
    record["missing_range"] = str(
        str(count_missing_range) +
        " (" + str(percentage_missing_range) + "%)"
    )
    record["unreportable_low"] = str(
        str(count_reportability_low) +
        " (" + str(percentage_reportability_low) + "%)"
    )
    record["unreportable_high"] = str(
        str(count_reportability_high) +
        " (" + str(percentage_reportability_high) + "%)"
    )
    record["detectable"] = str(
        str(count_detectable) +
        " (" + str(percentage_detectable) + "%)"
    )
    record["undetectable"] = str(
        str(count_undetectable) +
        " (" + str(percentage_undetectable) + "%)"
    )
    record["undetectable_low"] = str(
        str(count_undetectable_low) +
        " (" + str(percentage_undetectable_low) + "%)"
    )
    record["undetectable_high"] = str(
        str(count_undetectable_high) +
        " (" + str(percentage_undetectable_high) + "%)"
    )
    # Return information.
    return record


def organize_cohorts_phenotypes_hormones_missingness(
    table=None,
    report=None,
):
    """
    Organizes a summary table about missingness of hormone measurements in
    cohorts.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Copy information.
    table = table.copy(deep=True)
    # Define measurements.
    measurements = list()
    measurements.append("albumin")
    measurements.append("steroid_globulin")
    measurements.append("cholesterol")
    measurements.append("vitamin_d")
    measurements.append("testosterone")
    measurements.append("oestradiol")
    # Define cohorts for description.
    records_cohorts = ukb_strat.stratify_phenotype_cohorts_regression(
        table=table,
    )
    # Collect summary records and construct table.
    records = list()
    # Iterate on cohorts.
    for collection_cohort in records_cohorts:
        # Iterate on hormones.
        for measurement in measurements:
            # Determine column names.
            detection = str(str(measurement) + "_detection")
            missingness_range = str(str(measurement) + "_missingness_range")
            reportability_limit = str(str(measurement) + "_reportability_limit")
            # Organize information in record.
            record = organize_cohort_hormone_missingness_record(
                name_cohort=collection_cohort["name"],
                name_variable=measurement,
                column_variable=measurement,
                column_detection=detection,
                column_missingness_range=missingness_range,
                column_reportability_limit=reportability_limit,
                table=collection_cohort["table"],
            )
            # Collect records.
            records.append(record)
            pass
        pass
    # Organize table.
    table_missingness = pandas.DataFrame(data=records)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_cohorts_phenotypes_hormones_missingness()")
        utility.print_terminal_partition(level=3)
        print(table_missingness)
        pass
    # Return information.
    return table_missingness


##########
# Plot


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

    # TODO: (TCW, 3 August 2021) Consider histograms for alcohol variables...
    # "alcohol_frequency", "alcohol_previous", "alcohol_drinks_monthly",
    # "alcoholism"

    # Copy information.
    table = table.copy(deep=True)
    # Prepare table to summarize phenotype variables across cohorts and models.
    # These cohorts and models are simple and do not include multiple covariates
    # for genetic analyses.
    # Collect records of information about each cohort and model.
    records = list()
    records_novel = (
        ukb_strat.stratify_set_primary_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)
    records_novel = (
        ukb_strat.stratify_set_alcohol_sex_menopause_age(
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
            "body",
            "month",
            "albumin", "albumin_imputation",
            "steroid_globulin", "steroid_globulin_imputation",
            "oestradiol", "oestradiol_imputation",
            "oestradiol_bioavailable", "oestradiol_bioavailable_imputation",
            "oestradiol_free", "oestradiol_free_imputation",
            "testosterone", "testosterone_imputation",
            "testosterone_bioavailable", "testosterone_bioavailable_imputation",
            "testosterone_free", "testosterone_free_imputation",
            "vitamin_d", "vitamin_d_imputation",
            "alcohol_frequency",
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


def organize_phenotypes_plots_dot_trajectory_assessment(
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
        ukb_strat.stratify_set_primary_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)
    records_novel = (
        ukb_strat.stratify_set_alcohol_sex_menopause_age(
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
            "alcohol_frequency",
        ]
        # Iterate on phenotypes.
        for phenotype in phenotypes:
            # Plot hormone trajectories across month of assessment.
            # Chart.
            name_label = str(str(record["name"]) + "_" + str(phenotype))
            name_plot = str(name_label + "_dot_assessment")
            pail[name_plot] = plot_variable_means_dot_trajectory(
                label_title=name_label,
                column_phenotype=phenotype,
                column_trajectory="month",
                threshold_trajectory=13,
                title_abscissa="month of assessment and blood draw",
                title_ordinate="mean concentration (95% C.I.)",
                table=record["table"],
            )
            pass
        pass
    # Return information.
    return pail


def organize_phenotypes_plots_dot_trajectory_menstruation(
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
        ukb_strat.stratify_set_primary_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)
    records_novel = (
        ukb_strat.stratify_set_female_hormone_alteration(
            table=table,
        )
    )
    records.extend(records_novel)
    records_novel = (
        ukb_strat.stratify_set_alcohol_sex_menopause_age(
            table=table,
        )
    )
    records.extend(records_novel)

    # Collect information for plots.
    pail = dict()
    # Iterate on tables for cohorts and models.
    for record in records:
        # Determine whether the current cohort and model table is relevant to
        # the plot of phenotypes across days of menstrual cycle.
        if (record["menstruation"]):
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
                "alcohol_frequency",
            ]
            # Iterate on phenotypes.
            for phenotype in phenotypes:
                # Plot hormone trajectories across days of the menstrual cycle.
                # Chart.
                name_label = str(str(record["name"]) + "_" + str(phenotype))
                name_plot = str(name_label + "_dot_menstruation")
                pail[name_plot] = plot_variable_means_dot_trajectory(
                    label_title=name_label,
                    column_phenotype=phenotype,
                    column_trajectory="menstruation_days_threshold",
                    threshold_trajectory=30,
                    title_abscissa="days of menstrual cycle",
                    title_ordinate="mean concentration (95% C.I.)",
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


# TODO: I think

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
    #organize_report_contingency_table_stratification_by_missingness(
    #    table=table,
    #    report=report,
    #)

    # Organize report summary for missingness in biochemical measurements.
    table_missingness = organize_cohorts_phenotypes_hormones_missingness(
        table=table,
        report=report,
    )

    if False:

        ##########

        # Prepare table to summarize phenotype variables across cohorts and models.
        # These cohorts and models are simple and do not include multiple covariates
        # for genetic analyses.
        pail_phenotypes = (
            ukb_strat.stratify_cohorts_models_phenotypes_sets(
            table=table,
        ))
        # Collect summary records and construct table.
        records = list()
        for collection in pail_phenotypes:
            records_cohort = organize_cohort_model_variables_summary_long_records(
                name=collection["name"],
                cohort_model=collection["cohort_model"],
                category=collection["category"],
                table=collection["table"],
            )
            records.extend(records_cohort)
        # Organize table.
        table_phenotypes = pandas.DataFrame(data=records)

        ##########

        if (genotype_cohorts):
            # Read source information from file.
            table_kinship_pairs = ukb_strat.read_source_table_kinship_pairs(
                path_dock=path_dock,
                report=report,
            )
            # Prepare table to summarize phenotype variables across cohorts and models
            # for genetic analyses.
            if (set == "sex_hormones"):
                pail_genotypes = (
                    ukb_strat.stratify_genotype_cohorts_linear_set_sex_hormones(
                        table=table,
                        table_kinship_pairs=table_kinship_pairs,
                        report=report,
                ))
            elif (set == "bipolar_disorder_body"):
                pail_genotypes = (
                    ukb_strat.stratify_cohorts_genotypes_set_bipolar_body(
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
                record = organize_cohort_model_variables_summary_wide_record(
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
    pail["table_cohorts_measurements_missingness"] = table_missingness
    #pail["table_summary_cohorts_models_phenotypes"] = table_phenotypes
    #pail["table_summary_cohorts_models_genotypes"] = table_genotypes
    # Return information.
    return pail


# Plots


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
    pail["dot_trajectory_assessment"] = (
        organize_phenotypes_plots_dot_trajectory_assessment(
            table=table,
    ))
    pail["dot_trajectory_menstruation"] = (
        organize_phenotypes_plots_dot_trajectory_menstruation(
            table=table,
    ))

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: execute_plot_cohorts_models_phenotypes()")
        utility.print_terminal_partition(level=3)

    # Return information.
    return pail





#
