
"""
...

This module collects and organizes information about estimations of SNP
heritabilities and genetic correlations.

Method for estimation of heritability from GWAS summary statistics was linkage
disequilibrium score (LDSC) regression in tool LD SCore
(https://github.com/bulik/ldsc).

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


# Standard

import sys
#print(sys.path)
import os
import math
import statistics
import pickle
import copy
import random
import itertools
import time

# Relevant

import numpy
import pandas
import scipy.stats

# Custom
import promiscuity.utility as utility

###############################################################################
# Functionality

# TODO: TCW; 14 July 2022
# TODO: For convenience in preparing the table(s) for the article, I'd like to change the format of the summaries, especially for heritability.
# TODO: Refer to notes on the heritability table for the article.

##########
# Initialization


def initialize_directories(
    restore=None,
    path_dock=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        restore (bool): whether to remove previously existing directories
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
    paths["stratification"] = os.path.join(
        path_dock, "stratification_2022-04-09"
    )
    paths["heritability"] = os.path.join(path_dock, "heritability")
    paths["genetic_correlation"] = os.path.join(
        path_dock, "genetic_correlation"
    )

    paths["collection"] = os.path.join(path_dock, "collection")
    paths["collection_stratification"] = os.path.join(
        path_dock, "collection", "stratification"
    )
    paths["collection_heritability"] = os.path.join(
        path_dock, "collection", "heritability"
    )
    paths["collection_correlation"] = os.path.join(
        path_dock, "collection", "correlation"
    )
    paths["collection_correlation_combination"] = os.path.join(
        path_dock, "collection", "correlation_combination"
    )
    paths["correlation_format_split"] = os.path.join(
        path_dock, "collection", "correlation_format_split"
    )

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["collection"])
    # Initialize directories.
    utility.create_directories(
        path=paths["collection"]
    )
    utility.create_directories(
        path=paths["collection_stratification"]
    )
    utility.create_directories(
        path=paths["collection_heritability"]
    )
    utility.create_directories(
        path=paths["collection_correlation"]
    )
    utility.create_directories(
        path=paths["collection_correlation_combination"]
    )
    utility.create_directories(
        path=paths["correlation_format_split"]
    )
    # Return information.
    return paths


##########
# Read stratification table sample sizes


def read_extract_stratification_design_study_sample_count(
    file_name=None,
    file_prefix=None,
    file_suffix=None,
    design=None,
    path_parent_directory=None,
):
    """
    Reads and extracts information from log of LDSC for heritability estimation
    from GWAS summary statistics.

    arguments:
        file_name (str): name of a file
        file_prefix (str): file name prefix to recognize relevant files and to
            extract identifier
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        design (str): name of design for set of studies
        path_parent_directory (str): path to source parent directory
            for stratification tables

    raises:

    returns:
        (dict): information about estimation of a phenotype's heritability

    """

    # Extract study's identifier.
    study = str(
        file_name.replace(str(file_prefix), "").replace(str(file_suffix), "")
    )
    # Define path to file.
    path_file = os.path.join(
        path_parent_directory, file_name
    )
    # Read table from file.
    table = pandas.read_csv(
        path_file,
        sep="\t",
        header=0,
        #dtype="string",
    )
    # Read count of samples' records in table.
    count_samples = int(table.shape[0])
    # Collect information.
    record = dict()
    record["design"] = design
    record["study"] = study
    record["count_samples"] = count_samples
    # Return information.
    return record


def read_collect_organize_stratification_design(
    design=None,
    file_prefix=None,
    file_suffix=None,
    path_parent_directory=None,
    report=None,
):
    """
    Reads, collects, and organizes the counts of samples for each stratification
    table.

    arguments:
        design (str): name of design for set of studies
        file_prefix (str): file name prefix to recognize relevant files and to
            extract identifier
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        path_parent_directory (str): path to source parent directory
            for stratification tables
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of counts of samples in stratification
           tables for individual studies

    """

    # Collect names of child files within parent directory.
    files_names = utility.extract_directory_file_names(path=path_parent_directory)
    files_relevant = list(filter(
        lambda content: (
            (str(file_prefix) in content) and
            (str(file_suffix) in content)
        ), files_names
    ))
    records = list()
    for file_name in files_relevant:
        record = read_extract_stratification_design_study_sample_count(
            file_name=file_name,
            file_prefix=file_prefix,
            file_suffix=file_suffix,
            design=design,
            path_parent_directory=path_parent_directory,
        )
        records.append(record)
        pass
    # Organize table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    columns = [
        "design",
        "study",
        "count_samples",
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    table.sort_values(
        by=["design", "study",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Return information.
    return table


def read_collect_organize_stratification_designs_studies(
    path_parent_directory=None,
    report=None,
):
    """
    Reads, collects, and organizes the counts of samples for each stratification
    table.

    Within "stratification" parent directory, first degree child directories
    correspond to "design" names and child files within those directories
    correspond to "study" names.
    Create separate tables for each "design" directory.
    Create separate table records for each "study" file.

    arguments:
        path_parent_directory (str): path to source parent directory
            for stratification tables
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of Pandas data frames

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "report: read_collect_organize_stratification_sample_counts()"
        )

    # Collect names of container directories within the parent directory.
    child_directories = utility.extract_subdirectory_names(
        path=path_parent_directory
    )
    # Iterate on design directories.
    pail = dict()
    for child_directory in child_directories:
        path_child_directory = os.path.join(
            path_parent_directory, child_directory,
        )
        pail[child_directory] = read_collect_organize_stratification_design(
            design=child_directory,
            file_prefix="table_",
            file_suffix=".tsv",
            path_parent_directory=path_child_directory,
            report=report,
        )
        # Report.
        if report:
            print(pail[child_directory])
        pass
    # Return information.
    return pail


##########
# Read heritability estimates


def read_extract_heritability_design_study_detail(
    file_name=None,
    file_prefix=None,
    file_suffix=None,
    design=None,
    study=None,
    path_parent_directory=None,
):
    """
    Reads and extracts information from log of LDSC for heritability estimation
    from GWAS summary statistics.

    arguments:
        file_name (str): name of a file
        file_prefix (str): file name prefix to recognize relevant files and to
            extract identifier
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        design (str): name of design for set of studies
        study (str): name of study
        path_parent_directory (str): path to source parent directory
            for stratification tables

    raises:

    returns:
        (dict): information about estimation of a phenotype's heritability

    """

    # Define path to file.
    path_file = os.path.join(
        path_parent_directory, file_name
    )
    # Initialize variables for extraction.
    variants = float("nan")
    heritability = float("nan")
    heritability_error = float("nan")
    confidence_95_low = float("nan")
    confidence_95_high = float("nan")
    #confidence_95 = float("nan")
    ratio = float("nan")
    ratio_error = float("nan")
    # Read relevant lines from file.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=22,
        stop=30,
    )
    # Extract information from lines.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    prefix_heritability = "Total Observed scale h2: "
    prefix_ratio = "Ratio: "
    for line in lines:
        if prefix_variants in line:
            variants = int(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_heritability in line:
            content = line.replace(prefix_heritability, "")
            contents = content.split(" (")
            heritability_test = contents[0]
            if (not "NA" in heritability_test):
                heritability = float(contents[0])
                heritability_error = float(contents[1].replace(")", ""))
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_ratio in line)
        ):
            content = line.replace(prefix_ratio, "")
            contents = content.split(" (")
            ratio_test = contents[0]
            if (not "NA" in ratio_test):
                ratio = float(contents[0])
                ratio_error = float(
                    contents[1].replace(")", "")
                )
            pass
        pass
    # Organize information.
    if (
        (not pandas.isna(heritability)) and
        (not pandas.isna(heritability_error))
    ):
        confidence_95_low = (heritability - (1.96 * heritability_error))
        confidence_95_high = (heritability + (1.96 * heritability_error))
        pass
    confidence_95 = str(
        str(round(confidence_95_low, 4)) + " ... " +
        str(round(confidence_95_high, 4))
    )
    summary_error = str(
        str(heritability) + " (" + str(round(heritability_error, 4)) + ")"
    )
    summary_interval = str(
        "(h2: " + str(heritability) + "; 95% CI: " + str(confidence_95) + ")"
    )
    # Collect information.
    record = dict()
    record["design"] = design
    record["study"] = study
    record["summary_error"] = summary_error
    record["summary_interval"] = summary_interval
    record["variants"] = variants
    record["heritability"] = heritability
    record["standard_error"] = heritability_error
    record["confidence_95_low"] = confidence_95_low
    record["confidence_95_high"] = confidence_95_high
    record["confidence_95_range"] = confidence_95
    record["ratio"] = ratio
    record["ratio_standard_error"] = ratio_error
    # Return information.
    return record


def read_collect_organize_heritability_by_many_files_per_directory(
    file_suffix=None,
    path_source_directory=None,
):
    """
    Reads, collects, and organizes heritability estimates across phenotypes.

    This function is written to collect heritability information from multiple
    files within the same directory.
    This function needs to be updated. (TCW 9 September 2021)

    arguments:
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        path_source_directory (str): path to source directory that contains
            files with heritability estimations for phenotypes

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Collect names of files for metabolites' heritabilities.
    files = utility.extract_directory_file_names(path=path_source_directory)
    files_relevant = list(filter(
        lambda content: (str(file_suffix) in content), files
    ))
    records = list()
    for file in files_relevant:
        record = read_extract_phenotype_heritability(
            file=file,
            file_suffix=file_suffix,
            path_source_directory=path_source_directory,
        )
        records.append(record)
        pass
    # Organize heritability table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["heritability"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


def read_collect_organize_heritability_design(
    design=None,
    file_name=None,
    file_prefix=None,
    file_suffix=None,
    path_parent_directory=None,
    report=None,
):
    """
    Reads, collects, and organizes heritability estimates across studies that
    correspond to a single design.

    arguments:
        design (str): name of design for set of studies
        file_name (str): name of a file
        file_prefix (str): file name prefix to recognize relevant files and to
            extract identifier
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        path_parent_directory (str): path to source parent directory
            for stratification tables
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of details about heritability estimates for
            studies within design

    """

    # Collect names of child directories within the parent directory.
    child_directories = utility.extract_subdirectory_names(
        path=path_parent_directory
    )
    # Filter child directories to those that contain relevant files.
    child_directories_relevant = list()
    for child_directory in child_directories:
        path_file = os.path.join(
            path_parent_directory, child_directory, file_name,
        )
        if os.path.isfile(path_file):
            child_directories_relevant.append(child_directory)
    # Collect and organize details for heritability estimates across studies.
    records = list()
    for child_directory in child_directories_relevant:
        path_child_directory = os.path.join(
            path_parent_directory, child_directory,
        )
        record = read_extract_heritability_design_study_detail(
            file_name=file_name,
            file_prefix=file_prefix,
            file_suffix=file_suffix,
            design=design,
            study=child_directory,
            path_parent_directory=path_child_directory,
        )
        records.append(record)
        pass
    # Organize table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    columns = [
        "design",
        "study",
        "summary_error", "summary_interval",
        "variants",
        "heritability", "standard_error", "confidence_95_range",
        "confidence_95_low", "confidence_95_high",
        "ratio",
        "ratio_standard_error",
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    table.sort_values(
        by=["design", "study",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Return information.
    return table


def read_collect_organize_heritability_designs_studies(
    path_parent_directory=None,
    report=None,
):
    """
    Reads, collects, and organizes the details on estimates of heritability.

    Within "heritability" parent directory, first degree child directories
    correspond to "design" names and second degree child directories correspond
    to "study" names.
    Create separate tables for each "design" directory.
    Create separate table records for each "study" directory.

    arguments:
        path_parent_directory (str): path to source parent directory for
            heritability estimate reports
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of Pandas data frames

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "report: read_collect_organize_heritability_designs_studies()"
        )

    # Collect names of container directories within the parent directory.
    child_directories = utility.extract_subdirectory_names(
        path=path_parent_directory
    )
    # Iterate on design directories.
    pail = dict()
    for child_directory in child_directories:
        path_child_directory = os.path.join(
            path_parent_directory, child_directory,
        )
        pail[child_directory] = read_collect_organize_heritability_design(
            design=child_directory,
            file_name="heritability_report.log",
            file_prefix="",
            file_suffix="",
            path_parent_directory=path_child_directory,
            report=report,
        )
        # Report.
        if report:
            print(pail[child_directory])
        pass
    # Return information.
    return pail


##########
# Read genetic correlation estimates

# TODO: TCW; 20 July 2022
# report genetic correlation estimates and their 95% CIs to 4 decimal places



# I think this is obsolete
def define_cohort_names_to_extract_from_table_names():
    """
    Defines names of cohorts for extraction from names of tables.

    I do not know that I want to implement this because it will decrease the
    versatility of the collection procedure.

    arguments:

    raises:

    returns:
        (list<str>): names of cohorts

    """

    # Collect names of cohorts.
    cohorts = list()
    cohorts.append("female_male_priority_male")
    cohorts.append("female")
    cohorts.append("female_menstruation_regular")
    cohorts.append("female_premenopause")
    cohorts.append("female_perimenopause")
    cohorts.append("female_postmenopause")
    cohorts.append("male")
    cohorts.append("male_age_low")
    cohorts.append("male_age_middle")
    cohorts.append("male_age_high")

    # Return information
    return cohorts


def read_extract_correlation_design_study_pair_detail(
    file_name=None,
    file_prefix=None,
    file_suffix=None,
    design=None,
    study_primary=None,
    study_secondary=None,
    path_parent_directory=None,
):
    """
    Reads and extracts information from log of LDSC for genetic correlation
    estimation from GWAS summary statistics.

    arguments:
        file_name (str): name of a file
        file_prefix (str): file name prefix to recognize relevant files and to
            extract identifier
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        design (str): name of design for set of studies
        study_primary (str): name of primary study
        study_secondary (str): name of secondary study
        path_parent_directory (str): path to source parent directory
            for stratification tables

    raises:

    returns:
        (dict): information about estimation of genetic correlation for a pair
            of studies
    """

    # Define path to file.
    path_file = os.path.join(
        path_parent_directory, file_name
    )
    # Initialize variables.
    variants = float("nan")
    correlation = float("nan")
    correlation_error = float("nan")
    interval_95 = float("nan")
    confidence_95_low = float("nan")
    confidence_95_high = float("nan")
    confidence_95_not_one = float("nan")
    confidence_95_not_zero = float("nan")
    correlation_absolute = float("nan")
    probability = float("nan")
    # Read relevant lines from file.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=25,
        stop=57,
    )
    # Initialize flag.
    missing_correlation = False
    # Extract information from lines.
    prefix_variants = ""
    suffix_variants = " SNPs with valid alleles."
    prefix_correlation = "Genetic Correlation: "
    prefix_probability = "P: "
    for line in lines:
        if suffix_variants in line:
            variants = float(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_correlation in line:
            content = line.replace(prefix_correlation, "")
            contents = content.split(" (")
            correlation_test = contents[0]
            if (not "nan" in correlation_test):
                correlation = float(contents[0])
                correlation_absolute = math.fabs(correlation)
                correlation_error = float(contents[1].replace(")", ""))
            else:
                # Main report has a missing value ("nan") for correlation.
                # Importantly, the report also has a prefix for genetic
                # correlation.
                missing_correlation = True
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_probability in line)
        ):
            probability = float(line.replace(prefix_probability, ""))
            pass
        pass

    # Determine whether main report has a missing value ("nan") correlation.
    #if missing_correlation:
    if False:
        # Initialize flag.
        summary_table = False
        # Determine whether report has a summary table.
        # Keep index of prefix title.
        # Read relevant lines from file.
        lines = utility.read_file_text_lines(
            path_file=path_file,
            start=50,
            stop=250,
        )
        index = 50
        index_title = float("nan")
        count_lines = len(lines)
        for line in lines:
            if "Summary of Genetic Correlation Results" in line:
                summary_table = True
                index_title = index
            index += 1
        if summary_table:
            # Read values from bottom table in report.
            table = pandas.read_csv(
                path_file,
                sep="\s+",
                header=(index_title + 2),
                #skip_blank_lines=True,
            )
            print(path_file)
            print(table)
    # Organize information.
    if (
        (not pandas.isna(correlation)) and
        (not pandas.isna(correlation_error))
    ):
        interval_95 = float(1.96 * correlation_error)
        confidence_95_low = (correlation - interval_95)
        confidence_95_high = (correlation + interval_95)
        if (
            (confidence_95_low < 1.0) and (confidence_95_high < 1.0)
        ):
            confidence_95_not_one = 1
        else:
            confidence_95_not_one = 0
        if (
            ((confidence_95_low > 0) and (confidence_95_high > 0)) or
            ((confidence_95_low < 0) and (confidence_95_high < 0))
        ):
            confidence_95_not_zero = 1
        else:
            confidence_95_not_zero = 0
        pass
    confidence_95 = str(
        str(round(confidence_95_low, 5)) + " ... " +
        str(round(confidence_95_high, 5))
    )
    summary = str(
        "(rg: " + str(correlation) + "; 95% CI: " + str(confidence_95) + ")"
    )
    # Collect information.
    record = dict()
    record["design"] = design
    record["study_primary"] = study_primary
    record["study_secondary"] = study_secondary
    record["summary"] = summary
    record["variants"] = variants
    record["correlation"] = correlation
    record["standard_error"] = correlation_error
    record["interval_95"] = interval_95
    record["confidence_95_low"] = confidence_95_low
    record["confidence_95_high"] = confidence_95_high
    record["confidence_95_range"] = confidence_95
    record["confidence_95_not_one"] = confidence_95_not_one
    record["confidence_95_not_zero"] = confidence_95_not_zero
    record["correlation_absolute"] = correlation_absolute
    record["probability"] = probability
    # Return information.
    return record


def read_collect_organize_correlation_design_pairs(
    design=None,
    file_name=None,
    file_prefix=None,
    file_suffix=None,
    path_parent_directory=None,
    report=None,
):
    """
    Reads, collects, and organizes genetic correlation estimates across pairs of
    studies that correspond to a single design.

    arguments:
        design (str): name of design for set of studies
        file_name (str): name of a file
        file_prefix (str): file name prefix to recognize relevant files and to
            extract identifier
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        path_parent_directory (str): path to source parent directory
            for stratification tables
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of details about heritability estimates for
            studies within design

    """

    # Collect names of primary child directories within the parent directory.
    child_primary_directories = utility.extract_subdirectory_names(
        path=path_parent_directory
    )
    # Iterate on primary child directories.
    records = list()
    for child_primary_directory in child_primary_directories:
        path_child_primary_directory = os.path.join(
            path_parent_directory, child_primary_directory,
        )
        # Collect names of secondary child directories within the primary child
        # directory.
        child_secondary_directories = utility.extract_subdirectory_names(
            path=path_child_primary_directory
        )
        # Filter secondary child directories to those that contain relevant
        # files.
        child_secondary_directories_relevant = list()
        for child_secondary_directory in child_secondary_directories:
            path_file = os.path.join(
                path_parent_directory, child_primary_directory,
                child_secondary_directory, file_name,
            )
            if os.path.isfile(path_file):
                child_secondary_directories_relevant.append(
                    child_secondary_directory
                )
        # Collect and organize details for genetic correlation estimates across
        # pairs of studies.
        for child_secondary_directory in child_secondary_directories_relevant:
            path_child_secondary_directory = os.path.join(
                path_parent_directory, child_primary_directory,
                child_secondary_directory,
            )
            record = read_extract_correlation_design_study_pair_detail(
                file_name=file_name,
                file_prefix=file_prefix,
                file_suffix=file_suffix,
                design=design,
                study_primary=child_primary_directory,
                study_secondary=child_secondary_directory,
                path_parent_directory=path_child_secondary_directory,
            )
            records.append(record)
            pass

    # Organize table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    columns = [
        "design",
        "study_primary",
        "study_secondary",
        "summary",
        "variants",
        "correlation", "standard_error", "interval_95", "confidence_95_range",
        "probability",
        "confidence_95_low", "confidence_95_high",
        "confidence_95_not_one", "confidence_95_not_zero",
        "correlation_absolute",
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    table.sort_values(
        by=["design", "study_primary", "study_secondary"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Return information.
    return table


def read_collect_organize_correlation_designs_study_pairs(
    path_parent_directory=None,
    report=None,
):
    """
    Reads, collects, and organizes the details on estimates of genetic
    correlations between pairs of primary and secondary phenotypes.

    Within "genetic_correlation" parent directory, first degree child
    directories correspond to "design" names.

    Within those "design" parent directories, first degree child directories
    correspond to "primary" study names and second degree child directories
    correspond to "secondary" study names.
    Create separate tables for each "design" directory.
    Create separate table records for each "study" directory.

    arguments:
        path_parent_directory (str): path to source parent directory for
            heritability estimate reports
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "report: read_collect_organize_correlation_designs_study_pairs()"
        )

    # Collect names of container directories within the parent directory.
    child_directories = utility.extract_subdirectory_names(
        path=path_parent_directory
    )
    # Iterate on design directories.
    pail = dict()
    for child_directory in child_directories:
        path_child_directory = os.path.join(
            path_parent_directory, child_directory,
        )
        pail[child_directory] = read_collect_organize_correlation_design_pairs(
            design=child_directory,
            file_name="correlation_report.log",
            file_prefix="",
            file_suffix="",
            path_parent_directory=path_child_directory,
            report=report,
        )
        # Report.
        if report:
            print(pail[child_directory])
        pass
    # Return information.
    return pail


##########
# Organize summary information about genetic correlations


def read_source_gwas_primary_phenotypes(
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
        (dict): information about primary studies for genetic correlations

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "gwas_format_primary_secondary_phenotypes",
        "table_gwas_primary_phenotypes.tsv"
    )
    # Read file.
    table = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype={
            "gwas_study": "string",
            "primary_phenotype": "string",
            "dependence_type": "string",
        },
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table.set_index(
        "gwas_study",
        append=False,
        drop=True,
        inplace=True
    )
    pail_primaries = table.to_dict(
        orient="index",
    )

    # Return information.
    return pail_primaries


def read_source_gwas_primary_secondary_missing_fill(
    path_dock=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_dock, "parameters", "uk_biobank",
        "gwas_format_primary_secondary_phenotypes",
        "table_gwas_primary_secondary_missing_fill.tsv"
    )
    # Read file.
    table = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype={
            "study_primary": "string",
            "study_secondary": "string",
            "correlation": "float",
            "standard_error": "float",
            "interval_95": "float",
            "probability": "float",
        },
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Return information.
    return table


def define_search_string_extraction_parameters():
    """
    Defines parameters for multiple two-tier searches and extractions on longer
    main character strings within a single column of a table.

    Example main string: 'female_premenopause_unadjust_oestradiol_bioavailable'
    Example 'cohort' tier 1 search string: 'female_premenopause'
    Example 'cohort' tier 2 search string: 'female'
    Example 'phenotype' tier 1 search string: 'oestradiol_bioavailable'
    Example 'phenotype' tier 2 search string: 'oestradiol'

    arguments:

    raises:

    returns:
        (dict<dict<list<str>>>): collection of tier 1 and tier 2 shorter, search
            character strings across multiple variables of interest

    """

    pail = dict()
    pail["cohort"] = dict()
    pail["cohort"]["search_1"] = [
        "female_menstruation_regular",
        "female_premenopause",
        "female_perimenopause",
        "female_postmenopause",
        "male_age_low",
        "male_age_middle",
        "male_age_high",
    ]
    pail["cohort"]["search_2"] = [
        "female",
    ]
    pail["cohort"]["search_3"] = [
        "male",
    ]
    pail["model"] = dict() # variable 'model' indicates context and adjustment
    pail["model"]["search_1"] = [
        "joint_1",
        "unadjust",
    ]
    pail["model"]["search_2"] = []
    pail["model"]["search_3"] = []
    pail["phenotype_secondary"] = dict()
    pail["phenotype_secondary"]["search_1"] = [
        "albumin_imputation",
        "steroid_globulin_imputation_log",
        "oestradiol_detection",
        "oestradiol_log",
        "oestradiol_imputation_log",
        "oestradiol_bioavailable_imputation_log",
        "oestradiol_free_imputation_log",
        "testosterone_detection",
        "testosterone_log",
        "testosterone_imputation_log",
        "testosterone_bioavailable_imputation_log",
        "testosterone_free_imputation_log",
    ]
    pail["phenotype_secondary"]["search_2"] = [
        "oestradiol_imputation",
        "testosterone_imputation",
    ]
    pail["phenotype_secondary"]["search_3"] = []

    # Return.
    return pail


def combine_organize_correlation_tables(
    pail_correlation_tables=None,
    path_dock=None,
    report=None,
):
    """
    Combines and organizes information from estimates of genetic correlations.

    arguments:
        pail_correlation_tables (dict<object>): collection of Pandas data-frame
            tables for estimates of genetic correlations
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Read parameter table of information to fill in missing values for missing
    # records.
    table_missing_fill = read_source_gwas_primary_secondary_missing_fill(
        path_dock=path_dock,
    )
    # Append extra records to the genetic correlation tables to fill in missing
    # values for incomplete comparisons.
    tables_combination = list(pail_correlation_tables.values())
    tables_combination.append(table_missing_fill)
    # Concatenate separate tables for estimates of genetic correlations.
    table_combination = pandas.concat(
        tables_combination,
        axis="index",
        join="outer",
        ignore_index=True,
    )
    # Define parameters for extraction of information from character string
    # names.
    pail_extractions = define_search_string_extraction_parameters()
    # Extract information from character string names.
    table_extraction = (
        utility.drive_extract_search_strings_from_table_columns_main_strings(
            pail_extractions=pail_extractions,
            temporary_column_prefix="correlation_extraction",
            column_source="study_secondary",
            table=table_combination,
            report=report,
    ))
    # Organize information.
    pail = dict()
    pail["combination"] = table_combination
    pail["extraction"] = table_extraction
    # Return information.
    return pail


def select_split_genetic_correlation_table_by_primary_phenotypes(
    column_phenotype_primary=None,
    phenotypes_primary=None,
    table=None,
    report=None,
):
    """
    Select relevant primary phenotypes and split genetic correlation records
    by these primary phenotypes.

    arguments:
        column_phenotype_primary (str): name of table's column for name of
            primary phenotype
        phenotypes_primary (list<str>): names of primary phenotypes to include
        table (object): Pandas data frame of summary information from multiple
            regressions
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables

    """

    # Copy information.
    table = table.copy(deep=True)
    # Filter records by primary phenotype.
    table = table.loc[
        (table[column_phenotype_primary].isin(phenotypes_primary)), :
    ]
    # Stratify records in the table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table.set_index(
        [column_phenotype_primary],
        append=False,
        drop=True,
        inplace=True
    )
    tables_stratification = table.groupby(
        level=column_phenotype_primary,
        axis="index",
    )
    # Collect stratification tables for each group.
    pail_tables = dict()
    for name, table_stratification in tables_stratification:
        table_stratification.reset_index(
            level=None,
            inplace=True,
            drop=False,
        )
        # Collect table.
        pail_tables[name] = table_stratification
        pass
    # Return information.
    return pail_tables


def adjust_format_split_genetic_correlation_tables(
    table=None,
    path_dock=None,
    report=None,
):
    """
    Combines and organizes information from estimates of genetic correlations.

    arguments:
        pail_correlation_tables (dict<object>): collection of Pandas data-frame
            tables for estimates of genetic correlations
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables

    """

    # Read parameter table of information about GWAS primary phenotypes.
    pail_primaries = read_source_gwas_primary_phenotypes(
        path_dock=path_dock,
    )

    # Copy information.
    table = table.copy(deep=True)

    # Filter records by primary phenotype.
    # Select records for primary phenotypes with translations.
    table = table.loc[
        (table["study_primary"].isin(list(pail_primaries.keys()))), :
    ]
    # Primary phenotype.
    table["phenotype_primary"] = table.apply(
        lambda row: pail_primaries[row["study_primary"]]["phenotype"],
        axis="columns", # apply function to each row
    )
    # Dependence type.
    table["dependence_type"] = table.apply(
        lambda row: pail_primaries[row["study_primary"]]["dependence_type"],
        axis="columns", # apply function to each row
    )
    # Model context.
    table["model_context"] = "joint"
    # Model adjustment.
    table["model_adjustment"] = table.apply(
        lambda row:
            "unadjust" if (str(row["model"]) == "unadjust") else "adjust",
        axis="columns", # apply function to each row
    )
    # Variable.
    table["variable"] = table["phenotype_secondary"]
    # Create fill records with missing values.

    # Split table records by specific primary phenotypes.
    pail = select_split_genetic_correlation_table_by_primary_phenotypes(
        column_phenotype_primary="phenotype_primary",
        phenotypes_primary=[
            "alcohol_dependence",
            "alcohol_dependence_genotype",
            "alcohol_quantity",
            "depression",
            "bipolar",
            "bipolar_type_1",
            "bipolar_type_2",
            "schizophrenia",
        ],
        table=table,
        report=report,
    )
    # Return information.
    return pail


##########
# Driver


def read_collect_organize_source(
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Each primary dictionary (such as "pail_stratification") has keys that
    correspond to names of design tables with records for individual studies.

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Stratification table sample counts.
    pail_stratification = read_collect_organize_stratification_designs_studies(
        path_parent_directory=paths["stratification"],
        report=report,
    )
    # Heritability estimates.
    pail_heritability = read_collect_organize_heritability_designs_studies(
        path_parent_directory=paths["heritability"],
        report=report,
    )
    # Genetic correlation estimates.
    pail_correlation = read_collect_organize_correlation_designs_study_pairs(
        path_parent_directory=paths["genetic_correlation"],
        report=report,
    )
    # Collect information.
    pail = dict()
    pail["stratification"] = pail_stratification
    pail["heritability"] = pail_heritability
    pail["correlation"] = pail_correlation
    return pail


##########
# Write


def write_product_study_table(
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
        path_parent, str("table_" + name + ".tsv")
    )
    # Write information to file.
    information.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=False,
        na_rep="nan",
    )
    pass


def write_product_tables(
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
        write_product_study_table(
            name=name,
            information=information[name],
            path_parent=path_parent,
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

    # Cohort tables in PLINK format.
    write_product_tables(
        information=information["stratification"],
        path_parent=paths["collection_stratification"],
    )
    write_product_tables(
        information=information["heritability"],
        path_parent=paths["collection_heritability"],
    )
    write_product_tables(
        information=information["correlation"],
        path_parent=paths["collection_correlation"],
    )
    write_product_tables(
        information=information["correlation_combination"],
        path_parent=paths["collection_correlation_combination"],
    )
    write_product_tables(
        information=information["correlation_format_split"],
        path_parent=paths["correlation_format_split"],
    )
    pass




############SCRAP###################


# TODO: TCW 9 September 2021
# TODO: MAYBE still useful
def read_collect_primary_secondaries_genetic_correlations_by_files(
    file_suffix=None,
    path_source_directory=None,
):
    """
    Reads and collects estimations of genetic correlation between phenotype and
        metabolites.

    arguments:
        file_suffix (str): file name suffix to recognize relevant files and to
            extract identifier
        path_source_directory (str): path to source parent directory for files
            with genetic correlation estimations for metabolites

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Collect names of files for metabolites' heritabilities.
    files = utility.extract_directory_file_names(path=path_source_directory)
    files_relevant = list(filter(
        lambda content: (str(file_suffix) in content), files
    ))
    records = list()
    for file in files_relevant:
        record = read_extract_phenotypes_genetic_correlation(
            file=file,
            file_suffix=file_suffix,
            path_source_directory=path_source_directory,
        )
        records.append(record)
        pass
    # Organize heritability table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["correlation_probability"],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


# TODO: TCW 9 September 2021
# TODO: PROBABLY obsolete
def read_collect_primary_secondaries_genetic_correlations_by_directories(
    file=None,
    path_parent_directory=None,
):
    """
    Reads and collects estimations of genetic correlation between phenotype and
        metabolites.

    arguments:
        file (str): name of a file
        path_parent_directory (str): path to source parent directory for files
            with heritability estimations for metabolites

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Collect names of child directories within the parent directory.
    child_directories = utility.extract_subdirectory_names(
        path=path_parent_directory
    )
    # Filter child directories to those that contain relevant files.
    child_directories_relevant = list()
    for child_directory in child_directories:
        path_file = os.path.join(
            path_parent_directory, child_directory, file,
        )
        if os.path.isfile(path_file):
            child_directories_relevant.append(child_directory)
    # Collect and organize phenotype genetic correlations.
    records = list()
    for child_directory in child_directories_relevant:
        path_child_directory = os.path.join(
            path_parent_directory, child_directory,
        )
        record = read_extract_phenotypes_genetic_correlation(
            file=file,
            file_suffix="",
            path_source_directory=path_child_directory,
        )
        # Set name of cohort and hormone.
        record["identifier"] = child_directory
        records.append(record)
        pass
    # Organize heritability table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["correlation_probability"],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


# TODO: TCW 9 September 2021
# TODO: obsolete...
def read_source_hierarchy(
    primary_study=None,
    secondary_study=None,
    paths=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        primary_study (str): identifier of primary study, consisting of a single
            GWAS
        secondary_study (str): identifier of secondary study, consisting of
            multiple GWASes
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Reference table.
    path_table_cohort_model_phenotype_reference = os.path.join(
        paths["dock"], "parameters", "sexy_alcohol",
        "table_cohort_model_phenotype_reference.tsv"
    )
    table_cohort_model_phenotype_reference = pandas.read_csv(
        path_table_cohort_model_phenotype_reference,
        sep="\t",
        header=0,
        dtype={
            "identifier": "string",
            "cohort": "string",
            "cohort_sort": "int32",
            "hormone": "string",
            "hormone_sort": "int32",
            "unadjust": "int32",
        },
    )

    # Primary phenotype heritability.
    primary_heritability = read_extract_phenotype_heritability(
        file="heritability_report.log",
        file_suffix="",
        path_source_directory=paths["heritability_studies"][primary_study],
    )

    # Secondary sample counts table.
    table_secondary_samples_counts = (
        read_collect_cohorts_sample_counts_by_files(
            file_prefix="table_",
            file_suffix=".tsv",
            path_source_directory=(
                paths["stratification_studies"][secondary_study]
            ),
    ))

    # Secondary phenotypes heritability table.
    table_secondary_heritability = (
        read_collect_phenotypes_heritabilities_by_directories(
            file="heritability_report.log",
            path_parent_directory=(
                paths["heritability_studies"][secondary_study]
            ),
    ))
    # Genetic correlations between primary and secondary phenotypes.
    table_correlations = (
        read_collect_primary_secondaries_genetic_correlations_by_directories(
            file="correlation_report.log",
            path_parent_directory=(
                paths["correlation_studies"][primary_study][secondary_study]
            ),
    ))

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_cohort_hormone_reference)
        print(primary_heritability)
        print(table_secondary_heritability)
        print(table_correlations)
        utility.print_terminal_partition(level=2)
    # Compile and return information.
    return {
        "table_cohort_model_phenotype_reference": (
            table_cohort_model_phenotype_reference
        ),
        "primary_heritability": primary_heritability,
        "table_secondary_samples_counts": table_secondary_samples_counts,
        "table_secondary_heritability": table_secondary_heritability,
        "table_correlations": table_correlations,
    }


# TODO: "read_collect_organize_heritability_design_tables"
# TODO: "read_collect_organize_heritability_study_records"
# TODO: "read_collect_organize_correlation_design_tables"
# TODO: "read_collect_organize_correlation_study_records"

##########
# Summary
# TODO: 9 September 2021
# TODO: SOME of this functionality MIGHT still be useful...
# TODO: BUT keep it VERSATILE


def organize_cohort_model_phenotype_reference_table(
    table=None,
    identifier=None,
    keep_columns=None,
):
    """
    Organizes information about general attributes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        identifier (str): name of column for identifier of cohort and hormone
        keep_columns (list<str>): names of columns to keep

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy data.
    table = table.copy(deep=True)
    # Translate column names.
    translations = dict()
    translations[identifier] = "identifier"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Select relevant columns.
    table = table.loc[
        :, table.columns.isin(keep_columns)
    ]
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True
    )
    #table["identity"].astype("float")
    #table["identity"] = pandas.to_numeric(
    #    table["identity"],
    #    errors="coerce", # force any invalid values to missing or null
    #    downcast="float",
    #)
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


def organize_table_adjust_unadjust_models(
    table=None,
):
    """
    Organizes summary table with side by side columns for adjusted and
    unadjusted regression models.

    arguments:
        table (object): Pandas data frame of genetic correlations

    raises:

    returns:
        (object): Pandas data frame of phenotypes' heritabilities and genetic
            correlations

    """

    # Copy information.
    table = table.copy(deep=True)
    # Reset index.
    table.reset_index(
        level=None,
        inplace=True
    )
    # Determine combination identifier of cohort and hormone.
    table["cohort_hormone"] = table.apply(
        lambda row:
            str(str(row["cohort"]) + "_" + str(row["hormone"])),
        axis="columns", # apply across rows
    )
    # Now split the table by adjusted or unadjusted regression models.
    table_adjust = table.copy(deep=True)
    table_unadjust = table.copy(deep=True)
    table_adjust = table_adjust.loc[
        (table_adjust["unadjust"] == 0), :
    ]
    table_unadjust = table_unadjust.loc[
        (table_unadjust["unadjust"] == 1), :
    ]
    # Set new index.
    table_adjust["cohort_hormone"].astype("string")
    table_unadjust["cohort_hormone"].astype("string")
    table_adjust.set_index(
        "cohort_hormone",
        drop=True,
        inplace=True,
    )
    table_unadjust.set_index(
        "cohort_hormone",
        drop=True,
        inplace=True,
    )
    # Select relevant columns to simplify merge.
    table_adjust.drop(
        labels=[
            "unadjust",
        ],
        axis="columns",
        inplace=True
    )
    table_unadjust = table_unadjust.loc[
        :, table_unadjust.columns.isin([
            "cohort_hormone", "identifier",
            "heritability_summary",
            "heritability", "heritability_standard_error",
            "heritability_confidence_95",
            "correlation_summary",
            "correlation", "correlation_standard_error",
            "correlation_confidence_95",
            "correlation_probability",
        ])
    ]
    # Merge tables for adjusted or unadjusted regression models.
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table_adjust.merge(
        table_unadjust,
        how="outer",
        left_on="cohort_hormone",
        right_on="cohort_hormone",
        suffixes=("_adjust", "_unadjust"),
    )
    # Return information.
    return table_merge


def combine_organize_phenotypes_hierarchy_summary_table(
    table_cohort_model_phenotype_reference=None,
    primary_heritability=None,
    table_secondary_samples_counts=None,
    table_secondary_heritability=None,
    table_correlations=None,
    threshold_secondary_heritability=None,
    threshold_false_discovery_rate=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table_cohort_model_phenotype_reference (object): Pandas data frame of
            names and information about cohorts, regression models, and
            phenotypes in study
        primary_heritability (dict): information about estimation of a
            phenotype's heritability
        table_secondary_samples_counts (object): Pandas data frame of counts of
            sample records in each cohort and model
        table_secondary_heritability (object): Pandas data frame of
            metabolites' heritability estimates
        table_correlations (object): Pandas data frame of genetic correlations
        threshold_secondary_heritability (float): threshold for heritability of
            secondary phenotypes, or missing null
        threshold_false_discovery_rate (float): threshold for false discovery
            rate
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotypes' heritabilities and genetic
            correlations

    """

    # Organize reference table.
    table_cohort_model_phenotype_reference = (
        organize_cohort_model_phenotype_reference_table(
            table=table_cohort_model_phenotype_reference,
            identifier="identifier",
            keep_columns=[
                "identifier", "cohort", "cohort_sort",
                "hormone", "hormone_sort", "unadjust"
            ],
    ))

    # Merge tables for references and heritabilities.
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge_one = table_cohort_model_phenotype_reference.merge(
        table_secondary_samples_counts,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_reference", "_counts"),
    )

    # Merge tables for references and heritabilities.
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge_two = table_merge_one.merge(
        table_secondary_heritability,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_counts", "_heritability"),
    )

    # Merge tables for heritabilities and correlations.
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge_three = table_merge_two.merge(
        table_correlations,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_heritability", "_correlation"),
    )
    # Organize cohorts, hormones, and regression models.
    table = organize_table_adjust_unadjust_models(
        table=table_merge_three,
    )

    # Introduce columns for phenotype heritability.
    table["primary_heritability"] = primary_heritability["heritability"]
    table["primary_heritability_error"] = (
        primary_heritability["heritability_standard_error"]
    )
    # Select table rows for secondary phenotypes with valid heritability
    # estimates.
    # Only filter if threshold is not missing or null.
    if (not math.isnan(threshold_secondary_heritability)):
        table = table.loc[
            (
                table["heritability_adjust"] >= threshold_secondary_heritability
            ), :
        ]

    # Calculate False Discovery Rates (FDRs).
    if False:
        table = utility.calculate_table_false_discovery_rates(
            threshold=threshold_false_discovery_rate,
            probability="correlation_probability_adjust",
            discovery="correlation_discovery_adjust",
            significance="correlation_significance_adjust",
            table=table,
        )
    # Reset index.
    table.reset_index(
        level=None,
        inplace=True
    )
    # Sort table rows.
    table.sort_values(
        by=["cohort_sort", "hormone_sort"],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )

    # Sort table columns.
    columns_sequence = [
        "cohort",
        "hormone",
        "count_samples",
        "heritability_summary_adjust",
        "heritability_adjust", "heritability_standard_error_adjust",
        "heritability_confidence_95_adjust",
        "heritability_summary_unadjust",
        "heritability_unadjust", "heritability_standard_error_unadjust",
        "heritability_confidence_95_unadjust",
        "correlation_summary_adjust",
        "correlation_adjust", "correlation_standard_error_adjust",
        "correlation_probability_adjust",
        #"correlation_discovery_adjust",
        "correlation_confidence_95_adjust",
        "correlation_summary_unadjust",
        "correlation_unadjust", "correlation_standard_error_unadjust",
        "correlation_probability_unadjust",
        "correlation_confidence_95_unadjust",
        "heritability_ratio",
        "heritability_ratio_standard_error",
        "heritability_variants",
        "correlation_absolute",
        #"correlation_significance_adjust",
        "correlation_variants",
        "cohort_hormone",
        "identifier_adjust", "identifier_unadjust",
        "primary_heritability",
        "primary_heritability_error",
    ]
    table = table[[*columns_sequence]]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("combine_organize_phenotype_metabolites_summary_table()")
        print(table)
    # Return information.
    return table


##########
# Driver


def drive_collection_report_hierarchy_studies(
    primary_study=None,
    secondary_study=None,
    path_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Collect results from primary and secondary studies.
    Primary phenotype is the same across comparisons to multiple secondary
    studies.
    Secondary studies vary by cohort, regression model, and phenotype (such as
    hormone).

    arguments:
        primary_study (str): identifier of primary study, consisting of a single
            GWAS
        secondary_study (str): identifier of secondary study, consisting of
            multiple GWASes
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: drive_collection_report_hierarchy_studies()")
        print(primary_study)
        print(secondary_study)

    # Initialize directories.
    paths = initialize_directories_hierarchy(
        primary_study=primary_study,
        secondary_study=secondary_study,
        restore=False,
        path_dock=path_dock,
    )

    # Read source information from file.
    source = read_source_hierarchy(
        primary_study=primary_study,
        secondary_study=secondary_study,
        paths=paths,
        report=False,
    )
    print(source["table_secondary_samples_counts"])

    # Organize summary table.
    table_summary = combine_organize_phenotypes_hierarchy_summary_table(
        table_cohort_model_phenotype_reference=(
            source["table_cohort_model_phenotype_reference"]
        ),
        primary_heritability=source["primary_heritability"],
        table_secondary_samples_counts=source["table_secondary_samples_counts"],
        table_secondary_heritability=source["table_secondary_heritability"],
        table_correlations=source["table_correlations"],
        threshold_secondary_heritability=float("nan"),
        threshold_false_discovery_rate=0.05,
        report=report,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=5)
        print(table_summary)

    # Collect information.
    information = dict()
    information["table_summary"] = table_summary
    # Write product information to file.
    write_product(
        primary_study=primary_study,
        secondary_study=secondary_study,
        paths=paths,
        information=information
    )

    pass

# TODO:
def drive_collection_report_pair_studies(
    path_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Collect results from pairs of primary and secondary studies.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: drive_collection_report_pair_studies()")

    # Initialize directories.
    paths = initialize_directories_pair(
        restore=False,
        path_dock=path_dock,
    )

    # Read source information from file.
    source = read_source_pair(
        paths=paths,
        report=False,
    )
    #print(source["table_secondary_samples_counts"])

    # Report.
    if report:
        utility.print_terminal_partition(level=5)
        print(table_summary)


    pass


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
    pail_source = read_collect_organize_source(
        paths=paths,
        report=True,
    )
    # Organize information from Genetic Correlation estimates for
    # further analyses.
    pail_source["correlation_combination"] = (
        combine_organize_correlation_tables(
            pail_correlation_tables=pail_source["correlation"],
            path_dock=path_dock,
            report=True,
    ))
    # Adjust format of information about Genetic Correlation estimates.
    pail_source["correlation_format_split"] = (
        adjust_format_split_genetic_correlation_tables(
            table=pail_source["correlation_combination"]["extraction"],
            path_dock=path_dock,
            report=True,
    ))
    # Write product information to file.
    write_product(
        paths=paths,
        information=pail_source,
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
