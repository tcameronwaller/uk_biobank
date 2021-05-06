
"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

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

# Relevant

import numpy
import scipy.stats
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"

# Custom
import promiscuity.utility as utility
import promiscuity.plot as plot

###############################################################################
# Functionality


##########
# Initialization


def initialize_directories_cohorts(
    path_parent=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        path_parent (str): path to parent directory

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect information.
    paths = dict()
    sexes = ["female", "male",]
    alcoholisms = [
        "alcoholism_1", "alcoholism_2", "alcoholism_3", "alcoholism_4",
    ]
    groups = ["all", "case", "control"]
    for sex in sexes:
        paths[sex] = dict()
        for alcoholism in alcoholisms:
            paths[sex][alcoholism] = dict()
            for group in groups:
                paths[sex][alcoholism][group] = os.path.join(
                    path_parent, "cohorts", sex, alcoholism, group
                )
                # Initialize directories.
                utility.create_directories(
                    path=paths[sex][alcoholism][group]
                )
                pass
            pass
        pass
    # Return information.
    return paths


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
    paths["quality"] = os.path.join(
        path_dock, "organization", "quality"
    )
    paths["cohorts"] = os.path.join(
        path_dock, "organization", "cohorts"
    )

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["organization"])
    # Initialize directories.
    utility.create_directories(
        path=paths["organization"]
    )
    utility.create_directories(
        path=paths["quality"]
    )
    utility.create_directories(
        path=paths["cohorts"]
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
    path_table_assembly = os.path.join(
        path_dock, "assembly", "table_phenotypes.pickle"
    )
    path_table_ukb_samples = os.path.join(
        path_dock, "access", "ukb46237_imp_chr21_v3_s487320.sample"
    )
    # Read information from file.
    table_assembly = pandas.read_pickle(
        path_table_assembly
    )
    table_ukb_samples = pandas.read_csv(
        path_table_ukb_samples,
        sep="\s+",
        header=0,
        dtype="string",
    )
    # Compile and return information.
    return {
        "table_assembly": table_assembly,
        "table_ukb_samples": table_ukb_samples,
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


def convert_table_columns_variables_types_float(
    columns=None,
    table=None,
):
    """
    Converts data variable types.

    arguments:
        columns (list<str>): names of columns
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    table = table.copy(deep=True)
    # Convert data variable types.
    for column in columns:
        table[column] = pandas.to_numeric(
            table[column],
            errors="coerce", # force any invalid values to missing or null
            downcast="float",
        )
    # Return information.
    return table


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
    table_type = convert_table_columns_variables_types_float(
        columns=columns_match,
        table=table,
    )
    # Return information.
    return table_type


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
        (dict): collection of information about phenotype variables

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
    # Return information.
    return table


##########
# Sex, age, and body


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


def organize_sex_age_body_variables(
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

    # Convert variable types.
    columns_type = [
        "31-0.0", "22001-0.0", "21022-0.0", "21001-0.0"
    ]
    table = convert_table_columns_variables_types_float(
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
    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=["body_mass_index"],
        table=table,
    )
    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "31-0.0",
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
        "age", "body_mass_index", "body_mass_index_log",
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
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


##########
# Sex hormones
# Review:
# 11 March 2021: TCW verified formulas for estimation of free oestradiol
# 11 March 2021: TCW verified formulas for estimation of free testosterone
# 11 March 2021: TCW verified UK Biobank fields and their codings.

# TODO: implement bioavailable oestradiol and testosterone
# TODO: also implement the natural log definition of bioavailable testosterone...
# Chung, Pathology Informatics, 2017 (PubMed:28828199)

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
    # Return information.
    return table

# TODO: check accuracy of formula implementation
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

# TODO: check accuracy of formula implementation
def calculate_estimation_bioavailable_testosterone(
    testosterone_free=None,
    albumin=None,
    factors_concentration=None,
    associations=None,
):
    """
    Calculates an estimation of bioavailable testosterone (not bound to steroid
    hormone binding globulin).

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

# TODO: check accuracy of formula implementation
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

# TODO: check accuracy of formula implementation
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

    # Copy data.
    table = table.copy(deep=True)
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
    table = convert_hormone_concentration_units_moles_per_liter(
        table=table,
        factors_concentration=factors_concentration,
    )
    # Define association constants for calculation of free hormones.
    # units: L/mol
    associations = dict()
    associations["shbg_test"] = 5.97E8 #
    associations["shbg_oest"] = 3.14E8 #
    associations["albumin_test"] = 4.06E4 #
    associations["albumin_oest"] = 4.21E4 #
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
    # Convert variable types.
    columns_hormones = [
        "albumin", "steroid_globulin",
        "oestradiol", "oestradiol_free", "oestradiol_bioavailable",
        "testosterone", "testosterone_free", "testosterone_bioavailable",
    ]
    table = convert_table_columns_variables_types_float(
        columns=columns_hormones,
        table=table,
    )
    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=columns_hormones,
        table=table,
    )
    # Remove columns for variables that are not necessary anymore.
    # Pandas drop throws error if column names do not exist.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=[
            "30600-0.0", "30830-0.0", "30850-0.0", "30800-0.0",
        ],
        axis="columns",
        inplace=True
    )
    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        #"eid",
        "IID",
        "albumin", "albumin_log", "steroid_globulin", "steroid_globulin_log",
        "oestradiol", "oestradiol_log",
        "oestradiol_free", "oestradiol_free_log",
        "oestradiol_bioavailable", "oestradiol_bioavailable_log",
        "testosterone", "testosterone_log",
        "testosterone_free", "testosterone_free_log",
        "testosterone_bioavailable", "testosterone_bioavailable_log",
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    # Organize information for bioavailable testosterone.
    # Select records with valid (non-null) values of relevant variables.
    # Exclude missing values first to avoid interpretation of "None" as False.
    table_bioavailability = select_valid_records_all_specific_variables(
        names=[
            "eid", "IID",
            "sex", "sex_text",# "age", "body_mass_index",
            "testosterone", "steroid_globulin", "albumin",
        ],
        prefixes=[],#["genotype_pc_",],
        table=table_clean,
        drop_columns=True,
        report=False,
    )
    table_bioavailability_female = table_bioavailability.loc[
        (table_bioavailability["sex_text"] == "female"), :
    ]
    table_bioavailability_male = table_bioavailability.loc[
        (table_bioavailability["sex_text"] == "male"), :
    ]
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
        # Valid information for bioavailable testosterone.
        print("Valid information for bioavailable testosterone: ")
        print("... requires testosterone, steroid globulin, and albumin...")
        print(
            "Count females: " + str(table_bioavailability_female.shape[0])
        )
        print(
            "Count males: " + str(table_bioavailability_male.shape[0])
        )
        utility.print_terminal_partition(level=3)
        organize_report_column_pair_correlations(
            column_one="testosterone",
            column_two="testosterone_free",
            table=table_report,
        )
        organize_report_column_pair_correlations(
            column_one="oestradiol",
            column_two="oestradiol_free",
            table=table_report,
        )
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    # Return information.
    return pail


def organize_plot_cohort_sex_hormone_variable_distributions(
    prefix=None,
    bins=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        prefix (str): prefix for cohort and figure name
        bins (int): count of bins for histogram
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): collection of information about figures

    """

    # Collect information for plots.
    pail = dict()
    # Iterate on relevant variables.
    columns = [
        "oestradiol", "oestradiol_log",
        "oestradiol_free", "oestradiol_free_log",
        "testosterone", "testosterone_log",
        "testosterone_free", "testosterone_free_log",
        "steroid_globulin", "steroid_globulin_log",
        "albumin", "albumin_log",
        "age",
    ]
    for column in columns:
        if len(str(prefix)) > 0:
            name = str(prefix + "_" + column)
        else:
            name = column
        pail[name] = plot_variable_values_histogram(
            name=name,
            array=table[column].dropna().to_numpy(),
            bins=bins,
        )
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

    Data-Field "3700": "time since last menstrual period"
    UK Biobank data coding "100291" for variable field "3700".
    "days": 0 - 365
    "do not know": -1
    "prefer not to answer": -3

    Note: "Please count from the first day of your last menstrual period."

    Accommodate inexact float values.

    arguments:
        field_3700 (float): UK Biobank field 3700, count of days since previous
            menstruation (menstrual period)

    raises:

    returns:
        (bool): interpretation value

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


def interpret_menopause_hysterectomy(
    field_2724=None,
):
    """
    Intepret UK Biobank's coding for field 2724.

    Only interpret whether the variable indicates that the person had a
    hysterectomy.
    Do not propagate missing values since this variable primarily concerns a
    different subject.

    Data-Field "2724": "Had menopause"
    UK Biobank data coding "100579" for variable field "2724".
    "no": 0
    "yes": 1
    "not sure - had a hysterectomy": 2
    "not sure - other reason": 3
    "prefer not to answer": -3

    Accommodate inexact float values.

    arguments:
        field_2724 (float): UK Biobank field 2724, whether person has
            experienced menopause

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
            value = 0
        elif (0.5 <= field_2724 and field_2724 < 1.5):
            # 1: "yes"
            value = 0
        elif (1.5 <= field_2724 and field_2724 < 2.5):
            # 2: "not sure - had a hysterectomy"
            value = 1
        elif (2.5 <= field_2724 and field_2724 < 3.5):
            # 3: "not sure - other reason"
            value = 0
        elif (-3.5 <= field_2724 and field_2724 < -2.5):
            # -3: "prefer not to answer"
            value = 0
        else:
            # uninterpretable
            value = 0
    else:
        # null
        value = 0
    # Return.
    return value


def interpret_menopause_natural(
    field_2724=None,
):
    """
    Intepret UK Biobank's coding for field 2724.

    Only interpret whether the variable indicates that the person had a
    natural menopause, not involving hysterectomy.

    Data-Field "2724": "Had menopause"
    UK Biobank data coding "100579" for variable field "2724".
    "no": 0
    "yes": 1
    "not sure - had a hysterectomy": 2
    "not sure - other reason": 3
    "prefer not to answer": -3

    Accommodate inexact float values.

    Interpret missing or null values as False.
    The main justification is to avoid loss of persons from cohorts.
    This interpretation simplifies the definitions of "menopause_binary" and
    "menopause_ordinal".
    Only a "True" value of "menopause_natural" asserts postmenopause definition.
    A "False" (or null) value of "menopause_natural" relies on other variables
    for definition, including "oophorectomy", "menstruation_days", and "age".

    arguments:
        field_2724 (float): UK Biobank field 2724, whether person has
            experienced menopause

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
            value = 0
        elif (0.5 <= field_2724 and field_2724 < 1.5):
            # 1: "yes"
            value = 1
        elif (1.5 <= field_2724 and field_2724 < 2.5):
            # 2: "not sure - had a hysterectomy"
            value = float("nan")
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
        (bool): interpretation value

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
        (bool): interpretation value

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
        (bool): interpretation value

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
    field_3700=None,
):
    """
    Determine count of days since previous menstruation (menstrual period).

    arguments:
        sex_text (str): textual representation of sex selection
        field_3700 (float): UK Biobank field 3700, days since previous
            menstruation (menstrual period)

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret days since previous menstruation.
    menstruation_days = interpret_menstruation_days(
        field_3700=field_3700,
    )
    # Comparison.
    # Only define days since previous menstruation for female persons.
    if (
        (sex_text == "female")
    ):
        # Female.
        value = menstruation_days
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
    menopause_natural = interpret_menopause_natural(
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
    menopause_natural = interpret_menopause_natural(
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

    This definition is permissive of missing values in some variables.

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
    menopause_natural = interpret_menopause_natural(
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
                (menstruation_days >= threshold_menstruation_days_post)
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
        # Determine perimenopause
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
                    (menstruation_days >= threshold_menstruation_days_pre) and
                    (menstruation_days < threshold_menstruation_days_post)
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
                    (menstruation_days < threshold_menstruation_days_pre)
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
    menopause_natural = interpret_menopause_natural(
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


def determine_female_menstruation_phase(
    sex_text=None,
    menopause_ordinal=None,
    menstruation_days=None,
    threshold_premenopause=None,
    threshold_perimenopause=None,
):
    """
    Determine a female person's categorical menstruation phase.

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


def determine_female_pregnancy(
    sex_text=None,
    menopause_binary=None,
    field_3140=None,
):
    """
    Determine whether female persons were pregnant.

    This definition of pregnancy uses a definition of menopause
    (pre, peri, post) to avoid missing information and increase confidence.

    arguments:
        sex_text (str): textual representation of sex selection
        menopause_binary (float): binary logical representation of whether
            person has experienced menopause
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
    if (
        (sex_text == "female")
    ):
        if (pandas.isna(pregnancy)):
            # The variable has an uncertain value.
            # Determine whether person has experienced menopause.
            if (
                (not pandas.isna(menopause_binary)) and
                (menopause_binary == 1)
            ):
                value = 0
            else:
                value = float("nan")
        else:
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


def organize_report_cohort_variables_summaries_record(
    cohort=None,
    category=None,
    table=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        cohort (str): name of cohort
        category (str): name of category for cohort, to facilitate table sorts
            and comparisons
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for cohort.
    record = dict()
    record["combination"] = str(category + "_" + cohort)
    record["category"] = str(category)
    record["cohort"] = str(cohort)
    record["cohort_count"] = int(table.shape[0])
    # Collect information for general columns.
    columns = [
        "age",
        "menstruation_days",
        "albumin", "albumin_log", "steroid_globulin", "steroid_globulin_log",
        "oestradiol", "oestradiol_log",
        "oestradiol_free", "oestradiol_free_log",
        "testosterone", "testosterone_log",
        "testosterone_free", "testosterone_free_log",
    ]
    # Iterate on relevant columns.
    # Collect information for record.
    for column in columns:
        array = copy.deepcopy(table[column].dropna().to_numpy())
        # Determine count of valid values.
        count = int(array.size)
        if (count > 5):
            # Determine mean, median, standard deviation, and standard error of
            # values in array.
            mean = numpy.nanmean(array)
            median = numpy.nanmedian(array)
            standard_deviation = numpy.nanstd(array)
            standard_error = scipy.stats.sem(array)
        else:
            mean = float("nan")
            median = float("nan")
            standard_deviation = float("nan")
            standard_error = float("nan")
        # Collect information for record.
        record[str(column + "_count")] = str(count)
        record[str(column + "_mean")] = str(round(mean, 2))
        record[str(column + "_median")] = str(round(median, 2))
        record[str(column + "_stdev")] = str(round(standard_deviation, 2))
        record[str(column + "_stderr")] = str(round(standard_error, 2))
        pass
    # Return information.
    return record


def organize_report_female_male_cohorts_variables(
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
    # Stratify cohorts and collect tables.
    cohorts = list()

    # Sex

    cohort = dict()
    cohort["category"] = "sex"
    cohort["name"] = "female"
    cohort["table"] = table.loc[
        (table["sex_text"] == "female"), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "sex"
    cohort["name"] = "male"
    cohort["table"] = table.loc[
        (table["sex_text"] == "male"), :
    ]
    cohorts.append(cohort)

    # Pregnancy

    cohort = dict()
    cohort["category"] = "pregnancy"
    cohort["name"] = "female_pregnancy_definite"
    cohort["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["3140-0.0"] == 1)), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "pregnancy"
    cohort["name"] = "female_pregnancy_unsure"
    cohort["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["3140-0.0"] == 2)), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "pregnancy"
    cohort["name"] = "female_pregnancy_yes"
    cohort["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["pregnancy"] == 1)), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "pregnancy"
    cohort["name"] = "female_pregnancy_no"
    cohort["table"] = table.loc[
        ((table["sex_text"] == "female") & (table["pregnancy"] == 0)), :
    ]
    cohorts.append(cohort)

    # Menopause

    cohort = dict()
    cohort["category"] = "menopause"
    cohort["name"] = "female_menopause_unsure"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["2724-0.0"] == 3)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause"
    cohort["name"] = "female_menopause_blank"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["2724-0.0"] == -3)
        ), :
    ]
    cohorts.append(cohort)

    # Menopause binary

    cohort = dict()
    cohort["category"] = "menopause_binary"
    cohort["name"] = "female_premenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause_binary"
    cohort["name"] = "female_postmenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    # Menopause binary strict

    cohort = dict()
    cohort["category"] = "menopause_binary_strict"
    cohort["name"] = "female_premenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary_strict"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause_binary_strict"
    cohort["name"] = "female_postmenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_binary_strict"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    # Menopause ordinal

    cohort = dict()
    cohort["category"] = "menopause_ordinal"
    cohort["name"] = "female_premenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause_ordinal"
    cohort["name"] = "female_perimenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause_ordinal"
    cohort["name"] = "female_postmenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2)
        ), :
    ]
    cohorts.append(cohort)

    # Menopause ordinal strict

    cohort = dict()
    cohort["category"] = "menopause_ordinal_strict"
    cohort["name"] = "female_premenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal_strict"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause_ordinal_strict"
    cohort["name"] = "female_perimenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal_strict"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menopause_ordinal_strict"
    cohort["name"] = "female_postmenopause"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal_strict"] == 2)
        ), :
    ]
    cohorts.append(cohort)

    # Menstruation phase

    cohort = dict()
    cohort["category"] = "menstruation_phase"
    cohort["name"] = "female_menstruation_follicular"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_phase"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menstruation_phase"
    cohort["name"] = "female_menstruation_luteal"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menstruation_phase"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menstruation_phase"
    cohort["name"] = "female_premenopause_menstruation_follicular"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["menstruation_phase"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menstruation_phase"
    cohort["name"] = "female_premenopause_menstruation_luteal"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["menstruation_phase"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menstruation_phase"
    cohort["name"] = "female_perimenopause_menstruation_follicular"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["menstruation_phase"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "menstruation_phase"
    cohort["name"] = "female_perimenopause_menstruation_luteal"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["menstruation_phase"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    # Hormone-alteration therapies

    cohort = dict()
    cohort["category"] = "oral_contraception"
    cohort["name"] = "female_oral_contraception_yes"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["oral_contraception"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "oral_contraception"
    cohort["name"] = "female_oral_contraception_no"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["oral_contraception"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_replacement"
    cohort["name"] = "female_hormone_replacement_yes"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_replacement"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_replacement"
    cohort["name"] = "female_hormone_replacement_no"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_replacement"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_hormone_alteration_yes"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_hormone_alteration_no"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_premenopause_hormone_alteration_yes"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_premenopause_hormone_alteration_no"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 0) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_perimenopause_hormone_alteration_yes"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_perimenopause_hormone_alteration_no"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 1) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_postmenopause_hormone_alteration_yes"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["hormone_alteration"] == 1)
        ), :
    ]
    cohorts.append(cohort)

    cohort = dict()
    cohort["category"] = "hormone_alteration"
    cohort["name"] = "female_postmenopause_hormone_alteration_no"
    cohort["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["pregnancy"] == 0) &
            (table["menopause_ordinal"] == 2) &
            (table["hormone_alteration"] == 0)
        ), :
    ]
    cohorts.append(cohort)

    # Collect records for cohorts.
    records = list()
    for cohort in cohorts:
        record = organize_report_cohort_variables_summaries_record(
            cohort=cohort["name"],
            category=cohort["category"],
            table=cohort["table"],
        )
        records.append(record)

    # Organize table.
    table_summary = pandas.DataFrame(data=records)
    # Return information
    return table_summary


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
        "3700-0.0",
        "2784-0.0", "2794-0.0", "2804-0.0",
        "2814-0.0", "3536-0.0", "3546-0.0",
    ]
    table = convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )

    # Determine count of days since person's last menstrual period.
    table["menstruation_days"] = table.apply(
        lambda row:
            determine_female_menstruation_days(
                sex_text=row["sex_text"],
                field_3700=row["3700-0.0"],
            ),
        axis="columns", # apply across rows
    )
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
                menstruation_days=row["menstruation_days"],
                threshold_menstruation_days_pre=60, # threshold in days
                threshold_menstruation_days_post=365, # threshold in days
                field_2724=row["2724-0.0"],
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
    # Determine female persons' categorical menstruation phase.
    table["menstruation_phase"] = table.apply(
        lambda row:
            determine_female_menstruation_phase(
                sex_text=row["sex_text"],
                menopause_ordinal=row["menopause_ordinal"],
                menstruation_days=row["menstruation_days"],
                threshold_premenopause=13, # threshold in days
                threshold_perimenopause=15, # threshold in days
            ),
        axis="columns", # apply across rows
    )
    # Determine whether female persons were pregnant.
    table["pregnancy"] = table.apply(
        lambda row:
            determine_female_pregnancy(
                sex_text=row["sex_text"],
                menopause_binary=row["menopause_binary"],
                field_3140=row["3140-0.0"],
            ),
        axis="columns", # apply across rows
    )
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
            "3700-0.0",
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
        "menstruation_days",
        "hysterectomy", "oophorectomy", "hysterectomy_or_oophorectomy",
        "menopause_binary", "menopause_binary_strict",
        "menopause_ordinal", "menopause_ordinal_strict",
        "menstruation_phase",
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
    table_report_summary = organize_report_female_male_cohorts_variables(
        table=table,
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
        print(table_report_summary)
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["table_clean"] = table_clean
    pail["table_report"] = table_report
    pail["table_report_summary"] = table_report_summary
    # Return information.
    return pail


##########
# Neuroticism


def organize_neuroticism_variables(
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
        "neuroticism"
    ]
    table = convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Transform variables' values to normalize distributions.
    table = utility.transform_normalize_table_continuous_ratio_variables(
        columns=["neuroticism"],
        table=table,
    )
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


def determine_alcohol_consumption_frequency(
    frequency=None,
):
    """
    Translate information from UK Biobank about whether person
    never consumes any alcohol.

    Accommodate inexact float values.

    arguments:
        frequency (float): frequency of alcohol consumption, UK Biobank field
            1558

    raises:

    returns:
        (float): ordinal representation of person's frequency of alcohol
            consumption

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (
        (not pandas.isna(frequency)) and
        (0.5 <= frequency and frequency < 6.5)
    ):
        # The variable has a valid value.
        if (5.5 <= frequency and frequency < 6.5):
            # "never"
            alcohol_frequency = 0
        elif (4.5 <= frequency and frequency < 5.5):
            # "special occasions only"
            alcohol_frequency = 1
        elif (3.5 <= frequency and frequency < 4.5):
            # "one to three times a month"
            alcohol_frequency = 2
        elif (2.5 <= frequency and frequency < 3.5):
            # "once or twice a week"
            alcohol_frequency = 3
        elif (1.5 <= frequency and frequency < 2.5):
            # "three or four times a week"
            alcohol_frequency = 4
        elif (0.5 <= frequency and frequency < 1.5):
            # "daily or almost daily"
            alcohol_frequency = 5
    else:
        # "prefer not to answer" or null
        alcohol_frequency = float("nan")
    # Return information.
    return alcohol_frequency


def determine_current_alcohol_consumption(
    consumer=None,
    alcohol_frequency=None,
):
    """
    Translate information from UK Biobank about whether person consumes alcohol
    currently.

    Accommodate inexact float values.

    arguments:
        consumer (float): status of alcohol consumption, UK Biobank field
            20117
        alcohol_frequency (float): ordinal representation of person's frequency
            of alcohol consumption, derivation of UK Biobank field 1558

    raises:

    returns:
        (float): binary representation of whether person consumes alcohol
            currently

    """

    # Interpret consumer status.
    if (
        (not pandas.isna(consumer)) and
        (-0.5 <= consumer and consumer < 2.5)
    ):
        # The variable has a valid value.
        if (1.5 <= consumer and consumer < 2.5):
            # "current"
            consumer_boolean = True
        else:
            # "never", or "previous"
            consumer_boolean = False
    else:
        # "prefer not to answer" or null
        consumer_boolean = float("nan")
    # Interpret consumption frequency.
    if (not math.isnan(alcohol_frequency)):
        # The variable has a valid value.
        if (0.5 <= alcohol_frequency and alcohol_frequency < 5.5):
            # The person consumes alcohol currently.
            frequency_boolean = True
        else:
            frequency_boolean = False
    else:
        frequency_boolean = float("nan")
    # Integrate information from multiple variables.
    if (
        (not math.isnan(consumer_boolean)) or
        (not math.isnan(frequency_boolean))
    ):
        if (consumer_boolean or frequency_boolean):
            alcohol_current = 1
        else:
            alcohol_current = 0
    else:
        alcohol_current = float("nan")
    # Return information.
    return alcohol_current


def determine_previous_alcohol_consumption(
    consumer=None,
    former=None,
    comparison=None,
    alcohol_current=None,
):
    """
    Translate information from UK Biobank about whether person
    never consumes any alcohol.

    Accommodate inexact float values.

    arguments:
        consumer (float): status of alcohol consumption, UK Biobank field
            20117
        former (float): former alcohol consumption, UK Biobank field
            3731
        comparison (float): comparison to previous alcohol consumption, UK
            Biobank field 1628
        alcohol_current (float): binary representation of whether person
            consumes alcohol currently

    raises:

    returns:
        (float): binary representation of whether person consumed alcohol
            previously

    """

    # Interpret consumer.
    if (
        (not pandas.isna(consumer)) and
        (-0.5 <= consumer and consumer < 2.5)
    ):
        # The variable has a valid value.
        if (0.5 <= consumer and consumer < 2.5):
            # "current" or "previous"
            consumer_boolean = True
        else:
            # "never"
            consumer_boolean = False
    else:
        # "prefer not to answer" or null
        consumer_boolean = float("nan")
    # Interpret former.
    if (
        (not pandas.isna(former)) and
        (-0.5 <= former and former < 1.5)
    ):
        # The variable has a valid value.
        if (0.5 <= former and former < 1.5):
            # "yes"
            former_boolean = True
        else:
            # "no"
            former_boolean = False
    else:
        # "prefer not to answer" or null
        former_boolean = float("nan")
    # Interpret comparison.
    if (
        (not pandas.isna(comparison)) and
        (-0.5 <= former and former < 3.5)
    ):
        # "less nowadays", "about the same", or "more nowadays"
        comparison_boolean = True
    else:
        # "do not know", "prefer not to answer", or null
        comparison_boolean = float("nan")
    # Interpret alcohol current.
    if (not math.isnan(alcohol_current)):
        # The variable has a valid value.
        if (0.5 <= alcohol_current and alcohol_current < 1.5):
            # The person consumes alcohol currently.
            alcohol_current_boolean = True
        else:
            alcohol_current_boolean = False
    else:
        alcohol_current_boolean = float("nan")
    # Integrate information from multiple variables.
    if (
        (not math.isnan(consumer_boolean)) or
        (not math.isnan(former_boolean)) or
        (not math.isnan(comparison_boolean)) or
        (not math.isnan(alcohol_current_boolean))
    ):
        if (
            consumer_boolean or
            former_boolean or
            comparison_boolean or
            alcohol_current_boolean
        ):
            alcohol_previous = 1
        else:
            alcohol_previous = 0
    else:
        alcohol_previous = float("nan")
    # Return information.
    return alcohol_previous


# TODO: code alcohol_none as True / False / None
def determine_alcohol_none(
    alcohol_frequency=None,
    alcohol_current=None,
    alcohol_previous=None,
):
    """
    Translate information from UK Biobank about whether person
    never consumes any alcohol.

    Accommodate inexact float values.

    arguments:
        alcohol_frequency (float): ordinal representation of person's frequency
            of alcohol consumption, derivation of UK Biobank field 1558
        alcohol_current (float): binary representation of whether person
            consumes alcohol currently
        alcohol_previous (float): binary representation of whether person
            consumed alcohol previously

    raises:

    returns:
        (bool): whether person never consumes any alcohol

    """

    # Interpret consumption frequency.
    if (not math.isnan(alcohol_frequency)):
        # The variable has a valid value.
        if (-0.5 <= alcohol_frequency and alcohol_frequency < 0.5):
            # The person never consumes alcohol currently.
            frequency_boolean = True
        else:
            frequency_boolean = False
    else:
        frequency_boolean = float("nan")
    # Interpret current consumption.
    if (not math.isnan(alcohol_current)):
        # The variable has a valid value.
        if (-0.5 <= alcohol_current and alcohol_current < 0.5):
            # The person never consumes alcohol currently.
            current_boolean = True
        else:
            current_boolean = False
    else:
        current_boolean = float("nan")
    # Interpret previous consumption.
    if (not math.isnan(alcohol_previous)):
        # The variable has a valid value.
        if (-0.5 <= alcohol_previous and alcohol_previous < 0.5):
            # The person never consumed alcohol previously.
            previous_boolean = True
        else:
            previous_boolean = False
    else:
        previous_boolean = float("nan")
    # Integrate information from multiple variables.
    if (
        (not math.isnan(frequency_boolean)) and
        (not math.isnan(current_boolean)) and
        (not math.isnan(previous_boolean))
    ):
        if (
            frequency_boolean and
            current_boolean and
            previous_boolean
        ):
            # Person has never consumed alcohol, previously or currently.
            alcohol_none = 1
        else:
            # Person has consumed alcohol.
            alcohol_none = 0
    else:
        alcohol_none = float("nan")
    # Return information.
    return alcohol_none


def organize_alcohol_consumption_frequency_variables(
    table=None,
    report=None,
):
    """
    Organizes information about previous and current alcohol consumption.

    "frequency", field 1558: "Alcohol intake frequency"
    UK Biobank data coding 100402 for variable field 1558.
    "daily or almost daily": 1
    "three or four times a week": 2
    "once or twice a week": 3
    "one to three times a month": 4
    "special occasions only": 5
    "never": 6
    "prefer not to answer": -3

    "former", field 3731: "Former alcohol drinker"
    Variable 3731 was only collected for persons who never consume any alcohol
    currently (UK Biobank variable 1558).
    UK Biobank data coding 100352 for variable field 3731.
    "yes": 1
    "no": 0
    "prefer not to answer": -3

    "consumer", field 20117: "Alcohol drinker status"
    Variable 20117 is a derivation of variables 1558 and 3731.
    UK Biobank data coding 90 for variable field 20117.
    "current": 2
    "previous": 1
    "never": 0
    "prefer not to answer": -3

    "comparison", field 1628: "Alcohol intake versus 10 years previously"
    Variable 1628 was only collected for persons who consume alcohol currently
    (UK Biobank variable 1558).
    UK Biobank data coding 100417 for variable field 1628.
    "more nowadays": 1
    "about the same": 2
    "less nowadays": 3
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
        "1558-0.0", "3731-0.0", "1628-0.0", "20117-0.0",
    ]
    table = convert_table_columns_variables_types_float(
        columns=columns_type,
        table=table,
    )
    # Determine person's frequency of alcohol consumption.
    table["alcohol_frequency"] = table.apply(
        lambda row:
            determine_alcohol_consumption_frequency(
                frequency=row["1558-0.0"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person consumes alcohol currently.
    table["alcohol_current"] = table.apply(
        lambda row:
            determine_current_alcohol_consumption(
                consumer=row["20117-0.0"],
                alcohol_frequency=row["alcohol_frequency"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person never consumes any alcohol.
    table["alcohol_previous"] = table.apply(
        lambda row:
            determine_previous_alcohol_consumption(
                consumer=row["20117-0.0"],
                former=row["3731-0.0"],
                comparison=row["1628-0.0"],
                alcohol_current=row["alcohol_current"],
            ),
        axis="columns", # apply across rows
    )
    # Determine whether person never consumes any alcohol, either currently or
    # previously.
    table["alcohol_none"] = table.apply(
        lambda row:
            determine_alcohol_none(
                alcohol_frequency=row["alcohol_frequency"],
                alcohol_current=row["alcohol_current"],
                alcohol_previous=row["alcohol_previous"]
            ),
        axis="columns", # apply across rows
    )
    # Remove columns for variables that are not necessary anymore.
    table_clean = table.copy(deep=True)
    table_clean.drop(
        labels=["1558-0.0", "3731-0.0", "1628-0.0", "20117-0.0",],
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
            "alcohol_previous",
            "alcohol_current",
            "alcohol_none",
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
    table = convert_table_columns_variables_types_float(
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
    pail_quantity = organize_alcohol_consumption_quantity_variables(
        table=pail_consumption["table_clean"],
        report=report,
    )
    # Collect information.
    pail = dict()
    pail["consumption"] = pail_consumption
    pail["quantity"] = pail_quantity
    # Return information.
    return pail


##########
# Alcohol AUDIT questionnaire


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
    table = convert_table_columns_variables_types_float(
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
    table = convert_table_columns_variables_types_float(
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
# Cohort selection: general


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
        inplace=True
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
    variables=None,
    prefixes=None,
    table=None,
):
    """
    Selects records for males with sex-specific criteria.

    arguments:
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
        inplace=True
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
    # Return information.
    return table


def select_records_by_sex_specific_valid_variables_values(
    female=None,
    female_pregnancy=None,
    female_menopause_binary=None,
    female_menopause_ordinal=None,
    female_variables=None,
    female_prefixes=None,
    male=None,
    male_variables=None,
    male_prefixes=None,
    table=None,
):
    """
    Selects records by sex and by sex-specific criteria and variables.

    arguments:
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
        male_variables (list<str>): names of columns for variables in which rows
            must have valid values to keep records for males
        male_prefixes (list<str>): prefixes of columns for variables in which
            rows must have valid values to keep records for males
        table (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    """

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
            variables=male_variables,
            prefixes=male_prefixes,
            table=table,
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
        drop=True,
    )
    table_collection.set_index(
        "eid",
        append=False,
        drop=True,
        inplace=True
    )
    # Return information.
    return table_collection


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


##########
# Cohort selection: sexes, hormones
# For GWAS of hormones


def select_organize_plink_cohorts_variables_by_sex_hormone(
    hormone=None,
    table=None,
):
    """
    Organizes tables for specific cohorts in format for GWAS in PLINK.

    Notice that different variables are relevant for females and males.
    Combination of tables for females and males introduces missing or null
    values for these sex-specific variables.
    Remove records with null values separately for females and males but not
    after combination.

    arguments:
        hormone (str): name of column for hormone variable
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Compile information.
    pail = dict()
    # Select and organize variables across cohorts.
    # Translate variable encodings and table format for analysis in PLINK.

    # Cohort: non-pregnant females and males together

    table_female_male = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=True,
            male_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                hormone,
            ],
            male_prefixes=["genotype_pc_",],
            table=table,
    ))
    pail[str("table_female_male_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_male,
    ))

    # Cohort: all non-pregnant females together

    table_female = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female,
    ))

    table_female_combination = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_", "menopause_hormone_category_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_combination_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_combination,
    ))

    # Cohort: premenopausal females by binary menopause definition

    table_female_premenopause_binary = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0,],
            female_menopause_ordinal=[0, 1, 2],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "menstruation_days",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_premenopause_binary_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_premenopause_binary,
    ))

    # Cohort: postmenopausal females by binary menopause definition

    table_female_postmenopause_binary = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_postmenopause_binary_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_postmenopause_binary,
    ))

    # Cohort: premenopausal females by ordinal menopause definition

    table_female_premenopause_ordinal = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "menstruation_days",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_premenopause_ordinal_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_premenopause_ordinal,
    ))

    # Cohort: perimenopausal females by ordinal menopause definition

    table_female_perimenopause_ordinal = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[1,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "menstruation_days",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_perimenopause_ordinal_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_perimenopause_ordinal,
    ))

    # Cohort: postmenopausal females by ordinal menopause definition

    table_female_postmenopause_ordinal = (
        select_records_by_sex_specific_valid_variables_values(
            female=True,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[2,],
            female_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                "pregnancy", "menopause_binary", "menopause_ordinal",
                "hormone_alteration",
                hormone,
            ],
            female_prefixes=["genotype_pc_",],
            male=False,
            male_variables=[],
            male_prefixes=[],
            table=table,
    ))
    pail[str("table_female_postmenopause_ordinal_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_female_postmenopause_ordinal,
    ))

    # Cohort: males

    table_male = (
        select_records_by_sex_specific_valid_variables_values(
            female=False,
            female_pregnancy=[0,],
            female_menopause_binary=[0, 1,],
            female_menopause_ordinal=[0, 1, 2,],
            female_variables=[],
            female_prefixes=[],
            male=True,
            male_variables=[
                "eid", "IID",
                "sex", "sex_text", "age", "body_mass_index_log",
                hormone,
            ],
            male_prefixes=["genotype_pc_",],
            table=table,
    ))
    pail[str("table_male_" + hormone)] = (
        organize_phenotype_covariate_table_plink_format(
            boolean_phenotypes=[],
            binary_phenotypes=[],
            continuous_variables=[hormone],
            remove_null_records=False,
            table=table_male,
    ))
    # Return information.
    return pail


def select_organize_plink_cohorts_by_sex_hormones(
    table=None,
    report=None,
):
    """
    Organizes tables for specific cohorts in format for GWAS in PLINK.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Compile information.
    pail = dict()
    # Select and organize variables across cohorts.

    pail_albumin = select_organize_plink_cohorts_variables_by_sex_hormone(
        hormone="albumin_log",
        table=table,
    )
    pail.update(pail_albumin)

    pail_steroid_globulin = (
        select_organize_plink_cohorts_variables_by_sex_hormone(
            hormone="steroid_globulin_log",
            table=table,
    ))
    pail.update(pail_steroid_globulin)

    pail_oestradiol = select_organize_plink_cohorts_variables_by_sex_hormone(
        hormone="oestradiol_log",
        table=table,
    )
    pail.update(pail_oestradiol)

    pail_oestradiol_free = (
        select_organize_plink_cohorts_variables_by_sex_hormone(
            hormone="oestradiol_free_log",
            table=table,
    ))
    pail.update(pail_oestradiol_free)

    pail_testosterone = select_organize_plink_cohorts_variables_by_sex_hormone(
        hormone="testosterone_log",
        table=table,
    )
    pail.update(pail_testosterone)

    pail_testosterone_free = (
        select_organize_plink_cohorts_variables_by_sex_hormone(
            hormone="testosterone_free_log",
            table=table,
    ))
    pail.update(pail_testosterone_free)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "report: select_organize_plink_cohorts_by_sex_hormones()"
        )
        for table_name in pail.keys():
            utility.print_terminal_partition(level=5)
            print(table_name)
            print(
                "Count records: " + str(pail[table_name].shape[0])
            )
    # Return information.
    return pail




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


def plot_variable_series_histogram(
    series=None,
    bins=None,
    file=None,
    path_directory=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        comparison (dict): information for chart
        path_directory (str): path for directory

    raises:

    returns:

    """

    # Specify path to directory and file.
    path_file = os.path.join(
        path_directory, file
    )

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Report.
    utility.print_terminal_partition(level=1)
    print(file)
    print(len(series))
    # Create figure.
    figure = plot.plot_distribution_histogram(
        series=series,
        name="",
        bin_method="count",
        bin_count=bins,
        label_bins="values",
        label_counts="counts of persons per bin",
        fonts=fonts,
        colors=colors,
        line=False,
        position=1,
        text="",
    )
    # Write figure.
    plot.write_figure_png(
        path=path_file,
        figure=figure
    )

    pass


def organize_plot_variable_histogram_summary_charts(
    table=None,
    paths=None,
):
    """
    Organizes information about alcoholism from the AUDIT-C questionnaire.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (dict): collection of information about quantity of alcohol consumption

    """

    # Copy data.
    table = table.copy(deep=True)

    # Specify directories and files.
    # Create figures.
    plot_variable_series_histogram(
        series=table["alcohol_frequency"].dropna().to_list(),
        bins=6,
        file="histogram_alcohol_frequency.png",
        path_directory=paths["plot"],
    )
    plot_variable_series_histogram(
        series=table["alcohol_previous"].dropna().to_list(),
        bins=4,
        file="histogram_alcohol_previous.png",
        path_directory=paths["plot"],
    )
    plot_variable_series_histogram(
        series=table["alcohol_drinks_monthly"].dropna().to_list(),
        bins=50,
        file="histogram_alcohol_drinks_monthly.png",
        path_directory=paths["plot"],
    )
    plot_variable_series_histogram(
        series=table["alcoholism"].dropna().to_list(),
        bins=15,
        file="histogram_alcoholism.png",
        path_directory=paths["plot"],
    )

    pass




##########
# Plot


def plot_variable_values_histogram(
    name=None,
    array=None,
    bins=None,
):
    """
    Plots charts from the analysis process.

    arguments:
        name (str): name for plot
        array (object): NumPy array of values to bin and plot in histogram
        bins (int): count of bins for histogram

    raises:

    returns:

    """

    # Collect information about plot.
    pail = dict()
    pail["name"] = name
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
    pail["figure"] = plot.plot_distribution_histogram(
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
        label_text=name, # ""
        label_report=True,
    )
    # Return.
    return pail


##########
# Write


def write_product_quality(
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
    path_table_auditc = os.path.join(
        path_parent, "table_auditc.tsv"
    )
    path_table_audit = os.path.join(
        path_parent, "table_audit.tsv"
    )
    path_table_diagnosis = os.path.join(
        path_parent, "table_diagnosis.tsv"
    )
    path_table_alcoholism = os.path.join(
        path_parent, "table_alcoholism.tsv"
    )
    # Write information to file.
    information["table_auditc"].to_csv(
        path_or_buf=path_table_auditc,
        sep="\t",
        header=True,
        index=False,
    )
    information["table_audit"].to_csv(
        path_or_buf=path_table_audit,
        sep="\t",
        header=True,
        index=False,
    )
    information["table_diagnosis"].to_csv(
        path_or_buf=path_table_diagnosis,
        sep="\t",
        header=True,
        index=False,
    )
    information["table_alcoholism"].to_csv(
        path_or_buf=path_table_alcoholism,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_cohort_table(
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
        path_parent, str(name + ".tsv")
    )
    # Write information to file.
    information.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_cohorts(
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
        write_product_cohort_table(
            name=name,
            information=information[name],
            path_parent=path_parent,
        )
    pass


def write_product_trial(
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
        path_parent, "table_phenotypes_covariates.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        path_parent, "table_phenotypes_covariates.tsv"
    )
    # Write information to file.
    information["table_phenotypes_covariates"].to_pickle(
        path_table_phenotypes
    )
    information["table_phenotypes_covariates"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=False,
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

    # Quality control reports.
    write_product_quality(
        information=information["quality"],
        path_parent=paths["quality"],
    )
    # Cohort tables in PLINK format.
    write_product_cohorts(
        information=information["cohorts"],
        path_parent=paths["cohorts"],
    )

    # Trial organization.
    if False:
        write_product_trial(
            information=information["trial"],
            path_parent=paths["trial"],
        )
    pass


###############################################################################
# Procedure


def execute_genotype_sex_age_body(
    table=None,
    selection=None,
    report=None,
):
    """
    Organizes information about persons' genotypes, sex, age, and body mass
    index across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        selection (str): type of table to select and return
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank

    """

    # Organize information about genotype principal components.
    table_genotype = organize_genotype_principal_component_variables(
        table=table,
        report=report,
    )
    # Organize information about general attributes.
    pail_sex_age_body = organize_sex_age_body_variables(
        table=table_genotype,
        report=report,
    )

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_genotype_sex_age_body()")
        utility.print_terminal_partition(level=3)
        print(pail_sex_age_body[selection])
    # Return information.
    return pail_sex_age_body[selection]


def execute_sex_hormones(
    table=None,
    selection=None,
    report=None,
):
    """
    Organizes information about persons' sex hormones across UK Biobank.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        selection (str): type of table to select and return
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank

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
        print(pail_hormone[selection])
    # Return information.
    return pail_hormone[selection]


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


def execute_plot_hormones(
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
        (object): Pandas data frame of phenotype variables across UK Biobank

    """

    # Copy information.
    table = table.copy(deep=True)
    # Collect information for plots.
    pail = dict()
    # All persons in UK Biobank.
    table_not_pregnant = table.loc[
        (table["pregnancy"] < 0.5), :
    ]
    pail = organize_plot_cohort_sex_hormone_variable_distributions(
        prefix="",
        bins=50,
        table=table_not_pregnant,
    )
    # Filter to females.
    table_female = table.loc[
        (table["sex_text"] == "female"), :
    ]
    table_female_not_pregnant = table_female.loc[
        (table_female["pregnancy"] == 0), :
    ]
    pail_female = organize_plot_cohort_sex_hormone_variable_distributions(
        prefix="female",
        bins=50,
        table=table_female_not_pregnant,
    )
    pail.update(pail_female)
    # Stratify not pregnant females by recent use of oral contraception.
    table_contraception_no = table_female_not_pregnant.loc[
        (table_female_not_pregnant["oral_contraception"] == 0), :
    ]
    pail_contraception_no = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="contraception_no",
            bins=50,
            table=table_contraception_no,
    ))
    pail.update(pail_contraception_no)
    table_contraception_yes = table_female_not_pregnant.loc[
        (table_female_not_pregnant["oral_contraception"] == 1), :
    ]
    pail_contraception_yes = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="contraception_yes",
            bins=50,
            table=table_contraception_yes,
    ))
    pail.update(pail_contraception_yes)
    # Stratify by pre-menopausal, not pregnant females by recent use of
    # hormone replacement therapy.
    table_replacement_no = table_female_not_pregnant.loc[
        (table_female_not_pregnant["hormone_replacement"] == 0), :
    ]
    pail_replacement_no = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="replacement_no",
            bins=50,
            table=table_replacement_no,
    ))
    pail.update(pail_replacement_no)
    table_replacement_yes = table_female_not_pregnant.loc[
        (table_female_not_pregnant["hormone_replacement"] == 1), :
    ]
    pail_replacement_yes = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="replacement_yes",
            bins=50,
            table=table_replacement_yes,
    ))
    pail.update(pail_replacement_yes)

    # Filter to pre-menopausal, not pregnant females.
    table_premenopause = table_female_not_pregnant.loc[
        (table_female_not_pregnant["menopause"] == 0), :
    ]
    pail_premenopause = organize_plot_cohort_sex_hormone_variable_distributions(
        prefix="pre-menopause",
        bins=50,
        table=table_premenopause,
    )
    pail.update(pail_premenopause)
    # Plot days since previous menstrual period for pre-menopausal females.
    pail_menstruation = dict()
    pail_menstruation["pre-menopause_menstruation_day"] = (
        plot_variable_values_histogram(
            name="pre-menopause_menstruation_day",
            array=table_premenopause["menstruation_day"].dropna().to_numpy(),
            bins=50,
    ))
    pail.update(pail_menstruation)
    # Filter to post-menopausal, not pregnant females.
    table_postmenopause = table_female_not_pregnant.loc[
        (table_female_not_pregnant["menopause"] == 1), :
    ]
    pail_postmenopause = organize_plot_cohort_sex_hormone_variable_distributions(
        prefix="post-menopause",
        bins=50,
        table=table_postmenopause,
    )
    pail.update(pail_postmenopause)

    # Any hormone alteration (oral contraception or hormone-replacement).
    table_premenopause_alteration_yes = table_premenopause.loc[
        (table_premenopause["hormone_alteration"] == 1), :
    ]
    pail_premenopause_alteration_yes = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="premenopause_hormone_alteration_yes",
            bins=50,
            table=table_premenopause_alteration_yes,
    ))
    pail.update(pail_premenopause_alteration_yes)
    table_premenopause_alteration_no = table_premenopause.loc[
        (table_premenopause["hormone_alteration"] == 0), :
    ]
    pail_premenopause_alteration_no = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="premenopause_hormone_alteration_no",
            bins=50,
            table=table_premenopause_alteration_no,
    ))
    pail.update(pail_premenopause_alteration_no)

    table_postmenopause_alteration_yes = table_postmenopause.loc[
        (table_postmenopause["hormone_alteration"] == 1), :
    ]
    pail_postmenopause_alteration_yes = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="postmenopause_hormone_alteration_yes",
            bins=50,
            table=table_postmenopause_alteration_yes,
    ))
    pail.update(pail_postmenopause_alteration_yes)
    table_postmenopause_alteration_no = table_postmenopause.loc[
        (table_postmenopause["hormone_alteration"] == 0), :
    ]
    pail_postmenopause_alteration_no = (
        organize_plot_cohort_sex_hormone_variable_distributions(
            prefix="postmenopause_hormone_alteration_no",
            bins=50,
            table=table_postmenopause_alteration_no,
    ))
    pail.update(pail_postmenopause_alteration_no)

    # Filter to pregnant females.
    table_pregnancy = table_female.loc[
        (table_female["pregnancy"] == 1), :
    ]
    pail_pregnancy = organize_plot_cohort_sex_hormone_variable_distributions(
        prefix="pregnancy",
        bins=10,
        table=table_pregnancy,
    )
    pail.update(pail_pregnancy)

    # Filter to males.
    table_male = table.loc[
        (table["sex_text"] == "male"), :
    ]
    pail_male = organize_plot_cohort_sex_hormone_variable_distributions(
        prefix="male",
        bins=50,
        table=table_male,
    )
    pail.update(pail_male)

    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_plot_hormones()")
        utility.print_terminal_partition(level=3)
    # Return information.
    return pail


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

    # Organize Alchol Use Disorders Identification Test (AUDIT) questionnaire
    # variables, including separate scores for AUDIT-Consumption (AUDIT-C) and
    # AUDIT-Problem (AUDIT-P) portions of the questionnaire.
    pail_audit = organize_alcohol_audit_questionnaire_variables(
        table=pail_alcohol_consumption["quantity"]["table_clean"],
        report=False,
    )
    #print(pail_audit["audit"]["table_clean"])

    if False:
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

    # Copy information.
    table_alcohol = pail_audit["audit"]["table_clean"].copy(deep=True)
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("Report from organize_alcohol_consumption()")
        utility.print_terminal_partition(level=3)
        print(table_alcohol)
    # Return information.
    return table_alcohol


def execute_mental_health(
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
        (object): Pandas data frame of phenotype variables across UK Biobank

    """

    # Organize information about neuroticism.
    pail_neuroticism = organize_neuroticism_variables(
        table=table,
        report=report,
    )

    # Copy information.
    table_mental_health = pail_neuroticism["table_clean"].copy(deep=True)
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("report: execute_genotype_sex_age_body()")
        utility.print_terminal_partition(level=3)
        print(table_mental_health)
    # Return information.
    return table_mental_health
