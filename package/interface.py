"""
Organizes parameters and manages execution of procedures within the 'uk_biobank'
package.

Execution of this subpackage, 'uk_biobank' occurs under the management of a
higher level package. Importation paths represent this hierarchy.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of UK_Biobank
    (https://github.com/tcameronwaller/uk_biobank/).

    UK_Biobank supports analyses on data from the U.K. Biobank.
    Copyright (C) 2022 Thomas Cameron Waller

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

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.
import uk_biobank.assembly
import uk_biobank.importation
import uk_biobank.organization
import uk_biobank.stratification
import uk_biobank.description
import uk_biobank.regression
import uk_biobank.collection

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Package's main subparser.


def define_subparser_main(subparsers=None):
    """
    Defines subparser and parameters.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to subparser

    """

    # Define description.
    description = uk_biobank.interface.define_description_main()
    # Define epilog.
    epilog = uk_biobank.interface.define_epilog_main()
    # Define parser.
    parser = subparsers.add_parser(
        name="uk_biobank",
        description=description,
        epilog=epilog,
        help="Help for 'uk_biobank' routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser.add_argument(
        "-path_dock", "--path_dock", dest="path_dock", type=str, required=True,
        help=(
            "Path to dock directory for source and product " +
            "directories and files."
        )
    )
    parser.add_argument(
        "-assembly", "--assembly", dest="assembly",
        action="store_true",
        help=(
            "Assemble phenotype information from UK Biobank."
        )
    )
    parser.add_argument(
        "-importation", "--importation", dest="importation",
        action="store_true",
        help=(
            "Import phenotype information from UK Biobank."
        )
    )
    parser.add_argument(
        "-organization", "--organization", dest="organization",
        action="store_true",
        help=(
            "Organize phenotype information from UK Biobank."
        )
    )
    parser.add_argument(
        "-stratification", "--stratification",
        dest="stratification",
        action="store_true",
        help=(
            "Stratification of cohorts and formatting tables of phenotypes " +
            "and covariates for genetic analyses (especially GWAS in PLINK2)."
        )
    )
    parser.add_argument(
        "-description", "--description",
        dest="description",
        action="store_true",
        help=(
            "Description of cohorts and phenotypes with summary statistics " +
            "and plots."
        )
    )
    parser.add_argument(
        "-regression", "--regression",
        dest="regression",
        action="store_true",
        help=(
            "Regression analyses of phenotypes within cohorts."
        )
    )
    parser.add_argument(
        "-collection", "--collection",
        dest="collection",
        action="store_true",
        help=(
            "Collection and summary of reports from genotypic analyses."
        )
    )
    # Define behavior.
    parser.set_defaults(func=uk_biobank.interface.evaluate_parameters_main)
    # Return parser.
    return parser


def define_description_main():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Package's main procedure

        Do stuff.

        --------------------------------------------------
    """)
    return description


def define_epilog_main():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        main routine

        --------------------------------------------------
        additional notes...


        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def evaluate_parameters_main(arguments):
    """
    Evaluates parameters for model procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to main routine ...")
    # Execute procedure.
    if arguments.assembly:
        # Report status.
        print("... executing 'assembly' procedure ...")
        # Execute procedure.
        uk_biobank.assembly.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.importation:
        # Report status.
        print("... executing 'importation' procedure ...")
        # Execute procedure.
        uk_biobank.importation.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.organization:
        # Report status.
        print("... executing 'organization' procedure ...")
        # Execute procedure.
        uk_biobank.organization.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.stratification:
        # Report status.
        print("... executing 'stratification' procedure ...")
        # Execute procedure.
        uk_biobank.stratification.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.description:
        # Report status.
        print("... executing 'description' procedure ...")
        # Execute procedure.
        uk_biobank.description.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.regression:
        # Report status.
        print("... executing 'regression' procedure ...")
        # Execute procedure.
        uk_biobank.regression.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.collection:
        # Report status.
        print("... executing 'collection' procedure ...")
        # Execute procedure.
        uk_biobank.collection.execute_procedure(
            path_dock=arguments.path_dock
        )
    pass



###############################################################################
# Procedure



#
