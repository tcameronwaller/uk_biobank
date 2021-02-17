#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script accesses local information about persons and phenotypes in
# the UK Biobank.
# argument 1: path to local file with variable identifiers to access
# --- format is text file with new line delimiters
# argument 2: path to local directory in which to organize product files
###########################################################################
###########################################################################
###########################################################################

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "The script accesses local information about persons and phenotypes in "
echo "the UK Biobank."
echo "version: 1"
echo "----------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"

###########################################################################
# Organize arguments.
path_variables=$1
path_destination=$2

###########################################################################
# Organize paths.
cd ~/paths
path_ukb_phenotype_data=$(<"./ukbiobank_phenotype.txt")
path_exclusion="$path_ukb_phenotype_data/exclude.csv"
path_ukb_parameter=$(<"./ukbiobank_parameter.txt")
path_identifier_pairs="$path_ukb_parameter/link.file.csv"
path_ukb_tools=$(<"./ukbiobank_tools.txt")

###########################################################################
# Copy auxiliary files.

# Copy UK Biobank phenotype variables to destination directory.
cp $path_variables "$path_destination/uk_biobank_access_variables.txt"
cp $path_variables "$path_destination/variables.txt"
# Copy table of identifier pairs to destination directory.
cp $path_identifier_pairs "$path_destination/table_identifier_pairs.csv"
# Copy table of exclusion identifiers to destination directory.
cp $path_exclusion "${path_destination}/list_exclusion_identifiers.txt"

###########################################################################
# Access variables from each phenotype data release of UK Biobank.
# Access names of all current phenotype data releases.
# UK Biobank conversion tool can only accommodate path strings 1-64 characters.
# https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide
# https://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.3.pdf
cd $path_destination
accessions=($(ls -d $path_ukb_phenotype_data/ukb*))
for i in "${accessions[@]}"; do
    echo $i
    dir=`basename $i`
    # Remove log file to avoid error.
    rm $path_ukb_phenotype_data/$dir/$dir.log
    # Convert data to text file with comma ("csv") or tab ("txt") delimiters.
    # The option ("txt") for tab delimiters seems to be malfunctional.
    # The tab delimiter format has different counts of columns in header and
    # body rows.
    $path_ukb_tools/ukbconv $path_ukb_phenotype_data/$dir/$dir.enc_ukb csv -i"./variables.txt" -o"./$dir.raw"
    # Rename product file.
    #mv "./$dir.raw.csv" "./$dir.raw.tsv"
    # Remove log file to avoid error.
    rm $path_ukb_phenotype_data/$dir/$dir.log
done

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "Done."
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
