#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

################################################################################
################################################################################
################################################################################
# Notes
# 1. execute relevant procedures from Python program within virtual environment
################################################################################
################################################################################
################################################################################

################################################################################
# Organize argument variables.

path_environment=${1} # full path to Python virtual environment
path_uk_biobank=${2} # full path to parent directory for uk_biobank package
path_dock=${3} # full path to parent directory for source and product files

################################################################################
# Organize paths.

# Initialize dock parent directory.
if [ ! -d $path_dock ]; then
    mkdir -p $path_dock
fi

################################################################################
# Report.

echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "----------"
echo "----------"
echo "----------"
echo "Preparing to execute main procedures from 'uk_biobank' Python package."
echo "Modules 'assembly' and 'organization'."
echo "From repository 'uk_biobank' by 'tcameronwaller' on GitHub."
echo "I'll hopefully make this message look more fancy and official at some point... :/"
echo "----------"
echo "----------"
echo "----------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"
echo "--------------------------------------------------"

################################################################################
# Activate Python Virtual Environment.
source "${path_environment}/bin/activate"
echo "----------"
echo "confirm activation of Python virtual environment..."
which python3
sleep 5s

################################################################################
# Execute procedure.

# Echo each command to console.
set -x

# Execute Python procedure(s).
python3 ${path_uk_biobank}/uk_biobank/interface.py main --path_dock $path_dock --assembly
python3 ${path_uk_biobank}/uk_biobank/interface.py main --path_dock $path_dock --organization

################################################################################
# Deactivate Python Virtual Environment.
deactivate
echo "----------"
echo "confirm deactivation of Python virtual environment..."
which python3
