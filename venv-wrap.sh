#!/bin/sh

# wrap script for execution in the venv

if [[ $# -ne 2 ]] ; then
    echo "Usage: $0 venv-dir script" 1>&2
    exit 1
fi
venv=$1
script=$2

cat > $venv/bin/$script <<EOF
#!${venv}/bin/python3
import os
os.environ['PLATFORM_MODEL'] = '$venv/models/reduced/RandomForestClassifier'
EOF
cat >> $venv/bin/$script < scripts/$script.py
chmod +x $venv/bin/$script
