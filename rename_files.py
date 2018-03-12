#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rename files based on a csv file with old and new names.

@ V.R.Marcelino
Created on Feb 27 2018
"""

import csv
import pandas as pd
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--old_new_file', help='The path to the file containing old and new names', required=True)

args = parser.parse_args()
rename_file = args.old_new_file

#rename_file = "SraRunTable.csv"

with open(rename_file) as f:
    for line in f:
        old_name = line.split(",")[0]
        new_name = line.split(",")[1]
        os.rename(old_name, new_name)

