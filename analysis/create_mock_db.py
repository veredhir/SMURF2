#!/usr/bin/env

import argparse
import os
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np

USAGE=\
"""create_mock_db.py"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--idb', required=True, type="string",
                        help="path to input database [default = %default]",
                        default="/home/vered/EMIRGE/data/reference_db")
    parser.add_argument('--odb', required=True, type="string",
                        help="path to output database [default = %default]",
                        default="/home/vered/EMIRGE/data/s26_mock/reference_mock_db")
    parser.add_argument('-g', "--groups_to_ids_path", required=False, type="string",
                        help="path to Gari's csv file, converting group id to DB IDs [default = %default]",
                        default="/home/vered/EMIRGE/data/reference_db/Header_uni_forVered.csv")
    parser.add_argument('-n', "--groups_to_delete", required=True, type=int,
                        help="Number of groups to delete [default = %default]",
                        default=1)
    parser.add_argument('-m', "--mock_results_path", required=False, type="string",
                        help="path to SMURF reconstruction results (.mat format) [default = %default]",
                        default="/home/vered/EMIRGE/data/s26_mock/SUMRF_samples26.mat")

    (options, args) = parser.parse_args()

    # if os.pathoptions.odb

