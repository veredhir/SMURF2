
from Bio import SeqIO
import numpy as np
from smurf2_headers import *
import os
import scipy.io as sio
from fig2_evaluation_smurf import *

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
# from pandas.tools.plotting import table
from optparse import OptionParser, OptionGroup
import glob
import sys

def find_diff(seq1, seq2):
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            print("{}: {}!={}".format(i, seq1[i], seq2[i]))

SIMILAR = 'similar'
CONTAIN = 'contain'
RECALL = 'Recall'
PRECISION = 'Precision'

VALIDATION_THRESHOLD_THRESHOLD = 0.00001
EXPECTED_RES_FILE_NAME = "expected_res.csv"
ACTUAL_RES_FILE_NAME = "emirge_smurf_WFalseSTrue.csv"


class Header():
    ref_id = 'ref'
    prior = 'prior'
    sequence = 'sequence'
    region = 'region'
    weight = 'weight'
    is_changed = 'is_changed'
    new_id = 'new_id'


def read_db_ind_map(path_to_db_ind="/home/vered/EMIRGE/data/reference_db/Header_uni_forVered.csv"):
    db_map = pd.read_csv(path_to_db_ind, index_col=None, header=None)
    db_map['fasta_id'] = db_map[0].apply(lambda x: str(int(x)))
    return db_map


def to_set(x):
    return set(x)


def find_match_smurf(smurf_output_path, ground_truth_path):
    smurf_results = convert_format(smurf_output_path)
    gt_df = get_expected_df(ground_truth_path)
    if gt_df['ref'].dtype is not float:
        gt_df['ref_set'] = gt_df['ref'].apply(lambda x: set([x]))
    else:
        gt_df['ref_set'] = gt_df['ref'].apply(lambda x: set([float(x.strip('#'))]))

    test = gt_df.merge(smurf_results, on='region', suffixes=("_gt", "_smurf"))
    test.drop_duplicates(['ref_smurf', 'ref_gt'])
    test['is_match'] = test['ref_set'] <= test['smurf ids']

    test = test[test['is_match']]

    print 'recall = ', test.drop_duplicates('ref_gt')['prior_gt'].sum()
    print 'precision = ', test.drop_duplicates('ref_smurf')['prior_smurf'].sum()

    print 3



def convert_format(mat_path):
    """
    :param path: path to 'final_results.csv produced by smurf2.py
    :return: df hold the final results
    """
    matstruct_contents = sio.loadmat(mat_path)

    bacteriaMetaGroups = matstruct_contents['bactMetaGroups']
    found_bacteria = matstruct_contents['found_bacteria']
    comb_vecs=bacteriaMetaGroups['comb_vec'][0]
    db_ind = bacteriaMetaGroups['db_ind'][0]
    frequency = found_bacteria['frequency'][0][0][0]
    df = pd.DataFrame()
    db_map = read_db_ind_map()

    print "db map length = ", len(db_map)

    for comb_vectors, ids, prior, index in zip(comb_vecs, db_ind, frequency, range(len(frequency))):

        db_id = set(db_map.loc[ids[0]-1, 'fasta_id'].astype(float).to_list())
        comb_vector = comb_vectors[0]
        for ix in range(len(comb_vector)):
            if comb_vector[ix] == 1:
                df = df.append({'sequence': '',
                                'region': int(ix+1),
                                'prior': prior,
                                'smurf ids': db_id,
                                'ref': index}, ignore_index=True)

    return df



def main(argv=sys.argv[1:]):
    """
    command line interface to emirge

    """
    parser = OptionParser("Convert SMURF results format to SMURF2 results format")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "")
    group_reqd.add_option("-s", "--smurf_input",
                          type="string",
                          help="path to smurf2 results .mat directory",
                          default="/home/vered/EMIRGE/cluster/cluster_smurf_gary/Results_files_0")
    group_reqd.add_option("-2", "--smurf2_input",
                          type="string",
                          help="path to smurf results .mat directory",
                          default="/home/vered/EMIRGE/cluster/cluster_smurf/0/python_results")
    group_reqd.add_option("-o", "--output",
                          type="string",
                          help="path to output directory",
                          default="/home/vered/EMIRGE/cluster/cluster_smurf")
    parser.add_option_group(group_reqd)
    (options, args) = parser.parse_args(argv)

    smurf_dir = options.smurf_input
    smurf2_dir = options.smurf2_input
    output_dir=options.output

    for file_name in os.listdir(smurf_dir):
        if not file_name.endswith('mat'):
            continue

        test_name = file_name.split("sample_mock_ix_")[1].split('_results.mat')[0]
        smurf_path = os.path.join(smurf_dir, file_name)
        ground_truth = os.path.join(smurf2_dir, 'test_' + test_name, EXPECTED_RES_FILE_NAME)

        find_match_smurf(smurf_path, ground_truth)



if __name__ == "__main__":
    main()



