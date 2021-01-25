import pandas as pd


import numpy as np
from smurf2_headers import *
from smurf2_utills import *
import os
import traceback
import pandas as pd
from optparse import OptionParser, OptionGroup
import sys
import matplotlib.pylab as pylab
params = {'legend.fontsize': 12,
          'figure.figsize': (15, 5),
         'axes.labelsize': 20,
         'axes.titlesize':20,
         'xtick.labelsize':15,
         'ytick.labelsize':15}
pylab.rcParams.update(params)
SIMILAR = 'Identical'
CONTAIN = 'All subsets'
RECALL = 'Recall'
PRECISION = 'Precision'

ALL='all'
FROM_DB ='Intra DB'
MODIFIED = 'Modified bacteria'
ORIGIN = 'Intra DB - origin'

VALIDATION_THRESHOLD_THRESHOLD = 0.00001
EXPECTED_RES_FILE_NAME = "expected_res.csv"
ACTUAL_RES_FILE_NAME = "emirge_smurf_WFalseSTrue.csv"
ACTUAL_RES_FILE_NAME_NEW = "SMURF2_results.csv"

class Header():
    ref_id = 'ref'
    prior = 'prior'
    sequence = 'sequence'
    region = 'region'
    weight = 'weight'
    is_changed = 'is_changed'
    is_changed_source = 'is_changed_source'
    all = [ref_id, prior, sequence, region, weight, is_changed, is_changed_source]
    default = [ref_id, prior, sequence, region, weight]



def validate_priors(df, threshold=VALIDATION_THRESHOLD_THRESHOLD):
    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])
    df.prior = df.prior.apply(lambda r: r / sum_prior)

    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])

    # logging.debug("sum of priors is {}".format(sum_prior))
    if abs(sum_prior - 1) > threshold:
        raise Exception("sum of prior is not 1")


def get_test_df(path):
    """
    :param path: path to 'final_results.csv produced by smurf2.py
    :return: df hold the final results
    """
    df = pd.read_csv(path, index_col=None)
    df = df.rename(columns = {'Sequence': Header.sequence,
                              HeadersFormat.Region: Header.region,
                              HeadersFormat.Priors: Header.prior,
                              'Unique_Reference_id': Header.ref_id})
    df = df.drop_duplicates([Header.ref_id, Header.region])
    df[Header.weight] = df.groupby(Header.ref_id)[Header.region].transform('count')

    df = df[Header.default]
    validate_priors(df)

    return df


def get_expected_df(path):
    """
    :param path: path to
    :return: path to expected_res.csv produced by mock_creator.py
    """
    df = pd.read_csv(path, index_col=None)
    df = df.rename(columns={'sequence': Header.sequence,
                            'region': Header.region,
                            'prior': Header.prior,
                            'id': Header.ref_id})
    df[Header.weight] = df.groupby(Header.ref_id)[Header.region].transform('count')
    df[Header.is_changed] = df[Header.ref_id].apply(lambda id: '#' in str(id))
    changed_id = df[df[Header.is_changed] == True][Header.ref_id].tolist()
    changed_id_source  = [id.replace('#', '') for id in changed_id]
    df[Header.is_changed_source] = df[Header.ref_id].apply(lambda id: id in changed_id_source)
    df = df[Header.all]
    validate_priors(df)

    return df


def get_cmd_arguments(argv = sys.argv[1:]):
    USAGE = \
        """usage: %prog WORKING_DIR [required_options] [options]

        Compare the expected results to the actual results of emirge_smurf
        """

    parser = OptionParser(USAGE)

    (options, args) = parser.parse_args(argv)

    if len(args) != 1:
        parser.error(
            "WORKING_DIR is required, and all options except DIR should have a flag associated with them (options without flags: %s)" % args)
        return

    return options, args


def find_dist(s1, s2):
    return sum([1 for b1, b2 in zip(s1, s2) if b1 != b2])


def to_set(x):
    return set(x)


def calc_recall(best_match, max_distance=0, is_change=None, matching_type='is_subset'):
    if is_change is not None:
        best_match = best_match[best_match.is_changed == is_change]
        sum_prior = best_match.drop_duplicates('ref_gt')['prior_gt'].sum()
        best_match['prior_gt'] = best_match['prior_gt']/sum_prior
    best_match = best_match[best_match.ref_to_ref_dist <= max_distance]
    best_match_for_recall = best_match.groupby(['ref_gt', 'prior_gt'])[matching_type].sum().reset_index()
    return (best_match_for_recall['prior_gt']*(best_match_for_recall[matching_type] > 0)).sum()


def calc_precision(best_match, max_distance=0, is_change=None, matching_type='is_subset'):
    if is_change is not None:
        best_match = best_match[best_match.is_changed == is_change]
        sum_prior = best_match.drop_duplicates('ref_test')['prior_test'].sum()
        best_match['prior_test'] = best_match['prior_test'] / sum_prior
    best_match = best_match[best_match.ref_to_ref_dist <= max_distance]
    best_match_for_precision = best_match.groupby(['ref_test', 'prior_test'])[matching_type].sum().reset_index()
    return (best_match_for_precision['prior_test']*(best_match_for_precision[matching_type] > 0)).sum()


def find_best_match(test_df, ground_truth_df):
    all_merge_dfs = []

    #  Simple merge exact:
    tmp = ground_truth_df.groupby(Header.ref_id)[Header.region].apply(to_set).reset_index()
    ground_truth_df = ground_truth_df.merge(tmp.rename(columns={Header.region: 'regions_set'}), on=Header.ref_id)

    tmp = test_df.groupby(Header.ref_id)[Header.region].apply(to_set).reset_index()
    test_df = test_df.merge(tmp.rename(columns={Header.region: 'regions_set'}), on=Header.ref_id)

    test_df['test_id_region'] = test_df.apply(lambda r: str(r[Header.ref_id]) + "_" + str(r[Header.region]), axis=1)
    ground_truth_df['gt_id_region'] = ground_truth_df.apply(lambda r: str(r[Header.ref_id]) + "_" + str(r[Header.region]), axis=1)
    merge = test_df.merge(ground_truth_df, on=[Header.region, Header.sequence], how="inner", suffixes=('_test', '_gt'))
    merge['dist'] = 0
    all_merge_dfs.append(merge)

    best_mapping = merge

    test_mapped_ids = best_mapping['test_id_region'].unique()
    unmapped_test = test_df[~test_df['test_id_region'].isin(test_mapped_ids)]

    if not unmapped_test.empty:
        non_similar_merge = unmapped_test.merge(ground_truth_df, on=Header.region, suffixes=('_test', '_gt'))
        non_similar_merge['dist'] = non_similar_merge.apply(lambda r: find_dist(r[Header.sequence + "_test"],
                                                                                r[Header.sequence + "_gt"]), axis=1)
        non_similar_merge = non_similar_merge[non_similar_merge['dist'] ==
                                              non_similar_merge.groupby('test_id_region')['dist'].transform(min)]
        best_mapping = pd.concat([best_mapping, non_similar_merge], ignore_index=True)

    gt_mapped_ids = best_mapping['gt_id_region'].unique()
    unmapped_gt = ground_truth_df[~ground_truth_df['gt_id_region'].isin(gt_mapped_ids)]

    if not unmapped_gt.empty:
        non_similar_merge = test_df.merge(unmapped_gt, on=Header.region, suffixes=('_test', '_gt'))
        non_similar_merge['dist'] = non_similar_merge.apply(lambda r: find_dist(r[Header.sequence + "_test"],
                                                                                r[Header.sequence + "_gt"]), axis=1)
        non_similar_merge = non_similar_merge[non_similar_merge['dist'] ==
                                              non_similar_merge.groupby('gt_id_region')['dist'].transform(min)]
        all_merge_dfs.append(non_similar_merge)

        best_mapping = pd.concat([best_mapping, non_similar_merge], ignore_index=True)

    tmp = best_mapping.groupby(['ref_gt', 'ref_test'])[Header.region].apply(to_set).reset_index()
    best_mapping = best_mapping.merge(tmp.rename(columns={Header.region: 'common_regions_set'}))

    best_mapping['is_subset'] = best_mapping['regions_set_test'] == best_mapping['common_regions_set']
    best_mapping['is_identical'] = best_mapping['regions_set_gt'] == best_mapping['common_regions_set']

    best_mapping['ref_to_ref_dist'] = best_mapping.groupby(['ref_gt', 'ref_test'])['dist'].transform('sum')
    best_match = best_mapping.drop_duplicates(['ref_gt', 'ref_test'])[[u'ref_test',
                                                                       u'prior_test',
                                                                       u'ref_gt',
                                                                       u'prior_gt',
                                                                       u'is_changed',
                                                                       u'is_changed_source',
                                                                       u'is_subset',
                                                                       u'is_identical',
                                                                       u'ref_to_ref_dist']]

    return best_match


def main():
    print 'test evaluation'
    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    working_dir = os.path.abspath(args[0])
    single_experiment_evaluation(working_dir)


def single_experiment_evaluation(test_path):
    experiments_df = None

    actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME_NEW)
    expected_path = os.path.join(test_path, EXPECTED_RES_FILE_NAME)

    test_df = get_test_df(actual_path)
    ground_truth_df = get_expected_df(expected_path)

    best_mapping = find_best_match(test_df, ground_truth_df)

    print "recall identical {}".format(calc_recall(best_mapping, matching_type='is_identical'))
    print "precision identical {}".format(calc_precision(best_mapping, matching_type='is_identical'))
    print "recall identical, intra-db {}".format(calc_recall(best_mapping, matching_type='is_identical', is_change=False))
    print "precision identical, intra-db {}".format(calc_precision(best_mapping, matching_type='is_identical', is_change=False))
    print "recall identical, modified {}".format(calc_recall(best_mapping, matching_type='is_identical', is_change=True))
    print "precision identical, modified {}".format(calc_precision(best_mapping, matching_type='is_identical', is_change=True))

    print "recall subset {}".format(calc_recall(best_mapping))
    print "precision subset {}".format(calc_precision(best_mapping))
    print "recall subset, intra-db {}".format(calc_recall(best_mapping, is_change=False))
    print "precision subset, intra-db {}".format(calc_precision(best_mapping, is_change=False))
    print "recall subset, modified {}".format(calc_recall(best_mapping, is_change=True))
    print "precision subset, modified {}".format(calc_precision(best_mapping, is_change=True))

    print 3

if __name__ == "__main__":
    main()




