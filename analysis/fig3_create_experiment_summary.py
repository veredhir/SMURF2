
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
IDENTICAL = 'Identical'
ALL_SUBSETS = 'All subsets'
RECALL = 'Recall'
PRECISION = 'Precision'

ALL = 'all'
INTRA_DB ='Intra DB'
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
    changed_id_source = [id.replace('#', '') for id in changed_id]
    df[Header.is_changed_source] = df[Header.ref_id].apply(lambda id: id in changed_id_source)
    df = df[Header.all]
    validate_priors(df)

    return df



def to_set(x):
    return set(x)


def calc_table(expected_df, actual_df):

    tmp = expected_df.groupby(Header.ref_id)[Header.region].apply(to_set).reset_index()
    expected_df = expected_df.merge(tmp.rename(columns={Header.region: 'regions_set'}), on=Header.ref_id)
    tmp = actual_df.groupby(Header.ref_id)[Header.region].apply(to_set).reset_index()
    actual_df = actual_df.merge(tmp.rename(columns={Header.region: 'regions_set'}), on=Header.ref_id)
    merged = expected_df.merge(actual_df, on=[Header.region, Header.sequence], how="inner", suffixes=('_e', '_a'))

    merged['common_regions'] = merged.groupby([Header.ref_id + "_e", Header.ref_id + '_a'])[Header.region].transform('count')


    # Todo: decide if it should be contains or overlap
    merged['is_contains'] = (merged[Header.weight + '_a'] == merged['common_regions'])
    merged['is_equals'] = (merged[Header.weight + '_e'] == merged['common_regions'])


    score_table_eqaul = pd.pivot_table(merged,
                                       values='is_equals',
                                       index=[Header.ref_id + '_e'],
                                       columns=[Header.ref_id + '_a'],
                                       aggfunc=np.sum,
                                       fill_value=0)
    score_table_contain = pd.pivot_table(merged,
                                         values='is_contains',
                                         index=[Header.ref_id + '_e'],
                                         columns=[Header.ref_id + '_a'],
                                         aggfunc=np.sum,
                                         fill_value=0)


    actual_priors_index = merged.drop_duplicates(Header.ref_id + '_a')[Header.ref_id + '_a']
    filtered_actual = actual_df[actual_df[Header.ref_id].isin(actual_priors_index)].drop_duplicates(Header.ref_id)
    actual_priors = filtered_actual[Header.prior]
    actual_priors.index= filtered_actual[Header.ref_id]

    expected_priors_index = merged.drop_duplicates(Header.ref_id + '_e')[Header.ref_id + '_e']
    filtered_expected = expected_df[expected_df[Header.ref_id].isin(expected_priors_index)].drop_duplicates(
        Header.ref_id)
    expected_priors = filtered_expected[Header.prior]
    expected_priors.index = filtered_expected[Header.ref_id]

    return score_table_contain, score_table_eqaul, actual_priors, expected_priors


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


class Statistics(object):
    def __init__(self, data):
        self.mean = np.mean(data)
        self.std = np.std(data)
        self.min = min(data)
        self.max = max(data)


class SingleExperimentData(object):
    def __init__(self, workdir, index, changed_bacterias, changed_bases, mix_size=500000):
        self.changed_bacterias=changed_bacterias
        self.changed_bases=changed_bases
        self.mix_size = mix_size
        self.mix_size_str = str(mix_size/1000) + 'k'
        self.index= index
        self.results_path = os.path.join(workdir, str(index), "python_results")

    def calc_evaluation(self, best_match):
        dict_list = []
        basic_dict = {'exp-index': self.index,
                      '# Reads': self.mix_size_str,
                      '# Changes per bacterium':  self.changed_bases, 
                      '# Bacteria changed': self.changed_bacterias}
        
        for match_type in [IDENTICAL, ALL_SUBSETS]:
            for bacteria_group in [MODIFIED, ALL, INTRA_DB]:
                curr_dict = basic_dict.copy()
                curr_dict.update({'score-type': match_type,
                                  'value-type': PRECISION,
                                  'value': calc_precision(best_match,
                                                          max_distance=0,
                                                          bacteria_group=bacteria_group,
                                                          matching_type=match_type),
                                  'bacteria': bacteria_group})
                dict_list.append(curr_dict.copy())
                curr_dict.update({'score-type': match_type,
                                  'value-type': RECALL,
                                  'value': calc_recall(best_match,
                                                       max_distance=0,
                                                       bacteria_group=bacteria_group,
                                                       matching_type=match_type),
                                  'bacteria': bacteria_group})
                dict_list.append(curr_dict.copy())
        return pd.DataFrame(dict_list)


def main():
    print 'test evaluation'
    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    working_dir = os.path.abspath(args[0])

    experiments={}
    experiments[0] = SingleExperimentData(working_dir, 0, 1, 1, 500000)
    experiments[1] = SingleExperimentData(working_dir, 1, 1, 10, 500000)
    experiments[2] = SingleExperimentData(working_dir, 2, 1, 20, 500000)
    experiments[3] = SingleExperimentData(working_dir, 3, 3, 1, 500000)
    experiments[4] = SingleExperimentData(working_dir, 4, 3, 10, 500000)
    experiments[5] = SingleExperimentData(working_dir, 5, 3, 20, 500000)
    experiments[6] = SingleExperimentData(working_dir, 6, 3, 50, 500000)
    experiments[7] = SingleExperimentData(working_dir, 7, 3, 100, 500000)
    experiments[8] = SingleExperimentData(working_dir, 8, 3, 10, 100000)
    experiments[9] = SingleExperimentData(working_dir, 9, 3, 10, 50000)
    experiments[10] = SingleExperimentData(working_dir, 10, 10, 10, 500000)
    experiments[11] = SingleExperimentData(working_dir, 11, 20, 10, 500000)

    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacterium', '# Bacteria changed', 'score-type',
                           'value-type', 'value', 'bacteria']
    experiments_df = pd.DataFrame(columns=experiments_df_cols)

    for experiment in experiments.values():
        test_dirs = os.listdir(experiment.results_path)
        for test_dir in test_dirs:
            test_path = os.path.join(experiment.results_path, test_dir)
            eval = single_experiment_evaluation(experiment, test_path)
            if eval is not None:
                experiments_df = experiments_df.append(eval)

        experiments_df.to_csv(os.path.join(working_dir, "experiments_df.csv"), index=False)



def find_dist(s1, s2):
    return sum([1 for b1, b2 in zip(s1, s2) if b1 != b2])


def to_set(x):
    return set(x)


def calc_recall(best_match, max_distance=0, bacteria_group=None, matching_type=ALL_SUBSETS):
    if bacteria_group is not None and (bacteria_group is not ALL):
        is_change = bacteria_group == MODIFIED
        best_match = best_match[best_match.is_changed == is_change]
        sum_prior = best_match.drop_duplicates('ref_gt')['prior_gt'].sum()
        best_match['prior_gt'] = best_match['prior_gt']/sum_prior
    best_match = best_match[best_match.ref_to_ref_dist <= max_distance]
    best_match_for_recall = best_match.groupby(['ref_gt', 'prior_gt'])[matching_type].sum().reset_index()
    return (best_match_for_recall['prior_gt']*(best_match_for_recall[matching_type] > 0)).sum()


def calc_precision(best_match, max_distance=0, bacteria_group=None, matching_type=ALL_SUBSETS):
    if bacteria_group is not None and (bacteria_group is not ALL):
        is_change = bacteria_group == MODIFIED
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

    best_mapping[ALL_SUBSETS] = best_mapping['regions_set_test'] == best_mapping['common_regions_set']
    best_mapping[IDENTICAL] = best_mapping['regions_set_gt'] == best_mapping['common_regions_set']

    best_mapping['ref_to_ref_dist'] = best_mapping.groupby(['ref_gt', 'ref_test'])['dist'].transform('sum')
    best_match = best_mapping.drop_duplicates(['ref_gt', 'ref_test'])[[u'ref_test',
                                                                       u'prior_test',
                                                                       u'ref_gt',
                                                                       u'prior_gt',
                                                                       u'is_changed',
                                                                       u'is_changed_source',
                                                                       ALL_SUBSETS,
                                                                       IDENTICAL,
                                                                       u'ref_to_ref_dist']]

    return best_match


def single_experiment_evaluation(experiment, test_path):
    experiments_df = None
    try:
        actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME)
        expected_path = os.path.join(test_path, EXPECTED_RES_FILE_NAME)
        if not os.path.isfile(actual_path):
            actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME_NEW)
            if not os.path.isfile(actual_path):
                return experiments_df

        test_df = get_test_df(actual_path)
        ground_truth_df = get_expected_df(expected_path)

        best_match = find_best_match(test_df, ground_truth_df)

        experiments_df = experiment.calc_evaluation(best_match)

    except Exception as ex:
        traceback.print_exc()
        print("Failed add data for test ")

    return experiments_df


if __name__ == "__main__":
    # main()
    experiment = SingleExperimentData("/home/vered/EMIRGE/data", 0, 3, 10, 500000)
    test_path = "/home/vered/EMIRGE/data/test3"
    single_experiment_evaluation(experiment, test_path)



