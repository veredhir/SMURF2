

import numpy as np
from smurf2_headers import *
from smurf2_utills import *
import os
import seaborn as sns

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
# from pandas.tools.plotting import table
from optparse import OptionParser, OptionGroup
import glob
import sys

SIMILAR = 'similar'
CONTAIN = 'contain'
RECALL = 'Recall'
PRECISION = 'Precision'

VALIDATION_THRESHOLD_THRESHOLD = 0.00001
EXPECTED_RES_FILE_NAME = "expected_res.csv"
ACTUAL_RES_FILE_NAME = "emirge_smurf_WFalseSFalse.csv"

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


def count_diff(seq1, seq2):
    return sum(0 if x == y else 1 for x, y in zip(seq1, seq2))


def validate_priors(df, threshold=VALIDATION_THRESHOLD_THRESHOLD):
    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])
    df.prior = df.prior.apply(lambda r: r / sum_prior)

    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])

    # logging.debug("sum of priors is {}".format(sum_prior))
    if abs(sum_prior - 1) > threshold:
        raise Exception("sum of prior is not 1")


def get_test_df(path, is_smurf2):
    """
    :param path: path to 'final_results.csv produced by smurf2.py
    :return: df hold the final results
    """
    df = pd.DataFrame.from_csv(path, index_col=None)
    df = df.rename(columns = {'Sequence': Header.sequence,
                              HeadersFormat.Region: Header.region,
                              HeadersFormat.Priors: Header.prior,
                              'Unique_Reference_id': Header.ref_id})
    if is_smurf2:
        df = df.drop_duplicates([Header.ref_id, Header.region])
    df[Header.weight] = df.groupby(Header.ref_id)[Header.region].transform('nunique')

    # Fixing gary's results:
    df[Header.sequence] = df[Header.sequence].apply(lambda s: s[:126] + s[-126:])


    df = df[Header.default]
    validate_priors(df)

    return df


def get_expected_df(path):
    """
    :param path: path to
    :return: path to expected_res.csv produced by mock_creator.py
    """
    df = pd.DataFrame.from_csv(path, index_col=None)
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


def is_expected_contain_actual(merged, actual_df, expected_id, actual_id):
    """
    :param merged:
    :param actual_df:
    :param expected_id:
    :param actual_id:
    :return:
    """
    common_regions = len(merged[(merged[Header.ref_id + '_e'] == expected_id) & (merged[Header.ref_id + '_a'] == actual_id)])
    actual_length = len(actual_df[actual_df[Header.ref_id] == actual_id])
    return common_regions == actual_length


def is_actual_equals_expected(merged, expected_df, expected_id, actual_id):
    """
    :param merged:
    :param actual_df:
    :param expected_id:
    :param actual_id:
    :return:
    """
    common_regions = len(merged[(merged[Header.ref_id + '_e'] == expected_id) & (merged[Header.ref_id + '_a'] == actual_id)])
    expected_length = len(expected_df[expected_df[Header.ref_id] == expected_id])
    return common_regions == expected_length



def calc_full_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    actual_df = get_test_df(actual_path, 'gary' not in actual_path)
    return calc_table(expected_df, actual_df)

def calc_changed_bacteria_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    expected_df = expected_df[expected_df[Header.is_changed] == True]
    actual_df = get_test_df(actual_path)
    score_table_contain, score_table_eqaul, actual_priors, expected_priors =  calc_table(expected_df, actual_df)
    sum_priors = sum(actual_priors)
    actual_priors = [prior/sum_priors for prior in actual_priors]
    sum_priors = sum(expected_priors)
    expected_priors = [prior / sum_priors for prior in expected_priors]
    return score_table_contain, score_table_eqaul, actual_priors, expected_priors

def calc_unchagned_bacteria_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    expected_df = expected_df[(expected_df[Header.is_changed] == False) & (expected_df[Header.is_changed_source] == False)]
    actual_df = get_test_df(actual_path)
    score_table_contain, score_table_eqaul, actual_priors, expected_priors = calc_table(expected_df, actual_df)
    sum_priors = sum(actual_priors)
    actual_priors = [prior / sum_priors for prior in actual_priors]
    sum_priors = sum(expected_priors)
    expected_priors = [prior / sum_priors for prior in expected_priors]
    return score_table_contain, score_table_eqaul, actual_priors, expected_priors

def calc_chagned_source_bacteria_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    if expected_df.size == 0:
        return None, None, None, None
    expected_df = expected_df[expected_df[Header.is_changed_source] == True]
    actual_df = get_test_df(actual_path)
    score_table_contain, score_table_eqaul, actual_priors, expected_priors = calc_table(expected_df, actual_df)
    sum_priors = sum(actual_priors)
    actual_priors = [prior / sum_priors for prior in actual_priors]
    sum_priors = sum(expected_priors)
    expected_priors = [prior / sum_priors for prior in expected_priors]
    return score_table_contain, score_table_eqaul, actual_priors, expected_priors



def calc_table_no_pcr_error(expected_df, actual_df):
    # remove the primers
    # print ("removing the primers for calculation")
    PRIMER_SIZE=20
    expected_df[Header.sequence] = expected_df[Header.sequence].apply(lambda s: s[PRIMER_SIZE:-1*PRIMER_SIZE])
    actual_df[Header.sequence] = actual_df[Header.sequence].apply(lambda s: s[PRIMER_SIZE:-1*PRIMER_SIZE])

    expected_df['is_N'] = expected_df[Header.sequence].apply(lambda s: 'N' in s)
    actual_df['is_N'] = actual_df[Header.sequence].apply(lambda s: 'N' in s)

    if (True in expected_df['is_N'].unique().tolist()) or (True in actual_df['is_N'].unique().tolist()) :
        print '***\n***\n'
    merged = expected_df.merge(actual_df, on=[Header.region, Header.sequence], how="inner", suffixes=('_e', '_a'))
    merged = merged.drop_duplicates(['ref_e', 'ref_a', 'region'])

    common_regions_df = merged.groupby([Header.ref_id + "_e", Header.ref_id + '_a'])[
        Header.region].count().reset_index()
    if common_regions_df.size != 0:
        common_regions_df = common_regions_df.merge(
                merged[[Header.ref_id + "_e", Header.ref_id + '_a', Header.weight + '_e', Header.weight + '_a']],
                on=[Header.ref_id + "_e", Header.ref_id + '_a'],
                how="inner")
        common_regions_df['is_contains'] = (common_regions_df[Header.weight + '_a'] == common_regions_df[Header.region])
        common_regions_df['is_equals'] = (common_regions_df[Header.weight + '_e'] == common_regions_df[Header.region])
    else:
        # print ("\n\n the original has no match!")
        common_regions_df = merged.copy()
        common_regions_df['is_contains'] = False
        common_regions_df['is_equals'] = False

    score_table_eqaul = pd.pivot_table(common_regions_df,
                                       values='is_equals',
                                       index=[Header.ref_id + '_e'],
                                       columns=[Header.ref_id + '_a'],
                                       aggfunc=np.sum,
                                       fill_value=0)
    score_table_contain = pd.pivot_table(common_regions_df,
                                         values='is_contains',
                                         index=[Header.ref_id + '_e'],
                                         columns=[Header.ref_id + '_a'],
                                         aggfunc=np.sum,
                                         fill_value=0)


    actual_priors_index = common_regions_df.drop_duplicates(Header.ref_id + '_a')[Header.ref_id + '_a']
    filtered_actual = actual_df[actual_df[Header.ref_id].isin(actual_priors_index)].drop_duplicates(Header.ref_id)
    actual_priors = filtered_actual[Header.prior]
    actual_priors.index= filtered_actual[Header.ref_id]

    expected_priors_index = common_regions_df.drop_duplicates(Header.ref_id + '_e')[Header.ref_id + '_e']
    filtered_expected = expected_df[expected_df[Header.ref_id].isin(expected_priors_index)].drop_duplicates(
        Header.ref_id)
    expected_priors = filtered_expected[Header.prior]
    expected_priors.index = filtered_expected[Header.ref_id]
    #
    # actual_priors = actual_df.drop_duplicates(Header.ref_id)[Header.prior]
    # actual_priors.index = actual_df.drop_duplicates(Header.ref_id)[Header.ref_id]
    # expected_priors = expected_df.drop_duplicates(Header.ref_id)[Header.prior]
    # expected_priors.index = expected_df.drop_duplicates(Header.ref_id)[Header.ref_id]

    return score_table_contain, score_table_eqaul, actual_priors, expected_priors


def calc_table(expected_df, actual_df):
    # remove the primers
    # print ("removing the primers for calculation")
    PRIMER_SIZE=20
    expected_df[Header.sequence] = expected_df[Header.sequence].apply(lambda s: s[PRIMER_SIZE:-1*PRIMER_SIZE])
    actual_df[Header.sequence] = actual_df[Header.sequence].apply(lambda s: s[PRIMER_SIZE:-1*PRIMER_SIZE])

    expected_df['is_N'] = expected_df[Header.sequence].apply(lambda s: 'N' in s)
    actual_df['is_N'] = actual_df[Header.sequence].apply(lambda s: 'N' in s)

    if (True in expected_df['is_N'].unique().tolist()) or (True in actual_df['is_N'].unique().tolist()) :
        print '***\n***\n'
    merged = expected_df.merge(actual_df, on=[Header.region], how="inner", suffixes=('_e', '_a'))
    merged['seq_diff'] = merged.apply(lambda r: count_diff(r[Header.sequence+'_e'], r[Header.sequence+'_a']), axis=1)
    merged = merged[merged['seq_diff'] < 3]
    merged.sort_values(by=['seq_diff'], inplace=True)
    merged = merged.drop_duplicates(['ref_e', 'ref_a', 'region'])


    # test_m = expected_df.merge(actual_df, on=[Header.region, Header.sequence], how="inner", suffixes=('_e', '_a'))
    # merged = test_m.drop_duplicates(['ref_e', 'ref_a', 'region'])
    #
    # merged['seq_diff'].hist(bins=252, grid=False, xlabelsize=12, ylabelsize=12)
    # plt.xlabel("Diff bases", fontsize=15)
    # plt.ylabel("Frequency", fontsize=15)
    # plt.show()

    common_regions_df = merged.groupby([Header.ref_id + "_e", Header.ref_id + '_a'])[
        Header.region].count().reset_index()
    if common_regions_df.size != 0:
        common_regions_df = common_regions_df.merge(
                merged[[Header.ref_id + "_e", Header.ref_id + '_a', Header.weight + '_e', Header.weight + '_a']],
                on=[Header.ref_id + "_e", Header.ref_id + '_a'],
                how="inner")
        common_regions_df['is_contains'] = (common_regions_df[Header.weight + '_a'] == common_regions_df[Header.region])
        common_regions_df['is_equals'] = (common_regions_df[Header.weight + '_e'] == common_regions_df[Header.region])
    else:
        # print ("\n\n the original has no match!")
        common_regions_df = merged.copy()
        common_regions_df['is_contains'] = False
        common_regions_df['is_equals'] = False

    score_table_eqaul = pd.pivot_table(common_regions_df,
                                       values='is_equals',
                                       index=[Header.ref_id + '_e'],
                                       columns=[Header.ref_id + '_a'],
                                       aggfunc=np.sum,
                                       fill_value=0)
    score_table_contain = pd.pivot_table(common_regions_df,
                                         values='is_contains',
                                         index=[Header.ref_id + '_e'],
                                         columns=[Header.ref_id + '_a'],
                                         aggfunc=np.sum,
                                         fill_value=0)


    actual_priors_index = common_regions_df.drop_duplicates(Header.ref_id + '_a')[Header.ref_id + '_a']
    filtered_actual = actual_df[actual_df[Header.ref_id].isin(actual_priors_index)].drop_duplicates(Header.ref_id)
    actual_priors = filtered_actual[Header.prior]
    actual_priors.index= filtered_actual[Header.ref_id]

    expected_priors_index = common_regions_df.drop_duplicates(Header.ref_id + '_e')[Header.ref_id + '_e']
    filtered_expected = expected_df[expected_df[Header.ref_id].isin(expected_priors_index)].drop_duplicates(
        Header.ref_id)
    expected_priors = filtered_expected[Header.prior]
    expected_priors.index = filtered_expected[Header.ref_id]
    #
    # actual_priors = actual_df.drop_duplicates(Header.ref_id)[Header.prior]
    # actual_priors.index = actual_df.drop_duplicates(Header.ref_id)[Header.ref_id]
    # expected_priors = expected_df.drop_duplicates(Header.ref_id)[Header.prior]
    # expected_priors.index = expected_df.drop_duplicates(Header.ref_id)[Header.ref_id]

    return score_table_contain, score_table_eqaul, actual_priors, expected_priors




def calc_recall(score_table, expected_priors):
    """
    Recall  = sum(I_Ei * f_Ei)
    :param score_table:
    :param expected_priors:
    :return: recall
    """
    # print score_table.head()
    expected_indicator = (score_table > 0).any(axis=1)
    recall = sum(expected_indicator.mul(expected_priors))
    # logging.info( "recall = {}".format(recall))
    # logging.info(score_table[expected_indicator == False].index)
    return recall


def calc_precision(score_table, actual_priors):
    """
    Precision = sum(I_Ai * f_Ai)
    :param score_table:
    :param actual_priors:
    :return: precision
    """
    actual_indicator = (score_table > 0).any()
    precision = sum(actual_indicator.mul(actual_priors))
    # logging.info( "precision = {}".format(precision))
    return precision

    # logging.info( score_table.transpose()[actual_indicator == False].index)


def get_cmd_arguments(argv = sys.argv[1:]):
    USAGE = \
        """usage: %prog WORKING_DIR [required_options] [options]

        Compare the expected results to the actual results of emirge_smurf
        """

    parser = OptionParser(USAGE)

    # group_opt = OptionGroup(parser, "Optional flags",
    #                          "These flags are all optional.")
    #
    # group_opt.add_option("-r", dest="reference_id",
    #                       type="string", default=None,
    #                       help="reference id to compare")
    #
    # parser.add_option_group(group_opt)

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
    def __init__(self, workdir, index, changed_bacterias, changed_bases, mix_size=500000, gary_path=None):
        self.changed_bacterias=changed_bacterias
        self.changed_bases=changed_bases
        self.precisions_similarity=[]
        self.recall_similarity = []
        self.precisions_contains = []
        self.recall_contains=[]
        self.recall_contains_changed = []
        self.recall_similarity_changed = []
        self.recall_contains_unchanged = []
        self.recall_similarity_unchanged = []
        self.precision_contains_changed = []
        self.precision_similarity_changed = []
        self.precision_contains_unchanged = []
        self.precision_similarity_unchanged = []
        self.mix_size = mix_size
        self.mix_size_str = str(mix_size/1000) + 'k'
        self.index= index
        self.results_path = os.path.join(workdir, str(index), "python_results")
        self.actual_path = self.results_path
        self.is_smurf2=True
        if gary_path is not None:
            self.actual_path = os.path.join(gary_path, str(index))
            self.is_smurf2 = False

    def add_data(self, precisions_similarity, recall_similarity, precisions_contains, recall_contains, df_cols):
        self.precisions_similarity.append(precisions_similarity)
        self.precisions_contains.append(precisions_contains)
        self.recall_contains.append(recall_contains)
        self.recall_similarity.append(recall_similarity)
        return self.get_df(df_cols, 'all', precisions_similarity, precisions_contains, recall_similarity, recall_contains)

    def add_bact_change_data(self, precision_similarity, precision_contain, recall_similarity, recall_contain, df_cols):
        self.precision_contains_changed.append(precision_contain)
        self.precision_similarity_changed.append(precision_similarity)
        self.recall_contains_changed.append(recall_contain)
        self.recall_similarity_changed.append(recall_similarity)
        return self.get_df(df_cols, 'changed', precision_similarity, precision_contain, recall_similarity, recall_contain)

    def add_bact_change_source_data(self, precision_similarity, precision_contain, recall_similarity, recall_contain, df_cols):
        # print('add mutated-source!!!! - index = {}, precision = {}, recall = {}'.format(self.index, precision_contain, recall_contain))
        return self.get_df(df_cols, 'unchanged-origin', precision_similarity, precision_contain, recall_similarity,
                           recall_contain)

    def add_bact_unchange_data(self, precision_similarity, precision_contain, recall_similarity, recall_contain, df_cols):
        self.precision_contains_unchanged.append(precision_contain)
        self.precision_similarity_unchanged.append(precision_similarity)
        self.recall_contains_unchanged.append(recall_contain)
        self.recall_similarity_unchanged.append(recall_similarity)
        return self.get_df(df_cols, 'unchanged', precision_similarity, precision_contain, recall_similarity, recall_contain)

    def get_df(self, cols, bacteira_type, precisions_similarity, recall_similarity, precisions_contains, recall_contains):
        df_data = [[self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, SIMILAR, RECALL, recall_similarity,
                    bacteira_type, self.is_smurf2],
                   [self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, CONTAIN, RECALL, recall_contains,
                    bacteira_type, self.is_smurf2],
                   [self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, SIMILAR, PRECISION, precisions_similarity,
                    bacteira_type, self.is_smurf2],
                   [self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, CONTAIN, PRECISION, precisions_contains,
                    bacteira_type, self.is_smurf2]]
        return pd.DataFrame(df_data, columns=cols)


    @staticmethod
    def get_exp_key(changed_bacterias, changed_bases, is_smurf2):
        return "bctr{}bases{}is_smurf2{}".format(changed_bacterias, changed_bases, is_smurf2)


    def get_my_exp_key(self):
        return SingleExperimentData.get_exp_key(self.changed_bacterias, self.changed_bases, self.is_smurf2)


    def is_changed_bacterias(self, changed_bacterias):
        if self.changed_bacterias == changed_bacterias:
            return True
        return False

    def is_changed_bases(self, changed_bases):
        if self.changed_bases == changed_bases:
            return True
        return False


def get_stats_precision(experiments, ids, is_similarity):
    if is_similarity:
        precision_mean = [np.mean(experiments[id].precisions_similarity) for id in ids]
        precision_std = [np.std(experiments[id].precisions_similarity) for id in ids]
    else:
        precision_mean = [np.mean(experiments[id].precisions_contains) for id in ids]
        precision_std = [np.std(experiments[id].precisions_contains) for id in ids]
    return precision_mean, precision_std


def get_stats_recall(experiments, ids, is_similarity):
    if is_similarity:
        recall_mean = [np.mean(experiments[id].recall_similarity) for id in ids]
        recall_std = [np.std(experiments[id].recall_similarity) for id in ids]
    else:
        recall_mean = [np.mean(experiments[id].recall_contains) for id in ids]
        recall_std = [np.std(experiments[id].recall_contains) for id in ids]

    return recall_mean, recall_std


def get_stats_recall_changed_bact(experiments, ids, is_similarity):
    if is_similarity:
        recall_mean = [np.mean(experiments[id].recall_similarity_changed) for id in ids]
        recall_std = [np.std(experiments[id].recall_similarity_changed) for id in ids]
    else:
        recall_mean = [np.mean(experiments[id].recall_contains_changed) for id in ids]
        recall_std = [np.std(experiments[id].recall_contains_changed) for id in ids]

    return recall_mean, recall_std


def get_stats_recall_unchanged_bact(experiments, ids, is_similarity):
    if is_similarity:
        recall_mean = [np.mean(experiments[id].recall_similarity_unchanged) for id in ids]
        recall_std = [np.std(experiments[id].recall_similarity_unchanged) for id in ids]
    else:
        recall_mean = [np.mean(experiments[id].recall_contains_unchanged) for id in ids]
        recall_std = [np.std(experiments[id].recall_contains_unchanged) for id in ids]

    return recall_mean, recall_std


def get_stats_precision_changed_bact(experiments, ids, is_similarity):
    if is_similarity:
        precision_mean = [np.mean(experiments[id].precision_similarity_changed) for id in ids]
        precision_std = [np.std(experiments[id].precision_similarity_changed) for id in ids]
    else:
        precision_mean = [np.mean(experiments[id].precision_contains_changed) for id in ids]
        precision_std = [np.std(experiments[id].precision_contains_changed) for id in ids]

    return precision_mean, precision_std


def get_stats_precision_unchanged_bact(experiments, ids, is_similarity):
    if is_similarity:
        precision_mean = [np.mean(experiments[id].precision_similarity_unchanged) for id in ids]
        precision_std = [np.std(experiments[id].precision_similarity_unchanged) for id in ids]
    else:
        precision_mean = [np.mean(experiments[id].precision_contains_unchanged) for id in ids]
        precision_std = [np.std(experiments[id].precision_contains_unchanged) for id in ids]

    return precision_mean, precision_std


def create_figure3_old(experiments, working_dir, use_similarity):
    # Figure #3: SMIRF performance
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9, 7))

    # 3-1
    relevant_experiments = [4, 8, 9]
    precision_mean, precision_std = get_stats_precision(experiments, relevant_experiments, use_similarity)
    recall_mean, recall_std = get_stats_recall(experiments, relevant_experiments, use_similarity)

    axes[0].set_title('Number of reads')
    axes[0].set_xlabel("# Reads")
    axes[0].set_ylabel("Evaluation")
    axes[0].errorbar(range(len(relevant_experiments)),
                     precision_mean, yerr=precision_std,
                     fmt='o', color='g', label='precision')
    axes[0].errorbar(range(len(relevant_experiments)), recall_mean, yerr=recall_std,
                 fmt='o', color='b', label='recall')

    axes[0].set_xticks(range(len(relevant_experiments)))
    axes[0].set_xticklabels([str(experiments[id].mix_size / 1000) + "k" for id in relevant_experiments])

    # 3-2
    relevant_experiments = [3, 4, 5, 6, 7]
    precision_mean, precision_std = get_stats_precision(experiments, relevant_experiments, use_similarity)
    recall_mean, recall_std = get_stats_recall(experiments, relevant_experiments, use_similarity)

    axes[1].set_title('Number of changes per bacteria')
    axes[1].set_xlabel("# Changes per bacteria")
    axes[1].set_ylabel("Evaluation")
    axes[1].errorbar(range(len(relevant_experiments)),
                     precision_mean, yerr=precision_std,
                     fmt='o', color='g', label='precision')
    axes[1].errorbar(range(len(relevant_experiments)), recall_mean, yerr=recall_std,
                 fmt='o', color='b', label='recall')

    axes[1].set_xticks(range(len(relevant_experiments)))
    axes[1].set_xticklabels([str(experiments[id].changed_bases) for id in relevant_experiments])

    # 3-3
    relevant_experiments = [1, 4]
    precision_mean, precision_std = get_stats_precision(experiments, relevant_experiments, use_similarity)
    recall_mean, recall_std = get_stats_recall(experiments, relevant_experiments, use_similarity)

    axes[2].set_title('Number of bacteria changed')
    axes[2].set_xlabel("# Bacteria changed")
    axes[2].set_ylabel("Evaluation")
    axes[2].errorbar(range(len(relevant_experiments)),
                     precision_mean, yerr=precision_std,
                     fmt='o', color='g', label='precision')
    axes[2].errorbar(range(len(relevant_experiments)), recall_mean,  yerr=recall_std,
                     fmt='o', color='b', label='recall')

    axes[2].set_xticks(range(len(relevant_experiments)))
    axes[2].set_xticklabels([str(experiments[id].changed_bacterias) for id in relevant_experiments])

    legend = axes[0].legend()

    plt.tight_layout()
    if use_similarity:
        plt.savefig(os.path.join(working_dir, "Figure3_{}.png".format("similarity")), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(working_dir, "Figure3_{}.png".format("contain")), bbox_inches='tight')


def create_figure3_changedVsunchanged(experiments, working_dir, use_similarity):
    # Figure #3: SMIRF performance
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(9, 9))

    # 3-1
    relevant_experiments = [4, 8, 9]
    recall_mean_changed, recall_std_changed = get_stats_recall_changed_bact(experiments, relevant_experiments, use_similarity)
    recall_mean_unchanged, recall_std_unchanged = get_stats_recall_unchanged_bact(experiments, relevant_experiments,
                                                                       use_similarity)

    axes[0, 0].set_title('Number of reads')
    axes[0, 0].set_xlabel("# Reads")
    axes[0, 0].set_ylabel("Recall")
    axes[0, 0].errorbar(range(len(relevant_experiments)),
                     recall_mean_changed, yerr=recall_std_changed,
                     fmt='o', color='g', label='changed bacterias')
    axes[0, 0].errorbar(range(len(relevant_experiments)), recall_mean_unchanged, yerr=recall_std_unchanged,
                 fmt='o', color='b', label='unchanged bacterias')

    axes[0, 0].set_xticks(range(len(relevant_experiments)))
    axes[0, 0].set_xticklabels([str(experiments[id].mix_size / 1000) + "k" for id in relevant_experiments])

    # 3-2
    relevant_experiments = [3, 4, 5, 6, 7]
    recall_mean_changed, recall_std_changed = get_stats_recall_changed_bact(experiments, relevant_experiments,
                                                                            use_similarity)
    recall_mean_unchanged, recall_std_unchanged = get_stats_recall_unchanged_bact(experiments, relevant_experiments,
                                                                                  use_similarity)

    axes[0,1].set_title('Changes per bacteria')
    axes[0,1].set_xlabel("# Changes per bacteria")
    axes[0,1].set_ylabel("Recall")
    axes[0,1].errorbar(range(len(relevant_experiments)),
                     recall_mean_changed, yerr=recall_std_changed,
                     fmt='o', color='g', label='changed bacterias')
    axes[0,1].errorbar(range(len(relevant_experiments)), recall_mean_unchanged, yerr=recall_std_unchanged,
                 fmt='o', color='b', label='unchanged bacterias')

    axes[0,1].set_xticks(range(len(relevant_experiments)))
    axes[0,1].set_xticklabels([str(experiments[id].changed_bases) for id in relevant_experiments])

    # 3-3
    relevant_experiments = [1, 4]
    recall_mean_changed, recall_std_changed = get_stats_recall_changed_bact(experiments, relevant_experiments,
                                                                            use_similarity)
    recall_mean_unchanged, recall_std_unchanged = get_stats_recall_unchanged_bact(experiments, relevant_experiments,
                                                                                  use_similarity)

    axes[0,2].set_title('Bacteria changed')
    axes[0,2].set_xlabel("# Bacteria changed")
    axes[0,2].set_ylabel("Recall")
    axes[0,2].errorbar(range(len(relevant_experiments)),
                     recall_mean_changed, yerr=recall_std_changed,
                     fmt='o', color='g', label='changed bacterias')
    axes[0,2].errorbar(range(len(relevant_experiments)), recall_mean_unchanged, yerr=recall_std_unchanged,
                 fmt='o', color='b', label='unchanged bacterias')

    axes[0,2].set_xticks(range(len(relevant_experiments)))
    axes[0,2].set_xticklabels([str(experiments[id].changed_bacterias) for id in relevant_experiments])

    # 3-1
    relevant_experiments = [4, 8, 9]
    Precision_mean_changed, Precision_std_changed = get_stats_precision_changed_bact(experiments, relevant_experiments,
                                                                            use_similarity)
    Precision_mean_unchanged, Precision_std_unchanged = get_stats_precision_unchanged_bact(experiments, relevant_experiments,
                                                                                  use_similarity)

    axes[1, 0].set_title('Number of reads')
    axes[1, 0].set_xlabel("# Reads")
    axes[1, 0].set_ylabel("Precision")
    axes[1, 0].errorbar(range(len(relevant_experiments)),
                        Precision_mean_changed, yerr=Precision_std_changed,
                        fmt='o', color='g', label='changed bacterias')
    axes[1, 0].errorbar(range(len(relevant_experiments)), Precision_mean_unchanged, yerr=Precision_std_unchanged,
                        fmt='o', color='b', label='unchanged bacterias')

    axes[1, 0].set_xticks(range(len(relevant_experiments)))
    axes[1, 0].set_xticklabels([str(experiments[id].mix_size / 1000) + "k" for id in relevant_experiments])

    # 3-2
    relevant_experiments = [3, 4, 5, 6, 7]
    Precision_mean_changed, Precision_std_changed = get_stats_precision_changed_bact(experiments, relevant_experiments,
                                                                            use_similarity)
    Precision_mean_unchanged, Precision_std_unchanged = get_stats_precision_unchanged_bact(experiments, relevant_experiments,
                                                                                  use_similarity)

    axes[1, 1].set_title('Changes per bacteria')
    axes[1, 1].set_xlabel("# Changes per bacteria")
    axes[1, 1].set_ylabel("Precision")
    axes[1, 1].errorbar(range(len(relevant_experiments)),
                        Precision_mean_changed, yerr=Precision_std_changed,
                        fmt='o', color='g', label='changed bacterias')
    axes[1, 1].errorbar(range(len(relevant_experiments)), Precision_mean_unchanged, yerr=Precision_std_unchanged,
                        fmt='o', color='b', label='unchanged bacterias')

    axes[1, 1].set_xticks(range(len(relevant_experiments)))
    axes[1, 1].set_xticklabels([str(experiments[id].changed_bases) for id in relevant_experiments])

    # 3-3
    relevant_experiments = [1, 4]
    Precision_mean_changed, Precision_std_changed = get_stats_precision_changed_bact(experiments, relevant_experiments,
                                                                            use_similarity)
    Precision_mean_unchanged, Precision_std_unchanged = get_stats_precision_unchanged_bact(experiments, relevant_experiments,
                                                                                  use_similarity)

    axes[1, 2].set_title('Bacteria changed')
    axes[1, 2].set_xlabel("# Bacteria changed")
    axes[1, 2].set_ylabel("Precision")
    axes[1, 2].errorbar(range(len(relevant_experiments)),
                        Precision_mean_changed, yerr=Precision_std_changed,
                        fmt='o', color='g', label='changed bacterias')
    axes[1, 2].errorbar(range(len(relevant_experiments)), Precision_mean_unchanged, yerr=Precision_std_unchanged,
                        fmt='o', color='b', label='unchanged bacterias')

    axes[1, 2].set_xticks(range(len(relevant_experiments)))
    axes[1, 2].set_xticklabels([str(experiments[id].changed_bacterias) for id in relevant_experiments])

    plt.tight_layout()
    if use_similarity:
        plt.savefig(os.path.join(working_dir, "Figure3_ChangedVsUnchanged_{}.png".format("similarity")), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(working_dir, "Figure3_ChangedVsUnchanged_{}.png".format("contain")), bbox_inches='tight')


def create_fig2(df, workdir, hue_order, exp_ix):
    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacteria', '# Bacterias changed', 'score-type',
                           'value-type', 'value', 'bacteria', 'is_smurf2']
    relevant_experiments = [exp_ix]
    #hue_order = [0, 1]
    df = df[df['exp-index'].isin(relevant_experiments)]
    df = df.replace({"is_smurf2": {True: 'SMURF2', False: 'SMURF'}})

    g = sns.FacetGrid(df, row="value-type", col="score-type")
    g = (g.map(sns.barplot, 'is_smurf2', 'value', 'bacteria', hue_order=hue_order, palette=sns.color_palette('colorblind'),
          order=list(np.unique(df['is_smurf2']))).add_legend())
    save_fig(g, "exp-{}".format(exp_ix), workdir)
    return True


def create_fig_smurf_vs_smurf2(df, workdir, hue_order):
    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacteria', '# Bacterias changed', 'score-type',
                           'value-type', 'value', 'bacteria', 'is_smurf2']
    relevant_experiments = [0, 1]

    df = df[df['exp-index'].isin(relevant_experiments)]
    df = df[df["score-type"] == 'contain']
    df = df.replace({"is_smurf2": {True: 'SMURF2', False: 'SMURF'},
                     'exp-index': {0: 'Intra DB mixture', 1: 'Combination Intra DB and exterior bacteria'}})

    g = sns.FacetGrid(df, row="value-type", col="exp-index", aspect=1.2)
    colors = [sns.color_palette('colorblind')[2], sns.color_palette('colorblind')[-1]]
    g = (g.map(sns.barplot, 'is_smurf2', 'value', hue_order=['SMURF','SMURF2' ], palette=colors))
    save_fig(g, "smurfVsSmurf2", workdir)

    df_stats = df.groupby(["value-type", "exp-index", 'is_smurf2'])['value'].describe()
    df_stats.to_csv(os.path.join(workdir, "statistics.csv"))

    return True


def save_fig(g, fig_name, workdir):
    for ax in g.axes.flat:
        _ = plt.setp(ax.get_yticklabels(), visible=True)
        _ = plt.setp(ax.get_xticklabels(), visible=True)

    g.set_ylabels('')
    g.set_xlabels('')
    for ax in g.axes.flat:
        plt.setp(ax.texts, text="")
    g.set_titles(template='{row_name}|{col_name}')

    for i, axes_row in enumerate(g.axes):
        for j, axes_col in enumerate(axes_row):
            row, col = axes_col.get_title().split('|')
            if i == 0:
                axes_col.set_title(col.strip())
            else:
                axes_col.set_title('')
            # axes_col.set_title('')

            if j == 0:
                axes_col.set_ylabel(row.strip())
    g.savefig(os.path.join(workdir, "Figure2_{}.png".format(fig_name)))
    g.savefig(os.path.join(workdir, "Figure2_{}.pdf".format(fig_name)))

def main():

    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    working_dir = os.path.abspath(args[0])


    try:
        experiments_df = pd.read_csv(os.path.join(working_dir, "experiments_df.csv"), index_col=None)

        hue_order = ['all']  # , 'unchanged-origin']

        create_fig_smurf_vs_smurf2(experiments_df, working_dir, hue_order)
        create_fig2(experiments_df, working_dir, hue_order, 0)
        create_fig2(experiments_df, working_dir, hue_order, 1)

    except Exception as ex:
        print(ex)
        #
        # os.chdir(working_dir)

        experiments={}
        experiments[0] = SingleExperimentData(working_dir, 0, 0, 0, 500000,  "/home/vered/EMIRGE/cluster/cluster_smurf_gary")
        experiments[1] = SingleExperimentData(working_dir, 1, 10, 20, 500000,  "/home/vered/EMIRGE/cluster/cluster_smurf_gary")
        experiments[2] = SingleExperimentData(working_dir, 0, 0, 0, 500000)
        experiments[3] = SingleExperimentData(working_dir, 1, 10, 20, 500000)

        experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacteria', '# Bacterias changed', 'score-type',
                               'value-type', 'value', 'bacteria', 'is_smurf2']
        experiments_df = pd.DataFrame(columns=experiments_df_cols)

        for experiment in experiments.values():
            results_dirs = os.listdir(experiment.results_path)
            for test_dir in results_dirs:
                try:
                    expected_dir_path = os.path.join(experiment.results_path, test_dir)
                    actual_dir_path = os.path.join(experiment.actual_path, test_dir)
                    actual_path = os.path.join(actual_dir_path, ACTUAL_RES_FILE_NAME)
                    expected_path = os.path.join(expected_dir_path, EXPECTED_RES_FILE_NAME)
                    # print (actual_path)
                    # print (expected_path)
                    if( not os.path.isfile(expected_path)) or (not os.path.isfile(actual_path)):
                        continue

                    contain_table, similarity_table, actual_priors, expected_priors = calc_full_table(expected_path, actual_path)
                    df = experiment.add_data(calc_precision(similarity_table, actual_priors),
                                        calc_recall(similarity_table, expected_priors),
                                        calc_precision(contain_table, actual_priors),
                                        calc_recall(contain_table, expected_priors), experiments_df_cols)

                    experiments_df = experiments_df.append(df)

                except Exception as ex:
                    print("Failed add data for test = {}, ex = {}".format(test_dir, ex))


        experiments_df.to_csv(os.path.join(working_dir, "experiments_df.csv"), index=False)







if __name__ == "__main__":
    main()



