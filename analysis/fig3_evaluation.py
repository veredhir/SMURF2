
import numpy as np
from smurf2_headers import *
from smurf2_utills import *
import os
import seaborn as sns

# import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
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
    df = pd.DataFrame.from_csv(path, index_col=None)
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
    actual_df = get_test_df(actual_path)
    return calc_table(expected_df, actual_df)

def calc_changed_bacteria_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    expected_df = expected_df[expected_df[Header.is_changed] == True]
    actual_df = get_test_df(actual_path)
    score_table_contain, score_table_eqaul, actual_priors, expected_priors =  calc_table(expected_df, actual_df)
    sum_priors = sum(actual_priors)
    actual_priors = [prior/sum_priors for prior in actual_priors]
    sum_priors = expected_df.drop_duplicates(Header.ref_id)[Header.prior].sum()
    expected_priors = [prior / sum_priors for prior in expected_priors]
    return score_table_contain, score_table_eqaul, actual_priors, expected_priors

def calc_unchagned_bacteria_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    expected_df = expected_df[(expected_df[Header.is_changed] == False) & (expected_df[Header.is_changed_source] == False)]
    actual_df = get_test_df(actual_path)
    score_table_contain, score_table_eqaul, actual_priors, expected_priors = calc_table(expected_df, actual_df)
    sum_priors = sum(actual_priors)
    actual_priors = [prior / sum_priors for prior in actual_priors]
    sum_priors = expected_df.drop_duplicates(Header.ref_id)[Header.prior].sum()
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
    sum_priors = expected_df.drop_duplicates(Header.ref_id)[Header.prior].sum()
    expected_priors = [prior / sum_priors for prior in expected_priors]
    return score_table_contain, score_table_eqaul, actual_priors, expected_priors


def calc_table(expected_df, actual_df):
    merged = expected_df.merge(actual_df, on=[Header.region, Header.sequence], how="inner", suffixes=('_e', '_a'))

    common_regions_df = merged.groupby([Header.ref_id + "_e", Header.ref_id + '_a'])[
        Header.region].count().reset_index()
    if common_regions_df.size != 0:
        common_regions_df = common_regions_df.merge(
                merged[[Header.ref_id + "_e", Header.ref_id + '_a', Header.weight + '_e', Header.weight + '_a']],
                on=[Header.ref_id + "_e", Header.ref_id + '_a'],
                how="inner")
        common_regions_df['is_contains'] = True #(common_regions_df[Header.weight + '_a'] >= common_regions_df[Header.region])
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
        return self.get_df(df_cols, MODIFIED, precision_similarity, precision_contain, recall_similarity, recall_contain)

    def add_bact_change_source_data(self, precision_similarity, precision_contain, recall_similarity, recall_contain, df_cols):
        # print('add mutated-source!!!! - index = {}, precision = {}, recall = {}'.format(self.index, precision_contain, recall_contain))
        return self.get_df(df_cols, ORIGIN, precision_similarity, precision_contain, recall_similarity,
                           recall_contain)

    def add_bact_unchange_data(self, precision_similarity, precision_contain, recall_similarity, recall_contain, df_cols):
        self.precision_contains_unchanged.append(precision_contain)
        self.precision_similarity_unchanged.append(precision_similarity)
        self.recall_contains_unchanged.append(recall_contain)
        self.recall_similarity_unchanged.append(recall_similarity)
        return self.get_df(df_cols, FROM_DB, precision_similarity, precision_contain, recall_similarity, recall_contain)

    def get_df(self, cols, bacteira_type, precisions_similarity, recall_similarity, precisions_contains, recall_contains):
        df_data = [[self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, SIMILAR, RECALL, recall_similarity,
                    bacteira_type],
                   [self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, CONTAIN, RECALL, recall_contains,
                    bacteira_type],
                   [self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, SIMILAR, PRECISION, precisions_similarity,
                    bacteira_type],
                   [self.index, self.mix_size_str, self.changed_bases, self.changed_bacterias, CONTAIN, PRECISION, precisions_contains,
                    bacteira_type]]
        return pd.DataFrame(df_data, columns=cols)


    @staticmethod
    def get_exp_key(changed_bacterias, changed_bases):
        return "bctr{}bases{}".format(changed_bacterias, changed_bases)


    def get_my_exp_key(self):
        return SingleExperimentData.get_exp_key(self.changed_bacterias, self.changed_bases)


    def is_changed_bacterias(self, changed_bacterias):
        if self.changed_bacterias == changed_bacterias:
            return True
        return False

    def is_changed_bases(self, changed_bases):
        if self.changed_bases == changed_bases:
            return True
        return False


def create_fig3_bases(df, workdir, hue_order):
    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacterium', '# Bacteria changed', 'score-type',
                           'value-type', 'value', 'bacteria']
    relevant_experiments = [3, 4, 5, 6, 7]
    df = df[df['exp-index'].isin(relevant_experiments)]

    g = sns.FacetGrid(df, row="value-type", col="score-type")
    g = (g.map(sns.barplot, '# Changes per bacterium', 'value', 'bacteria', palette=sns.color_palette('colorblind'),
          hue_order=hue_order, order=list(np.unique(df['# Changes per bacterium']))))
    g.axes[1][0].legend(bbox_to_anchor=(1.2, 1.1))
    save_fig(g, "bases", workdir)

    return True

def create_fig3_reads(df, workdir, hue_order):
    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacterium', '# Bacteria changed', 'score-type',
                           'value-type', 'value', 'bacteria']
    relevant_experiments = [4, 8, 9]
    df = df[df['exp-index'].isin(relevant_experiments)]

    g = sns.FacetGrid(df, row="value-type", col="score-type")
    g = (g.map(sns.barplot, '# Reads', 'value', 'bacteria', palette=sns.color_palette('colorblind'),
          hue_order=hue_order, order=list(np.unique(df['# Reads']))))
    g.axes[1][0].legend(bbox_to_anchor=(1.2, 1.1))
    save_fig(g, "reads", workdir)

    return True


def create_fig3_bacterias(df, workdir, hue_order):
    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacterium', '# Bacteria changed', 'score-type',
                           'value-type', 'value', 'bacteria']
    relevant_experiments = [1, 4]
    df = df[df['exp-index'].isin(relevant_experiments)]

    g = sns.FacetGrid(df, row="value-type", col="score-type")
    g = (g.map(sns.barplot, '# Bacteria changed', 'value', 'bacteria', palette=sns.color_palette('colorblind'),
          hue_order=hue_order, order=list(np.unique(df['# Bacteria changed']))))
    g.axes[1][0].legend(bbox_to_anchor=(1.2, 1.1))
    save_fig(g, "bact", workdir)

    return True

def create_fig_3_unify(df, workdir, hue_order, type):
    df = df[(df['value-type'] == type) & (df['score-type'] == 'All subsets')]
    df = df.rename(columns={'value' : type})
    bases_experiments = [3, 4, 5, 6, 7]
    df_bases = df[df['exp-index'].isin(bases_experiments) ]
    reads_experiments = [4, 8, 9]
    df_reads = df[df['exp-index'].isin(reads_experiments)]
    bact_experiments = [1, 4, 10, 11]
    df_bact = df[df['exp-index'].isin(bact_experiments)]

    f, axes = plt.subplots(1, 3,figsize=(15,7.5))
    for ax in axes:
        ax.set(ylim=(0, 1.1))

    g = sns.barplot(y=type, x="# Changes per bacterium", data=df_bases, hue='bacteria', ax=axes[0], palette=sns.color_palette('colorblind'), hue_order=hue_order)
    g.legend_.remove()
    # Define some hatches
    hatch = '\\'
    # Loop over the bars
    for i, thisbar in enumerate(g.patches):
        # Set stripes on the common settings bar:
        if i in [1, 6]:
            thisbar.set_hatch(hatch)
    g = sns.barplot(y=type, x="# Reads", data=df_reads, hue='bacteria', order = ['50k', '100k', '500k'], ax=axes[1], palette=sns.color_palette('colorblind'), hue_order=hue_order)
    g.legend_.remove()
    # Loop over the bars
    for i, thisbar in enumerate(g.patches):
        # Set stripes on the common settings bar:
        if i in [2, 5]:
            thisbar.set_hatch(hatch)
    df_bact.rename(columns={"# Bacteria changed": "# Modified bacteria"}, inplace=True)
    g = sns.barplot(y=type, x="# Modified bacteria", data=df_bact, hue='bacteria', ax=axes[2], palette=sns.color_palette('colorblind'), hue_order=hue_order)
    g.legend(bbox_to_anchor=(1, 1))
    # Loop over the bars
    for i, thisbar in enumerate(g.patches):
        # Set stripes on the common settings bar:
        if i in [1, 5]:
            thisbar.set_hatch(hatch)

    save_fig_unify(f, axes, type, workdir)

    data_path=os.path.join(workdir, "statistics_{}_xx.csv".format(type))
    for col, df in [("# Changes per bacterium", df_bases), ("# Modified bacteria", df_bact), ("# Reads", df_reads)]:
        df.groupby([col, 'bacteria'])[type].describe().to_csv(data_path.replace('xx', col.replace(' ', '_')))




def save_fig(g, fig_name, workdir):
    for ax in g.axes.flat:
        _ = plt.setp(ax.get_yticklabels(), visible=True)
        _ = plt.setp(ax.get_xticklabels(), visible=True)

    g.set_ylabels('')
    for ax in g.axes.flat:
        plt.setp(ax.texts, text="")
    g.set_titles(template='SMURF2 Evaluation')

    for i, axes_row in enumerate(g.axes):
        for j, axes_col in enumerate(axes_row):
            row, col = axes_col.get_title().split('|')
            if i == 0:
                axes_col.set_title(col.strip())
            else:
                axes_col.set_title('')

            if j == 0:
                axes_col.set_ylabel(row.strip())
    g.savefig(os.path.join(workdir, "Figure3_{}.png".format(fig_name)))
    g.savefig(os.path.join(workdir, "Figure3_{}.pdf".format(fig_name)))

def save_fig_unify(g, axes, fig_name, workdir):
    for ax in axes:
        _ = plt.setp(ax.get_yticklabels(), visible=True)
        _ = plt.setp(ax.get_xticklabels(), visible=True)

    for i in range(1, len(axes)):
        axes[i].set_ylabel('')
        axes[i].set_ylabel('')

    g.savefig(os.path.join(workdir, "Figure3_{}.png".format(fig_name)))
    g.savefig(os.path.join(workdir, "Figure3_{}.pdf".format(fig_name)))

def main():
    print 'test evaluation'
    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    working_dir = os.path.abspath(args[0])


    try:
        experiments_df = pd.read_csv(os.path.join(working_dir, "experiments_df.csv"), index_col=None)
        print experiments_df.columns
        experiments_df = experiments_df.rename(columns={'# Bacterias changed': '# Bacteria changed',
                                                        '# Changes per bacteria': '# Changes per bacterium'})
        print experiments_df.columns
        experiments_df = experiments_df[experiments_df['bacteria'] != 'all']

        # old - not needed any more
        # hue_rename = {'Intra DB - origin': FROM_DB,
        #               MODIFIED: MODIFIED,
        #               'Intra DB': ORIGIN}
        #
        # score_rename = {'contain': CONTAIN,
        #                 'similar': SIMILAR}

        # experiments_df['bacteria'] = experiments_df['bacteria'].apply(lambda b : hue_rename[b])
        # experiments_df['score-type'] = experiments_df['score-type'].apply(lambda s: score_rename[s])

        if 'mutated' in working_dir:
            hue_order = [FROM_DB, MODIFIED, ORIGIN]
        else:
            hue_order = [FROM_DB, MODIFIED]

        # create_fig3_bases(experiments_df, working_dir, hue_order)
        # create_fig3_reads(experiments_df, working_dir, hue_order)
        # create_fig3_bacterias(experiments_df, working_dir, hue_order)
        create_fig_3_unify(experiments_df, working_dir, hue_order, 'Precision')
        create_fig_3_unify(experiments_df, working_dir, hue_order, 'Recall')
    except Exception as ex:
        print("\n\nTEST!!!\n\n")
        print(ex)

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
                try:
                    test_path = os.path.join(experiment.results_path, test_dir)
                    actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME)
                    expected_path = os.path.join(test_path, EXPECTED_RES_FILE_NAME)
                    if not os.path.isfile(actual_path):
                        actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME_NEW)
                        if not os.path.isfile(actual_path):
                            continue

                    contain_table, similarity_table, actual_priors, expected_priors = calc_full_table(expected_path, actual_path)
                    df = experiment.add_data(calc_precision(similarity_table, actual_priors),
                                        calc_recall(similarity_table, expected_priors),
                                        calc_precision(contain_table, actual_priors),
                                        calc_recall(contain_table, expected_priors), experiments_df_cols)

                    experiments_df = experiments_df.append(df)

                    contain_table, similarity_table, actual_priors, expected_priors = calc_changed_bacteria_table(expected_path, actual_path)
                    df = experiment.add_bact_change_data(calc_precision(similarity_table, actual_priors),
                                                    calc_precision(contain_table, actual_priors),
                                                    calc_recall(similarity_table, expected_priors),
                                                    calc_recall(contain_table, expected_priors), experiments_df_cols)
                    experiments_df = experiments_df.append(df)

                    contain_table, similarity_table, actual_priors, expected_priors = calc_unchagned_bacteria_table(expected_path, actual_path)

                    df = experiment.add_bact_unchange_data(calc_precision(similarity_table, actual_priors),
                                                      calc_precision(contain_table, actual_priors),
                                                      calc_recall(similarity_table, expected_priors),
                                                      calc_recall(contain_table, expected_priors), experiments_df_cols)
                    experiments_df = experiments_df.append(df)

                    contain_table, similarity_table, actual_priors, expected_priors = calc_chagned_source_bacteria_table(expected_path,
                                                                                                                  actual_path)
                    if contain_table is None:
                        print("mutated source - contain table is None")
                    else:
                        df = experiment.add_bact_change_source_data(calc_precision(similarity_table, actual_priors),
                                                             calc_precision(contain_table, actual_priors),
                                                             calc_recall(similarity_table, expected_priors),
                                                             calc_recall(contain_table, expected_priors), experiments_df_cols)

                        experiments_df = experiments_df.append(df)


                except Exception as ex:
                    print("Failed add data for test = {}, ex = {}".format(test_dir, ex))

        experiments_df.to_csv(os.path.join(working_dir, "experiments_df.csv"), index=False)


def test():
    working_dir=""
    experiment = SingleExperimentData(working_dir, 10, 10, 10, 500000)

    experiments_df_cols = ['exp-index', '# Reads', '# Changes per bacterium', '# Bacteria changed', 'score-type',
                           'value-type', 'value', 'bacteria']
    experiments_df = pd.DataFrame(columns=experiments_df_cols)
    # experiment.results_path

    try:
        test_path = "/home/vered/EMIRGE/data/test_4"
        actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME)
        expected_path = os.path.join(test_path, EXPECTED_RES_FILE_NAME)
        if not os.path.isfile(actual_path):
            actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME_NEW)
            if not os.path.isfile(actual_path):
                return

        contain_table, similarity_table, actual_priors, expected_priors = calc_full_table(expected_path,
                                                                                          actual_path)
        df = experiment.add_data(0,
                                 0,
                                 calc_precision(contain_table, actual_priors),
                                 calc_recall(contain_table, expected_priors), experiments_df_cols)

        experiments_df = experiments_df.append(df)

        contain_table, similarity_table, actual_priors, expected_priors = calc_changed_bacteria_table(
            expected_path, actual_path)
        df = experiment.add_bact_change_data(0,
                                             calc_precision(contain_table, actual_priors),
                                             0,
                                             calc_recall(contain_table, expected_priors), experiments_df_cols)
        experiments_df = experiments_df.append(df)

        contain_table, similarity_table, actual_priors, expected_priors = calc_unchagned_bacteria_table(
            expected_path, actual_path)

        df = experiment.add_bact_unchange_data(calc_precision(similarity_table, actual_priors),
                                               calc_precision(contain_table, actual_priors),
                                               calc_recall(similarity_table, expected_priors),
                                               calc_recall(contain_table, expected_priors), experiments_df_cols)
        experiments_df = experiments_df.append(df)

        contain_table, similarity_table, actual_priors, expected_priors = calc_chagned_source_bacteria_table(
            expected_path,
            actual_path)
        if contain_table is None:
            print("mutated source - contain table is None")
        df = experiment.add_bact_change_source_data(calc_precision(similarity_table, actual_priors),
                                                    calc_precision(contain_table, actual_priors),
                                                    calc_recall(similarity_table, expected_priors),
                                                    calc_recall(contain_table, expected_priors),
                                                    experiments_df_cols)

        experiments_df = experiments_df.append(df)


    except Exception as ex:
        print("Failed add data for test ")

    # experiments_df.to_csv(os.path.join(working_dir, "experiments_df.csv"), index=False)


if __name__ == "__main__":
    main()
    # test()



