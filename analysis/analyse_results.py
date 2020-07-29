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


def calc_table(expected_df, actual_df):
    merged = expected_df.merge(actual_df, on=[Header.region, Header.sequence], how="inner", suffixes=('_e', '_a'))

    common_regions_df = merged.groupby([Header.ref_id + "_e", Header.ref_id + '_a'])[
        Header.region].count().reset_index()
    if common_regions_df.size != 0:
        common_regions_df = common_regions_df.merge(
                merged[[Header.ref_id + "_e", Header.ref_id + '_a', Header.weight + '_e', Header.weight + '_a']],
                on=[Header.ref_id + "_e", Header.ref_id + '_a'],
                how="inner")
        common_regions_df['is_contains'] = (common_regions_df[Header.weight + '_a'] == common_regions_df[Header.region])
    else:
        # print ("\n\n the original has no match!")
        common_regions_df = merged.copy()
        common_regions_df['is_contains'] = False

    return common_regions_df

#/home/vered/EMIRGE/cluster/cluster_default/test_0


def get_cmd_arguments(argv = sys.argv[1:]):
    USAGE = \
        """usage: %prog WORKING_DIR [required_options] [options]

        Compare the expected results to the actual results of emirge_smurf
        """

    parser = OptionParser(USAGE)
    group_opt = OptionGroup(parser, "Optional parameters",
                            "Configurable parameters")
    group_opt.add_option("-r", "--results_path",
                         type="string", default='/tmp/smurf2_analyse_results.csv',
                         help="""Results file path (default=%default)""")
    parser.add_option_group(group_opt)
    (options, args) = parser.parse_args(argv)

    if len(args) != 1:
        parser.error(
            "WORKING_DIR is required, and all options except DIR should have a flag associated with them (options without flags: %s)" % args)
        return

    return options, args


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



def get_single_test_data(dir_path):
    expected = os.path.join(dir_path, EXPECTED_RES_FILE_NAME)
    actual = os.path.join(dir_path, ACTUAL_RES_FILE_NAME)
    expected_df = get_expected_df(expected)
    actual_df = get_test_df(actual)
    contains =  calc_table(expected_df, actual_df)
    contains.to_csv(os.path.join(dir_path, 'contains.csv'))
    return contains

def main():
    print 'test evaluation'
    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    working_dir = os.path.abspath(args[0])

    dfs = []
    for i in range(1):
        df = get_single_test_data(working_dir)
        df['test_name'] = i
        dfs.append(df)

    all_data = pd.concat(dfs)
    all_data.to_csv(options.results_path)




if __name__ == "__main__":
    main()