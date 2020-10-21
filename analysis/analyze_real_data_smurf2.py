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

def main():
    print 'test evaluation'
    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()





if __name__ == "__main__":
    main()