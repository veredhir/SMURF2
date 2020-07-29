import scipy.io as sio
import pandas as pd
from optparse import OptionParser, OptionGroup
import sys
import pandas as pd

def read_db_ind_map(path_to_db_ind):
    db_map = pd.read_csv(path_to_db_ind, index_col=None, header=None)
    db_map['fasta_id'] = db_map[0].apply(lambda x: int(x))
    db_dict = {}
    db_dict_inv = {}
    for id, fasta_id in zip(range(0, len(db_map)), db_map['fasta_id'].tolist()):
        db_dict[fasta_id] = id
        db_dict_inv[id] = fasta_id
    return db_dict, db_dict_inv


if __name__ == "__main__":
    """
    command line interface to emirge

    """
    parser = OptionParser("Convert SMURF2 results format to SMURF results format")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "")
    group_reqd.add_option("-1", "--smurf_results_mat",
                          type="string",
                          help="path to smurf results mat")
    group_reqd.add_option("-2", "--smurf2_results_mat",
                          type="string",
                          help="path to smurf results mat after V2G")
    group_reqd.add_option("-r", "--results_path",
                          type="string", default="comparison.csv",
                          help="path to smurf results mat after V2G")
    group_reqd.add_option("-t", "--taxa_path",
                          type="string",
                          help="taxa_and_head file path")

    parser.add_option_group(group_reqd)
    argv = sys.argv[1:]
    (options, args) = parser.parse_args(argv)

    db_map, db_map_inv = read_db_ind_map(options.taxa_path)  # "/home/vered/EMIRGE/data/reference_db/Header_uni_forVered.csv"

    smurf_s26_output = options.smurf_results_mat  # '/home/vered/EMIRGE/data/s26_mock/SUMRF_samples26.mat'

    smurf_s26 = sio.loadmat(smurf_s26_output)

    smurf_group_ids = [i[0][0] for i in smurf_s26['z_gr']]
    smurf_group_freq = [i[0] for i in smurf_s26['z_freq']]

    smurf2_s26_output = options.smurf2_results_mat  # '/home/vered/EMIRGE/data/s26_mock/sample_s26_mock_results_smurf.mat'

    smurf2_s26 = sio.loadmat(smurf2_s26_output)

    smurf2_group_ids = smurf2_s26['bactMetaGroups']['db_ind'][0]
    smurf2_group_freq = smurf2_s26['found_bacteria']['frequency'][0][0][0]
    print smurf2_group_freq

    matches = []
    matches2 = []
    matches3 = []
    data = []
    for smurf_gr, index, smurf_freq in zip(smurf_group_ids, range(len(smurf_group_ids)), smurf_group_freq):
        match_group = False
        match_group2 = False
        match_group3 = False
        match_ids = []
        smurf2_freq = 0
        # print "*******\n\n"
        print "Test SMURF ID = {}: freq = {}".format(index, smurf_freq)
        for id1 in smurf_gr:
            for smuf2_gr, curr_smurf2_freq, smurf2_index in zip(smurf2_group_ids, smurf2_group_freq, range(len(smurf2_group_freq))):
                for id2 in smuf2_gr[0]:
                    if db_map_inv[id2] == id1:
                        match_group = True
                        match_ids.append(id2)
                        smurf2_freq += curr_smurf2_freq
                        print curr_smurf2_freq
                        data.append({'smurf index': index,
                                     'smurf2_index': smurf2_index,
                                     'smurf2_id': id2,
                                     'smurf_id': id1,
                                     'smurf_groups': smurf_gr,
                                     'smurf_freq': smurf_freq,
                                     'smurf2_groups': match_ids,
                                     'smurf2_freq': curr_smurf2_freq})
                    if id2 == id1:
                        match_group3 = True
                        break
        matches.append(match_group)
        matches2.append(match_group2)
        matches3.append(match_group3)
        if match_group == False:
            data.append({'smurf index': index,
                         'smurf2_index': -1,
                         'smurf2_id': -1,
                         'smurf_id': -1,
                         'smurf_groups': smurf_gr,
                         'smurf_freq': smurf_freq,
                         'smurf2_groups': -1,
                         'smurf2_freq': 0})

    print "matches = {}/{}".format(sum(matches), len(matches))
    print "matches2 = {}/{}".format(sum(matches2), len(matches))
    print "matches3 = {}/{}".format(sum(matches3), len(matches))
    data_df = pd.DataFrame.from_dict(data)
    data_df.to_csv(options.results_path)

    smurf2_res = []
    for smuf2_gr, curr_smurf2_freq, smurf2_index in zip(smurf2_group_ids, smurf2_group_freq, range(len(smurf2_group_freq))):
        smurf2_res.append({'smurf2_index': smurf2_index,
                           'smurf2_freq': curr_smurf2_freq})

    pd.DataFrame.from_dict(smurf2_res).to_csv("s2.csv")

    smurf_res = []
    for smurf_gr, index, smurf_freq in zip(smurf_group_ids, range(len(smurf_group_ids)), smurf_group_freq):
        smurf2_res.append({'smurf_index': index,
                           'smurf_freq': smurf_freq})

    pd.DataFrame.from_dict(smurf_res).to_csv("s1.csv")

    print matches