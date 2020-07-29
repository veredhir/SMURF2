import pandas as pd
from Bio import SeqIO, pairwise2
import os
import sys
from optparse import OptionParser, OptionGroup
import numpy as np
from Bio import pairwise2

# db_full_16S_path = '/home/veredhi/EMIRGE_SMURF/GreenGenes_201305_unique_up_to_3_ambiguous_16S.fasta'
# s26_results_path = '/home/veredhi/EMIRGE_SMURF/s26_results_16s/ssrRNAs'
#
# output_path=""
# s26_smurf2_results = '/home/veredhi/results/s26_mock_test2/tmp/workdir/emirge_smurf_WFalseSTrue.csv'
# db_full_16s_regions = '/home/veredhi/EMIRGE_SMURF/ReferenceFullDB'
# primers_path = ""

# local parameters:
s26_results_path = "/home/vered/EMIRGE/data/s26_mock/ssrRNAs/"
primers_path = "/home/vered/EMIRGE/data/primers.csv"
database_gg_regions = "/home/vered/EMIRGE/data/reference_db"
MAX_AMBIGUITY = 4

# HPC parameters:
# s26_results_path = "/home/veredhi/EMIRGE_SMURF/s26_zymo/ssrRNAs"
# primers_path = "/home/veredhi/EMIRGE_SMURF/Code/primers.csv"
# database_gg_regions = "/home/veredhi/EMIRGE_SMURF/ReferenceFullDB"
# MAX_AMBIGUITY = 4

ids_path_for_validation = "/home/vered/EMIRGE/data/s26_mock/db_ids_for_zymo_s26.csv"


match_ids = ['4489602', '4489603', '4385232.01']


class Primer(object):
    """
    initialize and store the primers for each amplified region
    Require initialization using primer.csv file
        The header contain the amplified regions
        Each row in the table contain single primer used to amplify the corresponding region.
    """

    def __init__(self):
        self.all = {}
        self.histogram = {}
        self.min_size = 100
        self.similarity_score_th = 0.7

    def init(self, primers_path):
        """
        Init the primers using the primers.csv file

        primers_path: string
                      path to primers.csv file:
                      The header contain the amplified regions
                      Each row in the table contain single primer used to amplify the corresponding region.
        """
        primers_df = pd.read_csv(primers_path, index_col=False)
        for region in primers_df.columns:
            self.add_primer_for_region(int(region), primers_df[region].fillna('').tolist())

    def add_primer_for_region(self, region, primers_list):
        """
        Add list of primers for region
        Initiate the primers histogram

        region: int
                The considered region
        primers_list: list of strings
                The list of string corresponding to the region
        """
        primers_list = filter(lambda x: x != '', primers_list)
        self.all.update({region: primers_list})
        for primer in primers_list:
            self.histogram.update({primer: 0})
            if len(primer) < self.min_size:
                self.min_size = len(primer)

    def get_region_for_16s(self, seq_id, sequence, region_size=324, read_len=145):
        regions_dict = []
        MAX_DIFF = 4
        for r, primers in self.all.iteritems():
            region_exists = False
            best_score = MAX_DIFF + 1
            curr_best_region_data = None
            for primer in primers:
                primer_len = len(primer)
                for i in range(0, len(sequence)-region_size):
                    match_score = sum([c1 != c2 for c1, c2 in zip(primer, sequence[i: i + primer_len])])
                    if match_score == 0:
                        region_exists = True
                        regions_dict.append({'Region': r,
                                             'Sequence': sequence[i:i+read_len] , #+ complete_sequence(sequence[i+region_size-read_len: i+region_size]),
                                             'seq_id': seq_id})
                        break
                    elif match_score < best_score:
                        best_score = match_score
                        curr_best_region_data = {'Region': r,
                                                 'Sequence': sequence[i:i + read_len],
                                                 # + complete_sequence(sequence[i+region_size-read_len: i+region_size]),
                                                 'seq_id': seq_id}
                if region_exists:
                    break
            if best_score <= MAX_DIFF and region_exists == False:
                print("test ", best_score)
                regions_dict.append(curr_best_region_data)

        return regions_dict

    def split_16s_to_regions(self, results_path):
        full_data = []
        bact_counter = 0
        for f in os.listdir(results_path):
            for s26_record in SeqIO.parse(os.path.join(results_path, f), "fasta"):
                bact_counter += 1
                res_seq = s26_record.seq.__str__()
                res_id = s26_record.id.__str__()
                full_data = full_data + self.get_region_for_16s(res_id, res_seq)
        print bact_counter
        return pd.DataFrame.from_dict(full_data)


def complete_sequence(seq):
    comp_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return "".join([comp_dict[b] for b in seq][::-1])


def get_16s_ids_for_26s_results(full_db_path, results_path):
    s26_ids_by_bacteria = {}
    for f in os.listdir(results_path):
        bacteria = []
        for s26_record in SeqIO.parse(os.path.join(results_path, f), "fasta"):
            print "bacteria = {}:".format(s26_record.id.__str__())
            res_seq = s26_record.seq.__str__()
            best_score = MAX_AMBIGUITY
            best_id = ''
            is_identical = False
            for db_record in SeqIO.parse(full_db_path, "fasta"):
                db_seq = db_record.seq.__str__()
                if db_seq == res_seq:
                    bacteria.append(db_record.id.__str__())
                    print db_record.id.__str__()
                    is_identical = True
                    break
                else:
                    curr_score = sum(c1!=c2 for c1,c2 in zip(db_seq,res_seq))
                    if curr_score < best_score:
                        best_score = curr_score
                        best_id = db_record.id.__str__()
            if not is_identical and best_score < MAX_AMBIGUITY:
                bacteria.append(best_id)
                print "record = {}, score = {}".format(best_id, best_score)
            print "bacteria = {} matches = {}".format(s26_record.id.__str__(), len(bacteria))
        s26_ids_by_bacteria[f] = bacteria
    return s26_ids_by_bacteria


class BacteriumMatch(object):
    def __init__(self, min_score):
        self.score = min_score
        self.db_ids = {}
        self.results = []
        for i in range(1, 6):
            self.db_ids[i] = []

    def test_for_match(self, score, db_id, region):
        if score < self.score:
            self.db_ids[region].append({'id': db_id,
                                        'score': score,
                                        'region': region})

    def get_ids_for_region(self, region):
        df = pd.DataFrame.from_dict(self.db_ids[region])
        if df.empty:
            return []
        return df['id'].tolist()

    def get_db_ids(self):
        df = pd.DataFrame.from_dict(self.db_ids[1])
        if df.empty:
            return []
        df = df.append(pd.DataFrame.from_dict(self.db_ids[2]))
        df = df.append(pd.DataFrame.from_dict(self.db_ids[3]))
        df = df.append(pd.DataFrame.from_dict(self.db_ids[4]))
        df = df.append(pd.DataFrame.from_dict(self.db_ids[5]))

        df['score'] = df.groupby('id')['score'].transform('sum')
        df['count'] = df.groupby('id')['score'].transform('count')

        df = df.drop_duplicates('id')

        df = df[df['count'] == 5]
        df_min_score = df['score'].min()
        db_ids = df[df['score'] == df_min_score]['id'].tolist()

        return db_ids


def get_bacteria_ids_from_database(region_db_path, bacteria_to_extract_df, bacteria_to_remove='', max_alignment_mismatch=5, read_len=145):
    bacteria_matches = {}
    min_score = max_alignment_mismatch

    bacteria_to_test = [bact for bact in bacteria_to_extract_df['seq_id'].unique() if bacteria_to_remove in bact]
    for bacterium_id in bacteria_to_test:
        bacteria_matches[bacterium_id] = BacteriumMatch(min_score)

    baseline_region = -1
    is_first_fasta = True
    for region_fasta in os.listdir(region_db_path):
        if region_fasta.endswith('fasta'):
            region = int(region_fasta.split('.')[0][-1])
            print region
            if is_first_fasta:
                baseline_region = region
                for record in SeqIO.parse(os.path.join(region_db_path, region_fasta), "fasta"):
                    db_sequence = record.seq.__str__()[:read_len]
                    title = record.id.__str__()
                    for bacterium_id in bacteria_to_test:
                        bacterium_seq = bacteria_to_extract_df[(bacteria_to_extract_df['seq_id'] == bacterium_id) & (bacteria_to_extract_df['Region'] == region)]['Sequence'].iloc[0]
                        alignments_score = sum([c1 != c2 for c1, c2 in zip(bacterium_seq,db_sequence)])
                        bacteria_matches[bacterium_id].test_for_match(alignments_score, title, region)
                is_first_fasta = False
            else:
                record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
                for bacterium_id in bacteria_matches:
                    bacterium_seq = bacteria_to_extract_df[(bacteria_to_extract_df['seq_id'] == bacterium_id) & (bacteria_to_extract_df['Region'] == region)]['Sequence'].iloc[0]
                    for db_id in bacteria_matches[bacterium_id].get_ids_for_region(baseline_region):
                        try:
                            record = record_dict[db_id]
                            db_sequence = record.seq.__str__()[:read_len]
                            alignments_score = sum([c1 != c2 for c1, c2 in zip(bacterium_seq,db_sequence)])
                            bacteria_matches[bacterium_id].test_for_match(alignments_score, db_id, region)
                        except Exception as ex:
                            print "error with id = {}".format(db_id)

    db_id_by_bact_dict = []
    for bacterium_id in bacteria_matches:
        print bacterium_id
        for db_id in bacteria_matches[bacterium_id].get_db_ids():
            db_id_by_bact_dict.append({"bacteria_name": bacterium_id,
                                       "db_id":         db_id})


    df_db_id_by_bact = pd.DataFrame.from_dict(db_id_by_bact_dict)
    return df_db_id_by_bact


def test_1_bacteria_ids_from_database(region_db_path, bacteria_to_extract_df, test_match_ids, max_alignment_mismatch=5, read_len=145):
    bacteria_matches = {}
    min_score = max_alignment_mismatch

    bacteria_to_test = [bact for bact in bacteria_to_extract_df['seq_id'].unique() if bacteria_to_remove in bact]
    for bacterium_id in bacteria_to_test:
        bacteria_matches[bacterium_id] = BacteriumMatch(min_score)

    baseline_region = -1
    is_first_fasta = True
    for region_fasta in os.listdir(region_db_path):
        if region_fasta.endswith('fasta'):
            region = int(region_fasta.split('.')[0][-1])
            print region
            if is_first_fasta:
                baseline_region = region
                record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
                for match_id in test_match_ids:
                    record = record_dict[match_id]
                    db_sequence = record.seq.__str__()[:read_len]
                    title = record.id.__str__()
                    for bacterium_id in bacteria_to_test:
                        bacterium_seq = bacteria_to_extract_df[(bacteria_to_extract_df['seq_id'] == bacterium_id) & (bacteria_to_extract_df['Region'] == region)]['Sequence'].iloc[0]
                        alignments_score = sum([c1 != c2 for c1, c2 in zip(bacterium_seq,db_sequence)])
                        bacteria_matches[bacterium_id].test_for_match(alignments_score, title, region)
                is_first_fasta = False
            else:
                record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
                for bacterium_id in bacteria_matches:
                    bacterium_seq = bacteria_to_extract_df[(bacteria_to_extract_df['seq_id'] == bacterium_id) & (bacteria_to_extract_df['Region'] == region)]['Sequence'].iloc[0]
                    for db_id in bacteria_matches[bacterium_id].get_ids_for_region(baseline_region):
                        try:
                            record = record_dict[db_id]
                            db_sequence = record.seq.__str__()[:read_len]
                            alignments_score = sum([c1 != c2 for c1, c2 in zip(bacterium_seq,db_sequence)])
                            bacteria_matches[bacterium_id].test_for_match(alignments_score, db_id, region)
                        except Exception as ex:
                            print "error with id = {}".format(db_id)

    db_id_by_bact_dict = []
    for bacterium_id in bacteria_matches:
        print bacterium_id
        for db_id in bacteria_matches[bacterium_id].get_db_ids():
            db_id_by_bact_dict.append({"bacteria_name": bacterium_id,
                                       "db_id":         db_id})


    df_db_id_by_bact = pd.DataFrame.from_dict(db_id_by_bact_dict)
    return df_db_id_by_bact


def test_bacteria_ids_from_database(region_db_path, bacteria_to_extract_df, test_match_ids, max_alignment_mismatch=5, read_len=145):
    bacteria_matches = {}
    min_score = max_alignment_mismatch

    bacteria_to_test = [bact for bact in bacteria_to_extract_df['seq_id'].unique() if bacteria_to_remove in bact]
    for bacterium_id in bacteria_to_test:
        bacteria_matches[bacterium_id] = BacteriumMatch(min_score)

    baseline_region = -1
    is_first_fasta = True
    for region_fasta in os.listdir(region_db_path):
        if region_fasta.endswith('fasta'):
            region = int(region_fasta.split('.')[0][-1])
            print region
            if is_first_fasta:
                baseline_region = region
                record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
                for match_id in test_match_ids:
                    record = record_dict[match_id]
                    db_sequence = record.seq.__str__()[:read_len]
                    title = record.id.__str__()
                    for bacterium_id in bacteria_to_test:
                        bacterium_seq = bacteria_to_extract_df[(bacteria_to_extract_df['seq_id'] == bacterium_id) & (bacteria_to_extract_df['Region'] == region)]['Sequence'].iloc[0]
                        alignments_score = sum([c1 != c2 for c1, c2 in zip(bacterium_seq,db_sequence)])
                        bacteria_matches[bacterium_id].test_for_match(alignments_score, title, region)
                is_first_fasta = False
            else:
                record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
                for bacterium_id in bacteria_matches:
                    bacterium_seq = bacteria_to_extract_df[(bacteria_to_extract_df['seq_id'] == bacterium_id) & (bacteria_to_extract_df['Region'] == region)]['Sequence'].iloc[0]
                    for db_id in bacteria_matches[bacterium_id].get_ids_for_region(baseline_region):
                        try:
                            record = record_dict[db_id]
                            db_sequence = record.seq.__str__()[:read_len]
                            alignments_score = sum([c1 != c2 for c1, c2 in zip(bacterium_seq,db_sequence)])
                            bacteria_matches[bacterium_id].test_for_match(alignments_score, db_id, region)
                        except Exception as ex:
                            print "error with id = {}".format(db_id)

    db_id_by_bact_dict = []
    for bacterium_id in bacteria_matches:
        print bacterium_id
        for db_id in bacteria_matches[bacterium_id].get_db_ids():
            db_id_by_bact_dict.append({"bacteria_name": bacterium_id,
                                       "db_id":         db_id})


    df_db_id_by_bact = pd.DataFrame.from_dict(db_id_by_bact_dict)
    return df_db_id_by_bact


def write_new_database(full_fasta_path, new_db_fasta_dir, db_id_by_bact_dict, bacteria_to_remove):
    new_db_name = "database_without_{}".format(bacteria_to_remove)
    new_db_path = os.path.join(new_db_fasta_dir, new_db_name)

    remove_ids=[]
    for ix, row in db_id_by_bact_dict.iterrows():
        if bacteria_to_remove in row['bacteria_name']:
            remove_ids.append(row["db_id"])

    remove_ids = [int(rid) if ((rid - int(rid)) == 0) else rid for rid in remove_ids]
    remove_ids = set([str(rid) for rid in remove_ids])

    print "{}: ids: {}".format(bacteria_to_remove, remove_ids)

    if not os.path.isdir(new_db_path):
        os.mkdir(new_db_path)
    for region_fasta in os.listdir(full_fasta_path):
        if region_fasta.endswith('fasta'):
            records_index = SeqIO.index(os.path.join(full_fasta_path, region_fasta),  "fasta")
            curr_region_keys = set(records_index.keys())
            new_region_keys = curr_region_keys - remove_ids
            new_records_index = (records_index[key] for key in new_region_keys)
            with open(os.path.join(new_db_path, region_fasta), "w") as new_fasta:
                SeqIO.write(new_records_index, new_fasta, "fasta")


def validate_ids(ids_path, zymo_seq_path, region_db_path, bacteria, read_len=145):
    zymo_seqs = pd.read_csv(zymo_seq_path)
    db_ids = pd.read_csv(ids_path)
    db_ids['bact'] = db_ids['bacteria_name'].apply(lambda x: True if bacteria in x else False)
    db_ids = db_ids[db_ids['bact'] == True]

    zymo_seqs = zymo_seqs[zymo_seqs['bact'].apply(lambda x: True if bacteria in x else False)]

    results = []

    for region_fasta in os.listdir(region_db_path):
        if region_fasta.endswith('fasta'):
            region = int(region_fasta.split('.')[0][-1])
            record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
            zymo_region = zymo_seqs[zymo_seqs['Region']==region]
            for _, z_row in zymo_region.iterrows():
                z_id = z_row['seq_id']
                z_seq = z_row['Sequence']
                curr_db_ids = db_ids[db_ids['bacteria_name']==z_id]['db_id'].tolist()
                for db_id in curr_db_ids:
                    try:
                        db_id = str(db_id if db_id - int(db_id) != 0 else int(db_id))
                        record = record_dict[db_id]
                        db_sequence = record.seq.__str__()[:read_len]
                        score = sum([c1 != c2 for c1, c2 in zip(z_seq,db_sequence)])
                        results.append({'region': region,
                                        'seq_id': z_id,
                                        'db_id': db_id,
                                        'score':score})
                    except Exception as ex:
                        print "failed with {}".format(db_id)


    pd.DataFrame.from_dict(results).to_csv(ids_path + "_validation.csv")







if __name__ == "__main__":
    parser = OptionParser("Convert SMURF2 results format to SMURF results format")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "")
    group_reqd.add_option("-o", "--output",
                      type="string", default="/home/vered/EMIRGE/data/s26_mock",
                      help="path to output new database location")
    group_reqd.add_option("-b", "--bacteria",
                          type="string", default="Lactobacillus_fermentum",
                          help="bacteria to remove")

    parser.add_option_group(group_reqd)
    (options, args) = parser.parse_args(sys.argv[1:])

    bacteria_to_remove = options.bacteria
    output_path = options.output
    db_ids_for_bacteria_path = os.path.join(output_path, "db_ids_for_bacteria_{}.csv".format(options.bacteria))
    print db_ids_for_bacteria_path

    zymo_regions_path = os.path.join(s26_results_path, 'zymo_regions.csv')
    print zymo_regions_path
    if not os.path.exists(zymo_regions_path):
        primers = Primer()
        primers.init(primers_path)
        results_df = primers.split_16s_to_regions(s26_results_path)
        # s26_ids_by_bacteria = get_16s_ids_for_26s_results(db_full_16S_path, s26_results_path)
        # results_df = compare_results_sequences(db_full_16s_regions, s26_ids_by_bacteria)

        results_df['bact'] = results_df['seq_id'].apply(lambda x: x.split('_16')[0])
        results_df.to_csv(zymo_regions_path)
    else:
        results_df = pd.read_csv(zymo_regions_path)

    validate_ids(ids_path_for_validation, zymo_regions_path, database_gg_regions, "Pseudomonas_aeruginosa")

    if not os.path.exists(db_ids_for_bacteria_path):
        db_ids_for_bacterium = get_bacteria_ids_from_database(database_gg_regions, results_df, bacteria_to_remove, max_alignment_mismatch=5, read_len=145)
        # db_ids_for_bacterium = test_1_bacteria_ids_from_database(database_gg_regions, results_df, match_ids, max_alignment_mismatch=5, read_len=145)
        db_ids_for_bacterium.to_csv(db_ids_for_bacteria_path)
    db_ids_for_bacterium = pd.read_csv(db_ids_for_bacteria_path)
    write_new_database(database_gg_regions, output_path, db_ids_for_bacterium, bacteria_to_remove)







