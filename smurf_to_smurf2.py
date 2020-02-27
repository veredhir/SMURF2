
from Bio import SeqIO
import numpy as np
from emirge_headers import *
from emirge_utills import *
import os
import seaborn as sns
import scipy.io as sio

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


def add_taxa_to_results(taxa_path, results_df):
    taxa_df = pd.read_csv(taxa_path, sep='	 ', index_col=False)
    rename_dict = {'Domain':'domain',
                   'Phylum ': 'phylum',
                   'Class ': 'class',
                   'Order ': 'order',
                   'Family ': 'family',
                   'Genus ': 'genus',
                   'Species': 'seq_names'}
    taxa_df = taxa_df.rename(rename_dict)

    results_df['Header'] = results_df[Header.ref_id].apply(lambda r: r if '#' not in r else r.split('#')[0])
    df = results_df.merge(taxa_df, on='Header')

    return df



class Header():
    ref_id = 'ref'
    prior = 'prior'
    sequence = 'sequence'
    region = 'region'
    weight = 'weight'
    is_changed = 'is_changed'
    new_id = 'new_id'


def write_bacterium_to_fasta(new_bacteria_path,
                             new_bacteria):
    """
    :param new_bacteria_path:
    :param new_bacteria: dictionary, keys are the ids and values are the sequences
    :return:
    """
    with open(new_bacteria_path, 'w') as mock_fasta1:
        for id in new_bacteria.keys():
                mock_fasta1.write(">{}\n".format(id))
                mock_fasta1.write(new_bacteria[id] + '\n')
                mock_fasta1.write('\n')

def validate_priors(df, threshold=VALIDATION_THRESHOLD_THRESHOLD):
    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])
    df.prior = df.prior.apply(lambda r: r / sum_prior)

    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])

    # logging.debug("sum of priors is {}".format(sum_prior))
    if abs(sum_prior - 1) > threshold:
        raise Exception("sum of prior is not 1")


def write_new_bacterium_to_fasta(new_bacteria_df,
                                 regions,
                                 fasta_dir):
    regions = range(1, regions+1)
    for region in regions:
        if region not in new_bacteria_df.columns:
            continue
        fasta_name = "new_bacteria_{}.fasta".format(region)
        fasta_file = os.path.join(os.path.join(fasta_dir, fasta_name))
        with open(fasta_file, 'w') as mock_fasta1:
            for _, row in new_bacteria_df.iterrows():
                id = row[Header.new_id]
                seq = row[region]
                if seq != '':
                    mock_fasta1.write(">{}\n".format(id))
                    mock_fasta1.write(seq + '\n')
                    mock_fasta1.write('\n')


def are_similar(seq1, db_seq):
    length = len(seq1)/2
    seq2 = db_seq[:length] + db_seq[-1*length:]
    if seq2 == seq1:
        return True
    if 'N' in seq2:
        N_pos = [pos for pos, char in enumerate(seq2) if char == 'N']
        for pos in N_pos:
            seq1[pos] = 'N'
        if seq2 == seq1:
            return True
    return False


class Found_bacteria():
    frequency=[]
    assigned_reads=[]

class BacteriaMetaGroup():
    def __init__(self, db_ind, comb_vector, new_ind, is_out_of_db):
        self.db_ind=db_ind
        self.comb_vector=comb_vector


def read_db_ind_map(path_to_db_ind="/home/vered/EMIRGE/data/reference_db/Header_uni_forVered.csv"):
    db_map = pd.read_csv(path_to_db_ind, index_col=None, header=None)
    db_map['fasta_id'] = db_map[0].apply(lambda x: str(int(x)))
    return db_map


def convert_format(mat_path, fasta_dir, output_dir, read_length=126, num_regions=5):
    """
    :param path: path to 'final_results.csv produced by emirge_smurf.py
    :return: df hold the final results
    """
    matstruct_contents = sio.loadmat(mat_path)

    bacteriaMetaGroups = matstruct_contents['bactMetaGroups']
    found_bacteria = matstruct_contents['found_bacteria']
    comb_vecs=bacteriaMetaGroups['comb_vec'][0]
    db_ind = bacteriaMetaGroups['db_ind'][0]
    frequency = found_bacteria['frequency'][0][0][0]
    df = pd.DataFrame(columns=['Sequence', HeadersFormat.Region, HeadersFormat.Priors, 'Unique_Reference_id'])
    db_map = read_db_ind_map()

    print "db map length = %d", len(db_map)

    for comb_vectors, ids, prior in zip(comb_vecs, db_ind, frequency):
        id = ids[0]
        comb_vector = comb_vectors[0]
        for ix in range(len(comb_vector)):
            if comb_vector[ix] == 1:
                df = df.append({'Sequence': '',
                            HeadersFormat.Region: ix+1,
                            HeadersFormat.Priors: prior,
                            'Unique_Reference_id': id}, ignore_index=True

                )
    new_df = pd.DataFrame()
    fasta_files = get_fasta_files(fasta_dir, num_regions)
    for fasta_file in fasta_files:
        records = SeqIO.index(fasta_file.path, "fasta")
        curr_region = fasta_file.region
        print "\n\n\n**\n**\ncurr region: "
        print curr_region, ": ", fasta_file.path
        for index, row in df.iterrows():
            if row[HeadersFormat.Region] != curr_region:
                continue
            ids = row['Unique_Reference_id']
            db_seq = extract_from_records(records, ids, db_map)
            if(len(db_seq) > 1):
                size=len(db_seq)-1
                new_df = new_df.append(pd.DataFrame.from_dict({HeadersFormat.Region:[curr_region]*size,
                                                               HeadersFormat.Priors:[row[HeadersFormat.Priors]]*size,
                                                               'Unique_Reference_id': [ids[0]]*size,
                                                               'Sequence': db_seq[1:]}))
            seq = db_seq[0]
            df.loc[index, 'Sequence'] = seq
            df.loc[index, 'Unique_Reference_id'] = ids[0]
        records.close()

    print ("original size = {} added ={}".format(len(df), len(new_df)))
    df = df.append(new_df)
    res_file_name = "emirge_smurf_WFalseSFalse.csv"
    smurf2_path = os.path.join(output_dir, res_file_name)
    df.to_csv(smurf2_path)


def extract_from_records(records, ids, db_map):
    res = []
    for id in ids:
        try:
            fasta_id = str(db_map.loc[id-1, 'fasta_id'])
            seq = records[fasta_id].seq
            db_seq = seq.__str__()
            res.append(db_seq[:126] + db_seq[-126:])
        except Exception as ex:
            for i in range(60):
                ix_end = ".%02d" % (i,)
                new_fasta_id = fasta_id + ix_end
                if new_fasta_id in records.keys():
                    seq = records[new_fasta_id].seq
                    db_seq = seq.__str__()
                    res.append( db_seq[:126] + db_seq[-126:])
    res = list(set(res))
    return res
    # return res[0]


class FastaFile(object):
    def __init__(self, directory, file_name, total_num_of_regions=None, region=None):
        self.path = os.path.join(directory, file_name)
        self.file_name = file_name
        self.region = region
        self.total_num_of_regions = total_num_of_regions

    def initialize(self):
        if self.region == None:
            if self.total_num_of_regions is None:
                raise Exception("Missing 'total number of regions'")
            self.region = self.get_region_from_file_name()

    def get_region_from_file_name(self):
        region = -1
        for i in range(1, self.total_num_of_regions + 1):
            if str(i) in self.file_name:
                region = i
        return region


def get_fasta_files(fasta_dir, total_num_of_regions):
    fasta_files = []
    os.chdir(fasta_dir)
    files = glob.glob("*.fasta")
    for file in files:
        fasta_file = FastaFile(fasta_dir, file, total_num_of_regions)
        fasta_file.initialize()
        fasta_files.append(fasta_file)
    return fasta_files


def main(argv = sys.argv[1:]):
    """
    command line interface to emirge

    """
    parser = OptionParser("Convert SMURF2 results format to SMURF results format")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "")
    group_reqd.add_option("-i", "--input",
                      type="string",
                      help="path to smurf results .mat directory")
    group_reqd.add_option("-f", "--fasta",
                          type="string",
                          help="path to regions fasta directory")
    group_reqd.add_option("-o", "--output",
                          type="string",
                          help="path to output directory")
    parser.add_option_group(group_reqd)
    (options, args) = parser.parse_args(argv)

    input_dir=options.input
    fasta_dir=options.fasta
    output_dir=options.output

    for input_mat in os.listdir(input_dir):
        if not input_mat.endswith('mat'):
            continue
        res_dir_path = output_dir
        # res_dir_name = "test_"  + input_mat.split("sample_mock_ix_")[1].split("_results")[0]
        # res_dir_path = os.path.join(output_dir, res_dir_name)
        # print(res_dir_path)
        # os.mkdir(res_dir_path)
        # input_mat = os.path.join(input_dir, input_mat)

        convert_format(os.path.join(input_dir,input_mat), fasta_dir, res_dir_path)


if __name__ == "__main__":
    main()



