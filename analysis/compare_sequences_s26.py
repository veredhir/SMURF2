import pandas as pd
from Bio import SeqIO, pairwise2
import os
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt

db_full_16S_path = '/home/veredhi/EMIRGE_SMURF/GreenGenes_201305_unique_up_to_3_ambiguous_16S.fasta'
s26_results_path = '/home/veredhi/EMIRGE_SMURF/s26_results_16s/ssrRNAs'

s26_smurf2_results = '/home/veredhi/results/s26_mock_test2/tmp/workdir/emirge_smurf_WFalseSTrue.csv'
db_full_16s_regions = '/home/veredhi/EMIRGE_SMURF/ReferenceFullDB'
primers_path = ""

s26_results_path = "/home/vered/EMIRGE/data/s26_mock/ssrRNAs/"
primers_path = "/home/vered/EMIRGE/data/primers.csv"
GG_database_path = "/home/vered/EMIRGE/data/reference_db"

s26_smurf2_results_dir = "/home/vered/EMIRGE/data/s26_mock/dada2"
s26_smurf2_results1_dir = "/home/vered/EMIRGE/data/s26_mock/dada2/s26_without_Pseudomonas_aeruginosa/"
s26_smurf2_results2_dir = "/home/vered/EMIRGE/data/s26_mock/dada2/s26_without_Listeria_monocytogenes/"
s26_smurf2_results3_dir = "/home/vered/EMIRGE/data/s26_mock/dada2/s26_without_Lactobacillus_fermentum/"

s26_smurf2_dir_list = [s26_smurf2_results_dir, s26_smurf2_results1_dir, s26_smurf2_results2_dir, s26_smurf2_results3_dir]


# s26_smurf2_results = os.path.join(s26_smurf2_results_dir, "SMURF2_results.csv")


MAX_AMBIGUITY = 4


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


def compare_results_sequences(region_db_path, results_ids):
    s26_res_sequences = []
    for region_fasta in os.listdir(region_db_path):
        if region_fasta.endswith('fasta'):
            region = region_fasta.split('.')[0][-1]
            db_record_dict = SeqIO.index(os.path.join(region_db_path, region_fasta), "fasta")
            for bacteria in results_ids:
                for b_id in results_ids[bacteria]:
                    seq = db_record_dict[b_id].seq.__str__()
                    curr_dict = {'Sequence': seq,
                                 'Reference_id': id,
                                 'bacteria': bacteria,
                                 'Region': region}
                    s26_res_sequences.append(curr_dict)
    return pd.DataFrame.from_dict(s26_res_sequences)


def find_distance_between_bacteria(smurf2_bact, zymo_bact):
    bact_diff = 0
    diff_regions = []
    regions_counter= 0
    almost_equal_diff = 0
    almost_equal_regions = []
    for region in smurf2_bact.Region.unique():
        smurf2_seq = smurf2_bact[smurf2_bact.Region == region]['Sequence'].iloc[0]
        zymo_seq = zymo_bact[zymo_bact.Region == region]['Sequence'].iloc[0]
        # curr_diff = sum([c1 != c2 for c1, c2 in zip(zymo_seq, smurf2_seq)])

        alignments = pairwise2.align.globalxx(zymo_seq, smurf2_seq)
        curr_diff = len(zymo_seq) - alignments[0][2]
        bact_diff += curr_diff
        if curr_diff != 0:
            diff_regions.append(region)
        if curr_diff <= MAX_AMBIGUITY:
            almost_equal_diff+=curr_diff
            almost_equal_regions.append(region)
        regions_counter += 1

    return bact_diff, diff_regions, regions_counter, almost_equal_regions, almost_equal_diff


def find_distance_between_results(matches_df, zemo_df, smurf2_df):
    MAX_DIST = 1000
    diff_results = []
    unique_match = matches_df.drop_duplicates(['seq_id', 'changed_reference_id'])
    print 'unique matches = ', len(unique_match)
    for uid in unique_s26_smurf2.changed_reference_id.unique():
        dist = MAX_DIST
        for zymo_bact in unique_match[unique_match['changed_reference_id'] == uid]['seq_id']:
            zymo_df_bact = zemo_df[zemo_df['seq_id'] == zymo_bact]
            smurf2_df_bact = smurf2_df[smurf2_df['changed_reference_id'] == uid]
            dist, diff_regions, regions_counter, almost_equal_regions, almost_equal_diff = find_distance_between_bacteria(smurf2_df_bact, zymo_df_bact)
            diff_results.append({'smurf2_id': uid,
                                 'zymo_id': zymo_bact,
                                 'dist': dist,
                                 'diff_regions': diff_regions,
                                 '#diff_regions': len(diff_regions),
                                 'amplified_regions': regions_counter,
                                 'almost_equal_regions': len(almost_equal_regions),
                                 'almost_equal_diff': almost_equal_diff})
        if dist == MAX_DIST:
            best_dist = MAX_DIST
            best_match = None
            for zymo_bact in zemo_df.seq_id.unique():
                zymo_df_bact = zemo_df[zemo_df['seq_id'] == zymo_bact]
                smurf2_df_bact = smurf2_df[smurf2_df['changed_reference_id'] == uid]
                dist, diff_regions, regions_counter, almost_equal_regions, almost_equal_diff = find_distance_between_bacteria(smurf2_df_bact, zymo_df_bact)
                if dist < best_dist:
                    best_dist = dist
                    best_match = {'smurf2_id': uid,
                                  'zymo_id': zymo_bact,
                                  'dist': dist,
                                  'diff_regions': diff_regions,
                                  '#diff_regions': len(diff_regions),
                                  'amplified_regions': regions_counter,
                                  'almost_equal_regions': len(almost_equal_regions),
                                  'almost_equal_diff': almost_equal_diff}
            diff_results.append(best_match)
    return pd.DataFrame.from_dict(diff_results)


def get_distance_form_GG_single(db_index, sequence, db_id, read_len=145):
    if db_id - int(db_id) == 0:
        db_id = int(db_id)
    db_seq = db_index[str(db_id)].seq.__str__()[:read_len]
    smurf2_seq = "".join(sequence[:read_len])
    alignments = pairwise2.align.globalxx(db_seq, smurf2_seq)
    curr_diff = read_len - alignments[0][2]
    return curr_diff


def get_distance_from_GG(smurf2_results, gg_path):
    gg_indices = {}
    for i in range(1, 6):
        gg_region_path = os.path.join(gg_path, "GreenGenes_region{}.fasta".format(i))
        gg_indices[i] = SeqIO.index(gg_region_path, "fasta")

    dist_from_gg = smurf2_results.apply(lambda r: get_distance_form_GG_single(gg_indices[int(r['Region'])], r['Sequence'], r['Reference_id']), axis=1)
    return dist_from_gg


if __name__ == "__main__":
    zymo_regions_path = os.path.join(s26_results_path, 'zymo_regions.csv')
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
    distances = []
    comparison_without_bacteria = os.path.join(s26_smurf2_results_dir, "comparison_without_bacteria.csv")
    if not os.path.exists(comparison_without_bacteria) or True:
        comparison_res = []
        for s26_smurf2_dir in s26_smurf2_dir_list:
            s26_smurf2_results = os.path.join(s26_smurf2_dir, "SMURF2_results.csv")
            result_comparisons_dist_path = os.path.join(s26_smurf2_dir, "results_comparison_with_GG_dist.csv")
            if not os.path.exists(result_comparisons_dist_path):
                s26_smurf2_df = pd.read_csv(s26_smurf2_results)
                s26_smurf2_df = s26_smurf2_df.sort_values('Reference_id')
                unique_s26_smurf2 = s26_smurf2_df.drop_duplicates(['Region', 'changed_reference_id'])
                unique_s26_smurf2['dist_from_GG'] = get_distance_from_GG(unique_s26_smurf2, GG_database_path)
                unique_s26_smurf2['dist_from_GG'] = unique_s26_smurf2.groupby('changed_reference_id')['dist_from_GG'].transform('sum')
                unique_s26_smurf2['Sequence'] = unique_s26_smurf2['Sequence'].apply(lambda x: "".join(x[:145]))
                test = unique_s26_smurf2.merge(results_df, on=['Region', u'Sequence'])
                unique_s26_smurf2['count'] = unique_s26_smurf2.groupby('changed_reference_id')['Region'].transform('count')
                smurf2 = []
                for uid in unique_s26_smurf2.changed_reference_id.unique():
                    curr_dict = {'changed_reference_id': uid,
                                 'Prior': unique_s26_smurf2[unique_s26_smurf2['changed_reference_id'] == uid].iloc[0]['Priors']}
                    test_id = test[test.changed_reference_id==uid]
                    test_bact_id = test_id.groupby(['seq_id'])['Region'].count().reset_index()
                    test_bact_id['bact'] = test_bact_id['seq_id'].apply(lambda x: x.split('_16')[0])
                    test_bact_id = test_bact_id[test_bact_id['Region'] == test_bact_id['Region'].max()].drop_duplicates('bact')
                    for _, row in test_bact_id.iterrows():
                        curr_dict['#regions'] = row.Region
                        curr_dict['bact'] = row.bact
                        smurf2.append(curr_dict.copy())
                    if len(test_bact_id) == 0:
                        curr_dict['#regions'] = 0
                        curr_dict['bact'] = ''
                        smurf2.append(curr_dict.copy())
                smurf2_df = pd.DataFrame.from_dict(smurf2)
                smurf2_df['matches'] = smurf2_df.groupby('changed_reference_id')['bact'].transform('count')
                smurf2_df['matches'] = smurf2_df['matches']*(smurf2_df.bact != '')
                smurf2_df.groupby('bact')['Prior'].sum()
                distance_df = find_distance_between_results(test, results_df, unique_s26_smurf2)
                distance_df = distance_df.merge(unique_s26_smurf2[['changed_reference_id', 'Priors', 'dist_from_GG']].drop_duplicates(), left_on='smurf2_id', right_on='changed_reference_id')
                distance_df = distance_df.sort_values('dist').drop_duplicates('smurf2_id', keep='last')
                distance_df['identical_regions'] = distance_df['amplified_regions'] - distance_df['#diff_regions']
                distance_df['bact'] = distance_df['zymo_id'].apply(lambda x: x.split('_16')[0])


                distance_df.sort_values('Priors', ascending=False)[['dist', 'dist_from_GG', u'smurf2_id','Priors', u'zymo_id', 'amplified_regions', u'identical_regions', 'almost_equal_regions']]\
                    .to_csv(result_comparisons_dist_path, index=False)
            else:
                distance_df = pd.read_csv(result_comparisons_dist_path)
                distance_df['bact'] = distance_df['zymo_id'].apply(lambda x: x.split('_16')[0])
            distance_df['TEST_ID'] = s26_smurf2_dir.split('/')[-2]
            distances.append(distance_df)
            summarized_results = distance_df[distance_df.almost_equal_regions != 0].groupby('bact')['Priors'].sum()
            print summarized_results
            comparison_res.append(summarized_results)

        distances_df = pd.concat(distances)
        print "distance_df[distance_df.almost_equal_regions != 0].Priors.sum() {}".format(distances_df[distances_df.almost_equal_regions != 0].Priors.sum() )
        print len(distances_df[distances_df.almost_equal_regions != 0])
        distances_df = distances_df[distances_df.almost_equal_regions != 0]
        distances_df['bact_prior'] = distances_df.groupby(['bact', 'TEST_ID'])['Priors'].transform('sum')
        distances_df['w_dist_GG_bact'] = distances_df['dist_from_GG']*distances_df['Priors'] / distances_df['bact_prior']
        distances_df['w_dist_bact'] = distances_df['dist']*distances_df['Priors'] / distances_df['bact_prior']
        distances_df['zymo_prior'] = distances_df.groupby(['zymo_id', 'TEST_ID'])['Priors'].transform('sum')
        distances_df['w_dist_GG_zymo'] = distances_df['dist_from_GG'] * distances_df['Priors'] / distances_df['zymo_prior']
        distances_df['w_dist_zymo'] = distances_df['dist'] * distances_df['Priors'] / distances_df['zymo_prior']
        dist_test = distances_df.groupby(['bact', 'TEST_ID'])['w_dist_GG_bact', 'w_dist_bact', 'w_dist_GG_zymo', 'w_dist_zymo'].sum()
        dist_test = dist_test.reset_index()

        dist_test2 = distances_df.groupby(['zymo_id', 'TEST_ID'])['w_dist_GG_zymo', 'w_dist_zymo', 'Priors'].sum()
        dist_test2.to_csv(os.path.join(s26_smurf2_results_dir, "dist_by_zymo_id_comparison.csv"))

        dist_test = distances_df.groupby(['bact', 'TEST_ID'])['w_dist_GG_bact', 'w_dist_bact', 'Priors'].sum()
        dist_test.to_csv(os.path.join(s26_smurf2_results_dir, "dist_by_bact_comparison.csv"))

        test = pd.DataFrame.from_dict(comparison_res)
        test.to_csv(comparison_without_bacteria)
    else:
        test = pd.read_csv(comparison_without_bacteria)



    test = test*100
    test['SMURF2_test'] = ['Full DB', 'without Pseudomonas_aeruginosa', 'without Listeria_monocytogenes',
                           'without Lactobacillus_fermentum']
    test['SMURF2_test'] = ['Full Greengenes DB',
                           'Greengenes without\nPseudomonas aeruginosa',
                           'Greengenes without\nListeria monocytogenes',
                           'Greengenes without\nLactobacillus fermentum']
    test.index = test.SMURF2_test
    test = test[[ u'Bacillus_subtilis', u'Enterococcus_faecalis', u'Escherichia_coli', u'Lactobacillus_fermentum',u'Listeria_monocytogenes', u'Pseudomonas_aeruginosa', u'Salmonella_enterica', u'Staphylococcus_aureus']]


    fig = plt.figure()
    ax = plt.subplot(111)
    test.plot(kind='bar', stacked=True, ax=ax)

    # ax.set_title('SMURF2 detection of out of the database bacteria', fontsize=12)
    ax.set_ylabel("Bacterium frequency (%)", fontsize=12)
    ax.set_xlabel("Database used for reconstruction", fontsize=12, y=-0.3)

    # ax.legend(bbox_to_anchor=(1.01, 1.0), loc='upper left')
    plt.xticks(rotation=0)

    # ax.set_xticks(rotation=0)

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              fancybox=True, shadow=True, ncol=4)
    plt.show()
    # plt.savefig("sample.jpg")

    print len(test)
    a=3


