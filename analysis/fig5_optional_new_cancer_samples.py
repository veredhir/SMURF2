import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import logging
import time
import swifter

"""
23-01-2021
Hiseq_4_S_RDB284_Melanoma
Hiseq_4_S_RDB285_Melanoma
Hiseq_8_S_BC96_Lung
Ln_Batch5_S_Deborah3_S4_Melanoma
"""


def time_it(method):
    # logging.info('Start {}'.format(method.__name__))
    def timed(*args, **kw):
        logging.info('Entered %r' % method.__name__)
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        t = te - ts
        if t > 60:
            print('Done %s: %3.2f minutes\n' % (method.__name__, t/60))
        else:
            print('Done %s: %3.2f sec\n' % (method.__name__, t))
        return result
    return timed



def read_taxonomy_to_csv(taxonomy_path):
    taxonomy_df = pd.read_csv(taxonomy_path, sep='\t', skiprows=1)
    return taxonomy_df

@time_it
def to_int_swifter(df):
    for base in ['A', 'C', 'G', 'T']:
        df[base].update(df[base].swifter.apply(int))

@time_it
def to_int(df):
    for base in ['A', 'C', 'G', 'T']:
        df[base].update(df[base].apply(int))

@time_it
def to_int2(df):
    for base in ['A', 'C', 'G', 'T']:
        df[base].update(df[base].apply(lambda i: int(i)))
@time_it
def merge_data_left(df1, df2):
    return df1.merge(df2, how='left')


@time_it
def merge_mapping_taxa(mapping, taxa_df):
    return mapping.merge(taxa_df[['Class ', 'id']], left_on='Reference_id', right_on='id', how='left')


@time_it
def is_in(df, ids):
    return df[df.Unique_Reference_id.isin(ids)]


@time_it
def update_int_base(df1, df2):
    for base in ['A', 'C', 'G', 'T']:
        df1[base].update(df1[base].apply(int))
        df2[base].update(df2[base].apply(int))
    return df1, df2

@time_it
def or_dna(df):
    a = np.bitwise_and(df['A_ref'], df['A_read'])
    c = np.bitwise_and(df['C_ref'], df['C_read'])
    g = np.bitwise_and(df['G_ref'], df['G_read'])
    t = np.bitwise_and(df['T_ref'], df['T_read'])
    df['Score_bit'] = np.bitwise_or(np.bitwise_or(a, c), np.bitwise_or(g, t))
    df['Score'] = 0
    return df
    # logging.info("id {} done binary".format(g_id))

@time_it
def count_bits(df):
    # for iBit in range(0, 2*140 + 8, 8):
    #     df['Score_8lsb'] = df['Score_bit'] % (2 ** 8)
    #     df['Score'] = df['Score'] + np.unpackbits(
    #         [df['Score_8lsb'].astype('uint8')], axis=0).sum(axis=0)
    #     df['Score_bit'] = np.right_shift(df['Score_bit'], 8)
    df['Score'] = df['Score_bit'].apply(lambda x: bin(x).count("1"))

    return df




@time_it
def main():

    output_dir = '/home/vered/EMIRGE/data/real_data_fastq'

    results_summary = "Vered_re_analysis_HardDecision_Freq_ResultsSummary_SPECIES_cutFreq0.txt"
    species_cut_freq = "Vered_re_analysis_HardDecision_GroupsHeaders_SPECIES_cutFreq0.mat"
    read_counts = "Vered_re_analysis_HardDecision_ReadCountStats.txt"

    read_count = pd.read_csv(os.path.join(output_dir, read_counts), sep='\t')
    match_read_cols = ["Number of reads mathced to DB{}".format(ix) for ix in range(1, 6)]
    read_count['assigned_reads'] = read_count[match_read_cols].sum(axis=1)
    read_count_cols_to_save = ['Sample name', u'Number of loaded reads', u'Number of long reads', u'Number of good reads', 'assigned_reads']
    read_count[read_count_cols_to_save].to_csv(os.path.join(output_dir, "assigned_reads_tumers.csv"))

    smurf2_taxonomy_Hiseq_4_S_RDB284 = pd.read_csv(os.path.join(output_dir,"taxonomy_smurf2_Hiseq_4_S_RDB284.txt"),
                                                   sep='\t',
                                                   skiprows=1).rename(columns={'workdir/': 'Hiseq_4_S_RDB284'})
    smurf2_taxonomy_Hiseq_4_S_RDB285 = pd.read_csv(os.path.join(output_dir,"taxonomy_smurf2_Hiseq_4_S_RDB285.txt"),
                                                   sep='\t',
                                                   skiprows=1).rename(columns={'workdir/': 'Hiseq_4_S_RDB285'})

    smurf_taxonomy = pd.read_csv(os.path.join(output_dir, results_summary),
                                                   sep='\t',
                                                   skiprows=1).rename(columns={'2015-11-17 - Hiseq-4/RDB284_CCTGTCAA':'Hiseq_4_S_RDB284',
                                                                               '2015-11-17 - Hiseq-4/RDB285_CGACACTT':'Hiseq_4_S_RDB285'})

    merge_cols = [u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species']

    smurf2_taxonomy = smurf2_taxonomy_Hiseq_4_S_RDB284.merge(smurf2_taxonomy_Hiseq_4_S_RDB285, on=merge_cols, how='outer')
    comparison = smurf2_taxonomy.merge(smurf_taxonomy[merge_cols + ['Hiseq_4_S_RDB284', 'Hiseq_4_S_RDB285']],
                                       on=merge_cols,
                                       how='outer',
                                       suffixes=('_s2', '_s'))

    comparison = comparison.fillna(0)

    Hiseq_4_S_RDB284_comparison = comparison[merge_cols + ['Hiseq_4_S_RDB284_s2', 'Hiseq_4_S_RDB284_s']]
    Hiseq_4_S_RDB284_comparison['is_exists'] = (Hiseq_4_S_RDB284_comparison['Hiseq_4_S_RDB284_s2']>0) | (Hiseq_4_S_RDB284_comparison['Hiseq_4_S_RDB284_s'] >0)
    Hiseq_4_S_RDB284_comparison = Hiseq_4_S_RDB284_comparison[Hiseq_4_S_RDB284_comparison['is_exists']]
    Hiseq_4_S_RDB284_comparison_family = Hiseq_4_S_RDB284_comparison.groupby([u'domain', u'phylum', u'class', u'order', u'family'])[['Hiseq_4_S_RDB284_s2', 'Hiseq_4_S_RDB284_s']].sum().reset_index()
    Hiseq_4_S_RDB284_comparison_family.sort_values('Hiseq_4_S_RDB284_s2', ascending=False)\
        .to_csv(os.path.join(output_dir, "Hiseq_4_S_RDB284_comparison_family.csv"))


    Hiseq_4_S_RDB285_comparison = comparison[merge_cols + ['Hiseq_4_S_RDB285_s2', 'Hiseq_4_S_RDB285_s']]
    Hiseq_4_S_RDB285_comparison['is_exists'] = (Hiseq_4_S_RDB285_comparison['Hiseq_4_S_RDB285_s2'] > 0) | (
                Hiseq_4_S_RDB285_comparison['Hiseq_4_S_RDB285_s'] > 0)
    Hiseq_4_S_RDB285_comparison = Hiseq_4_S_RDB285_comparison[Hiseq_4_S_RDB285_comparison['is_exists']]
    Hiseq_4_S_RDB285_comparison_family = \
        Hiseq_4_S_RDB285_comparison.groupby([u'domain', u'phylum', u'class', u'order', u'family'])[
        ['Hiseq_4_S_RDB285_s2', 'Hiseq_4_S_RDB285_s']].sum().reset_index()
    Hiseq_4_S_RDB285_comparison_family.sort_values('Hiseq_4_S_RDB285_s2', ascending=False)\
        .to_csv(os.path.join(output_dir, "Hiseq_4_S_RDB285_comparison_family.csv"))



    smurf_taxonomy_results = "/home/vered/EMIRGE/data/real_data_fastq/tree_for_Vered.csv"
    taxonomy_dir = "/home/vered/EMIRGE/data/real_data_fastq/04-15-R02N1_S34"
    tested_sample = '04-15-R02N1_S34'
    skiprows = 2
    sep=","


    smurf2_dfs = []

    test_name = taxonomy_dir.split("/")[-1]
    curr_taxonomy_path = os.path.join(taxonomy_dir, 'taxonomy_smurf2.txt')
    taxonomy_df = read_taxonomy_to_csv(curr_taxonomy_path)
    mapped_reads=0
    with open(curr_taxonomy_path) as f:
        first_line = f.readline()
        mapped_reads = int(first_line.split('\t')[-1])
    print "smurf2: {}, mapped reads={}".format(test_name, mapped_reads)
    taxonomy_df.rename(columns={'workdir/': test_name}, inplace=True)

    taxonomy_df.fillna(0)
    smurf_taxonomy_result = pd.read_csv(smurf_taxonomy_results, skiprows=skiprows, sep=sep)
    smurf_taxonomy_result.rename(columns={'seq_names': 'species'}, inplace=True)
    smurf_taxonomy_result.fillna(0)
    print smurf_taxonomy_result.columns



    test_smurf2 = ''
    test_smurf = ''
    for col in taxonomy_df.columns:
        if tested_sample.replace("/", "_") in col:
            smurf2_cols = [u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species', col]
            test_smurf2=col
            print tested_sample, ": bacteria smurf2 ", len(taxonomy_df[taxonomy_df[col] != 0])
            break

    ix=-1
    for col in smurf_taxonomy_result.columns:
        ix += 1
        if tested_sample in col:
            smurf_cols = [u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species', col]
            test_smurf=col
            print tested_sample, ": bacteria smurf ", len(smurf_taxonomy_result[smurf_taxonomy_result[col] != 0])
            break
    with open(smurf_taxonomy_results) as smurf_f:
        for i in range(skiprows):
            first_line = smurf_f.readline()
            print "SMURF mapped_reads = ", first_line.split(sep)[ix]


    results_comparison = taxonomy_df[smurf2_cols].merge(smurf_taxonomy_result[smurf_cols],
                                                                   on=[u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species'],
                                                                   how='outer')

    results_comparison.fillna(0, inplace=True)

    results_comparison = results_comparison[(results_comparison[test_smurf]>0) | (results_comparison[test_smurf2]>0)]

    results_comparison[tested_sample+"_diff"] = (results_comparison[test_smurf2] - results_comparison[test_smurf]).abs()
    print len(results_comparison)
    print "large difference % ", len(results_comparison[results_comparison[tested_sample+"_diff"] > 0.7*results_comparison[[test_smurf2, test_smurf]].max(axis=1)])

    print "appears in only one results: ", len(results_comparison[
                  results_comparison[tested_sample + "_diff"] == results_comparison[[test_smurf2, test_smurf]].max(axis=1)])
    print "average abs diff = {}".format(results_comparison[tested_sample+"_diff"].mean())

    results_comparison.to_csv(os.path.join(output_dir, tested_sample.replace("/", "_") + "_smurf_vs_smurf2.csv"))
    table_for_figure = results_comparison.groupby([u'domain', u'phylum', u'class'])[[test_smurf, test_smurf2]].sum().reset_index()
    table_for_figure = table_for_figure.groupby(['class'])[[test_smurf, test_smurf2]].sum().reset_index()
    table_for_figure['Other'] = table_for_figure['04-15-R02N1_S34'] < 0.01

    table_for_figure.rename(columns={'Lane6/04-15-R02N1_S34': 'SMURF', '04-15-R02N1_S34': 'SMURF2'}, inplace=True)
    table_for_figure = table_for_figure.append({'class': 'Other',
                                                'SMURF': table_for_figure[table_for_figure['Other']]['SMURF'].sum(),
                                                'SMURF2': table_for_figure[table_for_figure['Other']]['SMURF2'].sum(),
                                                'Other': False}, ignore_index=True)
    table_for_figure = table_for_figure[table_for_figure.Other == 0][['class', 'SMURF', 'SMURF2']]

    smurf = table_for_figure[['class', 'SMURF']]
    smurf['type'] = 'SMURF'
    smurf.rename(columns={'SMURF': 'Frequency'}, inplace=True)
    smurf2 = table_for_figure[['class', 'SMURF2']]
    smurf2['type'] = 'SMURF2'
    smurf2.rename(columns={'SMURF2': 'Frequency'}, inplace=True)

    table_for_figure = smurf.append(smurf2)

    print table_for_figure
    ax = plt.figure(figsize=(14, 14))
    sns.set(font_scale=1.5)
    g = sns.barplot(data=table_for_figure, x='class', y='Frequency',  hue='type', palette='Paired')
    plt.setp(g.get_xticklabels(), rotation=10)
    ax.tight_layout()

    handles, labels = g.get_legend_handles_labels()
    g.legend(handles=handles, labels=labels)

    plt.ylim(0, 1)

    ax.savefig(os.path.join(taxonomy_dir, "Figure5a_04-15-R02N1_S34.png"))
    ax.savefig(os.path.join(taxonomy_dir, "Figure5a_04-15-R02N1_S34.pdf"))

    # plt.show()

# Go to taxa for vered
# Extract list of all the Gamma ids
# Go to smurf2 output 04-15-R02N1_S34
# Extract the relevant ids (merge)
# For each read mapped to Gamma: find the distance from the DB (from Gamma in the DB)
    ts = time.time()
    taxa_path = "/home/vered/EMIRGE/data/taxa_and_header.csv"
    taxa_df = pd.read_csv(taxa_path, sep='\t ')
    taxa_df['id2'] = taxa_df.index

    taxa_index = pd.read_csv('/home/vered/EMIRGE/data/reference_db/Header_uni_forVered.csv', header=None, names=['id'])
    taxa_index['id2'] = taxa_index.index

    taxa_df = taxa_df.merge(taxa_index)

    base_dtypes={'A': int, 'C': int, 'G': int, 'T': int}

    initial_mapping = pd.read_csv(os.path.join(taxonomy_dir, 'mapping.csv'))
    reference = pd.read_csv(os.path.join(taxonomy_dir, 'reference_db.csv'))
    ref_id_to_unique = pd.read_csv(os.path.join(taxonomy_dir, 'unique_ref_id_to_ref_id.csv'))

    te = time.time()
    t = te - ts
    print("Done reads", t, " sec")
    print ""

    mapping = initial_mapping.merge(ref_id_to_unique, how='left')

    mapping2 = mapping.merge(taxa_df[['Class ', 'id']], left_on='Reference_id', right_on='id', how='left')
    mapping2.drop_duplicates('Reads_group_id').groupby('Class ')['Count'].sum()

    gamma_reads = mapping2[mapping2['Class '] == 'Gammaproteobacteria']
    gamma_ids = gamma_reads.Unique_Reference_id.unique().tolist()

    reference_gamma = reference[reference.Unique_Reference_id.isin(gamma_ids)]

    gamma_reads, reference_gamma = update_int_base(gamma_reads, reference_gamma)

    reads_and_refs_chunk_df = pd.DataFrame.merge(reference_gamma, gamma_reads, on=['Unique_Reference_id', 'Region'], how='right',
                                                 suffixes=('_ref', '_read'))

    te = time.time()
    t = te - ts
    print("Done merges", t, " sec")
    print ""



    reads_and_refs_chunk_df = or_dna(reads_and_refs_chunk_df)

    reads_and_refs_chunk_df=count_bits(reads_and_refs_chunk_df)

    gamma_dist = 'Gammaproteobacteria\'s reads: Distance from the database'
    reads_and_refs_chunk_df[gamma_dist] = 140*2 - reads_and_refs_chunk_df['Score']

    data_for_hist = reads_and_refs_chunk_df[['Count', gamma_dist]]

    data_for_hist_without_count = data_for_hist.copy()

    for _, row in data_for_hist.iterrows():
        if row.Count > 1:
            data_for_hist_without_count = data_for_hist_without_count.append([row]*int(row.Count-1), ignore_index=True)

    ax = plt.figure(figsize=(14, 14))
    data_for_hist_without_count[gamma_dist] = data_for_hist_without_count[gamma_dist].astype(int)
    g = sns.countplot(data=data_for_hist_without_count, x=gamma_dist, color=sns.color_palette()[0])
    g = plt.setp(g.get_xticklabels(), fontsize=15)

    plt.text(1.5, 15000, 'SMURF', fontsize=15)  # add te
    # plt.text(13.5, 15000, 'SMURF2', fontsize=9)  # add te
    plt.axvline(2.5, color=sns.color_palette()[1], linestyle='--')
    # plt.axvline(14.5, color=sns.color_palette()[2], linestyle='--')
    ax.tight_layout()
    ax.savefig(os.path.join(taxonomy_dir, "Figure5b_hist_04-15-R02N1_S34.png"))
    ax.savefig(os.path.join(taxonomy_dir, "Figure5b_hist_04-15-R02N1_S34.pdf"))

    taxa_header_path = '/home/vered/EMIRGE/data/taxa_and_header.csv'

    print "mean = {}".format(data_for_hist_without_count[gamma_dist].mean())
    print "median = {}".format(data_for_hist_without_count[gamma_dist].median())
    print "% less 2 = {}".format(data_for_hist_without_count[data_for_hist_without_count[gamma_dist]<3][gamma_dist].count()/float(data_for_hist_without_count[gamma_dist].count()))
    print "sum less 2 = {}".format(
        data_for_hist_without_count[data_for_hist_without_count[gamma_dist] < 3][gamma_dist].count())
    print
    te = time.time()
    t = te - ts
    print("Done ", t, " sec")
    print ""


    print ""

if __name__ == "__main__":
    main()







