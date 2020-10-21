

import pandas as pd
import os






def compare_single_bacteria(smurf2_comparison_results_dir, smurf2_comparison_results_without_dir, bacteria_id):
    smurf2_comparison_results = os.path.join(smurf2_comparison_results_dir, "results_comparison_with_GG_dist.csv")
    smurf2_comparison_results_without = os.path.join(smurf2_comparison_results_without_dir, "results_comparison_with_GG_dist.csv")

    smurf2_comparison_results_seq = os.path.join(smurf2_comparison_results_dir, "SMURF2_results.csv")
    smurf2_comparison_results_without_seq = os.path.join(smurf2_comparison_results_without_dir, "SMURF2_results.csv")

    bacteria_ids = [bacteria_id.format(i) for i in range(10)]
    df_full = pd.read_csv(smurf2_comparison_results)
    df_without = pd.read_csv(smurf2_comparison_results_without)

    df_full_bact = df_full[df_full['zymo_id'].isin(bacteria_ids)]
    df_without_bact = df_without[df_without['zymo_id'].isin(bacteria_ids)]

    df_seq_full = pd.read_csv(smurf2_comparison_results_seq).merge(df_full_bact, right_on='smurf2_id',
                                                                   left_on='changed_reference_id')
    df_seq_without = pd.read_csv(smurf2_comparison_results_without_seq).merge(df_without_bact, right_on='smurf2_id',
                                                                              left_on='changed_reference_id')

    common_db_ids = list(set(df_seq_full.Reference_id.to_list()).intersection(df_seq_without.Reference_id.to_list()))
    full_common_freq = df_seq_full[df_seq_full['Reference_id'].isin(common_db_ids)].drop_duplicates('smurf2_id')[
        'Priors_x'].sum()
    without_common_freq = \
    df_seq_without[df_seq_without['Reference_id'].isin(common_db_ids)].drop_duplicates('smurf2_id')[
        'Priors_x'].sum()

    print "frequency of the common DB ids: Full --> {}, without -->{}".format(full_common_freq, without_common_freq)
    print "frequency of bacteria: Full --> {}, without --> {}".format(df_full_bact.Priors.sum(),
                                                                      df_without_bact.Priors.sum())
    print "amount of mapped DB bacteria: Full --> {}, without --> {}".format(
        len(df_seq_full['Reference_id'].drop_duplicates()), len(df_seq_without['Reference_id'].drop_duplicates()))
    print "amount of found bacteria: Full --> {}, without --> {}".format(
        len(df_seq_full['changed_reference_id'].drop_duplicates()),
        len(df_seq_without['changed_reference_id'].drop_duplicates()))

    df_seq_full_pivot = df_seq_full.drop_duplicates(['changed_reference_id', 'Region']).pivot(
        index='changed_reference_id', columns='Region', values='Sequence').reset_index().fillna('12345')
    df_seq_without_pivot = df_seq_without.drop_duplicates(['changed_reference_id', 'Region']).pivot(
        index='changed_reference_id', columns='Region', values='Sequence').reset_index().fillna('12345')

    for i in range(1, 6, 1):
        if i not in df_seq_full_pivot.columns:
            df_seq_full_pivot[i] = "12345"
        if i not in df_seq_without_pivot.columns:
            df_seq_without_pivot[i] = "12345"

    common_seq = df_seq_without_pivot.merge(df_seq_full_pivot, on=range(1, 6, 1))
    common_seq_freq_full = df_seq_full[
        df_seq_full.changed_reference_id.isin(common_seq.changed_reference_id_y.unique().tolist())]
    common_seq_all_regions_freq_full = common_seq_freq_full.drop_duplicates('changed_reference_id').Priors_x.sum()
    common_seq_freq_without = df_seq_without[
        df_seq_without.changed_reference_id.isin(common_seq.changed_reference_id_x.unique().tolist())]
    common_seq_all_regions_freq_without = common_seq_freq_without.drop_duplicates('changed_reference_id').Priors_x.sum()

    print "Frequency of common Sequence over all regions: Full {}, without {}".format(common_seq_all_regions_freq_full,
                                                                                      common_seq_all_regions_freq_without)

    seq_full = df_seq_full.Sequence.unique().tolist()
    seq_without = df_seq_without.Sequence.unique().tolist()

    common_seq = list(set(seq_full).intersection(seq_without))

    print "Sequence: Full {}, without {}, common {}".format(len(seq_full), len(seq_without), len(common_seq))

    df_seq_full = pd.read_csv(smurf2_comparison_results_seq)
    df_seq_without = pd.read_csv(smurf2_comparison_results_without_seq)
    df_seq_full_pivot = df_seq_full.drop_duplicates(['changed_reference_id', 'Region']).pivot(
        index='changed_reference_id', columns='Region', values='Sequence').reset_index().fillna('qqq')
    df_seq_without_pivot = df_seq_without.drop_duplicates(['changed_reference_id', 'Region']).pivot(
        index='changed_reference_id', columns='Region', values='Sequence').reset_index().fillna('qqq')
    print ("Total amount of sequences full = {}, without = {}".
           format(len(df_seq_full.drop_duplicates('changed_reference_id')),
                  len(df_seq_without.drop_duplicates('changed_reference_id'))))

    for i in range(1, 6, 1):
        if i not in df_seq_full_pivot.columns:
            df_seq_full_pivot[i] = "qqq"
        if i not in df_seq_without_pivot.columns:
            df_seq_without_pivot[i] = "qqq"

    common_seq = df_seq_without_pivot.merge(df_seq_full_pivot, on=range(1, 6, 1))
    common_seq_freq_full = df_seq_full[
        df_seq_full.changed_reference_id.isin(common_seq.changed_reference_id_y.unique().tolist())]
    common_seq_all_regions_freq_full = common_seq_freq_full.drop_duplicates('changed_reference_id').Priors.sum()
    common_seq_freq_without = df_seq_without[
        df_seq_without.changed_reference_id.isin(common_seq.changed_reference_id_x.unique().tolist())]
    common_seq_all_regions_freq_without = common_seq_freq_without.drop_duplicates('changed_reference_id').Priors.sum()
    print df_seq_without.drop_duplicates('changed_reference_id').Priors.sum()

    print "Frequency of common Sequence over all regions: Full {}, without {}".format(common_seq_all_regions_freq_full,
                                                                                      common_seq_all_regions_freq_without)


if __name__ == "__main__":
    smurf2_comparison_results_dir = "/home/vered/EMIRGE/data/s26_mock/dada2"

    smurf2_comparison_results_without_dirs = ["/home/vered/EMIRGE/data/s26_mock/dada2/s26_without_Pseudomonas_aeruginosa",
                                              "/home/vered/EMIRGE/data/s26_mock/dada2/s26_without_Listeria_monocytogenes",
                                              "/home/vered/EMIRGE/data/s26_mock/dada2/s26_without_Lactobacillus_fermentum/"]

    bacteria_templates = ["Pseudomonas_aeruginosa_16S_{}", "Listeria_monocytogenes_16S_{}", 'Lactobacillus_fermentum_16S_{}']

    for bid, bdir in zip(bacteria_templates, smurf2_comparison_results_without_dirs):
        print "\n\nbacteria: {}\n".format(bid)
        compare_single_bacteria(smurf2_comparison_results_dir, bdir, bid)



    a=5