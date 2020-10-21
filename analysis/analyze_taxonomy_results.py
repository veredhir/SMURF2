import pandas as pd
import os

def read_taxonomy_to_csv(taxonomy_path):
    taxonomy_df = pd.read_csv(taxonomy_path, sep='\t', skiprows=1)
    return taxonomy_df

if __name__ == "__main__":

    output_dir = '/home/vered/EMIRGE/data/real_data_fastq'


    smurf_taxonomy_results = "/home/vered/EMIRGE/data/real_data_fastq/tree_for_Vered.csv"
    smurf2_taxonomy_results = [
                               "/home/vered/EMIRGE/data/real_data_fastq/04-15-R02N1_S34",
                               "/home/vered/EMIRGE/data/real_data_fastq/06-15-T15S1_S46"]
    # "/home/vered/EMIRGE/data/real_data_fastq/06-15-T03N1_S7",
    tested_samples = ['04-15-R02N1_S34', '06-15-T15S1_S46']
    skiprows = 2
    sep=","

    # smurf_taxonomy_results = "/home/vered/EMIRGE/data/real_data_fastq/tax_for_vered_based_ALL_TMB_Nov19_NoIncluded.csv"
    # smurf2_taxonomy_results = ["/home/vered/EMIRGE/data/real_data_fastq/Batch184_RDB92_GCGTATCA",
    #                            "/home/vered/EMIRGE/data/real_data_fastq/Batch380_RDB100_TAAGTGGC"]
    # skiprows=1
    # sep = ','
    # tested_samples = ['Batch184/RDB92_GCGTATCA', 'Batch380/RDB100_TAAGTGGC']


    # old and irrelevant Ravid SMURF results.
    # smurf_taxonomy_results = "/home/vered/EMIRGE/data/real_data_fastq/TMB_upload_HardDecision_Freq_ResultsSummary_SPECIES_cutFreq0.txt"
    # smurf2_taxonomy_results = ["/home/vered/EMIRGE/data/real_data_fastq/Batch184_RDB92_GCGTATCA",
    #                            "/home/vered/EMIRGE/data/real_data_fastq/Batch380_RDB100_TAAGTGGC"]
    # skiprows=1
    # sep = '\t'

    # tested_samples = ['Batch184/RDB92_GCGTATCA']


    smurf2_dfs = []
    for taxonomy_dir in smurf2_taxonomy_results:
        test_name = taxonomy_dir.split("/")[-1]
        curr_taxonomy_path = os.path.join(taxonomy_dir, 'taxonomy_smurf2.txt')
        taxonomy_df = read_taxonomy_to_csv(curr_taxonomy_path)
        mapped_reads=0
        with open(curr_taxonomy_path) as f:
            first_line = f.readline()
            mapped_reads = int(first_line.split('\t')[-1])
        print "smurf2: {}, mapped reads={}".format(test_name, mapped_reads)
        taxonomy_df.rename(columns={'workdir/': test_name}, inplace=True)
        smurf2_dfs += [taxonomy_df]
    smurf2_taxonomy_result = pd.merge(smurf2_dfs[0], smurf2_dfs[1],  on=[u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species'], how='outer')
    smurf2_taxonomy_result.fillna(0)
    smurf_taxonomy_result = pd.read_csv(smurf_taxonomy_results, skiprows=skiprows, sep=sep)
    smurf_taxonomy_result.rename(columns={'seq_names': 'species'}, inplace=True)
    smurf_taxonomy_result.fillna(0)
    print smurf_taxonomy_result.columns


    for test in tested_samples:
        test_smurf2 = ''
        test_smurf = ''
        for col in smurf2_taxonomy_result.columns:
            if test.replace("/", "_") in col:
                smurf2_cols = [u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species', col]
                test_smurf2=col
                print test, ": bacteria smurf2 ", len(smurf2_taxonomy_result[smurf2_taxonomy_result[col] != 0])
                break

        ix=-1
        for col in smurf_taxonomy_result.columns:
            ix += 1
            if test in col:
                smurf_cols = [u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species', col]
                test_smurf=col
                print test, ": bacteria smurf ", len(smurf_taxonomy_result[smurf_taxonomy_result[col] != 0])
                break
        with open(smurf_taxonomy_results) as smurf_f:
            for i in range(skiprows):
                first_line = smurf_f.readline()
                print "SMURF mapped_reads = ", first_line.split(sep)[ix]


        results_comparison = smurf2_taxonomy_result[smurf2_cols].merge(smurf_taxonomy_result[smurf_cols],
                                                                       on=[u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species'],
                                                                       how='outer')

        results_comparison.fillna(0, inplace=True)

        results_comparison = results_comparison[(results_comparison[test_smurf]>0) | (results_comparison[test_smurf2]>0)]

        results_comparison[test+"_diff"] = (results_comparison[test_smurf2] - results_comparison[test_smurf]).abs()
        print len(results_comparison)
        print "large difference % ", len(results_comparison[results_comparison[test+"_diff"] > 0.7*results_comparison[[test_smurf2, test_smurf]].max(axis=1)])

        print "appears in only one results: ", len(results_comparison[
                      results_comparison[test + "_diff"] == results_comparison[[test_smurf2, test_smurf]].max(axis=1)])
        print "average abs diff = {}".format(results_comparison[test+"_diff"].mean())

        results_comparison.to_csv(os.path.join(output_dir, test.replace("/", "_") + "_smurf_vs_smurf2.csv"))
        # print results_comparison.groupby([u'domain', u'phylum', u'class', u'order', u'family', u'genus', u'species'])[[test_smurf, test_smurf2]].sum()
        print results_comparison.groupby([u'domain', u'phylum', u'class', u'order'])[[test_smurf, test_smurf2]].sum()
        print results_comparison.groupby([u'domain', u'phylum', u'class'])[[test_smurf, test_smurf2]].sum().reset_index()
        print results_comparison.groupby([u'domain', u'phylum'])[[test_smurf, test_smurf2]].sum()

        print ""


    print ""






