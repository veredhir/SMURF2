
import pandas as pd
import os
from smurf2_iteration import Smurf2Iteration


def test_filter_similar_mapped_reference(test_path):
    s2 = Smurf2Iteration(test_path)
    mapping_df = pd.read_csv(s2.paths.mapping)
    references = pd.read_csv(s2.paths.reference)
    s2.filter_similar_mapped_reference(mapping_df, references)


def test_find_initial_weight(test_path):
    weight_as_detected_regions = True
    s2 = Smurf2Iteration(test_path)
    s2.update_reference_db_with_weight(weight_as_detected_regions)


def test_update_ref(test_path, read_len):
    s2 = Smurf2Iteration(test_path)

    s2.iteration_index = 3
    s2.read_len = read_len
    s2.th.min_coverage_for_split = 0.005
    s2.allow_split = True

    rename_dict = {}
    for ix in range(2*read_len):
        for base in ['A', 'C', 'G', 'T']:
            rename_dict['{}{}_prob_n'.format(ix, base)]= str(ix)

    prob_n_dict = {
        'A': pd.read_csv(os.path.join(test_path, "prob_A.csv")).rename(columns=rename_dict),
        'C': pd.read_csv(os.path.join(test_path, "prob_C.csv")).rename(columns=rename_dict),
        'G': pd.read_csv(os.path.join(test_path, "prob_G.csv")).rename(columns=rename_dict),
        'T': pd.read_csv(os.path.join(test_path, "prob_T.csv")).rename(columns=rename_dict)
    }
    s2.update_bacteria_DB(prob_n_dict)


if __name__ == '__main__':
    test_path1 = "/home/vered/EMIRGE/data/test4/iter.01"
    test_update_ref(test_path1, read_len=126)

    # test_find_initial_weight(test_path1)

    # test_filter_similar_mapped_reference(test_path1)