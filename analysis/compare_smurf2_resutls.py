import pandas as pd


def get_agg_data(df):
    return df.groupby(['exp-index', 'score-type', 'value-type', 'bacteria'])['value'].mean().reset_index()


path1 = ""
path2 = ""

old_res = pd.read_csv(path1)
new_res = pd.read_csv(path2)

agg_old = get_agg_data(old_res)
agg_new = get_agg_data(new_res)

print agg_new.merge(agg_old, on=['score-type', 'value-type', 'bacteria'], suffixes=('_new', '_old'))