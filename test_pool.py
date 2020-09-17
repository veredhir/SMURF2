import pandas as pd
import itertools
from multiprocessing import  Pool, cpu_count

df = pd.DataFrame.from_dict({'A': [1, 1, 1, 2, 3, 4, 5, 5], 'B': [1, 1, 1, 5, 5, 5, 5, 5]})
grpby = df.groupby('A')

z, k, l, m = "vered", 4, 6.7, pd.DataFrame({'K':[1, 1, 1], 'M':['test', 'tets', 'vered']})


def test_func(((id, grp), z, k, l, m)):
    return id + l

pool = Pool(1)
ret = pool.map(test_func, itertools.izip(grpby, itertools.repeat(z), itertools.repeat(k), itertools.repeat(l), itertools.repeat(m)))


print ret