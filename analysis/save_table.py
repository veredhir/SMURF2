import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import table # EDIT: see deprecation warnings below
pd.options.display.latex.repr=True
ax = plt.subplot(111, frame_on=False) # no visible frame
ax.xaxis.set_visible(False)  # hide the x axis
ax.yaxis.set_visible(False)  # hide the y axis
df = pd.read_csv('/home/vered/EMIRGE/data/real_data_fastq/table1_compare_assigned_reads.csv', header=[0])

table(ax, df)  # where df is your data frame
plt.show()
plt.savefig('mytable.png')