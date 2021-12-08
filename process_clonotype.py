import pandas as pd
import glob
import matplotlib.pyplot as plt

'''Tabulate and process clonotype output for each RNA-Seq sample.
"run_mixcr.sh" script/pipeline was used to generate clonotype info for each sample'''


clonotype_output_files = glob.glob('/media/mcs/easystore2/1000_genomes_RNA/*clones.tsv')

top_clone_read_fraction = {}
number_of_distinct_clones = {}
reads_of_top_clone = {}

for i in clonotype_output_files:
    sample_id = i.split('.')[0].split('/')[-1]
    x = pd.read_table(i, index_col=0)
    reads_of_top = x.iloc[0,0]
    total_reads = np.sum(x.iloc[:,0])
    top_fraction = reads_of_top / total_reads
    clones = len(x.index)

    reads_of_top_clone[sample_id] = reads_of_top
    top_clone_read_fraction[sample_id] = top_fraction
    number_of_distinct_clones[sample_id] = clones

df = pd.DataFrame.from_dict(reads_of_top_clone, orient='index')
df.columns = ['reads_of_top_clone']
df['top_clone_read_fraction'] = df.index.map(lambda x:top_clone_read_fraction[x])
df['number_of_distinct_clones'] = df.index.map(lambda x:number_of_distinct_clones[x])

df.to_csv('11_12_21_clone_data.tsv', sep='\t')

# Figure xx
fig = plt.figure()
ax = plt.subplot()
df = df.sort_values('reads_of_top_clone')
plt.bar(range(len(df.index)), df['reads_of_top_clone'])
plt.savefig('reads_of_top_clone_barplot.svg')

# Figure xx
fig = plt.figure()
ax = plt.subplot()
df = df.sort_values('number_of_distinct_clones')
plt.bar(range(len(df.index)), df['number_of_distinct_clones'])
plt.savefig('number_of_distinct_clones.svg')
