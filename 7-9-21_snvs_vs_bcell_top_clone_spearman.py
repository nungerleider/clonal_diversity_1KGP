import pandas as pd
import glob
from scipy.stats import spearmanr, pearsonr
import numpy as np
import matplotlib.pyplot as plt


# Two batches of files were downloaded, then later combined
files = glob.glob('ERR*/results/variants/variants.vcf.filtered.vaf0.1.try9.tsv') 
files2 = glob.glob('ERR*/results/variants/variants.vcf.dbsnp_filtered.vcf')

filter_string = '/results/variants/variants.vcf.filtered.vaf0.1.try9.tsv'
db_string = '/results/variants/variants.vcf.dbsnp_filtered.vcf'

# Sample IDs (strip path)
names = [i.split('/')[0] for i in files2] + [i.split('/')[0] for i in files]

number_variants = {}
for file_id in names:
    
    filter_df = pd.read_table(file_id + filter_string)
    snp_df = pd.read_table(file_id + db_string, header=None)
    filter_df.index = [f'{i}:{j}' for i,j in zip(filter_df['chromosome'], filter_df['position'])]
    snp_df.index = [f'{i}:{j}' for i,j in zip(snp_df[0], snp_df[1])]
    common = set(filter_df.index) & set(snp_df.index)
    snp_df = snp_df.loc[common]
    n_mut = len(snp_df.index)
    number_variants[file_id] = n_mut

filtered_df = pd.DataFrame.from_dict(number_variants,orient='index')

rna_directory = '/media/mcs/easystore1/1000_genomes_RNA/spreadsheets/'

#gene_expression_path = rna_directory + '1kgenomes_rna_gene_expression.tsv' # EBV genes in index have "EBV_" prefix
wgs_map_path = "/media/mcs/easystore2/1000_genomes_wgs/etc/sra_run_table.tsv" 
rna_map_path = rna_directory + 'rna_patient_map.tsv'

#gene_expression = pd.read_table(gene_expression_path, index_col=0)
wgs_map = pd.read_table(wgs_map_path)
rna_map = pd.read_table(rna_map_path)

wgs_dict = dict(zip(wgs_map['Run'], wgs_map['Sample Name']))
rna_dict = dict(zip(rna_map['RNA_barcode'], rna_map['sample_id']))

bams = glob.glob('/media/mcs/easystore/1000_genomes_wgs/*bams/*idx')  # Go to dir w/ bam.idx files or manual type here

ebv_dict = {}
ebv_cn_dict = {}
depth_dict = {}

human = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']

for bam in bams:
    df = pd.read_table(bam, index_col=0, header=None)
    total_reads = np.sum(df.loc[human,2])
    human_chromosomal_length = np.sum(df.loc[human,1])
    ebv_reads = np.sum(df.loc['gi|1386454078|gb|MG298921.1|', 2])
    ebv_length = df.loc['gi|1386454078|gb|MG298921.1|', 1]
    human_depth = total_reads / human_chromosomal_length
    human_depth *= 100  # Average read length is 100nt
    ebv_depth = 100* ebv_reads / ebv_length
    
    ebv_dict[bam.split('.')[0]] = ebv_depth
    ebv_cn_dict[bam.split('.')[0]] = ebv_depth / (human_depth * 2)
    depth_dict[bam.split('.')[0]] = human_depth
    
    print(bam)


df = pd.DataFrame.from_dict(depth_dict, orient='index')
df.index = [i.split('_')[0] for i in df.index]
filtered_df.columns = ['snv']
filtered_df['depth'] = df
filtered_df['sample'] = [wgs_dict[i] for i in filtered_df.index]
df['sample'] = [wgs_dict[i] for i in df.index]
snvs = filtered_df.groupby('sample').mean()
df = df.groupby('sample').mean()
snvs['depth'] = df[0]
snvs = snvs.dropna()
snvs[snvs['depth']>5]
snvs = snvs[snvs['depth']>5]

clone['sample' ] = clone.index.map(lambda x:rna_dict[x])
clone = clone.groupby('sample').mean()

snvs['top_clone'] = clone['top_clone_over_total_igh_counts']
snvs.to_csv('snvs_and_top_clone.tsv', sep='\t')
spearmanr(snvs['snv'], snvs['top_clone']) # Copied from terminal
snvs = snvs.dropna()
snvs.to_csv('snvs_and_top_clone.tsv',sep='\t')
spearmanr(snvs['snv'], snvs['top_clone'])  # Copied from terminal
snvs['clonality_index'] = clone['clonality_index']
spearmanr(snvs['snv'], snvs['cloneality_index']) # Copied from terminal
spearmanr(snvs['snv'], snvs['clonality_index']) # Copied from terminal
snvs['total_igh_counts'] = clone['total_igh_counts']
plt.scatter(snvs['top_clone'], snvs['snv'])

snvs2 = snvs[snvs['depth']>10]
spearmanr(snvs2['snv'], snvs2['top_clone']) # Copied from terminal
plt.scatter(snvs['top_clone'], snvs['snv'])

snvs3=snvs[snvs['depth']>15]
spearmanr(snvs3['snv'], snvs3['top_clone']) # Copied from terminal
plt.scatter(snvs3['top_clone'], snvs3['snv'])

snvs4 = filtered_df.reset_index().drop_duplicates('index')
snvs4.index = snvs4['sample']
snvs4['top_clone'] = clone['top_clone_over_total_igh_counts']
snvs4[snvs4['depth']>10]
snvs4['depth'] = df[0]
snvs4[snvs4['depth']>5]
snvs4 = snvs4[snvs4['depth']>5]
spearmanr(snvs4['snv'], snvs4['top_clone']) # Copied from terminal
snvs4.dropna()
snvs4 = snvs4.dropna()
spearmanr(snvs4['snv'], snvs4['top_clone']) # Copied from terminal
plt.scatter(snvs3['top_clone'], snvs3['snv'])
plt.scatter(snvs3['top_clone'], snvs3['snv'])
