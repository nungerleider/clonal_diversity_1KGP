import pandas as pd
import numpy as np
import glob

''' This script was used to determine mapped read coverage for each WGS sample from the 1000 Genomes Project
The coverage results were bimodal: samples either had < ~5 reads/base or > ~25 reads/base. For our mutation 
analysis, we required samples to have over 25 reads/base of coverage.

 fastq files were aligned using 'bwa mem -t 20' 
 bam indices were generated using 'samtools idxstats' '''

bams = glob.glob('/media/mcs/easystore2/1000_genomes_wgs/*bams/*idx')  # Go to dir w/ bam.idx files or manual type here

depth_dict = {}
human = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']

for bam in bams:
    df = pd.read_table(bam, index_col=0, header=None)

    # For each 'idxstats' output file:
    #
    # Index contains the chromosome name
    # First column contains the length of each chromosome
    # Second column contains mapped reads

    total_reads = np.sum(df.loc[human, 2])  # Total mapped reads
    human_chromosomal_length = np.sum(df.loc[human, 1])  # Total length of human genome (~3 billion bases)
    human_depth = total_reads / human_chromosomal_length
    human_depth *= 100  # Average read length is 100nt in this dataset
    depth_dict[bam.split('.')[0]] = human_depth
    
    print(bam)

df = pd.DataFrame.from_dict(depth_dict, orient='index')
df.columns = ['depth']
df.to_csv('sequencing_depth_1000_genomes_project.tsv', sep='\t')

