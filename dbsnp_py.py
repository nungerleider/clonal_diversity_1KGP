import pandas as pd
import glob
from scipy.stats import spearmanr, pearsonr

# Two sets of files run through the 'bedtools intersect' command with the dbSNP set of variants
# Overlapping sites were removed from analysis
files = glob.glob('ERR*/results/variants/variants.vcf.filtered.vaf0.1.try9.tsv') 
files2 = glob.glob('ERR*/results/variants/variants.vcf.dbsnp_filtered.vcf')

filter_string = '/results/variants/variants.vcf.filtered.vaf0.1.try9.tsv'
db_string = '/results/variants/variants.vcf.dbsnp_filtered.vcf'
names = [i.split('/')[0] for i in files2]  # Extract sample name from path

number_variants = {}
for file_id in names:
    filter_df = pd.read_table(file_id + filter_string)  # vcf file filtered using 'filter_vcy_1kgenomes.py
    snp_df = pd.read_table(file_id + db_string, header=None)

    filter_df.index = [f'{i}:{j}' for i,j in zip(filter_df['chromosome'], filter_df['position'])] # Set index to chr:pos to match following df
    snp_df.index = [f'{i}:{j}' for i,j in zip(snp_df[0], snp_df[1])] # Capture 'chr:pos' for each row (unique id for each loci) in index
    
    # Only loci that have passed both filteration steps are retained
    # Use set membership test to query loci that match both
    common = set(filter_df.index) & set(snp_df.index) 
    
    # extract 'common' (passed both tests) loci from snp_df
    snp_df = snp_df.loc[common]
    snp_df.to_csv(f'{file_id}.pass_dbsnp_and_filter_vcf_steps.vcf', sep='\t')

    # Count the number of filtered variants then store
    n_mut = len(snp_df.index) 
    number_variants[file_id] = n_mut

filtered_df = pd.DataFrame.from_dict(number_variants, orient='index') 

