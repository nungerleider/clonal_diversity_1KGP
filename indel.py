import pandas as pd
import sys

''' Using the variants.vcf output file from Strelka analysis, retain only insertions or deletions'''

fh = sys.argv[1]
vcf = pd.read_table(fh, comment='#', header=None)
vcf.columns = ['chromosome', 'position', '_','ref', 'alt', '_','pass_fail', 'info1', 'info2', 'info3']

classification = []
for ref, alt in zip(vcf['ref'], vcf['alt']):

    # if vcf file had entries like ref:'A', alt:'ATGG'  
    if len(ref) < len(alt):
        classification.append('insertion')
    
    # if vcf file had entries like ref:'A', alt:''  
    elif len(ref) > len(alt):
        classification.append('deletion')

    # if vcf files had entries like ref:'A', alt: 'T'
    else:
        classification.append('SNV')

vcf['mut_type'] = classification
vcf['alt'] = vcf['info3'].map(lambda x:int(x.split(':')[5].split(',')[1]))
vcf['ref'] = vcf['info3'].map(lambda x:int(x.split(':')[5].split(',')[0]))
vcf['total_depth'] = vcf['ref'] + vcf['alt']
vcf['var_allele_frac'] = vcf['alt'] / (vcf['alt'] + vcf['ref'])
vcf = vcf[vcf['pass_fail'] == 'PASS']
max_var_allele_frac = 1
min_depth = 8
vcf = vcf[(vcf['var_allele_frac'] <= max_var_allele_frac) & (vcf['total_depth'] >= min_depth)]
vcf = vcf[vcf['mut_type'] != 'SNV']
vcf.to_csv(fh +'.indels.no_max_var_allele.tsv',sep='\t')


