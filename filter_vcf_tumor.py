import pandas as pd
import sys
import numpy as np

''' Variants were called using Strelka v.2. Output 'variants.vcf' files were processed and filtered 
by this script'''


def get_depth(info3):
    '''Calculate the depth of coverage at each called SNV, extracting alternate and ref depth
     using the third 'info' column'''
    
    alleles = [int(i) for i in info3.split(':')[5].split(',')]
    return np.sum(alleles)


def get_variant_fraction(info3):
    '''Calculate the fraction of alleles that have an alternate base at each called SNV, 
    using the third 'info' column'''
    
    alleles = [int(i) for i in info3.split(':')[5].split(',')]
    return alleles[1] / np.sum([alleles[1], alleles[0]])

def filter_single_sub(vcf):
    '''With the entire vcf file as input, require that each entry has a single base substitution
    (filter out any indels or consecutive substitutions of bases/ substitutions > 1 nucleotide)'''

    status = []
    for ref, alt in zip(vcf['ref'], vcf['alt']):
        if len(ref) == len(alt) == 1:
            status.append('snv')
        else:
            status.append('other')                                  

    vcf['variant_status']  = status
    return vcf[vcf['variant_status'] == 'snv']


def collapse(vcf):
    '''The strand of base substitutions cannot be determined using standard WGS. Each variant classified according to the following rules:
    if the reference p
     
    C       T   
    |   >   |  was always called 'C>T'
    G       A

    T       A
    |   >   |  was always called 'T>A'
    A       T

    etc. Starting base was collapsed to C or T'''

    snvs = []
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for ref, alt in zip(vcf['ref'], vcf['alt']):
        if ref not in ['C', 'T']:
            ref = complement[ref]
            alt = complement[alt]
        snv = f'{ref}>{alt}'
        snvs.append(snv)
    vcf['snv'] = snvs
    return vcf
            


vcf_path = sys.argv[1] # variants.vcf output from strelka v.2

max_var_allele_frac = .1  # We required an alt:(alt+ref) ratio of <= 0.1 (to avoid SNPs being interpreted as de novo mutations)
min_depth = 20  # There needed to be sufficient evidence through aligned reads to qualify a SNV


vcf = pd.read_table(vcf_path, header=None, comment='#')
vcf.columns = ['chromosome', 'position', '_','ref', 'alt', '_','pass_fail', 'info1', 'info2', 'info3']  # Strelka output has 3 info columns
vcf = vcf[vcf['pass_fail'] == 'PASS']
vcf = filter_single_sub(vcf)
vcf['depth'] = vcf['info3'].map(get_depth)
vcf['var_allele_frac'] = vcf['info3'].map(get_variant_fraction)
vcf = vcf[vcf['depth'] >= min_depth]
vcf = collapse(vcf)

vcf.to_csv(vcf_path + f'.filtered.vaf{max_var_allele_frac}.tumor_filter_pipeline.tsv', sep='\t')
