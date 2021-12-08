#!/bin/bash

# Intersect Strelka v2 output VCF file with dbSNP vcf file ('chr' prefix
# was added to each chromosome name in the original dbSNP vcf file to harmonize 
# naming scheme with our genome build).

# Take the complement (i.e. the variant sites that don't coincide with known population
# SNPs)

# VCF file must be already sorted by chromosomal coordinate

# sed command to add 'chr': sed '/#/!s/^/chr//' dbSNP.vcf > chr_added_dbSNP.vcf

vcf_input=$1
dbsnp_input=/media/mcs/easystore1/dbsnp_chr_added.vcf.gz
output_dir=${vcf_input.vcf.gz}
vcf_no_dbsnp=${vcf_input%vcf}intersect_no_dbsnp.vcf

# Requires bgzip and bcftools to be installed and located in $PATH
bgzip $vcf_input
bcftools index $vcf_input
bcftools isec  -C -p $output_dir $vcf_input $dbsnp_input
mv ${output_dir}/sites.txt $vcf_no_dbsnp