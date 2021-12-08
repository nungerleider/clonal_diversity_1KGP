#!/bin/bash

### Run the variant calling pipeline through this script ###

SRR=$1  # Provide accession number (i.e. SRR123545)
fastq_dir=$2  # Where to store the downloaded fastq file(s)
log_dir="${SRR}/logs"
mkdir -p $log_dir

# Double check the internet connection. This way if this script is run within a loop, 
# if the internet cuts out, it will pause the loop rather than try to download all files (and fail)

while true; do
    ping www.google.com -c1
    if [ $? -eq 0 ];then
        break
    fi
    echo "sleeping"
    sleep 30
done

# Make sure to have the sra-toolkit installed and binaries located in a $PATH directory 
# This script was set up to download paired end files from the NCBI SRA. If files of interest are
# single end, remove the "--split-files" flag in the "fasterq-dump" command below

echo "Downloading $SRR ..."
fasterq-dump -O $fastq_dir --split-files $SRR 2>${log_dir}/${SRR}.fasterqdump_log.err 1>${log_dir}/${SRR}.fasterqdump_log.sdout
echo "Done"

# If the download succeeded, variable "$?" will == "0". Otherwise retry the download

attempts=1
while [[ $? -ne 0 ]]; do
    attempts=$(( attempts + 1 ))

    if [[ $attempts -gt 10 ]];then
        echo "Download failed 10 times - check internet connection and/or make sure accession number is valid"
        exit 1
    fi

    echo "retrying download. Attempt number ${attempts}"
    fasterq-dump -O $fastq_dir --split-files $SRR 2>${log_dir}/${SRR}.fasterqdump_log.err 1>${log_dir}/${SRR}.fasterqdump_log.sdout

done

# After download is successful, run the fastq files through strelka (change strelka path below to local path)
echo "Running alignment and strelka $SRR ..."
bash /home/mcs/Downloads/strelka-2.9.10.centos6_x86_64/bin/strelka.sh ${fastq_dir}/${SRR}_1.fastq 2> ${log_dir}/${SRR}.bwa_and_strelka.error 1> ${log_dir}/${SRR}.bwa_and_strelka.sdout
echo "Done"

# Strelka will produce a gzipped VCF file. This file needs to be unzipped prior to the 'filter_vcf.py' command
echo "gunzipping $SRR"
gunzip ~/Desktop/${SRR}/results/variants/variants.vcf.gz
echo "Done"

# Run all variants through our initial filtration step (filter_vcf.py)
echo "filtering below 0.1 VAF and >8 depth $SRR ..."
python /home/mcs/Desktop/filter_vcf.py ~/Desktop/${SRR}/results/variants/variants.vcf 2>${log_dir}/${SRR}.filter_vcf.err 1>${log_dir}/${SRR}.filter_vcf.sderr
echo "Done"


echo "filtering dbsnps ${SRR} ..."
bash /home/mcs/Desktop/bcftool_vcf_conversion.sh ~/Desktop/${SRR}/results/variants/variants.vcf 2>${log_dir}/${SRR}.filter_dbsnp.err 1>${log_dir}/${SRR}.filter_dbsnp.sdout
echo "Done"



