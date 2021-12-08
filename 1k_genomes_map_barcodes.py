import pandas as pd

'''Map sample id to associated fastq files to match RNA seq to 
Whole Genome Seq data. The original SRA run table contains a mapping 
of library barcode to sample id for all 1000 genomes sequencing files'''

wgs_map_path = "/media/mcs/easystore2/1000_genomes_wgs/etc/sra_run_table.tsv" 
rna_map_path = '/media/mcs/easystore1/1000_genomes_RNA/spreadsheets/rna_patient_map.tsv'

wgs_map = pd.read_table(wgs_map_path)
rna_map = pd.read_table(rna_map_path)

wgs_dict = dict(zip(wgs_map['Run'], wgs_map['Sample Name']))
rna_dict = dict(zip(rna_map['RNA_barcode'], rna_map['sample_id']))


