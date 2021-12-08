import requests
import pandas as pd
import glob
import sys

'''Given a final, filtered vcf file, use the ensembl server to query predicted impact (i.e. 'missense', 'intronic', etc)'''

## come back to this ##
files = glob.glob('all_ERR_folders/*/results/variants/*try10*') 


# Some of this code was borrowed from the ensembl REST API documentation
# https://rest.ensembl.org/documentation/info/vep_region_post
def get_variant(handle):

    server = "https://rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    
    with open(handle) as infile, open(handle + '.with_consequence_annotation.tsv', 'w') as outfile:
        next(infile)  # skip header line
        for line in infile:
            line2 = line.split()[1:9]

            # Ensembl annotation has no 'chr' chromosomal prefex (i.e. instead of 'chr1', they use '1')
            line2[0] = line2[0].replace('chr','')

            # The API requires these fields be marked with '.' if they aren't to be included in query
            line2[-3], line2[-2], line2[-1] = '.', '.', '.'

            # Convert list to string for POST submission 
            line2 =' '.join(line2)
            
            # Format string using REST API instructions
            data = '{"variants": ["%s"]}' %line2

            r = requests.post(server + ext, headers=headers, data=data)       
            decoded = r.json()

            # If postulated variant is located within an intron of a transcript and an exon of a different transcript, 
            # classify the variant according to the highest potential impact it could make (coding sequence > non-coding, etc)
            consequence = decoded[0]['most_severe_consequence']
            line.strip('\n') +=f'\t{consequence}'
            outfile.write(f'{line}\n')



def cleanup(handle):

    '''The REST API returned entries spanning two lines - this function merged the two lines into one 
    necessary for downstream processing'''

    infile = open(handle)
    d = {}
    x = 0
    for i in infile:
        i = i.strip('\n')
        if x % 2 == 0:
            ky = i
        else:
            d[ky] = i
        x+=1
    outfile = open(f'{handle}2.tsv', 'w')
    for key,val in d.items():
        outfile.write(f'{key}\t{val}\n')
    outfile.close()