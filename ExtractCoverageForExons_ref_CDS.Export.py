###Get the coverage of exons in a gff for a BAM alignment file
###Feels insane that this doesn't exist
###The problem with featureCounts is the counts rather than coverage - it makes it hard to set a threshold.
###Set a threshold for all exons in a given transcript
###Usage: python ExtractCoverageForExons_ref_CDS.Export.py [path_to_input_file.bam] [path_to_input_annotation.gff/gtf] [Metadata delimiter in gff/gtf] [Feature to count] [path_to_output_file.tsv]

import os
from pathlib import Path
import sys
import os
import glob
import csv
import argparse
from icecream import ic
import ntpath
import pysam
import numpy as np
from tqdm import tqdm
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("input_file", help="input_file")
parser.add_argument("input_gff", help="Input gff")
parser.add_argument("delim", help="Delimiter for use in key:value extraction of metadata. '=' in the NCBI ref, ' ' in the Stringtie ref")
parser.add_argument("feature", help="The target feature to extract coverage info for.",default=str("exon"))
parser.add_argument("output_filepath", help="Output filepath",default=str(""))

args = parser.parse_args()

input_file = args.input_file
input_gff = args.input_gff
output_filepath = args.output_filepath
delim = args.delim
feature = args.feature
ic(delim)
samfile = pysam.AlignmentFile(input_file, "rb")
ic(samfile)

gff_data = [line.strip() for line in open(input_gff, 'r')]
#ic(gff_data)
ic(gff_data[2].split('\t'))


output_data = []
for line in tqdm(gff_data):
        if line.startswith("#"):
                continue
        chr, start, stop, geneID, source = line.split('\t')[0], int(line.split('\t')[3]) ,int(line.split('\t')[4]), line.split('\t')[8],line.split('\t')[2]
        if source == 'transcript':
                continue
        elif source == feature:
                metadata = dict()
                for item in geneID.split(';'):
                        data_split = dict()
                        if len(item.split(delim)) > 1:
#                               ic(item.lstrip())
                                data_split[item.lstrip().split(delim)[0]] = item.lstrip().split(delim)[1]
                                metadata.update(data_split)
#               print("-----")
#               ic(metadata)
                for key,value in metadata.items():
#                       ic(key)
#                       ic(value)
                        metadata[key] = value.replace('"',"")
#                       ic(key)
#                       ic(value)

#               ic(metadata)
#               metadata = [elem.replace('"','') for elem in metadata]
#               ic(metadata)
#               ic(metadata_flat)
                try:
                        coverage = samfile.count_coverage(chr,start,stop)
                except:
                        print("Error in this line:")
                        print(line)
                gene_coverage_list = []
                for letter in coverage:
                        gene_coverage_list.append(letter.tolist())
                gene_coverage_combined = zip(gene_coverage_list[0],gene_coverage_list[1],gene_coverage_list[2],gene_coverage_list[3])
                gene_coverage_combined_sum = [x + y + z + a for (x,y,z,a) in gene_coverage_combined]
                exon_cov_avg = sum(gene_coverage_combined_sum) / len(gene_coverage_combined_sum)
                exon_cov_over0 = sum(x > 0 for x in gene_coverage_combined_sum) / len(gene_coverage_combined_sum)
                exon_cov_over1 = sum(x > 1 for x in gene_coverage_combined_sum) / len(gene_coverage_combined_sum)
#               ic(exon_cov_avg,exon_cov_over0,exon_cov_over1)
                data_out = {'chr' : chr, 'start' : start, 'stop' : stop, 'exon_cov_avg' : exon_cov_avg, 'exon_cov_over0' : exon_cov_over0, 'exon_cov_over1' : exon_cov_over1}
#               ic(data_out)
                data_out.update(metadata)
#               output_line = [chr,start,stop,exon_cov_avg,exon_cov_over0,exon_cov_over1]
#               ic(output_line)
#               output_line = output_line + metadata
#               ic(data_out)
#               ic(metadata)
#               output_line = [str(i) for i in output_line]
                output_data.append(data_out)
        else:
                continue
                #ic("Not an exon or a transcript?")
#ic(output_data)
df = pd.DataFrame(output_data)
ic(df)
ic(output_filepath)
output_file_name = ""
if output_filepath == "":
        print("No output filepath detected, using input filename with .bam replaced by suffix.tsv")
        output_file_name = input_file.replace(".bam",".") + os.path.basename(input_gff)[-35:]+ ".tsv"
else:
        output_file_name = output_filepath
ic(output_filepath)
df.to_csv(output_file_name,index=False,sep='\t')

exit()
