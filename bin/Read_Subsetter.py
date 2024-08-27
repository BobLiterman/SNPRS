#!/usr/bin/env python3
#
# Arguments: 
# -g/--genomesize: A group genome size estimate in bp (e.g. 3500000000 for a primate dataset)
# -c/--coverage: When sampling reads, sample genome to X coverage
# -r/--reads: Directory containing folders of reads
# -o/--output: Directory to output subset reads
# -e/--extension (OPTIONAL): Read file extension (Default: fastq.gz)
# --forward (OPTIONAL): Suffix for forward reads (Default: _1.fastq.gz)
# --reverse (OPTIONAL): Suffix for reverse reads (Default: _2.fastq.gz)
#
# Output: For an analysis with x taxa groups, a genome size estimate of g, and a final desired coverage of c, this script will subset reads from each taxon group down to (c*g)/x bases

import os
from os import path
import sys
from glob import glob
import subprocess
import pandas as pd
import re
import argparse
from pathlib import Path
import concurrent.futures

def countBases(index_row):
    index, row = index_row  # Unpack the tuple into index and row
    if "," in row['Reads']:
        reads = row['Reads'].split(",")
        cmd = "reformat.sh in="+reads[0]+" in2="+reads[1]+" out=/dev/null"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bases = int(re.search(r"(\d+) bases", result.stderr.decode()).group(1))
    else:
        cmd = "reformat.sh in="+row['Reads']+" out=/dev/null"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bases = int(re.search(r"(\d+) bases", result.stderr.decode()).group(1))
    return bases

def subsetBases(index_row):
    index, row = index_row  # Unpack the tuple into index and row
    if "," in row['Reads']:
        reads = row['Reads'].split(",")
        cmd = "reformat.sh in="+reads[0]+" in2="+reads[1]+" out="+subset_output_dir+"/"+row['Dataset']+"_GenomeSubset_1.fq.gz out2="+subset_output_dir+"/"+row['Dataset']+"_GenomeSubset_2.fq.gz samplereadstarget="+str(int(row['ToSample']))
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = "reformat.sh in="+row['Reads']+" out="+subset_output_dir+"/"+row['Dataset']+"_GenomeSubset.fq.gz samplereadstarget="+str(int(row['ToSample']))
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
# Process arguments
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-g','--genomesize',action='store',required=True,type=int)
my_parser.add_argument('-c','--coverage',action='store',default=10,type=float)
my_parser.add_argument('-r','--reads',action='store',required=True,type=str)
my_parser.add_argument('-o','--output',action='store',required=True,type=str)
my_parser.add_argument('-e','--extension',action='store',default="fastq.gz",type=str)
my_parser.add_argument('--forward',action='store',default="_1.fastq.gz",type=str)
my_parser.add_argument('--reverse',action='store',default="_2.fastq.gz",type=str)
args = my_parser.parse_args()

genomeSize = int(args.genomesize)
genomeCov = float(args.coverage)
readPath = os.path.normpath(str(args.reads))
outPath = os.path.normpath(str(args.output))
if args.extension[0] == ".":
    read_ext = args.extension[1:]
else:
    read_ext = args.extension
    
# Check/set path where trimmed read files live
if not path.isdir(readPath):
    sys.exit("Path to reads does not exist.")
trim_read_dir = readPath

# Check if output folder exists and is empty
if path.isdir(outPath):
    if len(os.listdir(outPath)) > 0:
        sys.exit("Output path not empty.")
else:
    os.mkdir(outPath)
subset_output_dir = outPath
subset_log_file = outPath+"/Subset_Log.txt"

with open(subset_log_file, 'w+') as log:
    log.write("SNPRS Read Subsetter Log\n")
    log.write("-------------------------------------------------------\n\n")
    log.write(f"\t- Genome size estimate: {genomeSize} bp\n")
    log.write(f"\t- Desired coverage: {genomeCov}X\n")
    log.write(f"\t- Read directory: {os.path.abspath(trim_read_dir)}\n")
    log.write(f"\t- Output directory: {os.path.abspath(subset_output_dir)}\n")
    log.write("\n-------------------------------------------------------\n\n")

# Set taxa folders to only those containing files ending in args.extension
trim_read_tax_dirs = sorted(list(set([os.path.dirname(f) for f in list(Path(trim_read_dir).rglob("*."+read_ext))])))

# Ensure there is at least one directory with read files
if len(trim_read_tax_dirs) == 0:
    with open(subset_log_file, 'a+') as log:
        log.write("ERROR: No files found in read directory with extension ."+read_ext)
    sys.exit("No files found in read directory with extension ."+read_ext)
else:
    # Calculate subset depth
    subsetDepth = int((genomeCov*genomeSize)/len(trim_read_tax_dirs))
    with open(subset_log_file, 'a+') as log:
        log.write(f"\t- Found {len(trim_read_tax_dirs)} taxa groups with read files\n")
        log.write("\t- Based on a genome size estimate of " + str(genomeSize) + " bp, and with " + str(len(trim_read_tax_dirs)) + " groups, the requested subset depth at " + str(genomeCov) + "X coverage is " + str(subsetDepth) + " bp per group\n")
        log.write("\n-------------------------------------------------------\n\n")

#Initialize Pandas DF to get base counts and lists of paired and single-end files
df = pd.DataFrame(columns=['Taxon','Dataset','Reads'],dtype=object)

#For each taxa directory...
for tax_dir in trim_read_tax_dirs:
    
    # Get taxon ID
    taxon_ID = path.basename(tax_dir)

    # Find read files 
    files = sorted(glob(tax_dir+"/*."+read_ext))

    # Get PE data
    left_files = [s for s in files if args.forward in s]
    right_files = [s for s in files if args.reverse in s]

    # Reset file names and filter out single-end files
    left_files = [x.replace(args.forward, '') for x in left_files]
    right_files = [x.replace(args.reverse, '') for x in right_files]
    paired_files = list(set(left_files).intersection(right_files))

    #Reset file names and filter out single-end files
    left_pairs = []
    right_pairs = []
    for pair in paired_files:
        left_pairs.append(pair+args.forward)
        right_pairs.append(pair+args.reverse)
    paired_files = sorted(left_pairs + right_pairs)
    single_end = [x for x in files if x not in paired_files]
    
    # Add read data to dataframe
    for i in range(0,len(paired_files),2):
        df = df.append({'Taxon':taxon_ID,
                        'Dataset':path.basename(paired_files[i]).replace(args.forward,''),
                        'Reads':paired_files[i]+","+paired_files[i+1]},ignore_index=True)
    for s in single_end:
        df = df.append({'Taxon':taxon_ID,
                        'Dataset':path.basename(s).replace("."+read_ext,''),
                        'Reads':s},ignore_index=True)

# Add a BaseCount column to df by sending parallel tasks to count_bases
with concurrent.futures.ThreadPoolExecutor() as executor:
    df['Basecount'] = list(executor.map(countBases, df.iterrows()))

df = df.reset_index(drop=True)
df["Basecount"] = pd.to_numeric(df["Basecount"])
df = df.sort_values(by=['Taxon'])
df["ToSample"] = 0
taxa_list = sorted(df.Taxon.unique())

for taxa in taxa_list:
    
    taxaDF = df[df['Taxon']==taxa]
    dataset_count = taxaDF.shape[0] # Count samples for taxa
    per_dataset = subsetDepth/dataset_count # Calculate bases needed per sample if evenly distributed
    s = pd.Series(taxaDF["Basecount"]) # Grab basecount column
    
    # Check if taxa has enough total coverage. If not, use all reads and warn user.
    if(taxaDF["Basecount"].sum() < subsetDepth):
        with open(subset_log_file, 'a+') as log:
            log.write("\t- NOTE: " + taxa + " samples only have a total of " + str(taxaDF["Basecount"].sum()) + " bp. All " + taxa + " reads will be used.")
            log.write("\n\t- WARNING: Having under the required read threshold may have negative consequences on composite genome assembly and site calling.\n")
        df.loc[df["Taxon"] == taxa, ["ToSample"]] = list(df[df["Taxon"] == taxa].Basecount)

    # If there is enough total data, check sample values to optimize sampling
    else:
        if((s < per_dataset).any()):
            subset_countdown =  subsetDepth # Create countdown variable

            while((s < per_dataset).any()):

                low_datasets = list(taxaDF[taxaDF.Basecount<per_dataset]["Dataset"]) # Find datasets with too few bases
                high_datasets = list(taxaDF[taxaDF.Basecount>=per_dataset]["Dataset"]) # Find datasets with enough bases (currently)

                sum_low = taxaDF[taxaDF.Basecount<per_dataset].Basecount.sum() # Sum the bases of the low datasets
                
                with open(subset_log_file, 'a+') as log:
                    log.write("\t- Not all " + taxa + " samples have enough for even read sampling (" + ",".join(low_datasets) + "). Will take all reads from data poor samples and compensate from larger samples...\n")
                df.loc[df.Dataset.isin(low_datasets), ["ToSample"]] = list(df[df.Dataset.isin(low_datasets)].Basecount) # Set low datasets to use all reads

                taxaDF = taxaDF[taxaDF.Dataset.isin(high_datasets)] # Reset taxa DF to remove low samples
                dataset_count = taxaDF.shape[0] # Count remaining datasets
                subset_countdown = subset_countdown - sum_low

                per_dataset = int(subset_countdown/dataset_count) # Reset per_dataset

                s = pd.Series(taxaDF.Basecount) # Recount basecounts for while loop

            df.loc[df.Dataset.isin(high_datasets), ["ToSample"]] = per_dataset
        else:
            df.loc[df.Taxon == taxa, ["ToSample"]] = per_dataset
    
with concurrent.futures.ThreadPoolExecutor() as executor:
    list(executor.map(subsetBases, df.iterrows()))

df.to_csv(subset_output_dir+"/Subset_Scheme.csv",index=False)
with open(subset_log_file, 'a+') as log:
    log.write("\n-------------------------------------------------------\n\n")
    log.write("\n\n\t- Subset scheme written to " + subset_output_dir + "/Subset_Scheme.csv\n")
    log.write("\n-------------------------------------------------------\n\n")
    
print(subset_output_dir)