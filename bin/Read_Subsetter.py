#!/usr/bin/env python3

import os
import sys
import re
import pandas as pd
import argparse
from pathlib import Path
import concurrent.futures
import subprocess
from glob import glob

def count_bases(reads):
    if "," in reads:
        cmd = f"reformat.sh in={reads.replace(',', ' in2=')} out=/dev/null"
    else:
        cmd = f"reformat.sh in={reads} out=/dev/null"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    match = re.search(r"(\d+) bases", result.stderr.decode())
    return int(match.group(1)) if match else 0

def subset_bases(args, row):
    index, row = row 
    reads = row['Reads'].split(',')
    if len(reads) > 1:
        cmd = (f"reformat.sh in={reads[0]} in2={reads[1]} "
               f"out={args.output}/{row['Dataset']}_GenomeSubset_1.fq.gz "
               f"out2={args.output}/{row['Dataset']}_GenomeSubset_2.fq.gz "
               f"samplereadstarget={int(row['ToSample'])}")
    else:
        cmd = (f"reformat.sh in={reads[0]} out={args.output}/{row['Dataset']}_GenomeSubset.fq.gz "
               f"samplereadstarget={int(row['ToSample'])}")
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genomesize', type=int, required=True, help="Genome size estimate in bp")
    parser.add_argument('-c', '--coverage', type=float, default=10, help="Desired coverage")
    parser.add_argument('-r', '--reads', type=str, required=True, help="Directory containing read files")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for subset reads")
    parser.add_argument('-e', '--extension', type=str, default="fastq.gz", help="Read file extension")
    parser.add_argument('--forward', type=str, default="_1.fastq.gz", help="Suffix for forward reads")
    parser.add_argument('--reverse', type=str, default="_2.fastq.gz", help="Suffix for reverse reads")
    return parser.parse_args()

def setup_log(args):
    log_file = f"{args.output}/Subset_Log.txt"
    with open(log_file, 'w') as log:
        log.write("SNPRS Read Subsetter Log\n")
        log.write("-------------------------------------------------------\n\n")
        log.write(f"\t- Genome size estimate: {args.genomesize} bp\n")
        log.write(f"\t- Desired coverage: {args.coverage}X\n")
        log.write(f"\t- Read directory: {os.path.abspath(args.reads)}\n")
        log.write(f"\t- Output directory: {os.path.abspath(args.output)}\n")
        log.write("\n-------------------------------------------------------\n\n")
    return log_file

def validate_directories(read_dir, output_dir):
    if not Path(read_dir).is_dir():
        sys.exit("Path to reads does not exist.")
    if Path(output_dir).is_dir() and any(Path(output_dir).iterdir()):
        sys.exit("Output path not empty.")
    elif not Path(output_dir).is_dir():
        os.makedirs(output_dir)

args = setup_args()
read_ext = args.extension.lstrip('.')
validate_directories(args.reads, args.output)
log_file = setup_log(args)

# Find taxa directories
taxa_dirs = sorted({Path(f).parent for f in Path(args.reads).rglob(f"*.{read_ext}")})
if not taxa_dirs:
    with open(log_file, 'a') as log:
        log.write(f"ERROR: No files found with extension .{read_ext}\n")
    sys.exit(f"No files found with extension .{read_ext}")

subset_depth = int((args.coverage * args.genomesize) / len(taxa_dirs))

with open(log_file, 'a') as log:
    log.write(f"\t- Found {len(taxa_dirs)} taxa groups\n")
    log.write(f"\t- Requested subset depth: {subset_depth} bp per group\n")

# Prepare DataFrame
df = pd.DataFrame(columns=['Taxon', 'Dataset', 'Reads'])

for tax_dir in taxa_dirs:
    taxon_ID = tax_dir.name
    files = sorted(glob(f"{tax_dir}/*.{read_ext}"))
    left_files = [f for f in files if args.forward in f]
    right_files = [f for f in files if args.reverse in f]

    paired_files = {Path(f).stem for f in left_files} & {Path(f).stem for f in right_files}
    left_pairs = [f"{tax_dir}/{stem}{args.forward}" for stem in paired_files]
    right_pairs = [f"{tax_dir}/{stem}{args.reverse}" for stem in paired_files]
    paired_files = sorted(left_pairs + right_pairs)
    single_end_files = [f for f in files if f not in paired_files]

    for left, right in zip(left_pairs, right_pairs):
        df = df.append({'Taxon': taxon_ID, 'Dataset': Path(left).stem, 'Reads': f"{left},{right}"}, ignore_index=True)
    for s in single_end_files:
        df = df.append({'Taxon': taxon_ID, 'Dataset': Path(s).stem, 'Reads': s}, ignore_index=True)

df['Basecount'] = df['Reads'].apply(count_bases)
df = df.sort_values(by=['Taxon']).reset_index(drop=True)
df['ToSample'] = 0

for taxa in df['Taxon'].unique():
    taxa_df = df[df['Taxon'] == taxa]
    total_bases = taxa_df['Basecount'].sum()
    if total_bases < subset_depth:
        with open(log_file, 'a') as log:
            log.write(f"\t- NOTE: {taxa} samples have only {total_bases} bp. All reads will be used.\n")
        df.loc[df['Taxon'] == taxa, 'ToSample'] = taxa_df['Basecount']
    else:
        per_dataset = subset_depth / len(taxa_df)
        low_datasets = taxa_df[taxa_df['Basecount'] < per_dataset]['Dataset'].tolist()
        high_datasets = taxa_df[taxa_df['Basecount'] >= per_dataset]['Dataset'].tolist()
        df.loc[df['Dataset'].isin(low_datasets), 'ToSample'] = df[df['Dataset'].isin(low_datasets)]['Basecount']
        taxa_df = taxa_df[taxa_df['Dataset'].isin(high_datasets)]
        per_dataset = int((subset_depth - df[df['Dataset'].isin(low_datasets)]['Basecount'].sum()) / len(taxa_df))
        df.loc[df['Dataset'].isin(high_datasets), 'ToSample'] = per_dataset

with concurrent.futures.ThreadPoolExecutor() as executor:
    list(executor.map(lambda row: subset_bases(args, row), df.iterrows()))

df.to_csv(f"{args.output}/Subset_Scheme.csv", index=False)

with open(log_file, 'a') as log:
    log.write("\n-------------------------------------------------------\n")
    log.write(f"\t- Finished!\n")
    log.write(f"\t- Subset scheme written to {args.output}/Subset_Scheme.csv\n")
    log.write("-------------------------------------------------------\n")

print(args.output)
