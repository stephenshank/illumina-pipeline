import os
import csv
import shutil
from pathlib import Path

from Bio import SeqIO
import pandas as pd
import yaml


def load_config():
    with open('config.yml') as config_file:
        config = yaml.safe_load(config_file)
    return config


def situate_basespace_data():
    config = load_config()
    csv_filepaths = list(
        Path(config['data_root_directory']).expanduser().rglob('*.csv')
    )
    with open(csv_filepaths[0]) as f:
        metadata_rows = list(csv.reader(f))
    data_index = [i for i, row in enumerate(metadata_rows) if row[0] == '[Data]'][0]
    raw_keys = [
        row[0]
        for row in metadata_rows[data_index+2: ]
        if row[0] != ''
    ]

    key_hash = {}
    # process keys
    for raw_key in raw_keys:
        if '_Rep1' in raw_key:
            run_key = raw_key.replace('_Rep1', '')
            replicate_key = 'replicate-1'
        elif '_Rep2' in raw_key:
            run_key = raw_key.replace('_Rep2', '')
            replicate_key = 'replicate-2'
        else:
            run_key = raw_key
            replicate_key = 'replicate-1'
        new_value = {
            'original_key': raw_key,
            'forward_count': 0,
            'forward_filepath': '',
            'reverse_count': 0,
            'reverse_filepath': ''
        }
        if run_key in key_hash:
            key_hash[run_key][replicate_key] = new_value
        else:
            key_hash[run_key]= {replicate_key: new_value}

    # match to files
    for run_key in key_hash.keys():
        run_key_token = run_key.lower().replace('_', '').replace('-', '')
        fastq_paths = Path(config['data_root_directory']).expanduser().rglob('*.fastq.gz')
        for fastq_path in fastq_paths:
            fastq_basename = os.path.basename(fastq_path)
            fastq_token = fastq_basename.lower().replace('_', '').replace('-', '')
            if run_key_token in fastq_token:    
                if 'rep1' in fastq_token:
                    replicate_key = 'replicate-1'
                elif 'rep2' in fastq_token:
                    replicate_key = 'replicate-2'
                else:
                    replicate_key = 'replicate-1'
                is_forward = fastq_basename.split('_')[-2] == 'R1'
                direction_key = 'forward' if is_forward else 'reverse'    
                count_key = f'{direction_key}_count'
                key_hash[run_key][replicate_key][count_key] += 1
                key_hash[run_key][replicate_key][f'{direction_key}_filepath'] = fastq_path
                if key_hash[run_key][replicate_key][count_key] > 1:
                    raise ValueError("Multiple FASTQs found for same token... aborting!!")
    
    # clean up
    keys_to_delete = []
    for key, value in key_hash.items():
        if value['replicate-1']['forward_count'] == 0 and value['replicate-1']['reverse_count'] == 0:
            keys_to_delete.append(key)
    for key in keys_to_delete:
        del key_hash[key]
    
    # move or copy
    for sample_key, sample_value in key_hash.items():
        for replicate_key, replicate_value in sample_value.items():
            action = 'mov'
            old_forward_path = replicate_value['forward_filepath']
            old_forward_name = os.path.basename(old_forward_path)
            new_forward_path = os.path.join('data', sample_key, replicate_key, 'forward.fastq.gz')
            print(f'{action}ing {old_forward_name} to {new_forward_path}')
            directory_path = os.path.dirname(new_forward_path)
            os.makedirs(directory_path, exist_ok=True)
            shutil.copy(old_forward_path, new_forward_path)

            old_reverse_path = replicate_value['reverse_filepath']
            old_reverse_name = os.path.basename(old_reverse_path)
            new_reverse_path = os.path.join('data', sample_key, replicate_key, 'reverse.fastq.gz')
            print(f'{action}ing {old_reverse_name} to {new_reverse_path}\n')
            directory_path = os.path.dirname(new_forward_path)
            os.makedirs(directory_path, exist_ok=True)
            shutil.copy(old_reverse_path, new_reverse_path)

    sample_keys_filepath = os.path.join('data', 'samples.txt')
    with open(sample_keys_filepath, 'w') as samples_file:
        samples_file.write('\n'.join(key_hash.keys()))


def genbank_to_gtf(gbk_file, gtf_file):
    with open(gtf_file, 'w') as gtf_out:
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract required fields
                    chrom = record.id
                    source = "GenBank"
                    feature_type = "CDS"
                    start = int(feature.location.start) + 1  # GTF uses 1-based indexing
                    end = int(feature.location.end)  # Already 1-based in GenBank
                    score = "."
                    strand = "+" if feature.location.strand == 1 else "-"
                    frame = str(feature.qualifiers.get('codon_start', ['0'])[0])
                    gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
                    transcript_id = feature.qualifiers.get('protein_id', ['unknown'])[0]

                    # Format the GTF line
                    gtf_line = f"{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id \"{gene_name}\"; transcript_id \"{transcript_id}\";\n"
                    gtf_out.write(gtf_line)


def clean_varscan(varscan_df):
    varscan_df['Frequency'] = [
        float(format_.split(':')[6][:-1]) / 100
        for format_ in varscan_df['Sample1']
    ]
    return varscan_df


def merge_varscan(clean_varscan_dfs):
    variant_dictionary = {}
    for idx, clean_varscan_df in enumerate(clean_varscan_dfs):
        sample_key = 'sample_%d' % (idx+1)
        for _, row in clean_varscan_df.iterrows():
            variant_key = (row['#CHROM'], row['POS'], row['ALT'])
            variant_data = {'Frequency': row['Frequency']}
            if variant_key in variant_dictionary:
                variant_dictionary[variant_key][sample_key] = variant_data
            else:
                variant_dictionary[variant_key] = {sample_key: variant_data}
    
    variant_list = []
    for variant_key, variant_value in variant_dictionary.items():
        for idx in range(len(clean_varscan_dfs)):
            sample_key = 'sample_%d' % (idx+1)
            if not sample_key in variant_value:
                variant_value[sample_key] = {'Frequency': 0}
        
        flattened = {}
        for idx, sample_attributes in enumerate(variant_value.values()):
            for attribute_key, attribute_value in sample_attributes.items():
                flattened[f"{attribute_key}_{idx+1}"] = attribute_value
        variant_list.append({
            'segment': variant_key[0],
            'position': variant_key[1],
            'alt': variant_key[2],
            **flattened
        })
    return pd.DataFrame(variant_list)


def merge_varscan_io(input_tsv_filepaths, output_tsv_filepath):
        dfs = [
            pd.read_csv(input_tsv_filepath, sep='\t')
            for input_tsv_filepath in input_tsv_filepaths
        ]
        merged_df = merge_varscan(dfs)
        merged_df.to_csv(output_tsv_filepath, sep='\t')


def assign_coverage_bucket(coverage):
    if coverage == 0:
        return '0x'
    elif 1 <= coverage <= 100:
        return '1-100x'
    elif 101 <= coverage <= 1000:
        return '100-1000x'
    else:
        return '1000x+'


def compute_coverage_categories(df):
    # Calculate the number of sites in each range
    df['num_sites'] = df['end'] - df['start']
    
    # Assign each coverage value to a bucket
    df['coverage_bucket'] = df['coverage'].apply(assign_coverage_bucket)
    
    # Group by segment and coverage bucket and sum the number of sites
    bucket_summary = df.groupby(
        ['segment', 'coverage_bucket']
    )['num_sites'].sum().reset_index(name='count')
    
    # Create all possible combinations of segments and coverage buckets
    all_buckets = pd.DataFrame({
        'coverage_bucket': ['0x', '1-100x', '100-1000x', '1000x+']
    })
    segments = df['segment'].unique()
    all_combinations = pd.MultiIndex.from_product(
        [segments, all_buckets['coverage_bucket']],
        names=['segment', 'coverage_bucket']
    ).to_frame(index=False)
    
    # Merge the combinations with the actual bucket summary and fill missing values with 0
    bucket_summary = pd.merge(
        all_combinations,
        bucket_summary,
        on=['segment', 'coverage_bucket'],
        how='left'
    ).fillna(0)
    
    # Pivot to the desired format
    desired_structure = bucket_summary.pivot_table(
        index='segment', 
        columns='coverage_bucket', 
        values='count', 
        fill_value=0
    ).reset_index()
   
    desired_structure.columns.name = None
    return desired_structure

def compute_coverage_categories_io(input_coverage, output_summary):
    coverage_df = pd.read_csv(input_coverage, sep='\t')
    coverage_summary_df = compute_coverage_categories(coverage_df)
    coverage_summary_df.to_csv(output_summary, sep='\t', index=False)


if __name__ == '__main__':
    situate_basespace_data()