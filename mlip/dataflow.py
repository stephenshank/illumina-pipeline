import os
import sys
import subprocess
import csv
import gzip
import shutil
import re
from pathlib import Path
import argparse
from collections import defaultdict
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import yaml


with open('config.yml') as config_file:
    config = yaml.safe_load(config_file)


def load_reference_dictionary(subtype):
    reference_dictionary = {}
    with open('references.tsv', 'r', newline='') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        for row in reader:
            if row['subtype'] == subtype:
                reference_dictionary[row['segment_key']] = row
    return reference_dictionary


def check_duplicates(lines):
    seen = {}
    saw_dupes = False
    for line_num, line in enumerate(lines, 1):
        id_ = line.strip()
        if id_ in seen:
            saw_dupes = True
            print(f"Duplicate ID '{id_}' found on lines {seen[id_]} and {line_num}")
        else:
            seen[id_] = line_num
    if saw_dupes:
        print('You have duplicate sequencing IDs! MLIP cannot proceed until this is rectified.')
        sys.exit(1)


def preprocess(id_filepath, seq_key='Seq'):
    with open(id_filepath, 'r') as f:
        lines = f.read().splitlines()
    check_duplicates(lines)

    sorted_seq_ids = sorted(lines, key=lambda x: x.lower())
    seq_key_pattern = re.compile(rf'_{seq_key}(\d+)')
    key_hash = {}


    fieldnames = ['SequencingId', 'SampleId', 'Replicate']
    os.makedirs('data', exist_ok=True)
    f = open('data/metadata.tsv', 'w', newline='')
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    # process keys
    for seq_id in sorted_seq_ids:
        match = seq_key_pattern.search(seq_id)
        found_match = False
        if match:
            sample_id = seq_id.replace(match.group(0), '').lower()
            found_match = True
        else:
            sample_id = ''
        writer.writerow({
            'SequencingId': seq_id,
            'SampleId': sample_id,
            'Replicate': ''
        })

    f.close()
    print('Metadata spreadsheet written to data/metadata.tsv.')
    print('Please edit, then run the flow step.')
    return


def preprocess_cli(args):
    if args.key:
        preprocess(args.file, args.key)
    else:
        preprocess(args.file)


def tokenize(identifier):
    return identifier.lower().replace('_', '').replace('-', '')


def is_forward_fastq_path(fastq_path):
    return fastq_path.split('_')[-2] == 'R1'


def convert_forward_to_reverse(fastq_path):
    dirname, filename = os.path.split(fastq_path)
    new_filename = filename[::-1].replace("1R", "2R", 1)[::-1]  # Reverse swap
    return os.path.join(dirname, new_filename)


def fastq_is_low_coverage(filepath, min_reads=50):
    try:
        with gzip.open(filepath, 'rt') as f:
            lines = 0
            for line in f:
                lines += 1
                if lines > min_reads * 4:
                    return False
            return (lines // 4) < min_reads
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return True


def flow(args):
    with open('data/metadata.tsv', 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    fastq_paths = [
        str(fastq_path)
        for fastq_path in Path(config['data_root_directory']).expanduser().rglob('*.fastq.gz')
        if is_forward_fastq_path(str(fastq_path))
    ]
    fastq_tokens = [tokenize(os.path.basename(fastq_path)) for fastq_path in fastq_paths]
    fastq_hash = {
        token: {
            'path': path,
            'reverse': convert_forward_to_reverse(path),
            'seen': False
        }
        for token, path in zip(fastq_tokens, fastq_paths)
    }
    counter = Counter()
    empties = []
    for row in rows:
        sample_id = row['SampleId']
        counter[sample_id] += 1
        sequencing_token = tokenize(row['SequencingId'])

        token_pattern = re.compile(rf'^{sequencing_token}(?!\d)')
        found_match = False
        for fastq_token in fastq_tokens:
            if token_pattern.search(fastq_token):
                if fastq_hash[fastq_token]['seen'] == True:
                    print('FATAL ERROR! Tokens clashed. Please contact Stephen.')
                    sys.exit(1)
                found_match = True
                fastq_hash[fastq_token]['seen'] = True

                directory_path = os.path.join('data', sample_id, f'sequencing-{counter[sample_id]}')
                if not args.dry_run:
                    os.makedirs(directory_path, exist_ok=True)

                old_forward_path = fastq_hash[fastq_token]['path']
                new_forward_path = os.path.join(directory_path, 'forward.fastq.gz')
                old_reverse_path = fastq_hash[fastq_token]['reverse']
                new_reverse_path = os.path.join(directory_path, 'reverse.fastq.gz')

                # Check if fastq files are empty; if so, alert user, log sample id, and skip moving
                if fastq_is_low_coverage(old_forward_path) or fastq_is_low_coverage(old_reverse_path):
                    print(f"Alert: Fastq file(s) for sample {sample_id} are empty. Skipping move.")
                    empties.append(sample_id)
                    continue

                if not args.dry_run:
                    shutil.copy(old_forward_path, new_forward_path)
                    shutil.copy(old_reverse_path, new_reverse_path)

                print(f"{sample_id}: sequencing experiment {counter[sample_id]}, replicate {row['Replicate']}")
                print(f'\tForward: {old_forward_path}')
                print(f'\tMoving to {new_forward_path}\n')
                print(f'\tReverse: {old_reverse_path}')
                print(f'\tMoving to {new_reverse_path}\n\n')
        if not found_match:
            print(f'Warning! Could not find a match for {row['SequencingId']}!')

    with open("data/empty.txt", "w") as empty_file:
        for empty in set(empties):
            empty_file.write(empty+ "\n")


def flow_cli(args):
    flow(args)


def command_line_interface():
    parser = argparse.ArgumentParser(description="MLIP dataflow")
    subparsers = parser.add_subparsers(dest="command", required=True)

    pre_parser = subparsers.add_parser("preprocess", help="Convert Sequencing ID list to metadata spreadsheet.")
    pre_parser.add_argument("-f", "--file", required=True, type=str, help="Path to the ID list")
    pre_parser.add_argument("-k", "--key", required=False, type=str, default="Seq", help="Key to denote sequencing experiment")
    pre_parser.set_defaults(func=preprocess_cli)

    flow_parser = subparsers.add_parser("flow", help="Move data based on metadata spreadsheet.")
    flow_parser.add_argument("-n", "--dry-run", required=False, action="store_true", default=False, help="Show how FASTQs will move, but don't actually move them")
    flow_parser.set_defaults(func=flow_cli)

    args = parser.parse_args()
    args.func(args)


def load_metadata_dictionary():
    f = open('data/metadata.tsv', 'r')
    reader = csv.DictReader(f, delimiter='\t')
    md_dict = defaultdict(lambda: defaultdict(list))
    counter = Counter()
    for row in reader:
        sample_id = row['SampleId']
        replicate = row['Replicate']
        counter[sample_id] += 1
        md_dict[sample_id][replicate].append(counter[sample_id])
    f.close()
    return md_dict


def samples_to_analyze():
    md = load_metadata_dictionary()
    with open('data/empty.txt') as f:
        empties = [line.strip() for line in f.readlines()]
    for empty in empties:
        del md[empty]
    return md.keys()


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
    return varscan_df.loc[:, ['#CHROM', 'POS', 'REF', 'ALT', 'Frequency']]


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


consensus_coverage = config['consensus_minimum_coverage']
variant_coverage = config['variants_minimum_coverage']
coverage_bucket_labels = [
    '0x',
    f'1-{consensus_coverage}x',
    f'{consensus_coverage}-{variant_coverage}x',
    f'{variant_coverage}-1000x',
    '1000x+'
]


def assign_coverage_bucket(coverage):
    if coverage == 0:
        return coverage_bucket_labels[0]
    elif 1 <= coverage <= consensus_coverage:
        return coverage_bucket_labels[1]
    elif consensus_coverage <= coverage <= variant_coverage:
        return coverage_bucket_labels[2]
    elif variant_coverage <= coverage <= 1000:
        return coverage_bucket_labels[3]
    else:
        return coverage_bucket_labels[4]


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
        'coverage_bucket': coverage_bucket_labels
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
    ).reset_index()[['segment'] + coverage_bucket_labels]
   
    desired_structure.columns.name = None
    return desired_structure

def compute_coverage_categories_io(input_coverage, output_summary):
    coverage_df = pd.read_csv(input_coverage, sep='\t')
    coverage_summary_df = compute_coverage_categories(coverage_df)
    coverage_summary_df.to_csv(output_summary, sep='\t', index=False)


def define_coding_regions(gtf_file):
    with open(gtf_file, "r") as gtf:
        coding_regions = {}

        for line in gtf:
            if line.strip("\n") != "":      # ignore blank lines (otherwise throws an index error)
                sequence_name = line.split("\t")[0]
                annotation_type = line.split("\t")[2]
                start = int(line.split("\t")[3]) - 1  # adding the -1 here for 0 indexing
                stop = int(line.split("\t")[4]) - 1    # adding the -1 here for 0 indexing
                gene_name = line.split("\t")[8]
                gene_name = gene_name.split(";")[0]
                gene_name = gene_name.replace("gene_id ","")
                gene_name = gene_name.replace("\"","")

                if annotation_type.lower() == "cds":
                    if sequence_name not in coding_regions:
                        coding_regions[sequence_name] = {}
                        coding_regions[sequence_name][gene_name] = [start, stop]
                    elif sequence_name in coding_regions and gene_name not in coding_regions[sequence_name]:
                        coding_regions[sequence_name][gene_name] = [start, stop]
                    elif gene_name in coding_regions[sequence_name]:
                        coding_regions[sequence_name][gene_name].extend([start, stop])

            # sort coding region coordinates so that they are always in the correct order
            for sequence_name in coding_regions:
                for gene in coding_regions[sequence_name]:
                    coding_regions[sequence_name][gene] = sorted(coding_regions[sequence_name][gene])

    return(coding_regions)


amino_acid_abbreviations = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "O": "Pyl", "S": "Ser", "U": "Sec", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "B": "Asx", "Z": "Glx", "U": "Sec",
    "X": "Xaa", "J": "Xle", "*": "Stop"
}


def slice_fastas(coding_regions, reference_sequence_path):
    reference_sequence = {}

    for seq in SeqIO.parse(reference_sequence_path, "fasta"):
        sequence_name = str(seq.id)
        sequence = str(seq.seq).lower()
        reference_sequence[sequence_name] = sequence

    transcripts = {}
    stop_codons = ["taa","tag","tga"]
    refseqname = 'A sequence'

    for c in coding_regions:
        for gene in coding_regions[c]:
            transcripts[gene] = ""
            coordinates = coding_regions[c][gene]   # define the coding regions for each gene

            for i in range(0,int(len(coordinates)),2):   # go along the coordinates in chunks of 2 at a time
                sequence_chunk = reference_sequence[c][coordinates[i]:coordinates[i+1]+1]
                transcripts[gene] += sequence_chunk     # append each piece of the transcript together 

    # loop through each transcript to make sure that it begins with a start codon and ends with a stop codon
    for t in transcripts:
        if transcripts[t][0:3] != "atg":
            print("WARNING! " + refseqname + " " + t + "does not contain a start codon!")
        if transcripts[t][-3:] not in stop_codons:
            print("WARNING! " + refseqname + " " + t + " does not contain a stop codon! These are the last 3 nucleotides: " + transcripts[t][-3:])

    return transcripts


def annotate_amino_acid_changes(coding_regions, transcripts, vcf, outfilename):
    with open(vcf, "r") as csvfile:
        with open(outfilename, "w") as outfile:
            to_write = ["sample","gene","reference_position","reference_allele","variant_allele","coding_region_change","synonymous/nonsynonymous","frequency(%)","frequency","\n"]
            to_write2 = "\t".join(to_write)
            outfile.write(to_write2)

        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            # ignore comment lines
            if "##" not in row[0] and "#CHROM" not in row[0]:
                sequence_name = row[0]
                site = int(row[1]) - 1     # to make this 0 indexed
                reference_allele = row[3].lower()
                alternative_allele = row[4].lower()
 
                # pull out the frequency using a string search
                SearchStr = r'.+\:([0-9]{1,2}\.{0,1}[0-9]{0,2}\%)'
                result = re.search(SearchStr,row[9])
                if result:
                    frequency = '\t'.join(result.groups())
                    frequency_decimal = frequency.replace("%","")
                    frequency_decimal = (float(frequency_decimal))/100
                    frequency_decimal = "%.4f" % frequency_decimal    # only include 4 numbers after decimal

                else:
                    frequency = "none reported"

                # figure out whether the SNP lies within a coding region:
                for gene in coding_regions[sequence_name]:
                    coordinates = coding_regions[sequence_name][gene]

                    # go through gene coordinates 2 at a time; this is for genes with multiple regions
                    for i in range(0,int(len(coordinates)),2):
                        if site >= coordinates[i] and site <= coordinates[i+1]:  # if site is within the gene

                            # determine the coding region site, depending on if there are 2 frames or 1
                            if len(coordinates) == 2:
                                cds_site = site - coordinates[i]
                            elif len(coordinates) == 4:
                                cds_site = (coordinates[1] - coordinates[0]) + (site - coordinates[2] + 1)
                            # now determine whether the site is in the 1st, 2nd, or 3rd codon position
                            # if SNP is in 1st position in codon:
                            aa_site = int(cds_site/3)+1

                            if float(cds_site) % 3 == 0:
                                codon = transcripts[gene][cds_site:cds_site+3]
                                variant_codon = alternative_allele + codon[1:3]
                                variant_aa = Seq(variant_codon).translate()
                                ref_codon = reference_allele + codon[1:3]
                                ref_aa = Seq(ref_codon).translate()

                            # if variant is in the middle of the codon:
                            elif float(cds_site - 1) % 3 == 0:
                                codon = transcripts[gene][cds_site-1:cds_site+2]
                                variant_codon = codon[0] + alternative_allele + codon[2]
                                variant_aa = Seq(variant_codon).translate()
                                ref_codon = codon[0] + reference_allele + codon[2]
                                ref_aa = Seq(ref_codon).translate()

                            # if the variant is in the 3rd codon position
                            elif float(cds_site -2) % 3 == 0:
                                codon = transcripts[gene][cds_site-2:cds_site+1]
                                variant_codon = codon[0:2] + alternative_allele
                                variant_aa = Seq(variant_codon).translate()
                                ref_codon = codon[0:2] + reference_allele
                                ref_aa = Seq(ref_codon).translate()


                            # return the amino acid changes, with single letters converted to 3-letter aa codes
                            if ref_aa == variant_aa:
                                syn_nonsyn = "synonymous"
                            elif ref_aa != variant_aa and variant_aa == "*":
                                syn_nonsyn = "stop_gained"
                            elif ref_aa != variant_aa:
                                syn_nonsyn = "nonsynonymous"

                            amino_acid_change = amino_acid_abbreviations[ref_aa] + str(aa_site) + amino_acid_abbreviations[variant_aa]
                            with open(outfilename, "a") as outfile:
                                output = [sequence_name,gene,str(site+1),reference_allele.upper(),alternative_allele.upper(),amino_acid_change,syn_nonsyn,frequency,frequency_decimal,"\n"]
                                output2 = "\t".join(output)
                                outfile.write(output2)


def coverage_summary(input_tsvs, output_tsv):
    dfs = []
    for tsv_path in input_tsvs:
        split_path = tsv_path.split('/')
        sample_id = split_path[1]
        replicate = split_path[2]
        df = pd.read_csv(tsv_path, sep='\t')
        df['sample_id'] = f'{sample_id}-{replicate}'
        dfs.append(df)
    full_df = pd.concat(dfs, ignore_index=True)
    full_df['total_coverage'] = full_df[coverage_bucket_labels[0]] + \
        full_df[coverage_bucket_labels[1]] + \
        full_df[coverage_bucket_labels[2]] + \
        full_df[coverage_bucket_labels[3]] + \
        full_df[coverage_bucket_labels[4]]

    full_df.to_csv(output_tsv, sep='\t', index=False)


def mask_low_coverage(fasta_file, pileup_file, cov_threshold, qual_threshold, output_fasta):
    """
    Masks bases in the consensus FASTA if effective coverage is below threshold.
 
    The pileup file may not contain rows for positions with zero coverage.
    We build a dictionary of effective coverage (counting only bases with quality >= qual_threshold)
    for each reported position. Then, for every base (1-indexed) in the FASTA, if no row exists
    (implying zero coverage) or if the effective depth is below cov_threshold, that base is masked to 'N'.
    """
    # Build effective coverage dict: {segment: {position: effective_depth}}
    eff_cov = defaultdict(dict)
    with open(pileup_file) as pf:
        for line in pf:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])  # positions are 1-indexed
            # Use mpileup format if available; else fallback to 5-column pileup
            if len(fields) >= 6:
                qualities = fields[5]
            elif len(fields) == 5:
                qualities = fields[4]
            else:
                qualities = ""
            # Count only bases with quality >= qual_threshold
            depth = sum(1 for q in qualities if (ord(q) - 33) >= qual_threshold)
            eff_cov[chrom][pos] = depth

    # Process each FASTA sequence
    seq_records = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq = list(str(record.seq))
        segment = record.id
        seg_cov = eff_cov.get(segment, {})  # If no pileup rows, assume zero coverage everywhere
        # Iterate through every base position (1-indexed)
        for i in range(len(seq)):
            pos = i + 1
            depth = seg_cov.get(pos, 0)
            if depth < cov_threshold:
                seq[i] = 'N'
        record.seq = ''.join(seq)
        seq_records.append(record)
    SeqIO.write(seq_records, output_fasta, 'fasta')


def get_bases_and_qualities(raw_bases, qual_str):
    """Parses a pileup raw base string and its quality string, returning paired lists."""
    bases, qualities = [], []
    i, q_index = 0, 0
    while i < len(raw_bases):
        c = raw_bases[i]
        if c == '^':
            i += 2  # skip the read-start marker and its mapping quality
            continue
        if c == '$':
            i += 1
            continue
        if c in '+-':
            i += 1
            num = ''
            while i < len(raw_bases) and raw_bases[i].isdigit():
                num += raw_bases[i]
                i += 1
            if num:
                i += int(num)
            continue
        bases.append(c)
        if q_index < len(qual_str):
            qualities.append(ord(qual_str[q_index]) - 33)
            q_index += 1
        else:
            qualities.append(0)
        i += 1
    return bases, qualities

def parse_pileup(pileup_file, qual_threshold):
    """
    Parses a samtools pileup file and, for each position,
    computes the majority base and effective coverage using only bases with quality >= qual_threshold.

    Returns:
        dict: {segment: {position: (majority_base, effective_coverage)}}
    """
    pileup_data = {}
    with open(pileup_file) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 5:
                continue
            segment = fields[0]
            pos = int(fields[1])
            ref_base = fields[2].upper()
            raw_bases = fields[4]
            qual_str = fields[5] if len(fields) >= 6 else ""

            bases, quals = get_bases_and_qualities(raw_bases, qual_str)
            filtered = []
            for b, q in zip(bases, quals):
                if q >= qual_threshold:
                    filtered.append(ref_base if b in ['.', ','] else b.upper())
            effective_cov = len(filtered)
            majority = Counter(filtered).most_common(1)[0][0] if filtered else 'N'

            if segment not in pileup_data:
                pileup_data[segment] = {}
            pileup_data[segment][pos] = (majority, effective_cov)
    return pileup_data

def check_consensus(fasta_file, pileup_data, coverage_threshold):
    """
    Compares a consensus FASTA with pileup-derived majority bases (from high-quality bases only).
    Flags positions where:
      - A base is present when effective coverage is below threshold (should be 'N').
      - An 'N' is present despite sufficient effective coverage.
      - The consensus base does not match the majority call.
    """
    errors = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        segment = record.id
        seq = record.seq.upper()
        if segment not in pileup_data:
            print(f"Warning: Segment {segment} not found in pileup.")
            continue
        for pos in range(1, len(seq) + 1):
            base = seq[pos - 1]
            majority, cov = pileup_data[segment].get(pos, ('N', 0))
            expected = 'N' if cov < coverage_threshold else majority
            if base != expected:
                if base == 'N' and cov >= coverage_threshold:
                    err_type = "Unexpected N"
                elif base != 'N' and cov < coverage_threshold:
                    err_type = "N expected but not present"
                else:
                    err_type = "Consensus base mismatch"
                errors.append((segment, pos, base, expected, cov, err_type))
    if errors:
        print("Discrepancies found:")
        print(f"{'Segment':<12}{'Pos':<8}{'Consensus':<12}{'Expected':<12}{'Coverage':<10}{'Issue'}")
        for seg, pos, base, exp, cov, issue in errors:
            print(f"{seg:<12}{pos:<8}{base:<12}{exp:<12}{cov:<10}{issue}")
    else:
        print("No discrepancies found.")
    return errors


def check_consensus_io(input_consensus, input_pileup, output_tsv, sample, replicate):
    coverage_threshold = config['consensus_minimum_coverage']
    quality_threshold = config['minimum_quality_score']
    pileup_data = parse_pileup(input_pileup, quality_threshold)
    output = check_consensus(input_consensus, pileup_data, coverage_threshold)
    headers = [
        "Sample",
        "Replicate",
        "Segment",
        "Position",
        "Consensus",
        "Expected",
        "Coverage",
        "Issue"
    ]
    with open(output_tsv, "w", newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        writer.writerows([
            (sample, replicate) + row
            for row in output
        ])


if __name__ == '__main__':
    command_line_interface()
