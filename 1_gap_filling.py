from Bio import SeqIO
import pandas as pd
import os
from utils import *
from Bio.Seq import Seq
import re

def parse_subject_id(subject_id):
    # Use regular expression to match subject_id in different formats
    pattern = r'''
        ^(.*?_)           # Chromosome part (including underscore)
        (\d+)_(\d+)_      # Two position numbers
        (left|right)_flank$  # Complete flank type
    '''
    match = re.match(pattern, subject_id, re.VERBOSE)
    
    if match:
        chrom = match.group(1).rstrip('_')  # Remove possible trailing underscore
        start = int(match.group(2))
        end = int(match.group(3))
        flank_type = f"{match.group(4)}_flank"
        return chrom, start, end, flank_type
    
    # Alternative solution: handle other possible formats
    parts = subject_id.rsplit('_', 4)  # Split 4 times from the right
    try:
        chrom_part = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        flank_type = '_'.join(parts[3:])
        return chrom_part, start, end, flank_type
    except (IndexError, ValueError):
        raise ValueError(f"Cannot parse subject_id: {subject_id}")

# Define input and output file names
input_filename = './BLAST/flanking_sequences.out'
output_filename = './BLAST/flanking_sequences_processed.out'
# Open the input file and read the content

with open(input_filename, 'r') as file:
    lines = file.readlines()

# Filter out lines containing '#'
filtered_lines = [line for line in lines if '#' not in line]

# Write the filtered content to the output file
with open(output_filename, 'w') as file:
    file.writelines(filtered_lines)


# Screening process based on BLAST results
df = pd.read_csv('./BLAST/flanking_sequences_processed.out', sep='\t', header=None)

df.columns = ['query id', 'subject id', 'identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score']

df_filter = df[df['identity'] > 98]

fasta_file = './5k/flanking_sequences.fa'
subject_seq_length = get_length_file(fasta_file)
df_filter = df_filter.merge(subject_seq_length, left_on='subject id', right_on='reference_id', how='left')
df_filter = df_filter.drop(columns=['reference_id'])
# Remove sequences shorter than 200bp
df_filter = df_filter[df_filter['reference_length'] > 200]
df_filter = df_filter[df_filter['alignment length'] > df_filter['reference_length']*0.4]
# df_filter = df_filter[df_filter['mismatches'] <= df_filter['alignment length']*0.95]
df_filter = df_filter[((df_filter['subject id'].str.contains('right_flank', na=False)) & ((df_filter['s start'] == 1) | (df_filter['s end'] == 1))) | ((df_filter['subject id'].str.contains('left_flank', na=False)) & ((df_filter['s start'] == df_filter['reference_length']) | (df_filter['s end'] == df_filter['reference_length'])))]


# Remove duplicate sequences
# Group by 'subject id' and 'query id' together, and filter out groups with quantity greater than 1
duplicated_groups = df_filter.groupby(['subject id', 'query id']).filter(lambda x: len(x) > 1)

# Extract 'subject id' of groups with quantity greater than 1
duplicated_subject_ids = duplicated_groups['subject id'].unique()

# Use these 'subject id' to filter out rows in df_filter
df_filter = df_filter[~df_filter['subject id'].isin(duplicated_subject_ids)]
df_filter['q start'] = df_filter['q start'].astype(int)
df_filter['q end'] = df_filter['q end'].astype(int)


# fasta_file = 'first_1000_reads.fa'

# subject_seq_length = get_length_file(fasta_file)


# with open("flanking_sequences.fa", "r") as input_file:
#     all_flanking_sequences = list(SeqIO.parse(input_file, "fasta"))

# Gap filling process

# Read all sequences in the fasta file
all_sequences = load_sequences('./Assembly/assembly.fasta')
all_flanking_sequences = load_sequences('./5k/flanking_sequences.fa')
flanking_sequences_id = get_sequences_id('./5k/flanking_sequences.fa')

# Perform bilateral matching first

all_flanking_sequences_origin = all_flanking_sequences


bilateral_matched_subject_ids = []
error_gaps = []
bilateral_matched_subject_len = 0
# Initialize result storage list
results = []
for i in range(len(df_filter)):
    print(i)
    # query_seq = extract_sequences('filter.asm.bp.p_ctg.fa', df_filter.iloc[i]['query id'])
    # query_seq = query_seq[0].seq
    # subject_seq = extract_sequences('flanking_sequences.fa', df_filter.iloc[i]['subject id'])
    # subject_seq = subject_seq[0].seq
    query_seq_id = df_filter.iloc[i]['query id']
    subject_seq_id = df_filter.iloc[i]['subject id']
    
    # Extract sequences from the dictionary
    query_seq = extract_sequence(all_sequences, query_seq_id)
    subject_seq = extract_sequence(all_flanking_sequences, subject_seq_id)
    # Determine if the read fills the gap
    match_flanking_id = get_match_flanking(df_filter.iloc[i]['subject id'], flanking_sequences_id)
    
    df_detection = df_filter[(df_filter['query id'] == df_filter.iloc[i]['query id']) & (df_filter['subject id'] == match_flanking_id)]

    if len(df_detection) > 0:
        
        # Determine if bilateral matching has been performed
        if df_filter.iloc[i]['subject id'] in bilateral_matched_subject_ids:
            continue
        subject_seq2 = extract_sequence(all_flanking_sequences, match_flanking_id)
        subject_seq2 = subject_seq2[0]
        
        bilateral_matched_subject_ids.append(df_filter.iloc[i]['subject id'])
        bilateral_matched_subject_ids.append(match_flanking_id)
        
        # Extract gap position
        if "left"  in subject_seq_id:
            chrom, _, gap_start, _ = parse_subject_id(subject_seq_id)
        else:
            chrom, gap_end, _, _ = parse_subject_id(subject_seq_id)
        if "left"  in match_flanking_id:
            chrom, _, gap_start, _ = parse_subject_id(match_flanking_id)
        else:
            chrom, gap_end, _, _ = parse_subject_id(match_flanking_id)
        gap_start = gap_start +1
        gap_end = gap_end - 1
        
        
        if (df_filter.iloc[i]['s start'] < df_filter.iloc[i]['s end']):
            
            if "left"  in df_filter.iloc[i]['subject id']:
                seq_start = df_filter.iloc[i]['q end'] + 1
                seq_end = df_detection.iloc[0]['q start']
            else:
                seq_start = df_detection.iloc[0]['q end'] + 1
                seq_end = df_filter.iloc[i]['q start']
            
            new_subject_seq = subject_seq + query_seq[seq_start:seq_end] + subject_seq2
            
            bilateral_matched_subject_len = bilateral_matched_subject_len + len(query_seq[seq_start:seq_end])
            
            filled_seq = query_seq[seq_start:seq_end]
            
            # Determine if an unexpected gap occurs
            if len(filled_seq) < 1:
                error_gaps.append(df_filter.iloc[i]['subject id'])
                error_gaps.append(match_flanking_id)

                filled_seq = query_seq[seq_start:seq_end]
                # Record results
                results.append({
                    'chromosome': chrom,
                    'flank_type': 'bilateral_error',
                    'start_pos': gap_start,
                    'end_pos': gap_end,
                    'filled_sequence': str(filled_seq)
                })
                continue
            
            
        if (df_filter.iloc[i]['s start'] > df_filter.iloc[i]['s end']):
            
            if "left"  in df_filter.iloc[i]['subject id']:
                seq_start = df_detection.iloc[0]['q end']
                seq_end = df_filter.iloc[i]['q start'] + 1
            else:
                seq_start = df_filter.iloc[i]['q end']
                seq_end = df_detection.iloc[0]['q start'] + 1
            
            # Get complementary strand
            complement_seq = (Seq(query_seq._data[seq_start:seq_end])).complement()
            
            # Determine if an unexpected gap occurs
            if len(complement_seq) < 1:
                error_gaps.append(df_filter.iloc[i]['subject id'])
                error_gaps.append(match_flanking_id)
                
                
                complement_seq = (Seq(query_seq._data[seq_end:seq_start])).complement()
                # Record results
                results.append({
                    'chromosome': chrom,
                    'flank_type': 'bilateral_error',
                    'start_pos': gap_start,
                    'end_pos': gap_end,
                    'filled_sequence': str(complement_seq)
                })
                continue

            # Reverse the complementary strand to get the reverse complementary strand
            reversed_complement_seq = complement_seq[::-1]
            
            new_subject_seq = subject_seq + reversed_complement_seq + subject_seq2
            
            bilateral_matched_subject_len = bilateral_matched_subject_len + len(reversed_complement_seq)
            
            filled_seq = reversed_complement_seq
            
        flag=0
        for seq_name in all_flanking_sequences:
            if seq_name == df_filter.iloc[i]['subject id'] or seq_name == match_flanking_id:
                all_flanking_sequences[seq_name] = new_subject_seq
                flag += 1
            if flag == 2:
                break

        # Record results
        results.append({
            'chromosome': chrom,
            'flank_type': 'bilateral',
            'start_pos': gap_start,
            'end_pos': gap_end,
            'filled_sequence': str(filled_seq)
        })
        
# Unilateral matching

df_filter2 = df_filter[~df_filter['subject id'].isin(bilateral_matched_subject_ids)]

df_filter2['q start'] = df_filter2['q start'].astype(int)
df_filter2['q end'] = df_filter2['q end'].astype(int)

unilateral_matched_subject_len = 0
unilateral_matched_subject_ids = []

for i in range(len(df_filter2)):
    print(i)
    query_seq_id = df_filter2.iloc[i]['query id']
    subject_seq_id = df_filter2.iloc[i]['subject id']
    
    # Extract sequences from the dictionary
    query_seq = extract_sequence(all_sequences, query_seq_id)
    subject_seq = extract_sequence(all_flanking_sequences, subject_seq_id)
    # Determine if the read fills the gap
    match_flanking_id = get_match_flanking(df_filter2.iloc[i]['subject id'], flanking_sequences_id)
    
    df_detection = df_filter2[(df_filter2['query id'] == df_filter2.iloc[i]['query id']) & (df_filter2['subject id'] == match_flanking_id)]
    
    if len(df_detection) == 0: 
        
        # Parse position information
        chrom, original_start, original_end, flank_type = parse_subject_id(subject_seq_id)
        
        # if df_filter.iloc[i]['subject id'] in unilateral_matched_subject_ids:
        #     continue
        unilateral_matched_subject_ids.append(df_filter.iloc[i]['subject id'])
    
        # Fill the read into the sequence
        if ('right_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] < df_filter2.iloc[i]['s end']):
            new_subject_seq = query_seq[:df_filter2.iloc[i]['q start']-1] + subject_seq
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(query_seq[:df_filter2.iloc[i]['q start']-1])
            
            filled_seq = query_seq[:df_filter2.iloc[i]['q start']-1]

            
        elif ('right_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] > df_filter2.iloc[i]['s end']):
            # Get complementary strand
            complement_seq = (query_seq[df_filter2.iloc[i]['q end']:]).complement()
            
            # Reverse the complementary strand to get the reverse complementary strand
            reversed_complement_seq = complement_seq[::-1]
            new_subject_seq = reversed_complement_seq + subject_seq
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(reversed_complement_seq)
            
            filled_seq = reversed_complement_seq
            
        elif ('left_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] < df_filter2.iloc[i]['s end']):
            new_subject_seq = subject_seq + query_seq[df_filter2.iloc[i]['q end']:]
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(query_seq[df_filter2.iloc[i]['q end']:])
            
            filled_seq = query_seq[df_filter2.iloc[i]['q end']:]
            
        elif ('left_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] > df_filter2.iloc[i]['s end']):
            # Get complementary strand
            complement_seq = (query_seq[:df_filter2.iloc[i]['q start']-1]).complement()

            # Reverse the complementary strand to get the reverse complementary strand
            reversed_complement_seq = complement_seq[::-1]
            new_subject_seq = subject_seq + reversed_complement_seq
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(reversed_complement_seq)
            
            filled_seq = reversed_complement_seq
            
        # Modify the filled sequence to flanking_sequences
        for seq_name in all_flanking_sequences:
            if (seq_name == df_filter2.iloc[i]['subject id']) & (len(new_subject_seq) > len(all_flanking_sequences[seq_name])):
                all_flanking_sequences[seq_name] = new_subject_seq
                break
        
        # Record results
        results.append({
            'chromosome': chrom,
            'flank_type': flank_type,
            'start_pos': original_start,
            'end_pos': original_end,
            'filled_sequence': str(filled_seq)
        })

matched_subject_len = unilateral_matched_subject_len + bilateral_matched_subject_len
bilateral_matched_subject_ids = pd.Series(bilateral_matched_subject_ids).unique()
unilateral_matched_subject_ids = pd.Series(unilateral_matched_subject_ids).unique()


# ===== Generate final table =====
df_result = pd.DataFrame(results)

# Adjust column order and handle null values
final_columns = ['chromosome', 'flank_type', 'start_pos', 'end_pos', 'filled_sequence']
df_final = df_result[final_columns].copy()

# Save results
df_final.to_csv('./5k/gap_filling_results.csv', index=False)