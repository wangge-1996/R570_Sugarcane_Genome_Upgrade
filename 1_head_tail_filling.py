from Bio import SeqIO
import pandas as pd
import os
from utils import *
from Bio.Seq import Seq

# New function: Parse new head/tail format
def parse_headtail_id(subject_id):
    """Parse formats like Chr10os3_tail or scaffold_117_head"""
    try:
        # Split once from the right to get the region type
        parts = subject_id.rsplit('_', 1)
        if len(parts) != 2 or parts[1] not in ['head', 'tail']:
            raise ValueError
        return parts[0], parts[1]
    except:
        raise ValueError(f"Cannot parse subject_id: {subject_id}")

# Define input and output file names
input_filename = './BLAST/head_tail_sequences.out'
output_filename = './BLAST/head_tail_sequences_processed.out'

# Open the input file and read the content
with open(input_filename, 'r') as file:
    lines = file.readlines()

# Filter out lines containing '#'
filtered_lines = [line for line in lines if '#' not in line]

# Write the filtered content to the output file
with open(output_filename, 'w') as file:
    file.writelines(filtered_lines)

# Screening process based on BLAST results
df = pd.read_csv('./BLAST/head_tail_sequences_processed.out', sep='\t', header=None)

df.columns = ['query id', 'subject id', 'identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score']

df_filter = df[df['identity'] > 98]

fasta_file = './5k/head_tail_sequences.fa'
subject_seq_length = get_length_file(fasta_file)
df_filter = df_filter.merge(subject_seq_length, left_on='subject id', right_on='reference_id', how='left')
df_filter = df_filter.drop(columns=['reference_id'])
# Remove sequences shorter than 200bp
df_filter = df_filter[df_filter['reference_length'] > 200]
df_filter = df_filter[df_filter['alignment length'] > df_filter['reference_length']*0.4]
# df_filter = df_filter[df_filter['mismatches'] <= df_filter['alignment length']*0.95]
df_filter = df_filter[((df_filter['subject id'].str.contains('head', na=False)) & ((df_filter['s start'] == 1) | (df_filter['s end'] == 1))) | ((df_filter['subject id'].str.contains('tail', na=False)) & ((df_filter['s start'] == df_filter['reference_length']) | (df_filter['s end'] == df_filter['reference_length'])))]

# Remove duplicate sequences
# Group by 'subject id' and 'query id' together, and filter out groups with quantity greater than 1
duplicated_groups = df_filter.groupby(['subject id', 'query id']).filter(lambda x: len(x) > 1)

# Extract 'subject id' of groups with quantity greater than 1
duplicated_subject_ids = duplicated_groups['subject id'].unique()

# Use these 'subject id' to filter out rows in df_filter
df_filter = df_filter[~df_filter['subject id'].isin(duplicated_subject_ids)]
df_filter['q start'] = df_filter['q start'].astype(int)
df_filter['q end'] = df_filter['q end'].astype(int)

# Extension process

# Read all sequences in the fasta file
all_sequences = load_sequences('./Assembly/assembly.fasta')
head_tail_sequences = load_sequences('./5k/head_tail_sequences.fa')

head_tail_matched_subject_ids = []
head_tail_sequences_len = 0

# Initialize result storage list
results = []

# ===== Main processing loop =====
for i in range(len(df_filter)):
    print(i)
    query_seq_id = df_filter.iloc[i]['query id']
    subject_seq_id = df_filter.iloc[i]['subject id']
    
    # Extract sequences from the dictionary
    query_seq = extract_sequence(all_sequences, query_seq_id)
    subject_seq = extract_sequence(head_tail_sequences, subject_seq_id)
    
    # if df_filter.iloc[i]['subject id'] in unilateral_matched_subject_ids:
    #     continue
    head_tail_matched_subject_ids.append(df_filter.iloc[i]['subject id'])
    
    # try:
    #     chrom, region_type = parse_headtail_id(subject_seq_id)
    # except ValueError:
    #     continue
    chrom, region_type = parse_headtail_id(subject_seq_id)
    
    # Fill the read into the sequence
    if ('head' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] < df_filter.iloc[i]['s end']):
        new_subject_seq = query_seq[:df_filter.iloc[i]['q start']-1] + subject_seq
        head_tail_sequences_len = head_tail_sequences_len + len(query_seq[:df_filter.iloc[i]['q start']-1])
        
        filled_seq = query_seq[:df_filter.iloc[i]['q start']-1]
        
    elif ('head' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] > df_filter.iloc[i]['s end']):
        # Get complementary strand
        complement_seq = (query_seq[df_filter.iloc[i]['q end']:]).complement()

        # Reverse the complementary strand to get the reverse complementary strand
        reversed_complement_seq = complement_seq[::-1]
        new_subject_seq = reversed_complement_seq + subject_seq
        head_tail_sequences_len = head_tail_sequences_len + len(reversed_complement_seq)
        
        filled_seq = reversed_complement_seq
        
    elif ('tail' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] < df_filter.iloc[i]['s end']):
        new_subject_seq = subject_seq + query_seq[df_filter.iloc[i]['q end']:]
        head_tail_sequences_len = head_tail_sequences_len + len(query_seq[df_filter.iloc[i]['q end']:])
        
        filled_seq = query_seq[df_filter.iloc[i]['q end']:]
        
    elif ('tail' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] > df_filter.iloc[i]['s end']):
        # Get complementary strand
        complement_seq = (query_seq[:df_filter.iloc[i]['q start']-1]).complement()

        # Reverse the complementary strand to get the reverse complementary strand
        reversed_complement_seq = complement_seq[::-1]
        
        new_subject_seq = subject_seq + reversed_complement_seq
        head_tail_sequences_len = head_tail_sequences_len + len(reversed_complement_seq)
        
        filled_seq = reversed_complement_seq
        
    # Modify the filled sequence to head_tail_sequences
    for seq_name in head_tail_sequences:
        if (seq_name == df_filter.iloc[i]['subject id']) & (len(new_subject_seq) > len(head_tail_sequences[seq_name])):
            head_tail_sequences[seq_name] = new_subject_seq
            break

    # Record results
    results.append({
        'chromosome': chrom,
        'region_type': region_type,
        'start_pos': 1,
        'end_pos': 1,
        'filled_sequence': str(filled_seq)
    })

head_tail_matched_subject_ids = pd.Series(head_tail_matched_subject_ids).unique()
# ===== Generate final table =====
df_result = pd.DataFrame(results)
df_result = df_result[['chromosome', 'region_type', 'start_pos', 'end_pos', 'filled_sequence']]

# Save results
df_result.to_csv('./5k/head_tail_extensions.csv', index=False)