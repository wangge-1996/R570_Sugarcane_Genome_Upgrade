from Bio import SeqIO
import re
import json
import pandas as pd
import os


def locate_gaps(fasta_file):
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for match in re.finditer('N+', str(record.seq)):
                start = match.start()
                end = match.end() - 1
                print(f"{record.id}\t.\tgap\t{start}\t{end}\t.\t.\tName=gap{match.start()+1};size={end-start+1}")

def find_gaps(genome_seq, min_gap_size=1):
    gaps = []
    current_gap_start = None
    for i in range(len(genome_seq)):
        if genome_seq[i] == 'N':
            if current_gap_start is None:
                current_gap_start = i
        else:
            if current_gap_start is not None:
                current_gap_end = i
                if (current_gap_end - current_gap_start + 1) >= min_gap_size:
                    gaps.append((current_gap_start, current_gap_end))
                current_gap_start = None
    # Check for a gap that ends at the end of the sequence
    if current_gap_start is not None:
        if (len(genome_seq) - current_gap_start + 1) >= min_gap_size:
            gaps.append((current_gap_start, len(genome_seq) - 1))
    return gaps

def extract_sequences(output_file, genome_seq, gaps, genome_record, flank_size=5000):
    with open(output_file, "a+") as output_file:
        for gap_start, gap_end in gaps:
            # Calculate the flanking regions
            left_start = max(0, gap_start - flank_size)
            left_end = gap_start-1
            right_start = gap_end
            right_end = min(len(genome_seq) - 1, gap_end + flank_size - 1)

            # Extract the left flanking sequence (Tips: Python's slices are left closed and right open.)
            left_flank_seq = ''
            for i in range(left_end, left_start - 1, -1):  # Traverse reversely from left_end to left_start
                if genome_seq[i] == 'N':
                    left_flank_seq = genome_seq[i + 1:left_end + 1]  # Include the base at position i
                    left_start = i + 1  # Update left_start to the next position of 'N'
                    break
            else:  # No 'N' found in the left flanking region
                left_flank_seq = genome_seq[left_start:left_end + 1]

            # Extract the right flanking sequence and detect 'N'
            right_flank_seq = ''
            for i in range(right_start, right_end + 1):
                if genome_seq[i] == 'N':
                    right_flank_seq = genome_seq[right_start:i]
                    right_end = i - 1
                    break
            else:  # No 'N' found in the right flanking region
                right_flank_seq = genome_seq[right_start:right_end + 1]

            output_file.write(f">{genome_record.id}_{left_start}_{left_end}_left_flank\n{left_flank_seq}\n")
            output_file.write(f">{genome_record.id}_{right_start}_{right_end}_right_flank\n{right_flank_seq}\n")

def extract_head_tail_sequences(fasta_file, flank_size=5000):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    extracted_sequences = {}

    for seq_id, seq_record in sequences.items():
        # Extract the head sequence
        head_seq = seq_record.seq[:flank_size]
        # Check if the head sequence contains 'N'
        if 'N' in head_seq:
            head_start = head_seq.find('N')  # Find the position of the first 'N'
            head_seq = head_seq[:head_start]  # Extract the sequence before 'N'

        # Extract the tail sequence
        tail_seq = seq_record.seq[-flank_size:]  # Extract the last 5000bp of the tail
        # Check if the tail sequence contains 'N'
        if 'N' in tail_seq:
            tail_start = tail_seq.rfind('N')  # Find the position of the last 'N'
            tail_seq = tail_seq[tail_start+1:]  # Extract the sequence after 'N'

        # Store the extracted sequences in a dictionary
        extracted_sequences[seq_id] = (head_seq, tail_seq)

    return extracted_sequences


def extract_first_1000_sequences(input_file, output_file, read_num):
    count = 0
    with open(input_file, "r") as input_handle:
        with open(output_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                if count < read_num:
                    SeqIO.write(record, output_handle, "fasta")
                    count += 1
                else:
                    break

def extract_sequence_info_to_file(fasta_file, output_file):
    total_length = 0
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        for record in SeqIO.parse(fasta, 'fasta'):
            sequence_name = record.id
            sequence_length = len(record.seq)
            output.write(f"Sequence Name: {sequence_name}, Length: {sequence_length} bp\n")
            total_length += sequence_length
        output.write(f"Total length of all sequences: {total_length} bp\n")


def get_match_flanking(seq_id, flanking_sequences_id):
# Find the position of A in the list
    index = flanking_sequences_id.index(seq_id)
    
    # Check the keywords contained in A and select the next or previous string according to the rules
    if "left" in seq_id:
        # If A contains "left", select the next string
        if index + 1 < len(flanking_sequences_id):
            return flanking_sequences_id[index + 1]
        else:
            print("No next string in the list.")
    elif "right" in seq_id:
        # If A contains "right", select the previous string
        if index > 0:
            return flanking_sequences_id[index - 1]
        else:
            print("No previous string in the list.")

# Extract all sequence id
def get_sequences_id(fasta_file):
    with open(fasta_file, "r") as fasta_file:
        # Initialize an empty list to store sequence names
        sequence_names = []
        # Read each line of the file
        for line in fasta_file:
            # Check if it is a sequence name line
            if line.startswith(">"):
                # Remove the ">" at the beginning of the line and add to the list
                sequence_names.append(line[1:].strip())
    return sequence_names

# Read all sequences in the fasta file at once and store them in a dictionary
def load_sequences(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            sequences[record.id] = record.seq
    return sequences

# Extract the sequence of a specific ID
def extract_sequence(sequences, seq_id):
    return sequences.get(seq_id, None)

# Save the sequence as a json file format, dedicated to alphafold3
def save_json(sequence, path):
    data = {
        "name": str(sequence.id),
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence":str(sequence.seq)
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }
    
    with open(path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# Get sequence length
def get_length_file(fasta_file):
    extracted_sequences=[]
    with open(fasta_file, 'r', encoding='utf-8') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            extracted_sequences.append([record.id,len(record.seq)])
    
    df = pd.DataFrame(extracted_sequences)
    df.columns = ['reference_id', 'reference_length']
    
    return df

# Process blast files
def process_blast_results(input_file, query_lengths, reference_lengths, output_file):

    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Filter out lines that contain '#'
    filtered_lines = [line for line in lines if '#' not in line]

    with open('./Results/blast_processed.out', 'w') as file:
        file.writelines(filtered_lines)

    # Read the BLAST result file
    df = pd.read_csv('./Results/blast_processed.out', sep='\t', header=None)
    # Take the first four columns 
    df= df.iloc[:,:4 ]
    df.columns = ['query_id', 'reference_id', 'similarity', 'alignment_length']
    
    # The sum of the weighted similarity level and comparison length is calculated, grouped by query_id and reference_id
    grouped = df.groupby(['query_id', 'reference_id']).apply(lambda x: pd.Series({
        'similarity': (x['similarity'] * x['alignment_length']).sum() / x['alignment_length'].sum(),
        'alignment_length': float(x['alignment_length'].sum())
    })).reset_index()
    
    # Read length information for query and reference
    query_lengths.columns = ['query_id', 'query_length']
    
    reference_lengths.columns = ['reference_id', 'reference_length']
    
    # Merge length information to the grouped data frame
    grouped = pd.merge(grouped, query_lengths, on='query_id', how='left')
    grouped = pd.merge(grouped, reference_lengths, on='reference_id', how='left')

    grouped['query_length'] = grouped['query_length'].astype(float)
    grouped['reference_length'] = grouped['reference_length'].astype(float)
    
    grouped['alignment_length/query_length'] = grouped['alignment_length'] / grouped['query_length']
    grouped['alignment_length/reference_length'] = grouped['alignment_length'] / grouped['reference_length']

    # Save the result to a new file
    grouped.to_csv(output_file, sep='\t', header=False, index=False)