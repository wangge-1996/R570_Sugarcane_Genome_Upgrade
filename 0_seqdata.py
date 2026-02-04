from Bio import SeqIO
from utils import *

# Extract sequences flanking the gaps
output_file = "./5k/flanking_sequences.fa"
# Load the genome sequence from the FASTA file
for genome_record in SeqIO.parse("SofficinarumxspontaneumR570_771_v2.0.fa", "fasta"):
    genome_seq = genome_record.seq

    # Find all gaps in the genome sequence
    gaps = find_gaps(genome_seq)
    extract_sequences(output_file, genome_seq, gaps, genome_record)


# Read the head and tail sequences from the fasta file
fasta_file = 'SofficinarumxspontaneumR570_771_v2.0.fa'
extracted_seqs = extract_head_tail_sequences(fasta_file)

# Write the extracted sequences to a new fasta file
with open('./5k/head_tail_sequences.fa', 'w') as output_handle:
    for seq_id, (head_seq, tail_seq) in extracted_seqs.items():
        # Write the head sequence
        output_handle.write(f">{seq_id}_head\n{head_seq}\n")
        # Write the tail sequence
        output_handle.write(f">{seq_id}_tail\n{tail_seq}\n")