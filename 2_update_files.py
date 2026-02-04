from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from collections import defaultdict
import io

# Initialize data structures
reference = SeqIO.to_dict(SeqIO.parse("SofficinarumxspontaneumR570_771_v2.0.fa", "fasta"))
df = pd.read_csv("./5k/final_merged_results.csv")
bam_records = []

# Preprocessing: Filter invalid data and group by chromosome
df = df[df['filled_sequence'].notna()].copy()
df['original_start'] = df['start_pos']
df['original_end'] = df['end_pos']

# Sort by chromosome and original position
df = df.sort_values(['chromosome', 'start_pos'])

# Read the original GFF file and store information
with open('SofficinarumxspontaneumR570_771_v2.1.gene.gff3', 'r') as gff_file:
    lines = []
    for line in gff_file:
        if not line.startswith('#') and len(line.strip().split('\t')) >= 9:
            lines.append(line)
    gff_content = io.StringIO(''.join(lines))
    gff_df = pd.read_csv(
        gff_content,
        sep='\t',
        header=None,
        names=['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )

# Initialize position shift mapping dictionary
position_shift = defaultdict(list)

# Iterate over and process each chromosome
for chrom, group in df.groupby('chromosome'):
    if chrom not in reference:
        continue
    
    # Sorting logic
    group = group.copy()
    group['sort_priority'] = group['type'].map({
        'head': 0,
        'tail': 2
    }).fillna(1)
    group = group.sort_values(
        by=['sort_priority', 'original_start'],
        ascending=[True, True]
    ).drop(columns=['sort_priority'])
    
    # Initialize chromosome sequence and offset
    chrom_seq = list(str(reference[chrom].seq))
    offset = 0
    chrom_length = len(chrom_seq)
    
    # Process each operation
    for _, row in group.iterrows():
        seq = str(row['filled_sequence'])
        orig_start = row['original_start']
        orig_end = row['original_end']
        op_type = row['type']
        
        adj_start = orig_start + offset
        adj_end = orig_end + offset if not pd.isna(orig_end) else None
        
        try:
            if op_type == 'head':
                chrom_seq = list(seq) + chrom_seq
                offset += len(seq)
                shift = len(seq)
                # Record offset of original position 0
                position_shift[chrom].append((0, shift))
                bam_type = 'head'
                new_start, new_end = 0, len(seq) - 1
                
            elif op_type == 'tail':
                chrom_seq += list(seq)
                bam_type = 'tail'
                new_start, new_end = chrom_length, chrom_length + len(seq) - 1
                # Tail operation does not affect other positions, no need to record shift
                
            elif op_type == 'bilateral':
                if adj_end < adj_start or adj_end > len(chrom_seq):
                    continue
                del chrom_seq[adj_start:adj_end + 1]
                chrom_seq[adj_start:adj_start] = list(seq)
                shift = (len(seq) - (orig_end - orig_start + 1))
                # Record offset of original start position
                position_shift[chrom].append((orig_start, shift))
                offset += shift
                bam_type = 'bilateral'
                new_start, new_end = adj_start, adj_start + len(seq) - 1
                
            elif op_type == 'left_flank':
                chrom_seq[adj_end + 1:adj_end + 1] = list(seq)
                shift = len(seq)
                # Record offset of original end position + 1
                position_shift[chrom].append((orig_end + 1, shift))
                offset += shift
                bam_type = 'left_flank'
                new_start, new_end = adj_end + 1, adj_end + len(seq)
                
            elif op_type == 'right_flank':
                chrom_seq[adj_start:adj_start] = list(seq)
                shift = len(seq)
                # Record offset of original start position
                position_shift[chrom].append((orig_start, shift))
                offset += shift
                bam_type = 'right_flank'
                new_start, new_end = adj_start, adj_start + len(seq) - 1
                
            elif op_type == 'bilateral_error':
                if adj_end < adj_start or adj_end > len(chrom_seq):
                    continue
                # Delete the original gap region
                del chrom_seq[adj_start:adj_end + 1]
                gap_length = adj_end - adj_start + 1
                shift_gap = -(orig_end - orig_start + 1)
                position_shift[chrom].append((orig_start, shift_gap))
                offset -= gap_length
                # Process left deletion
                if len(seq) > 0:
                    delete_length = len(seq)
                    left_delete_start = max(0, adj_start - delete_length)
                    left_delete_end = adj_start - 1
                    if left_delete_end >= left_delete_start:
                        del_length = left_delete_end - left_delete_start + 1
                        del chrom_seq[left_delete_start:left_delete_end + 1]
                        shift_left = -del_length
                        # Calculate original left deletion start position
                        orig_left_start = left_delete_start - offset + gap_length  # Position before offset adjustment
                        position_shift[chrom].append((orig_left_start, shift_left))
                        offset -= del_length
                bam_type = 'bilateral_error'
                new_start, new_end = adj_start, adj_start
                
            else:
                continue

            # Record BAM information
            bam_records.append({
                'chromosome': chrom,
                'start_pos': new_start + 1,
                'end_pos': new_end + 1,
                'type': bam_type
            })
            chrom_length = len(chrom_seq)
            
        except IndexError:
            continue
    
    reference[chrom].seq = Seq(''.join(chrom_seq))

# Save the updated reference genome
SeqIO.write(reference.values(), "./5k/SofficinarumxspontaneumR570_771_v2.1.update.fa", "fasta")

# Generate BAM file
bam_df = pd.DataFrame(bam_records)
bam_df.to_csv("./5k/insertion_positions.bed", sep='\t', index=False, header=False)

# Update GFF positions
for idx in gff_df.index:
    chrom = gff_df.loc[idx, 'chromosome']
    start = gff_df.loc[idx, 'start']
    end = gff_df.loc[idx, 'end']
    total_shift = 0
    # Iterate over all original position offset records for the chromosome
    for orig_pos, shift in position_shift.get(chrom, []):
        if start >= orig_pos:
            total_shift += shift
    gff_df.loc[idx, 'start'] += total_shift
    gff_df.loc[idx, 'end'] += total_shift

# Save the updated GFF
gff_df.to_csv('./5k/SofficinarumxspontaneumR570_771_v2.1.gene_update.gff3', sep='\t', index=False, header=False)

print("Processing completed!")