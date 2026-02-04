import pandas as pd
import re

def get_chromosome_order(chrom):
    """Handle special rules for chromosome sorting"""
    # Extract base number (supports formats like Chr1A)
    base_num = re.findall(r'\d+', chrom)
    base_num = int(base_num[0]) if base_num else 9999
    
    # Sort by letter suffix (supports comparison like Chr2A vs Chr2B)
    suffix = re.sub(r'\d+', '', chrom).lower()
    
    # Sort special chromosome types
    type_order = {
        'un': 10000,
        'scaffold': 20000,
        'mito': 30000
    }
    for k in type_order:
        if k in chrom.lower():
            return (type_order[k], chrom)
    
    return (base_num, suffix, chrom)

def create_sort_key(row):
    """Generate sort key: (chromosome order, type priority, start position)"""
    # Process chromosome sorting
    chrom = str(row['chromosome'])
    chrom_order = get_chromosome_order(chrom)
    
    # Set type priority
    type_priority = {
        'head': 0,
        'left_flank': 1,
        'bilateral': 2,
        'right_flank': 3,
        'tail': 4
    }
    entry_type = row['type']
    priority = type_priority.get(entry_type, 5)  # Unknown types are placed last
    
    # Get valid position
    position = row['start_pos'] if pd.notnull(row['start_pos']) else row['end_pos']
    
    return (chrom_order, priority, position)

# Read and merge data
df_gap = pd.read_csv('./5k/gap_filling_results.csv')
df_headtail = pd.read_csv('./5k/head_tail_extensions.csv')

# Unify column names
merged_df = pd.concat([
    df_gap.rename(columns={'flank_type': 'type'}),
    df_headtail.rename(columns={'region_type': 'type'})
], ignore_index=True)

# Define grouping columns
group_cols = ['chromosome', 'type', 'start_pos', 'end_pos']

# Step 1: Filter valid groups
def filter_groups(group):
    # Preprocessing: filter empty sequences
    clean_group = group.dropna(subset=['filled_sequence'])
    
    # Return empty DataFrame for empty groups directly
    if clean_group.empty:
        return pd.DataFrame()
    
    # Case 1: All sequences are exactly the same
    if clean_group['filled_sequence'].nunique() == 1:
        return clean_group.head(1)
    
    # Case 2: Unique sequence in the group
    elif len(clean_group) == 1:
        return clean_group
    
    # Case 3: Keep the longest sequence
    else:
        # Calculate sequence lengths
        seq_lengths = clean_group['filled_sequence'].str.len()
        
        # Get indices of the longest sequence (keep the first occurrence)
        max_length = seq_lengths.max()
        max_indices = seq_lengths[seq_lengths == max_length].index
        
        # Return the first longest sequence
        return clean_group.loc[[max_indices[0]]]

# Process deduplication for bilateral_error type
# Split data into bilateral_error and other types
bilateral_mask = merged_df['type'] == 'bilateral_error'
df_bilateral = merged_df[bilateral_mask].copy()
df_other = merged_df[~bilateral_mask].copy()

# Deduplicate by chromosome, start and end positions (keep the first occurrence)
df_bilateral = df_bilateral.drop_duplicates(
    subset=['chromosome', 'start_pos', 'end_pos'], 
    keep='first'
)

# Merge data back together
merged_df = pd.concat([df_bilateral, df_other], ignore_index=True)

# Apply group filtering
filtered_df = merged_df.groupby(group_cols, group_keys=False).apply(filter_groups)

# Step 2: Generate sort keys and sort the DataFrame
filtered_df['sort_key'] = filtered_df.apply(create_sort_key, axis=1)
sorted_df = filtered_df.sort_values(by='sort_key').drop(columns=['sort_key'])

# Save the final result
sorted_df.to_csv('./5k/final_merged_results.csv', index=False)