from Bio import SeqIO
import pandas as pd
import os
from utils import *
from Bio.Seq import Seq
import re

def parse_subject_id(subject_id):
    # 使用正则表达式匹配不同格式的subject_id
    pattern = r'''
        ^(.*?_)           # 染色体部分（包含下划线）
        (\d+)_(\d+)_      # 两个位置数字
        (left|right)_flank$  # 完整的flank类型
    '''
    match = re.match(pattern, subject_id, re.VERBOSE)
    
    if match:
        chrom = match.group(1).rstrip('_')  # 去除末尾可能的下划线
        start = int(match.group(2))
        end = int(match.group(3))
        flank_type = f"{match.group(4)}_flank"
        return chrom, start, end, flank_type
    
    # 备用方案：处理其他可能格式
    parts = subject_id.rsplit('_', 4)  # 从右侧分割4次
    try:
        chrom_part = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        flank_type = '_'.join(parts[3:])
        return chrom_part, start, end, flank_type
    except (IndexError, ValueError):
        raise ValueError(f"无法解析subject_id: {subject_id}")

# 定义输入和输出文件的名称
input_filename = './5k/flanking_sequences.out'
output_filename = './5k/flanking_sequences_processed.out'
# 打开输入文件，读取内容

with open(input_filename, 'r') as file:
    lines = file.readlines()

# 过滤掉包含 '#' 的行
filtered_lines = [line for line in lines if '#' not in line]

# 将过滤后的内容写入到输出文件
with open(output_filename, 'w') as file:
    file.writelines(filtered_lines)


# 依据 Blast 结果筛选过程
df = pd.read_csv('./5k/flanking_sequences_processed.out', sep='\t', header=None)

df.columns = ['query id', 'subject id', 'identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score']

df_filter = df[df['identity'] > 98]


fasta_file = './5k/flanking_sequences.fa'
subject_seq_length = get_length_file(fasta_file)
df_filter = df_filter.merge(subject_seq_length, left_on='subject id', right_on='reference_id', how='left')
df_filter = df_filter.drop(columns=['reference_id'])
# 去掉200bp以下的序列
df_filter = df_filter[df_filter['reference_length'] > 200]
df_filter = df_filter[df_filter['alignment length'] > df_filter['reference_length']*0.4]
# df_filter = df_filter[df_filter['mismatches'] <= df_filter['alignment length']*0.95]
df_filter = df_filter[((df_filter['subject id'].str.contains('right_flank', na=False)) & ((df_filter['s start'] == 1) | (df_filter['s end'] == 1))) | ((df_filter['subject id'].str.contains('left_flank', na=False)) & ((df_filter['s start'] == df_filter['reference_length']) | (df_filter['s end'] == df_filter['reference_length'])))]


# 去掉重复序列
# 根据 'subject id' 和 'query id' 共同分组，并筛选出数量大于1的组
duplicated_groups = df_filter.groupby(['subject id', 'query id']).filter(lambda x: len(x) > 1)

# 提取数量大于1的组的 'subject id'
duplicated_subject_ids = duplicated_groups['subject id'].unique()

# 使用这些 'subject id' 来筛选掉 df_filter 中的行
df_filter = df_filter[~df_filter['subject id'].isin(duplicated_subject_ids)]
df_filter['q start'] = df_filter['q start'].astype(int)
df_filter['q end'] = df_filter['q end'].astype(int)


# fasta_file = 'first_1000_reads.fa'

# subject_seq_length = get_length_file(fasta_file)


# with open("flanking_sequences.fa", "r") as input_file:
#     all_flanking_sequences = list(SeqIO.parse(input_file, "fasta"))

# 填补过程

# 读取fasta文件中的所有序列
all_sequences = load_sequences('./ZC/assembly.fasta')
all_flanking_sequences = load_sequences('./5k/flanking_sequences.fa')
flanking_sequences_id = get_sequences_id('./5k/flanking_sequences.fa')

# 先进行双边匹配

all_flanking_sequences_origin = all_flanking_sequences


bilateral_matched_subject_ids = []
error_gaps = []
bilateral_matched_subject_len = 0
# 初始化结果存储列表
results = []
for i in range(len(df_filter)):
    print(i)
    # query_seq = extract_sequences('filter.asm.bp.p_ctg.fa', df_filter.iloc[i]['query id'])
    # query_seq = query_seq[0].seq
    # subject_seq = extract_sequences('flanking_sequences.fa', df_filter.iloc[i]['subject id'])
    # subject_seq = subject_seq[0].seq
    query_seq_id = df_filter.iloc[i]['query id']
    subject_seq_id = df_filter.iloc[i]['subject id']
    
    # 从字典中提取序列
    query_seq = extract_sequence(all_sequences, query_seq_id)
    subject_seq = extract_sequence(all_flanking_sequences, subject_seq_id)
    # 判断该reads是否填补了gap
    match_flanking_id = get_match_flanking(df_filter.iloc[i]['subject id'], flanking_sequences_id)
    
    df_detection = df_filter[(df_filter['query id'] == df_filter.iloc[i]['query id']) & (df_filter['subject id'] == match_flanking_id)]

    if len(df_detection) > 0:
        
        # 判断是否已经进行了双边匹配
        if df_filter.iloc[i]['subject id'] in bilateral_matched_subject_ids:
            continue
        subject_seq2 = extract_sequence(all_flanking_sequences, match_flanking_id)
        subject_seq2 = subject_seq2[0]
        
        bilateral_matched_subject_ids.append(df_filter.iloc[i]['subject id'])
        bilateral_matched_subject_ids.append(match_flanking_id)
        
        #提取gap位置
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
            
            #判断是否出现了不应该出现的gap
            if len(filled_seq) < 1:
                error_gaps.append(df_filter.iloc[i]['subject id'])
                error_gaps.append(match_flanking_id)

                filled_seq = query_seq[seq_start:seq_end]
                # 记录结果
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
            
            # 获取互补链
            complement_seq = (Seq(query_seq._data[seq_start:seq_end])).complement()
            
            #判断是否出现了不应该出现的gap
            if len(complement_seq) < 1:
                error_gaps.append(df_filter.iloc[i]['subject id'])
                error_gaps.append(match_flanking_id)
                
                
                complement_seq = (Seq(query_seq._data[seq_end:seq_start])).complement()
                # 记录结果
                results.append({
                    'chromosome': chrom,
                    'flank_type': 'bilateral_error',
                    'start_pos': gap_start,
                    'end_pos': gap_end,
                    'filled_sequence': str(complement_seq)
                })
                continue

            # 反转互补链以获得反向互补链
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

        # 记录结果
        results.append({
            'chromosome': chrom,
            'flank_type': 'bilateral',
            'start_pos': gap_start,
            'end_pos': gap_end,
            'filled_sequence': str(filled_seq)
        })
        
# 单边匹配

df_filter2 = df_filter[~df_filter['subject id'].isin(bilateral_matched_subject_ids)]

df_filter2['q start'] = df_filter2['q start'].astype(int)
df_filter2['q end'] = df_filter2['q end'].astype(int)

unilateral_matched_subject_len = 0
unilateral_matched_subject_ids = []

for i in range(len(df_filter2)):
    print(i)
    query_seq_id = df_filter2.iloc[i]['query id']
    subject_seq_id = df_filter2.iloc[i]['subject id']
    
    # 从字典中提取序列
    query_seq = extract_sequence(all_sequences, query_seq_id)
    subject_seq = extract_sequence(all_flanking_sequences, subject_seq_id)
    # 判断该reads是否填补了gap
    match_flanking_id = get_match_flanking(df_filter2.iloc[i]['subject id'], flanking_sequences_id)
    
    df_detection = df_filter2[(df_filter2['query id'] == df_filter2.iloc[i]['query id']) & (df_filter2['subject id'] == match_flanking_id)]
    
    if len(df_detection) == 0: 
        
        # 解析位置信息
        chrom, original_start, original_end, flank_type = parse_subject_id(subject_seq_id)
        
        # if df_filter.iloc[i]['subject id'] in unilateral_matched_subject_ids:
        #     continue
        unilateral_matched_subject_ids.append(df_filter.iloc[i]['subject id'])
    
        # reads填补到seq中
        if ('right_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] < df_filter2.iloc[i]['s end']):
            new_subject_seq = query_seq[:df_filter2.iloc[i]['q start']-1] + subject_seq
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(query_seq[:df_filter2.iloc[i]['q start']-1])
            
            filled_seq = query_seq[:df_filter2.iloc[i]['q start']-1]

            
        elif ('right_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] > df_filter2.iloc[i]['s end']):
            # 获取互补链
            complement_seq = (query_seq[df_filter2.iloc[i]['q end']:]).complement()
            
            # 反转互补链以获得反向互补链
            reversed_complement_seq = complement_seq[::-1]
            new_subject_seq = reversed_complement_seq + subject_seq
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(reversed_complement_seq)
            
            filled_seq = reversed_complement_seq
            
        elif ('left_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] < df_filter2.iloc[i]['s end']):
            new_subject_seq = subject_seq + query_seq[df_filter2.iloc[i]['q end']:]
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(query_seq[df_filter2.iloc[i]['q end']:])
            
            filled_seq = query_seq[df_filter2.iloc[i]['q end']:]
            
        elif ('left_flank' in df_filter2.iloc[i]['subject id']) & (df_filter2.iloc[i]['s start'] > df_filter2.iloc[i]['s end']):
            # 获取互补链
            complement_seq = (query_seq[:df_filter2.iloc[i]['q start']-1]).complement()

            # 反转互补链以获得反向互补链
            reversed_complement_seq = complement_seq[::-1]
            new_subject_seq = subject_seq + reversed_complement_seq
            unilateral_matched_subject_len = unilateral_matched_subject_len + len(reversed_complement_seq)
            
            filled_seq = reversed_complement_seq
            
        # 修改填补后的seq到flanking_sequences
        for seq_name in all_flanking_sequences:
            if (seq_name == df_filter2.iloc[i]['subject id']) & (len(new_subject_seq) > len(all_flanking_sequences[seq_name])):
                all_flanking_sequences[seq_name] = new_subject_seq
                break
        
        # 记录结果
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


# ===== 生成最终表格 =====
df_result = pd.DataFrame(results)

# 调整列顺序和空值处理
final_columns = ['chromosome', 'flank_type', 'start_pos', 'end_pos', 'filled_sequence']
df_final = df_result[final_columns].copy()

# 保存结果
df_final.to_csv('./5k/gap_filling_results.csv', index=False)

