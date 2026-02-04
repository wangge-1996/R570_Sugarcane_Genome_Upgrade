from Bio import SeqIO
import pandas as pd
import os
from utils import *
from Bio.Seq import Seq

# 新增函数：解析新的head/tail格式
def parse_headtail_id(subject_id):
    """解析类似 Chr10os3_tail 或 scaffold_117_head 的格式"""
    try:
        # 使用右侧分割一次获取区域类型
        parts = subject_id.rsplit('_', 1)
        if len(parts) != 2 or parts[1] not in ['head', 'tail']:
            raise ValueError
        return parts[0], parts[1]
    except:
        raise ValueError(f"无法解析subject_id: {subject_id}")



# 定义输入和输出文件的名称
input_filename = './5k/head_tail_sequences.out'
output_filename = './5k/head_tail_sequences_processed.out'
# 打开输入文件，读取内容

with open(input_filename, 'r') as file:
    lines = file.readlines()

# 过滤掉包含 '#' 的行
filtered_lines = [line for line in lines if '#' not in line]

# 将过滤后的内容写入到输出文件
with open(output_filename, 'w') as file:
    file.writelines(filtered_lines)

# 依据 Blast 结果筛选过程
df = pd.read_csv('./5k/head_tail_sequences_processed.out', sep='\t', header=None)

df.columns = ['query id', 'subject id', 'identity', 'alignment length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'evalue', 'bit score']

df_filter = df[df['identity'] > 98]


fasta_file = './5k/head_tail_sequences.fa'
subject_seq_length = get_length_file(fasta_file)
df_filter = df_filter.merge(subject_seq_length, left_on='subject id', right_on='reference_id', how='left')
df_filter = df_filter.drop(columns=['reference_id'])
# 去掉200bp以下的序列
df_filter = df_filter[df_filter['reference_length'] > 200]
df_filter = df_filter[df_filter['alignment length'] > df_filter['reference_length']*0.4]
# df_filter = df_filter[df_filter['mismatches'] <= df_filter['alignment length']*0.95]
df_filter = df_filter[((df_filter['subject id'].str.contains('head', na=False)) & ((df_filter['s start'] == 1) | (df_filter['s end'] == 1))) | ((df_filter['subject id'].str.contains('tail', na=False)) & ((df_filter['s start'] == df_filter['reference_length']) | (df_filter['s end'] == df_filter['reference_length'])))]

# 去掉重复序列
# 根据 'subject id' 和 'query id' 共同分组，并筛选出数量大于1的组
duplicated_groups = df_filter.groupby(['subject id', 'query id']).filter(lambda x: len(x) > 1)

# 提取数量大于1的组的 'subject id'
duplicated_subject_ids = duplicated_groups['subject id'].unique()

# 使用这些 'subject id' 来筛选掉 df_filter 中的行
df_filter = df_filter[~df_filter['subject id'].isin(duplicated_subject_ids)]
df_filter['q start'] = df_filter['q start'].astype(int)
df_filter['q end'] = df_filter['q end'].astype(int)

# 延伸过程

# 读取fasta文件中的所有序列
all_sequences = load_sequences('./ZC/assembly.fasta')
head_tail_sequences = load_sequences('./5k/head_tail_sequences.fa')

head_tail_matched_subject_ids = []
head_tail_sequences_len = 0

# 初始化结果存储列表
results = []

# ===== 主处理循环 =====
for i in range(len(df_filter)):
    print(i)
    query_seq_id = df_filter.iloc[i]['query id']
    subject_seq_id = df_filter.iloc[i]['subject id']
    
    # 从字典中提取序列
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
    
    # reads填补到seq中
    if ('head' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] < df_filter.iloc[i]['s end']):
        new_subject_seq = query_seq[:df_filter.iloc[i]['q start']-1] + subject_seq
        head_tail_sequences_len = head_tail_sequences_len + len(query_seq[:df_filter.iloc[i]['q start']-1])
        
        filled_seq = query_seq[:df_filter.iloc[i]['q start']-1]
        
    elif ('head' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] > df_filter.iloc[i]['s end']):
        # 获取互补链
        complement_seq = (query_seq[df_filter.iloc[i]['q end']:]).complement()

        # 反转互补链以获得反向互补链
        reversed_complement_seq = complement_seq[::-1]
        new_subject_seq = reversed_complement_seq + subject_seq
        head_tail_sequences_len = head_tail_sequences_len + len(reversed_complement_seq)
        
        filled_seq = reversed_complement_seq
        
    elif ('tail' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] < df_filter.iloc[i]['s end']):
        new_subject_seq = subject_seq + query_seq[df_filter.iloc[i]['q end']:]
        head_tail_sequences_len = head_tail_sequences_len + len(query_seq[df_filter.iloc[i]['q end']:])
        
        filled_seq = query_seq[df_filter.iloc[i]['q end']:]
        
    elif ('tail' in df_filter.iloc[i]['subject id']) & (df_filter.iloc[i]['s start'] > df_filter.iloc[i]['s end']):
        # 获取互补链
        complement_seq = (query_seq[:df_filter.iloc[i]['q start']-1]).complement()

        # 反转互补链以获得反向互补链
        reversed_complement_seq = complement_seq[::-1]
        
        new_subject_seq = subject_seq + reversed_complement_seq
        head_tail_sequences_len = head_tail_sequences_len + len(reversed_complement_seq)
        
        filled_seq = reversed_complement_seq
        
    # 修改填补后的seq到flanking_sequences
    for seq_name in head_tail_sequences:
        if (seq_name == df_filter.iloc[i]['subject id']) & (len(new_subject_seq) > len(head_tail_sequences[seq_name])):
            head_tail_sequences[seq_name] = new_subject_seq
            break

    # 记录结果
    results.append({
        'chromosome': chrom,
        'region_type': region_type,
        'start_pos': 1,
        'end_pos': 1,
        'filled_sequence': str(filled_seq)
    })

head_tail_matched_subject_ids = pd.Series(head_tail_matched_subject_ids).unique()
# ===== 生成最终表格 =====
df_result = pd.DataFrame(results)
df_result = df_result[['chromosome', 'region_type', 'start_pos', 'end_pos', 'filled_sequence']]

# 保存结果
df_result.to_csv('./5k/head_tail_extensions.csv', index=False)


