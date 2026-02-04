# R570_Sugarcane_Genome_Upgrade
Code for the paper: Upgrading the R570 sugarcane genome with ultra-long Nanopore reads improves resistance-gene discovery
## Installation
### Getting started
```bash
git clone https://github.com/wangge-1996/R570_Sugarcane_Genome_Upgrade.git
cd R570_Sugarcane_Genome_Upgrade
```
### Install Dependencies
```yaml
dependencies:
  - python=3.9
  - biopython>=1.84
  - pandas>=2.2.3
  - numpy>=1.26.4
  - bedtools>=2.31.1
  - openpyxl>=3.1.5
  - pytables>=3.9.2
  - pip>=24.3.1
```
## Usage
Prior to running the pipeline, the resequencing data must be assembled into contigs. Flye, minimap2 and hifiasm are recommended for this purpose.
```bash
# Step 1: Extract sequences from reference genome
python 0_seqdata.py
```
The extracted sequences need to be aligned with the reference genome using BLAST.

```bash
# Step 2: Fill genome gaps and extend head/tail sequences
python 1_gap_filling.py
python 1_head_tail_filling.py
# Step 3: Merge analysis results
python 1_result_merge.py
# Step 4: Update reference genome and annotation files
python 2_update_files.py
```
## Citation
