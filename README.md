# R570_Sugarcane_Genome_Upgrade
Code for the paper: Upgrading the R570 sugarcane genome with ultra-long Nanopore reads improves resistance-gene discovery
# Installation

# Usage
Prior to running the pipeline, the resequencing data must be assembled into contigs. Flye, minimap2 and hifiasm are recommended for this purpose.
```bash
python 0_seqdata.py
```
The extracted sequences need to be aligned with the reference genome using BLAST.

```bash
python 1_gap_filling.py
python 1_head_tail_filling.py
python 1_result_merge.py
python 2_update_files.py
