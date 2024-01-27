#!/usr/bin/bash

scp /Users/onkiwong/Desktop/Year_3/common_core_transdisciplinary_project/exploratory_test_2/gisaid_hcov_19.fasta angel@hpc02.sbms.hku.hk:/home/angel/data/hcov_19

# mafft --auto input > output
"/home/angel/miniconda3/envs/conda_env/bin/mafft"  --auto --reorder "gisaid_hcov_19.fasta" > "gisaid_hcov_19_aligned.fasta"

############################################

# Biopython for summarising alignment results (Platform: Google Colab)

!pip install biopython
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import SeqIO

##### nucleotides sequence #####

hcov_19_aligned = AlignIO.read('gisaid_hcov_19_aligned.fasta', 'fasta')

for seq in hcov_19_aligned:
    print(seq)

display = 0
i = 0
while display < 10 and i < 10:
    characters = hcov_19_aligned[:, i]
    print(characters)
    display += 1
    i += 1
    
print(hcov_19_aligned)

length = hcov_19_aligned.get_alignment_length()
diff_pos = {}
i = 0

while i < length:
    characters = hcov_19_aligned[:, i]
    characters = [nt for nt in characters if nt != '-']
    
    if len(set(characters)) > 1:
        diff_pos[i] = characters
        
    i += 1

print(diff_pos)

for k,v in diff_pos.items():
  print(f'position: {k}: {v}')


import csv

with open('nt_mismatch.csv', 'w', newline = '') as f:

  # Create a csv writer object
  writer = csv.writer(f)

  # Write the header row
  writer.writerow(['Position', 'Nucleotides variations'])
  # Write one key-value tuple per row
  for k, v in diff_pos.items():
    writer.writerow([k, v])


##### Amino Acids sequence #####

hcov_19_translated = AlignIO.read('gisaid_hcov_19_aligned.translated.fasta', 'fasta')

for seq in hcov_19_translated:
    print(seq)

display = 0
i = 0
while display < 10 and i < 10:
    characters = hcov_19_translated[:, i]
    print(characters)
    display += 1
    i += 1
    
print(hcov_19_translated)

length_aa = hcov_19_translated.get_alignment_length()
diff_pos_aa = {}
i = 0

while i < length_aa:
    residue = hcov_19_translated[:, i]
    residue = [aa for aa in residue if aa != '-']
    
    if len(set(residue)) > 1:
        diff_pos_aa[i] = residue
        
    i += 1

print(diff_pos_aa)

for k,v in diff_pos_aa.items():
  print(f'position: {k}: {v}')

import csv

with open('AA_mismatch.csv', 'w', newline = '') as f:

  # Create a csv writer object
  writer = csv.writer(f)

  # Write the header row
  writer.writerow(['Position', 'Amino Acids variations'])
  # Write one key-value tuple per row
  for k, v in diff_pos_aa.items():
    writer.writerow([k, v])

