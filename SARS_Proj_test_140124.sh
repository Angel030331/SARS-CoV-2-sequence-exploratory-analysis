##### 14/01/2024 SARS_Proj_test #####

# 1. linux: using MAFFT to align
# Environmental setting
# /home/angel
    # data
        # SARS_project_test
            # raw_fasta
                # patient_sequence.fasta
                # ref_sequence.fasta
    # SARS_proj_test_140124

# 1. Pairwise/Multiple alignment with MAFFT

# codes
# test: 200, 300, 500, 1000 iterations

mkdir /home/angel/data/SARS_project_test/align
cd /home/angel/data/SARS_project_test/raw_fasta
mafft --localpair --maxiterate 200 --add patient_sequence.fasta ref_sequence.fasta > /home/angel/data/SARS_project_test/align/output_it_200.aln
mafft --localpair --maxiterate 300 --add patient_sequence.fasta ref_sequence.fasta > /home/angel/data/SARS_project_test/align/output_it_300.aln
mafft --localpair --maxiterate 500 --add patient_sequence.fasta ref_sequence.fasta > /home/angel/data/SARS_project_test/align/output_it_500.aln
mafft --localpair --maxiterate 1000 --add patient_sequence.fasta ref_sequence.fasta > /home/angel/data/SARS_project_test/align/output_it_1000.aln

# export
scp -r angel@hpc02.sbms.hku.hk:/home/angel/data/SARS_project_test/align /Users/onkiwong/Downloads

# 2. Biopython for downstream analysis
# jupyter notebook

# Environmental setting
conda install anaconda::biopython
conda info --env
print(os.getcwd())
for f in os.listdir("/Users/onkiwong/Desktop/Year_3/common_core_transdisciplinary_project/exploratory_test/data/align"):
    print(f)

# Packages
import os
import Bio

# Analysis
from Bio import AlignIO
file_list = ['output_it_200.aln','output_it_300.aln','output_it_500.aln','output_it_1000.aln']
for f in file_list:
    align = AlignIO.read(f, "fasta")
    print(align)

# Compare differences in genome sequences
from Bio import SeqIO

def seq_info(seq):
    print("ID: ", seq.id)
    print("Length: ", len(seq))
    sequence = seq.seq
    print()
    return sequence

def search_gap(sequence):
    gap = 0
    pos_gap = []
    for pos, nt in enumerate(sequence):
        if nt == '-':
            gap += 1
            pos_gap.append(pos)
    return gap, pos_gap

def SNP(seq_list):
    length = len(seq_list[0])
    pos_diff = []
    for pos in range(length):
        nt_pos = [s[pos] for s in seq_list]
        if all(nt == nt_pos[0] for nt in nt_pos):
            continue
        else:
            pos_diff.append(pos)           
    return len(pos_diff)


seq_list = []
align_200 = SeqIO.parse('output_it_200.aln', 'fasta')
for seq in align_200:
    sequence = seq_info(seq)
    seq_list.append(sequence)
    print(f'gap(count, position): {search_gap(sequence)}')
    print()
    print(f'single nucleotide differences in genome (including gap): {SNP(seq_list)}')
    print()
    

# Translate to amino acid sequence

from Bio import SeqIO
from Bio.Seq import Seq

align_200 = list(SeqIO.parse('output_it_200.aln', 'fasta'))
AA_seq_list = []

def seq_info(seq):
    print("ID: ", seq.id)
    print("Length: ", len(seq))
    sequence = Seq(str(seq.seq))
    AA_seq = sequence.ungap('-').translate()
    print(AA_seq)
    return sequence, AA_seq

def AA_diff(AA_seq_list):
    length = min(len(seq) for seq in AA_seq_list)
    pos_diff = []
    for pos in range(length):
        aa_pos = [s[pos] for s in AA_seq_list]
        if all(aa == aa_pos[0] for aa in aa_pos):
            continue
        else:
            pos_diff.append(pos)           
    return len(pos_diff)

for seq in align_200:
    sequence, AA_seq = seq_info(seq)
    AA_seq_list.append(AA_seq)
    print()
    print(f'AA sequence differences: {AA_diff(AA_seq_list)}')



