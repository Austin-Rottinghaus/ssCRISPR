# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 15:28:44 2022

@author: austr
"""

import gzip
from Bio import SeqIO

sequence_MG='GATGTGTACTGACTGAATGGCGCTGGCCTCACCGCCCGGAACCGGAACTTCGGTGCCAATGACATAGCTCAGTTGCTCACGCTGGCAATCTGTCGCCACACTTTCCGCAGCAAAGCAAAGCACAGCAGCTCGTTCCGCAACCGTTTCTGGTGCTAACGGTATGGGATCCCCCGCGCAGGACATTGACGCATCAAGATGAATTTTACTGAAGCCGGCACGAACATATTCCTTTACCAGCTC'
sequence_Nissle='GATGTGTACTGACTGAATGGCGCTGGCCTCACCGCCCGGAACCGGAACTTCGGTGCCAATGACATAGTTCAGTTGCTCACGCTGGCAATCAGTCGCCACACTTTCCGCCGCCAGACAAAGCACAGCGGCACGTTCAGCAACCGTTTCTGGCGCTAACGGTATGGAATCGTCAGCGCAGGACATTGACGCATCAAGATGAATTTTACTGAAGCCGGCACGAACATATGCCTTTACCAGCTC'
sequence_10B='GATGTGTACTGACTGAATGGCGCTGGCCTCACCGCCCGGAACCGGAACTTCGGTGCCAATGACATAGCTCAGTTGCTCACGCTGGCAATCTGTCGCCACACTTTCCGCAGCAAAGCAAAGCACAGCAGCTCGTTCCGCAACCGTTTCTGGTGCTAACGGTATGGGATCCCCCGCGCAGGACATTGACGCATCAAGATGAATTTTACTGAAGCCGGCACGAACATATTCCTATACCAGCTC'
sequence_BL='GATGTGTACTGACTGAATGGCGCTGGCCTCACCGCCCGGAACCGGAACTTCGGTACCAATGACATAGCTCAGTTGCTCACGCTGGCAATCTGTCGCCACACTTTCCGCAGCAAAGCAAAGCACAGCAGCTCGTTCCGCAACCGTTTCTGGTGCTAACGGTATGGGATCCCCCGCGCAGGACATTGACGCATCAAGATGAATTTTACTGAAGCCAGCACGAACATATGCTTTTACCAGCTC'

number_total=0
number_MG=0
number_Nissle=0
number_10B=0
number_BL=0
with gzip.open("620-1_R2_001.fastq.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        sequence=record.seq
        for n in range(8):
            if sequence[n]=='G' and sequence[n+1]=='A' and sequence[n+2]=='T' and sequence[n+3]=='G':
                sequence_new=sequence[n:n+240]
                number_total=number_total+1
                if sequence_new==sequence_MG:
                    number_MG=number_MG+1
                if sequence_new==sequence_Nissle:
                    number_Nissle=number_Nissle+1
                if sequence_new==sequence_10B:
                    number_10B=number_10B+1
                if sequence_new==sequence_BL:
                    number_BL=number_BL+1
                continue
print(number_total)
print(number_MG+number_Nissle+number_10B+number_BL)
print(number_10B)
print(number_Nissle)
print(number_MG)
print(number_BL)