# FYP
Script for FYP submission

Aligner pipelines script contains the codes for the different aligners used

Stringie (in order):
  1. Stringtie - Aligner.sh
  2. ballgown.r
  3. Stringtie-last.sh


CrypSplice:

run following command
python CrypSplice.py –C WT1.bed,WT2.bed,WT3.bed –T MT1.bed,MT2.bed,MT3.bed –G MM10 –F 10 –M 0.95 –P 6
  
Dependencies:
Python version 2.7
bedtools version 2.29.2

CrypSplice pipeline requires CrypSplice.py and BBTestP.r to be in the same directory
