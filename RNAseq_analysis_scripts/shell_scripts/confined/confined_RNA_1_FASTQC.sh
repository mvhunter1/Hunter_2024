#!/bin/bash

R_DIR=/data/white/miranda/confined_RNAseq
R1="*_R1_001.fastq.gz"
R2="*_R2_001.fastq.gz"

FASTQC=/data/white/software/FastQC/fastqc

for i in $(cat $R_DIR/Sample_List.txt); do
echo "Running $i"
cd $R_DIR/$i

bsub \
-n1 \
-R rusage[mem=24] \
-o confined_RNA_1_FASTQC_%J.stdout \
-eo confined_RNA_1_FASTQC_%J.stderr \
-W 36:00 \
"$FASTQC $R1 $R2"
done

