#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
R_DIR=/data/white/miranda/confined_RNAseq
RT1="*TRIM_1P.fastq.gz"
RT2="*TRIM_2P.fastq.gz"

FASTQC=/data/white/software/FastQC/fastqc

# *******************************************************************************

# submit bsub for each sample
# FASTQC post-trim

for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i

        bsub \
                -n1 \
                -R rusage[mem=24] \
                -o confined_RNA_3_FASTQC2_%J.stdout \
                -eo confined_RNA_3_FASTQC2_%J.stderr \
                -W 36:00 \
                "$FASTQC $RT1 $RT2"
done

