#!/bin/bash

# Alignment with bowtie2
# To be completed after trimming with TrimGalore/cutadapt

# *******************************************************************************

# Experiment Parameters:
R_DIR=/data/white/miranda/HMGB2_OE/Miranda_analyses/ATACseq/
R1="*_R1_001_val_1.fq.gz"
R2="*_R2_001_val_2.fq.gz"

samtools=/opt/common/CentOS_7/samtools/samtools-1.13/bin/samtools
v2=sample_aligned.bam

numthreads=8

# *******************************************************************************

for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i

        bsub \
              	-n8 \
                -R rusage[mem=8] \
                -o HMG_ATAC_3_sort_%J.stdout \
                -eo HMG_ATAC_3_sort_%J.stderr \
                -R "span[ptile=8]" \
                -W 24:00 \
                "$samtools sort -@ $numthreads -o sample_aligned_sorted.bam $v2"
done
