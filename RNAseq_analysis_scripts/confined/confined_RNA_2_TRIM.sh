#!/bin/bash

R_DIR=/data/white/miranda/confined_RNAseq
R1="*_R1_001.fastq.gz"
R2="*_R2_001.fastq.gz"

TRIMTAG="_TRIM"
EXTENSION=".fastq.gz"

TRIM=/data/white/software/Trimmomatic-0.36/trimmomatic-0.36.jar
TRIM_FASTA=/data/white/software/Trimmomatic-0.36/adapters/TruNex-PE-1.fa

for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i

        bsub \
                -n4 \
                -R rusage[mem=12] \
		-R "span[ptile=4]" \
                -o confined_RNA_2_TRIM_%J.stdout \
                -eo confined_RNA_2_TRIM_%J.stderr \
                -W 36:00 \
                "java -jar $TRIM \
                        PE \
                        -threads 4 \
                        -phred33 \
                        -trimlog trimPE.log \
                        $R1 $R2 \
                        -baseout $i$TRIMTAG$EXTENSION \
                        ILLUMINACLIP:$TRIM_FASTA:2:30:10:4:true \
                        LEADING:3 \
                        TRAILING:3 \
                        SLIDINGWINDOW:4:15 \
                        MINLEN:36"
done
