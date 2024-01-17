#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
R_DIR=/data/white/miranda/HMGB2_OE/RNAseq
R1="*TRIM_1P.fastq.gz"
R2="*TRIM_2P.fastq.gz"

STAR=/data/white/software/STAR/STAR-2.5.3a/bin/Linux_x86_64/STAR

# Human genome from Nate
GTF=/data/white/genomes/Human/GRCh38_E90/annotation/Homo_sapiens.GRCh38.90.gtf
STAR_IDX=/data/white/genomes/Human/GRCh38_E90/STAR_Index/

# *******************************************************************************

# submit bsub for each sample
# STAR alignment

for i in $(cat $R_DIR/Sample_List.txt); do
	echo "Running $i"
	cd $R_DIR/$i

	bsub \
		-n8 \
		-R rusage[mem=8] \
		-o HMG_RNA_4_STAR_%J.stdout \
		-eo HMG_RNA_4_STAR_%J.stderr \
		-R "span[ptile=8]" \
		-W 12:00 \
		"$STAR \
			--runThreadN 8 \
			--genomeDir $STAR_IDX \
			--outSAMtype BAM SortedByCoordinate \
			--readFilesCommand gunzip -c \
			--readFilesIn $R1 $R2 \
			--quantMode GeneCounts \
			--sjdbGTFfile $GTF"
done



