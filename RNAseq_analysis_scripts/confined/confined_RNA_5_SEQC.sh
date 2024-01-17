#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
R_DIR=/data/white/miranda/confined_RNAseq
BAMIN=Aligned.sortedByCoord.out.bam
BAMRG=Aligned.sortedByCoord.RG.out.bam

GTF_TID=/data/white/genomes/Human/GRCh38_E90/annotation/Homo_sapiens.GRCh38.90.TID.gtf

REF_LOC=/data/white/genomes/Human/GRCh38_E90/sequence/
REF=Homo_sapiens.GRCh38.dna.primary_assembly.fa

JAVA_1_8=/programs/x86_64-linux/java/jdk1.8.0_144/bin/java
PICARD=/data/white/software/Picard/picard.jar
SEQC=/data/white/software/RNA-SeQC/RNA-SeQC_v1.1.8.jar
SAMTOOLS=/data/white/software/samtools-1.6/samtools


# *******************************************************************************

# submit bsub for each sample
# SeQC prep and run

for i in $(cat $R_DIR/Sample_List.txt); do
	echo "Running $i"
	cd $R_DIR/$i

	bsub \
		-n1 \
		-R rusage[mem=18] \
		-o confined_RNA_5_SeQC_%J.stdout \
		-eo confined_RNA_5_SeQC_%J.stderr \
		-W 6:00 \
		"$JAVA_1_8 -jar $PICARD AddOrReplaceReadGroups I=$BAMIN O=$BAMRG LB=MVH PL=illumina PU=MVH SM=MVH

		$SAMTOOLS index -b $BAMRG

		java -Xmx16G -jar $SEQC \
                        -s \"$i|$BAMRG|confined_RNA\" \
                        -t $GTF_TID \
                        -r $REF_LOC/$REF \
                        -o SeQC"
done

