#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
R_DIR=/data/white/miranda/confined_RNAseq
READS=ReadsPerGene.out.tab

OUTDIR='/Users/hunterm/Downloads/'

# *******************************************************************************

# make folder and copy files

for i in $(cat $R_DIR/Sample_List.txt); do
	echo "Copying $i"
	rsync -va --progress --no-p lilac:$R_DIR/$i/$READS "$OUTDIR/$I"
done
