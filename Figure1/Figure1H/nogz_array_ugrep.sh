#!/bin/bash

module load gcc/11.2.0 rust/1.66.1

SEQS=$1 #use input to seqs file
FILES=($(cat $2))
outfolder=$3
ERR=$4 # Use input into call.
indexfile=$5

IFS=$'\n'

# Take unique file identifier (SRR, etc.) and name output to match metadata (Customized to filenames for this iteration)
name=$(echo ${FILES[$SLURM_ARRAY_TASK_ID]} | sed -nr "s/.*\/(.*\.fastq)/\1/p")

# ugrep call which finds reads matching the query sequences. Assigns read names for processing with the AGREPAssess.R script
if [ $ERR == 0 ]
	then
		ug -z -H --format='> %G_%a%_seq%m%~%O%~' -f $SEQS ${FILES[$SLURM_ARRAY_TASK_ID]} > ${outfolder}/${name}_temp.fa
	else
		ug -Z${ERR} -z -H --format='> %G_%a%_seq%m%~%O%~' -f $SEQS ${FILES[$SLURM_ARRAY_TASK_ID]} > ${outfolder}/${name}_temp.fa

fi

# Adds RT/TS classifier separated by _
awk 'NR==FNR{ar["> "NR"_"]=$1 "_" $2} NR>FNR{match($0,/^> [0-9]+_/,m); gsub(/^> [0-9]+_/, "> " ar[m[0]]"_", $0);print} ' $indexfile ${outfolder}/${name}_temp.fa > ${outfolder}/${name}.fa

# awk 'NR==FNR{ar["> "NR"_"]=$0} NR>FNR{match($0,/^> [0-9]+_/,m); gsub(/^> [0-9]+_/, "> " ar[m[0]]"_", $0);print} ' SequenceNames.txt ${outfolder}/${name}_temp.fa > ${outfolder}/${name}.fa

rm ${outfolder}/${name}_temp.fa

