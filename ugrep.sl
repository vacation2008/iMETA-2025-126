#!/bin/bash

#SBATCH --account=hlilab
#SBATCH --partition=standard
#SBATCH --time=6-12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH -J PCA_ug
#SBATCH -o outputs/ug_out_%a
#SBATCH -e outputs/ug_err_%a

seqs=$1
files=$2
outfolder=$3
error=$4
indexfile=$5

bash nogz_array_ugrep.sh $seqs $files $outfolder $error $indexfile $SLURM_ARRAY_TASK_ID
