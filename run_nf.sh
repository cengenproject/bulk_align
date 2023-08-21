#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=nf_wrap_bsn11
#SBATCH -c 1
#SBATCH --mem=500M
#SBATCH --time=05:10:00

module load Java/17.0.4

#nextflow run align.nf -with-dag dag10.dot -resume
nextflow run align.nf -resume


