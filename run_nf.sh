#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=nf_wrap_bsn11
#SBATCH -c 1
#SBATCH --mem=500M
#SBATCH --time=05:10:00

module load Java/17.0.4
module load miniconda # for nudup and fastq-screen

#nextflow run align.nf -with-dag dag10.dot -resume
nextflow run align.nf -resume -with-dag dag12.dot


