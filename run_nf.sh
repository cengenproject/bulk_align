#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=nf_wrap_bsn9
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=21-05:10:00

nextflow run align.nf -with-dag dag9.dot -resume


