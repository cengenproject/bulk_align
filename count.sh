#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=count
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

echo "Starting $(date)"

WS="289"
refdir="/gpfs/ycga/work/hammarlund/aw853/references/WS"$WS/
ref_gtf=$refdir/"c_elegans.PRJNA13758.WS289.canonical_geneset.gtf.gz"
bam_dir="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_bams"



mkdir -p results/counts
module load Subread

featureCounts -a $ref_gtf \
	-o results/counts/231102_bsn12_count \
	-p --countReadPairs \
	-C \
	-t exon \
	-g gene_id \
	-T $SLURM_CPUS_PER_TASK \
	$bam_dir/*.bam



echo "Done $(date)"


