#!/usr/bin/env nextflow



// Define inputs batch directory name and associated suffix for each

/*batch_names = [ '200312_D00306_1303': '',
				'200320_D00306_1305': '',
				'201022_D00306_1322': 't1',
				'201022_D00306_1323': 't2']*/

batch_names = [ 'alexis_test': 'bis']


// Locations and parameters
params.root_dir = '/SAY/standard/mh588-CC1100-MEDGEN/raw_fastq/bulk'
params.WSversion = 'WS277'
params.references_dir = '/gpfs/ycga/project/ysm/hammarlund/aw853/references'

/* ----------------------------------------------------------------------------
              No  modification should be required below this line
---------------------------------------------------------------------------- */


// update when this script is modified
script_version = "bsn5"





// --------    Define references    --------
ref_genome = 'c_elegans.PRJNA13758.'+params.WSversion+'.genomic.fasta'
ref_gtf = 'c_elegans.PRJNA13758.'+params.WSversion+'.canonical_geneset.gtf'

ref_gtf_index = params.references_dir + '/' + params.WSversion + '/' + ref_gtf
STAR_ref_dir = params.references_dir + "/" + params.WSversion + "/star_index_2-7-7a/"






// -------- Prepare channels for input and output --------
batch_dirs = batch_names.keySet().collect { "${params.root_dir}/${it}/Sample_*" }
sample_dir = Channel.fromPath( batch_dirs, type: 'dir' ).map {it ->
//	(it =~ "^/SAY/standard/mh588-CC1100-MEDGEN/raw_fastq/bulk/([0-9_D]+)/Sample_([A-Zef1-9]{2,4}r[0-9]{1,4})")[0]
	(it =~ "^/SAY/standard/mh588-CC1100-MEDGEN/raw_fastq/bulk/(alexis_test)/Sample_([A-Zef12]{2,4}r[0-9]{1,4})")[0]

}


// count nb of Samples in each batch
def batch_sizes = batch_names.keySet().collect {
	def count = 0
	new File("${params.root_dir}/${it}").traverse(nameFilter: ~".*Sample_.*", maxDepth: 0) {
		return count++
	}
	return count
}
assert batch_sizes.size() == batch_names.size()

// repeat suffix as for each sample in the batch
def sample_suffixes_list = []
for( i=0; i<batch_names.size(); i++){
    for(j=0; j<batch_sizes[i]; j++){
        sample_suffixes_list.add(batch_names.values()[i])
    }
}
sample_suffix = Channel.fromList(sample_suffixes_list)



process testEverything{
	
	output:
	val true into tested_successfully
	
	shell:
	'''
	
	# Check that SQL connection can be established.
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 --user=$username --password=$password --database=test_cengen --execute="exit" || exit 1
	
	# Check that NuDup installed and found
	nudup.py --version || exit 1
	
	# Check that references are available
	if [ -d !{STAR_ref_dir} ]
	then
		if [ -e !{STAR_ref_dir}/SA ]
		then
			if [ $(ls !{STAR_ref_dir} | wc -w) -lt 14 ]
			then
				echo "STAR index directory contains too few elements."
				exit 1
			elif [ $(du !{STAR_ref_dir}/SA | cut -f1) -lt $((500 * 1024)) ]
			then
				echo "STAR index smaller than 500 MB."
				exit 1
			fi
		else
			echo "STAR index directory does not contain index as expected."
			exit 1
		fi
	else
		echo "STAR index directory does not exist."
		exit 1
	fi
	
	'''
}


process align{

	cpus '10'
	time '14d'
	memory '15GB'
	
	publishDir '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/'+script_version, mode: 'copy', saveAs: { filename -> sample + sample_suffix + (filename =~ /.*(\.log|\.bam)$/)[0][1] }

	input:
	val tested_successfully
	tuple path(sample_dir), val(batch), val(sample) from sample_dir
	val sample_suffix
	
	
	output:
	path "dedup.sorted.dedup.bam"
	path "sample.log"
	
	shell:
	'''
	echo "##########################################################################################"
	echo "###################                    Operations                   ######################"
	echo "##########################################################################################"
	echo
	echo
	echo "------------------------------------------------------------------------------------------"
	echo "---------------------------------- Merging FASTQ files -----------------------------------"
	echo "------------------------------------------------------------------------------------------"

	# extract sample name (e.g. VCr12)
	sample=!{sample}!{sample_suffix}
	
	
	
	echo "Merging files !{sample_dir}" >> merge.log
	echo >> merge.log

	cat !{sample_dir}/!{sample}_*_R1_00?.fastq.gz > merged_R1.fastq.gz
	cat !{sample_dir}/!{sample}_*_R2_00?.fastq.gz > merged_I1.fastq.gz
	cat !{sample_dir}/!{sample}_*_R3_00?.fastq.gz > merged_R2.fastq.gz

	# Trim the index to 8 bp for NuDup
	zcat merged_I1.fastq.gz | awk 'NR%4 == 1 {print} NR%4 == 2 {print substr($0, 0, 8)} NR%4 == 3 {print} NR%4 == 0 {print substr($0,0,8)}' > trimmed_I1.fq
	
	# Check for problems in nb of reads
	cnt_r1=$(zcat merged_R1.fastq.gz | wc -l)
	cnt_r2=$(zcat merged_R2.fastq.gz | wc -l)
	cnt_i1=$(zcat merged_I1.fastq.gz | wc -l)

	echo "Merged reads: $cnt_r1" >> merge.log
	echo >> merge.log
	echo >> merge.log

	echo "------ Merged reads: $cnt_r1 ------"	
	echo


	echo "------------------------------------------------------------------------------------------"
	echo "------------------------------- BBduk clipping of adaptors -------------------------------"
	echo "------------------------------------------------------------------------------------------"
	#module load BBMap/37.17-foss-2016a-Java-1.8.0_121
	module load BBMap/38.90-GCCcore-10.2.0


	# weirdly, BBduk puts the log in stderr
	bbduk.sh in1=merged_R1.fastq.gz \
		 in2=merged_R2.fastq.gz \
		 out1=trimmed_R1.fastq.gz \
		 out2=trimmed_R2.fastq.gz \
		 ref="/ycga-gpfs/apps/hpc/software/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa" \
		 ktrim=r k=23 mink=11 hdist=1 tpe tbo >> bbduk.log 2>&1

	echo "------------------------------------------------------------------------------------------"
	echo "------------------------------------- STAR alignment -------------------------------------"
	echo "------------------------------------------------------------------------------------------"
	module -q swap BBMap STAR/2.7.7a-GCCcore-10.2.0

	star_version=$(STAR --version)

	echo "Using STAR version $(STAR --version)"
	echo "Using STAR version $(STAR --version)" >> star.log
	
	STAR --runThreadN $SLURM_CPUS_PER_TASK \
		--genomeDir !{STAR_ref_dir} \
		--readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix aligned_ \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard \
		--outFilterMatchNminOverLread 0.3 >> star.log

	echo "------------------------------------------------------------------------------------------"
	echo "------------------------------------- Deduplications -------------------------------------"
	echo "------------------------------------------------------------------------------------------"


	module -q swap STAR SAMtools/1.11-GCCcore-10.2.0
	
	nudup.py --paired-end \
			 -f trimmed_I1.fq \
			 --start 8 \
			 --length 8 \
			 -o dedup \
			 aligned_Aligned.sortedByCoord.out.bam >> nudup.log
	
	




	echo "##########################################################################################"
	echo "###################     Reading quality metrics and full log file   ######################"
	echo "##########################################################################################"
	
	
	# Initialize Full Log
	# Then go through each step, register warnings or errors at the beginning of log file,
	#then dump all diagnostics into the file (also upload to DB), and  finally all other logs from previous softwares.
	
	echo "Finished running alignment !{script_version} of sample $sample on $(date +%F)."
	echo
	echo "     -----     "
	echo
	
	
	
	# Nb of reads after merge 
	
	if [ $cnt_r1 -ne $cnt_r2 -o $cnt_r1 -ne $cnt_i1 ]
	then
		echo "Error: number of reads in the R1/R2/I1 files do not match." >> sample.log
	fi
	
	
	# ----    BBduk    ----
	bbduk_input=$(awk '/^Input:/{print $2}' bbduk.log)
	bbduk_removed=$(awk '/^Total Removed:/{print $3}' bbduk.log)
	bbduk_result=$(awk '/^Result:/{print $2}' bbduk.log)

	if [ $(( $bbduk_input - $bbduk_removed )) -ne $bbduk_result ]
	then
		echo "Warning, BBduk number of reads don't add up for sample $sample" >> sample.log
	fi

	
	if [ $cnt_r1 -ne $(( 2 * bbduk_input )) ]
	then
		echo "Warning, number of reads in merge output and in BBduk input do not match" >> sample.log
	fi

	
	# ----    STAR    ----
	star_finished=$(cat star.log | grep -A5 "..... started STAR run" | grep "..... finished successfully")

	if [ ${#star_finished} -lt 1 ]
	then
		echo "Warning, STAR did not finish succesfully" >> sample.log
	fi

	star_input=$(awk -F "\t" '/Number of input reads/{print $2}' aligned_Log.final.out)
	star_unique_maps=$(awk -F "\t" '/Uniquely mapped reads number/{print $2}' aligned_Log.final.out)
	star_multi_maps=$(awk -F "\t" '/Number of reads mapped to multiple loci/{print $2}' aligned_Log.final.out)

	if [ $bbduk_result -ne $(( 2 * star_input )) ]
	then
		echo "Warning, number of output reads from BBduk and input reads for STAR don't match" >> sample.log
	fi



	# ----    NuDup    ----
	nudup_input=$(tail -1 dedup_dup_log.txt | cut -f1)
	nudup_pos_dup=$(tail -1 dedup_dup_log.txt | cut -f3)
	nudup_umi_dup=$(tail -1 dedup_dup_log.txt | cut -f5)
	
	
	# ----    Final counts    ----
	#	read the index field in the output BAM file
	
	bam_stats=$(samtools flagstat dedup.sorted.dedup.bam)
	
	
	regex="([[:digit:]]+) \\+ 0 in total \\(QC-passed reads \\+ QC-failed reads\\)\\
([[:digit:]]+) \\+ 0 secondary\\
([[:digit:]]+) \\+ 0 supplementary\\
([[:digit:]]+) \\+ 0 duplicates\\
([[:digit:]]+) \\+ 0 mapped \\([[:digit:].]+% : N/A\\)\\
([[:digit:]]+) \\+ 0 paired in sequencing\\
[[:digit:]]+ \\+ 0 read1\\
[[:digit:]]+ \\+ 0 read2\\
([[:digit:]]+) \\+ 0 properly paired \\([[:digit:].]+% : N/A\\)\\
([[:digit:]]+) \\+ 0 with itself and mate mapped\\
0 \\+ 0 singletons \\([[:digit:].]+% : N/A\\)\\
0 \\+ 0 with mate mapped to a different chr\\
0 \\+ 0 with mate mapped to a different chr \\(mapQ>=5\\)"
	
	
	
	
	regex="([[:digit:]]+) \\+ 0 in total \\(QC-passed reads \\+ QC-failed reads\\) ([[:digit:]]+) \\+ 0 secondary ([[:digit:]]+) \\+ 0 supplementary ([[:digit:]]+) \\+ 0 duplicates ([[:digit:]]+) \\+ 0 mapped \\([[:digit:].]+% : N/A\\) ([[:digit:]]+) \\+ 0 paired in sequencing [[:digit:]]+ \\+ 0 read1 [[:digit:]]+ \\+ 0 read2 ([[:digit:]]+) \\+ 0 properly paired \\([[:digit:].]+% : N/A\\) ([[:digit:]]+) \\+ 0 with itself and mate mapped 0 \\+ 0 singletons \\([[:digit:].]+% : N/A\\) 0 \\+ 0 with mate mapped to a different chr 0 \\+ 0 with mate mapped to a different chr \\(mapQ>=5\\)"


	regex="([[:digit:]]+) \\+ 0 in total \\(QC-passed reads \\+ QC-failed reads\\)[\n ]+([[:digit:]]+) \\+ 0 secondary[\n ]+([[:digit:]]+) \\+ 0 supplementary[\n ]+([[:digit:]]+) \\+ 0 duplicates[\n ]+([[:digit:]]+) \\+ 0 mapped \\([[:digit:].]+% : N/A\\)[\n ]+([[:digit:]]+) \\+ 0 paired in sequencing[\n ]+[[:digit:]]+ \\+ 0 read1[\n ]+[[:digit:]]+ \\+ 0 read2[\n ]+([[:digit:]]+) \\+ 0 properly paired \\([[:digit:].]+% : N/A\\)[\n ]+([[:digit:]]+) \\+ 0 with itself and mate mapped[\n ]+0 \\+ 0 singletons \\([[:digit:].]+% : N/A\\)[\n ]+0 \\+ 0 with mate mapped to a different chr[\n ]+0 \\+ 0 with mate mapped to a different chr \\(mapQ>=5\\)"

	if [[ $bam_stats =~ $regex ]]
	then
		final_total=${BASH_REMATCH[1]}
		final_secondary=${BASH_REMATCH[2]}
		final_supplementary=${BASH_REMATCH[3]}
		final_duplicates=${BASH_REMATCH[4]}
		final_mapped=${BASH_REMATCH[5]}
		final_paired=${BASH_REMATCH[6]}
		final_properly_paired=${BASH_REMATCH[7]}
	else
		echo "ERROR -- Could not read flagstat correctly. Final statistics incorrect" >> sample.log
		echo  			>> sample.log
		echo "Content:" >> sample.log
		echo $bam_stats >> sample.log
		echo  			>> sample.log
	fi
	
	echo                                                       >> sample.log
	echo "Quality metrics"                                     >> sample.log
	echo                                                       >> sample.log
	echo "Merged reads: $cnt_r1"                               >> sample.log
	echo "  --"                                                >> sample.log
	echo "BBduk input: $bbduk_input"                           >> sample.log
	echo "BBduk removed: $bbduk_removed"                       >> sample.log
	echo "BBduk result: $bbduk_result"                         >> sample.log
	echo "  --"                                                >> sample.log
	echo "STAR input: $star_input"                             >> sample.log
	echo "STAR mapped:\t"$(($star_unique_maps+$star_multi_maps)) >>sample.log
	echo "STAR unique maps: $star_unique_maps"                 >> sample.log
	echo "STAR multi maps: $star_multi_maps"                   >> sample.log
	echo "  --"                                                >> sample.log
	echo "NuDup input: $nudup_input"                           >> sample.log
	echo "NuDup positional duplicates: $nudup_pos_dup"         >> sample.log
	echo "NuDup UMI duplicates: $nudup_umi_dup"                >> sample.log
	echo "  --"                                                >> sample.log
	echo "Final file total: $final_total"                      >> sample.log
	echo "Final file secondary: $final_secondary"              >> sample.log
	echo "Final file supplementary: $final_supplementary"      >> sample.log
	echo "Final file duplicates: $final_duplicates"            >> sample.log
	echo "Final file mapped: $final_mapped"                    >> sample.log
	echo "Final file paired: $final_paired"                    >> sample.log
	echo "Final file properly paired: $final_properly_paired"  >> sample.log
	echo 													   >> sample.log
	echo "##################################################"  >> sample.log
	echo                                                       >> sample.log
	

	# Fill log file
	echo "    ------------    Merge    ------------" >> sample.log
	cat merge.log >> sample.log
	echo "    ------------    BBduk    ------------" >> sample.log
	cat bbduk.log >> sample.log
	echo "    ------------    STAR    ------------" >> sample.log
	cat star.log >> sample.log


	echo "Creating database entry" >> sample.log

	

	if [ $(( nudup_input - nudup_umi_dup )) -gt 5000000 ]
	then
		quality_check=0
	else
		quality_check=1
	fi


	sql_command="INSERT INTO alig2 (alig_lib_id, alig_pid, alig_sequ_batch, alig_date, alig_lab, alig_by, alig_line, alig_neurons, alig_reference_type, alig_reference_source, alig_reference_version, alig_software, alig_software_version, alig_aligned_reads, alig_quality_check,
	alig_merged_reads, alig_bbduk_input, alig_bbduk_removed, alig_bbduk_result,
	alig_star_input,alig_star_unique, alig_star_multi, alig_nudup_input, alig_nudup_pos_dup, alig_nudup_umi_dup,
	alig_final_total, alig_final_secondary, alig_final_supplementary, alig_final_duplicates, alig_final_mapped, alig_final_paired, alig_final_properly_paired,
	alig_comments)
	VALUES ('"$sample"', '!{script_version}', '!{batch}', '"$(date +%F)"', 'MH', '"$(echo $USER)"', '', '', 'genome', 'Wormbase', '!{params.WSversion}', 'STAR', '"$star_version"', '"$(($star_unique_maps+$star_multi_maps))"', '', '"$cnt_r1"', '"$bbduk_input"', '"$bbduk_removed"', '"$bbduk_result"', '"$star_input"', '"$star_unique_maps"', '"$star_multi_maps"', '"$nudup_input"', '"$nudup_pos_dup"', '"$nudup_umi_dup"', '"$final_total"', '"$final_secondary"', '"$final_supplementary"', '"$final_duplicates"', '"$final_mapped"', '"$final_paired"', '"$final_properly_paired"', '');"
	
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=test_cengen \
		--execute="$sql_command"

	echo "Uploaded to database." >> sample.log
	
	echo "Finished."
	'''
}



