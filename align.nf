#!/usr/bin/env nextflow



// Define inputs batch directory name and associated suffix for each

batch_names = [ '200312_D00306_1303': '',
		'200320_D00306_1305': '',
		'201022_D00306_1322': 't1',
		'201022_D00306_1323': 't2',
		'210108_D00306_1332': '',
		'210622_D00306_1363': 't1',
		'210622_D00306_1364': 't2']

// test samples
//batch_names = [ '000000_D00000_0000': 'bis']


// Locations and parameters
params.root_dir = '/SAY/standard/mh588-CC1100-MEDGEN/raw_fastq/bulk'
params.WSversion = 'WS281'
params.references_dir = '/gpfs/ycga/project/ysm/hammarlund/aw853/references'

/* ----------------------------------------------------------------------------
              No  modification should be required below this line
---------------------------------------------------------------------------- */


// update when this script is modified
script_version = "bsn7"





// --------    Define references    --------
ref_genome = 'c_elegans.PRJNA13758.'+params.WSversion+'.genomic.fa'
ref_gtf = 'c_elegans.PRJNA13758.'+params.WSversion+'.canonical_geneset.gtf'

ref_gtf_index = params.references_dir + '/' + params.WSversion + '/' + ref_gtf
WS_ref_dir = params.references_dir + "/" + params.WSversion
STAR_index_dir = WS_ref_dir + "/star_index_2-7-7a/"





// -------- Prepare channels for input and output --------
batch_dirs = batch_names.keySet().collect { "${params.root_dir}/${it}/Sample_*" }
sample_dir_for_init = Channel.fromPath( batch_dirs, type: 'dir' )
	.map {it ->
		(it =~ "^/SAY/standard/mh588-CC1100-MEDGEN/raw_fastq/bulk/([0-9_D]+)/Sample_([A-Zef1-9]{2,4}r[0-9]{1,4})")[0]}


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

sample_suffix_for_init = Channel.fromList(sample_suffixes_list)



process testEverything{
	
	output:
	val true into tested_successfully

	module 'R'
	module 'STAR'
	
	shell:
	'''
	
	# Check that SQL connection can be established.
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 --user=$username --password=$password --database=test_cengen --execute="exit" || exit 1
	
	# Check that NuDup installed and found
	nudup.py --version || exit 1
	
	# Check that references are available
	if [ ! -d !{STAR_index_dir} ]
	then
		echo "STAR index directory does not exist, downloading GTF and genome with wbData."
		R -e 'wbData::wb_get_gtf_path("!{params.WSversion}", "!{WS_ref_dir}")'
		zcat !{WS_ref_dir}/!{ref_gtf}.gz > !{WS_ref_dir}/!{ref_gtf}
		R -e 'wbData::wb_get_genome_path("!{params.WSversion}", "!{WS_ref_dir}")'
		zcat !{WS_ref_dir}/!{ref_genome}.gz > !{WS_ref_dir}/!{ref_genome}
		echo
		echo "Creating index"
		STAR --runThreadN $SLURM_CPUS_PER_TASK \
			--runMode genomeGenerate \
			--genomeDir !{STAR_index_dir} \
			--genomeFastaFiles !{WS_ref_dir}/!{ref_genome} \
			--sjdbGTFfile !{WS_ref_dir}/!{ref_gtf} \
			--sjdbOverhang 149 \
			--genomeSAindexNbases 12
	fi
	if [ ! -f !{STAR_index_dir}/SA ]
	then
		echo "STAR index directory does not contain index as expected."
		exit 1
	else
		if [ $(ls !{STAR_index_dir} | wc -w) -lt 14 ]
		then
			echo "STAR index directory contains too few elements."
			exit 1
		elif [ $(du !{STAR_index_dir}/SA | cut -f1) -lt $((500 * 1024)) ]
		then
			echo "STAR index smaller than 500 MB."
			exit 1
		fi
	fi
	'''
}

process initialize_db_entry{
	cpus '1'
	time '10m'
	memory '100MB'
	cache false
	
	input:
	  val tested_successfully
	  tuple path(sample_dir), val(batch), val(sample) from sample_dir_for_init
	  val sample_suffix_for_init
	
	output:
	  env ready_to_process
	  tuple path(sample_dir), val(batch), val(sample), val(sample_suffix) into sample_dir_for_merge
	
	shell:
	'''
	
	source ~/.cengen_database.id
	
	
	# First check if sample was already succesfully processed
	sql_command="SELECT alig_completed FROM alig
		WHERE alig_lib_id		= '!{sample}!{sample_suffix_for_init}' AND
		    alig_pid		= '!{script_version}' AND
		    alig_sequ_batch	= '!{batch}';"
			
			
	existing=$(mysql --host=23.229.215.131 \
					--user=$username \
					--password=$password \
					--database=cengen \
					--execute="$sql_command" \
					--silent --raw)
	
	
	
	if [[ -z "$existing" ]]
	then
		
		ready_to_process=1
		
		echo "New sample, will process: !{sample}!{sample_suffix_for_init}, !{script_version}, !{batch}"
		
		sql_command="INSERT INTO alig
			SET alig_lib_id		= '!{sample}!{sample_suffix_for_init}',
				alig_pid		= '!{script_version}',
				alig_sequ_batch	= '!{batch}',
				alig_completed	= '0',
				alig_date		= '"$(date +%F)"',
				alig_lab		= 'MH',
				alig_by		= '"$(echo $USER)"',
				alig_reference_type = 'genome',
				alig_reference_source = 'Wormbase',
				alig_reference_version = '!{params.WSversion}',
				alig_software	   = 'STAR',
				alig_software_version  = '2.7.7a';"

		mysql --host=23.229.215.131 \
			--user=$username \
			--password=$password \
			--database=cengen \
			--execute="$sql_command"
		
	elif [[ $(echo "$existing" | wc -l) -gt 1 ]]
	then
		echo "Potential error: several entries already exist!"
		exit 1
	elif [[ $existing == 0]]
	then
		ready_to_process=1
		echo "Sample already exists, but unfinished. Will reprocess: !{sample}!{sample_suffix_for_init}, !{script_version}, !{batch}"
	else
		ready_to_process=0
		echo "Sample already processed (completed=$existing ), will be ignored: !{sample}!{sample_suffix_for_init}, !{script_version}, !{batch}"
	fi
	
	'''
}


process merge{
	cpus '1'
	time '2h'
	memory '15GB'

	input:
	  val ready_to_process
	  tuple path(sample_dir), val(batch), val(sample), val(sample_suffix) from sample_dir_for_merge
	
	when:
		ready_to_process
	
	output:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) into merged_paths
	  path('merge.log') into merge_log
	
	shell:
	'''
	echo "---------------------------------- Merging FASTQ files -----------------------------------"

	
	echo "Merging files !{sample_dir}" >> merge.log
	echo >> merge.log

	cat !{sample_dir}/!{sample}_*_R1_00?.fastq.gz > merged_R1.fastq.gz
	cat !{sample_dir}/!{sample}_*_R2_00?.fastq.gz > merged_I1.fastq.gz
	cat !{sample_dir}/!{sample}_*_R3_00?.fastq.gz > merged_R2.fastq.gz

	# Trim the index to 8 bp for NuDup
	zcat merged_I1.fastq.gz | awk 'NR%4 == 1 {print} NR%4 == 2 {print substr($0, 0, 8)} NR%4 == 3 {print} NR%4 == 0 {print substr($0,0,8)}' > trimmed_I1.fq
	
	cnt_r1=$(zcat merged_R1.fastq.gz | wc -l)

	echo "Merged reads: $cnt_r1" >> merge.log
	echo >> merge.log
	echo >> merge.log

	echo "------ Merged reads: $cnt_r1 ------"	
	echo
	'''
}


merged_paths.into{ merged_paths_to_align; merged_paths_for_qc; merged_paths_for_logs }

process qc{
	cpus '2'
	time '2h'
	memory '10GB'

	module 'FastQC'
	
	input:
      tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from merged_paths_for_qc

	output:
      tuple path('merged_R1_fastqc/summary.txt'), path('merged_R2_fastqc/summary.txt') into fastqc_summaries
	  path('merged_R1_fastqc.html')
	  path('merged_R2_fastqc.html')

	publishDir mode: 'copy',
			pattern: '*.html',
			path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/'+script_version+'_logs/',
			saveAs: {filename ->
					read_file = (filename =~ /merged_R([12])_fastqc/)[0][1]
					sample+sample_suffix+'_R'+read_file+'_fastqc.html'}

	shell:
        '''
        fastqc --extract --threads $SLURM_CPUS_PER_TASK merged_R*
        '''	

}
	

process align{

	cpus '10'
	time '14d'
	memory '15GB'
	
	errorStrategy 'ignore'
	
	
	input:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from merged_paths_to_align

	publishDir mode: 'copy',
		pattern: '*.bam',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_bams/',
		saveAs: { sample + sample_suffix + '.bam' }
	
	publishDir mode: 'copy',
                pattern: '*_SJ.out.tab',
                path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_junctions/',
                saveAs: { sample + sample_suffix + '.SJ.tab' }
	
	publishDir mode: 'copy',
                pattern: '*.bam.bai',
                path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_bams/',
                saveAs: { sample + sample_suffix + '.bam.bai' }
	
	publishDir mode: 'copy',
		pattern: '*_Log.final.out',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '_STAR_stats.log'}
	
	publishDir mode: 'copy',
		pattern: '*_Log.out',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '_STAR.log'}
	
	
	output:
	  path "aligned_Log.final.out" into star_align_log_stats
	  path "aligned_Log.out"
	  tuple path("bbduk.log"), path("star.log"), path("nudup.log"), path("dedup_dup_log.txt") into alignment_logs
	  path "dedup.sorted.dedup.bam" into aligned_bam
	  path "aligned_SJ.out.tab"
	
	shell:
	'''

	echo "------------------------------- BBduk clipping of adaptors -------------------------------"
	module load BBMap/38.90-GCCcore-10.2.0


	# weirdly, BBduk puts the log in stderr
	bbduk.sh in1=merged_R1.fastq.gz \
		 in2=merged_R2.fastq.gz \
		 out1=trimmed_R1.fastq.gz \
		 out2=trimmed_R2.fastq.gz \
		 ref="/ycga-gpfs/apps/hpc/software/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa" \
		 ktrim=r k=23 mink=11 hdist=1 tpe tbo >> bbduk.log 2>&1

	echo "------------------------------------- STAR alignment -------------------------------------"
	module -q swap BBMap STAR/2.7.7a-GCCcore-10.2.0


	echo "Using STAR version $(STAR --version)"
	echo "Using STAR version $(STAR --version)" >> star.log
	
	STAR --runThreadN $SLURM_CPUS_PER_TASK \
		--genomeDir !{STAR_index_dir} \
		--readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix aligned_ \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard \
		--outFilterMatchNminOverLread 0.3 >> star.log

	echo "------------------------------------- Deduplications -------------------------------------"

	module -q swap STAR SAMtools/1.11-GCCcore-10.2.0
	
	nudup.py --paired-end \
			 -f trimmed_I1.fq \
			 --start 8 \
			 --length 8 \
			 -o dedup \
			 aligned_Aligned.sortedByCoord.out.bam >> nudup.log

	echo "------------------------------------- BAM indexing --------------------------------------"
	samtools index -@ $SLURM_CPUS_PER_TASK dedup.sorted.dedup.bam
	
	'''
}

process save_logs{
	
	cpus '1'
	time '1d'
	memory '1GB'
	
	module 'SAMtools'
	
	input:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from merged_paths_for_logs
	  path('merge.log') from merge_log
	  tuple path('R1_fastqc'), path('R2_fastqc') from fastqc_summaries
	  tuple path("bbduk.log"), path("star.log"), path("nudup.log"), path("dedup_dup_log.txt") from alignment_logs
	  path "dedup.sorted.dedup.bam" from aligned_bam
	  path "aligned_Log.final.out" from star_align_log_stats
	
	publishDir mode: 'copy',
		pattern: 'sample.log',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '.log'}

	output:
	  path "sample.log"
	  

	shell:
	'''
	echo "##########################################################################################"
	echo "###################     Reading quality metrics and full log file   ######################"
	echo "##########################################################################################"
	
	
	# Initialize Full Log
	# Then go through each step, register warnings or errors at the beginning of log file,
	#then dump all diagnostics into the file (also upload to DB), and  finally all other logs from previous softwares.
	
	sample=!{sample}!{sample_suffix}
	
	echo "Finished running alignment !{script_version} of sample $sample on $(date +%F)."
	echo
	echo "     -----     "
	echo
	
	
	
	# Nb of reads after merge 
	cnt_r1=$(zcat merged_R1.fastq.gz | wc -l)
	cnt_r2=$(zcat merged_R2.fastq.gz | wc -l)
	cnt_i1=$(cat trimmed_I1.fq | wc -l)
	
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
		echo "Content:" 	>> sample.log
		echo $bam_stats 	>> sample.log
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

	
	# Get FASTQC results
	arr1=($(cat R1_fastqc | cut -f1))
	arr2=($(cat R2_fastqc | cut -f1))

	worst_of (){
		if [ $1 == "FAIL"  ] || [ $2 == "FAIL" ]
		then
			echo "FAIL"
		elif [ $1 == "WARN" ] || [ $2 == "WARN" ]
		then
			echo "WARN"
		else
			echo "PASS"
		fi
	}



	comments=""
	
	
	sql_command="UPDATE alig
	  SET	alig_completed			= '1',
		alig_aligned_reads	 	= '"$(($star_unique_maps+$star_multi_maps))"',
		alig_merged_reads		= '"$cnt_r1"',
		alig_fqc_basic_stats		= '$(worst_of ${arr1[0]} ${arr2[0]})',
		alig_fqc_per_base_quality	= '$(worst_of ${arr1[1]} ${arr2[1]})',
		alig_fqc_per_tile_quality	= '$(worst_of ${arr1[2]} ${arr2[2]})',
		alig_fqc_per_seq_scores		= '$(worst_of ${arr1[3]} ${arr2[3]})',
		alig_fqc_per_base_seq_content	= '$(worst_of ${arr1[4]} ${arr2[4]})',
		alig_fqc_per_base_gc		= '$(worst_of ${arr1[5]} ${arr2[5]})',
		alig_fqc_per_base_n		= '$(worst_of ${arr1[6]} ${arr2[6]})',
		alig_fqc_seq_length		= '$(worst_of ${arr1[7]} ${arr2[7]})',
		alig_fqc_duplication		= '$(worst_of ${arr1[8]} ${arr2[8]})',
		alig_fqc_overrepresented_seq	= '$(worst_of ${arr1[9]} ${arr2[9]})',
		alig_fqc_adapter_content	= '$(worst_of ${arr1[10]} ${arr2[10]})',
		alig_bbduk_input		= '"$bbduk_input"',
		alig_bbduk_removed		= '"$bbduk_removed"',
		alig_bbduk_result		= '"$bbduk_result"',
		alig_star_input			= '"$star_input"',
		alig_star_unique		= '"$star_unique_maps"',
		alig_star_multi			= '"$star_multi_maps"',
		alig_nudup_input		= '"$nudup_input"',
		alig_nudup_pos_dup		= '"$nudup_pos_dup"',
		alig_nudup_umi_dup		= '"$nudup_umi_dup"',
		alig_final_total		= '"$final_total"',
		alig_final_secondary		= '"$final_secondary"',
		alig_final_supplementary	= '"$final_supplementary"',
		alig_final_duplicates		= '"$final_duplicates"',
		alig_final_mapped		= '"$final_mapped"',
		alig_final_paired		= '"$final_paired"',
		alig_final_properly_paired	= '"$final_properly_paired"',
		alig_comments			= '"$comments"'
	  WHERE alig_lib_id     = '"$sample"' AND
		alig_pid        = '!{script_version}' AND
		alig_sequ_batch = '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=cengen \
		--execute="$sql_command"

	echo "Uploaded to database." >> sample.log
	
	echo "Finished."
	'''
}





