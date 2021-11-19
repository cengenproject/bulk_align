#!/usr/bin/env nextflow



// Define inputs batch directory name and associated suffix for each

batch_names = [ '200312_D00306_1303': '',
		'200320_D00306_1305': '',
		'201022_D00306_1322': 't1',
		'201022_D00306_1323': 't2',
		'210108_D00306_1332': '',
		'210622_D00306_1363': 't1',
		'210622_D00306_1364': 't2',
		'211108_D00306_1377': 't1',
		'211103_D00306_1375': 't2']

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
script_version = "bsn9"





// --------    Define references    --------
ref_genome = 'c_elegans.PRJNA13758.'+params.WSversion+'.genomic.fa'
ref_gtf = 'c_elegans.PRJNA13758.'+params.WSversion+'.canonical_geneset.gtf'

ref_gtf_index = params.references_dir + '/' + params.WSversion + '/' + ref_gtf
WS_ref_dir = params.references_dir + "/" + params.WSversion
STAR_index_dir = WS_ref_dir + "/star_index_2-7-7a/"





// -------- Prepare channels for input and output --------
batch_dirs = batch_names.keySet().collect { "${params.root_dir}/${it}/Sample_*" }
samples_dir = Channel.fromPath( batch_dirs, type: 'dir' )
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

samples_suffix = Channel.fromList(sample_suffixes_list)



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
	
	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3
	
	input:
	  val tested_successfully
	  tuple path(sample_dir), val(batch), val(sample) from samples_dir
	  val sample_suffix from samples_suffix
	
	output:
	  env ready_to_process
	  tuple path(sample_dir), val(batch), val(sample), val(sample_suffix) into samples_pinit
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	
	source ~/.cengen_database.id
	
	
	# First check if sample was already succesfully processed
	sql_command="SELECT alig_completed FROM alig
		WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
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
		
		echo "New sample, will process: !{sample}!{sample_suffix}, !{script_version}, !{batch}"
		
		sql_command="INSERT INTO alig
			SET alig_lib_id		= '!{sample}!{sample_suffix}',
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
	elif [ "$existing" -eq 0 ]
	then
		ready_to_process=1
		echo "Sample already exists, but unfinished. Will reprocess: !{sample}!{sample_suffix}, !{script_version}, !{batch}"
	else
		ready_to_process=0
		echo "Sample already processed (completed=$existing ), will be ignored: !{sample}!{sample_suffix}, !{script_version}, !{batch}"
	fi
	
	'''
}


process merge{
	cpus '1'
	time '2h'
	memory '15GB'

	input:
	  val ready_to_process
	  tuple path(sample_dir), val(batch), val(sample), val(sample_suffix) from samples_pinit
	
	when:
		ready_to_process == '1'
	
	output:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) into samples_pmerge
	  tuple path('merge.log'), env(cnt_r1) into logs_pmerge
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	echo "---------------------------------- Merging FASTQ files -----------------------------------"

	
	echo "Merging files !{sample_dir}" >> merge.log
	echo >> merge.log

	cat !{sample_dir}/!{sample}_*_R1_00?.fastq.gz > merged_R1.fastq.gz
	cat !{sample_dir}/!{sample}_*_R2_00?.fastq.gz > merged_I1.fastq.gz
	cat !{sample_dir}/!{sample}_*_R3_00?.fastq.gz > merged_R2.fastq.gz

	# Trim the index to 8 bp for NuDup
	zcat merged_I1.fastq.gz | awk 'NR%4 == 1 {print} NR%4 == 2 {print substr($0, 0, 8)} NR%4 == 3 {print} NR%4 == 0 {print substr($0,0,8)}' > trimmed_I1.fq
	
	

	# Check nb of reads after merge 
	cnt_r1=$(zcat merged_R1.fastq.gz | wc -l)
	cnt_r2=$(zcat merged_R2.fastq.gz | wc -l)
	cnt_i1=$(cat trimmed_I1.fq | wc -l)
	
	if [ $cnt_r1 -ne $cnt_r2 -o $cnt_r1 -ne $cnt_i1 ]
	then
		echo "Error: number of reads in the R1/R2/I1 files do not match." >> merge.log
		echo "R1: $cnt_r1" >> merge.log
		echo "R2: $cnt_r2" >> merge.log
		echo "I1: $cnt_i1" >> merge.log
		exit 1
	fi
	
	echo "Merged reads: $cnt_r1" >> merge.log
	echo >> merge.log
	echo >> merge.log
	
	'''
}



process qc{
	cpus '2'
	time '2h'
	memory '10GB'

	module 'FastQC'
	
	input:
      tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from samples_pmerge
	  tuple path('merge.log'), val(merge_count) from logs_pmerge

	output:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) into samples_pqc
	  tuple path('merge.log'), val(merge_count), path('merged_R1_fastqc/summary.txt'), path('merged_R2_fastqc/summary.txt') into logs_pqc
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
	touch !{sample}!{sample_suffix}
	
	fastqc --extract --threads $SLURM_CPUS_PER_TASK merged_R*
	'''	

}

process trim{
	cpus 1
	time '2h'
	memory '10GB'
	module 'BBMap/38.90-GCCcore-10.2.0'
	
	input:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from samples_pqc
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc') from logs_pqc
	
	
	output:
	  tuple path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) into samples_ptrim
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc'), path("bbduk.log") into logs_ptrim
	
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	# weirdly, BBduk puts the log in stderr
	
	bbduk.sh in1=merged_R1.fastq.gz \
		 in2=merged_R2.fastq.gz \
		 out1=trimmed_R1.fastq.gz \
		 out2=trimmed_R2.fastq.gz \
		 ref="/ycga-gpfs/apps/hpc/software/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa" \
		 ktrim=r k=23 mink=11 hdist=1 tpe tbo >> bbduk.log 2>&1
	'''
}

process align_first{
	cpus 10
	time '5h'
	memory '15GB'
	module 'STAR/2.7.7a-GCCcore-10.2.0'
	
	
	input:
	  tuple path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from samples_ptrim
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc'), path("bbduk.log") from logs_ptrim
	
	output:
	  tuple path("aligned_Aligned.sortedByCoord.out.bam"), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) into samples_palign_first
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc'), path("bbduk.log"), path("aligned_Log.final.out"), path("aligned_Log.out"), path("star.log") into logs_palign_first
	  path("aligned_SJ.out.tab")
	  path("aligned_Log*")
	
	
	publishDir mode: 'copy',
		pattern: '*_SJ.out.tab',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_junctions/',
		saveAs: { sample + sample_suffix + '.SJ.tab' }
	
	
	publishDir mode: 'copy',
		pattern: '*_Log.final.out',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '_STAR_stats.log'}
	
	publishDir mode: 'copy',
		pattern: '*_Log.out',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '_STAR.log'}
	
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
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
	
	'''	
}

process dedup{
	cpus 1
	time '3h'
	memory '10GB'
	module 'SAMtools/1.11-GCCcore-10.2.0'
	
	errorStrategy 'retry'
	maxRetries 5
	
	
	input:
	  tuple path("aligned.bam"), path('trimmed_I1.fq'), val(sample), val(batch), val(sample_suffix) from samples_palign_first
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc'), path("bbduk.log"), path("aligned_Log.final.out"), path("aligned_Log.out"), path("star.log") from logs_palign_first
	
	output:
	  tuple path("dedup.sorted.dedup.bam"), val(sample), val(batch), val(sample_suffix) into samples_pdedup
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc'), path("bbduk.log"), path("aligned_Log.final.out"), path("aligned_Log.out"), path("star.log"), path("nudup.log"),     path("dedup_dup_log.txt") into logs_pdedup
	  path "dedup.sorted.dedup.bam*"
	
	
	publishDir mode: 'copy',
		pattern: '*.bam',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_bams/',
		saveAs: { sample + sample_suffix + '.bam' }
	
	publishDir mode: 'copy',
                pattern: '*.bam.bai',
                path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_bams/',
                saveAs: { sample + sample_suffix + '.bam.bai' }
	
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	nudup.py --paired-end \
			 -f trimmed_I1.fq \
			 --start 8 \
			 --length 8 \
			 -o dedup \
			 aligned.bam >> nudup.log
			 
	if [[ $(cat dedup_dup_log.txt | cut -f1 | tail -1) -eq 0 ]]
	then
		echo "Empty input, NudDup may have done something strange, retrying."
		exit 1
	fi
	
	samtools index -@ $SLURM_CPUS_PER_TASK dedup.sorted.dedup.bam
	'''
}




process save_logs{
	cpus '1'
	time '1h'
	memory '1GB'
	
	errorStrategy 'retry'
	maxRetries 3

	module 'SAMtools'
	
	input:
	  tuple path("dedup.bam"), val(sample), val(batch), val(sample_suffix) from samples_pdedup
	  tuple path('merge.log'), val(merge_count), path('R1_fastqc'), path('R2_fastqc'), path("bbduk.log"), path("aligned_Log.final.out"), path("aligned_Log.out"), path("star.log"), path("nudup.log"),     path("dedup_dup_log.txt") from logs_pdedup
	
	output:
	  path "sample.log"
	
	publishDir mode: 'copy',
		pattern: 'sample.log',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '.log'}

	
	  

	shell:
	'''
	
	# Initialize Full Log
	# Then go through each step, register warnings or errors at the beginning of log file,
	#then dump all diagnostics into the file (also upload to DB), and  finally all other logs from previous softwares.
	
	sample=!{sample}!{sample_suffix}
	
	echo "Finished running alignment !{script_version} of sample $sample on $(date +%F)." > sample.log
	
	
	
	
	# register warnings
	alig_comments=""
	
	
	# ---- FASTQC ----
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
	
	# ---- BBduk ----
	bbduk_input=$(awk '/^Input:/{print $2}' bbduk.log)
	bbduk_removed=$(awk '/^Total Removed:/{print $3}' bbduk.log)
	bbduk_result=$(awk '/^Result:/{print $2}' bbduk.log)

	if [ $(( $bbduk_input - $bbduk_removed )) -ne $bbduk_result ]
	then
		alig_comments="$alig_comments -- BBduk number of reads don't add up"
	fi

	
	if [ !{merge_count} -ne $(( 2 * bbduk_input )) ]
	then
		alig_comments="$alig_comments -- number of reads in merge output and in BBduk input do not match"
	fi

	
	# ----    STAR    ----
	star_finished=$(cat star.log | grep -A5 "..... started STAR run" | grep "..... finished successfully")

	if [ ${#star_finished} -lt 1 ]
	then
		alig_comments="$alig_comments -- STAR did not finish successfully"
	fi

	star_input=$(awk -F "\t" '/Number of input reads/{print $2}' aligned_Log.final.out)
	star_unique_maps=$(awk -F "\t" '/Uniquely mapped reads number/{print $2}' aligned_Log.final.out)
	star_multi_maps=$(awk -F "\t" '/Number of reads mapped to multiple loci/{print $2}' aligned_Log.final.out)
	star_unal_too_many_loci=$(awk -F "\t" '/Number of reads mapped to too many loci/{print $2}' aligned_Log.final.out)
	star_unal_many_mismatches=$(awk -F "\t" '/Number of reads unmapped: too many mismatches/{print $2}' aligned_Log.final.out)
	star_unal_too_short=$(awk -F "\t" '/Number of reads unmapped: too short/{print $2}' aligned_Log.final.out)
	star_unal_other=$(awk -F "\t" '/Number of reads unmapped: other/{print $2}' aligned_Log.final.out)

	if [ $bbduk_result -ne $(( 2 * star_input )) ]
	then
		alig_comments="$alig_comments -- number of output reads from BBduk and input reads for STAR do not match"
	fi
	
	
	if [ $star_input -ne $(( $star_unique_maps + $star_multi_maps + $star_unal_too_many_loci + $star_unal_many_mismatches +$star_unal_too_short + $star_unal_other )) ]
	then
		alig_comments="$alig_comments -- number of STAR reads do not match"
	fi



	# ----    NuDup    ----
	nudup_input=$(tail -1 dedup_dup_log.txt | cut -f1)
	nudup_unaligned=$(tail -1 dedup_dup_log.txt | cut -f2)
	nudup_pos_dup=$(tail -1 dedup_dup_log.txt | cut -f3)
	nudup_umi_dup=$(tail -1 dedup_dup_log.txt | cut -f5)
	
	if [ $nudup_unaligned -ne $(( 2 * ($star_unal_many_mismatches + $star_unal_too_many_loci + $star_unal_too_short + $star_unal_other) )) ]
	then
		alig_comments="$alig_comments -- number of unaligned reads do not match between STAR and NuDup"
	fi
	
	
	
	# ----    Final counts    ----
	
	bam_stats=$(samtools flagstat dedup.bam)
	
	
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
		alig_comments="$alig_comments -- failed to read flagstats"
	fi
	
	echo                                                       >> sample.log
	echo "$alig_comments"                                      >> sample.log
	echo                                                       >> sample.log
	echo  "-------------"                                      >> sample.log
	echo                                                       >> sample.log
	echo "Quality metrics"                                     >> sample.log
	echo                                                       >> sample.log
	echo "Merged reads: !{merge_count}"                               >> sample.log
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
	cat merge.log                                    >> sample.log
	echo "    ------------  FASTQC R1  ------------" >> sample.log
	cat R1_fastqc                                    >> sample.log
	echo "    ------------  FASTQC R2  ------------" >> sample.log
	cat R2_fastqc                                    >> sample.log
	echo "    ------------    BBduk    ------------" >> sample.log
	cat bbduk.log                                    >> sample.log
	echo "    ------------     STAR    ------------" >> sample.log
	cat star.log                                     >> sample.log
	echo "    ------------    NuDup    ------------" >> sample.log
	cat nudup.log                                    >> sample.log
	echo                                             >> sample.log
	echo dedup_dup_log.txt                           >> sample.log
	echo                                             >> sample.log
	echo                                             >> sample.log





	echo "Creating database entry" >> sample.log

	
	
	



	
	
	
	sql_command="UPDATE alig
	  SET	alig_completed			= '1',
		alig_aligned_reads	 		= '"$(($star_unique_maps+$star_multi_maps))"',
		alig_quality_check			= '0',
		alig_merged_reads			= '!{merge_count}',
		alig_fqc_basic_stats		= '$(worst_of ${arr1[0]} ${arr2[0]})',
		alig_fqc_per_base_quality	= '$(worst_of ${arr1[1]} ${arr2[1]})',
		alig_fqc_per_tile_quality	= '$(worst_of ${arr1[2]} ${arr2[2]})',
		alig_fqc_per_seq_scores		= '$(worst_of ${arr1[3]} ${arr2[3]})',
		alig_fqc_per_base_seq_content	= '$(worst_of ${arr1[4]} ${arr2[4]})',
		alig_fqc_per_base_gc		= '$(worst_of ${arr1[5]} ${arr2[5]})',
		alig_fqc_per_base_n			= '$(worst_of ${arr1[6]} ${arr2[6]})',
		alig_fqc_seq_length			= '$(worst_of ${arr1[7]} ${arr2[7]})',
		alig_fqc_duplication		= '$(worst_of ${arr1[8]} ${arr2[8]})',
		alig_fqc_overrepresented_seq	= '$(worst_of ${arr1[9]} ${arr2[9]})',
		alig_fqc_adapter_content	= '$(worst_of ${arr1[10]} ${arr2[10]})',
		alig_bbduk_input			= '"$bbduk_input"',
		alig_bbduk_removed			= '"$bbduk_removed"',
		alig_bbduk_result			= '"$bbduk_result"',
		alig_star_input				= '"$star_input"',
		alig_star_unique			= '"$star_unique_maps"',
		alig_star_multi				= '"$star_multi_maps"',
		alig_nudup_input			= '"$nudup_input"',
		alig_nudup_pos_dup			= '"$nudup_pos_dup"',
		alig_nudup_umi_dup			= '"$nudup_umi_dup"',
		alig_final_total			= '"$final_total"',
		alig_final_secondary		= '"$final_secondary"',
		alig_final_supplementary	= '"$final_supplementary"',
		alig_final_duplicates		= '"$final_duplicates"',
		alig_final_mapped			= '"$final_mapped"',
		alig_final_paired			= '"$final_paired"',
		alig_final_properly_paired	= '"$final_properly_paired"',
		alig_comments				= '"$alig_comments"'
	  WHERE alig_lib_id		= '"$sample"' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
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





