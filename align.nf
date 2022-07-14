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
params.WSversion = 'WS284'
params.references_dir = '/gpfs/ycga/project/ysm/hammarlund/aw853/references'
params.fastq_screen_conf = '/home/aw853/project/references/2205a_fastq_screen.conf'

params.db_name = 'cengen'
params.db_table_name = 'alig_11'

/* ----------------------------------------------------------------------------
              No  modification should be required below this line
---------------------------------------------------------------------------- */


// update when this script is modified
script_version = "bsn11"





// --------    Define references    --------
ref_genome = 'c_elegans.PRJNA13758.'+params.WSversion+'.genomic.fa'
ref_gtf = 'c_elegans.PRJNA13758.'+params.WSversion+'.canonical_geneset.gtf'

ref_gtf_index = params.references_dir + '/' + params.WSversion + '/' + ref_gtf
WS_ref_dir = params.references_dir + "/" + params.WSversion

STAR_version = '2.7.9a'
STAR_index_dir = WS_ref_dir + '/' + STAR_version + '/'





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

	module 'R/4.2.0-foss-2020b'
	module 'STAR/'+STAR_version+'-GCCcore-10.2.0'
	
	shell:
	'''
	
	# Check that SQL connection can be established.
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 --user=$username --password=$password --database=!{params.db_name} --execute="exit" || exit 1
	
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
	  tuple path(sample_dir), path("sample.log"), val(batch), val(sample), val(sample_suffix) into ch_init
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	# Init text log
	echo "Starting sample !{sample}!{sample_suffix} from batch !{batch} with !{script_version} on $(date +%F)." > sample.log

	
	# Check if sample was already succesfully processed
	source ~/.cengen_database.id
	sql_command="SELECT alig_completed FROM !{params.db_table_name} WHERE alig_lib_id = '!{sample}!{sample_suffix}' AND alig_pid = '!{script_version}' AND alig_sequ_batch = '!{batch}';"
	echo "Executing SQL command: $sql_command"		
			
	existing=$(mysql --host=23.229.215.131 --user=$username --password=$password --database=!{params.db_name} --execute="$sql_command" --silent --raw)
	
echo "done"
	
	
	if [[ -z "$existing" ]]
	then
		ready_to_process=1
		
		echo "New sample, will process: !{sample}!{sample_suffix}, !{script_version}, !{batch}"
		echo "New sample, will process." >> sample.log
		
		sql_command="INSERT INTO !{params.db_table_name}
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
				alig_software_version  = '!{STAR_version}';"

		mysql --host=23.229.215.131 \
			--user=$username \
			--password=$password \
			--database=!{params.db_name} \
			--execute="$sql_command"
		
	elif [[ $(echo "$existing" | wc -l) -gt 1 ]]
	then
		echo "Potential error: several entries already exist!"
		exit 1
	elif [ "$existing" -eq 0 ]
	then
		ready_to_process=1
		echo "Sample already exists, but unfinished. Will reprocess: !{sample}!{sample_suffix}, !{script_version}, !{batch}"
		echo "Sample already exists, but unfinished. Will reprocess." >> sample.log
	else
		ready_to_process=0
		echo "Sample already processed (completed=$existing ), will be ignored: !{sample}!{sample_suffix}, !{script_version}, !{batch}"
		echo "Sample already processed (completed=$existing ), will be ignored." >> sample.log
	fi
	
	echo >> sample.log
		
	'''
}


process merge{
	cpus '1'
	time '2h'
	memory '15GB'

	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3
	
	input:
	  val ready_to_process
	  tuple path(sample_dir), path("sample.log"), val(batch), val(sample), val(sample_suffix) from ch_init
	
	when:
		ready_to_process == '1'
	
	output:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'),
		path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_merge
	
	shell:
	'''
	# create empty file useful for debugging purpose
	touch !{sample}!{sample_suffix}
	
	echo "---------------------------------- Merging FASTQ files -----------------------------------"

	
	echo "----- Merging files !{sample_dir} -----" >> sample.log
	echo >> sample.log

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
		echo "Error: number of reads in the R1/R2/I1 files do not match." >> sample.log
		echo "R1: $cnt_r1" >> sample.log
		echo "R2: $cnt_r2" >> sample.log
		echo "I1: $cnt_i1" >> sample.log
		exit 1
	fi
	
	echo "Merged reads: $cnt_r1" >> sample.log
	echo >> sample.log	
	
	## Update SQL database
	sql_command="UPDATE !{params.db_table_name}
	  SET alig_merged_reads = '$cnt_r1'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "Uploaded to database." >> sample.log
	echo >> sample.log
	
	'''
}



process fastqc_raw{
	cpus '2'
	time '2h'
	memory '10GB'

	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3

	module 'FastQC/0.11.9-Java-11'
	
	input:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_merge

	output:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_fastqc_raw
	  path("*.html")

	publishDir mode: 'copy',
			pattern: '*.html',
			path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/'+script_version+'_logs/',
			saveAs: {filename ->
					read_file = (filename =~ /merged_R([12])_fastqc/)[0][1]
					sample+sample_suffix+'_fastqc_raw_R'+read_file+'.html'}

	shell:
	'''
	# create empty file useful for debugging purpose
	touch !{sample}!{sample_suffix}
	
	echo "----- FastQC raw (pre-alignment) -----"
	echo "----- FASTQC raw (pre-alignment) -----" >> sample.log
	
	
	fastqc --extract --threads $SLURM_CPUS_PER_TASK merged_R*
	
	
	echo "    ~~~~~~~~~~~~  FASTQC raw R1  "    >> sample.log
	cat merged_R1_fastqc/summary.txt    >> sample.log
	echo "    ~~~~~~~~~~~~  FASTQC raw R2  "    >> sample.log
	cat merged_R2_fastqc/summary.txt    >> sample.log
	
	
	# Process results: read summary file as a bash array
	
	arr1=($(cat merged_R1_fastqc/summary.txt | cut -f1))
	arr2=($(cat merged_R2_fastqc/summary.txt | cut -f1))

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
	
	
	sql_command="UPDATE !{params.db_table_name}
	  SET	alig_fqcraw_basic_stats		 = '$(worst_of ${arr1[0]} ${arr2[0]})',
		alig_fqcraw_per_base_quality	 = '$(worst_of ${arr1[1]} ${arr2[1]})',
		alig_fqcraw_per_tile_quality	 = '$(worst_of ${arr1[2]} ${arr2[2]})',
		alig_fqcraw_per_seq_scores		 = '$(worst_of ${arr1[3]} ${arr2[3]})',
		alig_fqcraw_per_base_seq_content = '$(worst_of ${arr1[4]} ${arr2[4]})',
		alig_fqcraw_per_base_gc		     = '$(worst_of ${arr1[5]} ${arr2[5]})',
		alig_fqcraw_per_base_n			 = '$(worst_of ${arr1[6]} ${arr2[6]})',
		alig_fqcraw_seq_length			 = '$(worst_of ${arr1[7]} ${arr2[7]})',
		alig_fqcraw_duplication		     = '$(worst_of ${arr1[8]} ${arr2[8]})',
		alig_fqcraw_overrepresented_seq	 = '$(worst_of ${arr1[9]} ${arr2[9]})',
		alig_fqcraw_adapter_content	     = '$(worst_of ${arr1[10]} ${arr2[10]})'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "FastQC-raw uploaded to database." >> sample.log
	echo >> sample.log
	
	'''	
}

process fastq_screen{
	cpus '6'
	time '2h'
	memory '10GB'

	module 'FastQ_Screen/0.14.1-foss-2018b-Perl-5.28.0'
	
	input:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_fastqc_raw

	output:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_fastqscreen
	  tuple path("merged_R1_screen.txt"), path("merged_R2_screen.txt"), val(sample), val(batch), val(sample_suffix) into summaries_fastqscreen
	  path("*_screen.*")

	publishDir mode: 'copy',
			pattern: '*_screen_.*',
			path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/'+script_version+'_logs/',
			saveAs: {filename ->
					println("Before"+filename)
					read_file = (filename =~ /merged_R([12])_screen\.(html|png|txt)/)[0]
					println("res "+read_file)
					sample+sample_suffix+'_fastqscreen_R'+read_file[1]+'.'+read_file[2]}

	shell:
	'''
	# create empty file useful for debugging purpose
	touch !{sample}!{sample_suffix}
	
	
	fastq_screen \
	--aligner bowtie2 \
	--conf !{params.fastq_screen_conf} \
	--subset 100000 \
	--threads $SLURM_CPUS_PER_TASK \
	merged_R*.fastq.gz
	
	echo "-----  FastQ-Screen -----" >> sample.log
	echo "   ~~~~~ R1 screen " >> sample.log
	cat merged_R1_screen.txt >> sample.log
	echo "   ~~~~~ R2 screen " >> sample.log
	cat merged_R2_screen.txt >> sample.log
	echo >> sample.log
	
	'''
}

process db_fastqscreen{
	cpus '1'
	time '10min'
	memory '500MB'

	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3
	
	module 'R/4.2.0-foss-2020b'
	
	input:
	  tuple path("merged_R1_screen.txt"), path("merged_R2_screen.txt"), val(sample), val(batch), val(sample_suffix) from summaries_fastqscreen
	  
	shell:
	'''
	#!/usr/bin/env Rscript
	
	library(readr)
	
	
	parse_fqscr_summary <- function(in_path){
		
		fqs_out_str <- read_lines(in_path)
		
		fqs_out <- readr::read_tsv(in_path,
							   skip = 1,
							   n_max = length(fqs_out_str) - 4,
							   col_types = cols(
								 Genome = col_character(),
								 `#Reads_processed` = col_integer(),
								 `#Unmapped` = col_integer(),
								 `%Unmapped` = col_double(),
								 `#One_hit_one_genome` = col_integer(),
								 `%One_hit_one_genome` = col_double(),
								 `#Multiple_hits_one_genome` = col_integer(),
								 `%Multiple_hits_one_genome` = col_double(),
								 `#One_hit_multiple_genomes` = col_integer(),
								 `%One_hit_multiple_genomes` = col_double(),
								 Multiple_hits_multiple_genomes = col_integer(),
								 `%Multiple_hits_multiple_genomes` = col_double()
							   ))
		
		uniq_worm <- fqs_out$`%One_hit_one_genome`[fqs_out$Genome == "Worm"]
		
		any_worm <- fqs_out$`%One_hit_one_genome`[fqs_out$Genome == "Worm"] +
		  fqs_out$`%Multiple_hits_one_genome`[fqs_out$Genome == "Worm"] +
		  fqs_out$`%One_hit_multiple_genomes`[fqs_out$Genome == "Worm"] +
		  fqs_out$`%Multiple_hits_multiple_genomes`[fqs_out$Genome == "Worm"]
		
		is_nonworm <- setdiff(fqs_out$Genome,
		                      c("Worm", "Mitochondria","rRNA","Vectors","Adapters"))
		
		any_other <- fqs_out$`%One_hit_one_genome`[fqs_out$Genome %in% is_nonworm] +
		  fqs_out$`%Multiple_hits_one_genome`[fqs_out$Genome %in% is_nonworm] +
		  fqs_out$`%One_hit_multiple_genomes`[fqs_out$Genome %in% is_nonworm] +
		  fqs_out$`%Multiple_hits_multiple_genomes`[fqs_out$Genome %in% is_nonworm]
		
		no_match <- stringr::str_match(string = fqs_out_str[length(fqs_out_str)],
									   pattern = "^%Hit_no_genomes: ([0-9.]{3,6})$")[,2]
		
		mitochondria <- fqs_out$`%One_hit_one_genome`[fqs_out$Genome == "Mitochondria"] +
		  fqs_out$`%Multiple_hits_one_genome`[fqs_out$Genome == "Mitochondria"] +
		  fqs_out$`%One_hit_multiple_genomes`[fqs_out$Genome == "Mitochondria"] +
		  fqs_out$`%Multiple_hits_multiple_genomes`[fqs_out$Genome == "Mitochondria"]
		
		
		rrna <- fqs_out$`%One_hit_one_genome`[fqs_out$Genome == "rRNA"] +
		  fqs_out$`%Multiple_hits_one_genome`[fqs_out$Genome == "rRNA"] +
		  fqs_out$`%One_hit_multiple_genomes`[fqs_out$Genome == "rRNA"] +
		  fqs_out$`%Multiple_hits_multiple_genomes`[fqs_out$Genome == "rRNA"]
		
		list(uniq_worm = uniq_worm, any_worm = any_worm, any_other = max(any_other), no_match = no_match, mitochondria = mitochondria, rrna = rrna)
	}
	
	cat("Updating DB with FastQ-Screen results for !{sample}!{sample_suffix}\n")
	
	fqscr_R1 <- parse_fqscr_summary("merged_R1_screen.txt")
	fqscr_R2 <- parse_fqscr_summary("merged_R2_screen.txt")
	

	sql_command <- paste0("UPDATE !{params.db_table_name}
	  SET	alig_fqscr_uniq_worm  = '", min(fqscr_R1[["uniq_worm"]], fqscr_R2$uniq_worm), "',
		alig_fqscr_any_worm	      = '", min(fqscr_R1$any_worm, fqscr_R2$any_worm), "',
		alig_fqscr_any_other	  = '", max(fqscr_R1$any_other, fqscr_R2$any_other), "',
		alig_fqscr_no_match		  = '", max(fqscr_R1$no_match, fqscr_R2$no_match), "',
		alig_fqscr_mitochondria	  = '", max(fqscr_R1$mitochondria, fqscr_R2$mitochondria), "',
		alig_fqscr_rrna		      = '", max(fqscr_R1$rrna, fqscr_R2$rrna), "'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';")
	
	source("~/.cengen_database.id")
	
	con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   host = "23.229.215.131",
                   user = username,
                   password = password,
                   dbname = "!{params.db_name}")
	
	DBI::dbGetQuery(con, sql_command)
	DBI::dbDisconnect(con)
  
	cat("Database updated with FastQ-Screen summary.\n")
	'''
}


process trim{
	cpus 1
	time '2h'
	memory '10GB'

	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3

	module 'BBMap/38.90-GCCcore-10.2.0'
	
	input:
	  tuple path('merged_R1.fastq.gz'), path('merged_R2.fastq.gz'), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_fastqscreen
	
	output:
	  tuple path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_trim
	
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	echo "----- Trimming with BBduk -----"
	
	# weirdly, BBduk puts the log in stderr
	
	bbduk.sh in1=merged_R1.fastq.gz \
		 in2=merged_R2.fastq.gz \
		 out1=trimmed_R1.fastq.gz \
		 out2=trimmed_R2.fastq.gz \
		 ref="/ycga-gpfs/apps/hpc/software/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa" \
		 ktrim=r k=23 mink=11 hdist=1 tpe tbo >> bbduk.log 2>&1
	
	echo "----  BBduk trim  ----" >> sample.log
	cat bbduk.log >> sample.log
	echo >> sample.log
	
	# Extract summary values
	bbduk_input=$(awk '/^Input:/{print $2}' bbduk.log)
	bbduk_removed=$(awk '/^Total Removed:/{print $3}' bbduk.log)
	bbduk_result=$(awk '/^Result:/{print $2}' bbduk.log)

	if [ $(( $bbduk_input - $bbduk_removed )) -ne $bbduk_result ]
	then
		echo "WARNING: BBduk number of reads do not add up!!"    >> sample.log
		echo "bbduk_input=$bbduk_input"                          >> sample.log
		echo "bbduk_removed=$bbduk_removed"                      >> sample.log
		echo "bbduk_result=$bbduk_result (should be equal to input - removed)"    >> sample.log
	fi
	
	
	sql_command="UPDATE !{params.db_table_name}
	  SET	alig_bbduk_input		= '"$bbduk_input"',
		alig_bbduk_removed			= '"$bbduk_removed"',
		alig_bbduk_result			= '"$bbduk_result"'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "Uploaded BBduk trim summary to database." >> sample.log
	echo >> sample.log
	
	'''
}

process align{
	cpus 10
	time '5h'
	memory '15GB'

	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3
	
	module 'STAR/'+ STAR_version +'-GCCcore-10.2.0'
	
	input:
	  tuple path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_trim
	
	output:
	  tuple path("aligned_Aligned.sortedByCoord.out.bam"), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_align
	  path("*_SJ.out.tab")
	
	publishDir mode: 'copy',
		pattern: '*_SJ.out.tab',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_junctions/',
		saveAs: { sample + sample_suffix + '.SJ.tab' }
	
	
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
	
	
	echo "-----   STAR alignment   -----" >> sample.log
	echo "   ~~~~~~  STAR output log "    >> sample.log
	cat star.log                          >> sample.log
	echo "   ~~~~~~  STAR self log "      >> sample.log
	cat aligned_Log.out                   >> sample.log
	echo "   ~~~~~~  STAR final stats "   >> sample.log
	cat aligned_Log.final.out             >> sample.log
	
	
	# Extract summaries
	
	star_finished=$(cat star.log | grep -A5 "..... started STAR run" | grep "..... finished successfully")
	
	star_input=$(awk -F "\t" '/Number of input reads/{print $2}' aligned_Log.final.out)
	star_unique_maps=$(awk -F "\t" '/Uniquely mapped reads number/{print $2}' aligned_Log.final.out)
	star_multi_maps=$(awk -F "\t" '/Number of reads mapped to multiple loci/{print $2}' aligned_Log.final.out)
	star_unal_too_many_loci=$(awk -F "\t" '/Number of reads mapped to too many loci/{print $2}' aligned_Log.final.out)
	star_unal_many_mismatches=$(awk -F "\t" '/Number of reads unmapped: too many mismatches/{print $2}' aligned_Log.final.out)
	star_unal_too_short=$(awk -F "\t" '/Number of reads unmapped: too short/{print $2}' aligned_Log.final.out)
	star_unal_other=$(awk -F "\t" '/Number of reads unmapped: other/{print $2}' aligned_Log.final.out)
	
	
	star_finished_successfully=1
	
	if [ ${#star_finished} -lt 1 ]
	then
		echo "WARNING: STAR did not finish succesfully!!" >> sample.log
		star_finished_successfully=0
	fi
	
	
	sql_command="UPDATE !{params.db_table_name}
	  SET	alig_star_finished  = '$star_finished_successfully',
	    alig_star_input			= '$star_input',
		alig_star_unique		= '$star_unique_maps',
		alig_star_multi			= '$star_multi_maps',
		alig_star_unal_too_many_loci  = '$star_unal_too_many_loci',
		alig_star_unal_many_mismatches = '$star_unal_many_mismatches',
		alig_star_unal_too_short       = '$star_unal_too_short',
		alig_star_unal_other           = '$star_unal_other'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "Uploaded STAR summaries to database." >> sample.log
	echo >> sample.log
	
	'''	
}


process dedup{
	cpus 1
	time '5h'
	memory '10GB'
	module 'SAMtools/1.11-GCCcore-10.2.0'
	
	errorStrategy 'retry'
	maxRetries 15
	
	input:
	  tuple path("aligned.bam"), path('trimmed_I1.fq'),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_align
	
	output:
	  tuple path("dedup.sorted.dedup.bam"),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_dedup_for_qc
	  tuple val(sample), path("dedup.sorted.dedup.bam") into ch_dedup_for_merging
	
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	nudup.py --paired-end \
			 -f trimmed_I1.fq \
			 --start 8 \
			 --length 8 \
			 -o dedup \
			 aligned.bam >> nudup.log
	
	echo "-----  NuDup output log   -----" >> sample.log
	cat nudup.log >> sample.log
	echo "     ~~~~~~ NuDup self log " >> sample.log
	cat dedup_dup_log.txt >> sample.log
	
	
	if [[ $(cat dedup_dup_log.txt | cut -f1 | tail -1) -eq 0 ]]
	then
		echo "Empty input, NuDup may have done something strange, retrying."
		exit 1
	fi
	
	
	
	# Check results and summary
	nudup_input=$(tail -1 dedup_dup_log.txt | cut -f1)
	nudup_unaligned=$(tail -1 dedup_dup_log.txt | cut -f2)
	nudup_pos_dup=$(tail -1 dedup_dup_log.txt | cut -f3)
	nudup_umi_dup=$(tail -1 dedup_dup_log.txt | cut -f5)
	
	
	sql_command="UPDATE !{params.db_table_name}
	  SET	alig_nudup_input	= '$nudup_input',
		alig_nudup_pos_dup		= '$nudup_pos_dup',
		alig_nudup_umi_dup		= '$nudup_umi_dup',
		alig_nudup_unaligned    = '$nudup_unaligned'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "Uploaded to database." >> sample.log
	echo >> sample.log
	'''
}



ch_tuple_for_merging = ch_dedup_for_merging
			.groupTuple()

						

process export_bam{
	cpus '5'
	time '30min'
	memory '3GB'
	
	module 'SAMtools/1.13-GCCcore-10.2.0'

	input:
          tuple val(sample), path(bam_files, stageAs: 'dedup*.bam') from ch_tuple_for_merging
	
	output:
	  path "*.bam"
	  path "*.bam.bai"
	
	
	publishDir mode: 'copy',
		pattern: '*.bam',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_bams/',
		saveAs: { sample + '.bam' }
	
	publishDir mode: 'copy',
                pattern: '*.bam.bai',
                path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_bams/',
                saveAs: { sample + '.bam.bai' }
	
	shell:
	'''
	
	samtools merge -@ $SLURM_CPUS_PER_TASK -o !{sample}.bam dedup*.bam
	
	samtools index -@ $SLURM_CPUS_PER_TASK !{sample}.bam
	
	
	'''
}




process fastqc_post{
	cpus '2'
	time '2h'
	memory '10GB'

	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3

	module 'FastQC/0.11.9-Java-11'
	
		
	input:
	  tuple path("dedup.bam"),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_dedup_for_qc
	
	output:
	  tuple path("dedup.bam"),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) into ch_qc_post
	  path "*.html"
	
	
    publishDir mode: 'copy',
			pattern: '*.html',
			path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/'+script_version+'_logs/',
			saveAs: {filename -> sample+sample_suffix+'_fastqc_post.html'}
	
	
	shell:
	'''
	touch !{sample}!{sample_suffix}
	
	fastqc --extract --threads $SLURM_CPUS_PER_TASK dedup.bam
	
	
	echo "---------  FASTQC post  ---------"   >> sample.log
	cat dedup_fastqc/summary.txt               >> sample.log
	
	
	# Process results: read summary file as a bash array
	
	arr=($(cat dedup_fastqc/summary.txt | cut -f1))
	
	
	
	sql_command="UPDATE !{params.db_table_name}
	  SET	alig_fqcpost_basic_stats	  = '${arr[0]}',
		alig_fqcpost_per_base_quality	  = '${arr[1]}',
		alig_fqcpost_per_tile_quality	  = '${arr[2]}',
		alig_fqcpost_per_seq_scores		  = '${arr[3]}',
		alig_fqcpost_per_base_seq_content = '${arr[4]}',
		alig_fqcpost_per_base_gc		  = '${arr[5]}',
		alig_fqcpost_per_base_n			  = '${arr[6]}',
		alig_fqcpost_seq_length			  = '${arr[7]}',
		alig_fqcpost_duplication		  = '${arr[8]}',
		alig_fqcpost_overrepresented_seq  = '${arr[9]}',
		alig_fqcpost_adapter_content	  = '${arr[10]}'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "FastQC-post uploaded to database." >> sample.log
	echo >> sample.log
	'''	
}



process finalize{
	cpus '1'
	time '1h'
	memory '1GB'
	
	errorStrategy { sleep((task.attempt - 1) * 2000); return 'retry' }
	maxRetries 3

	module 'SAMtools/1.13-GCCcore-10.2.0'
	
	input:
	  tuple path("dedup.bam"),
			path("sample.log"), val(sample), val(batch), val(sample_suffix) from ch_qc_post
	
	output:
	  path "sample.log"
	  path "final.log"
	
	publishDir mode: 'copy',
		pattern: 'sample.log',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '.log'}

	publishDir mode: 'copy',
		pattern: 'final.log',
		path: '/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/' + script_version + '_logs/',
		saveAs: { sample + sample_suffix + '.final.log'}
	  

	shell:
	'''
	sample=!{sample}!{sample_suffix}
	
	# Read what has been uploaded to DB,
	# read the bam stats,
	# check that everything consistent, create a final log with errors and warnings
	
	
	# ----    Final counts    ----
	
	bam_stats=$(samtools flagstat dedup.bam)
		
	regex="([[:digit:]]+) \\+ 0 in total \\(QC-passed reads \\+ QC-failed reads\\)[[:space:]]+[[:digit:]]+ \\+ 0 primary[[:space:]]+([[:digit:]]+) \\+ 0 secondary[[:space:]]+([[:digit:]]+) \\+ 0 supplementary[[:space:]]+([[:digit:]]+) \\+ 0 duplicates[[:space:]]+[[:digit:]]+ \\+ 0 primary duplicates[[:space:]]+([[:digit:]]+) \\+ 0 mapped \\([[:digit:].]+% : N/A\\)[[:space:]]+[[:digit:]]+ \\+ 0 primary mapped \\([[:digit:].]+% : N/A\\)[[:space:]]+([[:digit:]]+) \\+ 0 paired in sequencing[[:space:]]+[[:digit:]]+ \\+ 0 read1[[:space:]]+[[:digit:]]+ \\+ 0 read2[[:space:]]+([[:digit:]]+) \\+ 0 properly paired \\([[:digit:].]+% : N/A\\)[[:space:]]+([[:digit:]]+) \\+ 0 with itself and mate mapped[[:space:]]+0 \\+ 0 singletons \\([[:digit:].]+% : N/A\\)[[:space:]]+0 \\+ 0 with mate mapped to a different chr[[:space:]]+0 \\+ 0 with mate mapped to a different chr \\(mapQ>=5\\)"

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
		echo  			    >> sample.log
		echo "Content:" 	>> sample.log
		echo $bam_stats 	>> sample.log
		echo  		    	>> sample.log
		alig_comments="$warnings -- failed to read flagstats"
	fi
	
	
	
	
	echo "Finishing !{script_version} of sample $sample on $(date +%F)." >> sample.log
	echo "Finishing !{script_version} of sample $sample on $(date +%F)." >> final.log
	
	
	# register warnings
	warnings=""
	
	
	sql_command="SELECT * FROM !{params.db_table_name}
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	
	
	# run command twice, once to get col names, once to get data row
	
	IFS=$'\t' colnames=($(mysql --host=23.229.215.131 \
	                  --user=$username \
					  --password=$password \
					  --database=!{params.db_name} \
					  --execute="$sql_command" \
					  --batch \
					  | head -1))
	IFS=$'\t' existing=($(mysql --host=23.229.215.131 -u${username} -p${password} !{params.db_name} -e "${sql_command}" --batch --raw --silent))
	
	if [[ ${#colnames[@]} -ne ${#existing[@]} ]]
	then
		warnings="$warnings -- Error retrieving database: number of columns does not match:"
		warnings="$warnings ~ Expected ${#colnames[@]}, got ${#existing[@]} instead."
	fi
	
	
	# Create associative array (i.e. named array/dict/hash)
	declare -A db_content
	for (( i=0; i < ${#colnames[@]}; i++ ))
	do
		echo "  $i: ${colnames[$i]} is ${existing[$i]}"
		db_content["${colnames[$i]}"]="${existing[$i]}"
	done
	
	
	
	# Check that all values consistent
	
	if [ $(( ${db_content['alig_bbduk_input']} - ${db_content['alig_bbduk_removed']} )) -ne ${db_content['alig_bbduk_result']} ]
	then
		warnings="$warnings -- BBduk number of reads do not add up"
	fi
	
	
	if [ ${db_content['alig_merged_reads']} -ne $(( 2 * ${db_content['alig_bbduk_input']} )) ]
	then
		warnings="$warnings -- number of reads in merge output and in BBduk input do not match"
	fi
	
	
	if [ ${db_content['alig_star_finished']} -lt 1 ]
	then
		warnings="$warnings -- STAR did not finish successfully"
	fi


	if [ ${db_content['alig_bbduk_result']} -ne $(( 2 * ${db_content['alig_star_input']} )) ]
	then
		warnings="$warnings -- number of output reads from BBduk and input reads for STAR do not match"
	fi
	
	
	if [ ${db_content['alig_star_input']} -ne $(( ${db_content['alig_star_unique']} + ${db_content['alig_star_multi']} + ${db_content['alig_star_unal_too_many_loci']} + ${db_content['alig_star_unal_many_mismatches']} +${db_content['alig_star_unal_too_short']} + ${db_content['alig_star_unal_other']} )) ]
	then
		warnings="$warnings -- number of STAR reads do not match"
	fi
	
	
	if [ ${db_content['alig_nudup_unaligned']} -ne $(( 2 * (${db_content['alig_star_unal_many_mismatches']} + ${db_content['alig_star_unal_too_many_loci']} + ${db_content['alig_star_unal_too_short']} + ${db_content['alig_star_unal_other']}) )) ]
	then
		warnings="$warnings -- number of unaligned reads do not match between STAR and NuDup"
	fi
	
	
	
	
	echo "##################################################"  >> sample.log
	echo                                                       >> sample.log
	echo "------------   FINAL RESULTS   -------------"        >> sample.log
	echo                                                       >> sample.log
	echo "$warnings"                                           >> sample.log
	echo                                                       >> sample.log
	echo  "-------------------------------------------"        >> sample.log
	echo                                                       >> sample.log
	echo "Quality metrics"                                     >> sample.log
	echo                                                       >> sample.log
	echo "Merged reads: ${db_content['alig_merged_reads']}"    >> sample.log
	echo "  --"                                                  >> sample.log
	echo "BBduk input: ${db_content['alig_bbduk_input']}"        >> sample.log
	echo "BBduk removed: ${db_content['alig_bbduk_removed']}"    >> sample.log
	echo "BBduk result: ${db_content['alig_bbduk_result']}"      >> sample.log
	echo "  --"                                                                                   >> sample.log
	echo "STAR input: ${db_content['alig_star_input']}"                                           >> sample.log
	echo "STAR mapped:\t"$((${db_content['alig_star_unique']}+${db_content['alig_star_multi']}))  >>sample.log
	echo "STAR unique maps: ${db_content['alig_star_unique']}"                                    >> sample.log
	echo "STAR multi maps: ${db_content['alig_star_multi']}"                                      >> sample.log
	echo "  --"                                                                     >> sample.log
	echo "NuDup input: ${db_content['alig_nudup_input']}"                           >> sample.log
	echo "NuDup positional duplicates: ${db_content['alig_nudup_pos_dup']}"         >> sample.log
	echo "NuDup UMI duplicates: ${db_content['alig_nudup_umi_dup']}"                >> sample.log
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
	
	
	# Create a warnings file, of size 0 when there was no problem
	echo $warnings > final.log
	
	if [[ -z $warnings ]]
	then
	    alig_quality_check=0
	else
	    alig_quality_check=1
	fi





	echo "Creating database entry" >> sample.log

		
	
	sql_command="UPDATE !{params.db_table_name}
	  SET	alig_completed			= '1',
		alig_quality_check			= '$alig_quality_check',
		alig_comments				= '$warnings',
		alig_final_total			= '"$final_total"',
		alig_final_secondary		= '"$final_secondary"',
		alig_final_supplementary	= '"$final_supplementary"',
		alig_final_duplicates		= '"$final_duplicates"',
		alig_final_mapped			= '"$final_mapped"',
		alig_final_paired			= '"$final_paired"',
		alig_final_properly_paired	= '"$final_properly_paired"'
	  WHERE alig_lib_id		= '!{sample}!{sample_suffix}' AND
		alig_pid			= '!{script_version}' AND
		alig_sequ_batch		= '!{batch}';"
	
	source ~/.cengen_database.id
	mysql --host=23.229.215.131 \
		--user=$username \
		--password=$password \
		--database=!{params.db_name} \
		--execute="$sql_command"

	echo "Uploaded to database." >> sample.log
	echo >> sample.log
	
	echo "Finished."
	'''
}





