# subsample FASTQ

# subsample R1
zcat PVMr122_AGAACCAG_L001_R1_001.fastq.gz | paste - - - - | shuf -n 1000 | tr '\t' '\n' > PVMr122_subsample_L001_R1_001.fastq

# Unzip the files to look in
zcat PVMr122_AGAACCAG_L001_R2_001.fastq.gz > PVMr122_R2.fq
zcat PVMr122_AGAACCAG_L001_R3_001.fastq.gz > PVMr122_R3.fq

# Find the read names and extract the corresponding reads
cat PVMr122_subsample_L001_R1_001.fastq | paste - - - - | cut -f1 -d' ' > to_search.txt
grep -A3 --no-group-separator -f to_search.txt PVMr122_R2.fq > PVMr122_subsample_L001_R2_001.fastq
grep -A3 --no-group-separator -f to_search.txt PVMr122_R3.fq > PVMr122_subsample_L001_R3_001.fastq


# Same with a second set
# subsample R1
zcat PVMr122_AGAACCAG_L001_R1_002.fastq.gz | paste - - - - | shuf -n 1010 | tr '\t' '\n' > PVMr122_subsample_L001_R1_002.fastq

# Unzip the files to look in
zcat PVMr122_AGAACCAG_L001_R2_002.fastq.gz > PVMr122_R2.fq
zcat PVMr122_AGAACCAG_L001_R3_002.fastq.gz > PVMr122_R3.fq

# Find the read names and extract the corresponding reads
cat PVMr122_subsample_L001_R1_002.fastq | paste - - - - | cut -f1 -d' ' > to_search.txt
grep -A3 --no-group-separator -f to_search.txt PVMr122_R2.fq > PVMr122_subsample_L001_R2_002.fastq
grep -A3 --no-group-separator -f to_search.txt PVMr122_R3.fq > PVMr122_subsample_L001_R3_002.fastq

wc -l PVMr122_subsample_L001_R*
