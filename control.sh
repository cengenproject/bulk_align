(base) [aw853@c14n11 bulk_align]$ samtools idxstats /SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn6_bams/IL2r95ter.bam
I       15072434        974     0
II      15279421        2       0
III     13783801        0       0
IV      17493829        0       0
V       20924180        0       0
X       17718942        0       0
MtDNA   13794   0       0
*       0       0       3366


(base) [aw853@c14n11 bulk_align]$ samtools idxstats /SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn6_bams/PVMr122ter.bam
I       15072434        224     0
II      15279421        2       0
III     13783801        2       0
IV      17493829        0       0
V       20924180        0       0
X       17718942        2       0
MtDNA   13794   2       0
*       0       0       3808



(base) [aw853@c15n01 bulk_align]$ cat /SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn6_junctions/IL2r95nondec.SJ.tab
I       14785644        14785698        2       2       1       1       0       17
I       14893743        14894575        2       2       1       0       1       12
I       14894644        14894692        2       2       1       0       1       20

(base) [aw853@c15n01 bulk_align]$ cat /SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn6_junctions/PVMr122nondec.SJ.tab
I       15001959        15002533        1       1       1       1       0       16











Quality metrics

Merged reads: 8040
  --
BBduk input: 4020
BBduk removed: 0
BBduk result: 4020
  --
STAR input: 2010
STAR mapped:    327
STAR unique maps: 245
STAR multi maps: 82
  --
NuDup input: 488
NuDup positional duplicates: 138
NuDup UMI duplicates: 0
  --
Final file total: 4342
Final file secondary: 322
Final file supplementary: 0
Final file duplicates: 0
Final file mapped: 976
Final file paired: 4020
Final file properly paired: 654

##################################################

    ------------    Merge    ------------
Merging files Sample_IL2r95

Merged reads: 8040


    ------------    BBduk    ------------
java -ea -Xmx6246m -Xms6246m -cp /gpfs/ycga/apps/hpc/software/BBMap/38.90-GCCcore-10.2.0/current/ jgi.BBDuk in1=merged_R1.fastq.gz in2=merged_R2.fastq.gz out1=trimmed_R1.fastq.gz out2=trimmed_R2.fastq.gz ref=/ycga-gpfs/apps/hpc/software/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
Executing jgi.BBDuk [in1=merged_R1.fastq.gz, in2=merged_R2.fastq.gz, out1=trimmed_R1.fastq.gz, out2=trimmed_R2.fastq.gz, ref=/ycga-gpfs/apps/hpc/software/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa, ktrim=r, k=23, mink=11, hdist=1, tpe, tbo]
Version 38.90

maskMiddle was disabled because useShortKmers=true
0.018 seconds.
Initial:
Memory: max=6549m, total=6549m, free=6531m, used=18m

Added 3987 kmers; time:         0.024 seconds.
Memory: max=6549m, total=6549m, free=6523m, used=26m

Input is being processed as paired
Started output streams: 0.012 seconds.
Processing time:                0.196 seconds.

Input:                          4020 reads              406020 bases.
KTrimmed:                       84 reads (2.09%)        2154 bases (0.53%)
Trimmed by overlap:             12 reads (0.30%)        298 bases (0.07%)
Total Removed:                  0 reads (0.00%)         2452 bases (0.60%)
Result:                         4020 reads (100.00%)    403568 bases (99.40%)

Time:                           0.235 seconds.
Reads Processed:        4020    17.13k reads/sec
Bases Processed:        406k    1.73m bases/sec
    ------------    STAR    ------------
Using STAR version 2.7.7a
Sep 23 10:27:17 ..... started STAR run
Sep 23 10:27:17 ..... loading genome
Sep 23 10:27:18 ..... started mapping
Sep 23 10:27:20 ..... finished mapping
Sep 23 10:27:20 ..... started sorting BAM
Sep 23 10:27:20 ..... finished successfully
Creating database entry
Uploaded to database.
