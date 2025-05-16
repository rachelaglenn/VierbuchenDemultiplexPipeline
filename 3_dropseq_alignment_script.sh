#!/bin/bash 
inputpath="/scratch-rz/users/dol24jn/vierbuchen/samples/C_P2-5-11-25-30/rawdata/"
fastq1="$1"
fastq2="$2"
samplename="$3"
outputpath="/scratch-rz/users/dol24jn/vierbuchen/results/${samplename}/"

picardlocation="/home/dol24jn/software/Drop-seq_tools-2.5.4/jar/lib/picard-2.26.10.jar"
dropseqlocation="/home/dol24jn/software/Drop-seq_tools-2.5.4/jar/dropseq.jar"
starlocation="~/software/STAR/STAR/bin/Linux_x86_64_static/STAR"
genomedir="/home/dol24jn/software/database/dropseq_GRCm38/STAR"
referencesequence="/home/dol24jn/software/database/dropseq_GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa"
geneintervals="/home/dol24jn/software/database/dropseq_GRCm38/Mus_musculus.GRCm38..rRNA.intervals"
annotationsfile="/home/dol24jn/software/database/dropseq_GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.refFlat"


mkdir "$outputpath"


#1. Convert paired reads of 10x in 1 single unaligned BAM file (R1 is barcode + umi R2 is bio read)

java -jar ${picardlocation} FastqToSam \
FASTQ= "${inputpath}${fastq1}" \
FASTQ2= "${inputpath}${fastq2}" \
OUTPUT= "${outputpath}${samplename}.bam" \
SAMPLE_NAME="${samplename}"

#check first lines of file
samtools view "${outputpath}${samplename}.bam" | head -10 > "${outputpath}${samplename}_validation.bam"

#2. create bam with barcode tagged reads

java -jar ${dropseqlocation} TagBamWithReadSequenceExtended \
--BARCODED_READ 1 \
--INPUT "${outputpath}${samplename}.bam" \
--OUTPUT "${outputpath}${samplename}_tagged_Cell.bam" \
--SUMMARY "${outputpath}${samplename}_tagged_Cell_summary.txt" \
--BASE_RANGE 1-16 \
--BASE_QUALITY 10 \
--DISCARD_READ False \
--TAG_NAME XC \
--NUM_BASES_BELOW_QUALITY 1 

#check first lines of file
samtools view "${outputpath}${samplename}_tagged_Cell.bam" | head -10 > "${outputpath}${samplename}_tagged_Cell_validation.bam"


java -jar ${dropseqlocation} TagBamWithReadSequenceExtended \
--BARCODED_READ 1 \
--INPUT "${outputpath}${samplename}_tagged_Cell.bam" \
--OUTPUT "${outputpath}${samplename}_tagged_CellMolecular.bam" \
--SUMMARY "${outputpath}${samplename}_tagged_CellMolecular_summary.txt" \
--BASE_RANGE 17-28 \
--BASE_QUALITY 10 \
--DISCARD_READ TRUE \
--TAG_NAME XM \
--NUM_BASES_BELOW_QUALITY 1 

#check first lines of file
samtools view "${outputpath}${samplename}_tagged_CellMolecular.bam" | head -10 > "${outputpath}${samplename}_tagged_CellMolecular_validation.bam"

#3. Filter low quality barcodes
java -jar ${dropseqlocation} FilterBam \
TAG_REJECT=XQ \
INPUT="${outputpath}${samplename}_tagged_CellMolecular.bam" \
OUTPUT="${outputpath}${samplename}_tagged_CellMolecular_filtered.bam" 

#check first lines of file
samtools view "${outputpath}${samplename}_tagged_CellMolecular_filtered.bam" | head -10 > "${outputpath}${samplename}_tagged_CellMolecular_filtered_validation.bam"


#Trimming
java -jar ${dropseqlocation} TrimStartingSequence \
INPUT="${outputpath}${samplename}_tagged_CellMolecular_filtered.bam" \
OUTPUT="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed.bam" \
OUTPUT_SUMMARY="${outputpath}${samplename}_tagged_CellMolecular_filtered_report.txt" \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTACATGGG \
MISMATCHES=0 \
NUM_BASES=5 


#check first lines of file
samtools view "${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed.bam" | head -10 > "${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed_validation.bam"

#Trimming2
java -jar ${dropseqlocation} PolyATrimmer \
INPUT="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed.bam" \
OUTPUT="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2.bam" \
OUTPUT_SUMMARY="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2_report.txt" \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=false


#check first lines of file
samtools view "${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2.bam" | head -10 > "${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2_validation.bam"

#convert back to fastq
java -jar ${picardlocation} SamToFastq \
INPUT="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2.bam" \
FASTQ="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2.fastq"

${starlocation} \
--genomeDir ${genomedir} \
--readFilesIn "${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2.fastq" \
--outFileNamePrefix "${outputpath}${samplename}_" \
--runThreadN 64

java -Xmx4g -jar ${picardlocation} SortSam \
I="${outputpath}${samplename}_Aligned.out.sam" \
O="${outputpath}${samplename}_Aligned.out.sorted.bam" \
SO=queryname


java -Xmx4g -jar ${picardlocation} MergeBamAlignment \
REFERENCE_SEQUENCE=${referencesequence} \
UNMAPPED_BAM="${outputpath}${samplename}_tagged_CellMolecular_filtered_trimmed2.bam" \
ALIGNED_BAM="${outputpath}${samplename}_Aligned.out.sorted.bam" \
OUTPUT="${outputpath}${samplename}_merged.bam" \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false


#check first lines of file
samtools view "${outputpath}${samplename}_merged.bam" | head -10 > "${outputpath}${samplename}_merged_validation.bam"


java -jar ${dropseqlocation} TagReadWithInterval \
I="${outputpath}${samplename}_merged.bam" \
O="${outputpath}${samplename}_merged_gene_tagged.bam" \
INTERVALS= ${geneintervals} \
TAG=XG


java -jar ${dropseqlocation} TagReadWithGeneFunction \
INPUT="${outputpath}${samplename}_merged_gene_tagged.bam" \
O="${outputpath}${samplename}_merged_gene_function_tagged.bam" \
ANNOTATIONS_FILE=${annotationsfile}


#check first lines of file
samtools view "${outputpath}${samplename}_merged_gene_function_tagged.bam" | head -10 > "${outputpath}${samplename}_merged_gene_function_tagged_validation.bam"

# Stage 5: bead repair

java -jar ${dropseqlocation} DetectBeadSubstitutionErrors \
INPUT="${outputpath}${samplename}_merged_gene_function_tagged.bam" \
OUTPUT="${outputpath}${samplename}_merged_repaired.bam" \
MIN_UMIS_PER_CELL=20 \
OUTPUT_REPORT="${outputpath}${samplename}_merged_repaired_report.txt"


java -jar ${dropseqlocation} DetectBeadSynthesisErrors \
INPUT="${outputpath}${samplename}_merged_repaired.bam" \
MIN_UMIS_PER_CELL=20 \
OUTPUT_STATS="${outputpath}${samplename}_merged_repaired2_error_stats.txt" \
SUMMARY="${outputpath}${samplename}_merged_repaired2_summary.txt" \
REPORT="${outputpath}${samplename}_merged_repaired2_report.txt" \
CREATE_INDEX=true \
OUTPUT="${outputpath}${samplename}_merged_repaired2.bam"

#Quantification

#Digital Gene Expression
java -jar ${dropseqlocation} DigitalExpression \
I="${outputpath}${samplename}_merged_repaired2.bam" \
O="${outputpath}${samplename}_counts.dge.txt.gz" \
SUMMARY="${outputpath}${samplename}_counts_summary.txt" \
NUM_CORE_BARCODES=30000


