# Somatic Mutation Pipeline Script
#
##!/bin/bash

# Setting directories before running the script
PROJECT_DIR="/home1/vineet/pupilBio"
TRIMMOMATIC_JAR="/home1/vineet/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/home1/vineet/tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
REF_DIR="/home1/vineet/pupilBio/reffile"
OUTPUT_DIR="$PROJECT_DIR/output"
GATK_PATH="/home1/vineet/tools/gatk-4.6.1.0/gatk"
MUTECT2_SUPPORTING_FILES="$PROJECT_DIR/supportingfile_mutect2"
DATA_SOURCES_PATH="/home1/vineet/pupilBio/funco/funcotator_dataSources.v1.8.hg38.20230908s"
THREADS=8

# Creating directory just in case if missed in between
mkdir -p $PROJECT_DIR/qc
mkdir -p $PROJECT_DIR/trimmo
mkdir -p $PROJECT_DIR/qc2
mkdir -p $PROJECT_DIR/bwa
mkdir -p $PROJECT_DIR/gatk-sort
mkdir -p $PROJECT_DIR/recalibration
mkdir -p $PROJECT_DIR/gatk-bqsr
mkdir -p $PROJECT_DIR/gatk-metrics
mkdir -p $PROJECT_DIR/mutect2
mkdir -p $PROJECT_DIR/results
mkdir -p $PROJECT_DIR/tmp

# Step 1: Quality of the raw reads were checked using tool- (FastQC)
fastqc $PROJECT_DIR/*.fastq.gz -o $PROJECT_DIR/qc
multiqc $PROJECT_DIR/qc -o $PROJECT_DIR/qc/multiqc_pretrim_report

# Step 2: Trimming of adapter and low quality sequences  with  Trimmomatic
#Trimming of Tumor sample
java -jar "$TRIMMOMATIC_JAR" PE -phred33 -threads "$THREADS" $PROJECT_DIR/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz $PROJECT_DIR/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz $PROJECT_DIR/trimmo/PA220KH-lib09-P19-Tumor_S2_L001_R1_001_paired.fastq.gz $PROJECT_DIR/trimmo/PA220KH-lib09-P19-Tumor_S2_L001_R1_001_unpaired.fastq.gz $PROJECT_DIR/trimmo/PA220KH-lib09-P19-Tumor_S2_L001_R2_001_paired.fastq.gz $PROJECT_DIR/trimmo/PA220KH-lib09-P19-Tumor_S2_L001_R2_001_unpaired.fastq.gz ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

# Trimming of Normal sample
java -jar "$TRIMMOMATIC_JAR" PE -phred33 -threads "$THREADS" $PROJECT_DIR/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz $PROJECT_DIR/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz $PROJECT_DIR/trimmo/PA221MH-lib09-P19-Norm_S1_L001_R1_001_paired.fastq.gz $PROJECT_DIR/trimmo/PA221MH-lib09-P19-Norm_S1_L001_R1_001_unpaired.fastq.gz $PROJECT_DIR/trimmo/PA221MH-lib09-P19-Norm_S1_L001_R2_001_paired.fastq.gz $PROJECT_DIR/trimmo/PA221MH-lib09-P19-Norm_S1_L001_R2_001_unpaired.fastq.gz ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:1 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

# Rechecking quality post trimming
fastqc $PROJECT_DIR/trimmo/*_paired.fastq.gz -o $PROJECT_DIR/qc2
multiqc $PROJECT_DIR/qc2 -o $PROJECT_DIR/qc2/multiqc_posttrim_report

# Step 3: Alignment of trimmed sample using BWA

# Downloading reference files i.e hg38 and index them
wget -P $REF_DIR https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip $REF_DIR/hg38.fa.gz
samtools faidx $REF_DIR/hg38.fa
gatk CreateSequenceDictionary -R $REF_DIR/hg38.fa -O $REF_DIR/hg38.dict
bwa index $REF_DIR/hg38.fa

# Alignment of Normal and Tumor sample using BWA
bwa mem -t $THREADS -R "@RG\tID:PA220KH-Tumor\tPL:ILLUMINA\tSM:PA220KH-Tumor" $REF_DIR/hg38.fa $PROJECT_DIR/trimmo/PA220KH-lib09-P19-Tumor_S2_L001_R1_001_paired.fastq.gz $PROJECT_DIR/trimmo/PA220KH-lib09-P19-Tumor_S2_L001_R2_001_paired.fastq.gz > $PROJECT_DIR/bwa/PA220KH-Tumor.paired.sam

bwa mem -t $THREADS -R "@RG\tID:PA221MH-Normal\tPL:ILLUMINA\tSM:PA221MH-Normal" $REF_DIR/hg38.fa $PROJECT_DIR/trimmo/PA221MH-lib09-P19-Norm_S1_L001_R1_001_paired.fastq.gz $PROJECT_DIR/trimmo/PA221MH-lib09-P19-Norm_S1_L001_R2_001_paired.fastq.gz > $PROJECT_DIR/bwa/PA221MH-Normal.paired.sam

# Step 4: Marking the  duplicates and sorting them using GATK

#  Marking duplicates helps identify PCR artifacts, ensuring downstream analyses are not biased by duplicate reads.
# Sorting ensures BAM files are ready for variant calling steps.
gatk MarkDuplicatesSpark -I $PROJECT_DIR/bwa/PA220KH-Tumor.paired.sam -O $PROJECT_DIR/gatk-sort/PA220KH-Tumor_sorted_dedup_reads.bam --spark-master local[$THREADS]

gatk MarkDuplicatesSpark -I $PROJECT_DIR/bwa/PA221MH-Normal.paired.sam -O $PROJECT_DIR/gatk-sort/PA221MH-Normal_sorted_dedup_reads.bam --spark-master local[$THREADS]

# Step 5: Base quality recalibration using GATK

#this step adjusts base quality scores to correct for systematic errors in sequencing data, improving the accuracy of variant calling.
gatk BaseRecalibrator -I $PROJECT_DIR/gatk-sort/PA220KH-Tumor_sorted_dedup_reads.bam -R $REF_DIR/hg38.fa --known-sites $REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf -O $PROJECT_DIR/recalibration/PA220KH-Tumor_recal_data.table

gatk BaseRecalibrator -I $PROJECT_DIR/gatk-sort/PA221MH-Normal_sorted_dedup_reads.bam -R $REF_DIR/hg38.fa --known-sites $REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf -O $PROJECT_DIR/recalibration/PA221MH-Normal_recal_data.table

#Correcting Base Quality Scores for Improved Accuracy
gatk ApplyBQSR -I $PROJECT_DIR/gatk-sort/PA220KH-Tumor_sorted_dedup_reads.bam -R $REF_DIR/hg38.fa --bqsr-recal-file $PROJECT_DIR/recalibration/PA220KH-Tumor_recal_data.table -O $PROJECT_DIR/gatk-bqsr/PA220KH-Tumor_sorted_dedup_bqsr_reads.bam

gatk ApplyBQSR -I $PROJECT_DIR/gatk-sort/PA221MH-Normal_sorted_dedup_reads.bam -R $REF_DIR/hg38.fa --bqsr-recal-file $PROJECT_DIR/recalibration/PA221MH-Normal_recal_data.table -O $PROJECT_DIR/gatk-bqsr/PA221MH-Normal_sorted_dedup_bqsr_reads.bam

# Step 6: Generating Alignment Quality Reports and Insert Size Analysis
gatk CollectAlignmentSummaryMetrics -R $REF_DIR/hg38.fa -I $PROJECT_DIR/gatk-bqsr/PA220KH-Tumor_sorted_dedup_bqsr_reads.bam -O $PROJECT_DIR/gatk-metrics/PA220KH-Tumor_alignment_metrics.txt

gatk CollectAlignmentSummaryMetrics -R $REF_DIR/hg38.fa -I $PROJECT_DIR/gatk-bqsr/PA221MH-Normal_sorted_dedup_bqsr_reads.bam -O $PROJECT_DIR/gatk-metrics/PA221MH-Normal_alignment_metrics.txt

# Step 7: Identifying Somatic Variants with Mutect2
gatk Mutect2 -R $REF_DIR/hg38.fa -I $PROJECT_DIR/gatk-bqsr/PA220KH-Tumor_sorted_dedup_bqsr_reads.bam -I $PROJECT_DIR/gatk-bqsr/PA221MH-Normal_sorted_dedup_bqsr_reads.bam -tumor PA220KH-Tumor -normal PA221MH-Normal --germline-resource $MUTECT2_SUPPORTING_FILES/af-only-gnomad.hg38.vcf.gz --panel-of-normals $MUTECT2_SUPPORTING_FILES/1000g_pon.hg38.vcf.gz -O $PROJECT_DIR/mutect2/PA220KH_somatic_variants_mutect2.vcf.gz --f1r2-tar-gz $PROJECT_DIR/mutect2/PA220KH_f1r2.tar.gz


# STEP 8:Checking for Cross-Sample Contamination

# Run GetPileupSummaries for Tumor sample
${GATK_PATH} GetPileupSummaries --java-options '-Xmx50G' --tmp-dir ${PROJECT_DIR}/tmp/ -I $PROJECT_DIR/gatk-bqsr/PA220KH-Tumor_sorted_dedup_bqsr_reads.bam -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list -O ${PROJECT_DIR}/results/PA220KH-Tumor_getpileupsummaries.table

# Run GetPileupSummaries for Normal sample
${GATK_PATH} GetPileupSummaries --java-options '-Xmx50G' --tmp-dir ${PROJECT_DIR}/tmp/ -I $PROJECT_DIR/gatk-bqsr/PA221MH-Normal_sorted_dedup_bqsr_reads.bam -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list -O ${PROJECT_DIR}/results/PA221MH-Normal_getpileupsummaries.table


# Calculate contamination
${GATK_PATH} CalculateContamination -I $PROJECT_DIR/results/PA220KH-T_getpileupsummaries.table -matched $PROJECT_DIR/results/PA221MH-N_getpileupsummaries.table -O $PROJECT_DIR/results/PA220KH_pair_calculatecontamination.table

echo "CalculateContamination completed successfully."

# Calculate background mutation rate for normal sample
echo "Calculating background mutation rate for the normal sample..."
awk '{if ($1 != "contig") sum_AF+=$5; count++} END {print "Average Background Mutation Rate:", sum_AF/count}' $PROJECT_DIR/results/PA221MH-N_getpileupsummaries.table > $PROJECT_DIR/results/PA221MH-N_background_mutation_rate.txt


# STEP 9:Analyzing Read Orientation Bias for Accurate Variant Calling

# Run LearnReadOrientationModel
${GATK_PATH} LearnReadOrientationModel -I $PROJECT_DIR/results/PA220KH_f1r2.tar.gz -O $PROJECT_DIR/results/read-orientation-model.tar.gz

echo "LearnReadOrientationModel completed successfully."

# STEP 10: Refining Variants Identified by Mutect2

# Run FilterMutectCalls
${GATK_PATH} FilterMutectCalls -V $PROJECT_DIR/results/PA220KH_somatic_variants_mutect2.vcf.gz -R ${REF_DIR}/hg38.fa --contamination-table $PROJECT_DIR/results/PA220KH_pair_calculatecontamination.table --ob-priors $PROJECT_DIR/results/read-orientation-model.tar.gz -O $PROJECT_DIR/results/PA220KH_somatic_variants_filtered_mutect2.vcf

echo "FilterMutectCalls completed successfully."

# STEP 11: Adding Functional Annotations to Variants

# Run the GATK Funcotator tool
${GATK_PATH} Funcotator --variant $PROJECT_DIR/results/PA220KH_somatic_variants_filtered_mutect2.vcf --reference $REF_DIR/hg38.fa --ref-version hg38 -data-sources-path ${DATA_SOURCES_PATH} --output $PROJECT_DIR/results/PA220KH_somatic_variants_functotated.vcf --output-file-format VCF


