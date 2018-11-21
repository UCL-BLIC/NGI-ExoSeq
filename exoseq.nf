#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

========================================================================================
               N G I - E X O S E Q    B E S T    P R A C T I C E
========================================================================================
 Exome sequencing pipeline based on NGI_exoseq with FastQC, Trim_Galore and MultiQC
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - preprocess (FASTQC, TrimGalore)
 - alignment (BWA)
 - recalibration (GATK)
 - realignment (GATK)
 - variant calling (GATK)
 - vairant evaluation (SnpEff)
 - report (MultiQC) 
*/

// Package version
version = '1.0'

// Help message
helpMessage = """
===============================================================================
NGI-ExoSeq : Exome/Targeted sequence capture best practice analysis v${version}
	     (with FastQC, TrimGalore and MultiQC steps)
===============================================================================

Usage: nextflow_exoseq  --reads '*_R{1,2}.fastq.gz' --genome GRCh37

*************************************************************************************************************************
NOTE: fastq files are expected to have the typical Illumina naming convention (Ex: SampleName_S1_L001_R1_001.fastq.gz)
	to make sure that lanes are merged correctly and read groups are written correctly:
   		@RG	ID:SampleName_S1_L001	SM:SampleName	PL:illumina	PU:SampleName_S1_L001
		@RG	ID:SampleName_S1_L002	SM:SampleName	PL:illumina	PU:SampleName_S1_L002
*************************************************************************************************************************

This is a typical usage where the required parameters (with no defaults) were
given. The all avialable paramaters are listed below based on category

Required parameters:
--reads                        Absolute path to project directory
--genome                       Name of iGenomes reference

Output:
--outdir                       Path where the results to be saved [Default: './results']

Options:
--singleEnd                    Specifies that the input is single end reads

Kit files:
--kit                          Kit used to prep samples [Default: 'agilent_v5']
--bait                         Absolute path to bait file
--target                       Absolute path to target file
--target_bed                   Absolute path to target bed file (snpEff compatible format)

Genome/Variation files:
--dbsnp                        Absolute path to dbsnp file
--hapmap                       Absolute path to hapmap file
--omni                         Absolute path to omni file
--gfasta                       Absolute path to genome fasta file
--bwa_index                    Absolute path to bwa genome index

Other options:
--project                      Uppnex project to user for SLURM executor
--email                        Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
-name                          Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

For more detailed information regarding the paramaters and usage refer to NGI_exoseq package
documentation at https:// github.com/SciLifeLab/NGI-ExoSeq
"""

// Variables and defaults

params.singleEnd = false
params.saveTrimmed = true
params.notrim = false
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

params.help = false
params.reads = false
//params.genome = 'GRCh37'   #- defined in config
params.clusterOptions = false
params.project = false
//params.kit = 'agilent_v5'   #- defined in config
params.bait = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].bait ?: false : false
params.target = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target ?: false : false
params.target_bed = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target_bed ?: false : false
params.dbsnp = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].dbsnp ?: false : false
params.hapmap = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].hapmap ?: false : false
params.omni = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].omni ?: false : false
params.gfasta = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta ?: false : false
params.bwa_index = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].bwa_index ?: false : false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Nextflow version check
nf_required_version = '0.25.0'
try {
    if( ! workflow.nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
        }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// Check blocks for ceratin required parameters, to see they are given and exists
if (!params.reads || !params.genome){
    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
}
if (!params.kitFiles[ params.kit ] && ['bait', 'target'].count{ params[it] } != 2){
    exit 1, "Kit '${params.kit}' is not available in pre-defined config, so " +
            "provide all kit specific files with option '--bait' and '--target'"
}
if (!params.metaFiles[ params.genome ] && ['gfasta', 'bwa_index', 'dbsnp', 'hapmap', 'omni'].count{ params[it] } != 5){
    exit 1, "Genome '${params.genome}' is not available in pre-defined config, so you need to provide all genome specific " +
            "files with options '--gfasta', '--bwa_index', '--dbsnp', '--hapmap', '--omni' and '--target'"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input read files
 * Collect fastq files on sample_lane basis from given project directory
 */
    
 Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimgalore }


/* 
 * ADDED FASTQC STEP
*/

process fastqc {
    tag "$sample"

    publishDir "${params.outdir}/${sample_name}/preprocess/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(sample), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    sample_name = sample - ~/(_S\d*)?(_L\d*)?$/
    """
    fastqc -q $reads
    """
}

/*
 * ADDED TRIM_GALORE
 */

if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
    trimgalore_fastqc_reports = []
} else {
    process trim_galore {
        tag "$sample"

        publishDir "${params.outdir}/${sample_name}/preprocess/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(sample), file(reads) from raw_reads_trimgalore

        output:
        set val(sample), file('*.fq.gz') into trimmed_reads_bwa, trimmed_reads_sam
        file '*trimming_report.txt' into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        sample_name = sample - ~/(_S\d*)?(_L\d*)?$/
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $reads
            """
        }
    }
}


// Align the reads individually for each lane using BWA
process bwaAlign {
    tag "$sample"

    input:
    set val(sample), file (fastq) from trimmed_reads_bwa

    output:
    set val(sample), file("${sample}_bwa.sam") into raw_aln_sam

    script:
    """
    bwa mem \\
        -t 8 \\
        -k 2 \\
        $params.bwa_index \\
        $fastq \\
            > ${sample}_bwa.sam
    """
}

// Create unmapped bam files from raw fastq
process fastqToSam {
    tag "$sample"

    input:
    set val(sample), file(fastq) from trimmed_reads_sam

    output:
    set val(sample), file("${sample}_unaligned.bam") into raw_unaln_bam

    script:
    sam_name = sample - ~/(_S\d*)?(_L\d*)?$/
    fc_lane = sample - ~/^(.*_)$/
    """
    picard FastqToSam \\
        FASTQ=${fastq[0]} \\
        FASTQ2=${fastq[1]} \\
        QUALITY_FORMAT=Standard \\
        OUTPUT=${sample}_unaligned.bam \\
        READ_GROUP_NAME=$fc_lane \\
        SAMPLE_NAME=$sam_name \\
        PLATFORM_UNIT=$fc_lane \\
        VALIDATION_STRINGENCY=SILENT \\
        PLATFORM=illumina \\
        TMP_DIR=tmp \\
        SORT_ORDER=queryname \\
        MIN_Q=0 \\
        MAX_Q=93 \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        VERBOSITY=INFO \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Collect alinged and unaligned files as tuple for each sample
raw_aln_sam
    .cross(raw_unaln_bam)
    .map{ it -> [it[0][0], it[0][1], it[1][1]] }
    .set{ all_sample_bam }

// Create mergerd bam for each sample on flowcell_lane basis
process mergeLaneBam {
    tag "$sample"

    input:
    set val(sample), file(aln_sam), file(unaln_bam) from all_sample_bam

    output:
    set val(sample), file("${sample}_merged.bam") into lanes_merged_bam

    script:
    """
    picard MergeBamAlignment \\
        UNMAPPED_BAM=$unaln_bam \\
        ALIGNED_BAM=$aln_sam \\
        OUTPUT=${sample}_merged.bam \\
        REFERENCE_SEQUENCE=$params.gfasta \\
        PAIRED_RUN=true \\
        TMP_DIR=tmp \\
        CLIP_ADAPTERS=true \\
        IS_BISULFITE_SEQUENCE=false \\
        ALIGNED_READS_ONLY=false \\
        SORT_ORDER=coordinate \\
        READ1_TRIM=0 \\
        READ2_TRIM=0 \\
        MAX_INSERTIONS_OR_DELETIONS=1 \\
        CLIP_OVERLAPPING_READS=true \\
        VERBOSITY=INFO \\
        QUIET=false \\
        VALIDATION_STRINGENCY=SILENT \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Sort the bam merged bam files for lane
process sortLanesBam {
    tag "$sample"

    input:
    set val(sample), file(lane_bam) from lanes_merged_bam

    output:
    set val({sample - ~/(_S\d+)?(_L\d*)?$/}), file("${sample}_sorted.bam") into lanes_sorted_bam

    script:
    """
    picard SortSam \\
        INPUT=$lane_bam \\
        OUTPUT=${sample}_sorted.bam \\
        VERBOSITY=INFO \\
        SORT_ORDER=coordinate \\
        TMP_DIR=tmp \\
        VALIDATION_STRINGENCY=SILENT \\
        MAX_RECORDS_IN_RAM=500000 \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Group the bam file for each sample
lanes_sorted_bam
    .groupTuple()
    .set{ lanes_sorted_bam_group }

// Merge all bam files from lanes to one bam per sample
process mergeSampleBam {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/alignment", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/sorted/, "raw_sorted")}

    input:
    set val(sample), file(sample_bam) from lanes_sorted_bam_group

    output:
    set val(sample), file("${sample}_sorted.bam") into samples_sorted_bam

    script:
    if (sample_bam.properties.target)
        """
        picard MergeSamFiles \\
            ${sample_bam.target.flatten{"INPUT=$it"}.join(' ')} \\
            OUTPUT=${sample}_sorted.bam \\
            SORT_ORDER=coordinate \\
            TMP_DIR=tmp \\
            VALIDATION_STRINGENCY=SILENT \\
            VERBOSITY=INFO \\
            ASSUME_SORTED=false \\
            MERGE_SEQUENCE_DICTIONARIES=false \\
            USE_THREADING=false \\
            QUIET=false \\
            COMPRESSION_LEVEL=5 \\
            MAX_RECORDS_IN_RAM=500000 \\
            CREATE_INDEX=false \\
            CREATE_MD5_FILE=false
        """
    else
        """
        cp $sample_bam ${sample}_sorted.bam
        """
}

// Mark duplicates for all merged samples bam files
process markDuplicate {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/metrics", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".dup_metrics") > 0 ? filename : null }

    input:
    set val(sample), file(sorted_bam) from samples_sorted_bam

    output:
    set val(sample), file("${sample}_markdup.bam") into samples_markdup_bam
    file "${sample}.dup_metrics" into dup_metric_files

    script:
    """
    picard MarkDuplicates \\
        INPUT=$sorted_bam \\
        OUTPUT=${sample}_markdup.bam \\
        METRICS_FILE=${sample}.dup_metrics \\
        TMP_DIR=tmp \\
        VALIDATION_STRINGENCY=SILENT \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=false \\
        MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \\
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \\
        SORTING_COLLECTION_SIZE_RATIO=0.25 \\
        READ_NAME_REGEX=\"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \\
        VERBOSITY=INFO \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        CREATE_INDEX=false \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Recalibrate the bam file with known variants
process recalibrate {
    tag "$sample"

    input:
    set val(sample), file(markdup_bam) from samples_markdup_bam

    output:
    set val(sample), file("${sample}_recal.bam"), file("${sample}_recal.bai") into samples_recal_bam
    file '.command.log' into gatk_base_recalibration_results

    script:
    """
    gatk -T BaseRecalibrator \\
        -I $markdup_bam \\
        -R $params.gfasta \\
        -o ${sample}_table.recal \\
        -cov ReadGroupCovariate \\
        -cov QualityScoreCovariate \\
        -cov CycleCovariate \\
        -cov ContextCovariate \\
        -U \\
        -OQ \\
        --default_platform illumina \\
        --knownSites $params.dbsnp \\
        -l INFO

    gatk -T PrintReads \\
        -BQSR ${sample}_table.recal \\
        -I $markdup_bam \\
        -R $params.gfasta \\
        -o ${sample}_recal.bam \\
        -baq RECALCULATE \\
        -U \\
        -OQ \\
        -l INFO
    """
}

// Realign the bam files based on known variants
process realign {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/alignment", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/realign/, "sorted_dupmarked_recalibrated_realigned")}

    input:
    set val(sample), file(recal_bam), file(recal_bam_ind) from samples_recal_bam

    output:
    set val(sample), file("${sample}_realign.bam"), file("${sample}_realign.bai") into bam_vcall, bam_phasing, bam_metrics
    file '.command.log' into gatk_base_realign_results

    script:
    """
    gatk -T RealignerTargetCreator \\
        -I $recal_bam \\
        -R $params.gfasta \\
        -o ${sample}_realign.intervals \\
        --known $params.dbsnp \\
        -l INFO

    gatk -T IndelRealigner \\
        -I $recal_bam \\
        -R $params.gfasta \\
        -targetIntervals ${sample}_realign.intervals \\
        -o ${sample}_realign.bam \\
        -l INFO
    """
}

// Calculate certain metrics
process calculateMetrics {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/metrics", mode: 'copy'

    input:
    set val(sample), file(aligned_bam), file(aligned_bam_ind) from bam_metrics

    output:
    file("*{metrics,pdf}") into metric_files

    script:
    """
    picard CollectAlignmentSummaryMetrics \\
        INPUT=$aligned_bam \\
        OUTPUT=${sample}.align_metrics \\
        REFERENCE_SEQUENCE=$params.gfasta \\
        VALIDATION_STRINGENCY=SILENT \\
        MAX_INSERT_SIZE=100000 \\
        ASSUME_SORTED=true \\
        IS_BISULFITE_SEQUENCED=false \\
        METRIC_ACCUMULATION_LEVEL="ALL_READS" \\
        STOP_AFTER=0 \\
        VERBOSITY=INFO \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''

    picard CollectInsertSizeMetrics \\
        HISTOGRAM_FILE=${sample}_insert.pdf \\
        INPUT=$aligned_bam \\
        OUTPUT=${sample}.insert_metrics \\
        VALIDATION_STRINGENCY=SILENT \\
        DEVIATIONS=10.0 \\
        MINIMUM_PCT=0.05 \\
        STOP_AFTER=0 \\
        METRIC_ACCUMULATION_LEVEL="ALL_READS" \\
        ASSUME_SORTED=true \\
        VERBOSITY=INFO \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''

    picard CollectHsMetrics \\
        BAIT_INTERVALS=$params.bait \\
        TARGET_INTERVALS=$params.target \\
        INPUT=$aligned_bam \\
        OUTPUT=${sample}.hs_metrics \\
        METRIC_ACCUMULATION_LEVEL="ALL_READS" \\
        VERBOSITY=INFO \\
        VALIDATION_STRINGENCY=SILENT \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Call variants
process variantCall {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/variants/, "raw_variants")}

    input:
    set val(sample), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(sample), file("${sample}_variants.vcf"), file("${sample}_variants.vcf.idx") into raw_variants

    script:
    """
    gatk -T HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -o ${sample}_variants.vcf \\
        --annotation HaplotypeScore \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --standard_min_confidence_threshold_for_calling 30.0 \\
        --dbsnp $params.dbsnp -l INFO
    """
}

// Select variants
process variantSelect {
    tag "$sample"

    input:
    set val(sample), file(raw_vcf), file(raw_vcf_idx) from raw_variants

    output:
    set val(sample), file("${sample}_snp.vcf"), file("${sample}_snp.vcf.idx") into raw_snp
    set val(sample), file("${sample}_indels.vcf"), file("${sample}_indels.vcf.idx") into raw_indels

    script:
    """
    gatk -T SelectVariants \\
        -R $params.gfasta \\
        --variant $raw_vcf \\
        --out ${sample}_snp.vcf \\
        --selectTypeToInclude SNP

    gatk -T SelectVariants \\
        -R $params.gfasta \\
        --variant $raw_vcf \\
        --out ${sample}_indels.vcf \\
        --selectTypeToInclude INDEL \\
        --selectTypeToInclude MIXED \\
        --selectTypeToInclude MNP \\
        --selectTypeToInclude SYMBOLIC \\
        --selectTypeToInclude NO_VARIATION
    """
}

// Filter SNP
process filterSnp {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(raw_snp), file(raw_snp_idx) from raw_snp

    output:
    set val(sample), file("${sample}_filtered_snp.vcf"), file("${sample}_filtered_snp.vcf.idx") into filtered_snp

    script:
    """
    gatk -T VariantRecalibrator \\
        -R $params.gfasta \\
        --input $raw_snp \\
        --maxGaussians 4 \\
        --recal_file ${sample}_snp.recal \\
        --tranches_file ${sample}_snp.tranches \\
        -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $params.hapmap \\
        -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $params.omni \\
        -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $params.dbsnp \\
        --mode SNP \\
        -an QD \\
        -an FS \\
        -an MQ

    gatk -T ApplyRecalibration \\
        -R $params.gfasta \\
        --out ${sample}_filtered_snp.vcf \\
        --input $raw_snp \\
        --mode SNP \\
        --tranches_file ${sample}_snp.tranches \\
        --recal_file ${sample}_snp.recal
    """
}

// Filter indels
process filterIndel {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(raw_indel), file(raw_indel_idx) from raw_indels

    output:
    set val(sample), file("${sample}_filtered_indels.vcf"), file("${sample}_filtered_indels.vcf.idx") into filtered_indels

    script:
    """
    gatk -T VariantFiltration \\
        -R $params.gfasta \\
        --variant $raw_indel \\
        --out ${sample}_filtered_indels.vcf \\
        --filterName GATKStandardQD \\
        --filterExpression "QD < 2.0" \\
        --filterName GATKStandardReadPosRankSum \\
        --filterExpression "ReadPosRankSum < -20.0" \\
        --filterName GATKStandardFS \\
        --filterExpression "FS > 200.0"
    """
}

// Group filted snp and indels for each sample
filtered_snp
    .cross(filtered_indels)
    .map{ it -> [it[0][0], it[0][1], it[0][2], it[1][1], it[1][2]] }
    .set{ variants_filtered }

// Combine filtered snp and indels for each sample
process combineVariants {
    tag "$sample"

    input:
    set val(sample), file(fsnp), file(fsnp_idx), file(findel), file(findel_idx) from variants_filtered

    output:
    set val(sample), file("${sample}_combined_variants.vcf"), file("${sample}_combined_variants.vcf.idx") into combined_variants

    script:
    """
    gatk -T CombineVariants \\
        -R $params.gfasta \\
        --out ${sample}_combined_variants.vcf \\
        --genotypemergeoption PRIORITIZE \\
        --variant:${sample}_SNP_filtered $fsnp \\
        --variant:${sample}_indels_filtered $findel \\
        --rod_priority_list ${sample}_SNP_filtered,${sample}_indels_filtered
    """
}


// Group filted bam and vcf for each sample for phasing
bam_phasing
    .cross(combined_variants)
    .map{ it -> [it[0][0], it[0][1], it[0][2], it[1][1], it[1][2]] }
    .set{ files_for_phasing }

// Indetifying haplotypes and create phasing between them
process haplotypePhasing {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(bam), file(bam_ind), file(vcf), file(vcf_ind) from files_for_phasing

    output:
    set val(sample), file("${sample}_combined_phased_variants.vcf"), file("${sample}_combined_phased_variants.vcf.idx") into vcf_eval, vcf_anno

    script:
    """
    gatk -T ReadBackedPhasing \\
        -R $params.gfasta \\
        -I $bam \\
        --variant $vcf \\
        --out ${sample}_combined_phased_variants.vcf
    """
}

// Evaluate variants
process variantEvaluate {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(phased_vcf), file(phased_vcf_ind) from vcf_eval

    output:
    file "${sample}_combined_phased_variants.eval"
    file "${sample}_combined_phased_variants.eval" into gatk_variant_eval_results

    script:
    """
    gatk -T VariantEval \\
        -R $params.gfasta \\
        --eval $phased_vcf \\
        --dbsnp $params.dbsnp \\
        -o ${sample}_combined_phased_variants.eval \\
        -L $params.target \\
        --doNotUseAllStandardModules \\
        --evalModule TiTvVariantEvaluator \\
        --evalModule CountVariants \\
        --evalModule CompOverlap \\
        --evalModule ValidationReport \\
        --stratificationModule Filter \\
        -l INFO
    """
}

////////////////////////////////////////////////////////////////////////////////////////////
//// Download the GRCh37.75 database manually as this will not let you write in miniconda
// 	cd /shared/ucl/depts/cancer/apps/miniconda/3/share/snpeff-4.3.1t-1
//	wget https://kent.dl.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip
//	unzip snpEff_v4_3_GRCh37.75.zip
//// ALSO: this might fail with "java.lang.OutOfMemoryError: GC overhead limit exceeded" errors
// 	(although t should be ok with more than 1 core: https://github.com/bcbio/bcbio-nextgen/issues/1730)
//	or with "java.lang.OutOfMemoryError: Java heap space" errors, 
//	so use directly the snpEff.jar with -Xmx
/////////////////////////////////////////////////////////////////////////////////////////////

// Annotate variants
process variantAnnotate {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(phased_vcf), file(phased_vcf_ind) from vcf_anno
    val snpeffDb from Channel.value(params.metaFiles[ params.genome ].snpeffDb)

    output:
    file "*.{vcf,idx,snpeff}"
    file 'SnpEffStats.csv' into snpeff_results

    script:
    """
#    snpEff \\
#        -c /shared/ucl/depts/cancer/apps/miniconda/3/share/snpeff-4.3.1t-1/snpEff.config \\
#        -i vcf \\
#        -o gatk \\
#        -o vcf \\
#        -filterInterval $params.target_bed GRCh37.75 $phased_vcf \\
#            > ${sample}_combined_phased_variants.snpeff

    java -jar -Xmx${task.memory.toGiga()}g  \\
	/shared/ucl/depts/cancer/apps/miniconda/3/share/snpeff-4.3.1t-1/snpEff.jar  \\
        -c /shared/ucl/depts/cancer/apps/miniconda/3/share/snpeff-4.3.1t-1/snpEff.config \\
        -csvStats SnpEffStats.csv \\
        -i vcf \\
        -o gatk \\
        -o vcf \\
        -filterInterval $params.target_bed ${snpeffDb} $phased_vcf \\
            > ${sample}_combined_phased_variants.snpeff

    gatk -T VariantAnnotator \\
        -R $params.gfasta \\
        -A SnpEff \\
        --variant $phased_vcf \\
        --snpEffFile ${sample}_combined_phased_variants.snpeff \\
        --out ${sample}_combined_phased_annotated_variants.vcf
    """
}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

NGI_exoseq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'NGI_exoseq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Genome']       = params.genome
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="




def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'NGI-exoseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'NGI-exoseq Workflow Summary'
    section_href: 'https://github.com/senthil10/NGI-ExoSeq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    picard MarkDuplicates --version &> v_picard.txt  || true
    gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('picard/*') from dup_metric_files.collect()
    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results.collect()
    file ('snpEff/*') from snpeff_results.collect()
    file ('gatk_variant_eval/*') from gatk_variant_eval_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


/*
 * Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[NGI_exoseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[NGI_exoseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[NGI_exoseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[NGI_exoseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NGI_exoseq] Pipeline Complete"

}

