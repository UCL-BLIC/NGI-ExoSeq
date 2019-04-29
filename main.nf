#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

===============================================================================================
                     N G I - E X O S E Q    B E S T    P R A C T I C E
===============================================================================================
 Exome seq pipeline based on NGI_exoseq with FastQC, Trim_Galore, MultiQC, QualiMap  and GATK-4
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - preprocess (FASTQC, TrimGalore)
 - alignment (BWA)
 - recalibration (GATK4)
 - variant calling (GATK4)
 - variant evaluation (GATK3/SnpEff)
 - report (MultiQC) 
*/

// Package version
version = '1.0'

// Help message
helpMessage = """
===============================================================================
NGI-ExoSeq : Exome/Targeted sequence capture best practice analysis v${version}
	     (with FastQC, TrimGalore, MultiQC, QualiMap and GATK-4 steps)
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
--skip_markduplicates          Skip picard MarkDuplicates (useful for amplicons)
--skip_recalibration           Skip BaseRecalibration step (useful for genomes/mouse strains with poor SNP annotation)
--save_dedupBam                Save dedup BAM (not saved by default) 

Kit files:
--kit                          Kit used to prep samples [Default: 'agilent_v5_hg19']
--bait                         Absolute path to bait file
--target                       Absolute path to target file
--target_bed                   Absolute path to target bed file (snpEff compatible format)

Genome/Variation files:
--dbsnp                        Absolute path to dbsnp file
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
params.skip_markduplicates = false
params.skip_recalibration = false
params.save_dedupBam = false
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
//params.kit = 'agilent_v5_hg19'   #- defined in config
params.bait = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].bait ?: false : false
params.target = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target ?: false : false
params.target_bed = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target_bed ?: false : false
params.dbsnp = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].dbsnp ?: false : false
params.gfasta = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta ?: false : false
params.gfasta_fai_ucsc = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta_fai_ucsc ?: false : false
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
if (!params.metaFiles[ params.genome ] && ['gfasta', 'bwa_index', 'dbsnp'].count{ params[it] } != 3){
    exit 1, "Genome '${params.genome}' is not available in pre-defined config, so you need to provide all genome specific " +
            "files with options '--gfasta', '--bwa_index' and '--dbsnp'"
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
    file '*_fastqc.{zip,html}' into fastqc_results,fastqc_results2

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
    trimgalore_results2 = []
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
        file '*trimming_report.txt' into trimgalore_results,trimgalore_results2
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
        -t ${task.cpus} \\
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
    set val({sample - ~/(_S\d*)?(_L\d*)?$/}), file("${sample}_sorted.bam") into lanes_sorted_bam

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
    set val(sample), file("${sample}_sorted.bam"), file("${sample}_sorted.bai")  into samples_sorted_bam, samples_sorted_bam_bigwig

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
            CREATE_INDEX=true \\
            CREATE_MD5_FILE=false
        """
    else
        """
	echo "No need to merge $sample_bam as it was a single lane, just create the index"
	cp $sample_bam ${sample}_sorted.bam
	picard BuildBamIndex INPUT=$sample_bam OUTPUT=${sample}_sorted.bai
        """

}

// Get bigwigs

process bigwigs {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/alignment", mode: 'copy'

    input:
    set val(sample), file(raw_bam), file(raw_bam_ind) from samples_sorted_bam_bigwig

    output:
    file '*.bw'

    script:
    fasta=params.gfasta
    fastafai="${fasta}.fai"
    fastafaiucsc=params.gfasta_fai_ucsc
    if(params.genome == 'GRCh38')
  	  """
  	  bedtools genomecov -bg -ibam $raw_bam -g $fastafai > ${sample}_raw.bdg
  	  LC_COLLATE=C sort -k1,1 -k2,2n ${sample}_raw.bdg > ${sample}.raw_sorted.bdg
  	  bedGraphToBigWig ${sample}.raw_sorted.bdg $fastafai ${sample}_raw_sorted.bw
  	  """
    else
  	  """
  	  bedtools genomecov -bg -ibam $raw_bam -g $fastafai > ${sample}_raw.bdg
  	  LC_COLLATE=C sort -k1,1 -k2,2n ${sample}_raw.bdg > ${sample}.raw_sorted.bdg
  	  perl -p -i -e 's/^/chr/g' ${sample}.raw_sorted.bdg
  	  perl -p -i -e 's/chrMT/chrM/g' ${sample}.raw_sorted.bdg
  	  bedGraphToBigWig ${sample}.raw_sorted.bdg $fastafaiucsc ${sample}_raw_sorted.bw
  	  """
}


if(!params.skip_markduplicates){

	// Mark duplicates for all merged samples bam files
	process markDuplicate {
	    tag "$sample"
	    publishDir "${params.outdir}/${sample}/metrics", mode: 'copy',
	        saveAs: { filename -> filename.indexOf(".dup_metrics") > 0 ? filename : null }
	
	    input:
	    set val(sample), file(sorted_bam), file(sorted_bam_ind) from samples_sorted_bam
	
	    output:
	    set val(sample), file("${sample}_markdup.bam"), file("${sample}_markdup.bai") into samples_markdup_bam
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
	        CREATE_INDEX=true \\
	        MAX_RECORDS_IN_RAM=500000 \\
	        CREATE_MD5_FILE=false \\
	        GA4GH_CLIENT_SECRETS=''
	    """
	}

	/*
	 * Recalibrate BAM file with known variants and BaseRecalibrator
	 *
	*/

	process recalibrate {
	    tag "${sample}"
            publishDir "${params.outdir}/${sample}/alignment", mode: 'copy',
                saveAs: { filename -> 
			if(params.save_dedupBam && filename.indexOf(".bam") > 0) filename
			else if (params.save_dedupBam && filename.indexOf(".bai") > 0) filename
			else null
		}

	    input:
	    set val(sample), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam
	
	    output:
	    set val(sample), file("${sample}_recal.bam"), file("${sample}_recal.bai") into bam_vcall, bam_phasing, bam_metrics, bam_qualimap
	    file '.command.log' into gatk_base_recalibration_results,gatk_base_recalibration_results2
	
	    script:
	    dbSNP = params.dbsnp ? "--known-sites $params.dbsnp" : ''
	    if(!params.skip_recalibration){
		    """
		    gatk BaseRecalibrator \\
		        -I $markdup_bam \\
		        -R $params.gfasta \\
		        -O ${sample}_table.recal \\
		        -L $params.target \\
		        -ip 100 \\
			$dbSNP \\
		        --verbosity INFO \\
		        --java-options -Xmx${task.memory.toGiga()}g
	
		    gatk ApplyBQSR \\
		        -R $params.gfasta \\
		        -I $markdup_bam \\
		        --bqsr-recal-file ${sample}_table.recal \\
		        -O ${sample}_recal.bam \\
		        -L $params.target \\
		        -ip 100 \\
		        --create-output-bam-index true \\
		        --java-options -Xmx${task.memory.toGiga()}g
	
		    """
	    }else{
		    """
	            echo "*** SKIPPING RECALIBRATION ***"
		    cp $markdup_bam ${sample}_recal.bam 
		    cp $markdup_bam_ind ${sample}_recal.bai
		    """
	    }
	}

}else{

	// Recalibrate the bam file with known variants
	process recalibrate_wo_markdup {
	    tag "${sample}"

	    input:
	    set val(sample), file(sorted_bam), file(sorted_bam_ind) from samples_sorted_bam
	
	    output:
	    set val(sample), file("${sample}_recal.bam"), file("${sample}_recal.bai") into bam_vcall, bam_phasing, bam_metrics, bam_qualimap
	    file '.command.log' into gatk_base_recalibration_results,gatk_base_recalibration_results2
		
	    script:
	    dbSNP = params.dbsnp ? "--known-sites $params.dbsnp" : ''
	    if(!params.skip_recalibration){
		    """
		    gatk BaseRecalibrator \\
		        -I $sorted_bam \\
		        -R $params.gfasta \\
		        -O ${sample}_table.recal \\
		        -L $params.target \\
		        -ip 100 \\
		        $dbSNP \\
		        --verbosity INFO \\
		        --java-options -Xmx${task.memory.toGiga()}g
	
		    gatk ApplyBQSR \\
		        -R $params.gfasta \\
		        -I $sorted_bam \\
		        --bqsr-recal-file ${sample}_table.recal \\
		        -O ${sample}_recal.bam \\
		        -L $params.target \\
		        -ip 100 \\
		        --create-output-bam-index true \\
		        --java-options -Xmx${task.memory.toGiga()}g	
		    """
	    }else{
		    """
	            echo "*** SKIPPING RECALIBRATION ***"
		    cp $sorted_bam ${sample}_recal.bam 
		    cp $sorted_bam_ind ${sample}_recal.bai
		    """
	    }

	}
}

/*
 * Determine quality metrics of mapped BAM files using QualiMap 2
 *
*/
process qualiMap {
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/qualimap", mode: 'copy'

    input:
    set val(sample), file(recal_bam), file(recal_bam_ind) from bam_qualimap

    output:
    file "${sample}" into qualimap_results,qualimap_results2
    file '.command.log' into qualimap_stdout

    script:
    gcref = ''
    if(params.genome == 'GRCh37') gcref = '-gd HUMAN'
    if(params.genome == 'GRCh38') gcref = '-gd HUMAN'
    if(params.genome == 'mm9_6J') gcref = '-gd MOUSE'
    """
    qualimap bamqc $gcref \\
    -bam $recal_bam \\
    -outdir ${sample} \\
    --skip-duplicated \\
    --collect-overlap-pairs \\
    --outside-stats \\
    -nt ${task.cpus} \\
    -gff ${params.target_bed} \\
    --java-mem-size=${task.memory.toGiga()}G \\
    """
}


/*
 *  Calculate certain metrics
*/

process calculateMetrics {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/metrics", mode: 'copy'

    input:
    set val(sample), file(recal_bam), file(recal_bam_ind) from bam_metrics

    output:
    file("*{metrics,pdf}") into metric_files
    file "*metrics" into metric_results,metric_results2

    script:
    """
    java -jar -Xmx${task.memory.toGiga()}g \\
        /shared/ucl/apps/picard-tools/2.18.9/picard.jar CollectAlignmentSummaryMetrics \\
        INPUT=$recal_bam \\
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

    java -jar -Xmx${task.memory.toGiga()}g \\
        /shared/ucl/apps/picard-tools/2.18.9/picard.jar CollectInsertSizeMetrics \\
        HISTOGRAM_FILE=${sample}_insert.pdf \\
        INPUT=$recal_bam \\
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

    java -jar -Xmx${task.memory.toGiga()}g \\
	/shared/ucl/apps/picard-tools/2.18.9/picard.jar CollectHsMetrics \\
        BAIT_INTERVALS=$params.bait \\
        TARGET_INTERVALS=$params.target \\
        INPUT=$recal_bam \\
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
    set val(sample), file(recal_bam), file(recal_bam_ind) from bam_vcall

    output:
    set val(sample), file("${sample}_variants.vcf"), file("${sample}_variants.vcf.idx") into raw_variants_gvcf, raw_variants

    script:
    dbSNP = params.dbsnp ? "--dbsnp $params.dbsnp" : ''
    """
    gatk HaplotypeCaller \\
        -I $recal_bam \\
        -R $params.gfasta \\
        -O ${sample}_variants.vcf \\
        -ERC GVCF \\
        -L $params.target \\
        --create-output-variant-index \\
        -ip 100 \\
	--annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        $dbSNP \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}


/*
 * Genotype generate GVCFs using GATK's GenotypeGVCFs
 * 
*/ 

process genotypegvcfs{
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy' 

    input:
    set val(sample), file(raw_vcf), file(raw_vcf_idx) from raw_variants_gvcf

    output:
    set val(sample), file("${sample}_gvcf.vcf"), file("${sample}_gvcf.vcf.idx") into raw_gvcfs

    script:
    dbSNP = params.dbsnp ? "--dbsnp $params.dbsnp" : ''
    """
    gatk GenotypeGVCFs \\
    -R $params.gfasta \\
    -L $params.target \\
    -ip 100 \\
    $dbSNP \\
    -V $raw_vcf \\
    -O ${sample}_gvcf.vcf 
    """
}

// Select variants
process variantSelect {
    tag "$sample"

    input:
    set val(sample), file(raw_gvcf), file(raw_gvcf_idx) from raw_gvcfs

    output:
    set val(sample), file("${sample}_snp.vcf"), file("${sample}_snp.vcf.idx") into raw_snp
    set val(sample), file("${sample}_indels.vcf"), file("${sample}_indels.vcf.idx") into raw_indels

    script:
    """
    gatk SelectVariants \\
        -R $params.gfasta \\
        -V $raw_gvcf \\
        -O ${sample}_snp.vcf \\
        --select-type-to-include SNP

    gatk SelectVariants \\
        -R $params.gfasta \\
        -V $raw_gvcf \\
        -O ${sample}_indels.vcf \\
        --select-type-to-include INDEL \\
        --select-type-to-include MIXED \\
        --select-type-to-include MNP \\
        --select-type-to-include SYMBOLIC \\
        --select-type-to-include NO_VARIATION
    """
}


// Hard-filter SNPs
// thresholds suggested in https://software.broadinstitute.org/gatk/documentation/article?id=11097
// more info here: https://software.broadinstitute.org/gatk/documentation/article.php?id=6925
process filterSnp {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(raw_snp), file(raw_snp_idx) from raw_snp

    output:
    set val(sample), file("${sample}_filtered_snp.vcf"), file("${sample}_filtered_snp.vcf.idx") into filtered_snp

    script:
    """
    gatk VariantFiltration \\
        -R $params.gfasta \\
        -V $raw_snp \\
        -O ${sample}_filtered_snp.vcf \\
        --filter-name GATKStandardQD \\
        --filter-expression "QD < 2.0" \\
        --filter-name GATKStandardMQ \\
        --filter-expression "MQ < 40.0" \\
        --filter-name GATKStandardFS \\
        --filter-expression "FS > 60.0" \\
        --filter-name GATKStandardReadPosRankSum \\
        --filter-expression "ReadPosRankSum < -8.0"
    """
}

// Hard-filter indels
process filterIndel {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(sample), file(raw_indel), file(raw_indel_idx) from raw_indels

    output:
    set val(sample), file("${sample}_filtered_indels.vcf"), file("${sample}_filtered_indels.vcf.idx") into filtered_indels

    script:
    """
    gatk VariantFiltration \\
        -R $params.gfasta \\
        -V $raw_indel \\
        -O ${sample}_filtered_indels.vcf \\
        --filter-name GATKStandardQD \\
        --filter-expression "QD < 2.0" \\
        --filter-name GATKStandardReadPosRankSum \\
        --filter-expression "ReadPosRankSum < -20.0" \\
        --filter-name GATKStandardFS \\
        --filter-expression "FS > 200.0"
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
    set val(sample), file("${sample}_combined_variants.vcf"), file("${sample}_combined_variants.vcf.idx") into vcf_eval, vcf_anno

    script:
    """
    picard MergeVcfs \\
        O=${sample}_combined_variants.vcf \\
        I=$fsnp \\
        I=$findel \\
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
    file "${sample}_combined_phased_variants.eval" into gatk_variant_eval_results,gatk_variant_eval_results2

    script:
    dbSNP = params.dbsnp ? "--dbsnp $params.dbsnp" : ''
    """
    gatk -T VariantEval \\
        -R $params.gfasta \\
        --eval $phased_vcf \\
        $dbSNP \\
        -o ${sample}_combined_phased_variants.eval \\
        -L $params.target \\
        -ip 100 \\
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
//   Mouse from https://sourceforge.net/projects/snpeff/files/databases/v4_3/
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
    file "*.{vcf,id,genes.txt}"
    file '*_SnpEffStats.csv' into snpeff_results,snpeff_results2

    script:
    """
    java -jar -Xmx${task.memory.toGiga()}g  \\
	/shared/ucl/depts/cancer/apps/miniconda3/share/snpeff-4.3.1t-1/snpEff.jar  \\
        -c /shared/ucl/depts/cancer/apps/miniconda3/share/snpeff-4.3.1t-1/snpEff.config \\
        -csvStats ${sample}_SnpEffStats.csv \\
        -i vcf \\
        -o gatk \\
        -o vcf \\
	${snpeffDb} $phased_vcf \\
            > ${sample}_combined_phased_annotated_variants.vcf
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
summary['Skip Markdup'] = params.skip_markduplicates ? 'Yes' : 'No'
summary['Skip Recalibration'] = params.skip_recalibration ? 'Yes' : 'No'
summary['Save dedup BAM'] = params.save_dedupBam ? 'Yes' : 'No'
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
    file 'software_versions_mqc.yaml' into software_versions_yaml,software_versions_yaml2

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    picard MarkDuplicates --version &> v_picard.txt  || true
    qualimap --help > v_qualimap.txt
    gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}



if(!params.skip_markduplicates){

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
	    file ('metrics/*') from metric_results.collect()
	    file ('qualimap/*') from qualimap_results.collect()
	    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results.collect()
	    file ('snpeff/*') from snpeff_results.collect()
	    file ('gatk_variant_eval/*') from gatk_variant_eval_results.collect()
	    file software_versions from software_versions_yaml.collect()
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
}else{
	process multiqc_wo_markdup {
	    publishDir "${params.outdir}/MultiQC", mode: 'copy'
	
	    input:
	    file multiqc_config
	    file ('fastqc/*') from fastqc_results2.collect()
	    file ('trimgalore/*') from trimgalore_results2.collect()
	    file ('metrics/*') from metric_results2.collect()
	    file ('qualimap/*') from qualimap_results2.collect()
	    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results2.collect()
	    file ('snpeff/*') from snpeff_results2.collect()
	    file ('gatk_variant_eval/*') from gatk_variant_eval_results2.collect()
	    file software_versions from software_versions_yaml2.collect()
	    file workflow_summary from create_workflow_summary(summary)
	
	    output:
	    file "*multiqc_report.html" into multiqc_report2
	    file "*_data"
	
	    script:
	    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
	    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
	    """
	    multiqc -f $rtitle $rfilename --config $multiqc_config .
	    """
	}
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

