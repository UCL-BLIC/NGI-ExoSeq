# nfcore/rnaseq Usage

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
module load blic-modules
module load nextflow

nextflow_rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh38
```

This will launch the pipeline with the `legion` or `myriad` configuration profile, depending on where you submit the job from.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Main Arguments

### `-profile`
This parameter is NOT necessary as the shortcut `nextflow_rnaseq` takes care of selecting the appropiate configuration profile. But just for your information, profiles are used to give 
configuration presets for different compute environments.

* `legion`
    * A generic configuration profile to be used with the UCL cluster legion
* `myriad`
    * A generic configuration profile to be used with the UCL cluster myriad

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

## Alignment tool
By default, the pipeline uses [BWA](https://github.com/lh3/bwa) to align the raw FastQ reads to the reference genome. BWA-MEM is the latest and it's generally recommended as it is faster and more accurate compared to the other available BWA algorithms.

## Reference Genomes

The pipeline config files come bundled with paths to the following reference genome assemblies: GRCh37, GRCh38, mm9_6J. The pipeline has this aprameter set up as `GRCh37` by default, you can change it using the `--genome` flag.

* Human
  * `--genome GRCh37`
  * `--genome GRCh38`
* Mouse
  * `--genome mm9_6J`

Note: for genome mm9_6J, the HaplotypeCaller and the BaseRecalibration steps are run without a prior set of variants (no `--dbSNP` flag).

### `--bwa_index`, `--gfasta`, `--gfasta_fai_ucsc` 
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--bwa_index '[path to BWA index]' \
--gfasta '[path to FastA reference]' \
--gfasta_fai_ucsc '[path to USCS FastA reference index - for bigwigs]'
```

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder.
These can then be used for future pipeline runs, reducing processing times.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveAlignedIntermediates`
As above, by default intermediate BAM files from the alignment will not be saved. The final BAM files created after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM files from BWA and sorting steps.

## Adapter Trimming
If specific additional trimming is required (for example, from additional tags),
you can use any of the following command line parameters. These affect the command
used to launch TrimGalore!

### `--clip_r1 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).

### `--clip_r2 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).

### `--three_prime_clip_r1 [int]`
Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.

### `--three_prime_clip_r2 [int]`
Instructs Trim Galore to remove bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.


## Capture Kits
### Installed kits:
* Human 
  * `agilent_v5_hg19` (Default)
  * `agilent_v7_hg19`
  * `agilent_v7_hg38`
  * `roche_SeqCap_Exome_v3_hg19`
  * `roche_SeqCap_Exome_v3_hg38`
  * `roche_SeqCap_MedExome_hg19`
  * `roche_SeqCap_MedExome_hg38`
* Mouse
  * `agilent_mouse_V1_mm9`

You can select the appropiate kit using the `--kit` argument. If the capture kit that was used to sequence your sample is not one of the above, you can prepare the kit following [these instructions](./kits.md).

Once you have the required `interval_list` and `bed` files, you can specify them using the foillowing parameters:

```bash
--bait                         Absolute path to bait file
--target                       Absolute path to target file
--target_bed                   Absolute path to target bed file
```

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits on UPPMAX with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--skip_markduplicates`
This will skip the removal of duplicates. Might be useful for amplicon data.

### `--skip_recalibration`
This will skip the BaseRecalibration step. This is only recommended for mouse strains for which the SNP annotation is poor. For mm9_6J, skiping the recalibration is *not* recommended, as B6/J does not have variants and therefore any difference with the reference genome should be considered a technical error.

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.
