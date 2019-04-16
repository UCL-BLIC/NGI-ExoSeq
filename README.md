# ![nf-core/ExoSeq](https://raw.githubusercontent.com/nf-core/Exoseq/master/docs/images/ExoSeq_logo.png)
[![Build Status](https://travis-ci.org/nf-core/ExoSeq.svg?branch=master)](https://travis-ci.org/nf-core/ExoSeq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/exoseq)

## Introduction

**UCL-BLIC/ENGI-exoseq** is a bioinformatics analysis pipeline that performs best-practice analysis pipeline for Exome Sequencing data.

The pipeline is built based on [GATK](https://software.broadinstitute.org/gatk/best-practices/) best practices using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. The main steps done by pipeline are the following (more information about the processes can be found [here](docs/processes.md)).

* Alignment
* Marking Duplicates
* Recalibration
* Realignment
* Variant Calling
* Variant Filtration
* Variant Evaluation
* Variant Annotation

## Documentation
The UCL-BLIC/NGI-exoseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Pipeline installation and configuration instructions](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
   * [Preparing custom exome capture kits](docs/kits.md)
3. [Output and how to interpret the results](docs/output.md)
4. [Troubleshooting](docs/troubleshooting.md)

