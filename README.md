# CPS Extractor Pipeline <!-- omit in toc -->

This pipeline is in the early stages of development and is not fully tested!

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-23.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/singularity/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/harryhung/gps-unified-pipeline)

The CPS extractor pipeline is a Nextflow pipeline designed for processing Streptococcus pneumoniae FASTA sequences to extract the capsular locus sequence (CPS) and check for disruptive mutations

The pipeline is designed to be easy to set up and use, and is suitable for use on local machines and high-performance computing (HPC) clusters alike. Once you have downloaded the necessary docker/singularity images the pipeline can be used offline unless you have changed the selection of any database or container image.

The development of this pipeline is part of the GPS Project ([Global Pneumococcal Sequencing Project](https://www.pneumogen.net/gps/)). 

&nbsp;
# Table of contents <!-- omit in toc -->
- [Workflow](#workflow)
- [Usage](#usage)
  - [Requirement](#requirement)
  - [Accepted Inputs](#accepted-inputs)
  - [Profile](#profile)
  - [Options](#options)




&nbsp;
# Usage
## Requirement
- A POSIX-compatible system (e.g. Linux, macOS, Windows with [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)) with Bash 3.2 or later
- Java 11 or later (up to 21) ([OpenJDK](https://openjdk.org/)/[Oracle Java](https://www.oracle.com/java/))
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/)
  - For Linux, [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/) or [Docker Engine](https://docs.docker.com/engine/) is recommended over [Docker Desktop for Linux](https://docs.docker.com/desktop/). The latter is known to cause permission issues when running the pipeline on Linux. 

## Accepted Inputs
- Only Illumina paired-end short reads are supported
- Each sample is expected to be a pair of raw reads following this file name pattern: 
  - `*.{fa,fasta}` 
    
  
## Profile
- By default, Docker is used as the container engine and all the processes are executed by the local machine. To change this, you could use Nextflow's built-in `-profile` option to switch to other available profiles
  > ℹ️ `-profile` is a built-in Nextflow option, it only has one leading `-`
  ```
  nextflow run . -profile [profile name]
  ```
- Available profiles: 
  | Profile Name | Details |
  | --- | --- |
  | `standard`<br> (Default) | Docker is used as the container engine. <br> Processes are executed locally. |
  | `singularity` |  Singularity is used as the container engine. <br> Processes are executed locally. |
  | `lsf` | **The pipeline should be launched from a LSF cluster head node with this profile.** <br>Singularity is used as the container engine. <br> Processes are submitted to your LSF cluster via `bsub` by the pipeline. <br> (Tested on Wellcome Sanger Institute farm5 LSF cluster only) |

## Options
  ```
  |Usage:
  |nextflow run . [option] [value]
  |
  |--bakta_threads [INT]           Threads used for bakta. Default: 4
  |--bakta_db [PATH]               Path to bakta database. Default: /data/pam/software/bakta/v5
  |--blastdb [PATH]                Path to blast database. Default: cps_blastdb
  |--input [PATH]                  Path to the input directory that contains the sequences to be processed {.fa,.fasta}. Default: input
  |--output [PATH]                 Path to the output directory that save the results. Default: output
  |--prodigal-training-file [PATH] Path to prodigal training file used in annotation. Default: all.trn
  |--version                       Alternative workflow for getting versions of pipeline, container images, tools and databases
  |--help                          Print this help message
  ```