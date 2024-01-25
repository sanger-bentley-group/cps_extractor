# CPS Extractor Pipeline <!-- omit in toc -->

This pipeline is in the early stages of development and is not fully tested!

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-23.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/singularity/)

The CPS extractor pipeline is a Nextflow pipeline designed for processing Streptococcus pneumoniae reads to extract the capsular locus sequence (CPS) and check for disruptive mutations

The pipeline is designed to be easy to set up and use, and is suitable for use on local machines and high-performance computing (HPC) clusters alike. Once you have downloaded the necessary docker/singularity images the pipeline can be used offline unless you have changed the selection of any database or container image.

The development of this pipeline is part of the GPS Project ([Global Pneumococcal Sequencing Project](https://www.pneumogen.net/gps/)). 

&nbsp;
# Table of contents <!-- omit in toc -->
- [Workflow](#workflow)
- [Usage](#usage)
  - [Requirements](#requirements)
  - [Accepted Inputs](#accepted-inputs)
  - [Profile](#profile)
  - [Options](#options)




&nbsp;
# Usage
## Requirements
- A POSIX-compatible system (e.g. Linux, macOS, Windows with [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)) with Bash 3.2 or later
- Java 11 or later (up to 21) ([OpenJDK](https://openjdk.org/)/[Oracle Java](https://www.oracle.com/java/))
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/)
  - For Linux, [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/) or [Docker Engine](https://docs.docker.com/engine/) is recommended over [Docker Desktop for Linux](https://docs.docker.com/desktop/). The latter is known to cause permission issues when running the pipeline on Linux. 

## Accepted Inputs
- Only Illumina paired-end short reads are supported
- Each sample is expected to be a pair of raw reads following this file name pattern: 
  - `*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}`
    - example 1: SampleName_R1_001.fastq.gz, SampleName_R2_001.fastq.gz
    - example 2: SampleName_1.fastq.gz, SampleName_2.fastq.gz
    - example 3: SampleName_R1.fq, SampleName_R2.fq
    


## Setup 
1. Clone the repository (if Git is installed on your system)
    ```
    git clone https://github.com/Oliver-Lorenz-dev/cps_nf.git
    ```
    or 
    
    Download and unzip/extract the [latest release](https://github.com/Oliver-Lorenz-dev/cps_nf/releases)
2. Go into the local directory of the pipeline and it is ready to use without installation (the directory name might be different)
    ```
    cd cps_nf
    ```
3. Run the database setup to download all required additional files and container images, so the pipeline can be used at any time with or without the Internet afterwards.
    > ⚠️ Docker or Singularity must be running, and an Internet connection is required.
    - Using Docker as the container engine
      ```
      ./run_cps_extractor --setup
      ```
    - Using Singularity as the container engine
      ```
      ./run_cps_extractor --setup -profile singularity
      ```

## Run
> ⚠️ Docker or Singularity must be running.
<!-- -->
> ℹ️ By default, Docker is used as the container engine and all the processes are executed by the local machine. See [Profile](#profile) for details on running the pipeline with Singularity or on a HPC cluster.
- You can run the pipeline without options. It will attempt to get the raw reads from the default location (i.e. `input` directory inside the `cps_nf` local directory)
  ```
  ./run_cps_extractor
  ```
- You can also specify the location of the raw reads by adding the `--input` option
  ```
  ./run_cps_extractor --input /path/to/raw-reads-directory
  ```

## Options
  ```
  |Usage:
  |./run_cps_extractor [option] [value]
  |
  |--input [PATH]                  Path to the input directory that contains reads to be processed. Default: ./input
  |--output [PATH]                 Path to the output directory that save the results. Default: output
  |--serotype [STR]                Serotype (if known). Default: None
  |--setup                         Alternative workflow for setting up the required databases.
  |--version                       Alternative workflow for getting versions of pipeline, container images, tools and databases
  |--help                          Print this help message
  ```


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
  | `lsf` | **The pipeline should be launched from a LSF cluster head node with this profile.** <br>Singularity is used as the container engine. <br> Processes are submitted to your LSF cluster via `bsub` by the pipeline. <br> (Tested on Wellcome Sanger Institute farm5 LSF cluster only) <br>

## Resume
- If the pipeline is interrupted mid-run, Nextflow's built-in `-resume` option can be used to resume the pipeline execution instead of starting from scratch again
- You should use the same command of the original run, only add `-resume` at the end (i.e. all pipeline options should be identical) 
  > ℹ️ `-resume` is a built-in Nextflow option, it only has one leading `-`
  - If the original command is
    ```
    ./run_cps_extractor --input /path/to/raw-reads-directory
    ```
  - The command to resume the pipeline execution should be
    ```
    ./run_cps_extractor --input /path/to/raw-reads-directory -resume
    ```

## Clean Up
- During the run of the pipeline, Nextflow generates a considerable amount of intermediate files
- If the run has been completed and you do not intend to use the `-resume` option or those intermediate files, you can remove the intermediate files using one of the following ways:
  - Run the included `clean_pipeline` script
    - It runs the commands in manual removal for you
    - It removes the `work` directory and log files within the `cps_nf` local directory
    ```
    ./clean_pipeline
    ```
  - Manual removal 
    - Remove the `work` directory and log files within the `cps_nf` local directory
    ```
    rm -rf work
    rm -rf .nextflow.log*
    ```
  - Run `nextflow clean` command
    - This built-in command cleans up cache and work directories
    - By default, it only cleans up the latest run
    - For details and available options of `nextflow clean`, refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html#clean)
    ```
    ./nextflow clean
    ```
    
&nbsp;
# Pipeline Options
- The tables below contain the available options that can be used when you run the pipeline
- Usage:
  ```
  ./run_pipeline [option] [value]
  ```
> ℹ️ To permanently change the value of an option, edit the `nextflow.config` file inside the `cps_nf` local directory.
<!-- -->
> ℹ️ `$projectDir` is a [Nextflow built-in implicit variables](https://www.nextflow.io/docs/latest/script.html?highlight=projectdir#implicit-variables), it is defined as the local directory of `gps-pipeline`.
<!-- -->
> ℹ️ Pipeline options are not built-in Nextflow options, they are lead with `--` instead of `-`

## Alternative Workflows
  | Option      | Values | Description |
-------------| --- | ---| --- |
  | `--setup`   | `true` or `false`<br />(Default: `false`) | Use alternative workflow for initialisation, which means downloading all required additional files and container images, and creating databases.<br />Can be enabled by including `--setup` without value. |
  | `--version` | `true` or `false`<br />(Default: `false`)| Use alternative workflow for showing versions of pipeline, container images, tools and databases.<br />Can be enabled by including `--version` without value.<br /> (This workflow pulls the required container images if they are not yet available locally) |
  | `--help`    | `true` or `false`<br />(Default: `false`)| Show help message.<br />Can be enabled by including `--help` without value. |

