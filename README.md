# TBSEQ.cat-nf

## Introduction

TBSEQ.cat-nf is for the classification and clustering and SNP barcoding of *Mycobacterium tuberculosis* genomes for surveillance of outbreaks. The pipeline consists of 4 sub-workflows. (1) Single genomes analysis (<code>--single</code>), (2) Pairwise genome comparisons (<code>--pairwise</code>), (3) Analysis summary (<code>--summary</code>), and (4) Cluster barcoding (<code>--barcoding</code>). 

This pipeline is build to be *additive*, meaning new genomes can be analyses and appended to the database (BBDD), facilitating the continuing monitoring of TB outbreaks. Because of this, it is important to ensure that sufficient storage is available for all the data. I am still in the process of optimising the whole process to ensure that the most essential files are retained while being consiousnt of storage. A databse consisting of 1,300 genomes will require ~ 600 GB (excluding the reads).

![image](png/pipeline.png)

#### <u>sub-wf 1: Single genomes analsysis</u> (<code>--single</code>)
The workflow partitions *Mycobacterium tuberculosis* complex r-eads for downstream processing (removing potentially contaminating reads), and utilizes TB specific tools [TB-Profiler](https://github.com/jodyphelan/TBProfiler) and [MTBSeq](https://github.com/ngs-fzb/MTBseq_source) for resistance profiling, lineage classification and variant calling. This will create a database (BBDD) of analysed genomes that will be utulised for downstream analysis. 

#### <u>sub-wf 2: Pairwise genome analysis</u> (<code>--pairwise</code>)  [WIP]
All genomes that have been processed through sub-wf 1 will be partitioned into their sub-lineages at two levels (e.g. Lineage4.2), as designated by TB-profile SNP lineage barcoding. Within these groups the genomes will complete the MTBSeq analysis begining with the joining of the SNP profiles (the slowest step), generating SNP alignments, and clustering the genomes at SNP distance of 5, 10 and 15 (these distances can be modified in the <code>nextflow.config</code>.)
#### <u>sub-wf 3: Analysis (<code>--summary</code>)</u> [WIP]

#### <u>sub-wf 4: Cluster SNP barcoding</u> (<code>--barcoding</code>) [WIP]
This is an experimental aspect of the workflow that aims to begin characterising individual SNPs that are designated uniquely to a particular 15 SNP distance cluster (--distance 15 in MTBseq). The plan with this sub-wf is to quickly being identifying which genomic cluster a particular genome may belonge to prior to SNP clustering with the goal of reducing computational resources and speeding up the analysis. All genomes as part of the sub-wf 1 will have their SNP profiles compared to the cluster barcode SNPs and pre allocated a preliminary cluster for clustering in sub-wf 2.

#### <u>sub-utilities:</u> [WIP]

#### Removing genomes from database
Occasionally you may find that a genome has gone through the entire workflow that should not have, or you sequence the genome more and now have more reads and want to re-analyse it. In these situations the genome needs to be completely erased from the database (all intermediate files, all records of this genome.). The sub-utility erase (<code>--erase</code>) performs this. You provide a csv file containing identifies of genomes you want remove, same as the sample <code>CSV</code> for running the workflow, and RutiSeq will remove all record of this genome.

## Installation

Requirements:
The following software needs to be available on your path
- Nextflow
  - On the IGTP HPC conda is available in /soft/bin/nextflow, just add the path to your <code>~/.bashrc</code>. 
- Singularity
  - On the IGTP HPC conda is available in /soft/bin/singularity, just add the path to your <code>~/.bashrc</code>. You will also need to specify a cache that nextflow stores all the images (this is useful if you are re-running the workflow too as it wont have to rebuild each time and just access it from the cache.)
- Miniconda3/Anaconda/Condaforge - however you want <code>conda</code> to be available

1. First you need to clone the github repository and pull all the necessary scripts
```{sh}
git clone https://github.com/phesketh-igtp/TBSEQ.cat-nf.git
```
1. Download a kaiju database and add full paths to the <code>nextflow.config</code> file
```{sh}
# Example for downloading the pre-build kaiju RefSeq-rn (bacteria/archaea and fungi only)
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_refseq_ref_2023-07-05.tgz
```
Modify the nextflow.config file to update the paths. You will need to modify <code>kaiju_fmi</code>, <code>kaiju_nodes</code>, and <code>kaiju_names</code> with the path of the kaiju database being utilised.

![image](png/figure1.png)
3. If using conda, you will need to create the conda enviornments Create all the necessay onda environe
```{sh}
```

## Usage

On the IGTP HPC, the nextflow scripts need to be submitted as a job, and from there nextflow will spawn all the jobs as part of its workflow. Bear in mind that when you run this a lot of your fair-use weight will be occupied, so you are unlikely to be able to run anything else for some time. 

To run the workflow test, run the following:
```{sh}
qsub -S /bin/bash -cwd -V -N nf-main \
        -o qsub-nf.out -l mem_free=6G  \
        /path/to/RutiSeq-nf/submit-nf.sh \
        /path/to/RutiSeq-nf/main.nf \
        --samplesheet /path/to/RutiSeq-nf/test/samples.hpc.csv \
        --outdir /path/to/RutiSeq-nf/RutiSeq-test \
        -profile igtp,singularity_on # this specifies that the job should be submitted to the IGTP HPC using singularity, can also designate it to use conda with 'conda_on'
```

## Inputs

**<u>Sample sheet:</u>** This must be a <code>CSV</code> file that contains the following information:
| detailed_name | alias | forward_path | reverse_path |
| ------------- | ----- | ------------ | ------------ |
| sample1_XXX-AAA- | 1-XXX-AAA_LX | /path/to/R1.fastq.gz | /path/to/R2.fastq.gz |
| sample1_XXX-AAA- | 2-XXX-AAA_LX | /path/to/R1.fastq.gz | /path/to/R2.fastq.gz |
  
  The alias MUST have the following structure - **[SampleID]_[LibraryID]**. The sampleID can be have information separated by hyphens ( - ), BUT the undescore ( \_ ) must be reserve by distinguishing between the SampleID and the LibraryID. This is to satisfy a data naming convention requires by MTBseq, which requires reads are name: **\[SampleID]\_[LibID]\_[\*]_[Direction].f(ast)q.gz**.

  Retaining the original name of your sample in the file is for your own records and wont be used in the pipeline.

  The compelte path of the sample sheet use used for the following flag <code>--samplesheet /path/to/samples.csv</code>. It is recommended that you keep a dated record of all the samples ran.

  If a sample in your sample sheet is duplicated, or that SampleID already exists in the database, then that sample will no be analyzed, and an alert will be produced informing you of this.

**<u>Metadata [optional]:</u>** Metadata, if provided, is only utilized at the summary step when the nexus files are generated for visualising the median-joining networks. If the following metadata is provided with the correct headers, then the nexus files will also contain relevant metadata.
'detection_date' is the estimated date for the onset of infection, or diagnosis. 'location' can be anything you want, the district of the patient, the hospital that performed the diagnosis, country of sample origin.

| detailed_name | alias | detection_date | location |
| ------------- | ----- | ------------ | ------------ |
| sample1_XXX-AAA- | 1-XXX-AAA_LX | 2024-01-01 | Hospital A |
| sample1_XXX-AAA- | 2-XXX-AAA_LX | 2022-01-02 | Hospital B |
