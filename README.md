# TBSEQ.cat-nf

## Introduction

TBSEQ.cat-nf is for the classification and clustering and SNP barcoding of *Mycobacterium tuberculosis* genomes for surveillance of outbreaks. The pipeline consists of four sub-workflows, with the default being to run all workflows (<code>--workflow full</code>). (1) Single genomes analysis (<code>--workflow single</code>), (2) Pairwise genome comparisons (<code>--workflow pairwise</code>), (3) Analysis summary (<code>--workflow summary</code>), and (4) Cluster barcoding (<code>--workflow barcoding</code>). 

This pipeline is build to be *additive*, meaning new genomes can be analyses and appended to the database (BBDD), facilitating the continuing monitoring of TB outbreaks. Because of this, it is important to ensure that sufficient storage is available for all the data. I am still in the process of optimising the whole process to ensure that the most essential files are retained while being consiousnt of storage. A databse consisting of 1,300 genomes will require ~ 600 GB (excluding the reads).

![image](png/pipeline.png)

#### <u>Workflow 1: Single genomes analsysis</u> (<code>--workflow single</code>)
The workflow partitions *Mycobacterium tuberculosis* complex r-eads for downstream processing (removing potentially contaminating reads), and utilizes TB specific tools [TB-Profiler](https://github.com/jodyphelan/TBProfiler) and [MTBSeq](https://github.com/ngs-fzb/MTBseq_source) for resistance profiling, lineage classification and variant calling. This will create a database (BBDD) of analysed genomes that will be utulised for downstream analysis.

If a genome genome in your sample sheet already exists in the BBDD, that genome will not be re-analysed. The only way to reanlayse the genome is to remove it completely from the BBDD, a sub-utility is under construction.

#### <u>Workflow 2: Pairwise genome analysis</u> (<code>--workflow pairwise</code>)  [WIP]
During this sub-workflow, all genomes that have been processed through sub-wf 1 will be partitioned into their sub-lineages at two levels for Lineage 4 genomes (the majority), and at one level for all other lineages (until the lineage has more than 600 genomes, which is when it will be split at level 2), as designated by TB-profile SNP lineage barcoding. Within these groups the genomes will complete the MTBSeq analysis (<code>--TBjoin</code>, <code>--TBamend</code>, and <code>--TBgroup</code>), which includes joining SNP profiles, generating SNP alignments, and clustering the genomes at SNP distances (default; **5, 10, 15** - but these distances can be modified in the <code>nextflow.config</code>).

#### <u>Workflow 3: Analysis (<code>--workflow summary</code>)</u> [WIP]
During this sub-workflow all the outputs from sub-wf 1 and 2 are compile by their essential information, all the raw outputs still remain in the BBDD for detail consideration. The principal results will be generated in a multi-sheet <code>XLSX</code>. This will include (i) Summary of genomes (classification, mapping statistics, genome coverage and quality), (ii-iii) resistance profiles (using both TBDB and the WHO designations), (iv) genomic clusters at designated SNP distances.

Additional outputs include PDFs of SNP phylogeny (ML tree generated with IQ-Tree) colored by cluster identity, and NEX files for upload into [PopArt](https://popart.maths.otago.ac.nz/) for visualisation.

#### <u>Workflow 4: Cluster SNP barcoding</u> (<code>--workflow barcoding</code>) [WIP]
This is an experimental aspect of the workflow that aims to begin characterizing individual SNPs that are designated uniquely to a particular 15 SNP distance cluster (<code>--distance 15</code> in MTBseq). The plan with this sub-wf is to quickly being identifying which genomic cluster a particular genome may belong to prior to SNP clustering with the goal of reducing computational resources and speeding up the analysis. All genomes as part of the sub-wf 1 will have their SNP profiles compared to the cluster barcode SNPs and pre allocated a preliminary cluster for clustering in sub-wf 2.

In this workflow, all genomes SNP profiles merged into a single VCF (grouped by lineage), and the SNP profiles of genomes belonging to the same cluster are compared to all other genomes within the same lineage, to calculate the F~TS~ value (fixation index) for each SNP within the cluster population. SNPs that fulfill the following criteria are classified as a cluster specific SNP:
- F~TS~ = 1
- Minimum of 20 reads in both strands (20X cov)
- Minimum quality of 20
- Not annotated as: *PE/PPE/PGRS*; *maturase*; *phage*; or *13E12 repeat family protein*
- Not located in insertion sequences
- Not within InDels or in high density regions (> 3 SNPs in 10 bp)

The script will generate a <code>BED</code> file for each unique SNP and its cluster designation and the genomes will re-assessed by this barcoding method. A report will be generated on the quality of the predictions, wether the clsuter asignment was no better than random chance, or if there is consistency compared to the MTBSeq cluster assignment. In future all genomes will ideally be grouped by their predicted clsuters, speeding up the rate by which a genome cluster is identified.

#### <u>sub-utilities:</u> [WIP]

#### Removing genomes from database
Occasionally you may find that a genome has gone through the entire workflow that should not have, or you sequence the genome more and now have more reads and want to re-analyse it. In these situations the genome needs to be completely erased from the database (all intermediate files, all records of this genome.). The sub-utility erase (<code>--workflow erase</code>) performs this. You provide a csv file containing identifies of genomes you want remove, same as the sample <code>CSV</code> for running the workflow, and RutiSeq will remove all record of this genome.


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
2. Download a kaiju database and add full paths to the <code>nextflow.config</code> file
```{sh}
# Example for downloading the pre-build kaiju RefSeq-rn (bacteria/archaea and fungi only)
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_refseq_ref_2023-07-05.tgz
```
Modify the nextflow.config file to update the paths. You will need to modify <code>kaiju_fmi</code>, <code>kaiju_nodes</code>, and <code>kaiju_names</code> with the path of the kaiju database being utilised.

![image](png/figure1.png)
3. If using conda, you will need to create the conda enviornments Create all the necessary conda environent.

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

If you are adding new data to an existing database generated with this pipeline, the <code>--outdir</code> MUST be be given the path to that database.

## Inputs

**<u>Sample sheet:</u>** This must be a <code>CSV</code> file that contains the following information:
| name | alias | forward_path | reverse_path |
| ------------- | ----- | ------------ | ------------ |
| sample1_XXX-AAA- | 1-XXX-AAA_LX | /path/to/R1.fastq.gz | /path/to/R2.fastq.gz |
| sample1_XXX-AAA- | 2-XXX-AAA_LX | /path/to/R1.fastq.gz | /path/to/R2.fastq.gz |

```{sh}
$ cat samples.csv
name,alias,forward_path,reverse_path
sample1_XXX-AAA,1-XXX-AAA_LX,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
sample1_XXX-AAA,2-XXX-AAA_LX,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
```
  
The alias MUST have the following structure - **[SampleID]_[LibraryID]**. The sampleID can be have information separated by hyphens ( - ), BUT the undescore ( \_ ) must be reserve by distinguishing between the SampleID and the LibraryID. This is to satisfy a data naming convention requires by MTBseq, which requires reads are name: **\[SampleID]\_[LibID]\_[\*]_[Direction].f(ast)q.gz**.

If you sampleID's happen to follow this structure, you can just have the same value in **name** and **alias**.

Retaining the original name of your sample in the file is for your own records and wont be used in the pipeline.

The complete path of the sample sheet use used for the following flag <code>--samplesheet /path/to/samples.csv</code>. It is recommended that you keep a dated record of all the samples ran.

If a sample in your sample sheet is duplicated, or that SampleID already exists in the database, then that sample will no be analyzed, and an alert will be produced informing you of this.

**<u>Metadata [optional]:</u>** Metadata, if provided, is only utilized at the summary step when the nexus files are generated for visualising the median-joining networks. If the following metadata is provided with the correct headers, then the nexus files will also contain relevant metadata.
'detection_date' is the estimated date for the onset of infection, or diagnosis. 'location' can be anything you want, the district of the patient, the hospital that performed the diagnosis, country of sample origin.

| detailed_name | alias | detection_date | location |
| ------------- | ----- | ------------ | ------------ |
| sample1_XXX-AAA- | 1-XXX-AAA_LX | 2024-01-01 | Hospital A |
| sample1_XXX-AAA- | 2-XXX-AAA_LX | 2022-01-02 | Hospital B |

# To do 
- Test <code>--workflow pairwise</code>
- Write <code>--workflow summay</code> and <code>--workflow barcode</code>