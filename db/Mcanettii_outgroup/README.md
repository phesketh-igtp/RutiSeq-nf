In this directory is the genome of *Mycobacterium canettii* (strain CIPT140010059), which is often used as an outgroup in phylogenetic analyses of the *Mycobacterium tuberculosis* complex (MTBC). This strain is considered to be one of the closest known relatives of the MTBC and is frequently employed to root phylogenetic trees and provide evolutionary context.

The genome sequence is provided in FASTA format, along with the outputs of an alignment against the reference strain H37Rv (NC_000962.3) using Snippy. Provided are both the aligned contigs (`NC_000962.aligned.fa`) and the variant calls (`NC_000962.vcf`). These files can be useful for comparative genomics and phylogenetic studies involving MTBC strains.

If you intent to use a different genome for the snippy, then run `snippy --ref <M.canettii_genome.gbk> --ctgs <your_contigs.fasta> --outdir . --prefix <reference> --force`, to generate the aligned fasta and vcf files for your reference.

Within the workflow, this genome is appended to the full consensus alignment from snippy-core and used as an outgroup in IQ-TREE phylogenetic analyses after identifying variant and invariant sites.

## References

If you are using this genome in your analyses, please ensure to cite the original publication.

Supply et al., 2013 "Genomic analysis of smooth tubercle bacilli provides insights into ancestry and pathoadaptation of Mycobacterium tuberculosis" (doi: https://doi.org/10.1038/ng.2744)