# Description of MTBc ancestral sequence

This directory contains the genome that a inferred ancestral sequence of the Mycobacterium tuberculosis complex (MTBC). I have little else but to assume this is the inferred ancestral sequence that is used a lot in the literature, but reporting of this in the originating paper and downward articles is often missing or incorrect. I had to quest a bit to find the right reference let alone the zenodo repository! The original article doesnt cite the sequence directly, and I cannot get access to the supplementary data from the original article, but I found the sequence on zenodo which cites the original article.

- Zendo link: https://zenodo.org/records/3497110
- Direct download link: https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1
- Original publication: Comas et al., 2010 "Human T cell epitopes of Mycobacterium tuberculosis are evolutionarily hyperconserved" (doi: https://doi.org/10.1038/ng.590)

> Genomic coordinates follow the reference strain H37Rv coordinates and are equivalent. 

## Directory contents

- `MTB_ancestor_reference.fasta`: The inferred ancestral sequence of the MTBC in FASTA format.
- `MTB_ancestor_reference.pos.gz`: The inferred ancestral sequence by positions.

## Sanity check against H37Rv

Frustratingly, when I have been using it in identify variant sites in Mtb clusters I have observed no difference between H37Rv (NC_000962.3) and this ancestral sequence, which made me suspicious wether I even had the right sequence (the sequence was initially given to me by a colleague). So I decided to verify this by comparing the ancestral sequence to H37Rv for my own sanity check.

I compared this reconstruction against the reference genome H37Rv using
`snippy --ref NC_000962.3 --ctgs MTB_ancestor_reference.fasta --outdir MTBc_ancestral_vs_H37Rv --cpus 4`.

The variants called are summarized below:

```
Snippy v4.6.0
   stats: biallelic
          no. left trimmed                      : 1
          no. right trimmed                     : 1
          no. left and right trimmed            : 3
          no. right trimmed and left aligned    : 0
          no. left aligned                      : 0

       total no. biallelic normalized           : 5

       multiallelic
          no. left trimmed                      : 0
          no. right trimmed                     : 0
          no. left and right trimmed            : 0
          no. right trimmed and left aligned    : 0
          no. left aligned                      : 0

       total no. multiallelic normalized        : 0

       total no. variants normalized            : 5
       total no. variants observed              : 1152
       total no. reference observed             : 0
```

# Reference

Comas, I., Chakravartti, J., Small, P. et al. Human T cell epitopes of Mycobacterium tuberculosis are evolutionarily hyperconserved. Nat Genet 42, 498–503 (2010). https://doi.org/10.1038/ng.590