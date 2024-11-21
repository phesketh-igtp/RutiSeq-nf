# phesketh-igtp/TBSEQ.cat-wf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of phesketh-igtp/TBSEQ.cat-wf, created with the [nf-core](https://nf-co.re/) template. Currently only performs a single sample analysis of genomes and does not perform the 

### `Added`

- The following modules that will make up the subworkflow <code>SINGLE_SUBWF</code>:
    - <code>MTBSEQ_SINGLE</code> -- currently doesnt work, deletes all the data when it runs - very very good
    - <code>TBPROFILER_DB</code>
    - <code>TBPROFILER_PROFILE_TBDB</code>
    - <code>TBPROFILER_PROFILE_WHO</code>
- Test data added, 4 samples sub-samples to 10% of the average number of reads (20K)

### `Fixed`

### `Dependencies`

### `Deprecated`
