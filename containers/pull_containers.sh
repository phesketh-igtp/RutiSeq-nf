#!/bin/bash

singularity pull --name TBPROFILER.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data

singularity pull --name MTBC_READ_QC.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f00cd356ee92f5211e5941beeb4bcab6abfb341e0e5fa7ace8c043406c13381/data

singularity pull --name MTBSEQ.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data

singularity pull --name SNP_PROFILING.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/608c097132a7de8e156c452f40ea3b3fea6bf0a35b6988e4b2fe74d91524303f/data
