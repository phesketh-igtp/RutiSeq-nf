#!/bin/bash

### 

eval "$(conda shell.bash hook)"

# Configure bash
set -e          # exit immediately on error
set -u          # exit immediately if using undefined variables
set -o pipefail # ensure bash pipelines return non-zero status if any of their command fails

# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal
# to Nextlflow so that all jobs launched by the pipeline will be cancelled too.
_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}

trap _term TERM

# limit the RAM that can be used by nextflow
export NXF_JVM_ARGS="-Xms2g -Xmx5g"

# Run the pipeline. The command uses the arguments passed to this script, e.g:
#
# $ qsub -S /bin/bash -cwd -V -N nf-main -o qsub-nf.out -l mem_free=6G submit-nf.sh main.nf --samplesheet test/samples.hpc.csv --outdir RutiSeq -profile igtp,singularity_on
# Convenience::
# $ rm -rf qsub-nf.out .nextflow* nf-main.e6272*; qsub -S /bin/bash -cwd -V -N nf-main -o qsub-nf.out -l mem_free=6G submit-nf.sh main.nf --samplesheet test/samples.hpc.csv --outdir RutiSeq -profile igtp,singularity_on; sleep 2s; tail -f qsub-nf.out
nextflow run "$@" -ansi-log false & pid=$!

echo -e "Running:       
                nextflow run "$@" -ansi-log false & pid=$!
"

# Wait for the pipeline to finish
echo "Waiting for ${pid}"
wait $pid

# Return 0 exit-status if everything went well
exit 0

### SUBMIT TO HPC:
####
#### qsub -S /bin/bash -cwd -V -N nf-main -o qsub-nf.out -l mem_free=6G submit-nf.sh /imppc/labs/emlab/share/GitHub/RutiSeq-nf/main.nf -profile igtp,singularity_on --samplesheet sample-sheet/[XXXX].csv --runID [XXX]
