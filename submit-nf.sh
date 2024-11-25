#!/bin/bash

### qsub -S /bin/bash -cwd -V -N nf-main -o qsub-nf.out -e qsub-nf.err -l h=sge-exec-12 -q d10imppcv3 -pe smp 4 -l mem_free=6G submit-nf.sh --samplesheet test/samples.csv --outdir RutiSeq


eval "$(conda shell.bash hook)"

export PATH=$PATH:/imppc/labs/emlab/phesketh/.local/bin/


# Configure bash
set -e          # exit immediately on error
set -u          # exit immidiately if using undefined variables
set -o pipefail # ensure bash pipelines return non-zero status if any of their command fails

# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal
# to Nextlflow so that all jobs launche by the pipeline will be cancelled too.
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
# $ sbatch submit_nf.sh nextflow/rnatoy -with-singularity
#
# will use "nextflow/rnatoy -with-singularity" as arguments
nextflow run "$@" --profile igtp -ansi-log false & pid=$!

echo -e "Running:       nextflow run "$@" --profile igtp -ansi-log false "

# Wait for the pipeline to finish
echo "Waiting for ${pid}"
wait $pid

# Return 0 exit-status if everything went well
exit 0