#!/bin/bash

## Thank you to @dmdmckeow for the original script and inspriration
## @author: Dean A McKeown and Poppy J Hesketh-Best
## @version: v1.0.0
## @description: This script is used to submit a Nextflow pipeline to 
##               an HPC cluster using qsub. It sets up the environment, 
##               runs the pipeline, and handles job cancellation.
## @changelog:
##   v1.0.0-2024-11-01 : Initial version

green='\033[32m';red='\033[31m';cyan='\033[36m';purple='\033[35m';nocolor='\033[m'

### 

eval "$(conda shell.bash hook)"

# Configure bash
set -e          # exit immediately on error
set -u          # exit immediately if using undefined variables
set -o pipefail # ensure bash pipelines return non-zero status if any of their command fails

# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal
# to Nextflow so that all jobs launched by the pipeline will be cancelled too.
_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}

trap _term TERM

# limit the RAM that can be used by nextflow
export NXF_JVM_ARGS="-Xms2g -Xmx5g"

# Run the pipeline. The command uses the arguments passed to this script, e.g:

nextflow run "$@" -ansi-log false & pid=$!

echo -e "Running: nextflow run "$@" -ansi-log false & pid=$!\n"

echo -e "${red}$(date +'%d/%m/%Y %H:%M:%S')${nocolor}	qsub -S /bin/bash -cwd -V -N nf-main -o qsub-nf.out -l mem_free=6G submit-nf.sh "$@"" >> submit-nf.log

# Wait for the pipeline to finish
echo "Waiting for ${pid}"
wait $pid

# Return 0 exit-status if everything went well
exit 0