process CHECK_EXISTING_OUTPUTS {
    tag "${sampleID}"
    cache false

    input:
        tuple val(sampleID), path(forward), path(reverse)
    
    output:
        tuple val(sampleID), path(forward), path(reverse), env(all_outputs_exist)

    script:
    """
    #!/bin/bash
    set -e

    tbprofiler_tbdb_exists=0
    mtbseq_exists=0
    mtbseq_snp_exists=0

    tbprofiler_result="${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt"
    tbprofiler_who_result="${params.outdir}/bbdd/tbprofiler/who-only/results/${sampleID}.results.txt"
    mtbseq_classification="${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab"
    mtbseq_statistics="${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab"
    mtbseq_snp_profile="${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}.gatk.vcf.gz"

    echo "Checking files for ${sampleID}:"
    echo "TBProfiler result: \$tbprofiler_result"
    echo "TBProfiler WHO result: \$tbprofiler_who_result"
    echo "MTBseq classification: \$mtbseq_classification"
    echo "MTBseq statistics: \$mtbseq_statistics"
    echo "MTBseq SNP profile: \$mtbseq_snp_profile"

    if [[ -f "\$tbprofiler_result" ]] && [[ -f "\$tbprofiler_who_result" ]]; then
        tbprofiler_tbdb_exists=1
        echo "TBProfiler files exist"
    else
        echo "TBProfiler files do not exist"
    fi

    if [[ -f "\$mtbseq_classification" ]] && [[ -f "\$mtbseq_statistics" ]]; then
        mtbseq_exists=1
        echo "MTBseq files exist"
    else
        echo "MTBseq files do not exist"
    fi

    if [[ -f "\$mtbseq_snp_profile" ]]; then
        mtbseq_snp_exists=1
        echo "MTBseq SNP profile exists"
    else
        echo "MTBseq SNP profile does not exist"
    fi

    if [[ \$tbprofiler_tbdb_exists -eq 1 ]] && [[ \$mtbseq_exists -eq 1 ]] && [[ \$mtbseq_snp_exists -eq 1 ]]; then
        all_outputs_exist=true
    else
        all_outputs_exist=false
    fi

    echo "Sample: ${sampleID}"
    echo "tbprofiler_tbdb_exists: \$tbprofiler_tbdb_exists"
    echo "mtbseq_exists: \$mtbseq_exists"
    echo "mtbseq_snp_exists: \$mtbseq_snp_exists"
    echo "all_outputs_exist: \$all_outputs_exist"
    """
}