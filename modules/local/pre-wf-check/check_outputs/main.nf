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

cat <<'EOF' > check_script.sh

#!/bin/bash
set -e

sampleID="\$1"
outdir="\$2"

all_outputs_exist=true

if [ ! -d "\${outdir}/bbdd/mtbseq/samples/\${sampleID}" ]; then
    all_outputs_exist=false
fi

if [ ! -f "\${outdir}/bbdd/mtbseq/samples/\${sampleID}/Classification/Strain_Classification.tab" ]; then
    all_outputs_exist=false
fi

if [ ! -f "\${outdir}/bbdd/mtbseq/samples/\${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab" ]; then
    all_outputs_exist=false
fi

if [ ! -f "\${outdir}/bbdd/mtbseq/samples/\${sampleID}/SNP-Profiles/\${sampleID}.gatk.vcf.gz" ]; then
    all_outputs_exist=false
fi

if [ ! -f "\${outdir}/bbdd/tbprofiler/results/\${sampleID}.results.txt" ]; then
    all_outputs_exist=false
fi

if [ ! -f "\${outdir}/bbdd/tbprofiler/who-only/results/\${sampleID}.results.txt" ]; then
    all_outputs_exist=false
fi

echo "\$all_outputs_exist"
EOF

chmod +x check_script.sh
all_outputs_exist=\$(./check_script.sh "${sampleID}" "${params.outdir}")

echo "Sample ${sampleID}: all_outputs_exist=\$all_outputs_exist"
    """
}