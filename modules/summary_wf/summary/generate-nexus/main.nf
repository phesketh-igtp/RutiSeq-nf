process GENERATE_NEXUS {

    conda params.snippy_env 

    tag "cluster: ${clusterID}"

    publishDir "${params.outDir}/results/${params.runID}/networks/", 
        mode: 'copy', 
        overwrite: true, 
        saveAs: { filename ->
                    if (filename.startsWith("fasta/") || filename.startsWith("positions/") || filename.startsWith("nexus/")) {
                        return filename
                    }
                    return null  // Don't publish anything else
                }

    input:
        tuple val(lineage), 
                val(clusterID),
                file(snp_fasta),
                file(snp_tab),
                file(clusters_tab)

    output:
        path("fasta/*"),        optional: true
        path("positions/*"),    optional: true
        path("nexus/*"),        optional: true

        tuple val(clusterID),
            path("fasta/${clusterID}_refseq.fasta"), emit: snp_fasta

        tuple val(lineage), 
                val(clusterID),
                file("fasta/${clusterID}_refseq.fasta"),
                file(clusters_tab),        emit: annotated_nexus_ch
        
        path("${clusterID}.handover.out"), emit: handover_out

    script:

        """
        # Create output directories
        mkdir -p nexus/ fasta/ positions/

        rm -rf ${params.outDir}/results/${params.runID}/networks/*/${clusterID}*

        # Get the list of genomes in the cluster
        grep -w ${clusterID} ${clusters_tab} \\
            | cut -f3 \\
            | awk -F'-' 'NF > 1 && \$2 != ""' \\
            > tmp.sampleIDs

        grep -f tmp.sampleIDs ${snp_fasta} \\
            | sed  's/>//g' \\
            | sort \\
            | uniq \\
            > genomes.list

        # Check the number of genomes matches
        if [ "\$(wc -l < genomes.list)" -ne "\$(wc -l < tmp.sampleIDs)" ]; then
            echo "ERROR: lengths are unequal — duplicated IDs creating issue" >&2
            echo "genomes.list: \$(wc -l < genomes.list)" >&2
            echo "tmp.sampleIDs: \$(wc -l < tmp.sampleIDs)" >&2
            exit 1
        fi

        # Separate fasta for the cluster
        seqkit grep -w 0 -f \\
            "genomes.list" \\
            ${snp_fasta} > \\
            ${clusterID}.fasta

        # Check the correct number of genomes were extracted
        if [ "\$(grep -c '>' ${clusterID}.fasta)" -ne "\$(wc -l < genomes.list)" ]; then
            echo "ERROR: fasta sequences extracted do not match expected number" >&2
            echo "Extracted: \$(grep -c '>' ${clusterID}.fasta)" >&2
            echo "Expected: \$(wc -l < genomes.list)" >&2
            exit 1
        fi
        
        # Try to run snp-sites and check for SNPs
        if snp-sites ${clusterID}.fasta > ${clusterID}.snpsites.fasta 2>/dev/null && [[ -s ${clusterID}.snpsites.fasta ]]; then
            echo "Processing ${clusterID} with SNPs detected"
            
            # Normal processing with SNPs
            snp-sites ${clusterID}.fasta -v \\
                | cut -f2 \\
                | sed '1,4d' \\
                > positions/${clusterID}_positions.tab

            # H37Rv reference
            awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$3}' \\
                positions/${clusterID}_positions.tab \\
                ${snp_tab} > ${clusterID}_tmp_refseq

            paste -s -d "" ${clusterID}_tmp_refseq \\
                | sed '1i >H37Rv' \\
                > ${clusterID}_H37Rv.fasta

            # Genomic positions
            awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$1}' \\
                positions/${clusterID}_positions.tab "${snp_tab}" \\
                > positions/${clusterID}_genomic_positions.tab

            # MTBC ancestor processing
            cp ${params.mtbc_ancestor_path} tmp.MTB_anc.pos.gz
            gunzip tmp.MTB_anc.pos.gz

            awk 'NR==FNR {pos[\$1]; next} FNR in pos {print \$3}' \\
                positions/${clusterID}_genomic_positions.tab tmp.MTB_anc.pos \\
                > ${clusterID}_tmp_MTB_anc

            paste -s -d '' ${clusterID}_tmp_MTB_anc \\
                | sed '1i >MTB_anc' \\
                > ${clusterID}_MTB_anc.fasta

            # Combine FASTA files
            cat ${clusterID}.snpsites.fasta \\
                ${clusterID}_H37Rv.fasta \\
                ${clusterID}_MTB_anc.fasta \\
                > fasta/${clusterID}_refseq.fasta

            # Convert to NEXUS
            seqret -osformat2 nexus \\
                -sequence fasta/${clusterID}_refseq.fasta \\
                -outseq nexus/${clusterID}_refseq.nex

        else
            echo "No SNPs detected for ${clusterID}, creating minimal output files"
            
            # Create empty/minimal files
            cp ${clusterID}.fasta fasta/${clusterID}_refseq.fasta
            touch ${clusterID}.snpsites.fasta
            touch positions/${clusterID}_positions.tab
            touch positions/${clusterID}_genomic_positions.tab
            touch nexus/${clusterID}_refseq.nex
        fi

        # Clean up
        rm -rf tmp.*

        # Simplify the name of the file
        #cat ${snp_tab} > ${clusterID}.snp.tab

        # Final checks (your existing code)
        # check if the nexus generation was successful
        if [[ ! -f positions/${clusterID}_genomic_positions.tab ]]; then 
            echo "No SNPs detected for ${clusterID}" > positions/${clusterID}_genomic_positions.tab
        fi

        if [[ ! -f positions/${clusterID}_genomic_positions.tab ]]; then 
            echo "No SNPs detected for ${clusterID}, nexus generation failed." > nexus/${clusterID}_refseq.nex
        fi

        echo > "${clusterID}.handover" > ${clusterID}.handover.out
        """

}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0.1
@description:
    This process generates a NEXUS file for each cluster of genomes
    using the SNPs from the reference genome and the MTB_anc
    (Valencian ancestor, Iñaki Comas 2013) as a reference. The resulting 
    NEXUS file is intended to be use with PopArt for visualising median-joining
    networks.
@changelog
    v1.0-2025-04-01: Initial version
    v1.0.1-2025-04-04: Moved the function out of a BASH script and into the Nextflow script block
                    Added argument for the publishDir to only publish relevant files, as this
                        was causing clashes with other processes in the workflow
                    Imporved handling for cluster with no SNPs detected (usually as they are identical)
        v1.0.2-2026-01-05: Added checks to ensure the number of genomes extracted matches expectations.
*/