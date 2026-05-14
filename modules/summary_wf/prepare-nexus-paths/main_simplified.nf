process GENERATE_NEXUS_PATHS{

    conda params.snippy_env 

    tag "${lineage}; t=${distance}"

    publishDir "${params.outDir}/results/networks/${lineage}/", 
        mode: 'copy', 
        overwrite: true

    errorStrategy 'ignore'

    input:
        tuple val(lineage),
            val(distance),
            path(snp_fasta),
            path(snp_tab),
            path(clusters_tab)

    output:
        path("nexus.tuple.csv"), emit: nexus_tuple
        path(clusters_tab),      emit: clusters_tab

    script:

    """
    ###############################################################
    # Part 1: Create samplesheet
    ###############################################################
    lin=\$(echo "${lineage}" | sed 's/lineage/L/g')

    # Identify the unique clusters
        grep -w "${lineage}" ${clusters_tab} \\
            | cut -f4 \\
            | sort \\
            | uniq > unique.clusters.list

    # Remove clusters that are smaller than 3 genomes
        while read clusterID; do
            count=\$(grep -c "\${clusterID}" "${clusters_tab}")
            if [ "\${count}" -ge 3 ]; then
                echo "\${clusterID}" >> final_clusters.list
            fi
        done < unique.clusters.list

    ###############################################################
    # Part 2: Parse over the samplesheet and generate nexus files
    ###############################################################
    # Create output directories
    mkdir -p nexus/ fasta/ positions/

    for clusterID in $(cat unique.clusters.list); do
        
        # Clean up old results
        rm -rf ${params.outDir}/results/${params.runID}/networks/*/\${clusterID}*

                # Get the list of genomes in the cluster
        grep -w \${clusterID} ${clusters_tab} \\
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
            \${clusterID}.fasta

        # Check the correct number of genomes were extracted
        if [ "\$(grep -c '>' \${clusterID}.fasta)" -ne "\$(wc -l < genomes.list)" ]; then
            echo "ERROR: fasta sequences extracted do not match expected number" >&2
            echo "Extracted: \$(grep -c '>' \${clusterID}.fasta)" >&2
            echo "Expected: \$(wc -l < genomes.list)" >&2
            exit 1
        fi
        
        # Try to run snp-sites and check for SNPs
        if snp-sites \${clusterID}.fasta > \${clusterID}.snpsites.fasta 2>/dev/null && [[ -s \${clusterID}.snpsites.fasta ]]; then
            echo "Processing \${clusterID} with SNPs detected"
            
            # Normal processing with SNPs
            snp-sites \${clusterID}.fasta -v \\
                | cut -f2 \\
                | sed '1,4d' \\
                > positions/\${clusterID}_positions.tab

            # H37Rv reference
            awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$3}' \\
                positions/\${clusterID}_positions.tab \\
                ${snp_tab} > \${clusterID}_tmp_refseq

            paste -s -d "" \${clusterID}_tmp_refseq \\
                | sed '1i >H37Rv' \\
                > \${clusterID}_H37Rv.fasta

            # Genomic positions
            awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$1}' \\
                positions/\${clusterID}_positions.tab "${snp_tab}" \\
                > positions/\${clusterID}_genomic_positions.tab

            # MTBC ancestor processing
            cp ${params.mtbc_ancestor_path} tmp.MTB_anc.pos.gz
            gunzip tmp.MTB_anc.pos.gz

            awk 'NR==FNR {pos[\$1]; next} FNR in pos {print \$3}' \\
                positions/\${clusterID}_genomic_positions.tab tmp.MTB_anc.pos \\
                > \${clusterID}_tmp_MTB_anc

            paste -s -d '' \${clusterID}_tmp_MTB_anc \\
                | sed '1i >MTB_anc' \\
                > \${clusterID}_MTB_anc.fasta

            # Combine FASTA files
            cat \${clusterID}.snpsites.fasta \\
                \${clusterID}_H37Rv.fasta \\
                \${clusterID}_MTB_anc.fasta \\
                > fasta/\${clusterID}_refseq.fasta

            # Convert to NEXUS
            seqret -osformat2 nexus \\
                -sequence fasta/\${clusterID}_refseq.fasta \\
                -outseq nexus/\${clusterID}_refseq.nex

        else
            echo "No SNPs detected for \${clusterID}, creating minimal output files"
            
            # Create empty/minimal files
            cp \${clusterID}.fasta fasta/\${clusterID}_refseq.fasta
            touch \${clusterID}.snpsites.fasta
            touch positions/\${clusterID}_positions.tab
            touch positions/\${clusterID}_genomic_positions.tab
            touch nexus/\${clusterID}_refseq.nex
        fi

        # Clean up
        rm -rf tmp.*

        # Simplify the name of the file
        #cat ${snp_tab} > \${clusterID}.snp.tab

        # Final checks (your existing code)
        # check if the nexus generation was successful
        if [[ ! -f positions/\${clusterID}_genomic_positions.tab ]]; then 
            echo "No SNPs detected for \${clusterID}" > positions/\${clusterID}_genomic_positions.tab
        fi

        if [[ ! -f positions/\${clusterID}_genomic_positions.tab ]]; then 
            echo "No SNPs detected for \${clusterID}, nexus generation failed." > nexus/\${clusterID}_refseq.nex
        fi

    done
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 2.0.0
@description:
    This process prepares the paths for the NEXUS files for each cluster.
    It generates a CSV file with the paths to the NEXUS files and the
    corresponding tab files. The CSV file is used as input for the
    GENERATE_NEXUS process.
    The tuple has the following format:
    ["lineage", "clusterID", "fasta_path", "tab_path"]
@changelog:
    v1.0.0-2025-04-01: Initial version
    v1.0.1-2025-04-04: Added filtering to remove clusters with less than 3 genomes
    v1.0.2-2025-12-01: Updated the paths for the new system.
    v2.0.0-2026-05-13: Merged two modules to have a single module that parses over all the 
                        clusterIDs and genertes nexus files

*/