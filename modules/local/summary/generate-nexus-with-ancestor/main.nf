process GENERATE_NEXUS_W_MRCA {

    conda params.snp_profiling_env 

    tag "cluster: ${clusterID}"

    array 100

    publishDir "${params.outdir}/results/networks/", mode: 'copy'

    input:
        tuple val(lineage), 
                val(clusterID),
                file(snp_fasta),
                file(snp_tab),
                file(ancestor)
                
        file(pairwise_clusters)

    output:
        tuple val(clusterID),
                path("nexus/${clusterID}_refseq_mrca.nex"),            emit: nexus_w_no_metadata
        path("fasta/*")
        path("positions/*")
        path("nexus/*")

    script:

        """
        # create the output and temporary directories
            mkdir -p nexus/ fasta/ positions/

        # create the list of genomes within the cluster
            grep "${clusterID}" ${pairwise_clusters} | cut -f1 > ${clusterID}.genomes.list

        #·················································································#

        # create cluster directory and split up fasta file in cluster fastas
            while IFS=";" read -r genome; do
                seqkit grep -w 0 -n -p \${genome} ${snp_fasta} >> ${clusterID}.fasta
            done < ${clusterID}.genomes.list

        # run snp-sites on the fastas
            snp-sites ${clusterID}.fasta > ${clusterID}.snpsites.fasta
            snp-sites ${clusterID}.fasta -v | cut -f2 \\
                | sed '1,4d' > positions/${clusterID}_positions.tab


        #·················································································#

        # Create the variant alignment for the NODE ancestor

            #for i in `cat positions/${clusterID}_positions.tab`; do 
            #    sed -n \${i}'p' ${ancestor} | cut -f4
            #done > ${clusterID}_node_anc
            
            awk 'NR==FNR {pos[\$1]; next} FNR in pos {print \$4}' positions/${clusterID}_positions.tab ${ancestor} > ${clusterID}_node_anc


            # convert the column in fasta
                paste -s -d "" ${clusterID}_node_anc \\
                    | sed "1i >MRCA" > ${clusterID}_MRCA.fasta

        #·················································································#

        # H37Rv variance positions 
            #for i in `cat positions/${clusterID}_positions.tab`; do 
            #    sed -n \$((i+2))'p' ${snp_tab} | cut -f3
            #done > ${clusterID}_tmp_refseq

            awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$3}' positions/${clusterID}_positions.tab ${snp_tab} > ${clusterID}_tmp_refseq

                # convert column into fasta
                paste -s -d "" ${clusterID}_tmp_refseq \\
                    | sed '1i >H37Rv' > ${clusterID}_H37Rv.fasta

        #·················································································#

        # Get genomic positions
            #while read -r position; do
            #    sed -n \$((position+2))'p' ${snp_tab} | cut -f 1; 
            #done < positions/${clusterID}_positions.tab > positions/${clusterID}_genomic_positions.tab

            awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$1}' positions/${clusterID}_positions.tab ${snp_tab} > positions/${clusterID}_genomic_positions.tab


        #·················································································#

        # Valencian ancestor (MTB_anc) variance positions
            cp ${params.mtbc_ancestor_path} ${lineage}.tmp.MTB_anc.pos.gz
            gunzip ${lineage}.tmp.MTB_anc.pos.gz

            for i in `cat positions/${clusterID}_genomic_positions.tab`; do 
                sed -n \${i}'p' ${lineage}.tmp.MTB_anc.pos | cut -f3
            done > ${clusterID}_tmp_MTB_anc

            # convert the column in fasta
                paste -s -d "" ${clusterID}_tmp_MTB_anc \\
                    | sed '1i >MTB_anc' > ${clusterID}_MTB_anc.fasta

        # remove the large tab file
            rm -rf ${lineage}.tmp.MTB_anc.pos.gz

        #·················································································#

        # Create final FASTA file
            cat ${clusterID}.snpsites.fasta \\
                ${clusterID}_H37Rv.fasta \\
                ${clusterID}_MTB_anc.fasta \\
                ${clusterID}_MRCA.fasta \\
                > fasta/${clusterID}_refseq_mrca.fasta

        # convert to nexus for visualisation
            seqret -osformat2 nexus -sequence fasta/${clusterID}_refseq_mrca.fasta \\
                -outseq nexus/${clusterID}_refseq_mrca.nex
        """

}