process GENERATE_NEXUS {

    conda params.snp_profiling_env 

    tag "${lineage}"

    publishDir "${params.outdir}/bbdd/results/networks/", mode: 'copy'

    input:
        path pairwise_clusters
        tuple val(lineage), 
                path(joint_dir), 
                path(amend_dir)

    output:
        path("nexus/*"),                emit: nexus_dir
        path("fasta/*")
        path("positions/*")      

    script:
        """

        ln -s ${amend_dir}/${lineage}_*_phylo_w${params.mtbseq_window}.fasta sequences.fasta
        ln -s ${amend_dir}/${lineage}_*_phylo_w${params.mtbseq_window}.tab positions.tab

        cut -f1,7 ${pairwise_clusters} > ${lineage}.clusters.samples.tsv
        cut -f7 ${pairwise_clusters} | sort | uniq > clusters.list

        mkdir -p nexus/; mkdir -p fasta/; mkdir -p positions/
    
        while IFS='\t' read -r genome clusterID; do
        
        # create cluster directory and split up fasta file in cluster fastas
            
            mkdir -p tmp.\${clusterID}/

            while IFS=";" read -r genome; do

            # 1. create cluster directory and split up fasta file in cluster fastas
                seqkit grep -w 0 -n -p \${genome} sequences.fasta >> fasta/tmp.\${clusterID}/\${clusterID}.fasta

            done < clusters.list

            # 2. run snp-sites on the fastas
            snp-sites nexus/tmp.\${clusterID}/\${clusterID}.fasta > nexus/tmp.\${clusterID}/\${clusterID}.snpsites.fas
            snp-sites nexus/tmp.\${clusterID}/\${clusterID}.fasta -v | cut -f2 | sed '1,4d' > nexus/tmp.\${clusterID}/\${clusterID}_positions

            # 3. obtain the reference positions (H37Rv) for the cluster positions
            for i in `cat nexus/tmp.\${clusterID}/\${clusterID}_positions`; do 
                sed -n $((i+2))'p' positions.tab | cut -f3
            done > nexus/tmp.\${clusterID}/\${clusterID}_tmp_refseq
            # convert column into fasta
                paste -s -d "" nexus/tmp.\${clusterID}/\${clusterID}_tmp_refseq | sed '1i >H37Rv' > nexus/tmp.\${clusterID}/\${clusterID}_H37Rv.fas
                
                
            # 4. Get genomic positions
            while read -r position; do
                sed -n $((position+2))'p' positions.tab | cut -f 1; 
            done < nexus/tmp.\${clusterID}/\${clusterID}_positions > nexus/tmp.\${clusterID}/\${clusterID}_genomic_positions

            # 5. obtain the reference positions using Valencian ancestor (MTB_anc) for the cluster positions
            cp ${params.mtbc_ancestor_path} ${lineage}.tmp.MTB_anc.pos.gz; gunzip ${lineage}.tmp.MTB_anc.pos.gz

            for i in `cat nexus/tmp.\${clusterID}/\${clusterID}_genomic_positions`; do 
                sed -n ${i}'p' nexus/tmp.\${clusterID}/MTB_anc.pos | cut -f3 # doesnt need to +2 as the tsv file has no header
            done > nexus/tmp.\${clusterID}/\${clusterID}_tmp_MTB_anc
            # convert the column in fasta
                paste -s -d "" nexus/tmp.\${clusterID}/\${clusterID}_tmp_MTB_anc | sed '1i >MTB_anc' > nexus/tmp.\${clusterID}/\${clusterID}_MTB_anc.fas
                rm nexus/tmp.\${clusterID}/MTB_anc.pos

            # 6. join the fasta files
            cat nexus/tmp.\${clusterID}/\${clusterID}.snpsites.fas \
                    nexus/tmp.\${clusterID}/\${clusterID}_H37Rv.fas \
                    nexus/tmp.\${clusterID}/\${clusterID}_MTB_anc.fas > nexus/tmp.\${clusterID}/\${clusterID}_refseq.fas

            # 7. save as nexus
            seqret -osformat2 nexus \
                -sequence nexus/tmp.\${clusterID}/\${clusterID}_refseq.fas \
                -outseq nexus/tmp.\${clusterID}/\${clusterID}_refseq.nexus

            done < 

            #····································································#
            #····································································#

            #·········· CLEANUP ··········#

            mv nexus/tmp.\${clusterID}/\${clusterID}_refseq.nexus nexus/\${clusterID}_refseq.nex
            mv nexus/tmp.\${clusterID}/\${clusterID}.fasta nexus/fasta/\${clusterID}.fasta
            gzip nexus/fasta/\${clusterID}.fasta
            mv nexus/tmp.\${clusterID}/\${clusterID}_genomic_positions nexus/positions/

            # rm the temporary directory
            rm -rf nexus/tmp* ${lineage}.tmp.MTB_anc.pos.gz

        """

}