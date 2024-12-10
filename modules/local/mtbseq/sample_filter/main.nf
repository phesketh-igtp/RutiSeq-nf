process MTBSEQ_SAMPLE_FILTER {

    tag "${runID}"

    publishDir "${params.outdir}/bbdd/mtbseq/", mode: 'copy'

    input:
        val min_cov
        val runID
        path tbprof_tbdb_res

    output:
        path "Mapping_and_Variant_Statistics.tab", emit: mtbseq_stats
        path "Strain_Classification.tab", emit: mtbseq_class
        path "mtbseq.paths.minCov.paths"
        tuple val(sampleID), path(called_dir), path(positions_dir), emit: mtbseq_join_paths

    script:

        """

        # Create paths for the MTBseq files needed for pairwise analysis

        ls ${params.outdir}/bbdd/mtbseq/samples > samples.list
        sed 's@^@${params.outdir}/bbdd/mtbseq/samples/@g' samples.list | sed 's@\$@/Called/@g' > samples.called
        sed 's@^@${params.outdir}/bbdd/mtbseq/samples/@g' samples.list | sed 's@\$@/Position_Tables/@g' > samples.pos
        paste -d ',' samples.list samples.called samples.pos > mtbseq.paths.txt
        rm samples.called samples.pos samples.list

        # Concatenate the statistics files and classification files by MTBseq

        echo "Date	SampleID	LibraryID	FullID	Total Reads	Mapped Reads	% Mapped Reads	Genome Size	Genome GC	(Any) Total Bases	% (Any) Total Bases	(Any) GC-Content	(Any) Coverage mean	(Any) Coverage median	(Unambiguous) Total Bases	% (Unambiguous) Total Bases	(Unambiguous) GC-Content	(Unambiguous) Coverage mean	(Unambiguous) Coverage median	SNPs	Deletions	Insertions	Uncovered	Substitutions (Including Stop Codons)" > Mapping_and_Variant_Statistics.tab
        for f in ${params.outdir}/bbdd/mtbseq/samples/*/Statistics/Mapping_and_Variant_Statistics.tab; do
            sed '1d' \$f | sed "s@'@@g" >> Mapping_and_Variant_Statistics.tab
        done

        echo "Date	SampleID	LibraryID	FullID	Homolka species	Homolka lineage	Homolka group	Quality	Coll lineage (branch)	Coll lineage_name (branch)	Coll quality (branch)	Coll lineage (easy)	Coll lineage_name (easy)	Coll quality (easy)	Beijing lineage (easy)	Beijing quality (easy)" > Strain_Classification.tab
        for f in ${params.outdir}/bbdd/mtbseq/samples/*/Classification/Strain_Classification.tab; do
            sed '1d' \$f | sed "s@'@@g" >> Strain_Classification.tab
        done


        # Filter list of samples to only those with >= minimum coverage

        awk -v min_cov="${min_cov}" '(NR>1 && \$19 >= min_cov) {print \$4}' Mapping_and_Variant_Statistics.tab > sample.minCov.list

        grep -f sample.minCov.list mtbseq.paths.txt > mtbseq.paths.minCov.paths


        # Create the tuple for mtbseq_join_paths
        
        while IFS=',' read -r sampleID called_dir positions_dir; do
            echo "\$sampleID,\$called_dir,\$positions_dir"
        done < mtbseq.paths.minCov.paths > mtbseq_join_paths.txt
        """

}