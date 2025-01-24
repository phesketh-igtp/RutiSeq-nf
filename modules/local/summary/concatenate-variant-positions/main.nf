process CONCATENATED_VARIANT_FILES {

    publishDir "${params.outdir}/results/snps/", mode: 'copy'

    input:
        path("tabular_vars_*.csv")
        path("tabular_var_counts_*.csv")

    output:
        path("variant-positions.csv"),          emit: variant_positions
        path("variant-positions.counts.csv")
        path("cleanup-handover"),               emit: cleanup_handover


    script:

        """
        # Concatenate all tabular_vars files, excluding the header
            echo "positions;Sample;allel;cluster;genome_type;ID;name;start;stop;frame;product;description;function;cogcats;status_region;status_function;type;region_number;function_number" > variant-positions.csv

            cat tabular_vars_*.csv | sed '1!{/^positions/d;}' >> variant-positions.csv

        # Concatenate all tabular_var_counts files, excluding the header
            echo "positions;allel;cluster;genome_type;freq;ID;name;start;stop;frame;product;description;function;cogcats;status_region;status_function;type;region_number;function_number" > variant-positions.counts.csv
            
            cat tabular_var_counts_*.csv | sed '1!{/^positions/d;}' >> variant-positions.counts.csv

            touch cleanup-handover
        """

}