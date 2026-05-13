process PLOT_MAIN_PHYLOGENY {

    tag "${lineage}"

    conda params.r_phylogeny_env

    publishDir "${params.outDir}/results/${params.runID}/phylogeny/quarto/", 
                mode: 'copy', overwrite: true, pattern: "*.{qmd,Rdata,sh}"
    publishDir "${params.outDir}/results/${params.runID}/phylogeny/", 
                mode: 'copy', overwrite: true, pattern: "*.html"

    input:
        tuple val(lineage), 
            path(contree, stageAs: "snp.contree"), 
            path(alignments, stageAs: "snp.aln.fasta")

        path(processed_clusters, stageAs: "clusters.tsv")
        path(unprocessed_clusters, stageAs: "unprocesses_clusters.tsv")

    output:
    // Main outputs
        path("${lineage}.contree.Rdata"), optional: true
        path("${lineage}_phylogeny.html"), optional: true
    
    // Timetree channel output
        tuple val(lineage),  
            path(alignments),
            path(processed_clusters),
            path("unprocesses_clusters.tsv"), emit: timetree_ch

    script:
        """
        #!/usr/bin/env python
        from __future__ import annotations

        import argparse
        import json
        import os
        from pathlib import Path
        from typing import Dict, Any, List, Tuple

        import polars as pl

        

        """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0
@description:
    Plot the main phylogeny for the analysis, using the ML tree and cluster heatmap if clusters are present.
    If no clusters are present, plot the ML tree without the cluster heatmap.
@changelog:
    2025-04-01: Initial version.
    2025-10-01: Remove static tree PDF and replaced with quarto report generation.
    2025-10-28: Corrected the Quarto report, that was not producing the MonoPhyl output. 
*/
