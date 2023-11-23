#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define input parameters
params.matrix = '/home/luigui/Documents/nextflow_rna_seq/rawdata/salmon.gene.counts.matrix'
params.metadata = '/home/luigui/Documents/nextflow_rna_seq/rawdata/metadata.tsv'
params.outdir = '/home/luigui/Documents/nextflow_rna_seq/results/de_plots'

log.info """\
        R N A S E Q  P I P E L I N E    
        ===================================
        count matrix    : ${params.matrix}
        metadata    :   ${params.metadata}
        outdir  : ${params.outdir}
        """
        .stripIndent()

workflow {
    salmonImport()
}

process salmonImport {

    """
    R --no-save &&
    library("DESeq2")
    
    """


}