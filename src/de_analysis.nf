#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define input parameters
params.matrix = '/home/luigui/Documents/nextflow_rna_seq/results/de_analysis/salmonMerge/salmon.gene.TMM.EXPR.matrix'
params.metadata = '/home/luigui/Documents/nextflow_rna_seq/rawdata/samplesFile.tsv'
params.outdir = '/home/luigui/Documents/nextflow_rna_seq/results/de_analysis'

log.info """\
        R N A S E Q  P I P E L I N E    
        ===================================
        matrix        : ${params.matrix}
        outdir       : ${params.outdir}
        metadata :   ${params.metadata}
        """
        .stripIndent()


process edgerAnalysis {
    publishDir "${params.outdir}", mode: 'move', overwrite: false
    
    input:
    path matrix
    path metadata

    output:
    path "edgeR"

    """
    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix $matrix --method edgeR --samples_file $metadata --output edgeR
    
    cd edgeR

    $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../$matrix -P 1 -C 2 --samples ../$metadata
    """
}

process deseq2Analysis {
    publishDir "${params.outdir}", mode: 'move', overwrite: false

    input:
    path matrix
    path metadata

    output:
    path deseq2 

    """
    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix $matrix --method DESeq2 --samples_file $metadata --output deseq2
    
    cd deseq2

    $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../$matrix -P 1 -C 2 --samples ../$metadata
    """
}

// process heatmaps {
//     input:
//     path analysis
//     path matrix
//     path metadata

//     output:
//     path "*"


//     """
//     cd $analysis && $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../$matrix -P 1 -C 2 --samples ../$metadata
//     """
// }

workflow {
    edgerAnalysis(params.matrix, params.metadata)
    deseq2Analysis(params.matrix, params.metadata)
    // heatmaps(edgerAnalysis.out, params.matrix, params.metadata)
    // heatmaps(deseq2Analysis.out, params.matrix, params.metadata)

    view { it.trim() }
}