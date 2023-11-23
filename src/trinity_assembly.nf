#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define input parameters
params.reads = '/home/luigui/Documents/nextflow_rna_seq/rawdata/*_R{1,2}.fastq.gz'
params.samplesFile = '/home/luigui/Documents/nextflow_rna_seq/rawdata/samplesFile.tsv'
params.outdir = '/home/luigui/Documents/nextflow_rna_seq/results/de_analysis'

log.info """\
        R N A S E Q  P I P E L I N E    
        ===================================
        reads        : ${params.reads}
        outdir       : ${params.outdir}
        samplesFile :   ${params.samplesFile}
        """
        .stripIndent()

// Create Trinity 'de novo' assembly
process trinityAssembly {
    publishDir "${params.outdir}/trinityAssembly"

    input:
    path samplesFile

    output:
    path "trinityOutput"

    """
    Trinity --seqType fq --samples_file ${params.samplesFile} --max_memory 20G --CPU 16 --output trinityOutput
    """
}

// Create Salmon index
process salmonIndex {
    publishDir "${params.outdir}/salmonIndex"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome/Trinity.fasta -i index
    """
}

// Quantification with Salmon
process salmonQuant {
    publishDir "${params.outdir}/salmonQuant"

    tag "$pair_id"

    input:
    path index
    tuple val(pair_id), path(reads)
    path transcriptome

    output:
    path pair_id

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process salmonMerge {
    publishDir "${params.outdir}/salmonMerge"

    input:
    path transcriptome

    output:
    path "*"

    """
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map $transcriptome/Trinity.fasta.gene_trans_map  --name_sample_by_basedir --quant_files /home/luigui/Documents/nextflow_rna_seq/rawdata/quant_files.tsv
    """
}

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )

    trinityAssembly(params.samplesFile)
    salmonIndex(trinityAssembly.out)
    salmonQuant(salmonIndex.out, read_pairs_ch, trinityAssembly.out)
    salmonMerge(trinityAssembly.out)
}
