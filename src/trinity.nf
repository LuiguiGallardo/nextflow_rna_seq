params.str = 'your name'
params.samplesList = '/home/luigui/Documents/nextflow_rna_seq/samplesFile.tsv'
params.trinityOutput = '/home/luigui/Documents/nextflow_rna_seq/results/trinity_out_dir_2/'


process trinityAssembly {

    output:
        path './results/trinity_out_dir/'

    """
    docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity --seqType fq --samples_file `pwd`/${params.samplesList} --max_memory 20G --CPU 16 --output `pwd`/${params.trinityOutput}
    """
}

workflow {
    trinityAssembly | view { it.trim() }
}