// Declare syntax version
nextflow.enable.dsl=2

process CNVkit{

    container = "${projectDir}/../singularity-images/depot.galaxyproject.org-singularity-cnvkit-0.9.9--pyhdfd78af_0.img"

    input:
    path tumor_bam
    path normal_bam
    path fasta
    path target
    path annotate

    output:
    path "*.bed"
    path "*.cnn"
    path "*.cnr"
    path "*.cns"
    path "*.pdf"
    path "*.png"

    script:
    """
    cnvkit.py batch \\
        $tumor_bam \\
        --normal $normal_bam \\
        --fasta $fasta \\
        --targets $target \\
        --annotate $annotate \\
        --processes ${params.threads_num} \\
        --method hybrid --diagram --scatter
    cp *.bed ${launchDir}/${params.outdir}/
    cp *.cnn ${launchDir}/${params.outdir}/
    cp *.cnr ${launchDir}/${params.outdir}/
    cp *.cns ${launchDir}/${params.outdir}/
    cp *.pdf ${launchDir}/${params.outdir}/
    cp *.png ${launchDir}/${params.outdir}/
    """
}

workflow{
    tumor_bam  = Channel.of(params.tumor_bam)
    normal_bam = Channel.of(params.normal_bam)
    fasta      = Channel.of(params.ref_fa)
    target     = Channel.of(params.bed_file)
    annotate   = Channel.of(params.anno_file)
    CNVkit(tumor_bam,normal_bam,fasta,target,annotate)
}
