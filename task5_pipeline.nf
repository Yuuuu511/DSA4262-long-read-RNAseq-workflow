#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Task 5 - Long-read RNA-seq workflow
 *
 * Steps:
 * 1. Minimap2 alignment
 * 2. SAM to BAM conversion with samtools
 * 3. Transcript discovery and quantification with Bambu
 * 4. QC with NanoPlot and samtools stats
 */

params.reads = '/home/ubuntu/task3/fastq/*.fastq.gz'
params.refFa = '/home/ubuntu/task3/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.refGtf = '/home/ubuntu/task3/reference/Homo_sapiens.GRCh38.91.gtf'

params.outdir = 'results'
params.use_annotation = true


/*
 * STEP 1
 * Align reads to the reference genome using Minimap2.
 *
 * directRNA samples use -uf -k14
 * cDNA samples use -k14
 */
process MINIMAP2_ALIGN {

    cpus 4

    input:
    tuple val(sample_id), val(read_type), path(reads)
    path refFa

    output:
    tuple val(sample_id), val(read_type), path(reads),
          path("${sample_id}.sam")

    script:

    def extra_args = read_type == 'directRNA' ? '-uf -k14' : '-k14'

    """
    minimap2 -ax splice ${extra_args} -t ${task.cpus} \
        ${refFa} ${reads} > ${sample_id}.sam
    """
}


/*
 * STEP 2
 * Convert SAM to BAM using samtools.
 */
process SAM_TO_BAM {

    publishDir "${params.outdir}/bam", mode: 'copy'

    cpus 4

    input:
    tuple val(sample_id), val(read_type), path(reads), path(reads_sam)

    output:
    tuple val(sample_id), val(read_type), path(reads),
          path("${sample_id}.bam")

    script:

    """
    samtools view -@ ${task.cpus} -b ${reads_sam} > ${sample_id}.bam
    """
}


/*
 * STEP 3
 * Transcript discovery and quantification using Bambu.
 *
 * Scenario 1:
 *   params.use_annotation = true
 *
 * Scenario 2:
 *   params.use_annotation = false
 */
process BAMBU {

    publishDir "${params.outdir}/bambu", mode: 'copy'

    cpus 4

    input:
    path refFa
    path refGtf
    path reads_bam

    output:
    path "counts_transcript.txt", emit: transcript_counts
    path "counts_gene.txt", emit: gene_counts
    path "extended_annotations.gtf", emit: extended_gtf

    script:

    def annotation_code = params.use_annotation ?
        "annotations <- prepareAnnotations('${refGtf}')\n" +
        "se <- bambu(reads = samples.bam, annotations = annotations, genome = '${refFa}', NDR = 1, ncore = ${task.cpus})" :
        "se <- bambu(reads = samples.bam, genome = '${refFa}', NDR = 1, ncore = ${task.cpus})"

    """
    #!/usr/bin/env Rscript --vanilla

    library(bambu)

    samples.bam <- list.files(
        ".",
        pattern = "\\\\.bam\$",
        full.names = TRUE
    )

    ${annotation_code}

    writeBambuOutput(se, path = "./")
    """
}


/*
 * STEP 4
 * Quality control.
 *
 * NanoPlot:
 *   raw long-read sequencing quality and read-length statistics.
 *
 * samtools stats:
 *   alignment-based statistics from the BAM file.
 *
 * The Bambu GTF is included as an input so that this QC process
 * starts only after Bambu has successfully completed.
 */
process QC {

    publishDir "${params.outdir}/qc", mode: 'copy'

    cpus 2

    input:
    tuple val(sample_id), val(read_type), path(reads), path(reads_bam)
    path bambu_gtf

    output:
    path "${sample_id}_nanoplot"
    path "${sample_id}_samtools_stats.txt"

    script:

    """
    mkdir -p ${sample_id}_nanoplot

    NanoPlot \
        --fastq ${reads} \
        --outdir ${sample_id}_nanoplot \
        --threads ${task.cpus}

    samtools stats ${reads_bam} > ${sample_id}_samtools_stats.txt
    """
}


workflow {

    /*
     * Create one entry for each FASTQ file and determine whether
     * it is direct RNA or cDNA from the filename.
     */
    reads_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { reads ->
            def sample_id = reads.baseName.replace('.fastq', '')
            def read_type = reads.name.contains('directRNA') ? 'directRNA' : 'cDNA'

            tuple(sample_id, read_type, reads)
        }

    ref_fa_ch = Channel.value(file(params.refFa))
    ref_gtf_ch = Channel.value(file(params.refGtf))

    /*
     * Alignment
     */
    MINIMAP2_ALIGN(
        reads_ch,
        ref_fa_ch
    )

    /*
     * SAM -> BAM
     */
    SAM_TO_BAM(
        MINIMAP2_ALIGN.out
    )

    /*
     * Collect all four BAM files so Bambu analyses the full dataset
     * together.
     */
    bam_files_ch = SAM_TO_BAM.out
        .map { sample_id, read_type, reads, bam -> bam }
        .collect()

    /*
     * Bambu
     */
    BAMBU(
        ref_fa_ch,
        ref_gtf_ch,
        bam_files_ch
    )

    /*
     * QC is deliberately placed after Bambu.
     * Converting the Bambu GTF output to a value channel allows
     * every sample to run through QC after Bambu has completed.
     */
    bambu_done_ch = BAMBU.out.extended_gtf

    QC(
        SAM_TO_BAM.out,
        bambu_done_ch
    )
}
