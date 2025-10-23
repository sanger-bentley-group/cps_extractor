process ASSEMBLY_UNICYCLER {
    label 'unicycler_container'
    label 'farm_high_fallible'

    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads), val(sero)
    val min_contig_length

    output:
    tuple val(sample_id), path(fasta), val(sero), path(reads), emit: assembly_ch

    script:
    read1="${reads[0]}"
    read2="${reads[1]}"
    fasta="${sample_id}.contigs.fasta"
    """
    unicycler -1 "$read1" -2 "$read2"  -t "`nproc`" -o results --mode bold
    mv results/assembly.fasta "${fasta}"
    """
}