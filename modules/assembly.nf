process ASSEMBLY_SHOVILL {
    label 'shovill_container'
    label 'farm_high_fallible'

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
    shovill --R1 "$read1" --R2 "$read2" --outdir results --cpus "`nproc`" --minlen "$min_contig_length" --force
    mv results/contigs.fa "${fasta}"
    """
}