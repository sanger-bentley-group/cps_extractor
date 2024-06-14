// Extract protein sequences for each gene
process EXTRACT_PROTEIN_SEQUENCES {

    label 'gffread_container'
    label 'farm_low'
    label 'farm_scratchless'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(annotation_file), path(mutation_file), val(reference), path(cps_sequence)

    output:
    tuple val(sample_id), path(annotation_file), path(mutation_file), val(reference), emit: results_ch
    tuple val(sample_id), path(proteins_file), emit: proteins_ch

    script:
    proteins_file="${sample_id}_proteins"
    """
    sed -i "1s/.*/>contig_1/" ${cps_sequence}
    gffread -g ${cps_sequence} -y ${sample_id}_proteins ${annotation_file} -F
    """
}

process CREATE_PROTEIN_FILES {
    publishDir "${params.output}/${sample_id}/proteins", mode: 'copy', overwrite: true, pattern: "*.fa"

    label 'bash_container'
    label 'farm_low'
    label 'farm_scratchless'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(proteins_file)

    output:
    path("*.fa"), emit: proteins_ch

    script:
    """
    create_amino_acid_files.sh ${proteins_file}
    """
}