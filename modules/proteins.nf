// Extract protein sequences for each gene
process EXTRACT_PROTEIN_SEQUENCES {

    label 'gffread_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bakta_results), path(cps_sequence), path(annotation_file), path(gb_file), val(reference)

    output:
    tuple val(sample_id), path(annotation_file), val(reference), emit: results_ch
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

    tag "$sample_id"

    input:
    tuple val(sample_id), path(proteins_file)

    output:
    path("*.fa"), emit: proteins_ch

    script:
    """
    create_amino_acid_files.sh ${proteins_file} ${sample_id}
    """
}

// Extract protein sequences for each gene
process CONCAT_PROTEIN_SEQUENCES {
    publishDir "${params.output}/proteins", mode: 'copy', overwrite: true, pattern: "*.fasta"

    label 'bash_container'
    label 'farm_low'

    tag "$protein_name"

    input:
    tuple val(protein_name), path(files)

    output:
    path("*.fasta")

    script:

    """
    # reheader the protein files for the concatenated file to include the sample (rather than have the annotation as the header)
    for f in *.fa
    do
      separator="-${protein_name}.fa"
      echo \${separator}
      sample=\$(echo \${f} | awk -v FS="\${separator}" '{ print \$1 }')
      sed -i "1s/.*/>\${sample}/" \${f}
    done

    cat *.fa > ${protein_name}s.fasta
    """
}