# Save container information from the workflow into a JSON file

find_image () {
    sort -u "$PROCESSES_CONTAINERS_LIST" | grep "$1" | sed -r "s/.+\t(.+)/\1/"
}

BAKTA=$(find_image bakta)
BASH=$(find_image bash)
BLAST=$(find_image blast)
CHECK_GENE_CONTENT=$(find_image check_gene_content)
CPS_EXTRACTOR_PYTHON=$(find_image cps_extractor_python)
GAP_FILLER=$(find_image gap_filler)
PANAROO=$(find_image panaroo)
SEROBA=$(find_image seroba)
SHOVILL=$(find_image shovill)
SNP_DISTS=$(find_image snp_dists)

add_container () {
    jq -n --arg container "$1" '.container = $container'
}

jq -n \
    --argjson bakta "$(add_container "$BAKTA")" \
    --argjson bash "$(add_container "$BASH")" \
    --argjson blast "$(add_container "$BLAST")" \
    --argjson check_gene_content "$(add_container "$CHECK_GENE_CONTENT")" \
    --argjson cps_extractor_python "$(add_container "$CPS_EXTRACTOR_PYTHON")" \
    --argjson gap_filler "$(add_container "$GAP_FILLER")" \
    --argjson panaroo "$(add_container "$PANAROO")" \
    --argjson seroba "$(add_container "$SEROBA")" \
    --argjson shovill "$(add_container "$SHOVILL")" \
    --argjson snp_dists "$(add_container "$SNP_DISTS")" \
    '$ARGS.named' > "$JSON_FILE"