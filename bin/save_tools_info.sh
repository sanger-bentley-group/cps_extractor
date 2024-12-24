# Save received tools versions into a JSON file

add_version () {
    jq -n --arg version "$1" '.version = $version'
}

jq -n \
    --argjson ariba "$(add_version "$ARIBA_VERSION")" \
    --argjson bakta "$(add_version "$BAKTA_VERSION")" \
    --argjson bcftools "$(add_version "$BCFTOOLS_VERSION")" \
    --argjson bedtools "$(add_version "$BEDTOOLS_VERSION")" \
    --argjson blast "$(add_version "$BLAST_VERSION")" \
    --argjson bwa "$(add_version "$BWA_VERSION")" \
    --argjson panaroo "$(add_version "$PANAROO_VERSION")" \
    --argjson python "$(add_version "$PYTHON_VERSION")" \
    --argjson samtools "$(add_version "$SAMTOOLS_VERSION")" \
    --argjson seroba "$(add_version "$SEROBA_VERSION")" \
    --argjson shovill "$(add_version "$SHOVILL_VERSION")" \
    --argjson snpdists "$(add_version "$SNPDISTS_VERSION")" \
    '$ARGS.named' > "$JSON_FILE"
