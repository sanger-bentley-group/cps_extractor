# Save received tools versions into a JSON file

add_version () {
    jq -n --arg version "$1" '.version = $version'
}

jq -n \
    --argjson bakta "$(add_version "$BAKTA_VERSION")" \
    --argjson bedtools "$(add_version "$BEDTOOLS_VERSION")" \
    --argjson blast "$(add_version "$BLAST_VERSION")" \
    --argjson bwa "$(add_version "$BWA_VERSION")" \
    --argjson python "$(add_version "$PYTHON_VERSION")" \
    --argjson samtools "$(add_version "$SAMTOOLS_VERSION")" \
    --argjson shovill "$(add_version "$SHOVILL_VERSION")" \
    '$ARGS.named' > "$JSON_FILE"
