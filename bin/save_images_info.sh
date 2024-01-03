# Extract containers information from nextflow.config and save into a JSON file

find_image () {
    grep -E "container\s?=" -B 1 "$NEXTFLOW_CONFIG" | grep -v -- "^--$" | paste - - | sort -u | grep "$1" | sed -r "s/.+container\s?=\s?'(.+)'/\1/"
}

BAKTA=$(find_image bakta)
BASH=$(find_image bash)
BLAST=$(find_image blast)
BWA=$(find_image bwa)
CPS_EXTRACTOR_PYTHON=$(find_image cps_extractor_python)
SAMTOOLS=$(find_image samtools)
SHOVILL=$(find_image shovill)

add_container () {
    jq -n --arg container "$1" '.container = $container'
}

jq -n \
    --argjson bakta "$(add_container "$BAKTA")" \
    --argjson bash "$(add_container "$BASH")" \
    --argjson blast "$(add_container "$BLAST")" \
    --argjson bwa "$(add_container "$BWA")" \
    --argjson cps_extractor_python "$(add_container "$CPS_EXTRACTOR_PYTHON")" \
    --argjson samtools "$(add_container "$SAMTOOLS")" \
    --argjson shovill "$(add_container "$SHOVILL")" \
    '$ARGS.named' > "$JSON_FILE"
