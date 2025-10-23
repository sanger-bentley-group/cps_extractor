# Run SeroBA to serotype samples

{
    seroba runSerotyping "${SEROBA_DB}" "$READ1" "$READ2" "$SAMPLE_ID" && SEROTYPE=$(tail -1 "${SAMPLE_ID}/pred.csv" | awk -F ',' '{ print $3 }' | awk -F "(" '{ print $NF }' | sed 's|)||g')
} || {
    SEROTYPE="NA"
}

echo Serotype > "$SEROTYPE_REPORT"
echo $SEROTYPE >> "$SEROTYPE_REPORT"