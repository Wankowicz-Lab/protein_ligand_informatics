# Written almost completely by microsoft copilot
#!/bin/bash

INPUT_FILE="pdb_ids.txt"
OUTPUT_DIR="../PDB_Data_cdk2/cdk2_expanded"
mkdir -p "$OUTPUT_DIR"

BASE_URL="https://pdb-redo.eu/db"

# Define unit cell parameter ranges
MIN_A=52
MAX_A=54
MIN_B=70
MAX_B=72
MIN_C=71
MAX_C=73

# Function to check if a float is within range
within_range() {
    local val=$1 min=$2 max=$3
    awk -v v="$val" -v min="$min" -v max="$max" 'BEGIN { exit !(v >= min && v <= max) }'
}

while read -r pdb_id; do
    [[ -z "$pdb_id" || "$pdb_id" =~ ^# ]] && continue
    pdb_id_lower=$(echo "$pdb_id" | tr '[:upper:]' '[:lower:]')

    # Metadata URL
    info_url="${BASE_URL}/${pdb_id_lower}/data.json"

    # Fetch metadata
    metadata=$(wget -q -O - "$info_url")
    if [[ -z "$metadata" ]]; then
        echo "Failed to fetch metadata for $pdb_id_lower"
        continue
    fi

    # Extract fields using jq
    resolution=$(echo "$metadata" | jq -r '.properties.RESOLUTION')
    spacegroup=$(echo "$metadata" | jq -r '.properties.SPACEGROUP')
    a=$(echo "$metadata" | jq -r '.properties.AAXIS')
    b=$(echo "$metadata" | jq -r '.properties.BAXIS')
    c=$(echo "$metadata" | jq -r '.properties.CAXIS')
    alpha=$(echo "$metadata" | jq -r '.properties.ALPHA')
    beta=$(echo "$metadata" | jq -r '.properties.BETA')
    gamma=$(echo "$metadata" | jq -r '.properties.GAMMA')

    # Apply filters
    if (( $(echo "$resolution > 2.0" | bc -l) )); then
        echo "Skipping $pdb_id_lower: resolution $resolution Å too high"
        continue
    fi

    if [[ "$spacegroup" != "P 21 21 21" ]]; then
        echo "Skipping $pdb_id_lower: space group $spacegroup not allowed"
        continue
    fi

    if (( $(echo "$alpha != 90 || $beta != 90 || $gamma != 90" | bc -l) )); then
        echo "Skipping $pdb_id_lower: angles not 90°"
        continue
    fi

    if ! within_range "$a" "$MIN_A" "$MAX_A" || \
       ! within_range "$b" "$MIN_B" "$MAX_B" || \
       ! within_range "$c" "$MIN_C" "$MAX_C"; then
        echo "Skipping $pdb_id_lower: unit cell dimensions out of range"
        continue
    fi

    # Download structure
    pdb_url="${BASE_URL}/${pdb_id_lower}/${pdb_id_lower}_final.pdb"
    echo "Downloading $pdb_id_lower (resolution: $resolution Å)"
    wget -q -O "${OUTPUT_DIR}/${pdb_id}_final.pdb" "$pdb_url"

done < "$INPUT_FILE"