#!/bin/bash

# 📁 Load directory configuration
source "$(dirname "$0")/../../config/config.sh"

# 📄 Read SRA IDs into an array
mapfile -t SRA_IDS < "$METADATA_DIR/sra_ids.txt"

# 🚨 Exit if file is empty
if [[ ${#SRA_IDS[@]} -eq 0 ]]; then
  echo "❌ No SRA IDs found in $METADATA_DIR/sra_ids.txt"
  exit 1
fi

# 🛠 Generate paired-end download tasks (READ 1 and 2)
> "$METADATA_DIR/joblist.txt"
for SRA in "${SRA_IDS[@]}"; do
  echo "$SRA 1" >> "$METADATA_DIR/joblist.txt"
  echo "$SRA 2" >> "$METADATA_DIR/joblist.txt"
done

echo "✅ Created joblist at: $METADATA_DIR/joblist.txt"
