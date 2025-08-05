#!/bin/bash
set -euo pipefail

SRA=$1
READ=$2

# 📁 Load paths and helper functions
source "$(dirname "$0")/../../config/config.sh"
source "$(dirname "$0")/../../config/utils.sh"

FILE="${SRA}_${READ}.fastq.gz"
DEST="$RAW_DIR/$FILE"

# 🔎 Detect ENA block (needed to form correct FTP URL)
BLOCK_FOUND=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRA&result=read_run&fields=fastq_ftp&format=tsv" |
              grep -v fastq_ftp | head -n1 | awk -F '/' '{print $(NF-2)}')

if [[ -z "$BLOCK_FOUND" ]]; then
  echo "🚫 ENA files not found for: $SRA"
  exit 1
fi

URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/${SRA:0:6}/${BLOCK_FOUND}/$SRA/$FILE"

# ✅ Skip if file already exists and is valid
if [[ -f "$DEST" ]] && gunzip -t "$DEST" &>/dev/null; then
  echo "✅ Already verified: $FILE"
  exit 0
fi

# 🔁 Retry download and verify up to 3 times
ATTEMPT=1
MAX_TRIES=3
while [[ $ATTEMPT -le $MAX_TRIES ]]; do
  echo "🌐 Attempt $ATTEMPT: Downloading $FILE..."
  if wget -q --show-progress --progress=bar:force -c "$URL" -O "$DEST"; then
    echo "📦 Downloaded: $FILE"
    echo "🔍 Verifying..."
    if gunzip -t "$DEST" &>/dev/null; then
      echo "✅ Verified: $FILE"
      exit 0
    else
      echo "❌ Corrupted: $FILE — retrying..."
      rm -f "$DEST"
    fi
  else
    echo "⚠️  wget failed: $FILE"
  fi
  ((ATTEMPT++))
  sleep 2
done

# ❌ Give up after 3 failed attempts
echo "🚫 Failed after $MAX_TRIES attempts: $FILE"
exit 1

#!/bin/bash
set -euo pipefail

SRA=$1
READ=$2

# 📁 Load paths and helper functions
source "$(dirname "$0")/../../config/config.sh"
source "$(dirname "$0")/../../config/utils.sh"

FILE="${SRA}_${READ}.fastq.gz"
DEST="$RAW_DIR/$FILE"

# 🔎 Detect ENA block (needed to form correct FTP URL)
BLOCK_FOUND=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRA&result=read_run&fields=fastq_ftp&format=tsv" |
              grep -v fastq_ftp | head -n1 | awk -F '/' '{print $(NF-2)}')

if [[ -z "$BLOCK_FOUND" ]]; then
  echo "🚫 ENA files not found for: $SRA"
  exit 1
fi

URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/${SRA:0:6}/${BLOCK_FOUND}/$SRA/$FILE"

# ✅ Skip if file already exists and is valid
if [[ -f "$DEST" ]] && gunzip -t "$DEST" &>/dev/null; then
  echo "✅ Already verified: $FILE"
  exit 0
fi

# 🔁 Retry download and verify up to 3 times
ATTEMPT=1
MAX_TRIES=3
while [[ $ATTEMPT -le $MAX_TRIES ]]; do
  echo "🌐 Attempt $ATTEMPT: Downloading $FILE..."
  if wget -q --show-progress --progress=bar:force -c "$URL" -O "$DEST"; then
    echo "📦 Downloaded: $FILE"
    echo "🔍 Verifying..."
    if gunzip -t "$DEST" &>/dev/null; then
      echo "✅ Verified: $FILE"
      exit 0
    else
      echo "❌ Corrupted: $FILE — retrying..."
      rm -f "$DEST"
    fi
  else
    echo "⚠️  wget failed: $FILE"
  fi
  ((ATTEMPT++))
  sleep 2
done

# ❌ Give up after 3 failed attempts
echo "🚫 Failed after $MAX_TRIES attempts: $FILE"
exit 1
