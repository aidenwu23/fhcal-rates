#!/usr/bin/env bash
# ./scripts/data_prep/download.sh 0 9

(
  set -euo pipefail

  START_INDEX="${1:-0}"
  END_INDEX="${2:-0}"

  LIST="data/metadata/file_list.txt"

  START_LINE=$((START_INDEX + 1))
  END_LINE=$((END_INDEX + 1))

  mapfile -t FILES < <(sed -n "${START_LINE},${END_LINE}p" "$LIST")

  for FILE in "${FILES[@]}"; do
    echo "Downloading $FILE"
    rucio download "$FILE" --dir data/raw </dev/null
  done
)
