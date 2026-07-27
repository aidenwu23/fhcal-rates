#!/usr/bin/env bash
# ./scripts/data_prep/download.sh --start 0 --end 9 --metadata data/metadata/bkg.txt --output data/bkg

(
  set -euo pipefail

  START_INDEX=0
  END_INDEX=0
  LIST="data/metadata/bkg.txt"
  OUTPUT_DIR=""

  while (($#)); do
    case "$1" in
      --start)
        START_INDEX="${2:?--start requires a value}"
        shift 2
        ;;
      --end)
        END_INDEX="${2:?--end requires a value}"
        shift 2
        ;;
      --metadata)
        LIST="${2:?--metadata requires a value}"
        shift 2
        ;;
      --output)
        OUTPUT_DIR="${2:?--output requires a value}"
        shift 2
        ;;
      *)
        echo "Unknown option: $1" >&2
        exit 2
        ;;
    esac
  done

  START_LINE=$((START_INDEX + 1))
  END_LINE=$((END_INDEX + 1))

  mapfile -t FILES < <(sed -n "${START_LINE},${END_LINE}p" "$LIST")

  for FILE in "${FILES[@]}"; do
    echo "Downloading $FILE"
    if [[ -n "$OUTPUT_DIR" ]]; then
      rucio download "$FILE" --dir "$OUTPUT_DIR" --no-subdir </dev/null

      RELATIVE_PATH="${FILE#*:}"
      NESTED_FILE="$OUTPUT_DIR/${RELATIVE_PATH#/}"
      OUTPUT_FILE="$OUTPUT_DIR/${FILE##*/}"

      if [[ -f "$NESTED_FILE" ]]; then
        mv "$NESTED_FILE" "$OUTPUT_FILE"
        find "$OUTPUT_DIR" -mindepth 1 -type d -empty -delete
      fi
    else
      rucio download "$FILE" --dir data/raw </dev/null
    fi
  done
)
