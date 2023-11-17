#!/bin/bash

set -euo pipefail

FULLMAPPING=$1
CHR=$2
START=$3
END=$4
OUTFILE=$5

awk -v chr="$CHR" -v start="$START" -v end="$END" 'BEGIN {FS=OFS=" "} ($1 == chr) && ($2 > start) && ($3 < end)' "$FULLMAPPING" > "$OUTFILE"

