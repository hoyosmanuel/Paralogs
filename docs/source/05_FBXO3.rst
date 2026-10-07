I. IDENTIFICAR LOS PARALOGOS DE FBXO3
=====================================

1) Script to extract FBXO3s from TOGA
--------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

  cat > ANALYSES/FBXO3/scripts/01_find_FBXO3_hits.sh <<'EOF'
  #!/bin/bash
  set -euo pipefail
  
  ROOT="/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  SPECIES_DIR="$ROOT/SPECIES"
  OUT="$ROOT/ANALYSES/FBXO3/intermediate/FBXO3_all_hits.tsv"
  
  #### ^ENST00000265651\.   👉🏽 Starts exactly with ENST00000265651.
  #### [0-9]+               👉🏽 Any transcript version number
  #### #FBXO3#              👉🏽 Followed by #FBXO3#
  
  TARGET='^ENST00000265651\.[0-9]+#FBXO3#'
  
  printf "Species\tScaffold\tStart\tEnd\tStrand\tBlocks\tProjection\n" > "$OUT"
  
  for dir in "$SPECIES_DIR"/*; do
  
      [[ -d "$dir" ]] || continue
  
      species=$(basename "$dir")
  
      bed="$dir/query_annotation.bed"
      bedgz="$dir/query_annotation.bed.gz"
  
      if [[ -f "$bed" ]]; then
  
          awk -v sp="$species" -v target="$TARGET" '
          BEGIN { OFS="\t" }
  
          $4 ~ target {
              print sp,$1,$2,$3,$6,$10,$4
          }
          ' "$bed" >> "$OUT"
  
      elif [[ -f "$bedgz" ]]; then
  
          gzip -cd "$bedgz" | awk -v sp="$species" -v target="$TARGET" '
          BEGIN { OFS="\t" }
  
          $4 ~ target {
              print sp,$1,$2,$3,$6,$10,$4
          }
          ' >> "$OUT"
  
      else
          echo "WARNING: no query_annotation.bed for $species" >&2
      fi
  
  done
  
  echo "Created:"
  echo "$OUT"
  
  echo
  echo "FBXO3 projections:"
  awk 'NR > 1' "$OUT" | wc -l
  
  echo
  echo "Species represented:"
  awk 'NR > 1 {print $1}' "$OUT" | sort -u | wc -l
  EOF
