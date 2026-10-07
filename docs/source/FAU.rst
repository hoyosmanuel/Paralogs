I. Install SLiM
================

1) Script to extract FAUs from TOGA
-----------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  
  cat > ANALYSES/FAU/scripts/01_find_FAU_hits.sh <<'EOF'
  #!/bin/bash
  set -euo pipefail
  
  ROOT="/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  SPECIES_DIR="$ROOT/SPECIES"
  OUT="$ROOT/ANALYSES/FAU/intermediate/FAU_all_hits.tsv"
  
  #### ^ENST00000529639\.   👉🏽 Starts exactly with ENST00000529639.
  #### [0-9]+               👉🏽 Any transcript version number
  #### #FAU#                👉🏽 Followed by #FAU#
  
  TARGET='^ENST00000529639\.[0-9]+#FAU#'
  
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
  echo "FAU projections:"
  awk 'NR > 1' "$OUT" | wc -l
  
  echo
  echo "Species represented:"
  awk 'NR > 1 {print $1}' "$OUT" | sort -u | wc -l
  EOF



2) Execute
-----------

.. code-block:: bash


  chmod +x ANALYSES/FAU/scripts/01_find_FAU_hits.sh
  ./ANALYSES/FAU/scripts/01_find_FAU_hits.sh

You should see something like this:

.. code-block::

  ./ANALYSES/FAU/scripts/01_find_FAU_hits.sh
  Created:
  /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/FAU/intermediate/FAU_all_hits.tsv
  
  FAU projections:
       118
  
  Species represented:
       103
  (base) manuelhoyos@MacBookPro bat_HTF_genomic_analysis %
