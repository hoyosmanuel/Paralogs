I. IDENTIFICAR LOS PARALOGOS DE TMPO
=====================================

1) Script to extract TMPOs from TOGA
------------------------------------

First, retrieve all TOGA projections annotated as TMPO across all
species. At this stage, no assumption is made about which projection
corresponds to the canonical copy or to a candidate paralog.

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

  mkdir -p \
  ANALYSES/TMPO/scripts \
  ANALYSES/TMPO/intermediate \
  ANALYSES/TMPO/results \
  ANALYSES/TMPO/sequences

  cat > ANALYSES/TMPO/scripts/01_find_TMPO_hits.sh <<'EOF'
  #!/bin/bash
  set -euo pipefail

  ROOT="/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  SPECIES_DIR="$ROOT/SPECIES"
  OUT="$ROOT/ANALYSES/TMPO/intermediate/TMPO_all_hits.tsv"

  # TOGA projection names have the general structure:
  #
  # ENSTXXXXXXXXXXX.version#GENE#projection
  #
  # Here we initially retrieve every projection annotated as TMPO,
  # without assuming a particular human transcript.

  TARGET='#TMPO#'

  printf "Species\tScaffold\tStart\tEnd\tStrand\tBlocks\tProjection\n" > "$OUT"

  for dir in "$SPECIES_DIR"/*; do

      [[ -d "$dir" ]] || continue

      species=$(basename "$dir")

      bed="$dir/query_annotation.bed"
      bedgz="$dir/query_annotation.bed.gz"

      if [[ -f "$bed" ]]; then

          awk -v sp="$species" -v target="$TARGET" '
          BEGIN { OFS="\t" }

          index($4, target) {
              print sp,$1,$2,$3,$6,$10,$4
          }
          ' "$bed" >> "$OUT"

      elif [[ -f "$bedgz" ]]; then

          gzip -cd "$bedgz" | awk -v sp="$species" -v target="$TARGET" '
          BEGIN { OFS="\t" }

          index($4, target) {
              print sp,$1,$2,$3,$6,$10,$4
          }
          ' >> "$OUT"

      else

          echo "WARNING: no query_annotation.bed for $species" >&2

      fi

  done

  echo
  echo "Created:"
  echo "$OUT"

  echo
  echo "TMPO projections:"
  awk 'NR > 1' "$OUT" | wc -l

  echo
  echo "Species represented:"
  awk 'NR > 1 {print $1}' "$OUT" | sort -u | wc -l

  echo
  echo "Human transcripts represented in TOGA projections:"
  awk -F '\t' 'NR > 1 {split($7,a,"#"); print a[1]}' "$OUT" \
  | sort \
  | uniq -c \
  | sort -nr

  echo
  echo "Block-count distribution:"
  awk -F '\t' 'NR > 1 {print $6}' "$OUT" \
  | sort -n \
  | uniq -c

  EOF
