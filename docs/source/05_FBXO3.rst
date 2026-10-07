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

Execute
~~~~~~~

.. code-block:: bash

  chmod +x ANALYSES/FBXO3/scripts/01_find_FBXO3_hits.sh
  ./ANALYSES/FBXO3/scripts/01_find_FBXO3_hits.sh


2) Automatically identify the canonical copy and the additional non-4-block copy of the same scaffold.
--------------------------------------------------------------------------------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

.. code-block:: bash

  cat > ANALYSES/FBXO3/scripts/02_identify_FBXO3_paralogs.py <<'EOF'
  #!/usr/bin/env python3
  
  # ============================================================
  # 02_identify_FBXO3_paralogs.py
  #
  # Objective:
  # Identify candidate FBXO3 paralog projections from
  # FBXO3_all_hits.tsv.
  #
  # Operational criterion:
  #
  #   1. Canonical FBXO3 is the unique TOGA projection
  #      with 11 blocks.
  #
  #   2. A candidate paralog is an additional FBXO3
  #      projection with a block count different from 11.
  #
  #   3. The additional projection must occur on the same
  #      scaffold as the canonical FBXO3 copy.
  #
  # NOTE:
  # Multiple non-canonical TOGA projections may correspond
  # to the same genomic locus. They are NOT collapsed here.
  # That will be checked in the next step.
  #
  # Input:
  #   ANALYSES/FBXO3/intermediate/FBXO3_all_hits.tsv
  #
  # Output:
  #   ANALYSES/FBXO3/intermediate/
  #       FBXO3_same_scaffold_candidates_raw.tsv
  #
  # ============================================================
  
  import csv
  from collections import defaultdict
  from pathlib import Path
  
  
  ROOT = Path(
      "/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  )
  
  INPUT = (
      ROOT
      / "ANALYSES"
      / "FBXO3"
      / "intermediate"
      / "FBXO3_all_hits.tsv"
  )
  
  OUTPUT = (
      ROOT
      / "ANALYSES"
      / "FBXO3"
      / "intermediate"
      / "FBXO3_same_scaffold_candidates_raw.tsv"
  )
  
  
  # ------------------------------------------------------------
  # Read FBXO3 projections and group by species
  # ------------------------------------------------------------
  
  hits_by_species = defaultdict(list)
  
  with INPUT.open() as infile:
  
      reader = csv.DictReader(
          infile,
          delimiter="\t"
      )
  
      for row in reader:
  
          row["Start"] = int(
              row["Start"]
          )
  
          row["End"] = int(
              row["End"]
          )
  
          row["Blocks"] = int(
              row["Blocks"]
          )
  
          hits_by_species[
              row["Species"]
          ].append(row)
  
  
  # ------------------------------------------------------------
  # Identify canonical FBXO3 and candidate paralog projections
  # ------------------------------------------------------------
  
  candidate_pairs = []
  
  
  for species, hits in sorted(
      hits_by_species.items()
  ):
  
      # Canonical FBXO3:
      # unique projection with exactly 11 blocks
  
      canonical_hits = [
          hit
          for hit in hits
          if hit["Blocks"] == 11
      ]
  
  
      if len(canonical_hits) != 1:
  
          print(
              f"WARNING: {species} has "
              f"{len(canonical_hits)} "
              f"11-block FBXO3 projections"
          )
  
          continue
  
  
      canonical = canonical_hits[0]
  
  
      # Examine all additional FBXO3 projections
  
      for hit in hits:
  
          if hit is canonical:
              continue
  
  
          if (
              hit["Blocks"] != 11
              and
              hit["Scaffold"]
              == canonical["Scaffold"]
          ):
  
              candidate_pairs.append(
                  {
                      "Species":
                          species,
  
                      "Scaffold":
                          hit["Scaffold"],
  
                      "Paralog_start":
                          hit["Start"],
  
                      "Paralog_end":
                          hit["End"],
  
                      "Paralog_strand":
                          hit["Strand"],
  
                      "Paralog_blocks":
                          hit["Blocks"],
  
                      "Paralog_projection":
                          hit["Projection"],
  
                      "Canonical_start":
                          canonical["Start"],
  
                      "Canonical_end":
                          canonical["End"],
  
                      "Canonical_strand":
                          canonical["Strand"],
  
                      "Canonical_projection":
                          canonical["Projection"],
                  }
              )
  
  
  # ------------------------------------------------------------
  # Write output
  # ------------------------------------------------------------
  
  fields = [
      "Species",
      "Scaffold",
      "Paralog_start",
      "Paralog_end",
      "Paralog_strand",
      "Paralog_blocks",
      "Paralog_projection",
      "Canonical_start",
      "Canonical_end",
      "Canonical_strand",
      "Canonical_projection",
  ]
  
  
  with OUTPUT.open(
      "w",
      newline=""
  ) as outfile:
  
      writer = csv.DictWriter(
          outfile,
          fieldnames=fields,
          delimiter="\t"
      )
  
      writer.writeheader()
  
      writer.writerows(
          candidate_pairs
      )
  
  
  # ------------------------------------------------------------
  # Summary
  # ------------------------------------------------------------
  
  species_with_candidates = {
      row["Species"]
      for row in candidate_pairs
  }
  
  
  print()
  print(f"Created: {OUTPUT}")
  
  print(
      f"Same-scaffold non-canonical "
      f"FBXO3 projections: "
      f"{len(candidate_pairs)}"
  )
  
  print(
      f"Species with candidate "
      f"FBXO3 paralogs: "
      f"{len(species_with_candidates)}"
  )
  EOF

Execute
~~~~~~~

.. code-block:: bash

  chmod +x ANALYSES/FBXO3/scripts/02_identify_FBXO3_paralogs.py
  python3 ANALYSES/FBXO3/scripts/02_identify_FBXO3_paralogs.py
