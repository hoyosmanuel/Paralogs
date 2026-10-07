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

Check
~~~~~~~

.. code-block:: bash

  column -t ANALYSES/FBXO3/intermediate/FBXO3_same_scaffold_candidates_raw.tsv

.. code-block::

  Species                       Scaffold            Paralog_start  Paralog_end  Paralog_strand  Paralog_blocks  Paralog_projection             Canonical_start  Canonical_end  Canonical_strand  Canonical_projection
  Furipterus_horrens            manual_scaffold_2   21367653       21374263     +               2               ENST00000265651.8#FBXO3#1544   21311325         21341748       +                 ENST00000265651.8#FBXO3#12
  Rhinolophus_affinis           manual_scaffold_10  24007134       24035799     +               4               ENST00000265651.8#FBXO3#14959  24028650         24062946       +                 ENST00000265651.8#FBXO3#10
  Rhinolophus_ferrumequinum     scaffold_m29_p_11   66435470       66442140     -               3               ENST00000265651.8#FBXO3#3193   66387236         66420715       -                 ENST00000265651.8#FBXO3#10
  Rhinolophus_foetidus          manual_scaffold_3   23142453       23170541     +               3               ENST00000265651.8#FBXO3#87521  23163248         23203534       +                 ENST00000265651.8#FBXO3#10
  Rhinolophus_hipposideros      OZ077427            23133777       23140949     +               3               ENST00000265651.8#FBXO3#4854   23157144         23189003       +                 ENST00000265651.8#FBXO3#11
  Rhinolophus_pearsonii         LG06                23301452       23308437     +               3               ENST00000265651.8#FBXO3#45504  23324429         23358662       +                 ENST00000265651.8#FBXO3#10
  Rhinolophus_pearsonii         LG06                23301452       23308437     +               3               ENST00000265651.8#FBXO3#62166  23324429         23358662       +                 ENST00000265651.8#FBXO3#10
  Rhinolophus_pearsonii         LG06                23305776       23308437     +               2               ENST00000265651.8#FBXO3#66679  23324429         23358662       +                 ENST00000265651.8#FBXO3#10
  Rhinolophus_perniger_lanosus  manual_scaffold_11  23126157       23130990     +               3               ENST00000265651.8#FBXO3#18956  23147124         23180555       +                 ENST00000265651.8#FBXO3#10
  Rhinolophus_sedulus           manual_scaffold_10  23488810       23496589     +               3               ENST00000265651.8#FBXO3#6681   23513107         23548693       +                 ENST00000265651.8#FBXO3#10
  Triaenops_persicus            manual_scaffold_3   22901609       22933074     +               4               ENST00000265651.8#FBXO3#12820  22924851         22956878       +                 ENST00000265651.8#FBXO3#10


Missing Rhinolophus sinicus gene while present in clade. Need to verify.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sabemos que el parálogo tiene estas dos secuencias:
seq1="TTGTTGCTATGTCAATCGAAGACTTTGCCAGCTCTCAAGTCATGATCCGCTATGGAGAAGACATTGCAAAAGATACTGGCTGATATCTGA"
seq2="GGAAGACAAAGCACAGAAGAATCAGTGCTGGAAATCTCTCTTCATTGATACTTACTCTGATGTAGGAAGATATATTGACCATTATGCTGCTGTTAAAAAGGCCTGGGCCGATCTCAAGAAATATTTGGAGCCCAGACCTCCTCAGATGACTTTGTCTCTGCAAG"

Así que podemos buscar la ocurrencia de mas de una vez de esas dos secuencias en el el mismo Scaffold
Entonces hay que crear un fasta sonda de prueba:

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  >FBXO3_block_A
  TTGTTGCTATGTCAATCGAAGACTTTGCCAGCTCTCAAGTCATGATCCGCTATGGAGAAGACATTGCAAAAGATACTGGCTGATATCTGA
  >FBXO3_block_B
  GGAAGACAAAGCACAGAAGAATCAGTGCTGGAAATCTCTCTTCATTGATACTTACTCTGATGTAGGAAGATATATTGACCATTATGCTGCTGTTAAAAAGGCCTGGGCCGATCTCAAGAAATATTTGGAGCCCAGACCTCCTCAGATGACTTTGTCTCTGCAAG
  EOF
