I. IDENTIFICAR LOS PARALOGOS DE FAU
===================================

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



Execute
~~~~~~~~

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


2) Automatically identify the canonical 4-block copy and the additional non-4-block copy of the same scaffold.
---------------------------------------------------------------------------------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

.. code-block:: bash

  cat > ANALYSES/FAU/scripts/02_identify_FAU_paralogs.py <<'EOF'
  #!/usr/bin/env python3
  
  # ============================================================
  # 02_identify_FAU_paralogs.py
  #
  # Objetivo:
  # Identificar candidate FAU paralogs a partir de FAU_all_hits.tsv.
  #
  # Criterio operacional:
  #
  #   1. La copia canónica de FAU es la proyección TOGA completa
  #      con 4 bloques.
  #
  #   2. Un candidate paralog es una proyección adicional de FAU
  #      con un número de bloques diferente de 4.
  #
  #   3. La copia adicional debe estar localizada en el mismo
  #      scaffold que la copia canónica.
  #
  # IMPORTANTE:
  # La orientación opuesta entre canonical y paralog NO se utiliza
  # como criterio de selección. Es una observación posterior.
  #
  # Input:
  # ANALYSES/FAU/intermediate/FAU_all_hits.tsv
  #
  # Output:
  # ANALYSES/FAU/intermediate/FAU_same_scaffold_paralogs.tsv
  # ============================================================
  
  import csv
  from collections import defaultdict
  from pathlib import Path
  
  
  # ------------------------------------------------------------
  # Directorio principal del proyecto
  # ------------------------------------------------------------
  ROOT = Path("/Volumes/Expansion/project3/bat_HTF_genomic_analysis")
  
  INPUT = ROOT / "ANALYSES/FAU/intermediate/FAU_all_hits.tsv"
  
  OUTPUT = ROOT / "ANALYSES/FAU/intermediate/FAU_same_scaffold_paralogs.tsv"
  
  
  # ------------------------------------------------------------
  # Leer todas las proyecciones FAU y agruparlas por especie
  # ------------------------------------------------------------
  hits_by_species = defaultdict(list)
  
  with INPUT.open() as infile:
  
      reader = csv.DictReader(infile, delimiter="\t")
  
      for row in reader:
  
          # Convertir las variables numéricas
          row["Start"] = int(row["Start"])
          row["End"] = int(row["End"])
          row["Blocks"] = int(row["Blocks"])
  
          hits_by_species[row["Species"]].append(row)
  
  
  # ------------------------------------------------------------
  # Identificar canonical FAU y candidate paralogs
  # ------------------------------------------------------------
  candidate_pairs = []
  
  for species, hits in sorted(hits_by_species.items()):
  
      # --------------------------------------------------------
      # Canonical FAU:
      # proyección completa con exactamente 4 bloques
      # --------------------------------------------------------
      canonical_hits = [
          hit for hit in hits
          if hit["Blocks"] == 4
      ]
  
  
      # --------------------------------------------------------
      # En este dataset esperamos una única copia canónica
      # de 4 bloques por especie.
      #
      # Si no ocurre, imprimimos una advertencia y evitamos
      # hacer una clasificación automática dudosa.
      # --------------------------------------------------------
      if len(canonical_hits) != 1:
  
          print(
              f"WARNING: {species} has "
              f"{len(canonical_hits)} four-block FAU projections"
          )
  
          continue
  
  
      canonical = canonical_hits[0]
  
  
      # --------------------------------------------------------
      # Examinar todas las otras proyecciones FAU de la especie
      # --------------------------------------------------------
      for hit in hits:
  
          # Ignorar la propia copia canónica
          if hit is canonical:
              continue
  
  
          # ----------------------------------------------------
          # Candidate paralog:
          #
          # - debe ser una copia no-canónica
          # - debe tener != 4 bloques
          # - debe encontrarse en el mismo scaffold
          # ----------------------------------------------------
          if (
              hit["Blocks"] != 4
              and hit["Scaffold"] == canonical["Scaffold"]
          ):
  
              candidate_pairs.append(
                  {
                      "Species": species,
                      "Scaffold": hit["Scaffold"],
  
                      "Paralog_start": hit["Start"],
                      "Paralog_end": hit["End"],
                      "Paralog_strand": hit["Strand"],
                      "Paralog_blocks": hit["Blocks"],
  
                      "Canonical_start": canonical["Start"],
                      "Canonical_end": canonical["End"],
                      "Canonical_strand": canonical["Strand"],
                  }
              )
  
  
  # ------------------------------------------------------------
  # Escribir tabla final
  # ------------------------------------------------------------
  fields = [
      "Species",
      "Scaffold",
      "Paralog_start",
      "Paralog_end",
      "Paralog_strand",
      "Paralog_blocks",
      "Canonical_start",
      "Canonical_end",
      "Canonical_strand",
  ]
  
  with OUTPUT.open("w", newline="") as outfile:
  
      writer = csv.DictWriter(
          outfile,
          fieldnames=fields,
          delimiter="\t"
      )
  
      writer.writeheader()
      writer.writerows(candidate_pairs)
  
  
  # ------------------------------------------------------------
  # Resumen
  # ------------------------------------------------------------
  print()
  print(f"Created: {OUTPUT}")
  print(f"Candidate FAU paralogs: {len(candidate_pairs)}")
  EOF


Execute
~~~~~~~~

.. code-block:: bash

  chmod +x ANALYSES/FAU/scripts/02_identify_FAU_paralogs.py
  python ANALYSES/FAU/scripts/02_identify_FAU_paralogs.py

You should see something like this:

.. code-block::

  Created: /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/FAU/intermediate/FAU_same_scaffold_paralogs.tsv
  Candidate FAU paralogs: 15

Now we check if that sh*t worked:

.. code-block:: bash

  column -t ANALYSES/FAU/intermediate/FAU_same_scaffold_paralogs.tsv

.. code-block::

  Species                   Scaffold            Paralog_start  Paralog_end  Paralog_strand  Paralog_blocks  Canonical_start  Canonical_end  Canonical_strand
  Carollia_perspicillata    scaffold_4          119447942      119448261    -               2               119669674        119670778      +
  Cynopterus_sphinx         CM057383            143466352      143467136    +               3               143260443        143261494      -
  Eonycteris_spelaea        s4                  119662470      119662792    +               2               119439448        119440498      -
  Hipposideros_abae         HAP1_SUPER_7        6578634        6579422      -               3               6831441          6832501        +
  Hipposideros_armiger      LG08                125960229      125961015    +               3               125712363        125713424      -
  Hipposideros_caffer       HAP1_SUPER_8        4797044        4797837      -               3               5058275          5059335        +
  Hipposideros_jonesi       HAP1_SUPER_8        128305341      128306134    +               3               128033274        128034331      -
  Hipposideros_swinhoei     LG07                129900974      129901760    +               3               129638329        129639390      -
  Lonchorhina_inusitata     manual_scaffold_6   159124133      159124458    +               2               158902047        158903107      -
  Miniopterus_australis     manual_scaffold_11  4912507        4913416      -               3               5106650          5107723        +
  Miniopterus_natalensis    SUPER_11            70958412       70959318     +               3               70763943         70765016       -
  Miniopterus_schreibersii  OZ071104            68538622       68539528     +               3               68344261         68345316       -
  Rhinopoma_microphyllum    manual_scaffold_11  2593698        2594697      +               3               2352088          2353204        -
  Rhinopoma_muscatellum     manual_scaffold_13  2682494        2683279      +               3               2424464          2425571        -
  Rousettus_aegyptiacus     scaffold_m13_p_5    116490444      116490766    +               2               116269377        116270435      -
