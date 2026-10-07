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



3) Check if Charlie 8 is adjacent to the 3' end in the next 3 kb after the last exon of the paralog
----------------------------------------------------------------------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

.. code-block:: bash

  cat > ANALYSES/FAU/scripts/03_extract_FAU_3prime_3kb.sh <<'EOF'
  #!/bin/bash
  set -euo pipefail
  
  # ============================================================
  # 03_extract_FAU_3prime_3kb.sh
  #
  # Objective:
  # Extract the first 3 kb of genomic sequence immediately
  # downstream (3') of the terminal TOGA-projected exon/block
  # of each candidate FAU paralog.
  #
  # Expected architecture:
  #
  #   Paralog on + strand:
  #
  #       FAU ------>
  #                  |--- first 3 kb 3' ---|
  #
  #   Paralog on - strand:
  #
  #       |--- first 3 kb 3' ---|
  #                              <------ FAU
  #
  # In the BED12 TOGA projection:
  #
  #   + strand:
  #       paralog_end corresponds to the outer boundary of the
  #       terminal 3' block.
  #
  #   - strand:
  #       paralog_start corresponds to the outer boundary of the
  #       terminal 3' block.
  #
  # BED coordinates are 0-based, half-open.
  # samtools faidx coordinates are 1-based, inclusive.
  #
  # Input:
  #   intermediate/FAU_same_scaffold_paralogs.tsv
  #
  # Output:
  #   sequences/3prime_3kb/*.fa
  #   intermediate/FAU_3prime_3kb_regions.tsv
  # ============================================================
  
  
  ROOT="/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  
  INPUT="$ROOT/ANALYSES/FAU/intermediate/FAU_same_scaffold_paralogs.tsv"
  
  OUTDIR="$ROOT/ANALYSES/FAU/sequences/3prime_3kb"
  
  REGIONS="$ROOT/ANALYSES/FAU/intermediate/FAU_3prime_3kb_regions.tsv"
  
  
  # ------------------------------------------------------------
  # Check that samtools is available
  # ------------------------------------------------------------
  if ! command -v samtools >/dev/null 2>&1; then
      echo "ERROR: samtools is not available in the current environment." >&2
      echo "Activate samtools_env before running this script." >&2
      exit 1
  fi
  
  
  # ------------------------------------------------------------
  # Create output directory
  # ------------------------------------------------------------
  mkdir -p "$OUTDIR"
  
  
  # ------------------------------------------------------------
  # Metadata table describing every extracted region
  # ------------------------------------------------------------
  printf "Species\tScaffold\tParalog_start\tParalog_end\tStrand\tRegion_start\tRegion_end\tFASTA\n" \
      > "$REGIONS"
  
  
  # ------------------------------------------------------------
  # Read candidate paralogs
  # ------------------------------------------------------------
  tail -n +2 "$INPUT" | \
  while IFS=$'\t' read -r \
      species scaffold paralog_start paralog_end strand paralog_blocks \
      canonical_start canonical_end canonical_strand
  do
  
      echo "===== $species ====="
  
  
      # --------------------------------------------------------
      # Find the uncompressed genome FASTA.
      #
      # Ignore AppleDouble files generated by macOS (._filename).
      # --------------------------------------------------------
      fasta=$(find "$ROOT/SPECIES/$species" \
          -maxdepth 1 \
          -type f \
          -name "*.fa" \
          ! -name "._*" \
          -print)
  
  
      # --------------------------------------------------------
      # Require exactly one genome FASTA
      # --------------------------------------------------------
      fasta_count=$(printf '%s\n' "$fasta" | sed '/^$/d' | wc -l | tr -d ' ')
  
      if [[ "$fasta_count" -ne 1 ]]; then
          echo "ERROR: expected exactly one .fa for $species, found $fasta_count" >&2
          exit 1
      fi
  
  
      # --------------------------------------------------------
      # Create FASTA index if necessary
      # --------------------------------------------------------
      if [[ ! -f "${fasta}.fai" ]]; then
          echo "Creating FASTA index..."
          samtools faidx "$fasta"
      fi
  
  
      # --------------------------------------------------------
      # Determine first 3 kb immediately 3' of the FAU paralog.
      #
      # BED coordinates:
      #   start = 0-based
      #   end   = exclusive
      #
      # samtools coordinates:
      #   start = 1-based
      #   end   = inclusive
      # --------------------------------------------------------
  
      if [[ "$strand" == "+" ]]; then
  
          # BED interval:
          # [paralog_end, paralog_end + 3000)
          #
          # samtools:
          # paralog_end + 1 ... paralog_end + 3000
  
          region_start=$((paralog_end + 1))
          region_end=$((paralog_end + 3000))
  
  
      elif [[ "$strand" == "-" ]]; then
  
          # BED interval:
          # [paralog_start - 3000, paralog_start)
          #
          # Converted to samtools coordinates:
          # paralog_start - 2999 ... paralog_start
  
          region_start=$((paralog_start - 2999))
          region_end=$paralog_start
  
  
          # Prevent coordinates below 1
          if [[ "$region_start" -lt 1 ]]; then
              region_start=1
          fi
  
  
      else
          echo "ERROR: invalid strand '$strand' for $species" >&2
          exit 1
      fi
  
  
      # --------------------------------------------------------
      # Define samtools region
      # --------------------------------------------------------
      region="${scaffold}:${region_start}-${region_end}"
  
  
      # --------------------------------------------------------
      # Output FASTA
      # --------------------------------------------------------
      outfile="$OUTDIR/${species}_FAU_paralog_3prime_3kb.fa"
  
  
      echo "FASTA:  $(basename "$fasta")"
      echo "Region: $region"
      echo "Output: $(basename "$outfile")"
  
  
      # --------------------------------------------------------
      # Extract sequence.
      #
      # Sequence is kept in genomic orientation.
      # It is NOT reverse-complemented.
      # --------------------------------------------------------
      samtools faidx "$fasta" "$region" > "$outfile"
  
  
      # --------------------------------------------------------
      # Record extraction metadata
      # --------------------------------------------------------
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$species" \
          "$scaffold" \
          "$paralog_start" \
          "$paralog_end" \
          "$strand" \
          "$region_start" \
          "$region_end" \
          "$(basename "$fasta")" \
          >> "$REGIONS"
  
  done
  
  
  echo
  echo "========================================"
  echo "Finished"
  echo "========================================"
  
  echo
  echo "Extracted FASTA files:"
  find "$OUTDIR" -type f -name "*.fa" | wc -l
  
  echo
  echo "Region table:"
  echo "$REGIONS"
  
  EOF

Execute
~~~~~~~~

.. code-block:: bash

  chmod +x ANALYSES/FAU/scripts/03_extract_FAU_3prime_3kb.sh
  ./ANALYSES/FAU/scripts/03_extract_FAU_3prime_3kb.sh


Verify 
~~~~~~

.. code-block:: bash

  for f in ANALYSES/FAU/sequences/3prime_10kb/*.fa; do
      printf "%s\t" "$(basename "$f")"
      grep -v '^>' "$f" | tr -d '\n' | wc -c
  done

You should see something like this:

.. code-block::

  Carollia_perspicillata_FAU_paralog_3prime_3kb.fa	    3000
  Cynopterus_sphinx_FAU_paralog_3prime_3kb.fa	    3000
  Eonycteris_spelaea_FAU_paralog_3prime_3kb.fa	    3000
  Hipposideros_abae_FAU_paralog_3prime_3kb.fa	    3000
  Hipposideros_armiger_FAU_paralog_3prime_3kb.fa	    3000
  Hipposideros_caffer_FAU_paralog_3prime_3kb.fa	    3000
  Hipposideros_jonesi_FAU_paralog_3prime_3kb.fa	    3000
  Hipposideros_swinhoei_FAU_paralog_3prime_3kb.fa	    3000
  Lonchorhina_inusitata_FAU_paralog_3prime_3kb.fa	    3000
  Miniopterus_australis_FAU_paralog_3prime_3kb.fa	    3000
  Miniopterus_natalensis_FAU_paralog_3prime_3kb.fa	    3000
  Miniopterus_schreibersii_FAU_paralog_3prime_3kb.fa	    3000
  Rhinopoma_microphyllum_FAU_paralog_3prime_3kb.fa	    3000
  Rhinopoma_muscatellum_FAU_paralog_3prime_3kb.fa	    3000
  Rousettus_aegyptiacus_FAU_paralog_3prime_3kb.fa	    3000



4) Crear un FASTA de los terminales para enviar a CENSOR
---------------------------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

.. code-block:: bash

  cat > ANALYSES/FAU/scripts/04_make_FAU_CENSOR_fasta.sh <<'EOF'
  #!/bin/bash
  set -euo pipefail
  
  # ============================================================
  # 04_make_FAU_CENSOR_fasta.sh
  #
  # Objective:
  # Combine the 15 FAU paralog 3' 3-kb genomic regions into
  # a single multi-FASTA file for submission to CENSOR/Repbase.
  #
  # Each FASTA header contains:
  #   - species
  #   - scaffold
  #   - genomic coordinates
  #   - paralog strand
  #
  # Input:
  #   intermediate/FAU_3prime_3kb_regions.tsv
  #   sequences/3prime_3kb/*.fa
  #
  # Output:
  #   sequences/FAU_15_paralogs_3prime_10kb_CENSOR.fa
  # ============================================================
  
  ROOT="/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  
  REGIONS="$ROOT/ANALYSES/FAU/intermediate/FAU_3prime_3kb_regions.tsv"
  
  SEQDIR="$ROOT/ANALYSES/FAU/sequences/3prime_3kb"
  
  OUT="$ROOT/ANALYSES/FAU/sequences/FAU_15_paralogs_3prime_3kb_CENSOR.fa"
  
  
  # Start with an empty output file
  > "$OUT"
  
  
  # Read extraction metadata
  tail -n +2 "$REGIONS" | \
  while IFS=$'\t' read -r \
      species scaffold paralog_start paralog_end strand region_start region_end fasta
  do
  
      infile="$SEQDIR/${species}_FAU_paralog_3prime_3kb.fa"
  
      if [[ ! -s "$infile" ]]; then
          echo "ERROR: missing sequence for $species" >&2
          exit 1
      fi
  
  
      # Informative FASTA header
      printf ">%s|%s:%s-%s|FAU_paralog_3prime_3kb|strand=%s\n" \
          "$species" \
          "$scaffold" \
          "$region_start" \
          "$region_end" \
          "$strand" \
          >> "$OUT"
  
  
      # Copy sequence, excluding original samtools header
      grep -v '^>' "$infile" >> "$OUT"
  
  done
  
  
  echo "Created:"
  echo "$OUT"
  
  echo
  echo "Sequences in multi-FASTA:"
  grep -c '^>' "$OUT"
  EOF

Execute
~~~~~~~

.. code-block::

  chmod +x ANALYSES/FAU/scripts/04_make_FAU_CENSOR_fasta.sh
  ./ANALYSES/FAU/scripts/04_make_FAU_CENSOR_fasta.sh


You should see this:

.. code-block::

  Created:
  /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/FAU/sequences/FAU_15_paralogs_3prime_3kb_CENSOR.fa


Clean the headers (or this thing will not work)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

  awk '
  /^>/ {
      h=$0
      sub(/^>/,"",h)
      split(h,a,"|")
      gsub(/_/,"",a[1])
      print ">" a[1]
      next
  }
  {print}
  ' ANALYSES/FAU/sequences/FAU_15_paralogs_3prime_3kb_CENSOR.fa \
  > ANALYSES/FAU/sequences/FAU_15_paralogs_3prime_3kb_CENSOR_simpleheaders.fa


5) Send the FASTA to CENSOR
-----------------------------------

Open this webpage:
https://www.girinst.org/censor/

you should se this:

https://www.girinst.org/cgi-bin/censor/show_results.cgi?id=40507&lib=root


At this point, we confirmed that all 15 paralogs have the TE immediately adjacent to the 3′ region of interest. We can therefore extract the FAU paralog exons, align the sequences from all 15 species, and assess how well conserved the paralog is relative to the parental FAU gene. We can then map the alignment results onto the phylogenetic tree and visualize the pattern of conservation across species.


6) Extraer los exones
----------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  mkdir -p ANALYSES/FAU/sequences/paralog_projected_blocks
  nano ANALYSES/FAU/scripts/06_extract_FAU_paralog_blocks.py

.. code-block:: python

  #!/usr/bin/env python3
  
  from pathlib import Path
  from collections import defaultdict
  import csv
  import subprocess
  import sys
  
  
  # ============================================================
  # Paths
  # ============================================================
  
  ROOT = Path("/Volumes/Expansion/project3/bat_HTF_genomic_analysis")
  
  FAU_DIR = ROOT / "ANALYSES" / "FAU"
  
  INPUT_TSV = (
      FAU_DIR
      / "intermediate"
      / "FAU_same_scaffold_paralogs.tsv"
  )
  
  OUTDIR = (
      FAU_DIR
      / "sequences"
      / "paralog_projected_blocks"
  )
  
  METADATA_OUT = (
      FAU_DIR
      / "intermediate"
      / "FAU_paralog_projected_blocks.tsv"
  )
  
  OUTDIR.mkdir(parents=True, exist_ok=True)
  
  
  # ============================================================
  # Utility functions
  # ============================================================
  
  def reverse_complement(seq):
      """
      Reverse-complement a DNA sequence, including IUPAC ambiguity codes.
      """
      table = str.maketrans(
          "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
          "TGCAYRKMSWVHDBNtgcayrkmswvhdbn"
      )
  
      return seq.translate(table)[::-1]
  
  
  def find_genome_fasta(species_dir):
      """
      Find the genome FASTA in the top level of the species directory.
  
      If more than one FASTA file is present, use the largest one.
      """
  
      candidates = []
  
      for pattern in ("*.fa", "*.fasta", "*.fna"):
          candidates.extend(
              species_dir.glob(pattern)
          )
  
      candidates = [
          p
          for p in candidates
          if p.is_file()
      ]
  
      if not candidates:
          raise FileNotFoundError(
              f"No genome FASTA found in {species_dir}"
          )
  
      return max(
          candidates,
          key=lambda p: p.stat().st_size
      )
  
  
  def ensure_fai(fasta):
      """
      Create a samtools FASTA index if it does not already exist.
      """
  
      fai = Path(
          str(fasta) + ".fai"
      )
  
      if not fai.exists():
  
          print(
              f"Indexing {fasta.name} ...",
              file=sys.stderr
          )
  
          subprocess.run(
              [
                  "samtools",
                  "faidx",
                  str(fasta)
              ],
              check=True
          )
  
  
  def fetch_sequence(
      fasta,
      chrom,
      start0,
      end0
  ):
      """
      Extract sequence with samtools faidx.
  
      Input coordinates use BED convention:
  
          start0 = 0-based inclusive
          end0   = 0-based exclusive
  
      samtools faidx uses:
  
          1-based inclusive coordinates.
      """
  
      start1 = start0 + 1
      end1 = end0
  
      region = (
          f"{chrom}:"
          f"{start1}-"
          f"{end1}"
      )
  
      result = subprocess.run(
          [
              "samtools",
              "faidx",
              str(fasta),
              region
          ],
          check=True,
          capture_output=True,
          text=True
      )
  
      lines = (
          result.stdout
          .strip()
          .splitlines()
      )
  
      if not lines:
          raise RuntimeError(
              f"No sequence returned for "
              f"{fasta} {region}"
          )
  
      seq = "".join(
          lines[1:]
      ).upper()
  
      expected_length = (
          end0 - start0
      )
  
      if len(seq) != expected_length:
  
          raise RuntimeError(
              f"Length mismatch for {region}: "
              f"expected {expected_length}, "
              f"got {len(seq)}"
          )
  
      return seq
  
  
  def find_projection_bed_record(
      bedfile,
      scaffold,
      start,
      end,
      strand
  ):
      """
      Find the exact FAU paralog projection in query_annotation.bed
      using scaffold, start, end, strand, and the FAU projection label.
      """
  
      matches = []
  
      with open(bedfile) as fh:
  
          for line in fh:
  
              if (
                  not line.strip()
                  or line.startswith("#")
              ):
                  continue
  
              fields = (
                  line
                  .rstrip("\n")
                  .split("\t")
              )
  
              if len(fields) < 12:
                  continue
  
              chrom = fields[0]
  
              bed_start = int(
                  fields[1]
              )
  
              bed_end = int(
                  fields[2]
              )
  
              projection = fields[3]
  
              bed_strand = fields[5]
  
              if (
                  chrom == scaffold
                  and bed_start == start
                  and bed_end == end
                  and bed_strand == strand
                  and "#FAU#" in projection
              ):
  
                  matches.append(
                      fields
                  )
  
      if len(matches) != 1:
  
          raise RuntimeError(
              f"Expected exactly one BED record for "
              f"{scaffold}:{start}-{end}({strand}), "
              f"found {len(matches)}"
          )
  
      return matches[0]
  
  
  def write_fasta(
      path,
      records
  ):
      """
      Write FASTA records using 80 nt per line.
      """
  
      with open(path, "w") as out:
  
          for header, seq in records:
  
              out.write(
                  f">{header}\n"
              )
  
              for i in range(
                  0,
                  len(seq),
                  80
              ):
  
                  out.write(
                      seq[i:i + 80]
                      + "\n"
                  )
  
  
  # ============================================================
  # Main
  # ============================================================
  
  def main():
  
      if not INPUT_TSV.exists():
  
          raise FileNotFoundError(
              f"Missing input table: "
              f"{INPUT_TSV}"
          )
  
  
      records_by_block = defaultdict(
          list
      )
  
      metadata = []
  
      species_processed = []
  
  
      # ========================================================
      # Read the 15 FAU paralog candidates
      # ========================================================
  
      with open(INPUT_TSV) as fh:
  
          reader = csv.DictReader(
              fh,
              delimiter="\t"
          )
  
  
          required_columns = {
              "Species",
              "Scaffold",
              "Paralog_start",
              "Paralog_end",
              "Paralog_strand"
          }
  
  
          missing = (
              required_columns
              - set(
                  reader.fieldnames
                  or []
              )
          )
  
  
          if missing:
  
              raise RuntimeError(
                  "Missing required columns "
                  "in input TSV: "
                  + ", ".join(
                      sorted(missing)
                  )
              )
  
  
          # ====================================================
          # Process each species
          # ====================================================
  
          for row in reader:
  
              species = (
                  row["Species"]
              )
  
              scaffold = (
                  row["Scaffold"]
              )
  
              paralog_start = int(
                  row["Paralog_start"]
              )
  
              paralog_end = int(
                  row["Paralog_end"]
              )
  
              strand = (
                  row["Paralog_strand"]
              )
  
  
              species_processed.append(
                  species
              )
  
  
              species_dir = (
                  ROOT
                  / "SPECIES"
                  / species
              )
  
  
              bedfile = (
                  species_dir
                  / "query_annotation.bed"
              )
  
  
              if not species_dir.exists():
  
                  raise FileNotFoundError(
                      f"Missing species directory: "
                      f"{species_dir}"
                  )
  
  
              if not bedfile.exists():
  
                  raise FileNotFoundError(
                      f"Missing BED file: "
                      f"{bedfile}"
                  )
  
  
              # =================================================
              # Locate genome assembly
              # =================================================
  
              fasta = find_genome_fasta(
                  species_dir
              )
  
  
              ensure_fai(
                  fasta
              )
  
  
              print(
                  f"\n{species}",
                  file=sys.stderr
              )
  
  
              print(
                  f"  genome: "
                  f"{fasta.name}",
                  file=sys.stderr
              )
  
  
              print(
                  f"  paralog: "
                  f"{scaffold}:"
                  f"{paralog_start}-"
                  f"{paralog_end}"
                  f"({strand})",
                  file=sys.stderr
              )
  
  
              # =================================================
              # Find exact FAU projection in BED12
              # =================================================
  
              bed = find_projection_bed_record(
                  bedfile,
                  scaffold,
                  paralog_start,
                  paralog_end,
                  strand
              )
  
  
              chrom_start = int(
                  bed[1]
              )
  
  
              projection = (
                  bed[3]
              )
  
  
              block_count = int(
                  bed[9]
              )
  
  
              block_sizes = [
                  int(x)
                  for x in
                  bed[10]
                  .rstrip(",")
                  .split(",")
                  if x
              ]
  
  
              block_starts = [
                  int(x)
                  for x in
                  bed[11]
                  .rstrip(",")
                  .split(",")
                  if x
              ]
  
  
              # =================================================
              # Sanity checks
              # =================================================
  
              if (
                  len(block_sizes)
                  != block_count
              ):
  
                  raise RuntimeError(
                      f"{species}: "
                      f"blockSizes count "
                      f"({len(block_sizes)}) "
                      f"does not match "
                      f"blockCount "
                      f"({block_count})"
                  )
  
  
              if (
                  len(block_starts)
                  != block_count
              ):
  
                  raise RuntimeError(
                      f"{species}: "
                      f"blockStarts count "
                      f"({len(block_starts)}) "
                      f"does not match "
                      f"blockCount "
                      f"({block_count})"
                  )
  
  
              print(
                  f"  TOGA blocks: "
                  f"{block_count}",
                  file=sys.stderr
              )
  
  
              # =================================================
              # Build genomic blocks
              #
              # BED block order is always low genomic coordinate
              # to high genomic coordinate.
              # =================================================
  
              genomic_blocks = []
  
  
              for (
                  genomic_index,
                  (
                      rel_start,
                      size
                  )
              ) in enumerate(
                  zip(
                      block_starts,
                      block_sizes
                  ),
                  start=1
              ):
  
                  start0 = (
                      chrom_start
                      + rel_start
                  )
  
                  end0 = (
                      start0
                      + size
                  )
  
  
                  genomic_blocks.append(
                      {
                          "bed_block":
                              genomic_index,
  
                          "start0":
                              start0,
  
                          "end0":
                              end0,
  
                          "size":
                              size
                      }
                  )
  
  
              # =================================================
              # Convert BED genomic order to transcript order
              #
              # + strand:
              # genomic order = transcript order
              #
              # - strand:
              # genomic order must be reversed
              # =================================================
  
              if strand == "+":
  
                  transcript_blocks = (
                      genomic_blocks
                  )
  
  
              elif strand == "-":
  
                  transcript_blocks = list(
                      reversed(
                          genomic_blocks
                      )
                  )
  
  
              else:
  
                  raise RuntimeError(
                      f"{species}: "
                      f"unexpected strand "
                      f"value '{strand}'"
                  )
  
  
              # =================================================
              # Extract each projected block
              # =================================================
  
              for (
                  transcript_index,
                  block
              ) in enumerate(
                  transcript_blocks,
                  start=1
              ):
  
                  start0 = (
                      block["start0"]
                  )
  
                  end0 = (
                      block["end0"]
                  )
  
  
                  seq = fetch_sequence(
                      fasta,
                      scaffold,
                      start0,
                      end0
                  )
  
  
                  # =============================================
                  # Put sequence in transcript orientation
                  # =============================================
  
                  if strand == "-":
  
                      seq = reverse_complement(
                          seq
                      )
  
  
                  # =============================================
                  # Save sequence for multispecies FASTA
                  # =============================================
  
                  records_by_block[
                      transcript_index
                  ].append(
                      (
                          species,
                          seq
                      )
                  )
  
  
                  # =============================================
                  # Save metadata
                  # =============================================
  
                  metadata.append(
                      {
                          "Species":
                              species,
  
                          "Projection":
                              projection,
  
                          "Scaffold":
                              scaffold,
  
                          "Strand":
                              strand,
  
                          "Transcript_block":
                              transcript_index,
  
                          "BED_block":
                              block[
                                  "bed_block"
                              ],
  
                          "Start_0based":
                              start0,
  
                          "End_0based":
                              end0,
  
                          "Start_1based":
                              start0 + 1,
  
                          "End_1based":
                              end0,
  
                          "Length":
                              len(seq),
  
                          "Genome_FASTA":
                              fasta.name
                      }
                  )
  
  
                  print(
                      f"    transcript block "
                      f"{transcript_index}: "
                      f"{start0 + 1}-"
                      f"{end0} "
                      f"({len(seq)} bp)",
                      file=sys.stderr
                  )
  
  
      # ========================================================
      # Global sanity checks
      # ========================================================
  
      number_species = len(
          species_processed
      )
  
  
      if (
          len(
              set(species_processed)
          )
          != number_species
      ):
  
          raise RuntimeError(
              "Duplicate species found in "
              "FAU_same_scaffold_paralogs.tsv"
          )
  
  
      if number_species != 15:
  
          print(
              f"\nWARNING: expected "
              f"15 FAU paralog candidates, "
              f"but found "
              f"{number_species}.",
              file=sys.stderr
          )
  
  
      # ========================================================
      # Write one multispecies FASTA per projected block
      # ========================================================
  
      for block_number in sorted(
          records_by_block
      ):
  
          outfile = (
              OUTDIR
              / (
                  f"FAU_paralog_"
                  f"projected_block"
                  f"{block_number}.fa"
              )
          )
  
  
          write_fasta(
              outfile,
              records_by_block[
                  block_number
              ]
          )
  
  
      # ========================================================
      # Write metadata table
      # ========================================================
  
      fieldnames = [
          "Species",
          "Projection",
          "Scaffold",
          "Strand",
          "Transcript_block",
          "BED_block",
          "Start_0based",
          "End_0based",
          "Start_1based",
          "End_1based",
          "Length",
          "Genome_FASTA"
      ]
  
  
      with open(
          METADATA_OUT,
          "w",
          newline=""
      ) as out:
  
          writer = csv.DictWriter(
              out,
              fieldnames=fieldnames,
              delimiter="\t"
          )
  
          writer.writeheader()
  
          writer.writerows(
              metadata
          )
  
  
      # ========================================================
      # Final summary
      # ========================================================
  
      print(
          "\n"
          "========================================"
      )
  
      print(
          "FAU paralog projected-block extraction"
      )
  
      print(
          "========================================"
      )
  
  
      print(
          f"\nSpecies processed: "
          f"{number_species}"
      )
  
  
      total_blocks = 0
  
  
      for block_number in sorted(
          records_by_block
      ):
  
          records = (
              records_by_block[
                  block_number
              ]
          )
  
  
          lengths = [
              len(seq)
              for _, seq
              in records
          ]
  
  
          total_blocks += (
              len(records)
          )
  
  
          length_string = ",".join(
              str(length)
              for length
              in lengths
          )
  
  
          print(
              f"\nProjected block "
              f"{block_number}:"
          )
  
  
          print(
              f"  sequences: "
              f"{len(records)}"
          )
  
  
          print(
              f"  min length: "
              f"{min(lengths)} bp"
          )
  
  
          print(
              f"  max length: "
              f"{max(lengths)} bp"
          )
  
  
          print(
              f"  lengths: "
              f"{length_string}"
          )
  
  
      print(
          f"\nTotal projected blocks "
          f"extracted: "
          f"{total_blocks}"
      )
  
  
      print(
          "\nFASTA files:"
      )
  
  
      print(
          f"  {OUTDIR}"
      )
  
  
      print(
          "\nMetadata:"
      )
  
  
      print(
          f"  {METADATA_OUT}"
      )
  
  
      print(
          "\nDone."
      )
  
  
  # ============================================================
  # Run
  # ============================================================
  
  if __name__ == "__main__":
      main()


Execute
~~~~~~~~

.. code-block:: bash

  chmod +x ANALYSES/FAU/scripts/06_extract_FAU_paralog_blocks.py
  python3 ANALYSES/FAU/scripts/06_extract_FAU_paralog_blocks.py


.. code-block::

Cynopterus_sphinx:        ATGCAGCTGTTTGTCCGCGCCCGGGAGCTGCACACTCTTGAAGTGACCAGCTTCGAGACAGTTGTCCAGATCAAA
Rousettus_aegyptiacus:    ATGCAGCTGTTTGTCCGTGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTCGAGACAGTTGCCCAGATCAAG
Eonycteris_spelaea:       ATGCAGCTGTTTGTCCGTGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTCGAGACAGTTGCCCAGATCAAG
Rhinopoma_microphyllum:   ATGCAGCTGTTTGTACGCGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCTTGGAGACAGTTGCCCAGATCAAG
Rhinopoma_muscatellum:    ATGCAGCTGTTTGTACGCGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCTTGGAGACAGTTGCCCAGATCAAG
Hipposideros_abae:        ATGCAGCTGTTTGTCCGCTCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACGGTTGCCCAGATCAAG
Hipposideros_armiger:     ATGCAGCTGTTTGTCCGCTCCCGGGAGCTGCACACTCTTGAAGTGAGCGGCCTGGAGACAGTTGCCCAGATCAAG
Hipposideros_caffer:      ATGCAGCTGTTTGTCCGCTCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACGGTTGCCCAGATCAAG
Hipposideros_jonesi:      ATGCAGCTGTTTGTCCGCTCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACAGTTGCCCAGATCAAG
Hipposideros_swinhoei:    ATGCAGCTGTTTGTCCGATCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACAGTTGCACAGATCAAG
Miniopterus_australis:    ATGCAGCTGTTTGTCCGCGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACAGTTGCCCAGATCAAG
Miniopterus_natalensis:   ATGCAGCTGTTTGTCCGCGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACAGTTGCCCAGATCAAG
Miniopterus_schreibersii: ATGCAGCTGTTTGTCCGCGCCCGGGAGCTGCACACTCTTGAAGTGACCGGCCTGGAGACAGTTGCCCAGATCAAG
Lonchorhina_inusitata:    ATGCAGC--TCTGTCGGCTCC-GGGAGCTGCACACTCTTAAAGGGACAGACCTGGAGACAGTTGCCCAAATCAAA
Carollia_perspicillata:   ---CAGC--TCTGTCCGCGCC-GGGAGCTGCACACTCTTGAAGTGACCGGC-TGGAGACAGCTGCCCAAATCAAA
