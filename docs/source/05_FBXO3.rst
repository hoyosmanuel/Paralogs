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
  cat > ANALYSES/FBXO3/sequences/FBXO3_test_blocks.fa <<'EOF'
  >FBXO3_block_A
  TTGTTGCTATGTCAATCGAAGACTTTGCCAGCTCTCAAGTCATGATCCGCTATGGAGAAGACATTGCAAAAGATACTGGCTGATATCTGA
  >FBXO3_block_B
  GGAAGACAAAGCACAGAAGAATCAGTGCTGGAAATCTCTCTTCATTGATACTTACTCTGATGTAGGAAGATATATTGACCATTATGCTGCTGTTAAAAAGGCCTGGGCCGATCTCAAGAAATATTTGGAGCCCAGACCTCCTCAGATGACTTTGTCTCTGCAAG
  EOF

Hagamos un index del scafold donde deberia estar el parálogo de FBXO3:
Este index está en una especie donde sabemos por el código anterior que esl paralogo está presente entonces este es un control prositivo:

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  samtools faidx \
  SPECIES/Rhinolophus_hipposideros/mRhiHip1.1.hap1.fa \
  OZ077427 \
  > ANALYSES/FBXO3/sequences/Rhinolophus_hipposideros_OZ077427.fa

HAcemos una pequeña base de datos de blast

mkdir -p ANALYSES/FBXO3/intermediate/blast_db

.. code-block:: bash

  makeblastdb \
  -in ANALYSES/FBXO3/sequences/Rhinolophus_hipposideros_OZ077427.fa \
  -dbtype nucl \
  -out ANALYSES/FBXO3/intermediate/blast_db/RhiHip_OZ077427

Buscamos las secuencias:

.. code-block:: bash

  blastn \
  -query ANALYSES/FBXO3/sequences/FBXO3_test_blocks.fa \
  -db ANALYSES/FBXO3/intermediate/blast_db/RhiHip_OZ077427 \
  -task blastn-short \
  -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

Building a new DB, current time: 10/07/2026 17:11:21
New DB name:   /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/FBXO3/intermediate/blast_db/RhiHip_OZ077427
New DB title:  ANALYSES/FBXO3/sequences/Rhinolophus_hipposideros_OZ077427.fa
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 1 sequences in 0.737371 seconds.

Ahí están las secuencias:

.. code-block::

  FBXO3_block_A	OZ077427	95.556	90	4	0	1	90	23161427	23161516	3.18e-35	147
  FBXO3_block_A	OZ077427	94.444	90	5	0	1	90	23138058	23138147	7.74e-33	139
  FBXO3_block_B	OZ077427	96.951	164	5	0	1	164	23140786	23140949	1.08e-76	285
  FBXO3_block_B	OZ077427	96.324	136	5	0	1	136	23163669	23163804	5.53e-60	230

Ahora que sabemos que funciona, lo unico único que tenemos que hacer es extender esto a las otras especies

.. code-block:: bash


  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  
  cat > ANALYSES/FBXO3/scripts/03_blast_FBXO3_blocks_all_species.py <<'EOF'
  #!/usr/bin/env python3
  
  import csv
  import subprocess
  import sys
  from collections import defaultdict
  from pathlib import Path
  
  
  # ============================================================
  # 03_blast_FBXO3_blocks_all_species.py
  #
  # Objective:
  #
  # For every species:
  #
  #   1. Identify the canonical FBXO3 projection
  #      (unique 11-block TOGA projection).
  #
  #   2. Determine the scaffold containing canonical FBXO3.
  #
  #   3. Extract that entire scaffold from the genome.
  #
  #   4. BLAST the two FBXO3 diagnostic sequences
  #      against that scaffold.
  #
  # This allows detection of additional FBXO3-like copies
  # even when TOGA did not annotate them as FBXO3.
  #
  # Outputs:
  #
  # ANALYSES/FBXO3/sequences/canonical_scaffolds/
  #
  # ANALYSES/FBXO3/results/
  #   FBXO3_blocks_all_species_blast.tsv
  #
  # ANALYSES/FBXO3/results/
  #   FBXO3_blocks_all_species_summary.tsv
  #
  # ============================================================
  
  
  ROOT = Path(
      "/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
  )
  
  SPECIES_DIR = ROOT / "SPECIES"
  
  HITS_FILE = (
      ROOT
      / "ANALYSES"
      / "FBXO3"
      / "intermediate"
      / "FBXO3_all_hits.tsv"
  )
  
  QUERY = (
      ROOT
      / "ANALYSES"
      / "FBXO3"
      / "sequences"
      / "FBXO3_test_blocks.fa"
  )
  
  SCAFFOLD_DIR = (
      ROOT
      / "ANALYSES"
      / "FBXO3"
      / "sequences"
      / "canonical_scaffolds"
  )
  
  RESULTS_DIR = (
      ROOT
      / "ANALYSES"
      / "FBXO3"
      / "results"
  )
  
  RAW_OUT = (
      RESULTS_DIR
      / "FBXO3_blocks_all_species_blast.tsv"
  )
  
  SUMMARY_OUT = (
      RESULTS_DIR
      / "FBXO3_blocks_all_species_summary.tsv"
  )
  
  
  SCAFFOLD_DIR.mkdir(
      parents=True,
      exist_ok=True
  )
  
  RESULTS_DIR.mkdir(
      parents=True,
      exist_ok=True
  )
  
  
  
  # ============================================================
  # Check required programs
  # ============================================================
  
  import shutil
  
  for program in ["samtools", "blastn"]:
  
      if shutil.which(program) is None:
  
          sys.exit(
              f"ERROR: {program} was not found "
              f"in the current environment."
          )
  
  
  # ============================================================
  # Find genome FASTA
  #
  # Same general strategy used in the FAU workflow:
  # choose the largest top-level FASTA in each species folder.
  # ============================================================
  
  def find_genome_fasta(species_dir):
  
      candidates = []
  
      for pattern in [
          "*.fa",
          "*.fasta",
          "*.fna"
      ]:
  
          candidates.extend(
              species_dir.glob(pattern)
          )
  
      if not candidates:
  
          return None
  
      return max(
          candidates,
          key=lambda p: p.stat().st_size
      )
  
  
  # ============================================================
  # Ensure FASTA index
  # ============================================================
  
  def ensure_fai(genome):
  
      fai = Path(
          str(genome) + ".fai"
      )
  
      if not fai.exists():
  
          print(
              f"  Indexing genome: "
              f"{genome.name}"
          )
  
          subprocess.run(
              [
                  "samtools",
                  "faidx",
                  str(genome)
              ],
              check=True
          )
  
  
  # ============================================================
  # Read FBXO3 TOGA hits
  # ============================================================
  
  hits_by_species = defaultdict(list)
  
  with HITS_FILE.open() as infile:
  
      reader = csv.DictReader(
          infile,
          delimiter="\t"
      )
  
      for row in reader:
  
          row["Blocks"] = int(
              row["Blocks"]
          )
  
          hits_by_species[
              row["Species"]
          ].append(row)
  
  
  # ============================================================
  # Identify canonical FBXO3 scaffold for every species
  #
  # Canonical operational definition:
  # unique 11-block projection.
  # ============================================================
  
  canonical_by_species = {}
  
  
  for species, hits in sorted(
      hits_by_species.items()
  ):
  
      canonical = [
          h
          for h in hits
          if h["Blocks"] == 11
      ]
  
      if len(canonical) != 1:
  
          print(
              f"WARNING: {species}: "
              f"{len(canonical)} "
              f"11-block projections."
          )
  
          continue
  
      canonical_by_species[
          species
      ] = canonical[0]
  
  
  print()
  print(
      "Species with unique canonical FBXO3:",
      len(canonical_by_species)
  )
  print()
  
  
  # ============================================================
  # BLAST output columns
  # ============================================================
  
  raw_fields = [
      "Species",
      "Canonical_scaffold",
      "Canonical_start",
      "Canonical_end",
      "Canonical_strand",
      "Query",
      "Subject",
      "Percent_identity",
      "Alignment_length",
      "Query_length",
      "Query_start",
      "Query_end",
      "Subject_start",
      "Subject_end",
      "Evalue",
      "Bitscore",
      "Query_coverage"
  ]
  
  
  raw_rows = []
  
  summary = []
  
  
  # ============================================================
  # Process every species
  # ============================================================
  
  total_species = len(
      canonical_by_species
  )
  
  
  for number, species in enumerate(
      sorted(canonical_by_species),
      start=1
  ):
  
      canonical = canonical_by_species[
          species
      ]
  
      scaffold = canonical[
          "Scaffold"
      ]
  
      species_dir = (
          SPECIES_DIR
          / species
      )
  
  
      print(
          f"[{number}/{total_species}] "
          f"{species}"
      )
  
      print(
          f"  canonical scaffold: "
          f"{scaffold}"
      )
  
  
      # --------------------------------------------------------
      # Find genome
      # --------------------------------------------------------
  
      genome = find_genome_fasta(
          species_dir
      )
  
      if genome is None:
  
          print(
              "  WARNING: genome FASTA "
              "not found."
          )
  
          summary.append(
              {
                  "Species": species,
                  "Canonical_scaffold": scaffold,
                  "Block_A_hits": "NA",
                  "Block_B_hits": "NA",
                  "Status": "NO_GENOME"
              }
          )
  
          continue
  
  
      print(
          f"  genome: "
          f"{genome.name}"
      )
  
  
      # --------------------------------------------------------
      # Index genome if needed
      # --------------------------------------------------------
  
      ensure_fai(
          genome
      )
  
  
      # --------------------------------------------------------
      # Extract canonical FBXO3 scaffold
      # --------------------------------------------------------
  
      scaffold_fa = (
          SCAFFOLD_DIR
          / f"{species}.fa"
      )
  
  
      with scaffold_fa.open(
          "w"
      ) as outfile:
  
          result = subprocess.run(
              [
                  "samtools",
                  "faidx",
                  str(genome),
                  scaffold
              ],
              stdout=outfile,
              stderr=subprocess.PIPE,
              text=True
          )
  
  
      if result.returncode != 0:
  
          print(
              f"  WARNING: could not "
              f"extract scaffold {scaffold}"
          )
  
          print(
              result.stderr.strip()
          )
  
          summary.append(
              {
                  "Species": species,
                  "Canonical_scaffold": scaffold,
                  "Block_A_hits": "NA",
                  "Block_B_hits": "NA",
                  "Status": "SCAFFOLD_EXTRACTION_FAILED"
              }
          )
  
          continue
  
  
      # --------------------------------------------------------
      # BLAST the two query sequences
      #
      # word_size 7 gives more sensitivity than default blastn
      # for these relatively short sequences and allows
      # detection of diverged paralog copies.
      #
      # No identity cutoff is imposed here.
      # We want to inspect the complete evidence first.
      # --------------------------------------------------------
  
      blast_command = [
          "blastn",
  
          "-query",
          str(QUERY),
  
          "-subject",
          str(scaffold_fa),
  
          "-task",
          "blastn",
  
          "-word_size",
          "7",
  
          "-dust",
          "no",
  
          "-evalue",
          "1e-5",
  
          "-outfmt",
          (
              "6 "
              "qseqid "
              "sseqid "
              "pident "
              "length "
              "qlen "
              "qstart "
              "qend "
              "sstart "
              "send "
              "evalue "
              "bitscore "
              "qcovhsp"
          )
      ]
  
  
      result = subprocess.run(
          blast_command,
          capture_output=True,
          text=True,
          check=True
      )
  
  
      A_hits = 0
      B_hits = 0
  
  
      for line in result.stdout.splitlines():
  
          if not line.strip():
  
              continue
  
          fields = line.split("\t")
  
          (
              query,
              subject,
              pident,
              length,
              qlen,
              qstart,
              qend,
              sstart,
              send,
              evalue,
              bitscore,
              qcov
          ) = fields
  
  
          if query == "FBXO3_block_A":
  
              A_hits += 1
  
  
          elif query == "FBXO3_block_B":
  
              B_hits += 1
  
  
          raw_rows.append(
              {
                  "Species":
                      species,
  
                  "Canonical_scaffold":
                      scaffold,
  
                  "Canonical_start":
                      canonical["Start"],
  
                  "Canonical_end":
                      canonical["End"],
  
                  "Canonical_strand":
                      canonical["Strand"],
  
                  "Query":
                      query,
  
                  "Subject":
                      subject,
  
                  "Percent_identity":
                      pident,
  
                  "Alignment_length":
                      length,
  
                  "Query_length":
                      qlen,
  
                  "Query_start":
                      qstart,
  
                  "Query_end":
                      qend,
  
                  "Subject_start":
                      sstart,
  
                  "Subject_end":
                      send,
  
                  "Evalue":
                      evalue,
  
                  "Bitscore":
                      bitscore,
  
                  "Query_coverage":
                      qcov
              }
          )
  
  
      # --------------------------------------------------------
      # Simple descriptive summary
      # --------------------------------------------------------
  
      if (
          A_hits >= 2
          and
          B_hits >= 2
      ):
  
          status = (
              "MULTIPLE_HITS_BOTH_BLOCKS"
          )
  
      elif (
          A_hits >= 2
          or
          B_hits >= 2
      ):
  
          status = (
              "MULTIPLE_HITS_ONE_BLOCK"
          )
  
      elif (
          A_hits >= 1
          and
          B_hits >= 1
      ):
  
          status = (
              "SINGLE_PATTERN"
          )
  
      else:
  
          status = (
              "INCOMPLETE_OR_NO_HITS"
          )
  
  
      summary.append(
          {
              "Species":
                  species,
  
              "Canonical_scaffold":
                  scaffold,
  
              "Block_A_hits":
                  A_hits,
  
              "Block_B_hits":
                  B_hits,
  
              "Status":
                  status
          }
      )
  
  
      print(
          f"  Block A hits: "
          f"{A_hits}"
      )
  
      print(
          f"  Block B hits: "
          f"{B_hits}"
      )
  
      print(
          f"  {status}"
      )
  
      print()
  
  
  # ============================================================
  # Write raw BLAST hits
  # ============================================================
  
  with RAW_OUT.open(
      "w",
      newline=""
  ) as outfile:
  
      writer = csv.DictWriter(
          outfile,
          fieldnames=raw_fields,
          delimiter="\t"
      )
  
      writer.writeheader()
  
      writer.writerows(
          raw_rows
      )
  
  
  # ============================================================
  # Write summary
  # ============================================================
  
  summary_fields = [
      "Species",
      "Canonical_scaffold",
      "Block_A_hits",
      "Block_B_hits",
      "Status"
  ]
  
  
  with SUMMARY_OUT.open(
      "w",
      newline=""
  ) as outfile:
  
      writer = csv.DictWriter(
          outfile,
          fieldnames=summary_fields,
          delimiter="\t"
      )
  
      writer.writeheader()
  
      writer.writerows(
          summary
      )
  
  
  # ============================================================
  # Final report
  # ============================================================
  
  print()
  print(
      "========================================"
  )
  
  print(
      "FBXO3 scaffold-wide BLAST complete"
  )
  
  print(
      "========================================"
  )
  
  print()
  
  print(
      f"Raw hits:"
  )
  
  print(
      RAW_OUT
  )
  
  print()
  
  print(
      f"Summary:"
  )
  
  print(
      SUMMARY_OUT
  )
  
  print()
  
  print(
      "Species with >=2 hits for BOTH blocks:"
  )
  
  for row in summary:
  
      if (
          row["Status"]
          == "MULTIPLE_HITS_BOTH_BLOCKS"
      ):
  
          print(
              f"  {row['Species']} "
              f"A={row['Block_A_hits']} "
              f"B={row['Block_B_hits']}"
          )


.. code-block:: bash

  chmod +x ANALYSES/FBXO3/scripts/03_blast_FBXO3_blocks_all_species.py
  python3 ANALYSES/FBXO3/scripts/03_blast_FBXO3_blocks_all_species.py



Results:
~~~~~~~~~
.. code-block:: 

  ========================================
  FBXO3 scaffold-wide BLAST complete
  ========================================
  
  Raw hits:
  /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/FBXO3/results/FBXO3_blocks_all_species_blast.tsv
  
  Summary:
  /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/FBXO3/results/FBXO3_blocks_all_species_summary.tsv
  
  Species with >=2 hits for BOTH blocks:
    Hipposideros_abae A=2 B=2
    Hipposideros_armiger A=2 B=2
    Hipposideros_caffer A=2 B=2
    Hipposideros_jonesi A=2 B=2
    Hipposideros_larvatus A=2 B=2
    Hipposideros_swinhoei A=2 B=2
    Rhinolophus_affinis A=2 B=2
    Rhinolophus_ferrumequinum A=2 B=2
    Rhinolophus_foetidus A=2 B=2
    Rhinolophus_hipposideros A=2 B=2
    Rhinolophus_pearsonii A=2 B=2
    Rhinolophus_perniger_lanosus A=2 B=2
    Rhinolophus_sinicus A=2 B=2
    Rhinolophus_trifoliatus A=2 B=2
    Triaenops_persicus A=2 B=2
