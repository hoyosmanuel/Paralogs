I. IDENTIFICAR EL PARÁLOGO DE ZKSCAN7
=====================================

The ZKSCAN7 chimeric architecture was initially identified in
Rhinolophus hipposideros. The architecture contains host-derived
ZKSCAN7 exons followed by sequence derived from the transposable
element RhiSin-1.24.

To identify homologous genomic architectures independently of gene
annotation, exons 1–5 are used as host-gene probes. The TE-derived
exon 6 is analyzed separately as an RhiSin-1.24 probe.

1. Create the ZKSCAN7 host exon probes
---------------------------------------

Create an independent working directory for the ZKSCAN7 analysis

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  mkdir -p ANALYSES/ZKSCAN7/{intermediate,results,scripts,sequences}
  cd ANALYSES/ZKSCAN7
  

Create the host-exon FASTA

.. code-block:: bash

  nano sequences/ZKSCAN7_host_exon_probes.fa

Add the following sequences

.. code-block:: bash

  >ZKSCAN7_exon1_match
  AGTGAAGCTTCTCGCACTTCCGGCTCGGCCATTGTTTTCGCTTCTGCCGGCCGCGGTGGTGACCGTCGGCCCGACATCGCCGCCTTGGTGCCCCAGTTGGAGGACCGGTGGCGGGGAGAGAGGCGCTGTCGGGGCCGAGGTCGTAAGGCTGGAGCCACCGCCGCCCACGGCCGATCCCACGCTTCGCGCCTCCGGCCAGCGGGCGCGGTCCCTGTGATGCTCGGAAGCCGCCCGGTGCAGAGACTTCTCTTTGTCCCATTCCCCTGACCAGCCGTGCCTGGCGCGGAGCCGGACACCTGCGGGTACGGCGATGAACAGGGCGGACGTGGCCAG
  
  >ZKSCAN7_exon2_match
  GTGACACTCTCATCTGCCCACACCCCCAACTCTGAAAAACAGCTCCTGAGACCTGAATTACAGATTATCATTGGACAAACTGCTTTTAATCAAACAGACTGCTTCAGGAACAGGTTCTTTGGGGTCACAATGGCCACCCAAGGCAGGGGAACTTTAGGCCTCATCCCGAGGGCTGCTGTTCTCCAGAAGCAGGAGGAGCACCTGACTGTGAAGCAGGAGCCGTGCCTGACTGTGAAGCAGGAGCCGGGGAGCCAGACCTGGGGGCAGGGCTGCAGTCTCCAAAAGAGCCACCCTCCTGTCTGTGAAATTTTCCGGCTGCACTTCAGGCAGTTATGTTACCATGAGATGTCTGGGCCCCAAGAGGCGCTGAGCAGGCTCCGGGAGCTCTGTCACTGGTGGCTCATGCCCGAGGTGCACACCAAAGAGCAGATCCTGGAGCTGCTGGTGCTGGAGCAGTTCCTGAGCATCCTGCCTGCGGAGCTCCGGACCTGGGTGCAGCTGCATCACCCTGAGAGTGGCGAGCAGGCTGTGGCTGCGGTAGAGGATTTCCAGAGACACCTCGGTGGACCAGGAGAG
  
  >ZKSCAN7_exon3_match
  GTTTCAGCTCCAGCACAGGAACAGGAAATGCATTTGGAGGAGATGACAACCCTCGGTGCAACAGAGGAATCTCCTCCTACCTTACCCCTCAGTGAGGGTTCAGCCCCTGGAGCCCACCTTGAGCCTCCTCGTGTCCCAGGGGCACACCACCTCCCCAATGGGCACTCTG
  
  >ZKSCAN7_exon4_match
  CTCATGGTGCTTCCCGGGTGTCTGCCCTTCCCCGAGCAAGGCAGTTAGGAGACCAAGCAGTGGCACCTGTGCTTCGGATGGTCAGGCCCCAG
  
  >ZKSCAN7_exon5_match
  GGCTCTGCAGCATATGGGTTCCCGTCTGTGGACTGTACTGAGAAGAAGGGCAGAGGTCTAGCACGGAGTCGGAGAGGCCTGAACCGGAACATCAGGCCGGAAACGTACCGCAGCATGGCCTCGTTGG

Save the file and confirm that all five host probes are present

.. code-block:: bash

  grep '^>' sequences/ZKSCAN7_host_exon_probes.fa

The expected output is

.. code-block:: bash

  >ZKSCAN7_exon1_match
  >ZKSCAN7_exon2_match
  >ZKSCAN7_exon3_match
  >ZKSCAN7_exon4_match
  >ZKSCAN7_exon5_match


2) Create the RhiSin-1.24 TE-derived exon probe
-----------------------------------------------

Exon 6 corresponds to sequence derived from the transposable element
RhiSin-1.24. Because this sequence is TE-derived, it is analyzed
separately from the host-derived ZKSCAN7 exons.

From the ZKSCAN7 analysis directory

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/ZKSCAN7

Create the TE-derived probe FASTA

.. code-block:: bash

    nano sequences/ZKSCAN7_RhiSin-1.24_probe.fa

Add the following sequence

.. code-block:: bash

    >RhiSin-1.24_match
    acaagatgagtgcaaaaagaaagagctattccgtcgagtacaagaaaggaatcgtggaggactctcggggcaagaatcttacggctttctgcaaagagaagaagttggatctccgaatggtccgaaaatggcgagccgagtacgataacctcagtcaacaggtagacgaaggaaatgctaagaagcgcaagtgtggatcaggtcggcaaccattgtttcttgagctggaagacatcatctgtgagtagattgctgaaaggagagcaaagcctttggttgtgcacagggctgatattcaagcatttgtccttgcaatggcaccacagttagaaatatcccccgaaaaattcaaagcatcacaacactggctggatggcttccttccacaatatgaactgtctctaagaagatctaccacactgttctagctggaagatactgaagttattaaacatgcacttgcattcaagtccttcgttgatggcctcgatttttctaaatgccaactccaacatgattgctatgtatgaaactacagtgtttatgggccaaggatctcaaatgacaattgattagaggggtgcctcgtcaatctgctacattccctccactggttacgaaagtgcacgtgttatctgtattttggcaattcatctggttggaaagagagccccacctctaatcatcactaagggcaagaaagataaggttgaatgtgtttcaggcatttatgttcttgaaaccgaaaaagcctgctgcacacaagcagttaataaggaagtgggtcaatttaatgctgccacttgttttgcgaggtggccagagaggtctgctagtcttggattcagccagcactcactgcgctaaagacatgaagaactcccttgcagagaaaagaatagatcaaataatgattccctcaggaatgattgcctatctccagactcttgatactgcaataagcaagccattcaaggaacatttgctcatggatatccatgactacattgaaaatagaatggagagaaatcagtgtggaaactttgtgaagcctaatctgtaagaggttatgacttgggtgaagaattcatgagataaaaatcactgacagctgtgttgccaatgcactacgagcaggctacatggacaagaagtgctcatttcaggagagctccgttgctggacatgagaggttggggccaatgattttacaggaaatggagtccaaagaaattcaagccggaattcgtagtttgcagagttatgacaatgttccagaagaggatgacatgactatTTGAATAAATGTAGAATTTTGTACATGAATAAATGAATAAATGTAGATTGTtgtacatgaataaatgtagattgttgtacgtg

Remove accidental whitespace from the sequence

.. code-block:: bash

    sed -i '' '/^>/! s/[[:space:]]//g' \
    sequences/ZKSCAN7_RhiSin-1.24_probe.fa

Confirm the FASTA header

.. code-block:: bash

    head -1 sequences/ZKSCAN7_RhiSin-1.24_probe.fa

The expected output is

.. code-block:: bash

    >RhiSin-1.24_match

Finally, verify the probe length

.. code-block:: bash

    awk '
    /^>/ {
        if (name) print name, length(seq)
        name=$0
        seq=""
        next
    }
    {
        seq=seq $0
    }
    END {
        print name, length(seq)
    }
    ' sequences/ZKSCAN7_RhiSin-1.24_probe.fa



3) Test the ZKSCAN7 host and RhiSin-1.24 probes in Rhinolophus hipposideros
----------------------------------------------------------------------------

Use *Rhinolophus hipposideros* as the positive reference species.

From the ZKSCAN7 analysis directory

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/ZKSCAN7

Create the BLAST database directory

.. code-block:: bash

    mkdir -p intermediate/blast_db

Build a nucleotide BLAST database from the *Rhinolophus hipposideros*
genome

.. code-block:: bash

    makeblastdb \
    -in ../../SPECIES/Rhinolophus_hipposideros/*.fa \
    -dbtype nucl \
    -out intermediate/blast_db/Rhinolophus_hipposideros

Search the five host-derived ZKSCAN7 exons

.. code-block:: bash

    blastn \
    -query sequences/ZKSCAN7_host_exon_probes.fa \
    -db intermediate/blast_db/Rhinolophus_hipposideros \
    -task blastn \
    -word_size 7 \
    -dust no \
    -evalue 1e-10 \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
    > intermediate/ZKSCAN7_exon_matches_Rhinolophus_hipposideros.tsv

Search the RhiSin-1.24-derived sequence independently

.. code-block:: bash

    blastn \
    -query sequences/ZKSCAN7_RhiSin-1.24_probe.fa \
    -db intermediate/blast_db/Rhinolophus_hipposideros \
    -task blastn \
    -word_size 7 \
    -dust no \
    -evalue 1e-20 \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
    > intermediate/ZKSCAN7_RhiSin-1.24_matches_Rhinolophus_hipposideros.tsv

Combine both sets of matches

.. code-block:: bash

    cat \
    intermediate/ZKSCAN7_exon_matches_Rhinolophus_hipposideros.tsv \
    intermediate/ZKSCAN7_RhiSin-1.24_matches_Rhinolophus_hipposideros.tsv \
    > intermediate/ZKSCAN7_all_matches_Rhinolophus_hipposideros.tsv

Display the results grouped first by scaffold and then in the biological
probe order Exon1–Exon5 followed by RhiSin-1.24

.. code-block:: bash

    awk 'BEGIN{OFS="\t"}
    {
        if ($1=="ZKSCAN7_exon1_match")          order=1
        else if ($1=="ZKSCAN7_exon2_match")     order=2
        else if ($1=="ZKSCAN7_exon3_match")     order=3
        else if ($1=="ZKSCAN7_exon4_match")     order=4
        else if ($1=="ZKSCAN7_exon5_match")     order=5
        else if ($1=="RhiSin-1.24_match")        order=6
        else                                    order=99

        print $0, order
    }' intermediate/ZKSCAN7_all_matches_Rhinolophus_hipposideros.tsv \
    | sort -k2,2 -k13,13n \
    | cut -f1-12 \
    | column -t


ver:

.. code-block:: bash

  awk '$2=="OZ077431" {
      start=($8<$9?$8:$9)
      end=($8>$9?$8:$9)
  
      if (end >= 58250000 && start <= 58280000)
          print
  }' intermediate/ZKSCAN7_all_matches_Rhinolophus_hipposideros.tsv \
  | column -t



4) Search for the ZKSCAN7–RhiSin-1.24 architecture across species
------------------------------------------------------------------

Create the general detection script


.. code-block:: bash

    nano scripts/02_find_ZKSCAN7_chimeric_signature.py

Paste the following script


.. code-block:: python

  #!/usr/bin/env python3
  
  import sys
  import subprocess
  from pathlib import Path
  from collections import defaultdict
  
  
  # ============================================================
  # PARAMETERS
  # ============================================================
  
  MIN_EXON_IDENTITY = 70.0
  MIN_EXON_QCOV = 70.0
  
  MIN_TE_IDENTITY = 70.0
  MIN_TE_QCOV = 70.0
  
  MIN_EXONS_PER_LOCUS = 4
  
  # Maximum gap between neighboring host-exon hits
  MAX_EXON_GAP = 20000
  
  # Positive R. hipposideros control = 6651 bp
  MAX_TE_DISTANCE = 10000
  
  
  # ============================================================
  # PATHS
  # ============================================================
  
  SCRIPT = Path(__file__).resolve()
  ROOT = SCRIPT.parent.parent
  PROJECT = ROOT.parent.parent
  
  SPECIES_ROOT = PROJECT / "SPECIES"
  
  EXON_PROBES = ROOT / "sequences" / "ZKSCAN7_host_exon_probes.fa"
  TE_PROBE = ROOT / "sequences" / "ZKSCAN7_RhiSin-1.24_probe.fa"
  
  INTERMEDIATE = ROOT / "intermediate"
  RESULTS = ROOT / "results"
  
  
  # ============================================================
  # FUNCTIONS
  # ============================================================
  
  def run(cmd):
      subprocess.run(cmd, check=True)
  
  
  def strand_from_coordinates(sstart, send):
      return "+" if sstart <= send else "-"
  
  
  def interval(sstart, send):
      return min(sstart, send), max(sstart, send)
  
  
  def exon_number(query):
      try:
          return int(
              query.replace("ZKSCAN7_exon", "")
                   .replace("_match", "")
          )
      except Exception:
          return None
  
  
  def read_blast(path):
  
      hits = []
  
      with open(path) as fh:
  
          for line in fh:
  
              if not line.strip():
                  continue
  
              fields = line.rstrip().split("\t")
  
              if len(fields) != 12:
                  continue
  
              (
                  qseqid,
                  sseqid,
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
  
              sstart = int(sstart)
              send = int(send)
  
              start, end = interval(sstart, send)
  
              hits.append({
                  "query": qseqid,
                  "scaffold": sseqid,
                  "identity": float(pident),
                  "length": int(length),
                  "qlen": int(qlen),
                  "qstart": int(qstart),
                  "qend": int(qend),
                  "sstart": sstart,
                  "send": send,
                  "start": start,
                  "end": end,
                  "strand": strand_from_coordinates(sstart, send),
                  "evalue": evalue,
                  "bitscore": float(bitscore),
                  "qcov": float(qcov)
              })
  
      return hits
  
  
  # ============================================================
  # INPUT
  # ============================================================
  
  if len(sys.argv) != 2:
      sys.exit(
          "Usage: python3 scripts/02_find_ZKSCAN7_chimeric_signature.py "
          "<Species_name>"
      )
  
  species = sys.argv[1]
  
  species_dir = SPECIES_ROOT / species
  
  if not species_dir.exists():
      sys.exit(f"Species directory not found: {species_dir}")
  
  
  # ============================================================
  # FIND GENOME
  # ============================================================
  
  genomes = list(species_dir.glob("*.fa"))
  
  if len(genomes) == 0:
      sys.exit(f"No *.fa genome found for {species}")
  
  if len(genomes) > 1:
      sys.exit(
          f"More than one *.fa genome found for {species}: "
          + ", ".join(x.name for x in genomes)
      )
  
  genome = genomes[0]
  
  
  # ============================================================
  # DIRECTORIES
  # ============================================================
  
  species_intermediate = INTERMEDIATE / species
  species_intermediate.mkdir(parents=True, exist_ok=True)
  
  db_dir = INTERMEDIATE / "blast_db" / species
  db_dir.parent.mkdir(parents=True, exist_ok=True)
  
  RESULTS.mkdir(parents=True, exist_ok=True)
  
  db_prefix = str(db_dir)
  
  
  # ============================================================
  # BUILD BLAST DATABASE
  # ============================================================
  
  existing_db_files = list(
      db_dir.parent.glob(db_dir.name + ".*")
  )
  
  if not existing_db_files:
  
      print(f"Building BLAST database for {species}...")
  
      run([
          "makeblastdb",
          "-in", str(genome),
          "-dbtype", "nucl",
          "-out", db_prefix
      ])
  
  else:
  
      print(f"Reusing ZKSCAN7 BLAST database for {species}")
  
  
  # ============================================================
  # BLAST HOST EXONS
  # ============================================================
  
  exon_out = (
      species_intermediate /
      f"ZKSCAN7_exon_matches_{species}.tsv"
  )
  
  run([
      "blastn",
      "-query", str(EXON_PROBES),
      "-db", db_prefix,
      "-task", "blastn",
      "-word_size", "7",
      "-dust", "no",
      "-evalue", "1e-10",
      "-outfmt",
      "6 qseqid sseqid pident length qlen qstart qend "
      "sstart send evalue bitscore qcovhsp",
      "-out", str(exon_out)
  ])
  
  
  # ============================================================
  # BLAST RHISIN-1.24
  # ============================================================
  
  te_out = (
      species_intermediate /
      f"ZKSCAN7_RhiSin_matches_{species}.tsv"
  )
  
  run([
      "blastn",
      "-query", str(TE_PROBE),
      "-db", db_prefix,
      "-task", "blastn",
      "-word_size", "7",
      "-dust", "no",
      "-evalue", "1e-20",
      "-outfmt",
      "6 qseqid sseqid pident length qlen qstart qend "
      "sstart send evalue bitscore qcovhsp",
      "-out", str(te_out)
  ])
  
  
  # ============================================================
  # FILTER EXON HITS
  # ============================================================
  
  exon_hits = []
  
  for h in read_blast(exon_out):
  
      n = exon_number(h["query"])
  
      if n is None:
          continue
  
      if h["identity"] < MIN_EXON_IDENTITY:
          continue
  
      if h["qcov"] < MIN_EXON_QCOV:
          continue
  
      h["exon"] = n
  
      exon_hits.append(h)
  
  
  # ============================================================
  # FILTER TE HITS
  # ============================================================
  
  te_hits = [
      h for h in read_blast(te_out)
      if h["identity"] >= MIN_TE_IDENTITY
      and h["qcov"] >= MIN_TE_QCOV
  ]
  
  
  # ============================================================
  # GROUP HOST HITS BY SCAFFOLD AND STRAND
  # ============================================================
  
  groups = defaultdict(list)
  
  for h in exon_hits:
      groups[(h["scaffold"], h["strand"])].append(h)
  
  
  candidate_loci = []
  
  
  # ============================================================
  # RECONSTRUCT ZKSCAN7-LIKE LOCI
  # ============================================================
  
  for (scaffold, strand), hits in groups.items():
  
      hits = sorted(hits, key=lambda x: x["start"])
  
      clusters = []
      current = []
  
      for hit in hits:
  
          if not current:
              current = [hit]
              continue
  
          previous_end = max(x["end"] for x in current)
  
          gap = hit["start"] - previous_end
  
          if gap <= MAX_EXON_GAP:
              current.append(hit)
  
          else:
              clusters.append(current)
              current = [hit]
  
      if current:
          clusters.append(current)
  
  
      # ========================================================
      # EVALUATE EACH CLUSTER
      # ========================================================
  
      for cluster in clusters:
  
          best_by_exon = {}
  
          for h in cluster:
  
              n = h["exon"]
  
              if (
                  n not in best_by_exon
                  or h["bitscore"] > best_by_exon[n]["bitscore"]
              ):
                  best_by_exon[n] = h
  
  
          if len(best_by_exon) < MIN_EXONS_PER_LOCUS:
              continue
  
  
          selected = [
              best_by_exon[n]
              for n in sorted(best_by_exon)
          ]
  
  
          # ----------------------------------------------------
          # Verify transcriptional exon order
          # ----------------------------------------------------
  
          positions = [
              (h["start"] + h["end"]) / 2
              for h in selected
          ]
  
          if strand == "+":
  
              ordered = all(
                  positions[i] < positions[i + 1]
                  for i in range(len(positions) - 1)
              )
  
          else:
  
              ordered = all(
                  positions[i] > positions[i + 1]
                  for i in range(len(positions) - 1)
              )
  
  
          if not ordered:
              continue
  
  
          locus_start = min(
              h["start"] for h in selected
          )
  
          locus_end = max(
              h["end"] for h in selected
          )
  
  
          highest_exon = max(best_by_exon)
  
          terminal_exon = best_by_exon[highest_exon]
  
  
          # ----------------------------------------------------
          # Biological 3' end
          # ----------------------------------------------------
  
          if strand == "+":
              three_prime = terminal_exon["end"]
          else:
              three_prime = terminal_exon["start"]
  
  
          # ----------------------------------------------------
          # Find RhiSin-1.24 toward biological 3'
          # ----------------------------------------------------
  
          nearby_te = []
  
          for te in te_hits:
  
              if te["scaffold"] != scaffold:
                  continue
  
              if te["strand"] != strand:
                  continue
  
  
              if strand == "+":
  
                  distance = te["start"] - three_prime
  
                  if 0 <= distance <= MAX_TE_DISTANCE:
                      nearby_te.append(
                          (distance, te)
                      )
  
              else:
  
                  distance = three_prime - te["end"]
  
                  if 0 <= distance <= MAX_TE_DISTANCE:
                      nearby_te.append(
                          (distance, te)
                      )
  
  
          nearby_te.sort(
              key=lambda x: (
                  x[0],
                  -x[1]["bitscore"]
              )
          )
  
  
          # ----------------------------------------------------
          # CLASSIFY
          # ----------------------------------------------------
  
          exons_found = sorted(best_by_exon)
  
          if nearby_te:
  
              distance, te = nearby_te[0]
  
              if exons_found == [1, 2, 3, 4, 5]:
                  status = "ZKSCAN7_CHIMERIC_COMPLETE"
              else:
                  status = "ZKSCAN7_CHIMERIC_PARTIAL"
  
              te_start = te["start"]
              te_end = te["end"]
              te_identity = te["identity"]
              te_qcov = te["qcov"]
  
          else:
  
              status = "ZKSCAN7_LIKE_NO_3PRIME_RHISIN"
  
              distance = ""
              te_start = ""
              te_end = ""
              te_identity = ""
              te_qcov = ""
  
  
          candidate_loci.append({
  
              "species": species,
              "genome": genome.name,
              "scaffold": scaffold,
              "strand": strand,
  
              "locus_start": locus_start,
              "locus_end": locus_end,
  
              "exons": ",".join(
                  map(str, exons_found)
              ),
  
              "n_exons": len(exons_found),
  
              "three_prime": three_prime,
  
              "te_start": te_start,
              "te_end": te_end,
  
              "te_identity": te_identity,
              "te_qcov": te_qcov,
  
              "te_distance": distance,
  
              "status": status
          })
  
  
  # ============================================================
  # WRITE RESULTS
  # ============================================================
  
  result_file = (
      RESULTS /
      f"ZKSCAN7_signature_{species}.tsv"
  )
  
  header = [
  
      "Species",
      "Genome",
      "Scaffold",
      "Strand",
  
      "ZKSCAN7_start",
      "ZKSCAN7_end",
  
      "Exons_found",
      "Number_of_exons",
  
      "Biological_3prime",
  
      "RhiSin_start",
      "RhiSin_end",
  
      "RhiSin_identity",
      "RhiSin_qcov",
  
      "RhiSin_distance_3prime",
  
      "Status"
  ]
  
  
  with open(result_file, "w") as out:
  
      out.write(
          "\t".join(header) + "\n"
      )
  
      for x in sorted(
          candidate_loci,
          key=lambda z: (
              z["scaffold"],
              z["locus_start"]
          )
      ):
  
          row = [
  
              x["species"],
              x["genome"],
              x["scaffold"],
              x["strand"],
  
              x["locus_start"],
              x["locus_end"],
  
              x["exons"],
              x["n_exons"],
  
              x["three_prime"],
  
              x["te_start"],
              x["te_end"],
  
              x["te_identity"],
              x["te_qcov"],
  
              x["te_distance"],
  
              x["status"]
          ]
  
          out.write(
              "\t".join(map(str, row)) + "\n"
          )
  
  
  # ============================================================
  # HUMAN-READABLE SUMMARY
  # ============================================================
  
  print()
  print("=" * 72)
  print(f"ZKSCAN7 analysis: {species}")
  print("=" * 72)
  
  
  if not candidate_loci:
  
      print("No coherent ZKSCAN7-like locus detected.")
  
  
  else:
  
      for x in candidate_loci:
  
          print()
  
          print(
              f"{x['scaffold']}  "
              f"{x['locus_start']}-{x['locus_end']}  "
              f"strand {x['strand']}"
          )
  
          print(
              f"Exons: {x['exons']} "
              f"({x['n_exons']} exons)"
          )
  
  
          if x["status"].startswith(
              "ZKSCAN7_CHIMERIC"
          ):
  
              print(
                  f"RhiSin-1.24: "
                  f"{x['te_start']}-{x['te_end']}"
              )
  
              print(
                  f"Identity: "
                  f"{x['te_identity']:.3f}%"
              )
  
              print(
                  f"Coverage: "
                  f"{x['te_qcov']:.1f}%"
              )
  
              print(
                  f"3' distance: "
                  f"{x['te_distance']} bp"
              )
  
              print(
                  f"-> {x['status']}"
              )
  
          else:
  
              print(
                  "RhiSin-1.24: none within "
                  "3' search window"
              )
  
              print(
                  "-> ZKSCAN7_LIKE_NO_3PRIME_RHISIN"
              )
  
  
  print()
  print(f"Results: {result_file}")


Save the file and make it executable

.. code-block:: bash

    chmod +x scripts/02_find_ZKSCAN7_chimeric_signature.py

Before running the complete batch, the script can be checked against the
positive reference

.. code-block:: bash

    python3 scripts/02_find_ZKSCAN7_chimeric_signature.py \
    Rhinolophus_hipposideros

The expected positive locus should occur on

.. code-block:: 

    OZ077431

.. code-block:: bash

with exons 1–5 and a downstream RhiSin-1.24 match approximately 6.65 kb
from exon 5.

Once the positive control is recovered, create the batch-log directory

.. code-block:: bash

    mkdir -p results/batch_logs

    rm -f results/batch_logs/*.failed

.. code-block:: bash

Run the analysis across all species using three simultaneous processes

.. code-block:: bash

    find ../../SPECIES -mindepth 1 -maxdepth 1 -type d -print0 \
    | xargs -0 -n1 -P 3 sh -c '
        d="$1"
        species=$(basename "$d")

        echo "START: $species"

        if python3 scripts/02_find_ZKSCAN7_chimeric_signature.py "$species" \
            > "results/batch_logs/${species}.log" 2>&1
        then
            echo "OK:    $species"
        else
            echo "FAILED: $species"
            touch "results/batch_logs/${species}.failed"
        fi
    ' sh

