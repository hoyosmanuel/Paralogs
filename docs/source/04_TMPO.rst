I. IDENTIFICAR LOS PARALOGOS DE TMPO
=====================================

1) Script to extract TMPOs from TOGA
------------------------------------

Vamoa a usar tres evidencias el exon1, el exon 2 y el TE:

TMPO exon 1
TMPO exon 2
TE localizado hacia el 3' biológico del exon 2

La firma mínima va a ser:

.. code-block::

  Exon1 match
      +
  Exon2 match
      +
  SPIN_MI y/o MicTal-1.345 match
  en ≤5 kb hacia el 3' biológico del Exon2
      ↓
  TMPO_PARALOG_SIGNATURE_FOUND

Entramos al directorio

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO

Crear las sondas exónicas de TMPO

.. code-block:: bash

	cat > sequences/TMPO_exon_probes.fa <<'EOF'
	>TMPO_exon1_match
	ATGGCGGAGTTCCTGCAAGACCCCTGGGTCCTGACGACGGAGAAGCTGAAGAGGGAGCTGGTGGCCCACAATGTGCCGCTGCCGGCGGGGCAGCCGCGCAAAGACGTGTACGTGCAGCTCTACCTGCAGCACCTCACGGCGCGCAACCGGCCGCCGCTCGCCGCCGGCGCCGACAGCCAAGCGCCCCCGGACTTGTCCGGCGACGAGGAGCGCAAGCCCAGCCCCGTCCTCACCGCCCCGGTCCTCGCCGCCCCCGTCCTCGCCGCCGCGGGCCGCAGCCGCAGG
	>TMPO_exon2_match
	AAAGCCACAAAGAAAACTGATAAACCCAGACTAGAAGATAAAGATGATTTAGATGTAACTGAACTCACTAATGAAGCTCTTCTGAATCAGCTTGTGAAATATGGGGTGAATCCTGGCCCTATTGTGG
	>TMPO_exon3_match
	AGAGGAAGGGAGAGGGATAGAGAGTTAGAAACATCGATGAGAGAGATCGATCAGCTGCCTCCTGCACGCGCCCCAGTGGGGATGTGCCCGCAACCAAG
	EOF

Comprobamos:

.. code-block:: bash

  grep '^>' sequences/TMPO_exon_probes.fa

Crear las secuencias de los TEs

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO

  cat > sequences/MicTal-1.345_probe.fa <<'EOF'
  >MicTal-1.345_match
  TATGTTGTGAAATATTTATTTCCAATTACAAATAAATGAAATTCTGCCTTGAAGCATGGTGTAGCATGGGTAACAGTCTTAATATAACAACAGTAAACTACATTGCTACATTTTGGGACAGTATGCGTAAAAAGCCAACGTCAGTGCGAAGGTTGAGATTGCTTCTTTCTCACCAGATCAGAAATTCTCGGGGCAGTCTTTGCAAGAGCACAGCGAAGATCATCTTCCACAACAAGTCTCTTTCTGTATTTAGACTTGACAACCAACAAGCTGGAAAACCCTGTTTCACAAAGATGTGTTGTTGGAAATGGAAGAAGAAGGCGCAACACAGTCTCAGACAGAACAGGATACGACTGCAAATACTTGATCCAGAATTTTGTAACTGGCATTGTAGAGAAATCAGTTCTTGCTGCATCATTGTTAGTAAGTTCAATGAACTCCTCTTGTGCTGTCTCAGGAACATCTTCAACTTTGACTGTGAATGGGCTTCTAACAAGTGCAAATGGGGTGTCAGGTAGATTTGGAAAGTACCATGAAATTTCGCTTGCAAGCATGTGCAAGTGTTCTTTCACAGAAATTATCATCGTAATCCTTTTGTGTGGTTCAATGTCATGTTCCTCAAAGAAAGCAGATAAGGTAGGCAACATGTAAATTTTATTTTCATCCAATTTTTTTTGCCAAAGCTGAAGTTTCATTTGGAATGATTGGATCTTTTCAGATAAATCCAGACATGTTGCATTTGGTCCTTGCAGTGATAAATTTAATTCATGCAAGATGGCAAAGATGTCAACCAAATAAGCTATTTTCTGCAGTTCACTTTTATCGCTGAAAAGTGCTTCGAACTGTGGTTTTGCTTTCTGATTAAAAAATGTTTTGAGTTCATCACGAAGCTCAAAAACACGTTTTAAAATTTTTCCTCTCGACAACCATCTCGCTTCAGTGTGAAATAGCAGAGCATTGTTTGGCGCATCCAACTCATTGCACAGATGTGAAAACAGTCGACTGTTTAAAGTGCTCGCCTTTATAAAATTGACACAACTTATGACGCTTTTCATTACTTCTTGTAACTCTTGTGGCAGCGTCTTCCTTGCTAATATTTGCCCATGAATCATACAGTGAGTTCCGATGACTTTTGGTGACTCATTCAGTACCAAGCGTTGAAATCCAGATAGGCATCCTAGCGTAGCTGGAGCACCATCTGTGCAAACACCACAAACCTTTTCCCAAGAGATCTTATGCTCTTTCAGGAAAGAACCAACTGTGTCAAATACATCAAGGGCAGTAGCTGTCGTTTTAAGAGGTTTACAAAAAAGAAACTCATCCTTAAAGTCGCCATCATTAATATACCTCACGTAAACCAGTAACTGTGAACAGTTTGCAACGTCTGTAGATTCATCAAGCTGGATACTAAATATTGGAAGCGGAGCAGATTTAATTTCCTGGATTACCTGATTAAGAATATCAGCAGACATGTCATGTATTCTTCTGCGGATAGTGTCGTTAGATAAGGAAATTGCACTCAATTTCGTAACAAATTCGTCTCCGATCATAACTCGAACAATGTCTTTGGCCTCTGGCAACAGTAAATCCTTAGCAATGGTGTGAGGTTTCATAGCTCTGGCGATTCTGATTGCCACCAAATATGAAGCTTCAACGGCTGCTACGTTTTCTTTTTGGTACTTTCCACCAGTGTCAAGTCTGGCTTTCTTGAGTCCATCAGCTTTGCTTCTGAAATAGCTGGTATCCTTGCCGGCAAAGCTTGGATGCTTGCTATCAAAATGGCGTTTTAGTTTGTTCGGCTTCATAGATTCGGATGATAGAACGTCACAACAAATAACACATTGAGGTTTCTCAATTCCTTCTGTAATTATTGATGTGAAACCATACTGTAAAAAATCCTCACTATATTTTCTTTTCTTAACGTTCTTCTCTACACGGTCCATTGTCATGGGATGGTTTG
  EOF

  cat > sequences/SPIN_MI_probe.fa <<'EOF'
  >SPIN_MI_match
  CTGAAGAACGCACTAAAACAGCGATTCTCAATCTGTGGGTTGCGACCTTTTAGGGGTCAAACGACCCTCTCACAGGGGTCGCCTAAGACCATCGGAAAACACATATTTTCAGTAGTCTTAGGAACCCAGACACCGTTCCTCTATCCGTCTCCAGGCGGGTCTGCCCACATGCAGATACATGAGTACCTGGCGTGATGACATCATCGCACCAACCCCGTCACATACACCCTGTACAAATACAGGTGTATGTGACAGGGTTGGCGCCATAATGTACTTATGTGGACCAGTTACACATGTGTAGAGAGCAGCTACTCTGTAGAAAGCAGTTGCTGTGTTGAAAGCAGCTACTCTGTTAAAAGCAGCTACTGTGTTGAAAGCAGCTACTGTGTTGAAAGCAGCAGTATTGGAGGTAAAACGACACTTCGTGAATTATAATTACTGGGTAAATGTAAAATCATATAATGTAAAATCAACTACTGCAAAAAAACCACATCTGTACCATGGGAACTTAATCTGGATGCTGATCGGTCTTTTTATATTCAGCTGTGGTTGATGTGAATACTGCCCCCTTGTGATAGTAACAGGTATGTAAAAAACAACAACACAGAATGGTAAATCATAAGAAACTTTAAATGAACTGTATTGATTGAACTATGCCATGTATCATCTTTTGTATTATTAAAGCTATTGTTATATATTATTTTCATTAGCAAACCATCCCATGAGAATGGACCGTGTAGAGAAGAACGTTAAGAAAAGAAAATATAGTGAGGATTTTTTACAGTATGGTTTCACATCAATAATTACAGAAGGAATTGAGAAACCTCAATGTGTTATTTGTTGTGACGTTCTATCATCCGAATCTATGAAGCCGAACAAACTAAAACGCCATTTTGATAGCAAGCATCCGAGCTTTGCCGGCAAGGATACCAACTATTTTAGAAGCAAAGCTGATGGACTCAAGAAAGCCAGACTTGACACTGGTGGAAAGTACCAAAAAGAAAACGTAGCAGCCGTTGAAGCTTCATATTTGGTGGCACTCAGAATCGCCAGAGCTATGAAACCTCACACCATTGCTAAGGATTTACTGTTGCCAGAGGCCAAAGACATTGTTCGAGTTATGATCGGAGATGAATTTGTTACGAAATTGAGTGCAATTTCCTTATCTAATGACACTATCCGCAGAAGAATACATGACATGTCTGCTGATATTCTTAATCAGGTAATCCAGGAAATTAAATCTGCTCCACTTCCAATATTTAGTATCCAGCTTGATGAATCTACAGACGTTGCAAACTGTTCACAGTTACTGGTTTACGTAAGGTATATTAATGATGGCGACTTTAAGGATGAGTTTCTTTTTTGCAAACCCCTTAAAACGACAGCTACTGCCCTTGATGTATTTGACACAGTTGGTTCTTTTCTGAAAGAGCATAAGATCTCTTGGGAAAAGGTTTGTGGTGTTTGCACAGATGGTGCTCCAGCTACGCTAGGATGTCTATCTGGATTTCAACGTTTGGTACTGAATGAGTCACCAAAAGTCATCGGAACTCACTGTATGATTCATGGGCAAATATTAGCAAGGAAGACGCTGCCACAAGAGTTACAAGAAGTAATGAAAAGCGTCATAAGTTGTGTCAATTTTGTAAAGGCGAGCACTTTAAACAGTCGACTGTTTTCGCAACTGTGCAATGAGTTGGATGTGCCAAACAATGCTTTGCTATTTCACACTGAAGCGAGATGGTTGTCGAGAGGAAAAATTTTAAAACGTGTTTTTGAGCTTCGTGATGAACTCAAAACATTTTTTAATCAGAAAGCAAAACCACAGTTCGAAGCACTTTTCAGCGATAAAAGTGAACTGCAGAAAATAGCTTACTTGGTTGACATCTTTGCCATCTTGCATGAGTTAAATTTATCACTGCAAGGACCAAATGCAACATGCCTGGATTTGTCTGAAAAGATCCAATCATTCCAAATGAAACTTCAGCTTTGGCAAAAAAAATTGGATGAAAATAAAATTTACATGTTGCCTACCTTATCTGCTTTCTTTGAGGAACATGACATTGAACCACACAAAAGGATTACGATGATAATTTCTGTGAAAGAACACTTGCACATGCTTGCAAGCGAAATTTCATGGTACTTTCCAAATCTACCTGACACCCCATTTGCACTTGTTAGAAGCCCATTCACAGTCAAAGTTGAAGATGTTCCTGAGAAAGCACAAGAGGAGTTCATTGAACTTACTAACAACGATGCAGCGAGAACTGATTTCTCTACAATGCCCGTTACAAAATTCTGGATCAAGTATTTGCAGTCATATCCTGTTCTGTCCGAGACTGTGTTGCGCCTTCTTCTTCCATTTCCAACAACATATCTTTGTGAAACAGGGTTTTCCAGCTTGTTGGTTATCAAGTCTAAATACAGAAAGAGACTTGTTGTGGAAGATGATCTTCGTTGTGCTCTTGCAAAGACTGCCCCGAGAATTTCTGATCTGGTGAGAAAGAAGCAATCTCAACCTTCGCACTGACGTTGGCTTTTTACGCATACTGTCACAAAATGTAGCAATGTAGTTTACTGTTGTTATATTAAGACTGTTACCCATGCTTCAAGACAAAATTTCATTTATTTGTAATTAGAAATAAATATTTCACAATATATAATTACATATTGTTTTTGTGATTAATCACTGTGCTTTAATTATGTTCAATTTTTAACAGTGACAATACATTCTGCATATCAGATATTTACGATTCATAACAGTAGCAAAATTAAGTTATGAAGTAGCAACGAAAATAATTTTATGGTTGGGGGTGACCACAACATGAGGAACTGTGTTAAAGGGTCGCGGCATGAGGAAGGCTGAGAACCACTGCCTTAA AATTTATAAGC
  EOF

  #Borrar cualquier posible espacio
  sed -i '' '/^>/! s/[[:space:]]//g' sequences/SPIN_MI_probe.fa
  sed -i '' '/^>/! s/[[:space:]]//g' sequences/MicTal-1.345_probe.fa


Comprobar las longitudes:

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
  ' sequences/SPIN_MI_probe.fa sequences/MicTal-1.345_probe.fa

Juntar las dos secuencias:

.. code-block:: bash

  cat sequences/SPIN_MI_probe.fa \
      sequences/MicTal-1.345_probe.fa \
      > sequences/TMPO_TE_probes.fa

hacer el script:

.. code-block:: bash

	cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO

.. code-block:: bash


	#!/usr/bin/env python3
	
	from pathlib import Path
	import argparse
	import csv
	import subprocess
	import sys
	
	
	# ============================================================
	# TMPO PARALOG SIGNATURE FINDER
	# ============================================================
	#
	# Genome-wide strategy:
	#
	#   TMPO exon 1 match
	#          +
	#   TMPO exon 2 match
	#          +
	#   TE match <= 5 kb toward the biological 3' end of exon 2
	#
	# TE probes:
	#
	#   SPIN_MI_match
	#   MicTal-1.345_match
	#
	# Exon 3 is NOT used in the genome-wide search because its
	# short / low-complexity sequence produces excessive matches.
	#
	# IMPORTANT:
	#
	# A BLAST match to a TE probe identifies a candidate TE region.
	# Final TE-family identity should subsequently be confirmed
	# using the RepeatMasker annotation overlapping that region.
	# ============================================================
	
	
	# ============================================================
	# PATHS
	# ============================================================
	
	ROOT = Path(__file__).resolve().parent.parent
	PROJECT = ROOT.parent.parent
	
	SPECIES_ROOT = PROJECT / "SPECIES"
	
	EXON_PROBES = (
	    ROOT /
	    "sequences" /
	    "TMPO_exon_probes.fa"
	)
	
	TE_PROBES = (
	    ROOT /
	    "sequences" /
	    "TMPO_TE_probes.fa"
	)
	
	DB_DIR = (
	    ROOT /
	    "intermediate" /
	    "blast_db"
	)
	
	
	# ============================================================
	# PARAMETERS
	# ============================================================
	
	# Exon 1 and Exon 2 must occur within this genomic distance.
	MAX_EXON_SPAN = 50000
	
	# TE must occur within this distance toward the biological
	# 3' end of Exon 2.
	MAX_TE_DISTANCE = 5000
	
	# BLAST filters for TMPO exons.
	MIN_EXON_IDENTITY = 80.0
	MIN_EXON_QCOV = 80.0
	
	# BLAST filters for TE probes.
	MIN_TE_IDENTITY = 80.0
	MIN_TE_QCOV = 80.0
	
	
	# ============================================================
	# ARGUMENTS
	# ============================================================
	
	parser = argparse.ArgumentParser(
	    description=(
	        "Find TMPO-like loci and identify candidate TMPO "
	        "paralog signatures using Exon1 + Exon2 + 3-prime TE."
	    )
	)
	
	parser.add_argument(
	    "species",
	    help="Species directory name, e.g. Myotis_daubentonii"
	)
	
	args = parser.parse_args()
	
	species = args.species
	species_dir = SPECIES_ROOT / species
	
	
	# ============================================================
	# BASIC INPUT CHECKS
	# ============================================================
	
	if not species_dir.exists():
	    sys.exit(
	        f"\nERROR: species directory does not exist:\n"
	        f"{species_dir}\n"
	    )
	
	if not EXON_PROBES.exists():
	    sys.exit(
	        f"\nERROR: exon probe FASTA does not exist:\n"
	        f"{EXON_PROBES}\n"
	    )
	
	if not TE_PROBES.exists():
	    sys.exit(
	        f"\nERROR: TE probe FASTA does not exist:\n"
	        f"{TE_PROBES}\n"
	    )
	
	
	# ============================================================
	# FASTA READER
	# ============================================================
	
	def read_fasta(path):
	
	    records = {}
	    header = None
	    sequence = []
	
	    with path.open() as fh:
	
	        for line in fh:
	
	            line = line.strip()
	
	            if not line:
	                continue
	
	            if line.startswith(">"):
	
	                if header is not None:
	                    records[header] = "".join(sequence).upper()
	
	                header = line[1:].split()[0]
	                sequence = []
	
	            else:
	                sequence.append(line)
	
	    if header is not None:
	        records[header] = "".join(sequence).upper()
	
	    return records
	
	
	# ============================================================
	# CHECK PROBE CONTENT
	# ============================================================
	
	exon_records = read_fasta(EXON_PROBES)
	te_records = read_fasta(TE_PROBES)
	
	
	required_exons = {
	    "TMPO_exon1_match",
	    "TMPO_exon2_match",
	}
	
	required_tes = {
	    "SPIN_MI_match",
	    "MicTal-1.345_match",
	}
	
	
	missing_exons = (
	    required_exons -
	    set(exon_records)
	)
	
	missing_tes = (
	    required_tes -
	    set(te_records)
	)
	
	
	if missing_exons:
	
	    sys.exit(
	        "\nERROR: missing required exon probes:\n  "
	        +
	        "\n  ".join(sorted(missing_exons))
	        +
	        "\n"
	    )
	
	
	if missing_tes:
	
	    sys.exit(
	        "\nERROR: missing required TE probes:\n  "
	        +
	        "\n  ".join(sorted(missing_tes))
	        +
	        "\n"
	    )
	
	
	# Reject accidental characters in probe sequences.
	for name, sequence in {
	    **exon_records,
	    **te_records
	}.items():
	
	    bad = set(sequence) - set("ACGTN")
	
	    if bad:
	
	        sys.exit(
	            f"\nERROR: probe {name} contains invalid "
	            f"characters: {sorted(bad)}\n"
	        )
	
	
	# ============================================================
	# FIND GENOME FASTA
	# ============================================================
	
	fastas = [
	    p for p in species_dir.glob("*.fa")
	    if not p.name.startswith("._")
	]
	
	
	if len(fastas) == 0:
	
	    sys.exit(
	        f"\nERROR: no uncompressed genome .fa found in:\n"
	        f"{species_dir}\n"
	    )
	
	
	if len(fastas) > 1:
	
	    print(
	        "\nERROR: multiple uncompressed genome FASTAs found:\n"
	    )
	
	    for fasta in fastas:
	        print(f"  {fasta.name}")
	
	    sys.exit(
	        "\nSelect the appropriate genome manually.\n"
	    )
	
	
	genome = fastas[0]
	
	
	# ============================================================
	# WORKING PATHS
	# ============================================================
	
	species_work = (
	    ROOT /
	    "intermediate" /
	    species
	)
	
	species_work.mkdir(
	    parents=True,
	    exist_ok=True
	)
	
	DB_DIR.mkdir(
	    parents=True,
	    exist_ok=True
	)
	
	
	db_prefix = (
	    DB_DIR /
	    species
	)
	
	
	# Exon 3 is deliberately omitted here.
	core_exon_probes = (
	    species_work /
	    "TMPO_core_exon_probes.fa"
	)
	
	
	exon_blast = (
	    species_work /
	    f"TMPO_exon_matches_{species}.tsv"
	)
	
	
	te_blast = (
	    species_work /
	    f"TMPO_TE_matches_{species}.tsv"
	)
	
	
	result = (
	    ROOT /
	    "results" /
	    f"TMPO_signature_{species}.tsv"
	)
	
	
	# ============================================================
	# CREATE EXON 1 + EXON 2 FASTA
	# ============================================================
	
	with core_exon_probes.open("w") as out:
	
	    for name in [
	        "TMPO_exon1_match",
	        "TMPO_exon2_match",
	    ]:
	
	        out.write(
	            f">{name}\n"
	        )
	
	        out.write(
	            exon_records[name] +
	            "\n"
	        )
	
	
	# ============================================================
	# REPORT PROBE LENGTHS
	# ============================================================
	
	print("\nProbe lengths:")
	
	print(
	    f'  TMPO exon 1:    '
	    f'{len(exon_records["TMPO_exon1_match"])} bp'
	)
	
	print(
	    f'  TMPO exon 2:    '
	    f'{len(exon_records["TMPO_exon2_match"])} bp'
	)
	
	print(
	    f'  SPIN_MI:        '
	    f'{len(te_records["SPIN_MI_match"])} bp'
	)
	
	print(
	    f'  MicTal-1.345:   '
	    f'{len(te_records["MicTal-1.345_match"])} bp'
	)
	
	
	# ============================================================
	# BUILD BLAST DATABASE IF NECESSARY
	# ============================================================
	
	db_files = list(
	    DB_DIR.glob(
	        species + ".*"
	    )
	)
	
	
	if not db_files:
	
	    print(
	        f"\nCreating BLAST database for {species}..."
	    )
	
	    subprocess.run(
	        [
	            "makeblastdb",
	
	            "-in",
	            str(genome),
	
	            "-dbtype",
	            "nucl",
	
	            "-out",
	            str(db_prefix),
	        ],
	        check=True
	    )
	
	else:
	
	    print(
	        f"\nUsing existing BLAST database for {species}"
	    )
	
	
	# ============================================================
	# BLAST SETTINGS
	# ============================================================
	
	OUTFMT = (
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
	
	
	def run_blast(
	    query,
	    output,
	    evalue
	):
	
	    command = [
	        "blastn",
	
	        "-query",
	        str(query),
	
	        "-db",
	        str(db_prefix),
	
	        "-task",
	        "blastn",
	
	        "-word_size",
	        "7",
	
	        "-dust",
	        "no",
	
	        "-evalue",
	        str(evalue),
	
	        "-outfmt",
	        OUTFMT,
	    ]
	
	    with output.open("w") as out:
	
	        subprocess.run(
	            command,
	            stdout=out,
	            check=True
	        )
	
	
	# ============================================================
	# RUN BLAST
	# ============================================================
	
	print(
	    "\nSearching TMPO Exon 1 and Exon 2 matches..."
	)
	
	run_blast(
	    core_exon_probes,
	    exon_blast,
	    "1e-10"
	)
	
	
	print(
	    "Searching SPIN_MI and MicTal-1.345 matches..."
	)
	
	run_blast(
	    TE_PROBES,
	    te_blast,
	    "1e-20"
	)
	
	
	# ============================================================
	# READ BLAST OUTPUT
	# ============================================================
	
	def read_blast(path):
	
	    hits = []
	
	    with path.open() as fh:
	
	        reader = csv.reader(
	            fh,
	            delimiter="\t"
	        )
	
	        for row in reader:
	
	            # Our outfmt contains exactly 12 columns.
	            if len(row) != 12:
	                continue
	
	            sstart = int(row[7])
	            send = int(row[8])
	
	            hits.append(
	                {
	                    "query": row[0],
	
	                    "scaffold": row[1],
	
	                    "pident": float(row[2]),
	
	                    "length": int(row[3]),
	
	                    "qlen": int(row[4]),
	
	                    "qstart": int(row[5]),
	
	                    "qend": int(row[6]),
	
	                    "start": min(
	                        sstart,
	                        send
	                    ),
	
	                    "end": max(
	                        sstart,
	                        send
	                    ),
	
	                    "strand": (
	                        "+"
	                        if sstart < send
	                        else "-"
	                    ),
	
	                    "evalue": float(row[9]),
	
	                    "bitscore": float(row[10]),
	
	                    "qcov": float(row[11]),
	                }
	            )
	
	    return hits
	
	
	exon_hits = read_blast(
	    exon_blast
	)
	
	te_hits = read_blast(
	    te_blast
	)
	
	
	# ============================================================
	# FILTER EXON MATCHES
	# ============================================================
	
	exon1_hits = [
	
	    h for h in exon_hits
	
	    if (
	        h["query"]
	        ==
	        "TMPO_exon1_match"
	
	        and
	
	        h["pident"]
	        >=
	        MIN_EXON_IDENTITY
	
	        and
	
	        h["qcov"]
	        >=
	        MIN_EXON_QCOV
	    )
	]
	
	
	exon2_hits = [
	
	    h for h in exon_hits
	
	    if (
	        h["query"]
	        ==
	        "TMPO_exon2_match"
	
	        and
	
	        h["pident"]
	        >=
	        MIN_EXON_IDENTITY
	
	        and
	
	        h["qcov"]
	        >=
	        MIN_EXON_QCOV
	    )
	]
	
	
	# ============================================================
	# RECONSTRUCT TMPO-LIKE LOCI
	# ============================================================
	
	candidate_loci = []
	
	
	for exon1 in exon1_hits:
	
	    for exon2 in exon2_hits:
	
	        # Same scaffold.
	        if (
	            exon1["scaffold"]
	            !=
	            exon2["scaffold"]
	        ):
	            continue
	
	
	        # Same orientation.
	        if (
	            exon1["strand"]
	            !=
	            exon2["strand"]
	        ):
	            continue
	
	
	        strand = exon1["strand"]
	
	
	        # ----------------------------------------------------
	        # Correct exon order
	        # ----------------------------------------------------
	
	        if strand == "+":
	
	            # transcription ->
	            #
	            # Exon1 ... Exon2
	            if (
	                exon1["start"]
	                >=
	                exon2["start"]
	            ):
	                continue
	
	        else:
	
	            # transcription <-
	            #
	            # Exon2 ... Exon1
	            #
	            # Exon1 therefore has larger genomic coordinates.
	            if (
	                exon1["start"]
	                <=
	                exon2["start"]
	            ):
	                continue
	
	
	        locus_start = min(
	            exon1["start"],
	            exon2["start"]
	        )
	
	        locus_end = max(
	            exon1["end"],
	            exon2["end"]
	        )
	
	        locus_span = (
	            locus_end
	            -
	            locus_start
	            +
	            1
	        )
	
	
	        if (
	            locus_span
	            >
	            MAX_EXON_SPAN
	        ):
	            continue
	
	
	        candidate_loci.append(
	            {
	                "scaffold":
	                    exon1["scaffold"],
	
	                "strand":
	                    strand,
	
	                "start":
	                    locus_start,
	
	                "end":
	                    locus_end,
	
	                "span":
	                    locus_span,
	
	                "exon1":
	                    exon1,
	
	                "exon2":
	                    exon2,
	
	                "pair_score":
	                    (
	                        exon1["bitscore"]
	                        +
	                        exon2["bitscore"]
	                    ),
	            }
	        )
	
	
	# ============================================================
	# COLLAPSE DUPLICATE LOCI
	# ============================================================
	#
	# Keep the best Exon1 + Exon2 combination per scaffold/strand.
	# ============================================================
	
	best_loci = {}
	
	
	for locus in candidate_loci:
	
	    key = (
	        locus["scaffold"],
	        locus["strand"],
	    )
	
	
	    if key not in best_loci:
	
	        best_loci[key] = locus
	
	    else:
	
	        if (
	            locus["pair_score"]
	            >
	            best_loci[key]["pair_score"]
	        ):
	
	            best_loci[key] = locus
	
	
	loci = list(
	    best_loci.values()
	)
	
	
	# ============================================================
	# FILTER TE MATCHES
	# ============================================================
	
	te_hits = [
	
	    h for h in te_hits
	
	    if (
	        h["pident"]
	        >=
	        MIN_TE_IDENTITY
	
	        and
	
	        h["qcov"]
	        >=
	        MIN_TE_QCOV
	    )
	]
	
	
	# ============================================================
	# FIND TE MATCHES AT BIOLOGICAL 3' OF EXON 2
	# ============================================================
	
	results = []
	
	
	for locus in loci:
	
	    exon2 = locus["exon2"]
	
	    valid_te_matches = []
	
	
	    for te in te_hits:
	
	        if (
	            te["scaffold"]
	            !=
	            locus["scaffold"]
	        ):
	            continue
	
	
	        # ----------------------------------------------------
	        # PLUS STRAND
	        #
	        # Exon1 -> Exon2 -> TE
	        # ----------------------------------------------------
	
	        if locus["strand"] == "+":
	
	            if (
	                te["start"]
	                <
	                exon2["end"]
	            ):
	                continue
	
	            distance = (
	                te["start"]
	                -
	                exon2["end"]
	            )
	
	
	        # ----------------------------------------------------
	        # MINUS STRAND
	        #
	        # TE <- Exon2 <- Exon1
	        # ----------------------------------------------------
	
	        else:
	
	            if (
	                te["end"]
	                >
	                exon2["start"]
	            ):
	                continue
	
	            distance = (
	                exon2["start"]
	                -
	                te["end"]
	            )
	
	
	        if (
	            distance
	            <=
	            MAX_TE_DISTANCE
	        ):
	
	            valid_te_matches.append(
	                {
	                    **te,
	
	                    "distance":
	                        distance
	                }
	            )
	
	
	    # ========================================================
	    # KEEP BEST MATCH FOR EACH TE PROBE
	    # ========================================================
	
	    best_by_probe = {}
	
	
	    for te in valid_te_matches:
	
	        probe = te["query"]
	
	
	        if probe not in best_by_probe:
	
	            best_by_probe[probe] = te
	
	            continue
	
	
	        old = best_by_probe[probe]
	
	
	        new_rank = (
	            te["distance"],
	            -te["qcov"],
	            -te["bitscore"],
	        )
	
	        old_rank = (
	            old["distance"],
	            -old["qcov"],
	            -old["bitscore"],
	        )
	
	
	        if new_rank < old_rank:
	
	            best_by_probe[probe] = te
	
	
	    selected_te_matches = list(
	        best_by_probe.values()
	    )
	
	
	    if selected_te_matches:
	
	        status = (
	            "TMPO_PARALOG_SIGNATURE_FOUND"
	        )
	
	    else:
	
	        status = (
	            "TMPO_LIKE_LOCUS_NO_3PRIME_TE"
	        )
	
	
	    results.append(
	        {
	            **locus,
	
	            "te_matches":
	                selected_te_matches,
	
	            "status":
	                status,
	        }
	    )
	
	
	# ============================================================
	# WRITE RESULT TABLE
	# ============================================================
	
	columns = [
	
	    "Species",
	
	    "Genome",
	
	    "Scaffold",
	
	    "Locus_start",
	
	    "Locus_end",
	
	    "Strand",
	
	    "Exon1_start",
	
	    "Exon1_end",
	
	    "Exon1_identity",
	
	    "Exon1_qcov",
	
	    "Exon2_start",
	
	    "Exon2_end",
	
	    "Exon2_identity",
	
	    "Exon2_qcov",
	
	    "TE_probe_matches",
	
	    "Status",
	]
	
	
	with result.open("w") as out:
	
	    writer = csv.writer(
	        out,
	        delimiter="\t"
	    )
	
	    writer.writerow(
	        columns
	    )
	
	
	    for r in sorted(
	        results,
	
	        key=lambda x: (
	            x["status"]
	            !=
	            "TMPO_PARALOG_SIGNATURE_FOUND",
	
	            x["scaffold"]
	        )
	    ):
	
	        e1 = r["exon1"]
	        e2 = r["exon2"]
	
	
	        summaries = []
	
	
	        for te in sorted(
	            r["te_matches"],
	            key=lambda x: x["query"]
	        ):
	
	            summaries.append(
	                (
	                    f'{te["query"]}:'
	                    f'{te["start"]}-{te["end"]}:'
	                    f'id={te["pident"]:.3f}:'
	                    f'cov={te["qcov"]:.1f}:'
	                    f'distance={te["distance"]}'
	                )
	            )
	
	
	        writer.writerow(
	            [
	                species,
	
	                genome.name,
	
	                r["scaffold"],
	
	                r["start"],
	
	                r["end"],
	
	                r["strand"],
	
	
	                e1["start"],
	
	                e1["end"],
	
	                f'{e1["pident"]:.3f}',
	
	                f'{e1["qcov"]:.1f}',
	
	
	                e2["start"],
	
	                e2["end"],
	
	                f'{e2["pident"]:.3f}',
	
	                f'{e2["qcov"]:.1f}',
	
	
	                (
	                    ";".join(summaries)
	                    if summaries
	                    else "."
	                ),
	
	                r["status"],
	            ]
	        )
	
	
	# ============================================================
	# HUMAN-READABLE REPORT
	# ============================================================
	
	print("\n")
	print("=" * 64)
	print(species)
	print("=" * 64)
	
	
	for r in sorted(
	    results,
	
	    key=lambda x: (
	        x["status"]
	        !=
	        "TMPO_PARALOG_SIGNATURE_FOUND",
	
	        x["scaffold"]
	    )
	):
	
	    e1 = r["exon1"]
	    e2 = r["exon2"]
	
	
	    if (
	        r["status"]
	        ==
	        "TMPO_PARALOG_SIGNATURE_FOUND"
	    ):
	
	        label = "TMPO PARALOG"
	
	    else:
	
	        label = "TMPO-LIKE LOCUS"
	
	
	    print(
	        f"\n[{label}]"
	    )
	
	    print(
	        f'Scaffold: {r["scaffold"]}'
	    )
	
	    print(
	        f'Strand:   {r["strand"]}'
	    )
	
	    print()
	
	
	    print(
	        f'Exon1    '
	        f'{e1["start"]}-{e1["end"]}    '
	        f'id={e1["pident"]:.3f}%   '
	        f'cov={e1["qcov"]:.0f}%'
	    )
	
	
	    print(
	        f'Exon2    '
	        f'{e2["start"]}-{e2["end"]}    '
	        f'id={e2["pident"]:.3f}%   '
	        f'cov={e2["qcov"]:.0f}%'
	    )
	
	
	    if r["te_matches"]:
	
	        print(
	            f"\nTE matches within "
	            f"{MAX_TE_DISTANCE} bp "
	            "of Exon2 biological 3':"
	        )
	
	
	        for te in sorted(
	            r["te_matches"],
	            key=lambda x: (
	                x["distance"],
	                x["query"]
	            )
	        ):
	
	            print()
	
	            print(
	                f'  {te["query"]}'
	            )
	
	            print(
	                f'    coordinates: '
	                f'{te["start"]}-{te["end"]}'
	            )
	
	            print(
	                f'    identity:    '
	                f'{te["pident"]:.3f}%'
	            )
	
	            print(
	                f'    coverage:    '
	                f'{te["qcov"]:.0f}%'
	            )
	
	            print(
	                f'    distance:    '
	                f'{te["distance"]} bp'
	            )
	
	
	    else:
	
	        print(
	            f"\nTE match: NOT FOUND within "
	            f"{MAX_TE_DISTANCE} bp "
	            "of Exon2 biological 3'"
	        )
	
	
	    print()
	
	    print(
	        f'-> {r["status"]}'
	    )
	
	
	# ============================================================
	# SUMMARY
	# ============================================================
	
	positive = [
	
	    r for r in results
	
	    if (
	        r["status"]
	        ==
	        "TMPO_PARALOG_SIGNATURE_FOUND"
	    )
	]
	
	
	print(
	    "\n" +
	    "=" * 64
	)
	
	print(
	    f"TMPO paralog signatures detected: "
	    f"{len(positive)}"
	)
	
	print(
	    "=" * 64
	)
	
	
	if positive:
	
	    print(
	        "\nPositive loci:"
	    )
	
	    for r in positive:
	
	        print(
	            f'  {r["scaffold"]}:'
	            f'{r["start"]}-'
	            f'{r["end"]} '
	            f'({r["strand"]})'
	        )
	
	
	print(
	    f"\nResult table:\n{result}"
	)
	
	print(
	    f"\nRaw Exon1/Exon2 matches:\n{exon_blast}"
	)
	
	print(
	    f"\nRaw TE-probe matches:\n{te_blast}"
	)

.. code-block:: bash

	python3 scripts/02_find_TMPO_paralog_signature.py Myotis_daubentonii


el resultado es este:

.. code-block:: bash

	================================================================
	Myotis_daubentonii
	================================================================
	
	[TMPO PARALOG]
	Scaffold: NC_081859
	Strand:   -
	
	Exon1    11870994-11871278    id=99.649%   cov=100%
	Exon2    11859775-11859901    id=99.213%   cov=100%
	
	TE matches within 5000 bp of Exon2 biological 3':
	
	  SPIN_MI_match
	    coordinates: 11856166-11859021
	    identity:    94.443%
	    coverage:    100%
	    distance:    754 bp
	
	  MicTal-1.345_match
	    coordinates: 11856407-11858366
	    identity:    98.418%
	    coverage:    100%
	    distance:    1409 bp
	
	-> TMPO_PARALOG_SIGNATURE_FOUND
	
	[TMPO-LIKE LOCUS]
	Scaffold: NC_081841
	Strand:   -
	
	Exon1    17767416-17767670    id=89.804%   cov=89%
	Exon2    17753160-17753286    id=96.850%   cov=100%
	
	TE match: NOT FOUND within 5000 bp of Exon2 biological 3'
	
	-> TMPO_LIKE_LOCUS_NO_3PRIME_TE
	
	================================================================
	TMPO paralog signatures detected: 1
	================================================================
	
	Positive loci:
	  NC_081859:11859775-11871278 (-)
	
	Result table:
	/Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO/results/TMPO_signature_Myotis_daubentonii.tsv
	
	Raw Exon1/Exon2 matches:
	/Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO/intermediate/Myotis_daubentonii/TMPO_exon_matches_Myotis_daubentonii.tsv
	
	Raw TE-probe matches:
	/Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO/intermediate/Myotis_daubentonii/TMPO_TE_matches_Myotis_daubentonii.tsv


luego hacemos esto para ver qué coincide con los TEs en la secuencia:

.. code-block:: bash

	awk '$1=="NC_081859" && $2 < 11859021 && $3 > 11856166' \
	../../SPECIES/Myotis_daubentonii/mMyoDau2.1.pri.bed
	


2) Vamos a hacerlo para todas las especies:
-------------------------------------------

.. code-block:: bash

	cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/TMPO

.. code-block:: bash

	mkdir -p results/batch_logs

.. code-block:: bash

	find ../../SPECIES -mindepth 1 -maxdepth 1 -type d -print0 \
	| xargs -0 -n1 -P 3 sh -c '
	    d="$1"
	    species=$(basename "$d")
	
	    echo "START: $species"
	
	    if python3 scripts/02_find_TMPO_paralog_signature.py "$species" \
	        > "results/batch_logs/${species}.log" 2>&1
	    then
	        echo "OK:    $species"
	    else
	        echo "FAILED: $species"
	        touch "results/batch_logs/${species}.failed"
	    fi
	' sh


Generemos un documento:

.. code-block:: bash

	echo "TSV generated:"
	find results -maxdepth 1 -name 'TMPO_signature_*.tsv' | wc -l
	
	echo "Failed species:"
	find results/batch_logs -name '*.failed' | wc -l
	
	echo "Failed names:"
	find results/batch_logs -name '*.failed' -exec basename {} .failed \;




3) Resultados:
-------------------------------------------

.. code-block:: bash

	grep -h -v '^Species' results/TMPO_signature_*.tsv \
	| awk -F'\t' '{count[$NF]++} END {for (x in count) print x, count[x]}'

Identifcar las positivas:

.. code-block:: bash

	grep -h 'TMPO_PARALOG_SIGNATURE_FOUND' \
	results/TMPO_signature_*.tsv \
	| cut -f1-6


.. code-block:: bash

	(base) manuelhoyos@manuel TMPO % grep -h -v '^Species' results/TMPO_signature_*.tsv \
	| awk -F'\t' '{count[$NF]++} END {for (x in count) print x, count[x]}'
	 98O_LIKE_LOCUS_NO_3PRIME_TE
	 15O_PARALOG_SIGNATURE_FOUND
	(base) manuelhoyos@manuel TMPO % grep -h 'TMPO_PARALOG_SIGNATURE_FOUND' \
	results/TMPO_signature_*.tsv \
	| cut -f1-6
	Myotis_auriculus	mMyoAur1.1.pri.fa	SUPER__22	6541114	6552997	+
	Myotis_californicus	mMyoCal1.2.pri.fa	SUPER__22	6253904	6264499	+
	Myotis_daubentonii	mMyoDau2.1.pri.fa	NC_081859	11859775	11871278	-
	Myotis_evotis	mMyoEvo1.2.pri.fa	SUPER__22	11742535	11753529	-
	Myotis_lucifugus	mMyoLuc2.1.pri.fa	SUPER__22	6072092	6082862	+
	Myotis_myotis	mMyoMyo1.6.pri.fa	scaffold_m19_p_21	11543234	11553690	-
	Myotis_mystacinus	mMyoMys1.1.hap1.fa	OZ065572	6337773	6351730	+
	Myotis_nigricans	mMyoNig1.1.pri.fa	manual_scaffold_21	5415907	5426788	+
	Myotis_occultus	mMyoOcc1.1.pri.fa	SUPER__21	7575053	7585395	+
	Myotis_pilosus	mMyoPil1.2.pri.fa	LG21	6036664	6045937	+
	Myotis_thysanodes	mMyoThy1.2.pri.fa	SUPER__22	7519122	7530375	+
	Myotis_velifer	mMyoVel1.1.pri.fa	SUPER__22	11603049	11614185	-
	Myotis_vivesi	mMyoViv1.2.pri.fa	SUPER_20	11833187	11847488	-
	Myotis_volans	mMyoVol1.2.pri.fa	SUPER__22	11802255	11810690	-
	Myotis_yumanensis	mMyoYum1.1.hap1.fa	JAPQVU010000023	11610152	11618239	-
	(base) manuelhoyos@manuel TMPO %
	
