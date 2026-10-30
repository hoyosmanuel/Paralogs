I. IDENTIFY THE AP1AR GENOMIC ARCHITECTURE
==========================================

1. Create the AP1AR analysis directory
--------------------------------------

Create an independent working directory for the AP1AR analysis.

.. code-block:: bash

   cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
   mkdir -p ANALYSES/AP1AR/{intermediate,results,scripts,sequences}
   cd ANALYSES/AP1AR


2. Create the AP1AR host-exon probes
------------------------------------

Exons 1–9 are used as host-derived AP1AR probes. Exon 10 is analyzed
separately because it corresponds to the TE-derived portion of the
candidate transcript.

Create the host-probe FASTA.

.. code-block:: bash

   cat > sequences/AP1AR_host_exon_probes.fa <<'EOF'
   >AP1AR_exon1_match
   GGCGGCGGAGCGGCCCCGCGGCGCGACTGAGGCCGAGGCGAGCACGGCTCCCTGCGTGCGGCTCAGCGGCCCGGAGCCCGCATTCCGGCGCCCCTTCCTCGCGCGGCGGCGACTGGGCGGCCGCGGGAGTAGCGCCGCCTTTGCCCAGAGTCCCGCGTGGCCCCGGCGCGAATGCGGCTGCGGCCGGGCCTGCAGCGGGCGCGGGGGCGCCGGGCGGCCGCCATGGGGAACTGCTGCTGGGCGCAGTGCTTCGGGCGGATCCGCAGGGAGACGGGGCGGCTGCCGCGTGGCGGAGG
   >AP1AR_exon2_match
   ATCCAAATATTTTAGAACATGCTCAAGAAGTGAGCATTTAACGATAGAG
   >AP1AR_exon3_match
   TTTGAGAATCTAGTAGAAAGTGATGAA
   >AP1AR_exon4_match
   GGGGAAAGCCCAGGAAGCAGTCACAG
   >AP1AR_exon5_match
   GCCTCTTACTGAGGAAGAAATAATTGACCTAAGAGAGGCATTATGATTCTATTGCTGAGAAACAGAAAGATCTTGATATGAAAATTCAAAAAGAG
   >AP1AR_exon6_match
   TTAGCCTTACAAGAAGAGAAGCTAAGACTAGAAGAGGAGGCTGTGTACGCTGCCCAGCGCGAAGCAGCCAGGGCAGCAAAGCAGCGAAAGCTCTTGGAG
   >AP1AR_exon7_match
   CAGGAGAGGCTGGGCGCCGGTCCTGGATGCCGCCCCTCCAGCGGAGCCGCATGCCACAG
   >AP1AR_exon8_match
   CTCAGGACCAGAAGACGACTTCGAATCTTGTTTGCGAAACATAAAGTCACAGTATGAAGTTTTTCAAAGTAGTA
   >AP1AR_exon9_match
   GACTCTCGTCAGATGCCACCGTTTTGACACCAAATACAGAAAGCAGTTGTGACTTAATGACCAAAACCAAATCGACGAGTGGAAATGATGACAGCACGTCCTTGGATCTGGAGTGGGAAGATGAGGAAG
   EOF


3. Create the AP1AR TE-derived exon probe
-----------------------------------------

Create exon 10 as an independent TE probe.

.. code-block:: bash

   cat > sequences/AP1AR_Mariner3_pKuh_probe.fa <<'EOF'
   >Mariner3_pKuh_match
   AAATGAGCGCAAAAAGAAAAAGCTACTCCGTGGAGTACAAGAAAGGAATCGTGGAAGACTCTTGGGGCCAAAATCTTACTGCTTTCTGCAAAGAGATGATGTTGAGTATCCGATTGGTCCGAAAATGGCGAGCAGAGTATGGTAACCTCAGTGAACAAGTGGACAAGGGAAATGCTAAAAAGCACAAGTGTGGATCAGGTCAGCAACCATTATTTTCCGAGCTGGAAGACCTGATCTGTGAATGGGTTGGTGACAGGAGAGCAAAGGCTTTGGTTGTGAATAGGGATCAGATTCAAGAATTTGCCCTTGCAATGGCACCACAGTTTGACATAGCCCCCAAAGACTTCAAAGCATCACAACACTGGCTGGATAATTTCCTTCAGCGAAGTGAACTGTCTTTAAGAAGATCTACAACATTGTTTAGGTTGGAAGATACTCAAGTTATTAAGGGAGCACTTGCATTCAAGTCCTTTATTGATGATGCTGACTTCTCTAAATACAATCTTTCCAACATGATTGCTATGGATGAAACTGCAGTGTTTATGGGCCAATCATCTCAAATAACAATTGAACAGCGGGGTGCCTCCTTAGCGTACATTCCCTCCACTGCTTACGAAAGTGCACGTGTTTCCTGTATTTTGGCAATTCGTCTGGATGGCACTAAAATCCCACCTTTAATCATTTCTAAGGGCAAGAAGGAAAAGATTGAACGTGTTTTGGGCGTTTTTATTCTTGAAACTGAAAAAGCCTGGGCCACACAAGCAGTTATAAGGAAGTGGGTCGATTTAATGCTGCCGCTTGTTATGCGAGGGAGCCAAAGAGGTCTGCTAGTCTGGGATGCAGCCAGCACTCACCGCGCTAAAGACATGAAGACATTCCTTCACGAAAGAAAAATAGATCAGATAATGATTCCCGCAGGAATGACGGCCTATCTCCAGACTCTTGATATTGCAATAAACAAGCCAAAGACAATTTGCGCATGGAAATTAATGACTACATTGAAAATAGAATGGAGAGAAATCAGCGTGGAAACTTTGTAAAACCTAAACTGCAAGAAGTTGTGACTTGGGTGAAGAATTCATGGGAGAAAATTACTGACAGTTGCATTGCAAATGCACTGCGAGCAGGCTACCTTGACAAGAAGTACTCATTTAAAGACAATGCTATTGCTAAACATGAGAGATTTGGTCCACTGATTCTGAAAGAAATGGAGTCACAAGAAATTCACCAGGAAATTCAGGAATTGGATTGTTATGATGATGTTCCAGAAGATGATGACATGATTGCCATTGAATAAATGAGTAAATTCTTGATACTGTAGGATTCACAGTTTGTCTGGTTTGTGTTGACAACTACCGTACTACCATACCATATACAGGTAATAAATTTCCTTTATTTTCATTT
   EOF


4. Verify the probe sequences
-----------------------------

Confirm the host-probe FASTA headers.

.. code-block:: bash

   grep '^>' sequences/AP1AR_host_exon_probes.fa

The expected output is::

   >AP1AR_exon1_match
   >AP1AR_exon2_match
   >AP1AR_exon3_match
   >AP1AR_exon4_match
   >AP1AR_exon5_match
   >AP1AR_exon6_match
   >AP1AR_exon7_match
   >AP1AR_exon8_match
   >AP1AR_exon9_match

Check the host-exon lengths.

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
   ' sequences/AP1AR_host_exon_probes.fa

The expected lengths are::

   >AP1AR_exon1_match 296
   >AP1AR_exon2_match 49
   >AP1AR_exon3_match 27
   >AP1AR_exon4_match 26
   >AP1AR_exon5_match 95
   >AP1AR_exon6_match 99
   >AP1AR_exon7_match 59
   >AP1AR_exon8_match 74
   >AP1AR_exon9_match 129

Check the TE-derived probe.

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
   ' sequences/AP1AR_Mariner3_pKuh_probe.fa

The expected result is::

   >Mariner3_pKuh_match 1404


5. Validar Plecotus auritus
-----------------------------

.. code-block:: bash

   cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/AP1AR
   mkdir -p intermediate/blast_db

.. code-block:: bash

   makeblastdb \
   -in ../../SPECIES/Plecotus_auritus/*.fa \
   -dbtype nucl \
   -out intermediate/blast_db/Plecotus_auritus

.. code-block:: bash

   blastn \
   -query sequences/AP1AR_host_exon_probes.fa \
   -db intermediate/blast_db/Plecotus_auritus \
   -task blastn \
   -word_size 7 \
   -dust no \
   -evalue 1e-5 \
   -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
   > intermediate/AP1AR_exon_matches_Plecotus_auritus.tsv

.. code-block:: bash
   
   blastn \
   -query sequences/AP1AR_Mariner3_pKuh_probe.fa \
   -db intermediate/blast_db/Plecotus_auritus \
   -task blastn \
   -word_size 7 \
   -dust no \
   -evalue 1e-20 \
   -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
   > intermediate/AP1AR_Mariner3_matches_Plecotus_auritus.tsv

.. code-block:: bash

   sort -k1,1 -k11,11nr intermediate/AP1AR_exon_matches_Plecotus_auritus.tsv \
   | awk '!seen[$1]++' \
   | column -t

.. code-block:: bash

   sort -k11,11nr intermediate/AP1AR_Mariner3_matches_Plecotus_auritus.tsv \
   | head -20 \
   | column -t



6.  Lanzar contra 100 especies
------------------------------

.. code-block:: bash

   nano scripts/02_find_AP1AR_chimeric_signature.py

.. code-block:: bash

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
   
   # Exons 3 and 4 are only 27 and 26 bp long
   MIN_SHORT_IDENTITY = 85.0
   MIN_SHORT_QCOV = 90.0
   
   MIN_TE_IDENTITY = 70.0
   MIN_TE_QCOV = 50.0
   
   MIN_EXONS_PER_LOCUS = 6
   MAX_EXON_GAP = 20000
   
   # Positive Plecotus auritus control = 11,467 bp
   MAX_TE_DISTANCE = 20000
   
   NORMAL_EXONS = {1, 2, 5, 6, 7, 8, 9}
   SHORT_EXONS = {3, 4}
   
   
   # ============================================================
   # PATHS
   # ============================================================
   
   SCRIPT = Path(__file__).resolve()
   ROOT = SCRIPT.parent.parent
   PROJECT = ROOT.parent.parent
   
   SPECIES_ROOT = PROJECT / "SPECIES"
   
   HOST_PROBES = ROOT / "sequences" / "AP1AR_host_exon_probes.fa"
   TE_PROBE = ROOT / "sequences" / "AP1AR_Mariner3_pKuh_probe.fa"
   
   INTERMEDIATE = ROOT / "intermediate"
   RESULTS = ROOT / "results"
   
   
   # ============================================================
   # FUNCTIONS
   # ============================================================
   
   def run(cmd):
       subprocess.run(cmd, check=True)
   
   
   def strand_from_coordinates(sstart, send):
       return "+" if sstart <= send else "-"
   
   
   def exon_number(query):
       try:
           return int(
               query.replace("AP1AR_exon", "")
                    .replace("_match", "")
           )
       except Exception:
           return None
   
   
   def read_fasta(path):
       records = {}
       name = None
       seq = []
   
       with open(path) as fh:
           for line in fh:
               line = line.strip()
   
               if not line:
                   continue
   
               if line.startswith(">"):
                   if name is not None:
                       records[name] = "".join(seq)
   
                   name = line[1:]
                   seq = []
   
               else:
                   seq.append(line)
   
       if name is not None:
           records[name] = "".join(seq)
   
       return records
   
   
   def write_probe_subset(records, exons, output):
   
       with open(output, "w") as out:
   
           for name, seq in records.items():
   
               n = exon_number(name)
   
               if n in exons:
                   out.write(f">{name}\n{seq}\n")
   
   
   def read_blast(path):
   
       hits = []
   
       with open(path) as fh:
   
           for line in fh:
   
               if not line.strip():
                   continue
   
               f = line.rstrip().split("\t")
   
               if len(f) != 12:
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
               ) = f
   
               sstart = int(sstart)
               send = int(send)
   
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
                   "start": min(sstart, send),
                   "end": max(sstart, send),
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
           "Usage: python3 scripts/02_find_AP1AR_chimeric_signature.py "
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
   # PREPARE NORMAL AND SHORT EXON PROBES
   # ============================================================
   
   records = read_fasta(HOST_PROBES)
   
   normal_probe = species_intermediate / "AP1AR_normal_exons.fa"
   short_probe = species_intermediate / "AP1AR_short_exons.fa"
   
   write_probe_subset(records, NORMAL_EXONS, normal_probe)
   write_probe_subset(records, SHORT_EXONS, short_probe)
   
   
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
   
       print(f"Reusing AP1AR BLAST database for {species}")
   
   
   # ============================================================
   # BLAST NORMAL HOST EXONS
   # ============================================================
   
   normal_out = (
       species_intermediate /
       f"AP1AR_normal_exon_matches_{species}.tsv"
   )
   
   run([
       "blastn",
       "-query", str(normal_probe),
       "-db", db_prefix,
       "-task", "blastn",
       "-word_size", "7",
       "-dust", "no",
       "-evalue", "1e-5",
       "-outfmt",
       "6 qseqid sseqid pident length qlen qstart qend "
       "sstart send evalue bitscore qcovhsp",
       "-out", str(normal_out)
   ])
   
   
   # ============================================================
   # BLAST SHORT EXONS 3 AND 4
   # ============================================================
   
   short_out = (
       species_intermediate /
       f"AP1AR_short_exon_matches_{species}.tsv"
   )
   
   run([
       "blastn",
       "-query", str(short_probe),
       "-db", db_prefix,
       "-task", "blastn-short",
       "-dust", "no",
       "-evalue", "100",
       "-outfmt",
       "6 qseqid sseqid pident length qlen qstart qend "
       "sstart send evalue bitscore qcovhsp",
       "-out", str(short_out)
   ])
   
   
   # ============================================================
   # BLAST TE
   # ============================================================
   
   te_out = (
       species_intermediate /
       f"AP1AR_TE_matches_{species}.tsv"
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
   # FILTER HOST EXON HITS
   # ============================================================
   
   exon_hits = []
   
   
   for h in read_blast(normal_out):
   
       n = exon_number(h["query"])
   
       if n is None:
           continue
   
       if h["identity"] < MIN_EXON_IDENTITY:
           continue
   
       if h["qcov"] < MIN_EXON_QCOV:
           continue
   
       h["exon"] = n
       exon_hits.append(h)
   
   
   for h in read_blast(short_out):
   
       n = exon_number(h["query"])
   
       if n is None:
           continue
   
       if h["identity"] < MIN_SHORT_IDENTITY:
           continue
   
       if h["qcov"] < MIN_SHORT_QCOV:
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
   # RECONSTRUCT AP1AR-LIKE LOCI
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
       # EVALUATE CLUSTERS
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
   
   
           exons_found = sorted(best_by_exon)
   
           if len(exons_found) < MIN_EXONS_PER_LOCUS:
               continue
   
           # Exon 9 is required because it defines the biological 3' end
           if 9 not in best_by_exon:
               continue
   
   
           selected = [
               best_by_exon[n]
               for n in exons_found
           ]
   
   
           # ----------------------------------------------------
           # Verify transcriptional order
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
   
   
           # ----------------------------------------------------
           # Exon 9 defines biological 3'
           # ----------------------------------------------------
   
           exon9 = best_by_exon[9]
   
           if strand == "+":
               three_prime = exon9["end"]
           else:
               three_prime = exon9["start"]
   
   
           # ----------------------------------------------------
           # Search TE toward biological 3'
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
   
           if nearby_te:
   
               distance, te = nearby_te[0]
   
               if exons_found == list(range(1, 10)):
                   status = "AP1AR_CHIMERIC_COMPLETE"
               else:
                   status = "AP1AR_CHIMERIC_PARTIAL"
   
               te_start = te["start"]
               te_end = te["end"]
               te_identity = te["identity"]
               te_qcov = te["qcov"]
   
           else:
   
               status = "AP1AR_LIKE_NO_3PRIME_TE"
   
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
       f"AP1AR_signature_{species}.tsv"
   )
   
   header = [
   
       "Species",
       "Genome",
       "Scaffold",
       "Strand",
   
       "AP1AR_start",
       "AP1AR_end",
   
       "Exons_found",
       "Number_of_exons",
   
       "Biological_3prime",
   
       "TE_start",
       "TE_end",
   
       "TE_identity",
       "TE_qcov",
   
       "TE_distance_3prime",
   
       "Status"
   ]
   
   
   with open(result_file, "w") as out:
   
       out.write("\t".join(header) + "\n")
   
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
   # SUMMARY
   # ============================================================
   
   print()
   print("=" * 72)
   print(f"AP1AR analysis: {species}")
   print("=" * 72)
   
   
   if not candidate_loci:
   
       print("No coherent AP1AR-like locus detected.")
   
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
   
           if x["status"].startswith("AP1AR_CHIMERIC"):
   
               print(
                   f"TE: {x['te_start']}-{x['te_end']}"
               )
   
               print(
                   f"Identity: {x['te_identity']:.3f}%"
               )
   
               print(
                   f"Coverage: {x['te_qcov']:.1f}%"
               )
   
               print(
                   f"3' distance: {x['te_distance']} bp"
               )
   
           print(f"-> {x['status']}")
   
   
   print()
   print(f"Results: {result_file}")

.. code-block:: bash

   chmod +x scripts/02_find_AP1AR_chimeric_signature.py

.. code-block:: bash

   python3 scripts/02_find_AP1AR_chimeric_signature.py \
   Plecotus_auritus


.. code-block:: bash

   mkdir -p results/batch_logs

   rm -f results/batch_logs/*.log 2>/dev/null
   rm -f results/batch_logs/*.failed 2>/dev/null

.. code-block:: bash

   find ../../SPECIES -mindepth 1 -maxdepth 1 -type d -print0 \
   | xargs -0 -n1 -P 3 sh -c '
       d="$1"
       species=$(basename "$d")
   
       echo "START: $species"
   
       if python3 scripts/02_find_AP1AR_chimeric_signature.py "$species" \
           > "results/batch_logs/${species}.log" 2>&1
       then
           echo "OK:    $species"
       else
           echo "FAILED: $species"
           touch "results/batch_logs/${species}.failed"
       fi
   ' sh

.. code-block:: bash

   grep -h -v '^Species' results/AP1AR_signature_*.tsv \
   | awk -F'\t' 'NF>1 {count[$NF]++} END {for (x in count) print x, count[x]}' \
   | sort

.. code-block:: bash

   grep -h 'AP1AR_CHIMERIC_' results/AP1AR_signature_*.tsv \
   | column -t
