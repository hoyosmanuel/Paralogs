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
