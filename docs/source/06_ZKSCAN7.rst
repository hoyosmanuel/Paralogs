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
