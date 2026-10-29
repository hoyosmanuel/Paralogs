I. IDENTIFICAR EL PARÁLOGO DE NCL
==================================

The NCL chimeric architecture occurs in a genomic copy located close
to the canonical NCL locus. In *Myotis myotis*, TOGA identifies the
candidate paralog as::

    ENST00000322723.9#NCL#165659

Because the canonical NCL locus occurs nearby and shares substantial
sequence similarity with the candidate paralog, TOGA annotation alone
is not used to identify the duplicated locus.

Instead, genomic exon sequences are used as probes to identify
NCL-like loci independently of annotation.

The TE-derived terminal exon is treated separately from the
NCL-derived exons. Exons 1–9 are therefore used to identify the host
gene architecture, whereas exon 10 is used as an independent probe
for the Charlie8-derived sequence

1) Create diagnostic NCL exon probes
------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/NCL
  
  nano sequences/NCL_host_exon_probes.fa

.. code-block:: bash

  >NCL_exon1_match
  GCCTCAGGCTCGAACTCTAGCTGGTAGCGAAACATCCGCGCCGGGGTCTTCCGTACCAGCGGCCCCCTGGGCATTGTTCGGTTCCTCTGGGCTCGTCTGCCCCCTCACCCGCCACACACCCACCGTGTAACCCACAAAGATAAGGGGCCTCCCTCCGTAAGCACGAACGTCCAGGCCCGTTTTCTGAGGAACGCGCGCAGGGCGCCGCGAGGAGGAGCGTGTGCGCGCCGTCCCGGGCGCATTCGAGGGCGTCCTGCTGCGGGAAGTCTCGTGTGATTAGCGCGGAGGTCTTGCGGTTTGGGGTTCTTGGTGGCGAG
  
  >NCL_exon2_match
  TCTGGTAAAAATCAAGGTGACCCCAAGAAAATGGCTTCTCCTCCAAAGGAGGTAGAAGAAGATAGTGAAGATGAGGAAGTGCCAGATGAAGAAGATGAAAGCAGTGGAGGAGAG
  
  >NCL_exon3_match
  GTTATCATTCCTCAGAAAAAAGGCAAAAAGGCTACCACAACTCCAGGAAAGAAGGCGATGGTTTCCCCAACAGAAAAGATTGCAGTTGCCACAGCAGGCAAGAAAACAGTTGTCACCCCTGGCAAAGAGGCAACAGCTATGCTGGCTCCTGGCAAAAAGGCAGCAGTTACACCAGCCAAAGCCATAGCAACCCCTGGTAAGAAGGAAGCAGCCACCCCAGACAAGGGAGCAGGGAATGCCAAGAATGCCAAGAAGGAAGATAGTGATGAGGAAGATGATGGCAACAGTGAAGAGGAGGATGAGGTTGAATTTGAGCCAGCAGTGATGAAAGCTACTGCTGCTGCCCCTTTCTTGGATGATGATGATGAAGAAGATG
  
  >NCL_exon4_match
  ACTCTGAAGAAGAACCTATGGAGATTACACCAGCCAAAGGCAAGAAAGCTCCTGTAAAAGCTGTTCTTGTGAAGGCCAAGAGCACAGAAAATGAAGACGAAGAGGAGGATGATGATGAAG
  
  >NCL_exon5_match
  ATCCTATTAAAAAAACACCTGATAAACAGAAGAAGAAAATGGCCAAACAGAAGGCTTCCCCTGAAACCAAAAAACAGAAATTAGAAG
  
  >NCL_exon6_match
  CTGAACCACCTGTAACTTTCAAACTGTTTGTTGGAAACCTGAATTTCAACAAAACTGCTGCTGAACTAAAAACTGGTCTCAGGGAATTTTTTGCTAAAAATGATCTCACAGTTGTTGGTGTCAGAGTCAGCTCGTCCAG
  
  >NCL_exon7_match
  GAGATTTGGTTATGTGAATTTTAAATCTGTTGACGATCAGGAAAAAGCCCTGGAACTCAGTGAGACAAAGGTCCTTGGGTATGAAATTAAACTGAAGAAACCAAAGGGAAAGGAAACAAAAAAAG
  
  >NCL_exon8_match
  GTCGAAATGCCAAAACACTTTTGATCAAGAACTTGCCTAGCAAGGTCACTGAGCATGAACTAAAGGAAGTGTTTGAGGATGCTTTTCAAATCAGATTTCTGAGCAATGATGGGATAAGTAAAAG
  
  >NCL_exon9_match
  AATTGCCTATATTGACTTCAAGTCACAAGCTGATGCAGAGAGAACCTTGGAAGAGAAGCAAGGGACAGAAATTGGAGGGCTTGCCATAGTTCTGGAGCATGTTGGAGAGAAAAATCAAGGCCAAGAAGAGAGGGCCAGGAAGAGCAGAAGTTGGAGAG

comprobamos que tenemos 9 sondas:

.. code-block:: bash

  grep '^>' sequences/NCL_host_exon_probes.fa


2) Create the Charlie8-derived exon probe
-----------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/NCL

  nano sequences/NCL_Charlie8_probe.fa


.. code-block:: bash

  >Charlie8_match
  gaaaacaaaagatagaacaagaaagccacacatttcaagaaaaatgggagagagcttatttctttgtggaggtaaagaatgttcctatgtgtttaatatgcaaagagagtctgtctgtgtcgaaagaatataacctaagatgccattatgagacaaaccacagtaggaaccttgttggatatacagaaaagatgcgtgataaaaaactcaatgaactgaaaaaaaggttgaaacttgaacatgatttgcttttgaatgtgaataaaataactgatgctgcaatgaaatgtagttacgtgttacgtgagaaaattgcccgggcatcaaaaccccttacagatggtgcgtttataaaggaatgtctactgagtgctgcagaaattttgtgccctgaacagagacaagtgtttgccaacataaatctcactggtaatattgctgagcagcatgttgccaatgtgactgacaacttacagaacaccttgcaagaaaaagttaaatcatttgtggctttctctattgcagctcatgaaagcacttatataaataacgccccccagttagctgtatttattcgcggtgttgatgaaacttttgatgtgactgaggaacttttggacatggtacccctgacgggcacaacatcgggaaatgacttatttttatgtgttgagaaaagtcttagaaaatttgatgtagactggtcaaaattagtaagtgttagcacagatggtaattctgcaacggttggtgtgaagcaacttgttacaaaacttaagtctaaagtgtcagggctttgcaaagacacagaacttaaatctgtgcatggcctctttctccaagaatccctttgtgctaaaaagttacaaatggatcatgtcatggacatagtaatttacaccataacgtggctgcattcccatggctcgacccacaggaagtttagtgctttgttcagtgagttagacacacagtatggacacctgacctacatggaactgaagtggctgagccgtggaatggtcctaaggcagttttttgaactgttggaagaaattgacattttcatgtcttctaaaggaaaatctgtgcctcagcttcgtaacaaagattgggtcaaagacctagcctttttggttgacattatgacttatgtaaacatgttggatgttgctctgcaagggcattcgcagatggttacacaaatgtaccatgcagttcgctctttcctagcaaaattgtgtctttgggaaactcatttggcaagaaacaatctggcccacttcccgacgctgagactagtttctgagaatgaaagtgatggactgaactacattcccaaaattaaagagttgaagactgaatttcaaaaaaggttctctgatttcaaactttatgaaaacgaacttatactgttcagctcacctttctcaattaatattaataatgtgaatgaagagctacaaatggaggttattgaattgcagcacaacactatattgaaaactaaatatgatgatgttggaatcccagaattcttcaaatatcttgggaacagttaccctaaatataaaaatcattgtgcaaagattctgtccatgtttggaagtacatatatttgtgagcaactattttctgttatgacacaaaataatactattgctcccagttaaaggattcagagatgaatttgctcctgcatattgctacgagaggtccccaacctgatattgactgatgggcattagtgagtaacatcccattttagatttcaaatctgaggagtaacaaatttaaaaactctgaataattatttctatttccaaattacatttttctaatgtcaacttaaagcacaatctctggcatcatttgtttgtaccttataacataaaataaagtttgattataccctggc

eliminar saltos:

.. code-block:: bash

  sed -i '' '/^>/! s/[[:space:]]//g' \
  sequences/NCL_Charlie8_probe.fa


3) Test the NCL exon probes in Myotis myotis
--------------------------------------------

.. code-block:: bash

  cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/NCL
  mkdir -p intermediate/blast_db
  makeblastdb \
  -in ../../SPECIES/Myotis_myotis/mMyoMyo1.6.pri.fa \
  -dbtype nucl \
  -out intermediate/blast_db/Myotis_myotis


