I. IDENTIFY THE KATNAL1–hAT1_mMyo GENOMIC ARCHITECTURE
======================================================

The KATNAL1 chimeric architecture was initially identified in
*Myotis daubentonii*.

The candidate transcript contains nine host-derived KATNAL1 exons
followed by a terminal exon derived from the transposable element
hAT1_mMyo::

    KATNAL1 exon 1
        |
    KATNAL1 exon 2
        |
    KATNAL1 exon 3
        |
    KATNAL1 exon 4
        |
    KATNAL1 exon 5
        |
    KATNAL1 exon 6
        |
    KATNAL1 exon 7
        |
    KATNAL1 exon 8
        |
    KATNAL1 exon 9
        |
    hAT1_mMyo-derived exon 10

To identify homologous genomic architectures independently of gene
annotation, KATNAL1 exons 1–9 are used as host-gene probes.

The TE-derived exon 10 is analyzed independently as an hAT1_mMyo probe.


1. Create the KATNAL1 analysis directory
-----------------------------------------

Create an independent working directory

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

    mkdir -p ANALYSES/KATNAL1/{intermediate,results,scripts,sequences}

    cd ANALYSES/KATNAL1

Confirm the working directory

.. code-block:: bash

    pwd

The expected path is

.. code-block:: bash

    /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/KATNAL1


2. Create the KATNAL1 host exon probes
---------------------------------------

Create a FASTA file containing KATNAL1 exons 1–9

.. code-block:: bash

    nano sequences/KATNAL1_host_exon_probes.fa

Use the following FASTA headers and paste the corresponding sequences
supplied from the *Myotis daubentonii* transcript

.. code-block:: bash

  >KATNAL1_exon1_match
  gtaaagaatttttgtgagtttatttgagccaaactgtcgacaattgccgggaagaagtttctcaatggattgggataatgctccagagaatggcagttttttcatctatttttatacaatgtaatcaaaggaggggacttaggtgggttacttgaaatccacgggtgataggagggagaaagcaaaggggcggggaggggagaacgtctgggattggctaaaaagtaaaatgatagacacatacttaggtgggtgaggatagttaacagtcaacaattaacatgataacaataacaatgaaggaatttgtagtatctgtcctggcgcccaagcccaagcacattaggttgtgccccaaggggtctggaaaaatggaagttacaagttacccagacatttcaagaatatgttatcatagatgcaaaaagacagataggcacagttaaggtaaaggttgaccttgtcagtgaagatactggcctaggacataactacccaccataactgcttttagttaaagtttatttttttatttaaaaatatatttttgttgacttcagagaggaagggagaggagggagacatagaaacatcagtgatgagagtgaatcattgatttgctgcctcctgcacaccccccactggggatcaagcccacaacctgggcgtgtgcccttgacaggaatcgaccccccggggaaccttcagtccacaggctggtgctctatccactgagctaaaccagccagggcatgtttcagactatcctttgtggttactttaggtctttttttggcaagactgccttgcagcttgtcaggtaagcatgtggcccctttggtccacaCCTgtttgattccccctcagggcacacgcctgggttgatggctggatccccagttgggggcatgcaggaggcagccaattgaagtttctctcttatcgatgtttctatctctgtctcccttcttctctctctaaaatcaacaacatatatatttttttaaatgcaagggAAACCCATGCTACTTGCAAGTTCTCAAAACTCCTCAATCAGCCTAAAAGCCGTATGTATTACATCTTTCTAAAAGAGTCCCTTCTCTCCCACTGCTCAGTAAATGGCAGGAAGGAGGTGTGAAAGGTTGCAACTGTGATTTCATTCGTGAGAGCCCAAACAGAGGTGGCTCAGACTGCAAAGCGGCGGAAGGGATCAGCCCGCCGGGTGCCCTCCTCTGATGACGCAGGTGTCCTCTTAGGATCCAGCCGCATTCGCCAGCTCAGCCTCGGAACAGGGCACGCTGCCACCGGGATGGCTGCTGGCAGGGCTCATTAAGCAGCCCTAAAGCGCAGTCACCTGTGCTCCTCAATGTCACCGTCTCCCACTCTCAAACACGCACCTTTGCACCTCGTTTCCACCGCATGCGGCCGAATTCCCAGGCTCCTCCCTGGAGCTGACTCCGAGAACCCCCTCACTGACACCACCGCGCAGCCCGCGCTGTCATGGCAACCGTACGCCGCACCACCGCCTCTTCACGCACGTCGTGAACTAACGGACTCCACAGCGGCGGAGCGGTGGCCGGAGCTCCACCCCCACCCCACCAGACCTTGACCGCATAACCAACAGTCGGCCAGCCTGCGGCACCGGG
  
  >KATNAL1_exon2_match
  GTCTCTGAAAGAAGATGAACTTGGCTGAGATTTGTGACAATGCAAAGAAAGGAAGAGAATATGCACTTCTTGGAAATTATGACTCATCAATGGTATATTACCAAGGGGTGATCCAGCAGATTCAGAGACATTGCCAGTCAGTTAGAGACCCAGCTATCAGAGGCAAATGGCAACAG
  
  >KATNAL1_exon3_match
  GTTCGGCAGGAATTATTGGAAGAATATGAACAAGTTAAAAGTATTGTCAACACTTTGGAAAGTTTTAAAATTGACAAGCCCCCAGATTTCCCTGTGTCCTGTCAAGACGAACCATTTAGAGACCCTGCTGTTTGGCCGCCCCCTGTACCTGCAGAACACAG
  
  >KATNAL1_exon4_match
  AGCTCCACCTCAGATAAGGCGTCCCAATCGTGAAGTAAAACCTCTGAGGAAAGAAATGCCAGGAGCAGGAGCCCGGGGACCTGTAGGTCGAGCACATCCTATATCGAAGAGTGAGAAGCCCTCTACAAGTAGGGATAAGGATTGTAGAGCAAGAGGGAGAGATGACAAG
  
  >KATNAL1_exon5_match
  GGAAGAAAAAATATGCAAGATGGTGTGAGTGATGGTGAAATTCCAAAATTCGATGCTGCTGGATATGATAAGGATCTGGTGGAAGCCCTTGAGAGAGATATTGTGTCTAGGAATCTTAGCATTCATTG
  
  >KATNAL1_exon6_match
  GGATGACATAGCAGATCTGGAAGAAGCTAAGAAGTTGCTAAGGGAAGCTGTGGTTCTTCCCATGTGGATGCCTGACTTCTTCAAAGGGATTAGAAGGCCCTGGAAG
  
  >KATNAL1_exon7_match
  GGTGTGCTGATGGTTGGACCCCCAGGCACTGGTAAGACCATGCTAGCTAAAGCTGTTGCCACTGAATGTGGCACAACATTCTTCAATGTTTCCTCCTCTACACTGACATCCAAATATAGAGGTGAATCTGAGAAATTAGTCCGTCTGTTGTTTGAAATG
  
  >KATNAL1_exon8_match
  GCTAGATTTTACGCCCCCACTACGATCTTCATTGATGAGATAGATTCTATCTGCAGTCGAAGAGGGACCTCTGATGAGCACGAAGCAAGTCGCAGAGTCAAGTCTGAGCTACTCATTCAGATGGATG
  
  >KATNAL1_exon9_match
  GTGTTGGAGGAGCTTTAGAGAATGATGATCCTTCCAAAATGGTTATGGTATTAGCTGCTACTAATTTCCCATGGGATATTGATGAAGCTTTGCGAAGAAGATTAGAGAAGAGGATTTATATACCTCTCCCAACAG


Remove accidental spaces and other whitespace from sequence lines

.. code-block:: bash

    sed -i '' '/^>/! s/[[:space:]]//g' \
    sequences/KATNAL1_host_exon_probes.fa

Confirm the FASTA headers

.. code-block:: bash

    grep '^>' sequences/KATNAL1_host_exon_probes.fa

The expected output is

.. code-block:: bash

    >KATNAL1_exon1_match
    >KATNAL1_exon2_match
    >KATNAL1_exon3_match
    >KATNAL1_exon4_match
    >KATNAL1_exon5_match
    >KATNAL1_exon6_match
    >KATNAL1_exon7_match
    >KATNAL1_exon8_match
    >KATNAL1_exon9_match


3. Create the hAT1_mMyo TE-derived exon probe
----------------------------------------------

Exon 10 is derived from the transposable element hAT1_mMyo and is
therefore analyzed independently from the host-derived KATNAL1 exons.

Create the TE-derived probe FASTA

.. code-block:: bash

    nano sequences/KATNAL1_hAT1_mMyo_probe.fa

Use the following header and paste the complete exon 10 sequence

.. code-block:: bash

    >hAT1_mMyo_match
    gtgccctgtttgacatggctacaaagaagaagcaaagaaaaacggaagatgaaaaccgtgaatttaaagttgaatggaccaagattttcgcgtttattcaaaacttaaatggccttccgacctgtcttatttgtcaggaaaaattggcgcataacaagaaatcaaatttggagagacactttacaactaaacacgtgtcatttagtactaactatcctgtcggtgatgcaaggaagaaagcagttgaggaacttcaaaagagtcaagaaaaatcaagttctgtatttaattactggatgcaatcttccaacaatgtcaatattgcaagttttgtggttagtcaagagattggtaagagaggaaaaccatatacaaacggcgaatacataaaaagttgttttataaatgcatctgaagagctatttcgggattttaagaacaaagcagatattctaaaaaaaattaaagagttaccattgtctgctaaaacagtgaaagatagaacagtcaaaatgtcttcaaatataaccaaccagcaagtcgaagatcttaaattggtttcagctttatcaatagcagttgacgagtcttgtgacataaatgatacagcgcaagtttcactttttgtgcaatttatttcatctacaggtcctaaagagaacttttaggattattgccactcaaaggtcaaacacgtggagaggatatagcaaatgctgtaattgagtgcattggaaaacatcatattccactcgataaaattgtctcaatttcaacagacgggtcgaaaagtatgaccggcgtaagaaacgggtttgttgctattttgaaagaaaaaattaatcatgagatacttacttaccattgcatcattcatcaagaagcgctttgtgggcagacatttccagaagaaatttgcaaagttatggaattggtgattaagattatcaactccattatAgctaaagctcttaatcatcgccagtttaaagaatttttagttgaaatggagagtgagtacgcggatcttctgctgcataacaaggtatgttggttgtcaagaggaaatgttttaaaaagtttcgcttccttattaacagaaattaaagcatttctccttgagaagggtgttcactatccagaactaacaaacgatcagtggatccaaaaattttatttcatggtggaggttacatctcatctaaatcagcttaatcgtaaattacaggagaaaggaaacacaattttttcgatgttagaagaagttattttattcaaaaacaaatgatctctttttgctcgagattttgaaagagaaactttgattcactttccaagcctgttgaaacatcgccaagaaaataattctgatattgacatttgctatttcaaaacaacagttttaaatatgagggaagcttttctcaacaggtttcaggaattcagaaacagcaaagtgacgttagcctttgttaaaaaacctctcaatgctactgtaacagaattaaatttttctccctttaatattaacattggcgactttgaaatgcaattgctggattttaaaaaacaaagaactgtggagctcgaaatttgaatgtctttgtgctgatatagaaatattggagaaaaacaaatgtgaacttagttcacagcacaagtggtctgctttgaaagacctggagaaggaagatattttgatttttaacacctggaatagtattcctgactcacataatcaactgaaaaaactagcgtttgctgttctttccttattcgggtctacatatttatgcgaacaatcattttcaagcatgaaccttatcaagagtaaattgagaagtcgtcttattgatgagaacctggaatcatgcctaaaattaaaaacaacaacatacaaacctgacttacccaaactttccaaggagatgcaaggccattgttcacattagtgtttgtcagtcattgtcataatttaattttatggatattgtaaagtaaattgtaaggaatttaattttatcaataaataatatcgttattaagtatgatatcaagttttattcaatgtccctatagtttaactaagacttaaaactttaagtagagtttattaagttaactttagggaacgttgtggagtgaaagaagatgtagtggagaggattgagaaaggtatgttgagatggtttggtcatacagagaggatgaaCgaaaggagacagacgaaacaagtatacaagacgagtatggatgggagagttggaaggggtcgacctcagcgaacatacctcaatcagactgaggatgtttttagcaaaggccagcttaggagtaccctaattaagttaataacaatgtacctacctatatagtttaagtttaaaaaatttggctctcaaaaaaaatttcaatcattgtactgttgatatttgactctgttgactaataagtttgctgaccactagtctagAGCTTTGGACAAAGCTCTCTCCCAGTCTTCATTGCTTTGTAGTTGTTGATTTTTTAAAACCTTTTAAAATGGATTCTTATATATTTTTTGACTATTCACTTTCACTGCCTATGTCTTTTATTTTTATAAACAGGTGGTATATCATATTGACAGTTCTCTAGAAATAATTTCACAAAGTGGCTTCTATTATAAACTGCCTTAGCTAAAGGTATCCTTTCCTTTTTTACATTTTTATTTTTATTTATttattttttaaatgagagcagttttaggtttacaataaaattgagaaggagatatagagatttcctatatactccctgccccccacacatgcacagcctcccccattatcaacatcactcaacagaggagtcaacttattattattattttaccaaaaatgaaccaacattggcacattgttatcaccccaaatccatagtttaccttaggttacattctGTGGGTTTAagccgtgggcaaaccacagcccgcgggccagatccggcccatttgaaatg

Remove accidental whitespace

.. code-block:: bash

    sed -i '' '/^>/! s/[[:space:]]//g' \
    sequences/KATNAL1_hAT1_mMyo_probe.fa

Confirm the header

.. code-block:: bash

    head -1 sequences/KATNAL1_hAT1_mMyo_probe.fa

The expected output is

.. code-block:: bash

    >hAT1_mMyo_match


4. Verify probe lengths
------------------------

Verify the lengths of all host-derived probes

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
    ' sequences/KATNAL1_host_exon_probes.fa

Verify the length of the hAT1_mMyo-derived exon

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
    ' sequences/KATNAL1_hAT1_mMyo_probe.fa
