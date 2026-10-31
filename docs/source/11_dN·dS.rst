Preparación de secuencias codificantes para análisis de dN/dS
==============================================================

Objetivo
--------

El objetivo de este análisis es reconstruir las regiones codificantes homólogas de
las arquitecturas quiméricas identificadas durante el análisis genómico final.

En esta etapa **no se vuelve a determinar qué especies presentan cada arquitectura**.
Ese trabajo ya fue realizado en los análisis individuales de ``FBXO3``, ``TMPO``,
``NCL``, ``ZKSCAN7``, ``KATNAL1`` y ``ZFP28``.

Aquí se utilizan esos loci finales para reconstruir, cuando sea posible:

* la CDS derivada del gen hospedero (``HOST CDS``);
* la CDS derivada del elemento transponible (``TE CDS``);
* la CDS quimérica completa (``CHIMERA CDS``);
* sus respectivas traducciones.

Las secuencias se reconstruyen literalmente a partir del genoma.

No se deben corregir artificialmente:

* frameshifts;
* codones de parada internos;
* deleciones;
* bases ambiguas;
* regiones ausentes.

El análisis de ``dN/dS`` solamente se realizará después de completar la
reconstrucción y el control de calidad de las secuencias.


I. Crear el espacio de trabajo
==============================

Todos los archivos de esta etapa se mantienen separados de los análisis utilizados
originalmente para identificar las arquitecturas.

.. code-block:: bash

   cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis

   mkdir -p ANALYSES/dNdS_FINAL/{config,scripts,sequences,results,logs}

   mkdir -p ANALYSES/dNdS_FINAL/reference_transcripts
   mkdir -p ANALYSES/dNdS_FINAL/reference_cds
   mkdir -p ANALYSES/dNdS_FINAL/reference_probes


La estructura principal queda:

::

   ANALYSES/dNdS_FINAL/
   ├── config/
   ├── scripts/
   ├── reference_transcripts/
   ├── reference_cds/
   ├── reference_probes/
   ├── sequences/
   ├── results/
   └── logs/


II. Definir el conjunto final de especies
=========================================

Para evitar que especies utilizadas durante exploraciones preliminares entren
accidentalmente en el análisis final, se creó un manifiesto explícito.

Este archivo contiene exclusivamente las especies incluidas en Supplementary File 1,
Table 3, para las seis arquitecturas que actualmente tienen al menos tres especies
potencialmente comparables.

.. code-block:: bash

   cat > ANALYSES/dNdS_FINAL/config/final_architecture_species.tsv <<'EOF'
   Architecture	Species	Table_status	RNA_support	PreQC_note
   FBXO3	Rhinolophus_affinis	Near-complete	NO	.
   FBXO3	Rhinolophus_ferrumequinum	Near-complete	YES	.
   FBXO3	Rhinolophus_hipposideros	Near-complete	YES	.
   FBXO3	Rhinolophus_pearsonii	Near-complete	NO	.
   FBXO3	Rhinolophus_sinicus	Near-complete	YES	.
   FBXO3	Rhinolophus_trifoliatus	Near-complete	NO	.
   TMPO	Myotis_auriculus	Complete	NO	.
   TMPO	Myotis_californicus	Complete	NO	.
   TMPO	Myotis_daubentonii	Complete	YES	.
   TMPO	Myotis_evotis	Complete	NO	.
   TMPO	Myotis_lucifugus	Complete	NO	.
   TMPO	Myotis_myotis	Complete	YES	.
   TMPO	Myotis_mystacinus	Complete	YES	.
   TMPO	Myotis_nigricans	Near-complete	NO	.
   TMPO	Myotis_occultus	Complete	NO	.
   TMPO	Myotis_pilosus	Partial_TE_recovery	NO	TE_RECOVERY_85_PERCENT
   TMPO	Myotis_thysanodes	Complete	NO	.
   TMPO	Myotis_velifer	Near-complete	NO	.
   TMPO	Myotis_vivesi	Near-complete	NO	.
   TMPO	Myotis_volans	Near-complete	NO	.
   TMPO	Myotis_yumanensis	Near-complete	NO	.
   NCL	Antrozous_pallidus	Complete	NO	.
   NCL	Corynorhinus_mexicanus	Complete	NO	.
   NCL	Corynorhinus_townsendii	Complete	NO	.
   NCL	Eptesicus_fuscus	Complete	NO	.
   NCL	Eptesicus_nilssonii	Complete	NO	.
   NCL	Lasiurus_ega	Complete	NO	.
   NCL	Myotis_auriculus	Complete	NO	.
   NCL	Myotis_californicus	Complete	NO	.
   NCL	Myotis_daubentonii	Complete	YES	.
   NCL	Myotis_evotis	Complete	NO	.
   NCL	Myotis_lucifugus	Complete	NO	.
   NCL	Myotis_myotis	Complete	YES	.
   NCL	Myotis_mystacinus	Complete	YES	.
   NCL	Myotis_nigricans	Complete	NO	.
   NCL	Myotis_occultus	Complete	NO	.
   NCL	Myotis_pilosus	Complete	NO	.
   NCL	Myotis_thysanodes	Complete	NO	.
   NCL	Myotis_velifer	Complete	NO	.
   NCL	Myotis_vivesi	Complete	NO	.
   NCL	Myotis_volans	Complete	NO	.
   NCL	Myotis_yumanensis	Complete	NO	.
   NCL	Pipistrellus_kuhlii	Complete	NO	AMBIGUOUS_LOCUS_ASSIGNMENT
   NCL	Plecotus_auritus	Complete	YES	.
   ZKSCAN7	Rhinolophus_affinis	Complete	NO	.
   ZKSCAN7	Rhinolophus_ferrumequinum	Complete	YES	.
   ZKSCAN7	Rhinolophus_foetidus	Complete	NO	.
   ZKSCAN7	Rhinolophus_hipposideros	Complete	YES	.
   ZKSCAN7	Rhinolophus_pearsonii	Complete	NO	.
   ZKSCAN7	Rhinolophus_perniger_lanosus	Complete	NO	.
   ZKSCAN7	Rhinolophus_trifoliatus	Complete	NO	.
   KATNAL1	Myotis_auriculus	Complete	NO	.
   KATNAL1	Myotis_californicus	Complete	NO	.
   KATNAL1	Myotis_daubentonii	Complete	YES	.
   KATNAL1	Myotis_evotis	Complete	NO	.
   KATNAL1	Myotis_myotis	Complete	NO	.
   KATNAL1	Myotis_mystacinus	Complete	NO	.
   KATNAL1	Myotis_volans	Complete	NO	.
   ZFP28	Cynopterus_sphinx	High-confidence_complete	NO	.
   ZFP28	Desmodus_rotundus	Exon7_uncertain	NO	UNCERTAIN_EXON7_BOUNDARY
   ZFP28	Eonycteris_spelaea	High-confidence_complete	YES	.
   ZFP28	Eumops_nanus	High-confidence_complete	NO	.
   ZFP28	Hipposideros_abae	High-confidence_complete	NO	.
   ZFP28	Hipposideros_armiger	High-confidence_complete	NO	.
   ZFP28	Hipposideros_caffer	High-confidence_complete	NO	.
   ZFP28	Hipposideros_jonesi	High-confidence_complete	NO	.
   ZFP28	Hipposideros_larvatus	High-confidence_complete	NO	.
   ZFP28	Hipposideros_swinhoei	High-confidence_complete	NO	.
   ZFP28	Hypsignathus_monstrosus	High-confidence_complete	NO	.
   ZFP28	Molossus_alvarezi	High-confidence_complete	NO	.
   ZFP28	Molossus_molossus	High-confidence_complete	NO	.
   ZFP28	Rhinolophus_affinis	High-confidence_complete	NO	.
   ZFP28	Rhinolophus_ferrumequinum	High-confidence_complete	NO	.
   ZFP28	Rhinolophus_foetidus	Exon7_degraded	NO	DEGRADED_EXON7_BOUNDARY
   ZFP28	Rhinolophus_hipposideros	High-confidence_complete	NO	.
   ZFP28	Rhinolophus_pearsonii	High-confidence_complete	NO	.
   ZFP28	Rhinolophus_perniger_lanosus	High-confidence_complete	NO	.
   ZFP28	Rhinolophus_trifoliatus	High-confidence_complete	NO	.
   ZFP28	Rhinopoma_microphyllum	High-confidence_complete	NO	.
   ZFP28	Rhinopoma_muscatellum	High-confidence_complete	NO	.
   ZFP28	Rousettus_aegyptiacus	High-confidence_complete	YES	.
   ZFP28	Tadarida_brasiliensis	High-confidence_complete	NO	.
   EOF


Las arquitecturas ``AP1AR``, ``ZNF24`` y ``FAU`` se dejan por fuera de esta primera
ronda porque actualmente no tienen al menos tres recuperaciones completas y
comparables.


III. Identificar los transcritos quiméricos de referencia
=========================================================

Para reconstruir correctamente una CDS es necesario conocer primero la estructura
del transcrito quimérico real.

Los transcritos con evidencia de RNA utilizados como referencia son:

.. list-table::
   :header-rows: 1

   * - Arquitectura
     - Especie
     - RNA-seq
     - Transcrito
   * - FBXO3
     - Rhinolophus ferrumequinum
     - SRR18652180
     - STRG.3111.1
   * - ZKSCAN7
     - Rhinolophus ferrumequinum
     - SRR18652180
     - STRG.4351.1
   * - TMPO
     - Myotis myotis
     - SRR18652219
     - STRG.7092.1
   * - NCL
     - Myotis myotis
     - SRR18652219
     - STRG.3566.1
   * - ZFP28
     - Rousettus aegyptiacus
     - SRR18652228
     - STRG.7996.1
   * - KATNAL1
     - Myotis daubentonii
     - ERR13757851
     - STRG.3154.1


Inicialmente se buscaron candidatos en los resultados transcriptómicos:

.. code-block:: bash

   grep -RniE 
   'FBXO3|ZKSCAN7|ZFP28|KATNAL1|NCL|TMPO' 
   /Volumes/Expansion/project3/RESULTS 
   /Volumes/Expansion/project3/TRANSCRIPTOME+ANNOTATION 
   /Volumes/Expansion/project3/TRANSCRIPTOME_SOLO 
   2>/dev/null 
   | grep -E 'STRG\.|transcript_id' \
   > ANALYSES/dNdS_FINAL/reference_transcripts/candidate_transcript_hits.txt


Advertencia importante
----------------------

Los identificadores ``STRG`` producidos por StringTie **no son globalmente únicos**.

Por ejemplo, un identificador como:

::

   STRG.3111.1

puede aparecer en distintas especies o corridas.

Por tanto, un transcrito solamente se considera identificado correctamente cuando
se conocen simultáneamente:

::

   especie
   corrida RNA-seq
   STRG ID


IV. Transcritos de referencia recuperados
=========================================

Se recuperaron localmente cinco transcritos completos:

::

   FBXO3/
       FBXO3_rFer_SRR18652180_STRG.3111.1.fa
       FBXO3_rFer_SRR18652180_STRG.3111.1.gtf

   ZKSCAN7/
       ZKSCAN7_rFer_SRR18652180_STRG.4351.1.fa
       ZKSCAN7_rFer_SRR18652180_STRG.4351.1.gtf

   TMPO/
       TMPO_mMyo_SRR18652219_STRG.7092.1.fa
       TMPO_mMyo_SRR18652219_STRG.7092.1.gtf

   NCL/
       NCL_mMyo_SRR18652219_STRG.3566.1.fa
       NCL_mMyo_SRR18652219_STRG.3566.1.gtf

   ZFP28/
       ZFP28_rAeg_SRR18652228_STRG.7996.1.fa
       ZFP28_rAeg_SRR18652228_STRG.7996.1.gtf


Las longitudes observadas fueron:

::

   FBXO3     1725 nt
   ZKSCAN7   2602 nt
   TMPO      2610 nt
   NCL       3478 nt
   ZFP28     2746 nt


V. Inventariar ORFs en los transcritos de referencia
=====================================================

Se creó el script:

::

   ANALYSES/dNdS_FINAL/scripts/01_inventory_reference_ORFs.py


Su función es examinar los transcritos recuperados y registrar todos los ORFs
suficientemente largos para inspección posterior.

El resultado se almacena en:

::

   ANALYSES/dNdS_FINAL/results/reference_ORF_inventory.tsv


La lógica general utilizada por el script es:

.. code-block:: python

   from pathlib import Path
   from Bio.Seq import Seq
   from Bio import SeqIO
   import csv

   ROOT = Path(
       "/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
   )

   INPUT_ROOT = (
       ROOT
       / "ANALYSES"
       / "dNdS_FINAL"
       / "reference_transcripts"
   )

   OUTPUT = (
       ROOT
       / "ANALYSES"
       / "dNdS_FINAL"
       / "results"
       / "reference_ORF_inventory.tsv"
   )

   MIN_AA = 50

   rows = []

   for fasta in sorted(INPUT_ROOT.rglob("*.fa")):

       for record in SeqIO.parse(fasta, "fasta"):

           sequence = str(record.seq).upper()

           for frame in range(3):

               for start in range(frame, len(sequence) - 2, 3):

                   if sequence[start:start + 3] != "ATG":
                       continue

                   for end in range(start + 3, len(sequence) - 2, 3):

                       codon = sequence[end:end + 3]

                       if codon in {"TAA", "TAG", "TGA"}:

                           cds = sequence[start:end + 3]

                           protein = str(
                               Seq(cds).translate()
                           )

                           aa_length = len(protein.rstrip("*"))

                           if aa_length >= MIN_AA:

                               rows.append(
                                   {
                                       "File": fasta.name,
                                       "Start": start + 1,
                                       "End": end + 3,
                                       "CDS_nt": len(cds),
                                       "Protein_aa": aa_length,
                                       "Protein": protein,
                                   }
                               )

                           break

   OUTPUT.parent.mkdir(
       parents=True,
       exist_ok=True
   )

   with OUTPUT.open(
       "w",
       newline=""
   ) as handle:

       fields = [
           "File",
           "Start",
           "End",
           "CDS_nt",
           "Protein_aa",
           "Protein",
       ]

       writer = csv.DictWriter(
           handle,
           fieldnames=fields,
           delimiter="\t"
       )

       writer.writeheader()
       writer.writerows(rows)


Ejecutar:

.. code-block:: bash

   python3 ANALYSES/dNdS_FINAL/scripts/01_inventory_reference_ORFs.py


El inventario produjo 96 ORFs de al menos 50 aminoácidos.

La ORF correcta para cada arquitectura se determinó comparando las proteínas
resultantes con las proteínas quiméricas identificadas originalmente.

Las ORFs seleccionadas fueron:

.. list-table::
   :header-rows: 1

   * - Arquitectura
     - Inicio
     - Final
     - CDS
     - Proteína
   * - FBXO3
     - 16
     - 1674
     - 1659 nt
     - 552 aa
   * - ZKSCAN7
     - 443
     - 1501
     - 1059 nt
     - 352 aa
   * - TMPO
     - 237
     - 2468
     - 2232 nt
     - 743 aa
   * - NCL
     - 348
     - 3233
     - 2886 nt
     - 961 aa
   * - ZFP28
     - 73
     - 1245
     - 1173 nt
     - 390 aa


VI. Definir dónde termina HOST y comienza TE
============================================

Para cada arquitectura se identificó el exon terminal derivado del TE.

El archivo utilizado para guardar esas posiciones es:

::

   ANALYSES/dNdS_FINAL/config/reference_coding_segments.tsv


Se creó con:

.. code-block:: bash

   cat > ANALYSES/dNdS_FINAL/config/reference_coding_segments.tsv <<'EOF'
   Architecture	Reference_species	Run	Transcript	ORF_start	ORF_end	Terminal_TE_exon_start	Host_CDS_start	Host_CDS_end	TE_CDS_start	TE_CDS_end	Protein_aa
   FBXO3	Rhinolophus_ferrumequinum	SRR18652180	STRG.3111.1	16	1674	374	16	373	374	1674	552
   ZKSCAN7	Rhinolophus_ferrumequinum	SRR18652180	STRG.4351.1	443	1501	1254	443	1253	1254	1501	352
   TMPO	Myotis_myotis	SRR18652219	STRG.7092.1	237	2468	649	237	648	649	2468	743
   NCL	Myotis_myotis	SRR18652219	STRG.3566.1	348	3233	1561	348	1560	1561	3233	961
   ZFP28	Rousettus_aegyptiacus	SRR18652228	STRG.7996.1	73	1245	884	73	883	884	1245	390
   EOF


VII. Extraer HOST CDS, TE CDS y CHIMERA CDS de las referencias
==============================================================

Se creó:

::

   ANALYSES/dNdS_FINAL/scripts/02_extract_reference_components.py


La lógica del script consiste en:

1. leer ``reference_coding_segments.tsv``;
2. localizar el FASTA correspondiente;
3. cortar la región HOST;
4. cortar la región TE;
5. concatenarlas para reconstruir la CDS;
6. traducir las tres secuencias;
7. escribir FASTA independientes.


Una implementación equivalente es:

.. code-block:: python

   from pathlib import Path
   from Bio import SeqIO
   from Bio.Seq import Seq
   import csv

   ROOT = Path(
       "/Volumes/Expansion/project3/bat_HTF_genomic_analysis"
   )

   CONFIG = (
       ROOT
       / "ANALYSES"
       / "dNdS_FINAL"
       / "config"
       / "reference_coding_segments.tsv"
   )

   TRANSCRIPTS = (
       ROOT
       / "ANALYSES"
       / "dNdS_FINAL"
       / "reference_transcripts"
   )

   OUTPUT = (
       ROOT
       / "ANALYSES"
       / "dNdS_FINAL"
       / "reference_cds"
   )

   with CONFIG.open() as handle:

       reader = csv.DictReader(
           handle,
           delimiter="\t"
       )

       for row in reader:

           architecture = row["Architecture"]

           candidates = list(
               (
                   TRANSCRIPTS
                   / architecture
               ).glob("*.fa")
           )

           if len(candidates) != 1:

               print(
                   f"{architecture}: "
                   f"expected one FASTA, found "
                   f"{len(candidates)}"
               )

               continue

           fasta = candidates[0]

           record = next(
               SeqIO.parse(
                   fasta,
                   "fasta"
               )
           )

           transcript = str(
               record.seq
           ).upper()

           host_start = int(
               row["Host_CDS_start"]
           )

           host_end = int(
               row["Host_CDS_end"]
           )

           te_start = int(
               row["TE_CDS_start"]
           )

           te_end = int(
               row["TE_CDS_end"]
           )

           host = transcript[
               host_start - 1:host_end
           ]

           te = transcript[
               te_start - 1:te_end
           ]

           chimera = host + te

           outdir = (
               OUTPUT
               / architecture
           )

           outdir.mkdir(
               parents=True,
               exist_ok=True
           )

           def write_fasta(
               filename,
               header,
               sequence
           ):

               with (
                   outdir
                   / filename
               ).open("w") as out:

                   out.write(
                       f">{header}\n"
                   )

                   out.write(
                       f"{sequence}\n"
                   )

           write_fasta(
               "host_cds.fa",
               f"{architecture}|HOST",
               host
           )

           write_fasta(
               "te_cds.fa",
               f"{architecture}|TE",
               te
           )

           write_fasta(
               "chimera_cds.fa",
               f"{architecture}|CHIMERA",
               chimera
           )

           write_fasta(
               "host_protein.fa",
               f"{architecture}|HOST",
               str(
                   Seq(host).translate()
               )
           )

           write_fasta(
               "te_protein.fa",
               f"{architecture}|TE",
               str(
                   Seq(te).translate()
               )
           )

           write_fasta(
               "chimera_protein.fa",
               f"{architecture}|CHIMERA",
               str(
                   Seq(chimera).translate()
               )
           )

           print(
               architecture,
               f"HOST={len(host)}",
               f"TE={len(te)}",
               f"CHIMERA={len(chimera)}",
           )


Ejecutar:

.. code-block:: bash

   python3 
   ANALYSES/dNdS_FINAL/scripts/02_extract_reference_components.py


Los resultados obtenidos fueron:

::

   FBXO3
       HOST     358 nt
       TE       1301 nt
       CHIMERA  1659 nt
       Protein  552 aa

   ZKSCAN7
       HOST     811 nt
       TE       248 nt
       CHIMERA  1059 nt
       Protein  352 aa

   TMPO
       HOST     412 nt
       TE       1820 nt
       CHIMERA  2232 nt
       Protein  743 aa

   NCL
       HOST     1213 nt
       TE       1673 nt
       CHIMERA  2886 nt
       Protein  961 aa

   ZFP28
       HOST     811 nt
       TE       362 nt
       CHIMERA  1173 nt
       Protein  390 aa


VIII. Caso especial de KATNAL1
==============================

El transcrito original de referencia es:

::

   Myotis daubentonii
   ERR13757851
   STRG.3154.1


El archivo transcriptómico original no se encontró actualmente en el árbol local.

Sin embargo, el GTF original conserva las coordenadas completas del transcrito.

El transcrito está localizado en:

::

   NC_081841
   strand = -


y contiene diez exones.

Como el transcrito está en strand negativo, el orden biológico 5' → 3' es el
inverso del orden de las coordenadas genómicas.

Primero localizar el genoma:

.. code-block:: bash

   ROOT="/Volumes/Expansion/project3"

   find "$ROOT/SPECIES" 
   -type f 
   -name 'mMyoDau2.1.pri.fa' 
   2>/dev/null


Definir la ruta:

.. code-block:: bash

   ROOT="/Volumes/Expansion/project3"

   PROJ="$ROOT/bat_HTF_genomic_analysis"

   OUT="$PROJ/ANALYSES/dNdS_FINAL/reference_transcripts/KATNAL1"

   mkdir -p "$OUT"

   GENOME=$(
       find "$ROOT/SPECIES" 
       -type f 
       -name 'mMyoDau2.1.pri.fa' 
       2>/dev/null 
       | head -1
   )

   echo "Genome: $GENOME"

   samtools faidx "$GENOME"


Crear el archivo con los exones en orden transcripcional:

.. code-block:: bash

   cat > "$OUT/KATNAL1_STRG.3154.1_exons.tsv" <<'EOF'
   1	NC_081841	104836882	104838527
   2	NC_081841	104825954	104826129
   3	NC_081841	104823090	104823250
   4	NC_081841	104785802	104785970
   5	NC_081841	104766400	104766527
   6	NC_081841	104764568	104764673
   7	NC_081841	104759185	104759343
   8	NC_081841	104758356	104758482
   9	NC_081841	104742463	104742597
   10	NC_081841	104694088	104697084
   EOF


Reconstruir el transcrito.

Como está en strand negativo, cada exon debe ser reverse-complementado.

.. code-block:: bash

   : > "$OUT/KATNAL1_mDau_ERR13757851_STRG.3154.1.fa"

   echo 
   ">STRG.3154.1 reconstructed_from_genome species=Myotis_daubentonii strand=-" \
   > "$OUT/KATNAL1_mDau_ERR13757851_STRG.3154.1.fa"

   while IFS=$'\t' read -r exon chr start end
   do

       samtools faidx \
       "$GENOME" \
       "${chr}:${start}-${end}" \
       | sed '1d' \
       | tr -d '\n' \
       | tr 'ACGTacgt' 'TGCAtgca' \
       | rev

   done 
   < "$OUT/KATNAL1_STRG.3154.1_exons.tsv" \
   >> "$OUT/KATNAL1_mDau_ERR13757851_STRG.3154.1.fa"

   echo \
   >> "$OUT/KATNAL1_mDau_ERR13757851_STRG.3154.1.fa"


Comprobar la longitud reconstruida:

.. code-block:: bash

   awk '
   !/^>/ {
       n += length($0)
   }
   END {
       print "KATNAL1 transcript length =", n, "nt"
   }
   ' 
   "$OUT/KATNAL1_mDau_ERR13757851_STRG.3154.1.fa"


La longitud esperada es:

::

   5804 nt


El exon terminal 3' de KATNAL1 es:

::

   NC_081841:104694088-104697084


Este exon coincide exactamente con la copia de ``hAT1_mMyo`` encontrada durante
el análisis genómico:

::

   hAT1_mMyo
   NC_081841
   104694088-104697084
   identity = 100 %
   query coverage = 100 %


IX. Caso FBXO3: definir exactamente sus exones codificantes
===========================================================

El transcrito de referencia es:

::

   Rhinolophus ferrumequinum
   SRR18652180
   STRG.3111.1


Su estructura exónica es:

::

   exon 1 = transcript 1-119
   exon 2 = transcript 120-209
   exon 3 = transcript 210-373
   exon 4 = transcript 374-1725


El exon 4 es derivado de ``RhiSin-1.24``.

La CDS quimérica es:

::

   transcript 16-1674


Por tanto:

::

   exon 1:
       1-15      = 5' UTR
       16-119    = HOST CDS

   exon 2:
       completo  = HOST CDS

   exon 3:
       completo  = HOST CDS

   exon 4:
       1-1301    = TE CDS
       1302-1352 = 3' UTR


La reconstrucción es:

::

   HOST CDS =
       exon1[16-119]
       +
       exon2
       +
       exon3

   TE CDS =
       exon4[1-1301]

   CHIMERA CDS =
       HOST CDS
       +
       TE CDS


X. Crear los probes exónicos de FBXO3
=====================================

Los cuatro exones de ``STRG.3111.1`` se guardaron como FASTA.

.. code-block:: bash

   mkdir -p 
   ANALYSES/dNdS_FINAL/reference_probes/FBXO3


.. code-block:: bash

   cat > 
   ANALYSES/dNdS_FINAL/reference_probes/FBXO3/FBXO3_rFer_STRG.3111.1_exons.fa 
   <<'EOF'
   >FBXO3_exon1|HOST|CDS=16-119
   ACCCGGCCAAGCAAGATGGCGGCCATGGATGTGGTGTCGGCACCTCTGACCCTAGAATTGTTGCCCACCGATACCCTGCTCCTCATCTTATCCTTCTTGGACTACCAGGACCTAATCAA
   >FBXO3_exon2|HOST|CDS=1-90
   TTGTTGCTATGTCAGTCGAAGACTTCACCAGCTCGCAAGCCATGATCCACTATGGAGAAGACATTGCAAAAAATACTGGCTGATATCTGA
   >FBXO3_exon3|HOST|CDS=1-164
   GGAAGAGAAAGCACAGAAGAATCAGTGCTGGAAATCTCTCTTCATTGATACTTACTCTGATGTAGGAAGATACATTGACCATTATGCTGCTATTAAAAAGGCCTGGGCCGATCTCAGGAAATATTTGGAGCCCAGACCTCCTCAGATGGCTTTGTCTCTGCAAG
   >FBXO3_exon4|TE|CDS=1-1301
   ACAAGATGAGTGCAAAAAGAAAGAGCTATTCCGTCGAGTACAAGAAAGGAATCGTGGAGGACTCCCGGGGCAAGAATCTTACGGCTTTCTGCAAAGAGAAGAAGTTGGATCTCCGAATGGTCCGAAAATGGAGAGCAGAGTACGATAACCTCAGTCAACAGGTGGACGAGGGAAATGCTAAGAAGCGCAAGTGTGGATCAGGTCGGCAACCATTATTTCCTGAGCTGGAAGACATCATCTGTGAATGGATTGCTGACAGGAGAGCAAAGGCTTTGGCTGTGCGCAGGGCTGATATTCAAGCATTTGCCCTTGCAATGGCACCACAGTTAGAAATATCCCCAGAAGAATTCAAAGCATCACAACACTGGCTAGATGGCTTCCTTCAGCGATACGAACTGTCTCTAAGATCGACAACACTGTTCAAGCTAGAAGATACTGAAGTTATTAAACGTGCACTTGCATTCAAGTCCTTTGTTGATGGCATTGACTTTTCTAAATACCAACTCTCCAACATGATTGCTATGGATGAAACTGCAGTGTTTATGGGCCAAGGATCTCAAACGACAGTTGATCATAGGGGTGCCTCGTCAATCTACATTCCCTCCACTGGTTACGAAAGTGCACGTGTTACCTGTATTTTGGCAATTCGTCTGGATGGAAAGAAAGCCCCACCTCTAATCATCACTAAGGGCAAGAACGATAAGGTTGAATGTGTTTCAGGCATTTATGTTCTTGAAACCGAAAAAGTTTGGTGCACACAAGCAGTTATAAGGAAGTGGGTCGATTTAATGCTGCCACTTGTTTTGCGAGGTGGCCAAAGAGGTCTGCTAGTCTGGGATTCAGCCAGCACTCACCGCGCTAAAGACATAAAGAACTTCCTTGCAGAGAGAAGAATAGATCAAATAATGATTCCCGCAGGAATGACTGCTTATCTCCAGGCTCTTGATGTTGCAATAAACAAGCCATTCAAGGATCATTTGCGCATGGAAATCAATGACTACATTGAAAATAGCATGGAGAGAAATCAGCACAGAAACTTTGTGAAGCCTAGCCTGCAAGAGGTCGTGACTTGGGTGAAGAATTCATGGGATAAAATCACTGACAGCTGTGTTGCCAATGCACTACGAGCAGGCTACATGGACAGGAAGTGCTCATTTAAGGAAAGCTCTATTGCTGGACATGAGAGGTTGGGGCCAATGGTTCTACAGGAAATGGAGTCGCAAGAAATTCAAGCCGGAATTCAGGGTTTGGAGAGTTATGACGATGTTCCAGAAGAAGATGACATGACTGTATTTGAATAAATGTAGATTTGAATGCCAGCATTCTGAATAAACCTGCTTTGCTTTCCACCA
   EOF


Comprobar las longitudes:

.. code-block:: bash

   awk '
   /^>/ {
       if (seq != "") {
           print name, length(seq)
       }

       name = $0
       seq = ""
       next
   }

   {
       seq = seq $0
   }

   END {
       print name, length(seq)
   }
   ' 
   ANALYSES/dNdS_FINAL/reference_probes/FBXO3/FBXO3_rFer_STRG.3111.1_exons.fa


Resultado esperado:

::

   exon1 = 119 nt
   exon2 = 90 nt
   exon3 = 164 nt
   exon4 = 1352 nt


XI. Validar los probes de FBXO3 contra la especie de referencia
===============================================================

Antes de utilizar los probes en otras especies se comprobó que recuperan exactamente
el locus conocido de ``Rhinolophus ferrumequinum``.

Crear directorio:

.. code-block:: bash

   mkdir -p 
   ANALYSES/dNdS_FINAL/test_FBXO3_rFer


Crear una base BLAST del scaffold de referencia:

.. code-block:: bash

   makeblastdb 
   -in ANALYSES/FBXO3/sequences/canonical_scaffolds/Rhinolophus_ferrumequinum.fa 
   -dbtype nucl 
   -out ANALYSES/dNdS_FINAL/test_FBXO3_rFer/rFer


Ejecutar BLAST:

.. code-block:: bash

   blastn 
   -query ANALYSES/dNdS_FINAL/reference_probes/FBXO3/FBXO3_rFer_STRG.3111.1_exons.fa 
   -db ANALYSES/dNdS_FINAL/test_FBXO3_rFer/rFer 
   -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovs' 
   -out ANALYSES/dNdS_FINAL/test_FBXO3_rFer/FBXO3_exons_vs_rFer.tsv


Inspeccionar:

.. code-block:: bash

   column -t 
   ANALYSES/dNdS_FINAL/test_FBXO3_rFer/FBXO3_exons_vs_rFer.tsv 
   | head -30


Los cuatro exones recuperaron exactamente el locus de referencia:

::

   exon1
   66442155-66442037

   exon2
   66437824-66437735

   exon3
   66435109-66434946

   exon4
   66433426-66432075


Todos los exones presentaron:

::

   query coverage = 100 %


en la copia de referencia.


XII. BLAST de los exones de FBXO3 en las seis especies finales
===============================================================

Aunque el locus paralógico ya había sido identificado previamente en
``03_FBXO3.rst``, los cuatro exones se utilizaron adicionalmente para recuperar sus
secuencias exactas.

Se creó:

::

   ANALYSES/dNdS_FINAL/scripts/03_blast_FBXO3_final_species.sh


Contenido:

.. code-block:: bash

   #!/bin/bash

   set -euo pipefail

   ROOT="/Volumes/Expansion/project3/bat_HTF_genomic_analysis"

   QUERY="$ROOT/ANALYSES/dNdS_FINAL/reference_probes/FBXO3/FBXO3_rFer_STRG.3111.1_exons.fa"

   OUT="$ROOT/ANALYSES/dNdS_FINAL/FBXO3_blast"

   species=(
       Rhinolophus_affinis
       Rhinolophus_ferrumequinum
       Rhinolophus_hipposideros
       Rhinolophus_pearsonii
       Rhinolophus_sinicus
       Rhinolophus_trifoliatus
   )

   mkdir -p "$OUT"

   for sp in "${species[@]}"
   do

       echo "===== $sp ====="

       genome="$ROOT/ANALYSES/FBXO3/sequences/canonical_scaffolds/${sp}.fa"

       db="$OUT/${sp}"

       makeblastdb \
       -in "$genome" \
       -dbtype nucl \
       -out "$db" \
       >/dev/null

       blastn \
       -query "$QUERY" \
       -db "$db" \
       -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovs' \
       -out "$OUT/${sp}_FBXO3_exons.tsv"

   done

   echo
   echo "DONE"


Ejecutar:

.. code-block:: bash

   chmod +x 
   ANALYSES/dNdS_FINAL/scripts/03_blast_FBXO3_final_species.sh

   bash 
   ANALYSES/dNdS_FINAL/scripts/03_blast_FBXO3_final_species.sh


XIII. No seleccionar simplemente el primer hit BLAST
====================================================

El análisis mostró inmediatamente que escoger el hit con mayor bitscore sería
incorrecto.

FBXO3 posee:

* una copia canónica;
* una copia paralógica quimérica.

Por tanto, algunos exones producen dos hits muy similares en el mismo scaffold.

Además, el exon 4 contiene secuencia derivada de un TE y puede producir numerosos
hits adicionales en el genoma.

Por esta razón, los resultados de BLAST deben cruzarse con las coordenadas del
parálogo identificadas previamente en el análisis de FBXO3.


Las coordenadas diagnósticas del análisis previo son:

.. list-table::
   :header-rows: 1

   * - Especie
     - Strand
     - Exon 2 / Block A
     - Exon 3 / Block B
   * - R. affinis
     - +
     - 24009684--24009773
     - 24012407--24012570
   * - R. ferrumequinum
     - -
     - 66437735--66437824
     - 66434946--66435109
   * - R. hipposideros
     - +
     - 23138058--23138147
     - 23140786--23140949
   * - R. pearsonii
     - +
     - 23305777--23305866
     - 23308274--23308437
   * - R. sinicus
     - +
     - 23096298--23096387
     - 23099038--23099201
   * - R. trifoliatus
     - -
     - 169436228--169436317
     - 169432523--169432684


Estos loci ya habían sido establecidos antes del análisis de ``dN/dS``.

Por tanto:

**el BLAST actual no se utiliza para volver a descubrir FBXO3**.

Se utiliza para recuperar las secuencias codificantes correspondientes al locus que
ya conocemos.


XIV. Inspeccionar solamente los hits cercanos al parálogo FBXO3
===============================================================

Para evitar mezclar hits de la copia canónica con hits del parálogo se utilizó una
ventana genómica alrededor del locus previamente identificado.

.. code-block:: python

   from pathlib import Path
   import csv

   root = Path(
       "ANALYSES/dNdS_FINAL/FBXO3_blast"
   )

   regions = {

       "Rhinolophus_affinis":
           (24006000, 24017000),

       "Rhinolophus_ferrumequinum":
           (66430000, 66443000),

       "Rhinolophus_hipposideros":
           (23132000, 23145000),

       "Rhinolophus_pearsonii":
           (23300000, 23313000),

       "Rhinolophus_sinicus":
           (23090000, 23104000),

       "Rhinolophus_trifoliatus":
           (169420000, 169440000),
   }

   for sp, (lo, hi) in regions.items():

       f = (
           root
           / f"{sp}_FBXO3_exons.tsv"
       )

       print()
       print(
           f"===== {sp} ====="
       )

       rows = []

       with open(f) as fh:

           reader = csv.reader(
               fh,
               delimiter="\t"
           )

           for x in reader:

               s1 = int(x[7])
               s2 = int(x[8])

               a = min(s1, s2)
               b = max(s1, s2)

               if (
                   a <= hi
                   and
                   b >= lo
               ):

                   rows.append(x)

       rows.sort(
           key=lambda x: (
               x[0],
               -float(x[10])
           )
       )

       for x in rows:

           print(
               x[0],
               x[1],
               x[2],
               x[5],
               x[6],
               x[7],
               x[8],
               sep="\t"
           )


Este código puede ejecutarse directamente desde Bash mediante:

.. code-block:: bash

   python3 - <<'PY'

   from pathlib import Path
   import csv

   root = Path(
       "ANALYSES/dNdS_FINAL/FBXO3_blast"
   )

   regions = {

       "Rhinolophus_affinis":
           (24006000, 24017000),

       "Rhinolophus_ferrumequinum":
           (66430000, 66443000),

       "Rhinolophus_hipposideros":
           (23132000, 23145000),

       "Rhinolophus_pearsonii":
           (23300000, 23313000),

       "Rhinolophus_sinicus":
           (23090000, 23104000),

       "Rhinolophus_trifoliatus":
           (169420000, 169440000),
   }

   for sp, (lo, hi) in regions.items():

       f = (
           root
           / f"{sp}_FBXO3_exons.tsv"
       )

       print()
       print(
           f"===== {sp} ====="
       )

       rows = []

       with open(f) as fh:

           reader = csv.reader(
               fh,
               delimiter="\t"
           )

           for x in reader:

               s1 = int(x[7])
               s2 = int(x[8])

               a = min(s1, s2)
               b = max(s1, s2)

               if (
                   a <= hi
                   and
                   b >= lo
               ):

                   rows.append(x)

       rows.sort(
           key=lambda x: (
               x[0],
               -float(x[10])
           )
       )

       for x in rows:

           print(
               x[0],
               x[1],
               x[2],
               x[5],
               x[6],
               x[7],
               x[8],
               sep="\t"
           )

   PY


XV. Estado actual de FBXO3
==========================

Los exones 2 y 3 ya están definidos de manera independiente por el análisis original
de FBXO3.

El exon 4 derivado del TE también puede localizarse utilizando el locus quimérico y
la posición de ``RhiSin-1.24``.

El punto que todavía debe resolverse cuidadosamente es el exon 1 en algunas especies.

Actualmente:

::

   Rhinolophus ferrumequinum
       exon1 identificado

   Rhinolophus hipposideros
       exon1 identificado

   Rhinolophus affinis
       candidato de exon1 identificado

   Rhinolophus pearsonii
       candidato de exon1 identificado

   Rhinolophus sinicus
       candidato de exon1 identificado

   Rhinolophus trifoliatus
       exon1 no recuperado con el BLAST actual


La ausencia de un hit de exon 1 en ``R. trifoliatus`` **no demuestra todavía que el
exon esté ausente**.

Puede representar:

* divergencia suficiente para impedir su recuperación con el probe actual;
* degradación del exon;
* pérdida real del exon.

Por ahora:

::

   no se sabe


y no se debe reconstruir artificialmente.


XVI. Archivos de entrada disponibles para las otras arquitecturas
=================================================================

Los análisis genómicos anteriores ya produjeron archivos que contienen los hits
exónicos y de TE.

Estos archivos deben reutilizarse.


TMPO
----

Ejemplo:

::

   ANALYSES/TMPO/intermediate/Myotis_myotis/
       TMPO_exon_matches_Myotis_myotis.tsv
       TMPO_TE_matches_Myotis_myotis.tsv


El TE correcto es:

::

   SPIN_MI


Ejemplo de hit principal:

::

   SPIN_MI_match
   scaffold_m19_p_21
   qcov = 100 %


NCL
---

Los archivos tienen la forma:

::

   NCL_exon_matches_<species>.tsv
   NCL_CHARLIE8_matches_<species>.tsv


El TE correspondiente es:

::

   CHARLIE8


ZKSCAN7
-------

Los archivos tienen la forma:

::

   ZKSCAN7_exon_matches_<species>.tsv
   ZKSCAN7_RhiSin_matches_<species>.tsv


El TE correspondiente es:

::

   RhiSin-1.24


KATNAL1
--------

Los archivos tienen la forma:

::

   KATNAL1_exon_matches_<species>.tsv
   KATNAL1_hAT1_matches_<species>.tsv


El TE correspondiente es:

::

   hAT1_mMyo


ZFP28
-----

Los archivos tienen la forma:

::

   ZFP28_exon_matches_<species>.tsv
   ZFP28_GOLEM_matches_<species>.tsv


El TE correspondiente es:

::

   GOLEM-1p


XVII. Principio para la reconstrucción multi-especie
====================================================

Para cada arquitectura, el procedimiento será:

::

   1. leer el manifiesto final de especies

   2. leer el locus ya identificado en el análisis genómico

   3. recuperar solamente los hits de ese locus

   4. determinar qué porción de cada exon corresponde a CDS

   5. extraer las secuencias genómicas

   6. orientarlas 5' → 3'

   7. concatenar únicamente segmentos codificantes

   8. generar:

       HOST CDS
       TE CDS
       CHIMERA CDS

   9. traducir

   10. realizar QC


XVIII. Estructura de salida
===========================

Para cada arquitectura se generará:

::

   ANALYSES/dNdS_FINAL/sequences/FBXO3/
   ANALYSES/dNdS_FINAL/sequences/TMPO/
   ANALYSES/dNdS_FINAL/sequences/NCL/
   ANALYSES/dNdS_FINAL/sequences/ZKSCAN7/
   ANALYSES/dNdS_FINAL/sequences/KATNAL1/
   ANALYSES/dNdS_FINAL/sequences/ZFP28/


Cada directorio contendrá:

::

   host_cds.all.fa
   host_cds.usable.fa

   te_cds.all.fa
   te_cds.usable.fa

   chimera_cds.all.fa
   chimera_cds.usable.fa

   host_protein.all.fa
   host_protein.usable.fa

   te_protein.all.fa
   te_protein.usable.fa

   chimera_protein.all.fa
   chimera_protein.usable.fa


Los archivos ``all`` conservarán las reconstrucciones literales.

Los archivos ``usable`` contendrán solamente las secuencias que pasen QC.


XIX. Control de calidad
=======================

Cada componente debe evaluarse independientemente.

Los estados de QC previstos son:

::

   PASS

   AMBIGUOUS_BASES

   INCOMPLETE_CODING_REGION

   MISSING_HOST_SEGMENT

   MISSING_TE_CODING_REGION

   FRAMESHIFT

   INTERNAL_STOP

   AMBIGUOUS_LOCUS_ASSIGNMENT

   UNCERTAIN_EXON_BOUNDARY

   MULTIPLE_EQUIVALENT_HITS

   ORIENTATION_CONFLICT

   REFERENCE_CDS_UNAVAILABLE


Una especie puede pasar QC para ``HOST`` y fallar para ``TE``.

Por tanto, una región TE degradada no debe eliminar automáticamente una secuencia
HOST válida.


El reporte final tendrá como mínimo:

::

   Architecture
   Species
   Component
   Table_status
   RNA_support
   Genome
   Scaffold
   Strand
   Segments_expected
   Segments_recovered
   Coding_nt_expected
   Coding_nt_recovered
   Ambiguous_bases
   Length_mod3
   Internal_stops
   Start_codon
   Terminal_stop
   QC_status
   Usable_for_comparison
   Exclusion_reason


XX. Regla para pasar al análisis de dN/dS
=========================================

Todavía no se realizará ningún análisis de selección.

Primero se generará:

::

   ANALYSES/dNdS_FINAL/results/usable_sequence_counts.tsv


con:

::

   Architecture
   Component
   n_table
   n_reconstructed
   n_usable
   prepare_for_selection


Solamente se continuará cuando:

::

   n_usable >= 3


El siguiente paso, después del QC, será:

::

   CDS reconstruidas
          |
          v
   traducción
          |
          v
   alineamiento de proteínas
          |
          v
   retrotraducción / alineamiento codon-aware
          |
          v
   comprobación de frame
          |
          v
   análisis dN/dS


En otras palabras:

**primero reconstruimos correctamente las secuencias; después analizamos selección.**
