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

.. code-block:: bash

    blastn \
    -query sequences/NCL_host_exon_probes.fa \
    -db intermediate/blast_db/Myotis_myotis \
    -task blastn \
    -word_size 7 \
    -dust no \
    -evalue 1e-10 \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
    > intermediate/NCL_exon_matches_Myotis_myotis.tsv

.. code-block:: bash

    blastn \
    -query sequences/NCL_Charlie8_probe.fa \
    -db intermediate/blast_db/Myotis_myotis \
    -task blastn \
    -word_size 7 \
    -dust no \
    -evalue 1e-20 \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
    > intermediate/NCL_Charlie8_matches_Myotis_myotis.tsv

.. code-block:: bash

    cat \
    intermediate/NCL_exon_matches_Myotis_myotis.tsv \
    intermediate/NCL_Charlie8_matches_Myotis_myotis.tsv \
    > intermediate/NCL_all_matches_Myotis_myotis.tsv

.. code-block:: bash

    awk 'BEGIN{OFS="\t"}
    {
        pos = ($8 < $9 ? $8 : $9)
        print $0, pos
    }' intermediate/NCL_all_matches_Myotis_myotis.tsv \
    | sort -k2,2 -k13,13n \
    | cut -f1-12 \
    | column -t


4) Automate detection of candidate NCL-Charlie8 loci for all species
--------------------------------------------------------------------

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/NCL

.. code-block:: bash

    grep -h -v '^Species' results/NCL_signature_*.tsv \
    | awk -F'\t' '{count[$NF]++} END {for (x in count) print x, count[x]}'

.. code-block:: bash

    grep -h 'NCL_CHIMERIC_CANDIDATE' \
    results/NCL_signature_*.tsv \
    | cut -f1-6,7,14,15



el resultado es este:

Los 23 candidatos fuertes fueron:
- Antrozous pallidus
- Corynorhinus mexicanus
- Corynorhinus townsendii
- Eptesicus fuscus
- Eptesicus nilssonii
- Lasiurus ega
- Myotis auriculus
- Myotis californicus
- Myotis daubentonii
- Myotis evotis
- Myotis lucifugus
- Myotis myotis
- Myotis mystacinus
- Myotis nigricans
- Myotis occultus
- Myotis pilosus
- Myotis thysanodes
- Myotis velifer
- Myotis vivesi
- Myotis volans
- Myotis yumanensis
- Pipistrellus kuhlii (pero hay varias copias de NCL no es el mejor)
- Plecotus auritus

    
II. AUDITORÍA REPRODUCIBLE DE NCL–CHARLIE8
===========================================

Propósito y alcance
-------------------

Esta sección reemplaza el registro preliminar de la auditoría. Documenta en orden
cronológico por qué se revisó la detección de NCL, cómo se verificó el control
*Vespertilio murinus*, cómo se extendió la búsqueda a los genomas disponibles y
cómo se seleccionaron las 28 especies candidatas para evaluar CHARLIE8.

**Estado al cierre de esta auditoría:** se generó
``results/NCL_Charlie8_full_length_coordinates.tsv`` con **28 filas**.
La extracción de los 28 FASTA y el análisis de codones de terminación **aún no
se habían confirmado como ejecutados**. Tampoco se había validado la cobertura
de cada uno de los nueve exones en las 28 especies ni la expresión de todas las
quimeras. No interpretar «9/9» como nueve exones completos confirmados.

Convenciones utilizadas
~~~~~~~~~~~~~~~~~~~~~~~

* Los comandos del análisis comparativo se ejecutan desde
  ``/Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/NCL``.
  Los pasos iniciales en *V. murinus* indican explícitamente cuándo se trabaja
  desde la raíz del proyecto.
* Las coordenadas de BLAST y de ``samtools faidx`` son **1-based inclusivas**.
  Las coordenadas BED de RepeatMasker se interpretan como **0-based,
  intervalo semiabierto**; para llevarlas a 1-based inclusivas: inicio BED + 1,
  fin BED sin modificación.
* En la hebra negativa, la dirección biológica 3′ corresponde a coordenadas
  genómicas **decrecientes**. No se debe confundir «aguas abajo» con aumento
  numérico de la coordenada.
* El umbral máximo de separación usado en la búsqueda es **25.000 bp**.
  Es un criterio operativo de proximidad, no prueba de empalme o de
  incorporación del TE al transcrito.
* La sonda de CHARLIE8 mide **1.918 bp**. La cobertura BLAST expresa
  recuperación **respecto a esta sonda**, no respecto al transposón ancestral
  completo.

Paso 1. Identificar el falso negativo que motivó la auditoría
------------------------------------------------------------

**Problema.** El análisis inicial de firmas NCL/CHARLIE8 no recuperó
*Vespertilio murinus* como se esperaba. En esa especie, el locus NCL candidato
está en la hebra negativa. Una lectura incorrecta de la orientación 3′ o una
proyección incompleta de exones podía explicar la ausencia.

**Objetivo.** Auditar independientemente las clasificaciones previas sin
presuponer que las especies sin anotación de CHARLIE8 son verdaderos negativos.

Desde la raíz del proyecto:

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
    ls ANALYSES/NCL/results/NCL_signature_*.tsv

El primer intento fue una auditoría basada en firmas NCL ya clasificadas y
anotaciones BED de RepeatMasker. Se ejecutó con el archivo
``ANALYSES/NCL/auditar_NCL.py``. El código original se conserva íntegramente
más abajo en el **Apéndice A**, para que el intento inicial sea rastreable.

.. code-block:: bash

    cd ANALYSES/NCL
    python3 auditar_NCL.py --root ../..

**Nota de reproducibilidad:** desde ``ANALYSES/NCL``, la raíz del proyecto está
**dos**, no tres, niveles arriba. El comando histórico documentado como
``--root ../../..`` no apunta a la raíz correcta si se ejecuta desde ese
directorio. La auditoría BED produce:

* ``results/NCL_CHARLIE8_independent_audit.tsv``
* ``results/NCL_CHARLIE8_IGV_checklist.tsv``

La auditoría clasifica proximidad 3′/5′ y solapamiento, pero **no demuestra**
que el TE constituya un exón ni identifica transcritos quiméricos. Además,
algunos genomas carecían de anotaciones RepeatMasker BED utilizables: ausencia
de anotación no equivale a ausencia del elemento.

Paso 2. Inspeccionar el paralogo NCL de Vespertilio murinus
---------------------------------------------------------

**Por qué.** La proyección de TOGA ``ENST00000678364.1#NCL#786044`` contenía
**ocho bloques**, de modo que era necesario comprobar el locus en secuencia
genómica sin depender de esa proyección.

Desde la **raíz** del proyecto se extrajeron ocho bloques previamente
registrados, respetando la hebra negativa:

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis
    mkdir -p ANALYSES/NCL/sequences/Vespertilio_murinus
    GENOME="SPECIES/Vespertilio_murinus/mVesMur1.1.pri.fa"
    OUT="ANALYSES/NCL/sequences/Vespertilio_murinus/NCL_786044_blocks.fa"
    : > "$OUT"

    while read -r exon start end; do
        seq=$(samtools faidx "$GENOME" "OZ004711:${start}-${end}" \
              | tail -n +2 | tr -d '\n' \
              | tr 'ACGTacgt' 'TGCAtgca' | rev)
        printf ">Vespertilio_murinus_NCL_block%s\n%s\n" \
               "$exon" "$seq" >> "$OUT"
    done <<'EOF'
    1 64537101 64537118
    2 64534976 64535083
    3 64532829 64533132
    4 64532275 64532445
    5 64531425 64531511
    6 64530361 64530499
    7 64528873 64528997
    8 64528112 64528235
    EOF

Esas son las coordenadas de los **bloques de la reconstrucción inicial**, no
las coordenadas finales de todos los alineamientos de exones encontrados por
BLAST. El resultado de ocho bloques **no** podía demostrar que faltara el
noveno exón.

Paso 3. Comprobar las nueve sondas NCL y buscar independientemente todos los exones
---------------------------------------------------------------------------------

**Por qué.** Había que determinar si el exón 9 estaba realmente ausente o si
TOGA lo había omitido.

Desde la raíz:

.. code-block:: bash

    grep '^>' ANALYSES/NCL/sequences/NCL_host_exon_probes.fa

El archivo contiene las nueve consultas ``NCL_exon1_match`` hasta
``NCL_exon9_match``. Se extrajo una región suficientemente amplia del scaffold
``OZ004711`` para buscar los exones de manera independiente:

.. code-block:: bash

    samtools faidx \
      SPECIES/Vespertilio_murinus/mVesMur1.1.pri.fa \
      OZ004711:64523000-64542000 \
      > ANALYSES/NCL/sequences/Vespertilio_murinus/NCL_786044_extended.fa

    blastn \
      -query ANALYSES/NCL/sequences/NCL_host_exon_probes.fa \
      -subject ANALYSES/NCL/sequences/Vespertilio_murinus/NCL_786044_extended.fa \
      -task blastn -word_size 7 -dust no -evalue 1e-3 \
      -outfmt '6 qseqid pident length qlen qstart qend sstart send evalue bitscore'

**Resultado observado:** BLAST recuperó **los nueve exones**, incluido
``NCL_exon9_match`` (157 de 158 posiciones de consulta, 98,089 % de identidad).
La ausencia del exón 9 en aquella proyección de TOGA era, por tanto, un
**falso negativo de la anotación/proyección inicial**, no evidencia de ausencia
genómica. Esto justificó volver a examinar el conjunto de especies usando
búsquedas independientes.

Paso 4. Localizar CHARLIE8 en dirección 3′ del exón 9
---------------------------------------------------

**Por qué.** La arquitectura buscada es ``NCL exones 1–9 → CHARLIE8`` en el
orden **transcripcional**. Para *V. murinus*, ese orden corresponde a posiciones
genómicas decrecientes.

La anotación local de RepeatMasker situaba una región CHARLIE8 en
``OZ004711:64521164-64523425`` y el exón 9 de NCL en
``OZ004711:64524093-64524249``. Se extrajeron, sobre la hebra negativa,
para inspección y comparación:

.. code-block:: bash

    GENOME="SPECIES/Vespertilio_murinus/mVesMur1.1.pri.fa"
    OUTDIR="ANALYSES/NCL/sequences/Vespertilio_murinus"

    samtools faidx -i "$GENOME" OZ004711:64521164-64523425 \
      | awk 'NR==1 {$0=">Vespertilio_murinus_CHARLIE8"} {print}' \
      > "$OUTDIR/CHARLIE8.fa"

    samtools faidx -i "$GENOME" OZ004711:64524093-64524249 \
      | awk 'NR==1 {$0=">Vespertilio_murinus_NCL_exon9"} {print}' \
      > "$OUTDIR/NCL_exon9.fa"

    samtools faidx -i "$GENOME" OZ004711:64523426-64524092 \
      | awk 'NR==1 {$0=">NCL_exon9_CHARLIE8_intervening_667bp"} {print}' \
      > "$OUTDIR/NCL_exon9_CHARLIE8_intervening.fa"

Los 667 bp corresponden a la **separación entre el intervalo RepeatMasker**
y el exón 9; no equivalen a la separación respecto al alineamiento de la
sonda de 1.918 bp. Este último alineamiento va de
``64521191–64523093``, y la separación BLAST respecto al exón 9 es
``64524093 - 64523093 - 1 = 999 bp``.

Para revisar anotaciones de exones sobre esa región se empleó también:

.. code-block:: bash

    awk -F '\t' '
    BEGIN {OFS="\t"}
    $1=="OZ004711" && $3=="exon" &&
    $4<=64537118 && $5>=64521164 {
        print FILENAME,$4,$5,$7,$9
    }' /Volumes/Expansion/IGV/Vespertilio_murinus/ERR1375784*_MM.gtf

**Evidencia adicional observada en el control:** se inspeccionaron
alineamientos RNA-seq con la unión entre el exón 9 y el exón derivado de
CHARLIE8, incluyendo un intrón anotado en los alineamientos de **999 nt**.
El registro de trabajo contabilizó 27 alineamientos que respaldaban esa unión
entre cinco BAM (4, 8, 2, 6 y 7, respectivamente). Esta observación se limita
al control y no debe extrapolarse a las otras 27 especies.

Paso 5. Preparar la búsqueda de secuencia independiente en los genomas disponibles
--------------------------------------------------------------------------------

**Por qué.** El fallo de TOGA y la cobertura incompleta de archivos BED
obligaban a buscar directamente las sondas de NCL y CHARLIE8 en los FASTA
genómicos.

A partir de aquí, trabajar desde **ANALYSES/NCL**:

.. code-block:: bash

    cd /Volumes/Expansion/project3/bat_HTF_genomic_analysis/ANALYSES/NCL
    mkdir -p intermediate/NCL_103_blast intermediate/NCL_Charlie8_blast \
             intermediate/blast_db results

    wc -l intermediate/NCL_genomes_to_search.txt
    grep '^>' sequences/NCL_host_exon_probes.fa
    grep '^>' sequences/NCL_Charlie8_probe.fa

**Entradas de esta etapa:**

* ``intermediate/NCL_genomes_to_search.txt``: listado de **101 rutas de FASTA
  seleccionados** previamente; sus rutas deben existir antes de ejecutar BLAST.
* ``sequences/NCL_host_exon_probes.fa``: nueve secuencias de consulta para exones.
* ``sequences/NCL_Charlie8_probe.fa``: sonda CHARLIE8 de **1.918 bp**.
* ``SPECIES/<Species>/``: directorios de ensamblajes genómicos.

Había **103 directorios de especies**, pero **101 genomas** disponibles para esta
búsqueda. Los directorios ``Hipposideros_cyclops`` y
``Rhinolophus_sedulus`` no aportaban los FASTA genómicos esperados; sus datos
TOGA estaban asociados a otros nombres taxonómicos. Para
``Miniopterus_natalensis`` se seleccionó el ensamblaje hap2. Las dos especies
sin FASTA no deben contarse como negativas.

**Dependencia:** el procedimiento de selección que produjo
``NCL_genomes_to_search.txt`` pertenece a una etapa anterior. Esta sección
conserva ese listado **como entrada fija**; no inventa un algoritmo nuevo de
selección de ensamblajes.

Paso 6. Ejecutar BLAST genómico para NCL y CHARLIE8
--------------------------------------------------

**Por qué.** Se buscó cada sonda por separado en los 101 genomas,
independientemente de TOGA y RepeatMasker. El siguiente bloque reconstruye
el procedimiento reproducible **a partir del listado de FASTA ya seleccionado**.
Si existen los TSV originales y no se desean sobreescribir, omitir este paso.

.. code-block:: bash

    while IFS= read -r genome; do
        [ -n "$genome" ] || continue
        [ -f "$genome" ] || { echo "Falta FASTA: $genome" >&2; exit 1; }
        species=$(basename "$(dirname "$genome")")
        db="intermediate/blast_db/$species"

        makeblastdb -in "$genome" -dbtype nucl -out "$db" >/dev/null || exit 1

        blastn -query sequences/NCL_host_exon_probes.fa \
          -db "$db" -task blastn -word_size 7 -dust no -evalue 1e-3 \
          -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
          > "intermediate/NCL_103_blast/${species}.tsv" || exit 1

        blastn -query sequences/NCL_Charlie8_probe.fa \
          -db "$db" -task blastn -word_size 7 -dust no -evalue 1e-3 \
          -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp' \
          > "intermediate/NCL_Charlie8_blast/${species}.tsv" || exit 1
    done < intermediate/NCL_genomes_to_search.txt

Se obtuvieron **101 archivos TSV de exones NCL** y **101 archivos TSV de
CHARLIE8**. En cada uno se usaron estas 12 columnas:

``qseqid, sseqid, pident, length, qlen, qstart, qend, sstart, send, evalue, bitscore, qcovhsp``.

Los porcentajes de cobertura e identidad en lo que sigue proceden de estas
salidas BLAST. **``qcovhsp`` corresponde a cada HSP**, no necesariamente a la
cobertura acumulada de varios HSP del mismo locus.

Paso 7. Buscar los pares exón 9–CHARLIE8 en orientación 3′
--------------------------------------------------------

**Por qué.** Antes de examinar los nueve exones, se identificaron posiciones
genómicas en las que un hit de CHARLIE8 quedara **después del exón 9** en la
misma hebra y el mismo scaffold, a una distancia de **0–25.000 bp**.

.. code-block:: bash

    printf 'Species\tScaffold\tStrand\tNCL_exon9_start\tNCL_exon9_end\tCharlie8_start\tCharlie8_end\tDistance_bp\tCharlie8_identity\tCharlie8_qcov\n' \
      > results/NCL_Charlie8_3prime_candidates.tsv

    for ncl in intermediate/NCL_103_blast/*.tsv; do
        species=$(basename "$ncl" .tsv)
        te="intermediate/NCL_Charlie8_blast/${species}.tsv"

        awk -v species="$species" '
        BEGIN {FS=OFS="\t"}
        FNR==NR {
            if ($1!="NCL_exon9_match") next
            strand=($8<$9 ? "+" : "-")
            s=($8<$9 ? $8 : $9)
            e=($8>$9 ? $8 : $9)
            key=$2 SUBSEP strand
            count[key]++
            exonStart[key,count[key]]=s
            exonEnd[key,count[key]]=e
            next
        }
        {
            strand=($8<$9 ? "+" : "-")
            s=($8<$9 ? $8 : $9)
            e=($8>$9 ? $8 : $9)
            key=$2 SUBSEP strand
            for (i=1;i<=count[key];i++) {
                if (strand=="+") distance=s-exonEnd[key,i]-1
                else distance=exonStart[key,i]-e-1
                if (distance>=0 && distance<=25000)
                    print species,$2,strand,exonStart[key,i],exonEnd[key,i],s,e,distance,$3,$12
            }
        }' "$ncl" "$te" >> results/NCL_Charlie8_3prime_candidates.tsv
    done

    wc -l results/NCL_Charlie8_3prime_candidates.tsv
    tail -n +2 results/NCL_Charlie8_3prime_candidates.tsv \
      | cut -f1 | sort -u | wc -l

**Resultado observado:** 74 líneas con encabezado = **73 pares candidatos**
entre **54 especies**. No eran 73 arquitecturas independientes: un locus podía
generar varios hits de CHARLIE8.

Paso 8. Examinar por qué algunas especies tenían varios pares
------------------------------------------------------------

**Por qué.** Varias filas por especie podían reflejar HSP repetitivos o
parciales, no varios parálogos completos.

.. code-block:: bash

    awk -F '\t' '
    NR>1 {count[$1]++}
    END {for (species in count) if (count[species]>1)
             print species,count[species]}
    ' results/NCL_Charlie8_3prime_candidates.tsv | sort

Se observaron **13 especies** con más de un par: por ejemplo,
``Lasiurus_ega`` (4) y ``Rhinopoma_microphyllum`` (5). En esta etapa no se
asignó a estos pares el estatus de parálogos independientes.

Paso 9. Exigir los nueve alineamientos NCL en el orden correcto
--------------------------------------------------------------

**Por qué.** Tener un exón 9 cercano a CHARLIE8 no era suficiente. Se buscó
para cada par una cadena ``exón 1 → 2 → … → 9`` del mismo scaffold y hebra.
Se aplicó un límite operativo de **25 kb entre exones consecutivos**.

**Corrección del código durante la auditoría:** una primera versión extraía
erróneamente el número del exón desde ``a[1]`` (que contenía ``NCL``). Otra
versión usaba una coordenada inicial incorrecta. Ambas fueron descartadas.
**El código siguiente es la versión corregida que produjo los resultados
reportados.**

.. code-block:: bash

    printf 'Species\tScaffold\tStrand\tNCL_exons\tCharlie8_distance_bp\n' \
      > results/NCL_9exons_Charlie8_positive.tsv

    for ncl in intermediate/NCL_103_blast/*.tsv; do
        species=$(basename "$ncl" .tsv)

        awk -F '\t' -v species="$species" '
        BEGIN {OFS="\t"}
        FNR==NR {
            if ($1 !~ /^NCL_exon[1-9]_match$/) next
            split($1,a,"_")
            exon=a[2]
            sub(/^exon/,"",exon)
            n++
            chr[n]=$2
            dir[n]=($8<$9 ? "+" : "-")
            ex[n]=exon+0
            st[n]=($8<$9 ? $8 : $9)
            en[n]=($8>$9 ? $8 : $9)
            next
        }
        FNR==1 {next}
        $1!=species {next}
        {
            chrom=$2
            strand=$3
            pos=(strand=="+" ? $4 : $5)
            found=1
            for (exon=8; exon>=1; exon--) {
                best=0
                bestgap=25001
                for (i=1; i<=n; i++) {
                    if (chr[i]!=chrom || dir[i]!=strand || ex[i]!=exon) continue
                    gap=(strand=="+" ? pos-en[i]-1 : st[i]-pos-1)
                    if (gap>=0 && gap<bestgap) {
                        best=i
                        bestgap=gap
                    }
                }
                if (!best) {found=0; break}
                pos=(strand=="+" ? st[best] : en[best])
            }
            if (found) print species,chrom,strand,"9/9",$8
        }' "$ncl" results/NCL_Charlie8_3prime_candidates.tsv
    done | sort -u >> results/NCL_9exons_Charlie8_positive.tsv

    column -t -s $'\t' results/NCL_9exons_Charlie8_positive.tsv

    tail -n +2 results/NCL_9exons_Charlie8_positive.tsv \
      | cut -f1 | sort -u | wc -l
    tail -n +2 results/NCL_9exons_Charlie8_positive.tsv | wc -l

**Resultado observado:** **28 especies** y **32 pares candidatos**.
*Vespertilio murinus* se recuperó con ``OZ004711``, hebra ``-``, ``9/9`` y
separación de **999 bp** respecto al hit de CHARLIE8.

**Limitación crítica:** este filtro verifica que existan alineamientos de los
nueve exones con orden, hebra y distancias compatibles. **No comprueba todavía
cobertura mínima, completitud codificante, unicidad del locus ni empalme
RNA-seq**. La selección codifica «9/9 hits ordenados», no «9 exones íntegros».

Paso 10. Recuperar la identidad y cobertura de CHARLIE8 del BLAST original
------------------------------------------------------------------------

**Por qué.** El archivo del paso 7 ya conservaba las dos métricas necesarias:
``Charlie8_identity`` (columna 9) y ``Charlie8_qcov`` (columna 10). No era
necesario repetir BLAST ni reconstruir estos porcentajes.

Para asociar los pares del paso 9 con sus métricas, se cruzaron especie,
scaffold, hebra **y distancia exón 9–TE** (esta última era la columna 5 del
archivo de positivos y columna 8 del archivo de pares):

.. code-block:: bash

    awk -F '\t' '
    BEGIN {OFS="\t"}
    NR==FNR {
        if (FNR>1) positive[$1 FS $2 FS $3 FS $5]=1
        next
    }
    FNR==1 {next}
    ($1 FS $2 FS $3 FS $8) in positive {
        print $1,$2,$3,$8,$10,$9
    }' results/NCL_9exons_Charlie8_positive.tsv \
       results/NCL_Charlie8_3prime_candidates.tsv \
    | sort -k1,1 -k4,4n \
    > results/NCL_Charlie8_recovery_metrics.tsv

    column -t results/NCL_Charlie8_recovery_metrics.tsv

Las seis columnas del archivo de salida son, en orden: especie, scaffold,
hebra, distancia (bp), ``qcovhsp`` (%), identidad (%).

**Resultado observado:** 32 filas. De ellas, **28 HSP** cubrían **99–100 %**
de la sonda de CHARLIE8: **24 al 100 %** y **4 al 99 %**. Las cuatro
coincidencias adicionales de ``Antrozous_pallidus`` (una) y
``Lasiurus_ega`` (tres) tenían coberturas de 4, 5, 4 y 10 %, y no se
consideraron recuperaciones casi completas.

Identidad entre los 28 hits con cobertura ≥99 %: **86,018–100,000 %**.
El control ``Vespertilio_murinus`` presentó **100 % de cobertura y 95,151 %
de identidad** respecto a la sonda. Las coberturas de 99–100 % no garantizan
identidad de ORF, ausencia de indels ni ausencia de codones de terminación.

Paso 11. Crear la tabla de coordenadas de los 28 CHARLIE8 casi completos
------------------------------------------------------------------------

**Por qué.** Antes de extraer secuencias, era necesario conservar las
coordenadas genómicas y la orientación de **cada HSP seleccionado**. Estas
coordenadas no estaban en la tabla reducida del paso 10, pero sí en el
archivo original del paso 7.

.. code-block:: bash

    awk -F '\t' '
    BEGIN {OFS="\t"}
    NR==FNR {
        if (FNR>1) positive[$1 FS $2 FS $3 FS $5]=1
        next
    }
    FNR>1 && $10>=99 &&
    ($1 FS $2 FS $3 FS $8) in positive {
        print $1,$2,$3,$6,$7,$9,$10
    }' results/NCL_9exons_Charlie8_positive.tsv \
       results/NCL_Charlie8_3prime_candidates.tsv \
    > results/NCL_Charlie8_full_length_coordinates.tsv

    wc -l results/NCL_Charlie8_full_length_coordinates.tsv

**Resultado observado en la terminal:**

.. code-block:: text

    28 results/NCL_Charlie8_full_length_coordinates.tsv

Las columnas son: ``Species``, ``Scaffold``, ``Strand``, ``Charlie8_start``,
``Charlie8_end``, ``Charlie8_identity``, ``Charlie8_qcov``. El archivo no
incluye encabezado, para facilitar su lectura desde ``bash``.

Paso 12. Estado de la auditoría y análisis pendientes
----------------------------------------------------

Al cierre de la conversación documentada, **el último resultado ejecutado y
confirmado es el archivo de 28 coordenadas del paso 11**.

Las siguientes acciones quedaron propuestas, **pero NO se deben presentar
como completadas**:

#. Extraer, desde los FASTA de las 28 especies, las regiones de CHARLIE8
   usando ``samtools faidx`` y aplicando la orientación inversa complementaria
   cuando corresponda. Mantener la asociación exacta entre especie, scaffold
   y coordenadas; verificar que se obtengan 28 FASTA no vacíos.
#. Evaluar los posibles indels y la estructura de los alineamientos, pues
   el porcentaje de cobertura BLAST por sí solo no garantiza una secuencia
   codificante continua.
#. Determinar el marco de lectura **a partir de la CDS y la unión NCL–CHARLIE8**.
   Traducir CHARLIE8 aislado en los tres marcos no permite distinguir por sí
   solo un *stop* prematuro del *stop* terminal fisiológico.
#. Inspeccionar la cobertura de cada uno de los nueve exones para confirmar la
   integridad de la arquitectura, y evaluar el soporte de RNA para cada
   especie cuando haya datos disponibles.

**Conclusión restringida a lo observado:** se identificaron **28 especies
candidatas** que contienen nueve hits NCL ordenados y uno de CHARLIE8 con
cobertura ≥99 % de la sonda, en orientación 3′ y separados por ≤25 kb. **No se
sabe aún** cuántas tienen una ORF fusionada intacta ni cuántas contienen
codones de terminación prematuros.

Apéndice A. Código original de la auditoría basada en BED
--------------------------------------------------------

Se conserva aquí el primer script histórico, sin reinterpretar sus banderas
como confirmación transcriptómica. Su dependencia de los archivos
``NCL_signature_*.tsv`` y de las anotaciones BED implica que **no sustituye**
la auditoría independiente por secuencias descrita en los pasos 3–11.

.. code-block:: python


    #!/usr/bin/env python3
    """Auditoría independiente de NCL/Charlie8 en todas las especies analizadas.
    No clasifica expresión ni demuestra una fusión transcrita.
    Uso: python3 auditar_NCL_CHARLIE8_103.py --root /Volumes/Expansion/project3/bat_HTF_genomic_analysis
    """
    import argparse
    import csv
    import re
    from pathlib import Path
    from collections import defaultdict
    
    TE = re.compile(r'charlie[\s_\-]*8', re.I)
    
    
    def get_signature_files(root):
        return sorted((root/'ANALYSES'/'NCL'/'results').glob('NCL_signature_*.tsv'))
    
    
    def read_loci(files):
        loci=[]; problems=[]
        for fp in files:
            with fp.open(errors='replace') as f:
                for lineno,line in enumerate(f,1):
                    if not line.strip() or line.lower().startswith('species\t'): continue
                    c=line.rstrip('\n').split('\t')
                    if len(c)<7:
                        problems.append((str(fp),lineno,'<7 columns'));continue
                    species,assembly,chrom,strand=c[:4]
                    try: start,end=int(c[4]),int(c[5])
                    except ValueError:
                        problems.append((str(fp),lineno,'invalid coordinates'));continue
                    if strand not in ('+','-') or start>=end:
                        problems.append((str(fp),lineno,'invalid strand/interval'));continue
                    status=c[-1].strip()  # The original classifier is the last TSV field
                    loci.append(dict(species=species,assembly=assembly,chrom=chrom,strand=strand,
                                     start=start,end=end,exons=c[6],original_status=status,
                                     original_file=fp.name,original_line=lineno))
        return loci,problems
    
    
    def bed_candidates(root,species):
        d=root/'SPECIES'/species
        if not d.exists(): return []
        # RepeatMasker annotations in this project are often named e.g. mVesMur1.1.pri.bed.
        # Do NOT use TOGA query_annotation.bed as evidence for transposable elements.
        return sorted(p for p in d.glob('*.bed') if p.is_file()
                      and 'query_annotation' not in p.name.lower()
                      and 'toga' not in p.name.lower()
                      and not p.name.startswith('._'))
    
    
    def parse_charlie_bed(paths):
        hits=defaultdict(list); read_files=[]; suspect=[]
        for path in paths:
            n=0; found=0
            try:
                with path.open(errors='replace') as f:
                    for line in f:
                        if not line.strip() or line.startswith(('#','track ','browser ')):continue
                        c=line.rstrip('\n').split('\t')
                        if len(c)<4: continue
                        # BED requires columns chrom,start0,end0,name.
                        try: s,e=int(c[1]),int(c[2])
                        except ValueError: continue
                        if e<=s: continue
                        n+=1
                        # Search annotation fields, rather than positions or scaffold IDs.
                        if not any(TE.search(field) for field in c[3:]):continue
                        found+=1
                        hits[c[0]].append((s+1,e,c[3],path.name))  # to 1-based inclusive
            except (OSError,UnicodeError) as exc:
                suspect.append((path.name,str(exc)));continue
            read_files.append((path.name,n,found))
        return hits,read_files,suspect
    
    
    def relation(locus, hit, window):
        hs,he,_,_=hit
        st,en=locus['start'],locus['end']
        # Report 3' proximity strictly using transcriptional orientation.
        # Overlap of TE with locus is reported separately and needs exon-level inspection.
        if hs<=en and he>=st:
            return 'OVERLAPS_NCL_INTERVAL',0
        if locus['strand']=='+':
            if hs>en: return 'DOWNSTREAM_3PRIME',hs-en-1
            return 'UPSTREAM_5PRIME',st-he-1
        else:
            if he<st: return 'DOWNSTREAM_3PRIME',st-he-1
            return 'UPSTREAM_5PRIME',hs-en-1
    
    
    def main():
        ap=argparse.ArgumentParser()
        ap.add_argument('--root',type=Path,required=True)
        ap.add_argument('--window',type=int,default=25000,
                        help='maximum distance outside NCL locus in its biological 3-prime direction (bp)')
        args=ap.parse_args(); root=args.root.expanduser().resolve()
        source=get_signature_files(root)
        if not source: raise SystemExit('ERROR: no NCL_signature_*.tsv in ANALYSES/NCL/results')
        loci,problems=read_loci(source)
        if not loci: raise SystemExit('ERROR: signature tables yielded no loci')
        resultdir=root/'ANALYSES'/'NCL'/'results'
        out=resultdir/'NCL_CHARLIE8_independent_audit.tsv'
        checklist=resultdir/'NCL_CHARLIE8_IGV_checklist.tsv'
        all_rows=[]; missing=[]; read_summary=[]
        for species in sorted({d['species'] for d in loci}):
            bedfiles=bed_candidates(root,species)
            if not bedfiles:
                missing.append(species)
                hitmap={}; valid=[];suspect=[]
            else:
                hitmap,valid,suspect=parse_charlie_bed(bedfiles)
                read_summary.extend((species,*x) for x in valid)
            usable=any(n>0 for _,n,_ in valid)
            for locus in (x for x in loci if x['species']==species):
                chrhits=hitmap.get(locus['chrom'],[])
                # Compare all Charlie8 hits on this scaffold, distinguishing 3' vs 5'.
                nearby=[]
                for hit in chrhits:
                    direction,distance=relation(locus,hit,args.window)
                    if distance<=args.window:
                        nearby.append((hit,direction,distance))
                nearby.sort(key=lambda z:(z[2], 0 if z[1]=='DOWNSTREAM_3PRIME' else 1))
                supported=[x for x in nearby if x[1] in ('DOWNSTREAM_3PRIME','OVERLAPS_NCL_INTERVAL')]
                best=supported[0] if supported else (nearby[0] if nearby else None)
                if not usable: flag='ANNOTATION_MISSING_OR_UNREADABLE'
                elif supported: flag='REVIEW_3PRIME_OR_OVERLAP'
                elif nearby: flag='ONLY_5PRIME_NEARBY'
                else: flag='NO_CHARLIE8_WITHIN_25KB_OF_NCL' if args.window==25000 else 'NO_CHARLIE8_WITHIN_WINDOW'
                hit,side,gap=(best if best is not None else (None,'',''))
                all_rows.append({
                    'Species':species,'Scaffold':locus['chrom'],'Strand':locus['strand'],
                    'NCL_start':locus['start'],'NCL_end':locus['end'],
                    'NCL_exons_detected':locus['exons'],
                    'Original_classification':locus['original_status'],
                    'RepeatMasker_BEDs':';'.join(p.name for p in bedfiles),
                    'BED_has_valid_records':int(usable),
                    'Charlie8_hits_on_scaffold':len(chrhits),
                    'Charlie8_hits_within_window':len(nearby),
                    'Best_Charlie8_start':hit[0] if hit else '',
                    'Best_Charlie8_end':hit[1] if hit else '',
                    'Best_Charlie8_name':hit[2] if hit else '',
                    'Best_Charlie8_track':hit[3] if hit else '',
                    'Relative_position':side,
                    'Separation_bp':gap,
                    'Audit_flag':flag,
                    'IGV_region':f"{locus['chrom']}:{max(1, locus['start']-args.window):,}-{locus['end']+args.window:,}".replace(',',''),
                    'Caution':'Proximity does not prove TE exon inclusion or RNA splicing'
                })
        cols=list(all_rows[0]);resultdir.mkdir(parents=True,exist_ok=True)
        with out.open('w',newline='') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=cols);w.writeheader();w.writerows(all_rows)
        checklist_rows=[x for x in all_rows if x['Audit_flag']=='REVIEW_3PRIME_OR_OVERLAP'
                        or x['Audit_flag']=='ANNOTATION_MISSING_OR_UNREADABLE'
                        or (('NCL_CHIMERIC_CANDIDATE' in x['Original_classification']) and x['Audit_flag']!='REVIEW_3PRIME_OR_OVERLAP')]
        with checklist.open('w',newline='') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=cols);w.writeheader();w.writerows(checklist_rows)
        print('Species with NCL signatures:',len({x['Species'] for x in all_rows}))
        print('NCL loci:',len(all_rows))
        for flag in sorted(set(x['Audit_flag'] for x in all_rows)):
            print(flag,':',sum(x['Audit_flag']==flag for x in all_rows))
        print('Species without usable local BED:',len(set(x['Species'] for x in all_rows if not x['BED_has_valid_records'])))
        print('Malformed signature lines:',len(problems))
        print('ALL LOCI:',out)
        print('MANUAL REVIEW:',checklist)
        print('IMPORTANT: if BED annotations are missing, NO negative result can be inferred.')
        print('IMPORTANT: inspect original signature columns and RNA junctions before changing the manuscript.')
        if problems:
            print('Malformed rows (first 5):',problems[:5])
    
    if __name__=='__main__': main()

