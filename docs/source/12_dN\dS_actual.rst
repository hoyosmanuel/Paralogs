HyPhy 3: Reanálisis de selección evolutiva en quimeras HOST–TE
================================================================

Objetivo
--------

Construir un procedimiento completamente reproducible mediante código
para analizar selección evolutiva en regiones codificantes quiméricas
de murciélagos, sin depender de la plataforma web Datamonkey.

Conservar los análisis anteriores como referencia histórica.


Proyecto
--------

Directorio principal::

    /Volumes/Expansion/project3/bat_HTF_genomic_analysis

Análisis original::

    ANALYSES/dNdS_FINAL/

Nuevo espacio de trabajo::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/


Datos biológicos
----------------

Cinco arquitecturas:

* FBXO3
* NCL
* TMPO
* ZFP28
* ZKSCAN7

Cada arquitectura fue dividida en tres componentes:

* HOST: región codificante derivada del gen hospedador.
* TE: región codificante derivada del elemento transponible.
* CHIMERA: región codificante quimérica completa.

Total: 15 conjuntos de alineamientos codónicos.

KATNAL1 fue excluido del análisis comparativo de la quimera
completa porque solamente una especie pasó el control de calidad.


Análisis anteriores
-------------------

Software verificado: HyPhy 2.5.103 (macOS ARM64).

BUSTED-MH:

* 15 análisis.
* Variación de tasas sinónimas activada.
* Sustituciones multinucleotídicas dobles y triples.
* Todas las ramas como ramas de prueba.
* Ningún análisis significativo a p < 0.05.

FEL-MH:

* 15 análisis.
* Variación de tasas sinónimas activada.
* Sustituciones dobles y triples.
* Tasas multinucleotídicas por sitio estimadas.
* 43 sitios TE con señal purificadora nominal (p <= 0.05).
* 2 sitios TE con señal diversificadora nominal (p <= 0.05).

Corrección Benjamini-Hochberg:

* Aplicada independientemente a cada alineamiento.
* Ningún codón conservó significancia a q <= 0.05.
* Segunda implementación de verificación: script 38.

Estos resultados no demuestran selección purificadora significativa
en codones individuales.

Los archivos JSON y LOG originales deben conservarse intactos.


Nueva etapa: HyPhy 3
--------------------

Comprobar la existencia de los archivos
---------------------------------------------

Se ejecutó::

    find ANALYSES/dNdS_FINAL/BUSTED_MH_inputs \
      -type f -name "*.BUSTED_MH.nex" | sort

Resultado:

Se identificaron 15 archivos NEXUS correspondientes a
cinco arquitecturas y tres componentes por arquitectura:

* FBXO3: HOST, TE, CHIMERA.
* NCL: HOST, TE, CHIMERA.
* TMPO: HOST, TE, CHIMERA.
* ZFP28: HOST, TE, CHIMERA.
* ZKSCAN7: HOST, TE, CHIMERA.


1.2. Inspeccionar dimensiones y árboles
---------------------------------------

Se ejecutó::

    for f in ANALYSES/dNdS_FINAL/BUSTED_MH_inputs/*/*.BUSTED_MH.nex; do
        echo "=== $(basename "$f") ==="
        grep -Ei 'DIMENSIONS|BEGIN TREES|TREE tree1' "$f" | cut -c 1-110
    done

Los 15 archivos contienen declaraciones DIMENSIONS,
bloques BEGIN TREES y definiciones TREE tree1.

Las dimensiones declaradas fueron:

.. list-table:: Dimensiones de los alineamientos
   :header-rows: 1

   * - Arquitectura
     - Especies
     - HOST (codones)
     - TE (codones)
     - CHIMERA (codones)
   * - FBXO3
     - 6
     - 119
     - 432
     - 552
   * - NCL
     - 15
     - 411
     - 556
     - 968
   * - TMPO
     - 15
     - 138
     - 605
     - 744
   * - ZFP28
     - 15
     - 270
     - 121
     - 392
   * - ZKSCAN7
     - 7
     - 278
     - 81
     - 360

Las longitudes en codones se calcularon dividiendo
el valor nchar entre tres.


1.3. Validar alineamientos y correspondencia de taxones
-------------------------------------------------------

Se utilizó Biopython para leer los archivos NEXUS,
extraer sus alineamientos y árboles filogenéticos,
comparar los identificadores de taxones y comprobar
la longitud de los alineamientos.

El comando ejecutado fue::

    python3 - <<'PY'
    from pathlib import Path
    from Bio import AlignIO, Phylo
    from io import StringIO

    root = Path("ANALYSES/dNdS_FINAL/BUSTED_MH_inputs")

    for archdir in sorted(p for p in root.iterdir() if p.is_dir()):
        trees = []
        for f in sorted(archdir.glob("*.BUSTED_MH.nex")):
            alignment = AlignIO.read(str(f), "nexus")
            tree = Phylo.read(str(f), "nexus")
            seq_ids = {r.id for r in alignment}
            tip_ids = {t.name for t in tree.get_terminals()}
            trees.append(Phylo._io.write if False else
                         tuple(sorted(tree.get_terminals(),
                                      key=lambda x: x.name)))
            assert seq_ids == tip_ids, f"Taxon mismatch: {f}"
            assert alignment.get_alignment_length() % 3 == 0, \
                f"Frame length: {f}"
            print(
                f"{f.name}: taxa={len(seq_ids)} "
                f"codons={alignment.get_alignment_length()//3} PASS"
            )
    PY

Resultado:

Los 15 archivos produjeron PASS.

Esto confirma que:

* Biopython puede leer los alineamientos y árboles.
* Los taxones del alineamiento coinciden con los terminales
  del árbol correspondiente.
* Todas las longitudes de alineamiento son divisibles por tres.
* Los números de especies y codones coinciden con las
  dimensiones previamente inspeccionadas.


2. Verificación de las topologías filogenéticas
===============================================

Objetivo
--------

Determinar si los archivos NEXUS correspondientes a HOST,
TE y CHIMERA utilizan exactamente la misma topología
filogenética dentro de cada arquitectura.

Esta comprobación es necesaria para garantizar que las
comparaciones entre componentes se realizan sobre árboles
con la misma estructura evolutiva.


Procedimiento
-------------

Se utilizó Biopython para leer los árboles incluidos en
los archivos NEXUS.

Para comparar las topologías se representó recursivamente
cada nodo mediante el conjunto ordenado de sus descendientes,
independientemente del orden de las ramas en el archivo.

Se ejecutó::

    python3 - <<'PY'
    from pathlib import Path
    from Bio import Phylo

    root = Path("ANALYSES/dNdS_FINAL/BUSTED_MH_inputs")

    for archdir in sorted(p for p in root.iterdir() if p.is_dir()):
        trees = {}

        for component in ("HOST", "TE", "CHIMERA"):
            f = archdir / f"{archdir.name}.{component}.BUSTED_MH.nex"
            tree = Phylo.read(str(f), "nexus")

            def topology(node):
                if node.is_terminal():
                    return node.name
                return tuple(sorted(
                    (topology(child) for child in node.clades),
                    key=str
                ))

            trees[component] = topology(tree.root)

        identical = (
            trees["HOST"] == trees["TE"] == trees["CHIMERA"]
        )

        print(
            f"{archdir.name}: "
            f"{'PASS' if identical else 'FAIL'}"
        )
    PY


Resultados
----------

La comparación produjo::

    FBXO3: PASS
    NCL: PASS
    TMPO: PASS
    ZFP28: PASS
    ZKSCAN7: PASS


Interpretación
--------------

Las cinco arquitecturas presentan topologías enraizadas
idénticas entre sus tres componentes.

Por tanto, se confirmó que HOST, TE y CHIMERA utilizan
la misma estructura filogenética dentro de cada
arquitectura.

Esta verificación no demuestra que las longitudes de las
ramas sean idénticas, ni que los árboles hayan sido
correctamente derivados del árbol maestro.


Estado
------

* Arquitecturas examinadas: 5.
* Comparaciones satisfactorias: 5.
* Diferencias topológicas detectadas: 0.
* Archivos originales modificados: 0.

Verificación de topologías: PASS.

Siguiente comprobación: identidad de las longitudes
de rama entre HOST, TE y CHIMERA.



4. Selección de arquitecturas para los análisis de dN/dS
=========================================================

Objetivo
--------

Documentar los criterios utilizados para seleccionar las
arquitecturas quiméricas que podían someterse a análisis
comparativos de selección evolutiva.

La identificación de una arquitectura quimérica no implica
automáticamente que existan suficientes secuencias codificantes
homólogas para realizar análisis de dN/dS.


Conjunto inicial
----------------

El manuscrito describe nueve arquitecturas quiméricas
HOST–TE de alta confianza.

Para la preparación de los análisis de selección se
consideraron inicialmente seis arquitecturas:

* FBXO3
* NCL
* TMPO
* ZFP28
* ZKSCAN7
* KATNAL1

Las arquitecturas FAU, ZNF24 y AP1AR quedaron fuera de
esta primera ronda debido a que no reunían al menos
tres recuperaciones codificantes completas y comparables,
según la documentación previa de reconstrucción de CDS.


Evaluación de KATNAL1
---------------------

La arquitectura KATNAL1–hAT1_mMyo fue recuperada
genómicamente en siete especies de Myotis.

Sin embargo, la reconstrucción literal de sus regiones
codificantes mostró que solamente Myotis daubentonii
conservaba una ORF quimérica intacta equivalente
a la referencia.

Las otras seis especies presentaban alteraciones del
marco de lectura y/o codones de parada internos.

Por tanto, KATNAL1 fue excluida del análisis comparativo
de dN/dS de la quimera completa.

Esta exclusión no implica ausencia de la arquitectura
genómica ni de evidencia transcriptómica.

Los resultados describen las secuencias de los
ensamblajes disponibles y no permiten descartar
errores locales de ensamblaje.


Conjunto final analizado
------------------------

Se conservaron cinco arquitecturas:

.. list-table:: Arquitecturas incluidas
   :header-rows: 1

   * - Arquitectura
     - Especies
     - Componentes
   * - FBXO3
     - 6
     - HOST, TE, CHIMERA
   * - NCL
     - 15
     - HOST, TE, CHIMERA
   * - TMPO
     - 15
     - HOST, TE, CHIMERA
   * - ZFP28
     - 15
     - HOST, TE, CHIMERA
   * - ZKSCAN7
     - 7
     - HOST, TE, CHIMERA

Cada arquitectura fue evaluada mediante tres
alineamientos codónicos independientes.

Total:

* 5 arquitecturas.
* 15 alineamientos codónicos.
* 15 archivos NEXUS originales.


Consideraciones metodológicas
-----------------------------

El número de especies disponibles no es equivalente
a la potencia estadística del análisis.

La ausencia de significancia estadística tampoco
demuestra evolución neutral.

La exclusión de una arquitectura del análisis de
dN/dS no constituye evidencia de ausencia de
función biológica o de exaptación.

La selección de las cinco arquitecturas corresponde
al conjunto de datos utilizado en los análisis
históricos de BUSTED-MH y FEL-MH.


Estado
------

* Arquitecturas identificadas en el manuscrito: 9.
* Arquitecturas consideradas inicialmente: 6.
* Arquitecturas incluidas finalmente: 5.
* Alineamientos codónicos disponibles: 15.

La justificación de KATNAL1 se encuentra documentada
en el análisis previo de reconstrucción de CDS.

La procedencia y el control de calidad de las demás
arquitecturas deberán vincularse con sus respectivos
archivos originales para completar la trazabilidad.


5. Verificación del árbol filogenético maestro
==============================================

Objetivo
--------

Comprobar la existencia y la estructura básica del árbol
filogenético maestro utilizado como referencia para los
análisis de selección evolutiva.

Archivo de referencia::

    ANALYSES/dNdS_FINAL/config/master_bat_tree.tre


Procedimiento
-------------

Se utilizó Biopython para leer el árbol en formato Newick
y comprobar su número de terminales, identificadores
únicos y disponibilidad de longitudes de rama.

Se ejecutó::

    python3 - <<'PY'
    from pathlib import Path
    from Bio import Phylo

    path = Path("ANALYSES/dNdS_FINAL/config/master_bat_tree.tre")

    print("Archivo existente:", path.is_file())

    if path.is_file():
        tree = Phylo.read(path, "newick")
        tips = tree.get_terminals()
        names = [tip.name for tip in tips]

        print("Número de terminales:", len(tips))
        print("Nombres únicos:", len(set(names)))
        print("Longitudes de rama ausentes:",
              sum(t.branch_length is None for t in tree.find_clades()))
        print("Primeros 10 taxones:", sorted(names)[:10])
    PY


Resultados
----------

La ejecución produjo::

    Archivo existente: True
    Número de terminales: 111
    Nombres únicos: 111
    Longitudes de rama ausentes: 0


Interpretación
--------------

El árbol maestro contiene 111 terminales con identificadores
únicos y longitudes de rama definidas para todos sus nodos.

La prueba confirma su disponibilidad y lectura correcta,
pero no demuestra todavía que los árboles incluidos en los
15 archivos NEXUS sean subárboles derivados de él.


Estado
------

* Existencia del archivo: PASS.
* Lectura mediante Biopython: PASS.
* Identificadores únicos: PASS.
* Longitudes de rama presentes: PASS.

Pendiente: comparar cada árbol NEXUS con el árbol maestro
podado al conjunto correspondiente de especies.



6. ACTUALIZACIÓN DEL 22–23 DE NOVIEMBRE DE 2026: NCL–CHARLIE8
================================================================

IMPORTANTE: esta sección es el estado de trabajo más reciente y tiene
precedencia sobre las cifras históricas de NCL y KATNAL1 anteriores.
Los análisis originales BUSTED-MH/FEL-MH y sus NEXUS siguen siendo
históricos; NO se han recalculado con estas CDS nuevas.

6.1. Decisión sobre especies para HyPhy 3
----------------------------------------

Las siguientes listas fueron propuestas por el investigador como conjuntos
CANDIDATOS, no como secuencias codificantes validadas:

* NCL–CHARLIE8: las 28 especies de la auditoría independiente de NCL,
  EXCLUYENDO ``Pipistrellus_kuhlii`` y ``Cistugo_seabrae``: 26 especies.
* TMPO (14): Myotis_myotis, Myotis_daubentonii, Myotis_mystacinus,
  Myotis_evotis, Myotis_thysanodes, Myotis_auriculus,
  Myotis_californicus, Myotis_occultus, Myotis_lucifugus,
  Myotis_velifer, Myotis_yumanensis, Myotis_nigricans,
  Myotis_vivesi, Myotis_volans.
* KATNAL1 (5): Myotis_myotis, Myotis_daubentonii,
  Myotis_mystacinus, Myotis_evotis, Myotis_volans.
* ZKSCAN7 (7): Rhinolophus_perniger, Rhinolophus_trifoliatus,
  Rhinolophus_foetidus, Rhinolophus_pearsonii, Rhinolophus_affinis,
  Rhinolophus_hipposideros, Rhinolophus_ferrumequinum.
* FBXO3 (3): Rhinolophus_sinicus, Rhinolophus_hipposideros,
  Rhinolophus_ferrumequinum.
* ZFP28–GOLEM-1p (22): Rousettus_aegyptiacus,
  Hypsignathus_monstrosus, Eonycteris_spelaea, Cynopterus_sphinx,
  Rhinopoma_muscatellum, Rhinopoma_microphyllum,
  Rhinolophus_ferrumequinum, Rhinolophus_hipposideros,
  Rhinolophus_affinis, Rhinolophus_pearsonii,
  Rhinolophus_trifoliatus, Rhinolophus_perniger,
  Hipposideros_caffer, Hipposideros_abae, Hipposideros_jonesi,
  Hipposideros_armiger, Hipposideros_swinhoii, Hipposideros_larvatus,
  Molossus_molossus, Molossus_alvarezi, Eumops_nanus,
  Tadarida_brasiliensis.

Total: 6 arquitecturas y 77 entradas especie–arquitectura. Los nombres
científicos deben cotejarse con los identificadores de los FASTA; la
ortografía aquí está normalizada provisionalmente, no validada contra
todos los ensamblajes. KATNAL1 era histórica y provisionalmente excluida;
su reingreso está propuesto, NO validado. Tres secuencias para FBXO3
suponen potencia estadística limitada. Todas requieren revisión de
ortología, límites codificantes y ausencia de paradas internas.

6.2. Listado EXACTO de 26 candidatos NCL
---------------------------------------

La consulta a
``ANALYSES/NCL/results/NCL_28_species_exon_TE_recovery.tsv``
produjo 26 especies tras excluir las dos citadas::

    Antrozous_pallidus
    Corynorhinus_mexicanus
    Corynorhinus_townsendii
    Eptesicus_fuscus
    Eptesicus_nilssonii
    Lasiurus_ega
    Myotis_auriculus
    Myotis_californicus
    Myotis_daubentonii
    Myotis_evotis
    Myotis_lucifugus
    Myotis_myotis
    Myotis_mystacinus
    Myotis_nigricans
    Myotis_occultus
    Myotis_pilosus
    Myotis_thysanodes
    Myotis_velifer
    Myotis_vivesi
    Myotis_volans
    Myotis_yumanensis
    Nyctalus_aviator
    Pipistrellus_nathusii
    Pipistrellus_pygmaeus
    Plecotus_auritus
    Vespertilio_murinus

En la reconstrucción histórica existían 23 registros (22 de las 26
seleccionadas + Pipistrellus_kuhlii); de los 22 seleccionados,
15 eran ``PASS`` y 7 fallaban QC. Faltaban cuatro CDS:
``Nyctalus_aviator``, ``Pipistrellus_nathusii``,
``Pipistrellus_pygmaeus``, ``Vespertilio_murinus``.

Las siete que FALLARON el QC anterior eran:

.. list-table:: Fallos del QC HISTÓRICO de NCL
   :header-rows: 1

   * - Especie
     - nt quimera antes
     - stops internos antes
     - estado
   * - Corynorhinus_mexicanus
     - 2855
     - 50
     - fallo
   * - Corynorhinus_townsendii
     - 2855
     - 49
     - fallo
   * - Eptesicus_fuscus
     - 2852
     - 51
     - fallo
   * - Eptesicus_nilssonii
     - 2855
     - 52
     - fallo
   * - Lasiurus_ega
     - 2852
     - 55
     - fallo
   * - Myotis_volans
     - 2885
     - 60
     - fallo
   * - Plecotus_auritus
     - 2855
     - 50
     - fallo

El indicador histórico ``FRAMESHIFT`` fue asignado por longitud
``%3 != 0`` de la reconstrucción; no constituye por sí solo prueba de
mutación frameshift en el genoma, ni los stops de una reconstrucción
incorrecta prueban pseudogenización.

6.3. Hallazgo fundamental: RNA-seq de Plecotus auritus
-----------------------------------------------------

En ``/Volumes/Expansion/IGV/Plecotus_auritus/`` existen los BAM,
BAI y GTF de cinco experimentos: ``ERR13757870`` a ``ERR13757874``.
Se comprobó el locus ``OY734026`` en hebra ``-``:

* intervalo NCL histórico: 32902868–32914225;
* exón 3 reconstruido: 32910222–32910549 (328 nt);
* exón 4 reconstruido previamente: 32909700–32909836 (137 nt);
* exón 5: 32909073–32909159 (87 nt);
* coincidencia CHARLIE8 histórica: 32899963–32901868.

Se contaron uniones de RNA-seq mediante la operación ``N`` de CIGAR,
con coordenadas (última base del bloque alineado a la izquierda,
primera base del bloque alineado a la derecha):

.. list-table:: Soporte de uniones en Plecotus auritus
   :header-rows: 1

   * - BAM
     - exón 4–5 (32909159,32909666)
     - NCL exón 9–CHARLIE8 (32901868,32902868)
   * - ERR13757870
     - 14
     - 5
   * - ERR13757871
     - 7
     - 1
   * - ERR13757872
     - 31
     - 14
   * - ERR13757873
     - 27
     - 5
   * - ERR13757874
     - 11
     - 4
   * - TOTAL
     - 90
     - 29

La unión previamente supuesta ``(32909159,32909700)`` produjo 0
alineamientos en las cinco muestras. En contraste, el inicio del
bloque alineado del exón 4 en coordenadas ascendentes es 32909666,
y la otra unión observada es ``(32909836,32910222)``. Esto respalda
para ese modelo un exón 4 de ``32909666–32909836`` (171 nt),
34 nt más largo que el reconstruido inicialmente. Hay soporte de
las uniones individuales, no secuenciación de una única molécula
completa de extremo a extremo. Estos conteos son alineamientos,
no necesariamente moléculas independientes. Cuatro BAM emitieron
advertencia de índice BAI anterior al BAM (fecha de modificación):
verificar índices si se necesita auditoría adicional.

Secuencia de 34 nt incorporada SOLO para la reconstrucción provisional::

    GTGGAGGTAGAGGTGGAGGTGGAGGAAGAGGAAG

QC antes: 2855 nt, módulo 3 = 2, 50 stops, sin stop terminal.
QC provisional después: 2889 nt, módulo 3 = 0, 0 stops internos,
ATG inicial y TAA terminal (stop en codón 963).

La evidencia de RNA apoya la arquitectura de transcrito NCL–CHARLIE8
y el límite corregido; NO demuestra proteína funcional ni que todas
las uniones se encuentren en una única isoforma completa.

6.4. Corrección por homología de otras cinco especies
---------------------------------------------------

En el exón 4 de cinco especies, el script histórico había escogido
motivos GT anteriores dentro de una región repetitiva. Se halló un
donador alternativo con contexto homólogo al de Plecotus::

    AGGAAG|GTAATTGA

NO hay RNA-seq confirmado para estas cinco especies. Son extensiones
genómicas provisionales basadas en homología y continuidad de ORF,
no validaciones transcriptómicas directas.

.. list-table:: Extensiones exón 4 verificadas EN MEMORIA
   :header-rows: 1

   * - Especie
     - Anterior
     - Extensión
     - Nueva longitud
     - CDS nueva nt
     - stops internos nueva
   * - Corynorhinus_mexicanus
     - 137
     - +34
     - 171
     - 2889
     - 0
   * - Corynorhinus_townsendii
     - 137
     - +34
     - 171
     - 2889
     - 0
   * - Eptesicus_fuscus
     - 134
     - +37
     - 171
     - 2889
     - 0
   * - Eptesicus_nilssonii
     - 137
     - +40
     - 177
     - 2895
     - 0
   * - Lasiurus_ega
     - 137
     - +40
     - 177
     - 2892
     - 0

Las cinco CDS provisionales tienen ``ATG`` inicial, ``TAA`` terminal,
longitud múltiplo de 3 y ninguna parada interna. Secuencias añadidas::

    Corynorhinus_mexicanus:   GTGGAGGTGGAGGTGGAGGTGGAGGAAGAGGAAG
    Corynorhinus_townsendii: GTGGAGGTGGAGGTGGAGGTGGAGGAAGAGGAAG
    Eptesicus_fuscus:        GTGGTGGAGGTGGAGGTGGAGGTGGAGGAAGAGGAAG
    Eptesicus_nilssonii:    GTGGAGGTGGAGGTGGAGGTGGAGGTGGAGGAAGAGGAAG
    Lasiurus_ega:           GTGGAGGTAGAAGTGGAGGTGGAGGTGGAGGAAGAGGAAG

IMPORTANTE: en ``Myotis_volans`` el exón 4 reconstruido (135 nt)
ya llega a ``ATGAAG|GTAATTGA``. Su problema era DISTINTO.

6.5. Corrección del inicio CDS de Myotis volans
----------------------------------------------

El exón 2 recuperado tenía 98 nt frente a 114 en la referencia
``Myotis_myotis``. Su homólogo ``ATGGCT`` comienza en la posición 15,
no en la posición 31. El script anterior recortaba FIJAMENTE los
primeros 30 nt de cada exón 2, dejando solamente 68 nt de la región
codificante, en lugar de los 84 nt homólogos desde ATG.

Al recortar 14 nt, en lugar de 30 nt, se recuperan 16 nt:

* reconstrucción anterior: 2885 nt, primer codón AGG,
  60 stops internos, codón terminal observado GTT;
* corregida: 2901 nt, primer codón ATG, cero stops internos,
  TAA terminal, stop en codón 967.

Esta corrección se basa en posición homóloga del ATG; no implica
evidencia transcriptómica directa de la isoforma de ``M. volans``.

6.6. ARCHIVOS NUEVOS YA ESCRITOS EN EL MAC — NO SOBREESCRIBIR
------------------------------------------------------------

La siguiente creación de archivos fue confirmada mediante salida
terminal del usuario (Mac; NO disponibles necesariamente en el entorno
de ChatGPT):

Directorio::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/

Archivos::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_22_chimera_cds_provisional.fa
    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_22_reconstruction_QC.tsv

La escritura terminó con::

    CDS saved: 22
    Original historical files were not modified.

Se conservaron 15 secuencias históricas ``PASS`` sin modificaciones;
se añadieron las siete correcciones indicadas anteriormente. Se excluyó
``Pipistrellus_kuhlii``. Los registros FASTA mantienen IDs como::

    NCL|Plecotus_auritus|OY734026|strand=-|CHIMERA

Se verificó posteriormente que en **las 22 secuencias**:

* el segmento TE terminal permanece idéntico a
  ``ANALYSES/dNdS_FINAL/sequences/NCL/te_cds.fa``;
* ``HOST_nt % 3 == 1``;
* el codón que atraviesa la unión HOST–TE es ``GGA``;
* las reconstrucciones contienen ATG inicial, 0 stops internos,
  TAA terminal y longitud divisible por 3.

Las regiones HOST tienen tamaños diversos (ej.: Myotis_myotis 1213 nt,
Plecotus_auritus 1216 nt, Myotis_volans 1228 nt); TE = 1673 nt
para estas 22 reconstrucciones. Mantener la fase es un control de
consistencia, NO prueba independiente de verdadera unión ni ortología.

PRIORIDAD: NO volver a ejecutar ``17_reconstruct_NCL_chimeras.py``
sobre rutas históricas, porque sobrescribiría las secuencias/QC originales.
NO sobrescribir los dos archivos provisionales nuevos; crear versiones
separadas si se hacen otras correcciones.

6.7. Cuatro CDS pendientes y evidencia de sus loci
--------------------------------------------------

``ANALYSES/NCL/results/NCL_9exons_Charlie8_positive.tsv`` devuelve:

.. list-table:: Cuatro candidatos sin CDS construida
   :header-rows: 1

   * - Especie
     - Scaffold BLAST
     - Hebra
     - Distancia exón9–TE
     - Identidad TE
     - Cobertura TE
   * - Nyctalus_aviator
     - CM072928
     - -
     - 994 bp
     - 95.929 %
     - 99 %
   * - Pipistrellus_nathusii
     - hap2_chr_8
     - -
     - 922 bp
     - 96.138 %
     - 99 %
   * - Pipistrellus_pygmaeus
     - OX465310
     - -
     - 999 bp
     - 95.720 %
     - 99 %
   * - Vespertilio_murinus
     - OZ004711
     - -
     - 999 bp
     - 95.151 %
     - 100 %

BLAST y genomas comprobados en el Mac:

* Nyctalus_aviator:
  ``SPECIES/Nyctalus_aviator/mNycAvi1.1.pri.fa``;
  TE HSP ``CM072928:27729146–27731050``.
* Pipistrellus_nathusii:
  únicamente se encontró ``SPECIES/Pipistrellus_nathusii/mPipNat1.1.hap1.fa``.
  El resultado BLAST se refiere a ``hap2_chr_8``; TE HSP en
  ``hap2_chr_8:27717136–27719040``. **EL USUARIO CONFIRMÓ QUE NO TIENE
  EL GENOMA HAP2** y una búsqueda de ``*PipNat*hap2*.fa`` bajo
  ``/Volumes/Expansion/project3`` no encontró nada. NO aplicar las
  coordenadas hap2 al FASTA hap1. Para incluir esta especie habría
  que obtener el ensamblaje hap2 exacto o rehacer la identificación
  de loci/exones/TE independientemente sobre hap1, con coordenadas
  nuevas y sin asumir equivalencia.
* Pipistrellus_pygmaeus:
  ``SPECIES/Pipistrellus_pygmaeus/mPipPyg1.1.pri.fa``;
  TE HSP ``OX465310:27958546–27960450``.
* Vespertilio_murinus:
  ``SPECIES/Vespertilio_murinus/mVesMur1.1.pri.fa``;
  TE HSP ``OZ004711:64521191–64523093``. Existe evidencia previa
  de 27 alineamientos RNA-seq para la unión NCL–CHARLIE8 en cinco
  BAM (4,8,2,6,7) y auditoría que corrigió un falso negativo de TOGA.

Las cuatro especies tienen 9 hits NCL en orden y TE próximo en 3′.
ESTO NO GARANTIZA nueve exones completos, CDS intacta o identidad
inequívoca de parálogo. En ``P. pygmaeus`` la primera auditoría
``NCL_all_loci_exon_CHARLIE8.tsv`` clasificó un locus de nueve
exones en OX465310 como sin CHARLIE8, pero la auditoría posterior
recuperó el TE 3′; no confundir clasificaciones de los dos métodos.

6.8. Estado científico y lista de pendientes
--------------------------------------------

HECHOS CONSTATADOS

* 26 especies NCL propuestas; 22 secuencias quiméricas provisionales
  guardadas y comprobadas sin stops internos; cuatro pendientes.
* Para Plecotus auritus hay datos de RNA-seq de la unión exón 4–5
  y la unión NCL–CHARLIE8; la corrección de 34 nt es informada por RNA.
* Para cinco especies adicionales las extensiones del exón 4 son
  inferencias de homología. Para Myotis volans el ATG homólogo
  corrige el recorte artificial de exón 2.
* El genoma hap2 de Pipistrellus nathusii NO está disponible
  según el usuario. Solo hay hap1 visible en SPECIES.

LO QUE **NO SE SABE TODAVÍA**

* Si todas las 22 reconstrucciones representan isoformas reales y
  parálogos ortólogos adecuados para pruebas codónicas.
* Si se obtendrán ORF completas para las cuatro especies pendientes.
* Si todas las uniones exónicas inferidas están respaldadas por RNA.
* Si el TE está completo para todas las especies en términos de CDS,
  pese a los resultados de BLAST; ni funcionalidad proteica.

PENDIENTES, EN ESTE ORDEN RECOMENDADO

#. Confirmar en disco los dos archivos provisionales y conservarlos.
#. Reconstruir Nyctalus aviator, Pipistrellus pygmaeus y
   Vespertilio murinus desde las coordenadas de sus FASTA, con la
   debida selección del locus y auditoría de exones/TE; comenzar por
   Vespertilio murinus debido a evidencia RNA previa.
#. Resolver Pipistrellus nathusii: localizar/obtener hap2 o rehacer
   el análisis de secuencias con hap1; no transferir coordenadas.
#. Crear resultados HOST, TE y CHIMERA para todas las CDS aceptadas,
   respetar codón de unión fase 1 (GGA), alinear por codones,
   comprobar ortología y filtros de QC independientes.
#. Comparar los árboles específicos con el maestro de 111 terminales
   debidamente podado; evaluar topología y longitudes de rama.
#. Auditar las otras cinco arquitecturas propuestas (incluida
   KATNAL1); sus números de especies son CANDIDATOS.
#. Solo entonces ejecutar nuevas pruebas HyPhy, conservando intactos
   los 15 NEXUS/resultados históricos y registrando versión/opciones.

NO ESTÁ EJECUTADO HyPhy 3 sobre las 22 CDS nuevas.

6.9. Punto exacto para reanudar mañana
-------------------------------------

El usuario acabó diciendo que **NO posee el genoma hap2 de
Pipistrellus nathusii**. No repetir esa búsqueda y no tratar de
extraer ``hap2_chr_8`` de ``mPipNat1.1.hap1.fa``.

Preferencia de trabajo: proporcionar UN SOLO COMANDO cada vez, en
macOS y de solo lectura cuando se esté inspeccionando; solicitar
salida terminal y continuar. Para acciones que creen archivos,
explicar previamente exactamente qué se escribirá y jamás
sobrescribir históricos. Distinguir hechos, inferencias y dudas.

**Siguiente análisis recomendado: Vespertilio murinus.**
Antes de generar una CDS, inspeccionar con un solo comando las filas
BLAST correspondientes a los exones 2–9 en
``ANALYSES/NCL/intermediate/NCL_103_blast/Vespertilio_murinus.tsv``
y la anotación previa del locus, y escoger cada exón por mismo scaffold,
hebra, orden, proximidad, cobertura y, si procede, RNA.

COMANDO INICIAL OPCIONAL (SOLO LECTURA), DESDE RAÍZ::

    awk -F '\t' '$1 ~ /^NCL_exon[2-9]_match$/ && $2=="OZ004711" \
      {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$12}' \
      ANALYSES/NCL/intermediate/NCL_103_blast/Vespertilio_murinus.tsv \
      | head -90

**No ejecutar ningún script histórico de reconstrucción con salidas
existentes sin redireccionarlo a un directorio nuevo.**
