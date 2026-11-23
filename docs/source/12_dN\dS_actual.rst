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


7. ACTUALIZACION DEL 23 DE NOVIEMBRE DE 2026: 25 CDS, ALINEAMIENTOS Y ARBOLES
=============================================================================

ESTA SECCION ES LA CONTINUACION MAS RECIENTE DE LA BITACORA.
Tiene precedencia sobre el estado de avance y los pendientes de la seccion 6,
particularmente 6.7-6.9. NO reemplaza la documentacion historica ni implica
que las inferencias de ortologia, expresion y seleccion esten terminadas.

7.1. Estado biologico exacto: 25 de 26 candidatos NCL
----------------------------------------------------

De los 26 candidatos definidos en 6.2 (28 recuperaciones iniciales menos
``Cistugo_seabrae`` y ``Pipistrellus_kuhlii``), se reconstruyeron 25 CDS
quimericas provisionales NCL-CHARLIE8, que pasaron los controles de ORF:
ATG inicial, TAA terminal y ausencia de codones de parada internos.
Estas comprobaciones NO prueban por si mismas expresion, ortologia inequivoca,
funcionalidad proteica ni homologia de todas las posiciones del alineamiento.

La vigesimosexta especie candidata, ``Pipistrellus_nathusii``, SIGUE PENDIENTE:
el BLAST corresponde a ``hap2_chr_8`` y en el Mac solo se dispone de un
ensamblaje hap1. NO transferir coordenadas entre haplotipos.

``Cistugo_seabrae`` se mantiene excluida de este conjunto de dN/dS por
recuperacion insuficientemente fiable del locus/exon 1; ``Pipistrellus_kuhlii``
por ambiguedad de asignacion entre copias NCL. Estas exclusiones son decisiones
de QC del conjunto de seleccion, no una negacion de hits genomicos.

Archivo CDS ACTUAL, provisional y NO SOBRESCRIBIR::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_chimera_cds_Antrozous_corrected.fa

Archivos historicos intermedios, tambien NO SOBRESCRIBIR::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_22_chimera_cds_provisional.fa
    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_22_reconstruction_QC.tsv
    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_23_chimera_cds_provisional.fa
    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_24_chimera_cds_provisional.fa
    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_chimera_cds_provisional.fa

**Importante:** esta lista de rutas refleja las escrituras informadas en el
terminal durante el trabajo. Los archivos residen en el Mac del investigador;
no se han vuelto a inspeccionar directamente desde este entorno de ChatGPT.

7.2. Nuevas reconstrucciones posteriores al conjunto de 22
---------------------------------------------------------

``Vespertilio_murinus`` (especie 23):

* Ensamblaje ``SPECIES/Vespertilio_murinus/mVesMur1.1.pri.fa``;
  scaffold ``OZ004711``, hebra negativa.
* RNA-seq en cinco BAM ``ERR13757845`` a ``ERR13757849``: 27
  alineamientos de junction NCL-CHARLIE8 (4, 8, 2, 6, 7 por BAM).
* Exones codificantes: exon 2 64534976-64535083 (78 nt desde el ATG);
  exon 3 64532829-64533132 (304 nt); exon 4 64532275-64532445
  (171 nt); exon 5 64531425-64531511 (87 nt);
  exon 6 64530361-64530499 (139 nt);
  exon 7 64528873-64528997 (125 nt);
  exon 8 64528112-64528235 (124 nt);
  exon 9 64524093-64524250 (158 nt).
* TE codificante reconstruido: ``OZ004711:64521421-64523093``
  (1673 nt); HOST = 1186 nt; quimera = 2859 nt, 952 aa sin stop.
* Se obtuvo alineamiento TE de 1673/1673 nt, sin gaps, con identidad
  de 96.294 % respecto a la region de referencia Myotis myotis usada
  para la reconstruccion. No confundir este porcentaje con la identidad
  del BLAST historico (95.151 %) recogida en la tabla de recuperacion.
* El apoyo RNA indica uniones de empalme observadas, NO prueba
  por si solo una molecula completa ni una proteina funcional.

``Pipistrellus_pygmaeus`` (especie 24):

* Ensamblaje ``SPECIES/Pipistrellus_pygmaeus/mPipPyg1.1.pri.fa``;
  scaffold ``OX465310``, hebra negativa. Esta reconstruccion se baso
  en evidencia genomica y homologia; NO se verifico RNA local aqui.
* Exon 2 27977258-27977371 (84 nt codificantes desde ATG);
  exon 3 27974562-27974865 (304 nt);
  exon 4 27974004-27974180 (177 nt);
  exon 5 27971487-27971573 (87 nt);
  exon 6 27970862-27971000 (139 nt);
  exon 7 27969347-27969471 (125 nt);
  exon 8 27968571-27968694 (124 nt);
  exon 9 27961450-27961607 (158 nt).
* TE 27958778-27960450 (1673 nt); HOST = 1198 nt;
  quimera = 2871 nt; 956 aa sin stop.
* Identidad del segmento TE reconstruido frente al de referencia:
  96.593 %, 1673/1673 nt sin gaps. BLAST historico: 95.720 %.
* El valor ``RNA support = YES`` de una tabla previa requiere
  trazabilidad de la evidencia original; no se ha ratificado en
  esta reconstruccion.

``Nyctalus_aviator`` (especie 25):

* Ensamblaje ``SPECIES/Nyctalus_aviator/mNycAvi1.1.pri.fa``;
  scaffold ``CM072928``, hebra negativa. Sin RNA-seq local verificado.
* Exon 2 codificante 27745830-27745913 (84 nt);
  exon 3 27743197-27743500 (304 nt);
  exon 4 27742636-27742812 (177 nt);
  exon 5 27742067-27742153 (87 nt);
  exon 6 27741470-27741608 (139 nt);
  exon 7 27739968-27740092 (125 nt);
  exon 8 27739191-27739314 (124 nt);
  exon 9 27732036-27732190 (155 nt).
* TE 27729378-27731050 (1673 nt); HOST = 1195 nt;
  quimera = 2868 nt; 955 aa sin stop.
* El exon 9 es 3 nt mas corto que el de P. pygmaeus (155 frente
  a 158 nt). Se prefirio el modelo por homologia/ORF frente a
  un candidato de 174 nt peor alineado.
* Identidad TE reconstruido frente a referencia: 96.772 %,
  1673/1673 nt sin gaps. BLAST historico: 95.929 %.
* Se comprobo posteriormente un gap de un codon completo en
  la columna MAFFT 425, dos codones antes de la union. El codon
  ``TGG`` (W) presente en M. myotis y P. pygmaeus no aparece
  en N. aviator en ese alineamiento. Compatible con la diferencia
  de 3 nt del exon 9; NO se ha demostrado el evento ancestral.

7.3. Correccion especifica de Antrozous pallidus
------------------------------------------------

Para ``Antrozous_pallidus`` se reviso el exon 4 en el scaffold
``manual_scaffold_1``, hebra negativa. La reconstruccion historica
incluia 126 nt (32434710-32434835); se propuso incluir otros
42 nt en 32434668-32434709, obteniendo exon 4 de 168 nt.

* Contexto donador corregido ``AGAGGAAG|GTAATTGA`` frente al
  antiguo ``TGAAGATG|GTGTGGAG``.
* Modelo previo quimera: 2835 nt; corregido: 2877 nt.
* HOST corregido: 1204 nt; sin stops internos, ATG inicial,
  TAA terminal, union hibrida GGA.
* Apoyo: homologia genomica y ORF; NO RNA-seq local comprobado.
* El exon 3 de Antrozous (322 nt) se parece en longitud y homologia
  al de Corynorhinus/Plecotus (328 nt), aunque difiere del de
  Myotis myotis (376 nt); no forzar equivalencia de longitud.

El archivo de 25 CDS anteriormente guardado
``NCL_25_chimera_cds_provisional.fa`` NO es el canon actualizado:
la version activa incorpora esta correccion y se llama
``NCL_25_chimera_cds_Antrozous_corrected.fa``.

7.4. Alineamiento proteico de las 25 quimeras
--------------------------------------------

Programas verificados en el entorno ``samtools_env`` del Mac:

* ``mafft`` instalado en
  ``/Users/manuelhoyos/miniconda3/envs/samtools_env/bin/mafft``.
* ``hyphy`` instalado en
  ``/Users/manuelhoyos/miniconda3/envs/samtools_env/bin/hyphy``.
* ``macse``, ``pal2nal.pl``, ``prank`` y ``muscle`` no se
  encontraron en PATH durante la comprobacion.

Archivo proteico guardado, extraido de las CDS y sin stop terminal::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_chimera_proteins_provisional.fa

Caracteristicas: 25 proteinas, 952-966 aa; ninguna parada interna.

Se ejecuto MAFFT 7.526 con estrategia L-INS-i::

    mafft --localpair --maxiterate 1000 --thread 2 INPUT > OUTPUT

Resultado: codigo de salida 0; convergencia de refinamiento.
Archivo guardado::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_chimera_proteins_MAFFT_LINSI.fa

QC comprobado mediante Biopython:

* 25 secuencias, 25 identificadores unicos;
  cada proteina original se conserva exactamente tras retirar gaps.
* Longitud del alineamiento: 983 columnas de aminoacidos.
* 921 columnas sin gaps; 62 columnas con al menos un gap.
* 35 columnas tienen gaps en al menos 10 especies;
  maximo 24 gaps en una sola columna.
* Se identificaron 12 bloques consecutivos con gaps:
  columnas 2, 9, 26, 40, 61-78, 110, 119-126,
  132-133, 141, 149-150, 189-213 y 425.
* Columna 40: gap en 24/25 especies.
* Columnas 61-78: gap en 10/25 especies cada una.
* Columnas 149-150: gaps en 15/25 especies cada una.
* Bloque 189-213: hasta 20 gaps en una columna y
  21 especies afectadas en alguna posicion del bloque.
* Las diferencias de HOST alrededor de las regiones repetitivas
  ricas en E/D (glutamato/aspartato), especialmente cerca de
  189-213, tienen homologia local potencialmente ambigua.

NO eliminar bloques automaticamente: evaluar sensibilidad de los
analisis con/sin columnas ambiguas, despues de justificar un criterio
reproducible. Una columna sin gaps no demuestra por si sola homologia.

7.5. Union HOST-CHARLIE8 y codon hibrido
---------------------------------------

En las 25 CDS actuales:

* La porcion genomica TE reconstruida mide 1673 nt.
* ``HOST_nt % 3 == 1``.
* El codon que atraviesa la union tiene un nucleotido HOST
  y dos nucleotidos TE y es ``GGA`` en las 25 especies.
* Ese codon ocupa en todas ellas la columna **427** del
  alineamiento MAFFT, sin gaps.
* Se examinaron las columnas 417-437: solo ``Nyctalus_aviator``
  presenta un gap cercano, en la columna 425.
* Para referencia, codones 424/425/426/427:

  - Myotis_myotis: AGT / TGG / AGA / GGA.
  - Pipistrellus_pygmaeus: CGT / TGG / AGA / GGA.
  - Nyctalus_aviator: CGG / --- / AGA / GGA.

Estos controles respaldan la coherencia del alineamiento en la union;
no prueban por si solos que todas las uniones se transcriban ni
que todos los loci sean ortologos.

7.6. Proyeccion a codones y separacion HOST / TE / CHIMERA
--------------------------------------------------------

Se proyecto el alineamiento proteico de MAFFT sobre las CDS originales
mediante Biopython: cada aminoacido se sustituyo por su codon original
y cada gap proteico por ``---``. Se excluyo el ``TAA`` terminal.
La traduccion de cada codon se cotejo con su aminoacido, y la
secuencia sin gaps se confronto con la CDS sin stop.

Archivo CHIMERA guardado::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_chimera_codons_MAFFT_LINSI.fa

QC de CHIMERA:

* 25 especies; 983 columnas codonicas, 2949 nt alineados.
* Codones invalidos: 0; codones stop internos: 0.
* Columnas codonicas con gaps: 62.
* Columna 427: ``GGA`` en todas las especies.
* Los gaps son tripletes completos; no se alteraron las CDS fuente.

Separacion guardada en dos FASTA adicionales::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_HOST_codons_MAFFT_LINSI.fa
    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/sequences/NCL/NCL_25_TE_codons_MAFFT_LINSI.fa

Criterio codonico explicito:

.. list-table:: Componentes NCL HyPhy 3 provisionales
   :header-rows: 1

   * - Componente
     - Columnas MAFFT (1-based)
     - Codones
     - Especies
   * - HOST
     - 1-426
     - 426
     - 25
   * - Union mixta (solo CHIMERA)
     - 427
     - 1
     - 25
   * - TE
     - 428-983
     - 556
     - 25
   * - CHIMERA
     - 1-983
     - 983
     - 25

Los codones completos de HOST y TE se analizaran separadamente,
excluyendo el codon GGA de ambos para no asignarlo enteramente
a un solo origen. La CHIMERA conserva GGA. La separacion y
validacion de tripletes se confirmaron en la salida terminal.
NO se ha ejecutado HyPhy 3 sobre estos nuevos alineamientos.

ADVERTENCIA: estos alineamientos son candidatos con integridad
codonica demostrada; la evaluacion de homologia ambigua, ortologia,
calidad del arbol y modelos estadisticos sigue pendiente.

7.7. Arboles de referencia transferidos a carpeta ``arbol``
---------------------------------------------------------

Directorio creado en el Mac::

    ANALYSES/dNdS_FINAL/HyPhy3_reanalysis/arbol/

Archivos presentes segun ``ls -l`` del usuario::

    Myotis_cactus_neutral_ML_RE-ROOTED.tre
    Myotis_cactus_neutral_ML_rooted.tre

El primero (10146 bytes) es **NEXUS de FigTree**, no Newick simple.
El segundo (3766 bytes) es **Newick**. Los archivos ``._*`` son
metadatos auxiliares de macOS, NO arboles, e ignorarlos.

Se leyeron correctamente con::

    Phylo.read("Myotis_cactus_neutral_ML_RE-ROOTED.tre", "nexus")
    Phylo.read("Myotis_cactus_neutral_ML_rooted.tre", "newick")

Auditoria confirmada por terminal:

* Ambos: 111 terminales, 111 identificadores unicos.
* Etiquetas compartidas: 111/111; diferencias: ninguna.
* Longitud total de ramas reportada: 6.88873 en ambos.
* ``rooted`` de Biopython: NEXUS ``True``; Newick ``False``.
  Esto es una propiedad de lectura/representacion, NO una prueba
  del enraizamiento biologico.
* Ambos presentan tres descendientes directos en la raiz
  de su representacion: grupos de 103, 7 y 1 taxones.
* El grupo de 1 contiene ``HLequAsi3`` (Equus asinus, por codigo
  inferido; confirmar nomenclatura con tabla de identificadores).
* Grupo de 7 incluye, entre los ejemplos, ``HLbosTau10``,
  ``HLeleMax2``, ``HLmanPen4``, ``HLneoVis2``, ``HLsorFum1``.
* Comparacion de topologias NO ENRAIZADAS mediante biparticiones:
  NEXUS 108; Newick 108; compartidas 108;
  exclusivas NEXUS 0; exclusivas Newick 0.
  Resultado: **topologias no enraizadas IDENTICAS**.

No se ha probado formalmente identidad de cada longitud de rama
individual, solo la longitud total. Tampoco se ha validado la
ubicacion biologica de la raiz. Ambos utilizan etiquetas codificadas
(p. ej. ``HLmyoMyo6``), no los identificadores ``Myotis_myotis``
de los FASTA NCL: el cotejo de nombres exactos sin diccionario
arroja 0 coincidencias y NO implica ausencia de las especies.

No utilizar indiscriminadamente el archivo NEXUS completo como
Newick: ``Phylo.parse(..., "newick")`` interpreto erroneamente
las instrucciones NEXUS/FigTree como 99 pseudoarboles.

7.8. Tabla del manuscrito: auditoria y cambios necesarios
--------------------------------------------------------

La tabla de 28 recuperaciones NCL proporcionada por el investigador
usa columnas ``Architecture``, ``Species``, ``Genus``,
``Age (Ma) of the common node (Ancestor)``, ``Bat family``,
``Architecture recovery``, ``TE recovery (%)``,
``TE identity (%)`` y ``RNA support``.

**No confundir** los hits de BLAST de esta tabla con la recuperacion
integra de CDS ni con el QC de los 25 alineamientos de HyPhy 3.

Correcciones/revisiones prioritarias:

#. ``Plecotus auritus``: cambiar ``RNA support`` de NO a YES:
   29 alineamientos para la union NCL-CHARLIE8 y 90 para exon 4-5
   (cinco BAM; comprobados). No implica validacion proteica.
#. ``Vespertilio murinus``: conservar YES, respaldado por 27
   alineamientos de union NCL-CHARLIE8 en cinco BAM.
#. ``Pipistrellus pygmaeus``: el YES original NO fue verificado
   mediante RNA local en la reconstruccion reciente. Marcar
   ``YES (source to verify)`` o ``not independently verified``
   hasta localizar y examinar la evidencia original; NO cambiar
   a NO afirmando ausencia demostrada de RNA.
#. ``Antrozous pallidus``: actualizar nota de recuperacion para
   reflejar correccion provisional del exon 4 de +42 nt;
   sin soporte RNA confirmado aqui.
#. ``Nyctalus aviator``: actualizar nota para reflejar CDS
   provisional intacta y exon 9 de 155 nt, con un codon gap
   respecto a otras especies, sin evidencia RNA confirmada.
#. ``Cistugo seabrae`` y ``Pipistrellus kuhlii``: conservar como
   hallazgos genomicos historicos pero agregar columna de
   elegibilidad: **excluidas de HyPhy 3 por QC de locus**.
#. ``Pipistrellus nathusii``: agregar nota **no incluida en las
   25 CDS actuales**: BLAST en hap2; falta ensamblaje hap2.
#. Los 25 taxones del nuevo alineamiento deben poder distinguirse
   claramente de los 28 hits historicos; 25/28 aceptados
   provisionalmente para alineamiento, 2 excluidos, 1 pendiente.
#. La identidad TE del BLAST original no debe reemplazarse por
   la identidad del segmento TE de la CDS reconstruida. Si se
   presenta esta ultima, introducir una nueva columna y definir
   referencia, longitud y metodo.
#. La cifra 36.37 Ma es un valor heredado de la tabla, no
   recalculado ni auditado durante esta etapa. No afirmar una
   edad nueva sin revisar arbol fechado/fuente original.
#. ``9/9 NCL exon matches`` solo indica hits; no garantiza
   que los nueve exones correspondan a una CDS valida.

En particular: ``P. pygmaeus`` puede tener evidencia transcriptomica
original documentada fuera del conjunto local; **no se sabe** con
la informacion disponible aqui si ese YES es correcto.

7.9. Punto EXACTO de reanudacion y normas de seguridad
-----------------------------------------------------

ESTADO AL CIERRE DEL 23 DE NOVIEMBRE DE 2026:

* 25 CDS NCL-CHARLIE8 provisionales reconstruidas y guardadas.
* Proteinas alineadas con MAFFT L-INS-i; 983 columnas.
* CHIMERA/HOST/TE codonicos separados y guardados para 25 especies.
* QC de integridad codonica aprobado; region HOST de baja
  complejidad todavia requiere analisis de sensibilidad.
* Dos arboles de 111 terminales en ``arbol/``: NEXUS y Newick,
  topologia no enraizada identica (108/108 splits).
* Todavia NO hay un arbol podado validado para las 25 especies.
* Todavia NO se han corrido nuevos tests HyPhy 3 sobre estas CDS.

SIGUIENTE PASO RECOMENDADO (uno solo, SOLO LECTURA):

#. Localizar el diccionario / mapeo fiable entre etiquetas de
   arbol ``HL...`` y especies binomiales de los FASTA NCL; comparar
   con identificadores de antiguos NEXUS NCL y los archivos de
   referencia del proyecto. No inferir a ciegas desde prefijos.

DESPUES (aun NO ejecutado):

#. Confirmar correspondencia exacta 25/25 y unicidad de etiquetas.
#. Podar una **copia nueva** del arbol escogido a 25 especies,
   conservar topologia y documentar longitudes de ramas y raiz.
#. Comparar topologia y longitud de ramas con arbol NCL historico.
#. Auditar regiones HOST ambiguas, decidir mascara de sensibilidad
   con criterios explicitos y reproducibles (sin borrado arbitrario).
#. Evaluar si el arbol de especies basado en sitios neutrales es
   apropiado para el test que se realizara; HyPhy puede reestimar
   longitudes de rama bajo su modelo de sustitucion.
#. Ejecutar HyPhy 3 en nuevo directorio de resultados y registrar
   version, parametros, arbol, datasets, JSON/LOG y correcciones
   multiples. No reinterpretar p nominal como q significativa.

REGLAS PRACTICAS:

* Un solo comando ejecutable por turno cuando se trabaje en Mac.
* Preferir auditorias de SOLO LECTURA antes de escribir.
* NO modificar ni sobrescribir FASTA, tablas, NEXUS, JSON y LOG
  historicos o provisionales previamente confirmados.
* Distinguir explicitamente: observacion / inferencia / desconocido.
* Mantener ``Pipistrellus_nathusii`` fuera de la matriz de 25
  hasta disponer de un locus y una CDS validos en el ensamblaje
  correcto. No convertir hits hap2 usando coordenadas hap1.
