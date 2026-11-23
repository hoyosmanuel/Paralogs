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
