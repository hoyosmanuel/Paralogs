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

Estado inicial:

* Directorio del proyecto comprobado.
* Espacio HyPhy3_reanalysis creado.
* Instalación de HyPhy 3 todavía no realizada.
* Compatibilidad y resultados todavía no verificados.

Objetivos:

1. Registrar versión y procedencia del software.
2. Verificar los archivos de entrada originales.
3. Ejecutar controles de calidad.
4. Reproducir primero los análisis piloto.
5. Comparar resultados de HyPhy 2 y HyPhy 3.
6. Ejecutar los análisis restantes después de validar el piloto.
7. Aplicar BH-FDR y comprobar los resultados.
8. Regenerar las tablas para el manuscrito.


Reglas de trabajo
-----------------

* Ejecutar un paso a la vez.
* Documentar comandos, salidas, errores y correcciones.
* No sobrescribir los análisis originales.
* Distinguir hechos comprobados de hipótesis e inferencias.
* No modificar las secuencias para obtener resultados deseados.
* Conservar scripts, parámetros, versiones y resultados.
* Actualizar este registro después de cada etapa verificada.


Registro cronológico
--------------------

2026-11-22

* Confirmado el directorio ANALYSES/dNdS_FINAL.
* Creado HyPhy3_reanalysis con documentation, scripts,
  results y validation.
* Preparada la documentación inicial del nuevo proyecto.
