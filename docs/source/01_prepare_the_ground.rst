I. Set Directories
==================

1) Set ANALYSES Directories: 
----------------------------------

.. code-block:: bash
  
  /Volumes/Expansion/project3/bat_HTF_genomic_analysis
  mkdir -p ANALYSES/{FAU,ZFP28,ZKSCAN7,FBXO3,NCL,TMPO,ZNF24,AP1AR,KATNAL1}

  for gene in FAU ZFP28 ZKSCAN7 FBXO3 NCL TMPO ZNF24 AP1AR KATNAL1; do
    mkdir -p ANALYSES/"$gene"/{scripts,intermediate,sequences,results}
  done

2) Set SPECIES Directories: 
----------------------------------

Aquí adentro de múltiples maneras hice carpetas de las especies de interés las cuales son:

::

  Antrozous_pallidus
  Artibeus_intermedius
  Artibeus_lituratus
  Aselliscus_stoliczkanus
  Brachyphylla_cavernarum
  Carollia_perspicillata
  Centurio_senex
  Choeroniscus_minor
  Cistugo_seabrae
  Corynorhinus_mexicanus
  Corynorhinus_townsendii
  Craseonycteris_thonglongyai
  Cynopterus_sphinx
  Desmodus_rotundus
  Diaemus_youngii
  Diphylla_ecaudata
  Ectophylla_alba
  Eonycteris_spelaea
  Eptesicus_fuscus
  Eptesicus_nilssonii
  Erophylla_bombifrons
  Eumops_nanus
  Furipterus_horrens
  Glossophaga_mutica
  Glossophaga_soricina
  Glyphonycteris_daviesi
  Hipposideros_abae
  Hipposideros_armiger
  Hipposideros_caffer
  Hipposideros_cyclops
  Hipposideros_jonesi
  Hipposideros_larvatus
  Hipposideros_swinhoei
  Hypsignathus_monstrosus
  Lasiurus_ega
  Leptonycteris_yerbabuenae
  Lionycteris_spurrelli
  Lonchorhina_inusitata
  Macrophyllum_macrophyllum
  Macrotus_waterhousii
  Megaderma_spasma
  Micronycteris_megalotis
  Miniopterus_australis
  Miniopterus_natalensis
  Miniopterus_schreibersii
  Molossus_alvarezi
  Molossus_molossus
  Molossus_nigricans
  Mops_condylurus
  Mormoops_megalophylla
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
  Mystacina_tuberculata
  Myzopoda_aurita
  Natalus_tumidirostris
  Noctilio_leporinus
  Nyctalus_aviator
  Nycteris_thebaica
  Phyllostomus_discolor
  Phyllostomus_hastatus
  Pipistrellus_kuhlii
  Pipistrellus_nathusii
  Pipistrellus_pygmaeus
  Platyrrhinus_guianensis
  Plecotus_auritus
  Rhinolophus_affinis
  Rhinolophus_ferrumequinum
  Rhinolophus_foetidus
  Rhinolophus_hipposideros
  Rhinolophus_pearsonii
  Rhinolophus_perniger_lanosus
  Rhinolophus_sedulus
  Rhinolophus_sinicus
  Rhinolophus_trifoliatus
  Rhinophylla_pumilio
  Rhinopoma_microphyllum
  Rhinopoma_muscatellum
  Rhynchonycteris_naso
  Rousettus_aegyptiacus
  Saccopteryx_bilineata
  Saccopteryx_leptura
  Tadarida_brasiliensis
  Taphozous_melanopogon
  Thyroptera_tricolor
  Trachops_cirrhosus
  Triaenops_persicus
  Trinycteris_nicefori
  Uroderma_convexum
  Vampyressa_thyone
  Vespertilio_murinus
