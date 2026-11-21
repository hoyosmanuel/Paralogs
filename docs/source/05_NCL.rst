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

    
II. AUDITORIA DE NCL
==================================

como no encontré a Vespertilio (que corre en sentido contrario) murinus muy seguramente cometí un error.
así que vamos a probar de nuevo.

.. code-block:: bash

    nano auditar_NCL.py

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



RUM!

.. code-block:: bash

    python3 auditar_NCL.py --root ../../..
