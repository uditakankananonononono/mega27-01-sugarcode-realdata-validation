"""Inventory of sugarcode-ai modules and bio/ libraries with validation status, derived from sweep/VERDICTS.md and benchmarks/.
Status is 'validated' only for names in the curated EVIDENCE map below, each tied to a VERDICTS.md section and a committed result file; everything else is 'not yet validated'."""
import os, re, sys, subprocess, csv
SC=sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai'
V=open('sweep/VERDICTS.md').read().lower(); BN=' '.join(os.listdir('benchmarks')).lower()
head=subprocess.run(['git','-C',SC,'rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip()
rows=[]
EVIDENCE={('module','chem_descriptors'):'benchmarks/results.json (VERDICTS table)',('module','chem_similarity'):'benchmarks/results.json (VERDICTS table)',
('module','mit_offtarget'):'VERDICTS table; sweep/crispr_vs_crispor.py',('module','crisprscan_score'):'VERDICTS table + CRISPRscan finding',
('module','rna_nussinov'):'VERDICTS table; sweep/rna_vs_vienna.py',('module','profile_hmm'):'benchmarks/profile_hmm_cv.json, profile_hmm_cv_local.json',
('module','codon_opt'):'benchmarks/sweep_codon_cai.json',('module','pgx_guidelines'):'benchmarks/sweep_pgx_cpic.json',
('module','neohunter'):'VERDICTS neohunter section',('module','docking_studio'):'benchmarks/sweep_docking.json',
('module','deepsplice'):'discovery/splice_region_vus/cv_results_me.json, maxent_cmp.json',
('bio','codon'):'benchmarks/sweep_codon_cai.json',('bio','structures'):'benchmarks/sweep_sasa_freesasa.json',('bio','primer'):'benchmarks/sweep_primer_primer3.json',
('bio','phylo'):'benchmarks/sweep_phylo_dendropy.json',('bio','orf'):'benchmarks/sweep_orf_orfipy.json',('bio','vcf'):'benchmarks/sweep_vcf_pysam.json',
('bio','de'):'benchmarks/sweep_de_pydeseq2.json',('bio','rnaseq'):'benchmarks/sweep_de_pydeseq2.json',('bio','proteinprops'):'benchmarks/sweep_proteinprops_restriction.json',('bio','restriction'):'benchmarks/sweep_proteinprops_restriction.json',('bio','stockholm'):'benchmarks/sweep_stockholm_biopython.json',('bio','gff'):'benchmarks/sweep_gff_gffutils.json',('bio','newick'):'benchmarks/sweep_newick_biophylo.json',('bio','sam'):'benchmarks/sweep_sam_pysam.json',('bio','fastq'):'benchmarks/sweep_fastq_biopython.json',('bio','pdb'):'benchmarks/sweep_pdb_gemmi.json',('bio','pileup'):'benchmarks/sweep_pileup_pysam.json',('bio','motif'):'benchmarks/sweep_motif_biopython.json',('bio','pwm'):'benchmarks/sweep_motif_biopython.json (via bio.motif)',('bio','align'):'benchmarks/sweep_align_parasail.json',('bio','bed'):'benchmarks/sweep_bed_bioframe.json',('bio','genbank'):'benchmarks/sweep_genbank_biopython.json',('bio','sequence'):'benchmarks/sweep_sequence_biopython.json',('bio','kmer'):'benchmarks/sweep_kmer_sourmash.json',('bio','gstats'):'benchmarks/sweep_gstats_scipy.json',('bio','uniprot'):'benchmarks/sweep_uniprot_hgnc.json',('bio','entrez'):'benchmarks/sweep_entrez_hgnc.json',('bio','reactome'):'benchmarks/sweep_reactome_download.json',('bio','chembl'):'benchmarks/sweep_chembl_order.json',('bio','gnomad'):'benchmarks/sweep_gnomad_constraint.json',('bio','splice'):'benchmarks/sweep_splice_clinvar.json',('module','cfd_offtarget'):'benchmarks/sweep_cfd_doench.json',('module','str_scope'):'benchmarks/sweep_str_pytrf.json',('module','acmg_bayesian'):'benchmarks/sweep_acmg_erepo.json',('module','openclinvar'):'benchmarks/sweep_openclinvar.json',('module','crispr_opt'):'benchmarks/sweep_crispr_ontarget_fcres.json',('module','qsar_bench'):'benchmarks/sweep_qsar_descriptors.json',('module','evofold_4d'):'benchmarks/sweep_evofold_anm.json',('module','structural_biophysics'):'benchmarks/sweep_vina_score.json; benchmarks/sweep_sasa_freesasa_recheck14.json',('module','molecule_eval'):'benchmarks/sweep_molecule_eval.json',('module','riboswitch'):'benchmarks/sweep_nussinov.json',('module','prime_design'):'benchmarks/sweep_prime_design.json',('module','promoter_lib'):'benchmarks/sweep_stability_promoter.json',('module','stability_ai'):'benchmarks/sweep_stability_promoter.json',('bio','fasta'):'benchmarks/sweep_fasta.json',('module','crisprater'):'benchmarks/sweep_crisprater.json',('module','omega_stats'):'benchmarks/sweep_omega_stats.json',('module','microbiome_exp'):'benchmarks/sweep_microbiome.json',('module','rna_decoder'):'benchmarks/sweep_rna_decoder.json',('module','evidence_mining'):'benchmarks/sweep_evidence_mining.json',('module','rarenet_ai'):'benchmarks/sweep_rarenet.json',('module','oncocircuit'):'benchmarks/sweep_oncocircuit.json'}
def loc(p):
    n=0
    for dp,_,fs in os.walk(p) if os.path.isdir(p) else [(os.path.dirname(p),[],[os.path.basename(p)])]:
        for f in fs:
            if f.endswith('.py'): n+=sum(1 for _ in open(os.path.join(dp,f),errors='ignore'))
    return n
tests=' '.join(open(os.path.join(dp,f),errors='ignore').read() for dp,_,fs in os.walk(SC+'/tests') for f in fs if f.endswith('.py'))
for kind,base in (('module',SC+'/src/sugarcode/modules'),('bio',SC+'/src/sugarcode/bio')):
    for n in sorted(os.listdir(base)):
        if n.startswith('_') or n=='data': continue
        name=n[:-3] if n.endswith('.py') else n
        if not (n.endswith('.py') or os.path.isdir(os.path.join(base,n))): continue
        pat=name.lower(); keys=[pat]+([ 'bio.'+pat, pat+'.py'] if kind=='bio' else [])
        ev=EVIDENCE.get((kind,name),'')
        nt=len(re.findall(r'(?:modules|bio)[./]'+re.escape(name)+r'\b',tests))
        rows.append((kind,name,loc(os.path.join(base,n)),nt,'validated' if ev else 'not yet validated',ev))
with open('manifests/module_inventory.tsv','w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['kind','name','python_lines','test_references','status','evidence']); w.writerows(rows)
print('sugarcode-ai',head,len(rows),sum(r[4]=='validated' for r in rows))
