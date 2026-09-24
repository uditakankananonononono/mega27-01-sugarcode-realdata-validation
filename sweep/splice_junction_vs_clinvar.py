"""bio.splice.junction_map vs ClinVar canonical splice-site SNVs (c.N+1/+2, c.N-1/-2), live NCBI E-utilities 2026-09-25.
For each variant the named c. position must be a junction in sugarcode's map and the reference base at the named intronic offset
in sugarcode's window must equal ClinVar's reference allele. ClinVar names use the MANE Select transcript."""
import sys, json, re
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import splice as S, entrez as E
GENES=['TP53','BRCA1','BRCA2','MLH1','MSH2','APC','PTEN','ATM','CDH1','PALB2','CFTR','NF1']
pat=re.compile(r'^(NM_\d+\.\d+)\((\w+)\):c\.(\d+)([+-])([12])([ACGT])>([ACGT])$')
res=[]; per={}
for g in GENES:
    jm=S.junction_map(g)
    if jm.get('status')!='ok': per[g]={'status':jm.get('status')}; continue
    ids=E.esearch('clinvar',f'{g}[gene] AND (splice_donor_variant[Molecular consequence] OR splice_acceptor_variant[Molecular consequence])',retmax=40)
    sm=E.esummary('clinvar',ids)
    ok=tot=junc=0; tx=set(); bad=[]
    for u,d in sm.items():
        t=d.get('title',''); m=pat.match(t.split(' ')[0])
        if not m: continue
        nm,gg,n,sign,k,ref,alt=m.groups(); n=int(n); k=int(k); tx.add(nm); tot+=1
        if sign=='+':
            w=jm['donors'].get(n); base=w[2+k] if w else None
        else:
            w=jm['acceptors'].get(n); base=w[14-k] if w else None  # window = 14 intronic + first exonic base (seq[c-15:c], 1-based c)
        junc+= w is not None; good= base==ref; ok+=good
        if not good and len(bad)<3: bad.append({'title':t.split(' ')[0],'window':w,'base':base})
    per[g]={'accession':jm['accession'],'clinvar_transcripts':sorted(tx),'n':tot,'junction_found':junc,'ref_base_ok':ok,'examples_bad':bad,'n_exons_map':len(jm['exons'])}
    print(g,per[g]['accession'],sorted(tx),tot,junc,ok,bad[:1])
out={'reference':'ClinVar (esearch/esummary, canonical +/-1,2 SNVs, MANE-named)','per_gene':per,
 'total':sum(v.get('n',0) for v in per.values()),'junction_found':sum(v.get('junction_found',0) for v in per.values()),'ref_base_ok':sum(v.get('ref_base_ok',0) for v in per.values())}
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else 'benchmarks/sweep_splice_clinvar.json','w'),indent=1)
print(out['total'],out['junction_found'],out['ref_base_ok'])
