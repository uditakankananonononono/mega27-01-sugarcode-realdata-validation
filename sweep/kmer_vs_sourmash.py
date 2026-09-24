"""bio.kmer vs sourmash on 6 individually fetched coronavirus genomes (NCBI nuccore FASTA, 2026-09-25):
SARS-CoV-2 NC_045512.2, RaTG13 MN996532.2, pangolin GX-P4L MT040333.1, SARS-CoV Tor2 NC_004718.3, bat ZC45 MG772933.1, MERS NC_019843.3.
Exact canonical k-mer sets (k=21) vs sourmash MinHash(scaled=1, all hashes kept); minimizer-sketch Jaccard estimates vs sourmash FracMinHash(scaled=100)."""
import sys, glob, json, itertools, os
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import kmer as K
import sourmash
seqs={os.path.basename(f)[:-3]:''.join(l.strip() for l in open(f) if not l.startswith('>')).upper() for f in sorted(glob.glob('data/kmer/*.fa'))}
k=21; out={'reference':'sourmash %s'%sourmash.VERSION,'k':k,'genomes':{},'pairs':[]}
mh1={}; mh100={}; sets={}
for n,s in seqs.items():
    a=sourmash.MinHash(n=0,ksize=k,scaled=1); a.add_sequence(s,force=True); mh1[n]=a
    b=sourmash.MinHash(n=0,ksize=k,scaled=100); b.add_sequence(s,force=True); mh100[n]=b
    sets[n]=K.kmer_set(s,k)
    out['genomes'][n]={'len':len(s),'distinct_kmers_sugarcode':len(sets[n]),'distinct_hashes_sourmash':len(a),'equal':len(sets[n])==len(a)}
err_sc=[]; err_sm=[]
for x,y in itertools.combinations(seqs,2):
    ex=K.compare_sequences(seqs[x],seqs[y],k); sm=mh1[x].jaccard(mh1[y])
    sk=K.compare_sketches(seqs[x],seqs[y],k,w=100)  # expected sketch density ~2/(w+1)
    est100=mh100[x].jaccard(mh100[y])
    out['pairs'].append({'a':x,'b':y,'jaccard_exact_sugarcode':ex['jaccard'],'jaccard_sourmash_scaled1':sm,'abs_diff':abs(ex['jaccard']-sm),
        'containment_a_in_b_sugarcode':ex['containment_a_in_b'],'containment_sourmash_scaled1':mh1[x].contained_by(mh1[y]),
        'mash_distance_sugarcode':K.mash_distance(ex['jaccard'],k),'minimizer_w100_estimate':sk['jaccard_estimate'],'minimizer_sketch_sizes':[sk['sketch_a'],sk['sketch_b']],
        'sourmash_scaled100_estimate':est100,'sourmash_sketch_sizes':[len(mh100[x]),len(mh100[y])]})
    err_sc.append(abs(sk['jaccard_estimate']-ex['jaccard'])); err_sm.append(abs(est100-ex['jaccard']))
p=out['pairs']
out['summary']={'n_pairs':len(p),'max_abs_diff_exact_jaccard':max(q['abs_diff'] for q in p),'max_abs_diff_containment':max(abs(q['containment_a_in_b_sugarcode']-q['containment_sourmash_scaled1']) for q in p),
 'distinct_counts_equal':sum(g['equal'] for g in out['genomes'].values()),'sketch_mae_sugarcode_minimizer_w100':sum(err_sc)/len(err_sc),'sketch_mae_sourmash_scaled100':sum(err_sm)/len(err_sm)}
json.dump(out,open('benchmarks/sweep_kmer_sourmash.json','w'),indent=1)
print(json.dumps(out['summary'],indent=1)); print(json.dumps(out['genomes']))
for q in p: print(q['a'],q['b'],round(q['jaccard_exact_sugarcode'],4),round(q['jaccard_sourmash_scaled1'],4),round(q['minimizer_w100_estimate'],4),round(q['sourmash_scaled100_estimate'],4))
