import json,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
B='../benchmarks/'; D='../discovery/splice_region_vus/'
rna=json.load(open(B+'sweep_rna.json')); rv=json.load(open(B+'sweep_rna_vienna_methods.json'))
fams=list(rna); m=['Nussinov','MFE','centroid','MEA g=1']; vals=[[rna[f]['nussinov_f1'],rv[f]['mfe'],rv[f]['cent'],rv[f]['mea1']] for f in fams]
fig,ax=plt.subplots(figsize=(6,3)); x=np.arange(len(fams))
for i,n in enumerate(m): ax.bar(x+i*0.2-0.3,[v[i] for v in vals],0.2,label=n)
ax.set_xticks(x); ax.set_xticklabels(fams); ax.set_ylabel('base-pair F1'); ax.legend(fontsize=7); ax.set_title('Rfam seed structures (n=25/family)',fontsize=9)
fig.tight_layout(); fig.savefig('figs/rna_f1.png',dpi=200)
cal=json.load(open(D+'cnn_calibration.json'))['thresholds']; t=list(cal)
fig,ax=plt.subplots(figsize=(6,3))
ax.plot([cal[k]['fpr'] for k in t],[cal[k]['sensitivity'] for k in t],'o-',label='sensitivity vs FPR')
for k in t: ax.annotate(k,(cal[k]['fpr'],cal[k]['sensitivity']),fontsize=7)
ax.set_xscale('log'); ax.set_xlabel('false-positive rate (gene-grouped CV)'); ax.set_ylabel('sensitivity'); ax.set_title('Splice CNN operating points',fontsize=9)
fig.tight_layout(); fig.savefig('figs/cnn_operating.png',dpi=200)
p=json.load(open(B+'sweep_pgx_cpic.json')); g=['CYP2C19','CYP2D6']
fig,ax=plt.subplots(figsize=(6,3)); x=np.arange(2)
ax.bar(x-0.2,[p['before_fix'][k]['agree']/p['before_fix'][k]['of'] for k in g],0.4,label='before fix')
ax.bar(x+0.2,[p[k]['agree']/p[k]['cpic_diplotypes'] for k in g],0.4,label='after fix')
ax.set_xticks(x); ax.set_xticklabels(g); ax.set_ylabel('fraction of all CPIC diplotypes\nmatched (uncovered = miss)',fontsize=8); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig('figs/pgx_cpic.png',dpi=200)
