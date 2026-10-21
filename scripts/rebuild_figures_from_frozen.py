#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description='Rebuild data-driven manuscript and SI figures from frozen result tables.')
ap.add_argument('--results', help='Directory containing compatible result CSVs. Defaults to expected/.')
a=ap.parse_args()
base=Path(a.results) if a.results else ROOT/'expected'
out=ROOT/'figures'/'rebuilt'; out.mkdir(parents=True,exist_ok=True)

# Figure 2: cohort benchmark
c=pd.read_csv(base/'cohort_benchmark.csv'); x=np.arange(len(c)); w=.25
fig,ax=plt.subplots(figsize=(10,4.8))
ax.bar(x-w,c.rrr_mean,w,yerr=c.rrr_sd,label='RRR',capsize=2)
ax.bar(x,c.pcr_mean,w,yerr=c.pcr_sd,label='PCR',capsize=2)
ax.bar(x+w,c.ridge_mean,w,yerr=c.ridge_sd,label='Ridge',capsize=2)
ax.axhline(0,lw=.8); ax.set_xticks(x,c.cohort,rotation=35,ha='right'); ax.set_ylabel('Held-out multivariate $R^2$'); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(out/'Figure2_rebuilt.png',dpi=300); fig.savefig(out/'Figure2_rebuilt.pdf'); plt.close(fig)

# Figure 3: sample-size calibration thresholds
r=pd.read_csv(base/'rank_thresholds.csv')
labels={35:'Kang',60:'Wandro',72:'Jacobs',77:'Erawijantari',176:'Franzosa',222:'He',230:'Wang',278:'Yachida'}
fig,ax=plt.subplots(figsize=(9,5.4))
ax.scatter(r['n'],100*r['per_direction_share'],s=55)
for _,z in r.iterrows():
    ax.annotate(labels.get(int(z['n']),str(int(z['n']))),(z['n'],100*z['per_direction_share']),xytext=(8,8),textcoords='offset points')
ax.set_xlabel('Rank-selection sample size'); ax.set_ylabel('Per-direction signal share needed for\nmedian selected rank ≥ 5 (%)')
ax2=ax.twinx(); lo,hi=ax.get_ylim(); ax2.set_ylim(lo/10,hi/10); ax2.set_ylabel('Total population signal on simulated grid')
fig.tight_layout(); fig.savefig(out/'Figure3_rebuilt.png',dpi=300); fig.savefig(out/'Figure3_rebuilt.pdf'); plt.close(fig)

# Figure 4: functional compression and seed robustness
n=pd.read_csv(base/'hcom2_nmf_threshold.csv'); s=pd.read_csv(base/'hcom2_support_strength.csv')
fig,axs=plt.subplots(1,2,figsize=(10.5,4.1))
axs[0].plot(n['rank'],n.p3_pearson,'o-',label='ordinary NMF'); axs[0].plot(s['latent_dimension'],s.p3_pearson,'o-',label='support-preserving strengths'); axs[0].axhline(.772953,ls='--',lw=1,label='full model'); axs[0].set(xlabel='Representation dimension',ylabel='Passage-3 Pearson'); axs[0].legend(frameon=False,fontsize=8)
axs[1].plot(n['rank'],n.seed_pass_fraction,'o-',label='ordinary NMF'); axs[1].plot(s['latent_dimension'],s.seed_pass_fraction,'o-',label='support-preserving strengths'); axs[1].set(xlabel='Representation dimension',ylabel='Fraction of seeds passing',ylim=(-.03,1.03)); axs[1].axhline(.5,ls=':',lw=.8)
fig.tight_layout(); fig.savefig(out/'Figure4_rebuilt.png',dpi=300); fig.savefig(out/'Figure4_rebuilt.pdf'); plt.close(fig)

# Figure 5: Ho support controls
h=pd.read_csv(base/'ho_support_controls.csv'); fig,ax=plt.subplots(figsize=(8.5,4.3)); xx=np.arange(len(h)); bars=ax.bar(xx,h.primary_error)
ax.set_xticks(xx,h.support_model,rotation=28,ha='right'); ax.set_ylabel('Primary author-style error')
# Mark rewire exploratory by hatch when found
for bar,name in zip(bars,h.support_model):
    if 'Degree-preserving' in name: bar.set_hatch('//')
fig.tight_layout(); fig.savefig(out/'Figure5_rebuilt.png',dpi=300); fig.savefig(out/'Figure5_rebuilt.pdf'); plt.close(fig)

# SI Fig S1: niche vs predictive dimension
q=pd.read_csv(base/'niche_dimensions.csv'); fig,ax=plt.subplots(figsize=(6,4.5)); ax.scatter(q.niche_mean,q.rrr_dimension)
for _,z in q.iterrows(): ax.annotate(z.cohort,(z.niche_mean,z.rrr_dimension),xytext=(5,4),textcoords='offset points',fontsize=8)
ax.set_xlabel('Effective niche dimension'); ax.set_ylabel('Median selected RRR dimension'); fig.tight_layout(); fig.savefig(out/'Supplementary_Figure_S1_rebuilt.png',dpi=300); fig.savefig(out/'Supplementary_Figure_S1_rebuilt.pdf'); plt.close(fig)

# SI Fig S2: ordinary NMF seed pass fraction at high ranks
fig,ax=plt.subplots(figsize=(6,4)); ax.plot(n['rank'],n.seed_pass_fraction,'o-'); ax.axvline(52,ls=':',lw=.8); ax.set(xlabel='Ordinary NMF rank',ylabel='Fraction of seeds satisfying criterion',ylim=(-.03,1.03)); fig.tight_layout(); fig.savefig(out/'Supplementary_Figure_S2_rebuilt.png',dpi=300); fig.savefig(out/'Supplementary_Figure_S2_rebuilt.pdf'); plt.close(fig)

# SI Fig S3: weak-edge ablation
wdf=pd.read_csv(base/'hcom2_weak_edge_ablation.csv'); y=pd.to_numeric(wdf.retained_improvement,errors='coerce'); fig,ax=plt.subplots(figsize=(6,4)); ax.plot(100*wdf.level,y,'o-'); ax.axhline(.95,ls='--',lw=.8); ax.set(xlabel='Weakest positive edges removed (%)',ylabel='Retained predictive improvement'); fig.tight_layout(); fig.savefig(out/'Supplementary_Figure_S3_rebuilt.png',dpi=300); fig.savefig(out/'Supplementary_Figure_S3_rebuilt.pdf'); plt.close(fig)
print('Wrote',out)
