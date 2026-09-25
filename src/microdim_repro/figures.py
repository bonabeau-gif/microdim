from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


def figure1(out):
    fig,ax=plt.subplots(figsize=(7.2,6.0));ax.axis('off')
    boxes=[(.16,.68,.68,.22,'Mechanistic consumption dimension','rank/subspace of uptake-driven competition\nin an explicit consumer-resource model\nTheorems 1–2'),
           (.16,.38,.68,.22,'Effective niche dimension','latent environmental axes organizing abundance\nunder a specified ecological model\nSrinivasan et al. (2026)'),
           (.16,.08,.68,.22,'Predictive dimension','rank or latent dimension selected because it\nimproves held-out prediction or compression\npaired-cohort benchmark')]
    for x,y,w,h,title,sub in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.01',fill=False,linewidth=1.5))
        ax.text(x+w/2,y+h*.67,title,ha='center',va='center',weight='bold',fontsize=13)
        ax.text(x+w/2,y+h*.36,sub,ha='center',va='center',fontsize=9)
    ax.annotate('',xy=(.5,.68),xytext=(.5,.60),arrowprops=dict(arrowstyle='->',lw=1.3))
    ax.text(.515,.635,'requires an identification bridge\nto uptake or interactions',va='center',fontsize=8)
    ax.annotate('',xy=(.5,.38),xytext=(.5,.30),arrowprops=dict(arrowstyle='->',lw=1.3))
    ax.text(.515,.335,'requires ecological generative\nassumptions and validation',va='center',fontsize=8)
    ax.text(.5,.015,'All three can be low; equality is an additional claim, not a consequence of low rank alone.',ha='center',fontsize=8)
    fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure2(reference_csv,out):
    d=pd.read_csv(reference_csv);models=['RRR','PCR','Ridge'];x=np.arange(len(d));w=.25
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    for j,m in enumerate(models):
        ax.bar(x+(j-1)*w,d[f'{m}_mean_R2'],w,yerr=d[f'{m}_sd_R2'],capsize=2,label=m)
    ax.axhline(0,linewidth=1);ax.set_xticks(x,d.cohort,rotation=40,ha='right');ax.set_ylabel('Outer-fold predictive R² (mean ± SD)');ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure3(primary_csv,out):
    d=pd.read_csv(primary_csv);fig,ax=plt.subplots(figsize=(7.2,5.0))
    for cohort in ['Wandro preterms','Franzosa IBD','Yachida CRC']:
        s=d[(d.cohort==cohort)&(d.true_rank==10)].sort_values('per_direction_variance_share')
        ax.plot(100*s.per_direction_variance_share,s.fraction_true_directions_retained,marker='o',label=f'{cohort} (rank-selection n≈{int(s.n_rank_selection.iloc[0])})')
    ax.set_xlabel('Variance share per true direction (%)');ax.set_ylabel('Mean fraction of true directions retained');ax.set_ylim(-.02,1.02);ax.legend(fontsize=8,frameon=False)
    fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure4(curve_csv,out):
    d=pd.read_csv(curve_csv);fig,ax=plt.subplots(figsize=(6.4,4.8));ax.plot(d.eta,d['median'],marker='o',label='median');ax.fill_between(d.eta,d.q025,d.q975,alpha=.25,label='2.5–97.5 percentile')
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Resource time-scale factor eta');ax.set_ylabel('Maximum log-abundance error');ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure_s1(primary_csv,out):
    d=pd.read_csv(primary_csv);d=d[d.true_rank==10]
    order=['Kang autism','Wandro preterms','Jacobs IBD families','Erawijantari gastric cancer','Franzosa IBD','He infants','Wang ESRD','Yachida CRC']
    mat=np.array([[d[(d.cohort==c)&(d.population_R2==r)].median_selected_rank.iloc[0] for r in [0.02,.05,.10,.15,.20,.30,.50]] for c in order])
    fig,ax=plt.subplots(figsize=(7.3,5.2));im=ax.imshow(mat,aspect='auto');ax.set_yticks(range(8),[f'{c.split()[0]} (n={int(d[d.cohort==c].n_rank_selection.iloc[0])})' for c in order]);ax.set_xticks(range(7),['.02','.05','.10','.15','.20','.30','.50']);ax.set_xlabel('Population R² of the rank-10 signal');ax.set_ylabel('Cohort-sized rank-selection dataset')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):ax.text(j,i,f'{mat[i,j]:g}',ha='center',va='center')
    fig.colorbar(im,ax=ax,label='Median selected rank');fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure_s2(spec_csv,out):
    d=pd.read_csv(spec_csv);fig,ax=plt.subplots(figsize=(6.4,4.8))
    for spec in ['equal','geometric']:
        s=d[d.spectrum==spec].sort_values('population_R2');ax.plot(s.population_R2,s.median_selected_rank,marker='o',label=spec)
    ax.axhline(10,ls='--',lw=1);ax.set_xlabel('Population variance share carried by rank-10 signal');ax.set_ylabel('Median selected rank');ax.legend(frameon=False);fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure_s3(micom_csv,out):
    d=pd.read_csv(micom_csv).sort_values('sample_id');x=np.arange(len(d));w=.25;fig,ax=plt.subplots(figsize=(7.2,4.8))
    for j,(col,label) in enumerate([('binary_participation_rank','Binary uptake'),('flux_participation_rank','Per-biomass flux'),('weighted_participation_rank','Abundance-weighted flux')]):
        ax.bar(x+(j-1)*w,d[col],w,label=label)
    ax.set_xticks(x,d.sample_id);ax.set_xlabel('Community');ax.set_ylabel('Participation-ratio effective dimension');ax.legend(frameon=False);fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)


def figure_s4(sens_csv,out):
    d=pd.read_csv(sens_csv);fig,ax=plt.subplots(figsize=(7.2,4.8))
    for sid,g in d.groupby('sample_id'):
        ax.plot(g.flux_cutoff,g.participation_rank,marker='o',label=sid)
    ax.set_xscale('log');ax.set_xlabel('Minimum absolute uptake flux retained');ax.set_ylabel('Per-biomass participation-ratio dimension');ax.legend(title='Community',frameon=False);fig.tight_layout();fig.savefig(out,dpi=300);plt.close(fig)
