#!/usr/bin/env python3
"""Reproduce PLOS ONE V8.7.6 Tables 1-5, data-derived main Figs 3-5,
and Supporting Figs S1-S6 from the frozen derived CSVs in S2 File.

Figs 1-2 are conceptual/design schematics and are therefore not generated from
numerical CSVs. All comparisons are descriptive; this script performs no
significance testing and does not treat pairwise distances as independent
experimental replicates.
"""
from pathlib import Path
import argparse, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

FAMILY_ORDER=['count_plus_initial','dynamic_scale','seriousness']
FAMILY_LABEL={'count_plus_initial':'Count + initial','dynamic_scale':'Dynamic scale','seriousness':'Seriousness'}
BLOCK_ORDER=['count_plus_initial','dynamic_scale','seriousness','temporal_path']
BLOCK_LABEL={'count_plus_initial':'Count + initial','dynamic_scale':'Dynamic scale','seriousness':'Seriousness','temporal_path':'Temporal path'}
ARMS=['P1','P4']
P1='#6B7280'; P4='#0072B2'; LINK='#CBD5E1'; TEAL='#2F7D6D'; PURPLE='#6A4C93'; GOLD='#B07A18'; GRID='#E5E7EB'; INK='#20262E'

def mean_arm(df,col): return df.groupby('arm')[col].mean().reindex(ARMS)
def popsd_arm(df,col): return df.groupby('arm')[col].agg(lambda x: np.std(x.to_numpy(float), ddof=0)).reindex(ARMS)
def clean(ax,axis='y'):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    if axis in ('y','both'): ax.grid(axis='y',color=GRID,lw=.7,alpha=.75,zorder=0)
    if axis in ('x','both'): ax.grid(axis='x',color=GRID,lw=.7,alpha=.75,zorder=0)

def seed_from_identity(x): return int(re.search(r'seed(\d+)',str(x)).group(1))
def pairify(df):
    x=df.copy()
    if 'seed_i' not in x.columns:
        x['seed_i']=x['identity_i'].map(seed_from_identity); x['seed_j']=x['identity_j'].map(seed_from_identity)
    x['seed_low']=x[['seed_i','seed_j']].min(axis=1).astype(int); x['seed_high']=x[['seed_i','seed_j']].max(axis=1).astype(int)
    x['pair_key']=list(zip(x.seed_low,x.seed_high)); return x

def table2(d):
    rows=[]
    loc=pd.read_csv(d/'G2_LOCATION_PAIRWISE.csv'); cv=pd.read_csv(d/'G2_MARGINAL_SCALE_CV.csv')
    pa=pd.read_csv(d/'G2_FACTOR_PRINCIPAL_ANGLES.csv'); prj=pd.read_csv(d/'G2_PROJECTION_DISTANCE.csv'); proc=pd.read_csv(d/'G2_PROCRUSTES_COMPARISON.csv')
    ba=pd.read_csv(d/'G2_BLOCK_DELTA_ALLOCATION.csv'); rho=pd.read_csv(d/'G2_DELTA_RHO_PAIRWISE.csv')
    def add(metric,vals,role): rows.append([metric,vals['P1'],vals['P4'],role])
    add('location L2',mean_arm(loc,'location_l2_distance'),'marginal reproducibility')
    add('marginal-scale CV',mean_arm(cv,'cv_marginal_sd'),'marginal reproducibility')
    for fam in FAMILY_ORDER: add(f'principal-angle sin_F: {fam}',mean_arm(pa[pa.family=='F_'+fam],'sin_F'),'local subspace geometry')
    for fam in FAMILY_ORDER: add(f'projection distance: {fam}',mean_arm(prj[prj.family=='F_'+fam],'projection_frobenius_distance'),'local subspace geometry')
    for fam in FAMILY_ORDER: add(f'Procrustes: {fam}',mean_arm(proc[proc.family=='F_'+fam],'raw_procrustes_distance'),'local subspace geometry')
    ba=ba.assign(abs_delta=ba.delta_allocation.abs())
    for block in BLOCK_ORDER: add(f'block allocation delta: {block}',mean_arm(ba[ba.semantic_block==block],'abs_delta'),'aggregate allocation reproducibility')
    add('delta rho',mean_arm(rho,'delta_rho'),'aggregate allocation reproducibility')
    return pd.DataFrame(rows,columns=['metric','P1','P4','interpretive role'])

def table3(d):
    metrics=[]
    def add(name,fn,col,filt=None):
        x=pd.read_csv(d/fn)
        if filt is not None: x=x.query(filt)
        m=mean_arm(x,col); s=popsd_arm(x,col); ch=(m['P4']/m['P1']-1)*100
        metrics.append([name,m['P1'],s['P1'],m['P4'],s['P4'],ch])
    add('Sigma relative Frobenius','G3_SIGMA_PAIRWISE.csv','relative_frobenius')
    add('correlation RMSE','G3_CORRELATION_CONSEQUENCE.csv','correlation_rmse')
    add('spectrum L2','G3_SPECTRUM_CONSEQUENCE.csv','spectrum_l2_distance')
    for b in BLOCK_ORDER: add(f'block relative Frobenius: {b}','G3_BLOCK_COVARIANCE_CONSEQUENCE.csv','block_relative_frobenius',f"block == '{b}'")
    add('B (AR1) relative Frobenius','G3_COMPONENT_PAIRWISE.csv','relative_frobenius_B')
    add('L (low-rank) relative Frobenius','G3_COMPONENT_PAIRWISE.csv','relative_frobenius_L')
    return pd.DataFrame(metrics,columns=['metric','P1 mean','P1 SD','P4 mean','P4 SD','Descriptive P4 vs P1 change (%)'])

def table5(d):
    x=pd.read_csv(d/'V8_P4_VS_P1_DESCRIPTIVE_COMPARISON.csv')
    specs=[('Stage I representation','location L2','location_l2_P4_reduction_pct'),('Stage I representation','principal-angle sin_F','principal_angle_sin_F_P4_reduction_pct'),('Stage I representation','projection Frobenius','projection_frobenius_P4_reduction_pct'),('Stage I representation','raw Procrustes','raw_procrustes_P4_reduction_pct'),('Stage II covariance','covariance relative Frobenius','covariance_reproducibility_distance_P4_reduction_pct'),('Stage II covariance','correlation RMSE','correlation_rmse_P4_reduction_pct'),('Stage II covariance','normalized-spectrum L2','normalized_spectrum_l2_P4_reduction_pct'),('Stage II component','L relative Frobenius','L_relative_frobenius_P4_reduction_pct'),('Known-target accuracy','mean-vector error','mean_error_to_target_P4_reduction_pct'),('Known-target accuracy','covariance-target error','covariance_error_to_target_P4_reduction_pct'),('Known-target accuracy','analytic Gaussian KL','analytic_gaussian_kl_P4_reduction_pct')]
    rows=[]
    for layer,name,col in specs:
        v=x[col].astype(float); rows.append([layer,name,int((v>0).sum()),float(v.min()),float(v.max())])
    return pd.DataFrame(rows,columns=['metric layer','metric','P4-lower cells (of 6)','min reduction (%)','max reduction (%)'])

def pair_panel(ax,df,col,title,ylabel):
    q=pairify(df); pv=q.pivot(index='pair_key',columns='arm',values=col); matched=pv.dropna()
    for j,(_,r) in enumerate(matched.iterrows()):
        jit=(j-(len(matched)-1)/2)*.014
        ax.plot([0+jit,1+jit],[r.P1,r.P4],color=LINK,lw=.85,zorder=1)
        ax.scatter(0+jit,r.P1,s=24,color=P1,edgecolor='white',lw=.3,zorder=3); ax.scatter(1+jit,r.P4,s=24,color=P4,edgecolor='white',lw=.3,zorder=3)
    if set(['P1','P4']).issubset(pv.columns):
        extra=pv[pv.P1.isna() & pv.P4.notna()].P4
        for j,v in enumerate(extra): ax.scatter(1+(j-(len(extra)-1)/2)*.018,v,s=27,facecolor='white',edgecolor=P4,lw=1,zorder=3)
    for arm,x,c in [('P1',0,P1),('P4',1,P4)]:
        vals=q[q.arm==arm][col].to_numpy(float); ax.errorbar(x,vals.mean(),yerr=vals.std(ddof=0),fmt='D',ms=5,mfc='white',mec=c,ecolor=c,capsize=2.5,lw=.9,zorder=5)
    ax.set_xlim(-.28,1.28); ax.set_xticks([0,1]); ax.set_xticklabels(['P1\n6 pairs','P4\n10 pairs']); ax.set_ylabel(ylabel); ax.set_title(title,loc='left',fontweight='bold'); clean(ax)

def family_panel(ax,df,col,title,ylabel):
    q=pairify(df); bases=np.arange(3)*1.55
    for fi,fam in enumerate(FAMILY_ORDER):
        qq=q[q.family=='F_'+fam]; pv=qq.pivot(index='pair_key',columns='arm',values=col); matched=pv.dropna()
        for j,(_,r) in enumerate(matched.iterrows()):
            jit=(j-(len(matched)-1)/2)*.008; x1=bases[fi]-.18+jit; x4=bases[fi]+.18+jit
            ax.plot([x1,x4],[r.P1,r.P4],color=LINK,lw=.8,zorder=1); ax.scatter(x1,r.P1,s=20,color=P1,zorder=3); ax.scatter(x4,r.P4,s=20,color=P4,zorder=3)
        extra=pv[pv.P1.isna() & pv.P4.notna()].P4 if set(['P1','P4']).issubset(pv.columns) else []
        for j,v in enumerate(extra): ax.scatter(bases[fi]+.18+(j-(len(extra)-1)/2)*.012,v,s=24,facecolor='white',edgecolor=P4,lw=1,zorder=3)
        for arm,x,c in [('P1',bases[fi]-.18,P1),('P4',bases[fi]+.18,P4)]:
            vals=qq[qq.arm==arm][col].to_numpy(float); ax.errorbar(x,vals.mean(),yerr=vals.std(ddof=0),fmt='D',ms=5,mfc='white',mec=c,ecolor=c,capsize=2.5,lw=.9,zorder=5)
    ax.set_xticks(bases); ax.set_xticklabels([FAMILY_LABEL[f] for f in FAMILY_ORDER]); ax.set_ylabel(ylabel); ax.set_title(title,loc='left',fontweight='bold'); clean(ax)

def fig3(d,o):
    pa=pd.read_csv(d/'G2_FACTOR_PRINCIPAL_ANGLES.csv'); prj=pd.read_csv(d/'G2_PROJECTION_DISTANCE.csv'); pro=pd.read_csv(d/'G2_PROCRUSTES_COMPARISON.csv'); loc=pd.read_csv(d/'G2_LOCATION_PAIRWISE.csv')
    fig,axs=plt.subplots(2,2,figsize=(10.8,7.4))
    family_panel(axs[0,0],pa,'sin_F','A  Principal-angle subspace discrepancy','Principal-angle sine norm')
    family_panel(axs[0,1],prj,'projection_frobenius_distance','B  Projection-matrix discrepancy','Projection Frobenius distance')
    family_panel(axs[1,0],pro,'raw_procrustes_distance','C  Procrustes discrepancy','Raw Procrustes distance')
    pair_panel(axs[1,1],loc,'location_l2_distance','D  Guide-location reproducibility','Location L2 distance')
    fig.tight_layout(); fig.savefig(o/'Fig3_reproduced.png',dpi=300); plt.close(fig)

def fig4(d,o):
    sig=pd.read_csv(d/'G3_SIGMA_PAIRWISE.csv'); cor=pd.read_csv(d/'G3_CORRELATION_CONSEQUENCE.csv'); sp=pd.read_csv(d/'G3_SPECTRUM_CONSEQUENCE.csv')
    fig,axs=plt.subplots(2,2,figsize=(10.8,7.4))
    pair_panel(axs[0,0],sig,'relative_frobenius','A  Covariance Frobenius distance','Relative Frobenius distance')
    pair_panel(axs[0,1],cor,'correlation_rmse','B  Correlation RMSE','Correlation RMSE')
    pair_panel(axs[1,0],sp,'spectrum_l2_distance','C  Spectrum discrepancy','Spectrum L2 distance')
    ax=axs[1,1]; rows=[('Covariance distance',sig,'relative_frobenius'),('Correlation RMSE',cor,'correlation_rmse'),('Spectrum L2',sp,'spectrum_l2_distance')]
    ys=np.arange(3)[::-1]
    for y,(lab,df,col) in zip(ys,rows):
        q=pairify(df); p1=q[q.arm=='P1'][col].to_numpy(float); p4all=q[q.arm=='P4'][col].to_numpy(float); p4m=q[(q.arm=='P4') & (q.seed_i!=4004) & (q.seed_j!=4004)][col].to_numpy(float)
        rall=p4all.mean()/p1.mean(); rmatch=p4m.mean()/p1.mean()
        ax.plot([rmatch,rall],[y,y],color=LINK,lw=1.2); ax.scatter(rmatch,y,s=44,color=TEAL); ax.scatter(rall,y,s=44,facecolor='white',edgecolor=TEAL,lw=1.2)
        ax.text(max(rmatch,rall)+.015,y,f'{rall:.2f}×',va='center',fontsize=8)
    ax.axvline(1,ls='--',lw=.9,color='#7C8794'); ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows]); ax.set_xlabel('P4 / P1 mean ratio'); ax.set_title('D  Descriptive P4/P1 ratios',loc='left',fontweight='bold'); clean(ax,'x')
    fig.tight_layout(); fig.savefig(o/'Fig4_reproduced.png',dpi=300); plt.close(fig)

def fig5(d,o):
    comp=pd.read_csv(d/'G3_COMPONENT_PAIRWISE.csv'); blk=pd.read_csv(d/'G3_BLOCK_COVARIANCE_CONSEQUENCE.csv'); decomp=pd.read_csv(d/'G3_COMPONENT_DECOMPOSITION_EXTENDED.csv')
    rng=np.random.default_rng(20260825); fig,axs=plt.subplots(2,2,figsize=(12,8))
    ax=axs[0,0]
    for ci,(lab,col) in enumerate([('B (AR1/diagonal)','relative_frobenius_B'),('L (low-rank)','relative_frobenius_L')]):
        for arm,c in [('P1',P1),('P4',P4)]:
            vals=comp[comp.arm==arm][col].to_numpy(float); x=ci+(-.12 if arm=='P1' else .12); ax.scatter(np.full(len(vals),x)+rng.normal(0,.02,len(vals)),vals,s=20,color=c,alpha=.8); ax.errorbar(x,vals.mean(),yerr=vals.std(ddof=0),fmt='D',mfc='white',mec=c,ecolor=c,capsize=4)
    ax.set_xticks([0,1]); ax.set_xticklabels(['B (AR1/diagonal)','L (low-rank)']); ax.set_ylabel('Relative Frobenius distance'); ax.set_title('A  Covariance-component reproducibility',loc='left',fontweight='bold'); clean(ax)
    ax=axs[0,1]
    for bi,b in enumerate(BLOCK_ORDER):
        for arm,c in [('P1',P1),('P4',P4)]:
            vals=blk[(blk.block==b)&(blk.arm==arm)].block_relative_frobenius.to_numpy(float); yy=bi+(-.08 if arm=='P1' else .08); ax.errorbar(vals.mean(),yy,xerr=vals.std(ddof=0),fmt='o',color=c,capsize=4)
    ax.set_yticks(range(4)); ax.set_yticklabels([BLOCK_LABEL[b] for b in BLOCK_ORDER]); ax.invert_yaxis(); ax.set_xlabel('Relative Frobenius distance (mean ± SD)'); ax.set_title('B  Semantic-block covariance consequence',loc='left',fontweight='bold'); clean(ax,'x')
    ax=axs[1,0]; xx=np.arange(len(decomp)); labs=[f"{a}-{int(seed)}" for a,seed in zip(decomp.arm,decomp.seed)]; ax.bar(xx,decomp.ar1_trace_share,label='AR1/diagonal-temporal'); ax.bar(xx,decomp.low_rank_trace_share,bottom=decomp.ar1_trace_share,label='Low-rank'); ax.set_ylim(0,1); ax.set_xticks(xx); ax.set_xticklabels(labs,rotation=45,ha='right'); ax.set_ylabel('Share of total covariance trace'); ax.set_title('C  Per-run covariance trace composition',loc='left',fontweight='bold'); ax.legend(frameon=False,ncol=2); clean(ax)
    ax=axs[1,1]
    for arm,c in [('P1',P1),('P4',P4)]:
        q=decomp[decomp.arm==arm]; ax.scatter(q.low_rank_trace_share,q.low_rank_frobenius_ratio,s=40,color=c,label=arm)
        for _,r in q.iterrows(): ax.annotate(str(int(r.seed)),(r.low_rank_trace_share,r.low_rank_frobenius_ratio),xytext=(3,3),textcoords='offset points',fontsize=7)
    ax.set_xlabel('Low-rank trace share'); ax.set_ylabel('Low-rank Frobenius ratio'); ax.set_title('D  Low-rank magnitude across evaluable V3 fits',loc='left',fontweight='bold'); ax.legend(frameon=False); clean(ax,'both')
    fig.tight_layout(); fig.savefig(o/'Fig5_reproduced.png',dpi=300); plt.close(fig)

def s1(d,o):
    x=pd.read_csv(d/'V8_P4_VS_P1_DESCRIPTIVE_COMPARISON.csv'); specs=[('A  Location L2','location_l2'),('B  Raw Procrustes','raw_procrustes'),('C  Covariance distance','covariance_reproducibility_distance'),('D  Correlation RMSE','correlation_rmse'),('E  Spectrum L2','normalized_spectrum_l2'),('F  Low-rank distance','L_relative_frobenius')]
    fig,axs=plt.subplots(2,3,figsize=(12,7.2)); ys=np.arange(len(x))[::-1]
    for ax,(title,base) in zip(axs.flat,specs):
        for yi,(_,r) in zip(ys,x.iterrows()):
            v1=float(r[base+'_P1']); v4=float(r[base+'_P4']); red=float(r[base+'_P4_reduction_pct']); ax.plot([v1,v4],[yi,yi],color=LINK,lw=1.2); ax.scatter(v1,yi,s=25,color=P1); ax.scatter(v4,yi,s=25,color=P4); ax.text(v1,yi,f' P4 {"↓" if red>=0 else "↑"}{abs(red):.1f}%',va='center',fontsize=7.5)
        ax.set_yticks(ys); ax.set_yticklabels(x.target_id); ax.set_title(title,loc='left',fontweight='bold'); ax.set_xlabel('Metric value'); clean(ax,'x')
    fig.tight_layout(); fig.savefig(o/'S1_Fig_reproduced.png',dpi=300); plt.close(fig)

def s2(d,o):
    x=pd.read_csv(d/'V8_CELL_SUMMARY_FLAT.csv'); metrics=[('Location L2','location_l2'),('Raw Procrustes','raw_procrustes'),('Covariance distance','covariance_reproducibility_distance'),('Correlation RMSE','correlation_rmse'),('Spectrum L2','normalized_spectrum_l2'),('Low-rank distance','L_relative_frobenius')]
    wide=[]
    for label,col in metrics:
        p=x.pivot(index='target_id',columns='particles',values=col); wide.append(np.log(p.P4/p.P1).rename(label))
    W=pd.concat(wide,axis=1); C=W.corr(method='spearman')
    fig,ax=plt.subplots(figsize=(7.2,6.4)); mask=np.triu(np.ones_like(C,dtype=bool),k=1); arr=C.to_numpy().copy(); arr[mask]=np.nan
    im=ax.imshow(arr,cmap='coolwarm',vmin=-1,vmax=1)
    ax.set_xticks(range(len(C))); ax.set_xticklabels(C.columns,rotation=30,ha='right'); ax.set_yticks(range(len(C))); ax.set_yticklabels(C.index)
    for i in range(len(C)):
        for j in range(len(C)):
            if not mask[i,j]: ax.text(j,i,f'{C.iloc[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(C.iloc[i,j])>.65 else INK)
    cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04); cb.set_label('Spearman ρ'); fig.tight_layout(); fig.savefig(o/'S2_Fig_reproduced.png',dpi=300); plt.close(fig)

def s3(d,o):
    x=pd.read_csv(d/'V8_P4_VS_P1_DESCRIPTIVE_COMPARISON.csv'); specs=[('Representation','Location L2','location_l2_P4_reduction_pct'),('Representation','Principal-angle sine norm','principal_angle_sin_F_P4_reduction_pct'),('Representation','Projection Frobenius','projection_frobenius_P4_reduction_pct'),('Representation','Raw Procrustes','raw_procrustes_P4_reduction_pct'),('Covariance consequence','Covariance distance','covariance_reproducibility_distance_P4_reduction_pct'),('Covariance consequence','Correlation RMSE','correlation_rmse_P4_reduction_pct'),('Covariance consequence','Spectrum L2','normalized_spectrum_l2_P4_reduction_pct'),('Covariance consequence','Low-rank distance','L_relative_frobenius_P4_reduction_pct'),('Known-target accuracy','Mean error to target','mean_error_to_target_P4_reduction_pct'),('Known-target accuracy','Covariance error to target','covariance_error_to_target_P4_reduction_pct'),('Known-target accuracy','Analytic Gaussian KL','analytic_gaussian_kl_P4_reduction_pct')]
    fig,ax=plt.subplots(figsize=(8.5,8.5)); y=len(specs)-1; ys=[]; labs=[]
    for grp,lab,col in specs:
        vals=x[col].astype(float).to_numpy(); ax.hlines(y,vals.min(),vals.max(),color=LINK,lw=1.5); ax.scatter(vals.mean(),y,s=45,color=P4); ax.text(vals.mean()+1.0,y,f'{vals.mean():.1f}',va='center',fontsize=8); ys.append(y); labs.append(lab); y-=1
    ax.axvline(0,color='#7C8794',ls='--',lw=.9); ax.set_yticks(ys); ax.set_yticklabels(labs); ax.set_xlabel('P4 reduction relative to P1 (%)'); clean(ax,'x'); fig.tight_layout(); fig.savefig(o/'S3_Fig_reproduced.png',dpi=300); plt.close(fig)

def s4(d,o):
    alloc=pd.read_csv(d/'G2_BLOCK_COVARIANCE_ALLOCATION.csv'); delta=pd.read_csv(d/'G2_BLOCK_DELTA_ALLOCATION.csv'); rho=pd.read_csv(d/'G2_DELTA_RHO_PAIRWISE.csv'); cv=pd.read_csv(d/'G2_MARGINAL_SCALE_CV.csv')
    fig,axs=plt.subplots(2,2,figsize=(11,7.2)); bases=np.arange(4)
    ax=axs[0,0]
    for i,b in enumerate(BLOCK_ORDER):
        for arm,xoff,c in [('P1',-.06,P1),('P4',.06,P4)]:
            vals=alloc[(alloc.semantic_block==b)&(alloc.arm==arm)].block_fraction_of_total_variance.to_numpy(float); ax.scatter(np.full(len(vals),i+xoff),vals,s=18,color=c,alpha=.75); ax.errorbar(i+xoff,vals.mean(),yerr=vals.std(ddof=0),fmt='D',mfc='white',mec=c,ecolor=c,capsize=2)
    ax.set_xticks(bases); ax.set_xticklabels([BLOCK_LABEL[b] for b in BLOCK_ORDER]); ax.set_ylabel('Fraction of total variance'); ax.set_title('A  Semantic-block variance allocation',loc='left',fontweight='bold'); clean(ax)
    ax=axs[0,1]
    for i,b in enumerate(BLOCK_ORDER):
        q=pairify(delta[delta.semantic_block==b]); pv=q.pivot(index='pair_key',columns='arm',values='delta_allocation'); matched=pv.dropna()
        for _,r in matched.iterrows(): ax.plot([i-.10,i+.10],[r.P1,r.P4],color=LINK,lw=.7)
        for arm,xoff,c in [('P1',-.10,P1),('P4',.10,P4)]:
            vals=q[q.arm==arm].delta_allocation.to_numpy(float); ax.scatter(np.full(len(vals),i+xoff),vals,s=16,color=c,alpha=.75); ax.errorbar(i+xoff,vals.mean(),yerr=vals.std(ddof=0),fmt='D',mfc='white',mec=c,ecolor=c,capsize=2)
    ax.axhline(0,color='#888',ls='--',lw=.7); ax.set_xticks(bases); ax.set_xticklabels([BLOCK_LABEL[b] for b in BLOCK_ORDER]); ax.set_ylabel('Delta allocation'); ax.set_title('B  Pairwise block-allocation difference',loc='left',fontweight='bold'); clean(ax)
    pair_panel(axs[1,0],rho,'delta_rho','C  Pairwise AR1-parameter difference','Pairwise |Δρ|')
    ax=axs[1,1]
    for i,b in enumerate(BLOCK_ORDER):
        for arm,xoff,c in [('P1',-.08,P1),('P4',.08,P4)]:
            vals=cv[(cv.semantic_block==b)&(cv.arm==arm)].cv_marginal_sd.to_numpy(float); # deterministic spread by rank
            offs=np.linspace(-.035,.035,len(vals)); ax.scatter(np.full(len(vals),i+xoff)+offs,vals,s=6,color=c,alpha=.22); ax.errorbar(i+xoff,vals.mean(),yerr=vals.std(ddof=0),fmt='D',mfc='white',mec=c,ecolor=c,capsize=2)
    ax.set_xticks(bases); ax.set_xticklabels([BLOCK_LABEL[b] for b in BLOCK_ORDER]); ax.set_ylabel('Marginal-scale CV'); ax.set_title('D  Latent-scale variability by semantic block',loc='left',fontweight='bold'); clean(ax)
    fig.tight_layout(); fig.savefig(o/'S4_Fig_reproduced.png',dpi=300); plt.close(fig)

def s5(d,o):
    x=pd.read_csv(d/'V8_CELL_SUMMARY_FLAT.csv'); groups={'Representation':[('Location L2','location_l2'),('Principal-angle sine','principal_angle_sin_F'),('Projection Frobenius','projection_frobenius'),('Raw Procrustes','raw_procrustes')],'Covariance consequence':[('Covariance distance','covariance_reproducibility_distance'),('Correlation RMSE','correlation_rmse'),('Spectrum L2','normalized_spectrum_l2'),('B relative Frobenius','B_relative_frobenius'),('L relative Frobenius','L_relative_frobenius')],'Known-target accuracy':[('Mean error to target','mean_error_to_target'),('Covariance error to target','covariance_error_to_target'),('Analytic Gaussian KL','analytic_gaussian_kl')]}; colors={'Representation':PURPLE,'Covariance consequence':TEAL,'Known-target accuracy':GOLD}
    rows=[]
    for grp,metrics in groups.items():
        for lab,col in metrics:
            p=x.pivot(index='target_id',columns='particles',values=col); rows.append((grp,lab,(p.P4/p.P1).dropna()))
    fig,ax=plt.subplots(figsize=(8.5,8.5)); y=len(rows)-1; ys=[]; labs=[]
    for grp,lab,s in rows:
        vals=s.to_numpy(float); ax.hlines(y,vals.min(),vals.max(),color=LINK,lw=1.5); ax.scatter(vals,np.full(vals.shape,y,dtype=float),s=16,color=colors[grp],alpha=.35); ax.scatter(vals.mean(),y,s=45,color=colors[grp]); ax.text(max(vals.max(),vals.mean())+.02,y,f'{vals.mean():.2f}×',va='center',fontsize=8); ys.append(y); labs.append(lab); y-=1
    ax.axvline(1,color='#7C8794',ls='--',lw=.9); ax.set_yticks(ys); ax.set_yticklabels(labs); ax.set_xlabel('P4 / P1 ratio across 6 targets'); clean(ax,'x'); fig.tight_layout(); fig.savefig(o/'S5_Fig_reproduced.png',dpi=300); plt.close(fig)

def s6(d,o):
    x=pd.read_csv(d/'V8_CELL_SUMMARY_FLAT.csv'); piv=x.pivot(index='target_id',columns='particles'); targets=list(piv.index); specs=[('A  Fitted low-rank trace share','fitted_low_rank_trace_share','Trace share'),('B  B relative Frobenius','B_relative_frobenius','B relative Frobenius'),('C  L relative Frobenius','L_relative_frobenius','L relative Frobenius')]
    fig,axs=plt.subplots(1,3,figsize=(11,4.0)); yy=np.arange(len(targets))[::-1]
    for ax,(title,col,xlab) in zip(axs,specs):
        p1=piv[(col,'P1')]; p4=piv[(col,'P4')]
        for y,t in zip(yy,targets):
            v1=float(p1.loc[t]); v4=float(p4.loc[t]); ax.plot([v1,v4],[y,y],color=LINK,lw=1); ax.scatter(v1,y,s=24,color=P1); ax.scatter(v4,y,s=24,color=P4)
        ax.set_yticks(yy); ax.set_yticklabels(targets); ax.set_xlabel(xlab); ax.set_title(title,loc='left',fontweight='bold'); clean(ax,'x')
    fig.tight_layout(); fig.savefig(o/'S6_Fig_reproduced.png',dpi=300); plt.close(fig)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data-dir',default='.'); ap.add_argument('--output-dir',default='reproduced_outputs'); a=ap.parse_args(); d=Path(a.data_dir).resolve(); o=Path(a.output_dir).resolve(); o.mkdir(parents=True,exist_ok=True)
    pd.read_csv(d/'V3_RUN_MATRIX.csv').to_csv(o/'Table1_run_matrix.csv',index=False); t2=table2(d); t2.to_csv(o/'Table2_stage_I.csv',index=False); t3=table3(d); t3.to_csv(o/'Table3_stage_II.csv',index=False); pd.read_csv(d/'G3_COMPONENT_DECOMPOSITION_EXTENDED.csv').to_csv(o/'Table4_component_spectrum.csv',index=False); t5=table5(d); t5.to_csv(o/'Table5_simulation_summary.csv',index=False)
    fig3(d,o); fig4(d,o); fig5(d,o); s1(d,o); s2(d,o); s3(d,o); s4(d,o); s5(d,o); s6(d,o)
    checks=[(float(t2.loc[t2.metric=='location L2','P1'].iloc[0]),2.04564,5e-5,'Table2 P1 location L2'),(float(t2.loc[t2.metric=='marginal-scale CV','P4'].iloc[0]),0.0323272,5e-7,'Table2 P4 marginal CV'),(float(t3.loc[t3.metric=='Sigma relative Frobenius','P1 mean'].iloc[0]),0.280254,5e-7,'Table3 P1 covariance'),(float(t3.loc[t3.metric=='correlation RMSE','P4 mean'].iloc[0]),0.00780648,5e-8,'Table3 P4 corr RMSE')]
    for got,exp,tol,label in checks:
        if abs(got-exp)>tol: raise RuntimeError(f'Frozen value check failed: {label}: {got} vs {exp}')
    expected=['Table1_run_matrix.csv','Table2_stage_I.csv','Table3_stage_II.csv','Table4_component_spectrum.csv','Table5_simulation_summary.csv','Fig3_reproduced.png','Fig4_reproduced.png','Fig5_reproduced.png','S1_Fig_reproduced.png','S2_Fig_reproduced.png','S3_Fig_reproduced.png','S4_Fig_reproduced.png','S5_Fig_reproduced.png','S6_Fig_reproduced.png']
    missing=[f for f in expected if not (o/f).exists()]
    if missing: raise RuntimeError('Missing outputs: '+str(missing))
    print('REPRODUCTION CHECK: PASS')
    print('OUTPUT COUNT:',len(expected))
if __name__=='__main__': main()
