#!/usr/bin/env python3
"""Display-only assembly of the RetinaPlex-State results into Figure 1-5."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'figures'/'final'; OUT.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
 'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,
 'legend.fontsize':6.5,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.right':False,
 'axes.spines.top':False,'axes.linewidth':.65,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
TRAITS=['AMD','POAG','RE']; TC={'AMD':'#5B7FA3','POAG':'#B56576','RE':'#C68A3A'}
CL=['Muller glia','RGC','amacrine cells','astrocytes','bipolar cells','cones','horizontal cells','microglia','rods']
CLL={'Muller glia':'Müller glia','RGC':'RGC','amacrine cells':'Amacrine','astrocytes':'Astrocytes','bipolar cells':'Bipolar','cones':'Cones','horizontal cells':'Horizontal','microglia':'Microglia','rods':'Rods'}
AL={'MG':'Müller','AC':'Amacrine','BC':'Bipolar','Cone':'Cone','Rod':'Rod','HC':'Horizontal','RGC':'RGC','Astrocyte':'Astrocyte','Microglia':'Microglia'}
def rd(x): return pd.read_csv(ROOT/x,sep='\t')
SUPERSCRIPT = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')
def format_q(q):
    q=float(q)
    if q < 0.01 and q > 0:
        exp=int(math.floor(math.log10(q)))
        mant=q/(10**exp)
        return f'{mant:.2f} × $10^{{{exp}}}$'
    return f'{q:.2f}'
def save(fig,n):
    fig.savefig(OUT/f'Figure{n}.svg',bbox_inches='tight')
    fig.savefig(OUT/f'Figure{n}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'Figure{n}.png',dpi=600,bbox_inches='tight')
    fig.savefig(OUT/f'Figure{n}.tiff',dpi=600,bbox_inches='tight')
    plt.close(fig)
def pl(ax,s): ax.text(-.12,1.04,s,transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
def grid(ax): ax.grid(axis='y',color='#E7E7E7',lw=.45,zorder=0); ax.tick_params(length=2.5,width=.55)
def box(ax,xy,wh,s,face='#EEF3F7',edge='#5B7FA3',fs=7):
    x,y=xy; w,h=wh; ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.02',transform=ax.transAxes,facecolor=face,edgecolor=edge,lw=.8)); ax.text(x+w/2,y+h/2,s,transform=ax.transAxes,ha='center',va='center',fontsize=fs)
def arrow(ax,p,q,c='#7A858F'): ax.add_patch(FancyArrowPatch(p,q,transform=ax.transAxes,arrowstyle='-|>',mutation_scale=8,lw=.75,color=c))

def fig1():
    rg=rd('results/phase1/pairwise_ldsc_rg.tsv'); re_id='RE_2026_EUR_COMPOSITE'
    rows=[rg[(rg.trait_1=='AMD_GCST003219')&(rg.trait_2=='POAG_GCST90011766')].iloc[0],rg[(rg.trait_1=='AMD_GCST003219')&(rg.trait_2==re_id)].iloc[0],rg[(rg.trait_1=='POAG_GCST90011766')&(rg.trait_2==re_id)].iloc[0]]
    fig=plt.figure(figsize=(7.2,4.75)); gs=fig.add_gridspec(2,2,height_ratios=[1.1,1],hspace=.48,wspace=.34)
    a=fig.add_subplot(gs[0,:]); a.axis('off'); pl(a,'a'); a.text(0,.98,'Three European ocular GWAS',transform=a.transAxes,va='top',fontsize=8,fontweight='bold')
    for x,s,f,e in [(.01,'AMD\nGCST003219','#EAF0F6','#5B7FA3'),(.23,'POAG\nGCST90011766','#F6E9EC','#B56576'),(.45,'European composite\nRE liability','#F8F0E4','#C68A3A')]: box(a,(x,.30),(.18,.38),s,f,e)
    box(a,(.72,.57),(.24,.22),'signed genome-wide\narchitecture','#F1F1F1','#68737D'); box(a,(.72,.23),(.24,.22),'retinal RNA-defined\ncellular context','#E9F3EE','#4F8A6F')
    for x in (.19,.41,.63): arrow(a,(x,.49),(.71,.67),'#68737D'); arrow(a,(x,.49),(.71,.34),'#4F8A6F')
    a.text(.5,.04,'HDL-L local resolution  →  donor-aware scDRS/MAGMA  →  Müller state and HRCA OCR boundaries',transform=a.transAxes,ha='center',fontsize=6.7,color='#4E5963')
    b=fig.add_subplot(gs[1,0]); pl(b,'b'); y=np.arange(3); v=np.array([r.rg for r in rows]); se=np.array([r.rg_SE for r in rows]); b.errorbar(v,y,xerr=1.96*se,fmt='o',color='#3F6D8C',ecolor='#3F6D8C',lw=1.1,capsize=2.5,ms=4.8,zorder=3); b.axvline(0,color='#9AA2A8',lw=.6); b.set_yticks(y,['AMD–POAG','AMD–RE','POAG–RE']); b.set_xlabel('Signed genome-wide genetic correlation (rᵍ)'); b.set_title('Partial and oppositely signed relationships',loc='left',pad=6); b.set_xlim(-.38,.20)
    for yy,r in zip(y,rows): b.text(.99,yy,f'FDR {r.FDR:.2g}',transform=b.get_yaxis_transform(),ha='right',va='center',fontsize=6)
    b.text(.02,1.01,'POAG–RE no-UKBB: rg = -0.1766, P = 3.20e-9',transform=b.transAxes,fontsize=5.8,color='#4E5963',va='bottom'); grid(b)
    c=fig.add_subplot(gs[1,1]); pl(c,'c'); pairs=['AMD–POAG','AMD–RE','POAG–RE']; interpretable=[107,8,6]; hotspots=[0,0,0]; x=np.arange(3); width=.34; c.bar(x-width/2,interpretable,width,color='#B8CADB',label='Interpretable regions'); c.bar(x+width/2,hotspots,width,color='#E3E7EA',label='FDR-supported hotspots'); c.set_xticks(x,pairs); c.set_ylabel('Count'); c.set_title('Local architecture at tested resolution',loc='left',pad=6); c.set_ylim(0,120); c.legend(frameon=False,fontsize=5.1,loc='upper right',handlelength=1.1); [c.text(xx-width/2,v+3,str(v),ha='center',fontsize=6.5,fontweight='bold') for xx,v in zip(x,interpretable)]; [c.text(xx+width/2,4,'0',ha='center',fontsize=6.5,fontweight='bold',color='#69747D') for xx in x]; grid(c); save(fig,1)

def fig2():
    d=rd('results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv'); m=rd('results/phase3a1/magma_broad_celltype_enrichment.tsv'); d['broad_class']=pd.Categorical(d.broad_class,categories=CL,ordered=True); mat=d.pivot(index='trait',columns='broad_class',values='mean_donor_delta').reindex(index=TRAITS,columns=CL); f=d.pivot(index='trait',columns='broad_class',values='global_broad_FDR').reindex(index=TRAITS,columns=CL)
    fig=plt.figure(figsize=(7.8,5.25)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.35],hspace=.72,wspace=.42); a=fig.add_subplot(gs[0,:]); pl(a,'a'); im=a.imshow(mat.values,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2.5)); a.set_xticks(range(9),[CLL[x] for x in CL]); [lab.set_rotation(35) or lab.set_rotation_mode('anchor') or lab.set_fontsize(6.2) for lab in a.get_xticklabels()]; a.set_yticks(range(3),TRAITS); a.set_title('Donor-aware broad retinal cellular effects',loc='left',pad=12)
    for i in range(3):
        for j in range(9):
            val=mat.iloc[i,j]; q=f.iloc[i,j]; a.text(j,i,f'{val:+.2f}\nFDR {format_q(q)}',ha='center',va='center',fontsize=5.9,color='white' if abs(val)>1.1 else '#1E2730')
    cb=fig.colorbar(im,ax=a,fraction=.02,pad=.015); cb.set_label('Mean within-donor effect',fontsize=6.5); cb.ax.tick_params(labelsize=6,length=2)
    donor=rd('results/phase3a1/donor_class_deltas.tsv'); mg=donor[donor.broad_class=='Muller glia']; b=fig.add_subplot(gs[1,0]); pl(b,'b')
    for i,t in enumerate(TRAITS):
        vv=mg.loc[mg.trait==t,'delta'].to_numpy(); b.scatter(np.full(len(vv),i)+np.linspace(-.13,.13,len(vv)),vv,s=19,color=TC[t],zorder=3,edgecolor='white',lw=.35); r=d[(d.trait==t)&(d.broad_class=='Muller glia')].iloc[0]; b.errorbar(i,r.mean_donor_delta,yerr=[[r.mean_donor_delta-r.CI95_lower],[r.CI95_upper-r.mean_donor_delta]],fmt='D',color='#252C32',ms=4,capsize=3,lw=1,zorder=4)
    b.axhline(0,color='#9AA2A8',lw=.6); b.set_xticks(range(3),TRAITS); b.set_ylabel('Within-donor effect'); b.set_title('Müller glia: n = 8 donors, all positive',loc='left',pad=6); b.set_ylim(-.2,2.4); b.text(.02,.03,'n = 8 donors\n8/8 positive for each trait\nDiamonds: mean and 95% CI',transform=b.transAxes,fontsize=6.2,color='#4E5963'); grid(b)
    c=fig.add_subplot(gs[1,1]); pl(c,'c'); ac=donor[(donor.trait=='RE')&(donor.broad_class.isin(['Muller glia','amacrine cells']))]
    for i,cl in enumerate(['Muller glia','amacrine cells']):
        vv=ac.loc[ac.broad_class==cl,'delta'].to_numpy(); c.scatter(np.full(len(vv),i)+np.linspace(-.12,.12,len(vv)),vv,s=20,color='#4F8A6F' if cl=='Muller glia' else '#C68A3A',edgecolor='white',lw=.35,zorder=3); r=d[(d.trait=='RE')&(d.broad_class==cl)].iloc[0]; c.errorbar(i,r.mean_donor_delta,yerr=[[r.mean_donor_delta-r.CI95_lower],[r.CI95_upper-r.mean_donor_delta]],fmt='D',color='#252C32',ms=4,capsize=3,lw=1,zorder=4)
    c.axhline(0,color='#9AA2A8',lw=.6); c.set_xticks([0,1],['Müller','Amacrine']); c.set_ylabel('Within-donor effect'); c.set_title('RE adds an amacrine component',loc='left',pad=6); c.set_ylim(.8,2.2); c.text(.02,.03,'n = 8 donors\nRE amacrine: +1.2530\nFDR 2.88 × $10^{-6}$; MAGMA 0.0261',transform=c.transAxes,fontsize=5.8,color='#4E5963'); grid(c); save(fig,2)

def fig3():
    rg=rd('results/phase1/pairwise_ldsc_rg.tsv'); p=rd('results/phase3a1/cross_trait_broad_profile_similarity.tsv'); st=rd('results/phase3a1/cross_trait_standardized_effect_profiles.tsv'); pairs=['AMD-POAG','AMD-RE','POAG-RE']; key={'AMD-POAG':('AMD_GCST003219','POAG_GCST90011766'),'AMD-RE':('AMD_GCST003219','RE_2026_EUR_COMPOSITE'),'POAG-RE':('POAG_GCST90011766','RE_2026_EUR_COMPOSITE')}; xx=[]; yy=[]
    for z in pairs:
        k=key[z]; r=rg[(rg.trait_1==k[0])&(rg.trait_2==k[1])].iloc[0]; q=p[(p.trait_a==z.split('-')[0])&(p.trait_b==z.split('-')[1])].iloc[0]; xx.append(r.rg); yy.append(q.spearman_rho)
    fig=plt.figure(figsize=(7.2,4.7)); gs=fig.add_gridspec(1,3,width_ratios=[1,1.35,.95],wspace=.52); a=fig.add_subplot(gs[0,0]); pl(a,'a')
    for z,x,y in zip(['AMD–POAG','AMD–RE','POAG–RE'],xx,yy): a.scatter(x,y,s=40 if z=='POAG–RE' else 28,color='#B56576' if z=='POAG–RE' else '#5B7FA3',zorder=3); a.text(x+.012,y+.015,z,fontsize=6,fontweight='bold' if z=='POAG–RE' else 'normal')
    a.axvline(0,color='#A5ADB3',lw=.55); a.axhline(.5,color='#D8DDE1',lw=.55,ls='--'); a.set_xlabel('Signed genome-wide rᵍ'); a.set_ylabel('Cellular-context similarity ρ'); a.set_title('Distinct cross-trait dimensions',loc='left',pad=6); a.set_xlim(-.22,.10); a.set_ylim(.25,.98); a.text(-.21,.315,'Bootstrap 95% CI\nAMD–POAG 0.5476–0.8167\nAMD–RE 0.1905–0.6833\nPOAG–RE 0.5714–0.9524',fontsize=4.8,color='#4E5963',va='top'); grid(a)
    b=fig.add_subplot(gs[0,1]); pl(b,'b'); q=st.set_index('broad_class').loc[CL]; ix=np.arange(9); b.plot(ix,q.POAG,color=TC['POAG'],marker='o',ms=3.5,lw=1.2,label='POAG'); b.plot(ix,q.RE,color=TC['RE'],marker='o',ms=3.5,lw=1.2,label='RE'); b.axhline(0,color='#A5ADB3',lw=.55); b.set_xticks(ix,[CLL[x] for x in CL]); [lab.set_rotation(42) or lab.set_rotation_mode('anchor') for lab in b.get_xticklabels()]; b.set_ylabel('Standardized class effect'); b.set_title('POAG–RE broad-cell profile',loc='left',pad=6); b.legend(loc='upper right',frameon=False,handlelength=1.5); b.text(.02,.03,'rᵍ = −0.1732\nρ = +0.8833\nbootstrap 95% CI\n0.5714–0.9524',transform=b.transAxes,fontsize=6.1,color='#4E5963'); grid(b)
    c=fig.add_subplot(gs[0,2]); pl(c,'c'); c.axis('off'); c.text(.5,.97,'Interpretive distinction',ha='center',va='top',fontsize=8,fontweight='bold'); box(c,(.08,.59),(.84,.22),'Signed genome-wide\nallelic architecture','#F6E9EC','#B56576'); box(c,(.08,.20),(.84,.22),'Unsigned / profile-level\ncellular context','#E9F3EE','#4F8A6F'); arrow(c,(.5,.57),(.5,.44)); c.text(.5,.49,'can diverge',ha='center',va='center',fontsize=6.5,color='#4E5963'); c.text(.5,.06,'Signed allelic direction ≠\ncellular-context localization',ha='center',va='bottom',fontsize=6.3,color='#4E5963'); save(fig,3)

def fig4():
    sm=rd('results/phase3a1/mapping_selection_sensitivity.tsv'); lo=rd('results/phase3a1/leave_one_donor_out.tsv'); d=rd('results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv'); m=rd('results/phase3a1/magma_broad_celltype_enrichment.tsv'); fig=plt.figure(figsize=(8.0,4.9)); gs=fig.add_gridspec(1,3,width_ratios=[1.1,1.1,1.35],wspace=.55); a=fig.add_subplot(gs[0,0]); pl(a,'a'); rules=['A_TOP_75_PERCENT_MAPPING','B_EXCLUDE_LOWEST_MAPPING','C_MEDIAN_DONOR_CLASS_SCORE']; labels=['Top 75%\nmapping','Exclude lowest\ndonor','Median donor\nscore']
    for i,t in enumerate(TRAITS):
        z=sm[(sm.trait==t)&(sm.broad_class=='Muller glia')]
        for j,rn in enumerate(rules): a.scatter(j+(i-1)*.13,z[z.sensitivity==rn].effect.iloc[0],color=TC[t],s=20,zorder=3)
    a.axhline(0,color='#9AA2A8',lw=.55); a.set_xticks(range(3),labels); a.set_xticklabels(labels,fontsize=6,rotation=15,ha='right',rotation_mode='anchor'); a.set_ylabel('Müller effect'); a.set_title('Mapping-selection sensitivity',loc='left',pad=6); a.legend(handles=[plt.Line2D([0],[0],marker='o',ls='',color=c,label=t,ms=4) for t,c in TC.items()],frameon=False,loc='upper left',fontsize=6); a.set_ylim(.7,2.3); grid(a)
    b=fig.add_subplot(gs[0,1]); pl(b,'b')
    for i,t in enumerate(TRAITS):
        z_full=lo[(lo.trait==t)&(lo.broad_class=='Muller glia')].leave_one_out_mean_delta.to_numpy(dtype=float); z=z_full[np.isfinite(z_full)]; assert len(z)==8; b.vlines(i,z.min(),z.max(),color=TC[t],lw=2,alpha=.55); b.scatter(np.full(len(z),i),z,s=17,color=TC[t],edgecolor='white',lw=.35,zorder=3); pri=d[(d.trait==t)&(d.broad_class=='Muller glia')].iloc[0].mean_donor_delta; b.scatter(i,pri,marker='D',s=23,color='#20262B',zorder=4)
    b.axhline(0,color='#9AA2A8',lw=.55); b.set_xticks(range(3),TRAITS); b.set_ylabel('Leave-one-donor-out effect'); b.set_title('Donor sensitivity',loc='left',pad=6); b.text(.02,.03,'Dots: each omitted donor\nDiamonds: primary mean',transform=b.transAxes,fontsize=6,color='#4E5963'); b.set_ylim(.8,2.3); grid(b)
    c=fig.add_subplot(gs[0,2]); pl(c,'c')
    for t in TRAITS:
        for cl in CL:
            r=d[(d.trait==t)&(d.broad_class==cl)].iloc[0]; q=m[(m.trait==t)&(m.broad_class==cl)].iloc[0]; mk='D' if r.global_broad_FDR<.05 else 'o'; c.scatter(r.mean_donor_delta,-math.log10(max(float(q.global_broad_FDR),1e-12)),color=TC[t],marker=mk,s=24 if mk=='D' else 13,alpha=.85,edgecolor='white',lw=.25)
    c.axhline(-math.log10(.05),color='#9AA2A8',lw=.55,ls='--'); c.axvline(0,color='#C4CBD0',lw=.55); c.set_xlabel('scDRS donor-aware effect'); c.set_ylabel('-log10 MAGMA broad FDR'); c.set_title('Orthogonal associative evidence',loc='left',pad=6); c.text(.03,.97,'Diamonds: scDRS Tier 1\nMAGMA is not causal validation',transform=c.transAxes,fontsize=6,color='#4E5963',va='top'); c.set_ylim(-.05,4); grid(c); save(fig,4)

def fig5():
    d=rd('results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv'); s=rd('results/phase3a2/muller_state_eligibility.tsv'); at=rd('results/phase3b/sldsc_retinal_ocr_enrichment.tsv'); at=at[at.cell_class.isin(list(AL))]; fig=plt.figure(figsize=(7.2,5.35)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.1],hspace=.72,wspace=.48)
    a=fig.add_subplot(gs[0,0]); pl(a,'a'); rows=[('AMD','Muller glia'),('POAG','Muller glia'),('RE','Muller glia'),('RE','amacrine cells')]; vv=[d[(d.trait==t)&(d.broad_class==cl)].iloc[0].mean_donor_delta for t,cl in rows]; a.barh(range(4),vv,color=[TC[t] if cl=='Muller glia' else '#C68A3A' for t,cl in rows],height=.58); a.axvline(0,color='#9AA2A8',lw=.55); a.set_yticks(range(4),['AMD Müller','POAG Müller','RE Müller','RE amacrine']); a.set_xlabel('Donor-aware effect'); a.set_title('Broad RNA tier: supported',loc='left',pad=6)
    for i,(v,(t,cl)) in enumerate(zip(vv,rows)): a.text(v+.05,i,f'FDR {d[(d.trait==t)&(d.broad_class==cl)].iloc[0].global_broad_FDR:.1e}',va='center',fontsize=5.8)
    a.set_xlim(0,2.3); grid(a)
    b=fig.add_subplot(gs[0,1]); pl(b,'b'); s=s.sort_values('state'); yy=np.arange(len(s)); b.barh(yy,s.donors_with_30_cells,color='#AEB7BE',height=.62); b.axvline(5,color='#B56576',ls='--',lw=.8); b.set_yticks(yy,s.state.str.replace('MG_state_','S',regex=False),fontsize=5.5); b.set_xlabel('Donors with ≥30 cells'); b.set_title('Müller state tier: unsupported',loc='left',pad=6); b.text(.03,.03,'11 Leiden states; 0 eligible\nAll states failed prespecified replication/composition criteria',transform=b.transAxes,fontsize=5.8,color='#4E5963'); b.set_xlim(0,8.6); grid(b)
    c=fig.add_subplot(gs[1,0]); pl(c,'c'); order=list(AL); pv=at.pivot(index='cell_class',columns='trait',values='global_primary_fdr').reindex(index=order,columns=TRAITS); im=c.imshow(-np.log10(pv.astype(float).clip(lower=1e-12)),aspect='auto',cmap='Greys',vmin=0,vmax=.5); c.set_yticks(range(9),[AL[x] for x in order]); c.set_xticks(range(3),TRAITS); c.set_title('Class-level OCR tier: unsupported',loc='left',pad=6); c.set_xlabel('Trait')
    for i in range(9):
        for j in range(3): c.text(j,i,f'{pv.iloc[i,j]:.2f}',ha='center',va='center',fontsize=5.7,color='white' if -math.log10(pv.iloc[i,j])>.30 else '#20252B')
    c.text(0,-.24,'No cell class reached FDR < 0.05; OCR resource is pooled',transform=c.transAxes,fontsize=5.8,color='#4E5963')
    e=fig.add_subplot(gs[1,1]); pl(e,'d'); e.axis('off'); e.text(.5,.97,'Evidence boundary',ha='center',va='top',fontsize=8,fontweight='bold'); box(e,(.10,.68),(.80,.16),'RNA-defined broad\ncellular convergence','#E9F3EE','#4F8A6F',7.3); box(e,(.10,.43),(.80,.16),'Discrete Müller state\nnot supported','#F2F2F2','#8B959C',7.3); box(e,(.10,.18),(.80,.16),'Class-level OCR enrichment\nnot supported','#F2F2F2','#8B959C',7.3); arrow(e,(.5,.66),(.5,.60)); arrow(e,(.5,.41),(.5,.35)); e.text(.5,.04,'Claim resolution stops at broad transcriptional context',ha='center',va='bottom',fontsize=6.2,color='#4E5963'); save(fig,5)

if __name__=='__main__':
    fig1(); fig2(); fig3(); fig4(); fig5(); print(f'Wrote Figure1-Figure5 assets to {OUT}')
