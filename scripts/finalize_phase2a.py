from pathlib import Path
import pandas as pd

root = Path('.')
out = root / 'results/phase2a'
x = pd.read_csv(out / 'hdll_local_rg_all.tsv', sep='\t')
reg = pd.read_csv(root / 'metadata/hdll_region_reference.tsv', sep='\t')
sens = pd.read_csv(out / 'hdll_noukbb_sensitivity.tsv', sep='\t')
priority = x[(x.bivariate_interpretable) & x.trait_pair.isin(['AMD-POAG','POAG-RE','AMD-RE'])].sort_values(['trait_pair','p']).groupby('trait_pair').head(10).copy()
priority = priority.merge(reg[['chr','region','region_id']], on=['chr','region'], how='left')
priority = priority.merge(sens[['trait_pair','chr','region','direction_stable']], on=['trait_pair','chr','region'], how='left')
priority['priority_rank'] = priority.groupby('trait_pair').cumcount() + 1
priority['noukbb_stable'] = priority['direction_stable'].fillna(False)
priority['lava_support'] = 'NO_ELIGIBLE_TARGET'
priority['global_leaveout_effect'] = 'NOT_RUN_NO_GLOBAL_OR_PAIR_FDR_TARGET'
priority['phase2b_priority'] = 'TIER3'
priority['warnings'] = priority['warning'].fillna('')
cols = ['priority_rank','chr','start','end','region_id','trait_pair','h2_trait1','h2_trait2','local_covariance','local_rg','rg_lower','rg_upper','p','pair_fdr','global_fdr','direction','noukbb_stable','lava_support','global_leaveout_effect','warnings','phase2b_priority']
priority[cols].to_csv(out / 'LOCAL_PLEIOTROPY_PRIORITY_TABLE.tsv', sep='\t', index=False)
pd.DataFrame(columns=['trait_pair','chr','region','start','end','hdll_local_rg','lava_local_rg','direction_agreement','inference_agreement','status']).to_csv(out / 'lava_targeted_sensitivity.tsv', sep='\t', index=False)
pd.DataFrame(columns=['trait_pair','chr','region','start','end','original_rg','leave_region_out_rg','change_rg','se','p','status']).to_csv(out / 'leave_region_out_global_rg.tsv', sep='\t', index=False)

def summary(pair):
    y = x[(x.trait_pair == pair) & x.bivariate_interpretable]
    return len(y), int((y.pair_fdr < .05).sum()), int((y.global_fdr < .05).sum()), int((y.local_rg > 0).sum()), int((y.local_rg < 0).sum())

lines = ['# RetinaPlex-State Phase 2A Local Architecture', '', 'Date: 2026-09-01', '', '## Verdict', '', '`DIFFUSE_ARCHITECTURE`', '', 'No region passed global or pair-specific FDR after local h2 eligibility filtering. The strong global POAG-RE relationship is interpreted as diffuse at the available regional resolution; AMD-POAG has no robust local hotspot, and AMD-RE has no robust opposing local signals.', '', '## Pairwise summary', '', '| Pair | Interpretable | Pair-FDR significant | Global-FDR significant | Positive | Negative |', '|---|---:|---:|---:|---:|---:|']
for pair in ['AMD-POAG','POAG-RE','AMD-RE']:
    a,b,c,d,e = summary(pair)
    lines.append(f'| {pair} | {a} | {b} | {c} | {d} | {e} |')
lines += ['', '## Sensitivity', '', f'No-UKBB sensitivity was completed for {len(sens)} selected top regions from AMD-POAG and POAG-RE. Directional stability was retained for {int(sens.direction_stable.sum())}/{len(sens)} evaluated results. No LAVA targets or leave-one-region-out targets met the locked selection rule.', '', '## Interpretation', '', '- POAG-RE global sharing is not represented by a small set of robust local HDL-L hotspots.', '- AMD-POAG remains suggestive globally without a reproducible local hotspot.', '- AMD-RE does not support robust canceling local pleiotropy.', '- The results do not justify Phase 2B SNP-level locus discovery.', '', '## Next recommended phase', '', '`D. polygenic/pathway analysis instead of locus discovery`', '', 'No cell-state or biological annotation analysis was started in Phase 2A.']
(root / 'reports/PHASE2A_LOCAL_ARCHITECTURE.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
print('finalized Phase 2A tables and report')
