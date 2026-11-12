// Post-process gene-set names: restore acronyms lost by sentence-casing (run after build_human_luad.R)
const fs = require('fs'); const f = process.argv[2] || '../site/data/human_luad_genesets.json'; const g = JSON.parse(fs.readFileSync(f));
const AC = 'DNA RNA UV KRAS MYC E2F G2M TNFA NFKB IL2 IL6 JAK STAT3 STAT5 TGF MTORC1 PI3K AKT MTOR P53 WNT ATP NADH NADPH MHC ER GTP GTPASE TCA ROS ECM HDL LDL VEGF EGFR ERBB ERBB2 FGFR PDGF NMDA AMPA GABA TRNA RRNA MRNA MIRNA SNRNA NCRNA LNCRNA SRP ERK MAPK MAP JNK NF KAPPAB IGF IGF1 HIF BMP SMAD TLR TP53 APC RHO RAC CDC42 II I III IV B T NK G1 S M G2 CD4 CD8 HLA UPR SUMO ISG15 NEDD8 PTEN BRCA1 RAS RAF GPCR CAMP CGMP IP3 COPII COPI ESCRT SARS COV HIV EGR1 TNF IFN IL1 IL10 IL12 IL17 IL4 IL13 FC FCERI FCGR PD1 CTLA4 NOTCH HH SHH CTNNB1 YAP TAZ AMPK SREBP PPAR PPARA RXR LXR FXR VDR GLI FOXO NRF2 KEAP1 ATF4 ATF6 IRE1 PERK XBP1'.split(' ');
const set = new Set(AC);
const fix = n => n.split(' ').map(w => { const U = w.toUpperCase(); if (U === 'DN') return 'down'; return set.has(U.replace(/[^A-Z0-9]/g, '')) ? U : w; }).join(' ').replace(/^./, c => c.toUpperCase());
for (const c of Object.values(g.catalogs)) { c.names = c.names.map(fix); if (!Array.isArray(c.members)) c.members = c.ids.map(id => c.members[id]); }
fs.writeFileSync(f, JSON.stringify(g)); console.log('prettified', f);
