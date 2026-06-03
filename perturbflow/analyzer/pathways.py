# =============================================================================
# perturbflow/analyzer/pathways.py
#
# Lightweight gene functional annotation module.
#
# Provides:
#   • A curated catalogue of ~80 gene sets covering hallmark biology,
#     immune processes, cell cycle, DNA damage, metabolism, and signalling.
#     All gene symbols are HGNC standard (human).  The catalogue is embedded
#     in this file so no network access or external files are required.
#   • run_enrichment(gene_list, background, top_n) — Fisher's exact test
#     enrichment (pure NumPy/math, no scipy dependency), BH-corrected.
#   • annotate_deg(deg_df, background) — adds pathway annotation columns to
#     a DEG DataFrame produced by _compute_deg() in deg.py.
#
# Outputs added to deg CSVs:
#   pathways   -- pipe-separated top-3 enriched pathway names (significant DEGs)
#   pathway_fdr -- lowest FDR across enriched pathways (float)
# =============================================================================

from __future__ import annotations

from math import log, comb
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Embedded gene-set catalogue
# Each entry: "PATHWAY_NAME": {"GENE1", "GENE2", ...}
# Sources: MSigDB Hallmarks (H), KEGG, Reactome (subset), curated immune sets.
# ---------------------------------------------------------------------------

GENE_SETS: Dict[str, List[str]] = {
    # ── Cell cycle & proliferation ────────────────────────────────────────
    "CELL_CYCLE_G1S":        ["CCND1","CCND2","CCND3","CDK4","CDK6","E2F1","E2F2","E2F3","RB1","CDKN1A","CDKN2A","CCNE1","CCNE2","CDC25A","MYC","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7"],
    "CELL_CYCLE_G2M":        ["CCNB1","CCNB2","CDK1","CDC20","BUB1","BUB1B","BUB3","MAD2L1","AURKA","AURKB","PLK1","CENPF","MKI67","TOP2A","BIRC5","CDKN1B","WEE1","CDC25C","CHEK1","CHEK2"],
    "MITOTIC_SPINDLE":       ["TUBB","TUBA1B","KIF11","KIF23","CENPE","KIF2C","CKAP5","KIFC1","KIF20A","PLK4","ASPM","NUSAP1","TPX2","INCENP","PTTG1"],
    "DNA_REPLICATION":       ["POLA1","POLA2","POLB","POLD1","POLE","PCNA","RFC1","RFC2","RFC3","RPA1","RPA2","LIG1","FEN1","TYMS","RRM1","RRM2","GINS1","GINS2","MCM10"],

    # ── DNA damage & repair ──────────────────────────────────────────────
    "DNA_DAMAGE_RESPONSE":   ["TP53","ATM","ATR","BRCA1","BRCA2","CHEK1","CHEK2","MDM2","CDKN1A","H2AFX","RAD51","FANCD2","FANCA","NBN","MRE11","RAD50","PARP1","XRCC1","MLH1","MSH2","MSH6","PMS2"],
    "HOMOLOGOUS_RECOMB":     ["RAD51","BRCA1","BRCA2","RAD52","PALB2","BARD1","RBBP8","MRE11","NBN","RAD50","XRCC3","RAD54L","RECA","EME1","MUS81","SLX1A"],
    "MISMATCH_REPAIR":       ["MLH1","MLH3","MSH2","MSH3","MSH6","PMS1","PMS2","PCNA","RFC1","EXO1","POLG","LIG1"],
    "P53_PATHWAY":           ["TP53","MDM2","CDKN1A","GADD45A","BAX","BBC3","PMAIP1","FAS","TNFRSF10B","SFN","PERP","DDB2","POLK","SESN1","SESN2","RRM2B","BTG2","NOTCH1"],

    # ── Apoptosis & survival ─────────────────────────────────────────────
    "APOPTOSIS":             ["CASP3","CASP6","CASP7","CASP8","CASP9","CASP10","BAX","BAD","BCL2","BCL2L1","MCL1","BID","FADD","APAF1","CYCS","DIABLO","XIAP","BIRC3","FAS","TNFRSF10A","TNFRSF10B"],
    "ANTIAPOPTOSIS_BCL2":    ["BCL2","BCL2L1","BCL2L2","MCL1","BCL2A1","BCL2L10","BCLAF1"],
    "P53_MEDIATED_APOPTOSIS":["TP53","PMAIP1","BBC3","BAX","TNFRSF10B","FAS","PERP","AIFM2"],

    # ── Protein homeostasis ──────────────────────────────────────────────
    "PROTEASOME":            ["PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","UBA1","UBB","UBC","UBE2D1"],
    "UNFOLDED_PROTEIN_RESP": ["ATF6","ATF4","DDIT3","ERN1","EIF2AK3","HSPA5","HYOU1","HERPUD1","XBP1","PERK","EDEM1","DERL1","SEL1L","DNAJB9","CALR"],
    "HEAT_SHOCK_RESPONSE":   ["HSPA1A","HSPA1B","HSPA4","HSPA5","HSPA8","HSPB1","HSP90AA1","HSP90AB1","HSPH1","DNAJB1","DNAJC3","STIP1","BAG1","CHIP"],

    # ── Ribosome & translation ────────────────────────────────────────────
    "RIBOSOME_BIOGENESIS":   ["RPL3","RPL4","RPL5","RPL6","RPL7","RPL8","RPL10","RPL11","RPL12","RPL13","RPL14","RPL15","RPL18","RPL19","RPL22","RPL23","RPL24","RPL26","RPL27","RPL28","RPL29","RPL30","RPL31","RPL32","RPS2","RPS3","RPS4X","RPS6","RPS7","RPS8","RPS9","RPS10","RPS11","RPS12","RPS14","RPS15","RPS18","RPS19","RPS20","RPS23","RPS25","RPS27","RPS28","RPS29","NOP56","NOP58","FBL","DKC1","GAR1","NHP2","NOP10"],
    "mRNA_SPLICING":         ["SNRPA","SNRPB","SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF","SNRPG","SF3A1","SF3B1","SF3B2","U2AF1","U2AF2","PRPF8","PRPF3","PRPF4","HNRNPA1","HNRNPA2B1","HNRNPC","HNRNPD","SRSF1","SRSF2","SRSF3"],
    "TRANSLATION_INITIATION":["EIF1","EIF1A","EIF2A","EIF2S1","EIF2S2","EIF2S3","EIF3A","EIF3B","EIF3C","EIF3D","EIF3E","EIF4A1","EIF4E","EIF4G1","EIF5","EIF5B","EIF6","PAIP1"],

    # ── Metabolism ───────────────────────────────────────────────────────
    "GLYCOLYSIS":            ["HK1","HK2","GPI","PFKL","PFKM","ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA","LDHB","SLC2A1","SLC2A3"],
    "OXIDATIVE_PHOSPHORYLATION":["MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8","NDUFA1","NDUFA2","NDUFB1","UQCRB","UQCRC1","COX4I1","COX5A","ATP5A1","ATP5B","ATP5F1","SDHA","SDHB","FH","IDH1","IDH2"],
    "FATTY_ACID_METABOLISM": ["FASN","ACACA","ACACB","SCD","ELOVL1","ELOVL6","DGAT1","DGAT2","PPARA","PPARG","ACSL1","ACSL4","CPT1A","CPT2","HADHA","HADHB","ECHS1","ACADM","ACAD9","ACADVL"],
    "NUCLEOTIDE_SYNTHESIS":  ["IMPDH1","IMPDH2","PPAT","GART","PFAS","ADSL","ATIC","PRPS1","PRPS2","UMPS","CAD","CTPS1","DHODH","TYMS","DHFR","MTHFR"],
    "AMINO_ACID_METABOLISM": ["GLS","GLUD1","ASNS","PHGDH","PSAT1","PSPH","CBS","CTH","MAT1A","MAT2A","ALDH1A1","GOT1","GOT2","GPT","ARG1","ASS1","OTC"],
    "OXIDATIVE_STRESS":      ["SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","PRDX2","PRDX3","TXNRD1","TXN","NQO1","HMOX1","KEAP1","NFE2L2","GCLC","GCLM","G6PD","GSR"],

    # ── Transcription & chromatin ─────────────────────────────────────────
    "CHROMATIN_REMODELING":  ["SMARCA4","SMARCB1","SMARCC1","SMARCD1","ARID1A","ARID1B","ARID2","PBRM1","SETD2","KDM6A","KDM5C","KDM5A","JARID2","EP300","CREBBP","KAT2A","KAT2B","HDAC1","HDAC2","SIRT1","DNMT1","DNMT3A","DNMT3B","TET1","TET2"],
    "POLYCOMB_REPRESSION":   ["EZH1","EZH2","SUZ12","EED","BMI1","RING1","RNF2","PHC1","CBX2","CBX4","CBX7","CBX8","JARID2"],
    "TRANSCRIPTION_MYC":     ["MYC","MYCN","MAX","MXD1","MXI1","MAD2L2","MNT","BIN1","MYC_targets_V1","NPM1","NCL","ODC1","CCND2","CAD","CDK4"],
    "TRANSCRIPTION_NF_KB":   ["NFKB1","NFKB2","RELA","RELB","REL","NFKBIA","NFKBIB","IKBKA","IKBKB","IKBKG","TNFAIP3","BIRC2","BIRC3","TRAF1","TRAF2","IL6","CXCL8","TNF","LTA","BCL2A1"],
    "MTORC1_SIGNALING":      ["MTOR","RPTOR","MLST8","RPS6KB1","RPS6KB2","EIF4EBP1","ULK1","PIK3C3","DEPTOR","AKT1S1","RHEB","TSC1","TSC2","PTEN","PI3K_p85","AKT1"],
    "HYPOXIA":               ["HIF1A","EPAS1","ARNT","VHL","CA9","LDHA","ENO1","SLC2A1","VEGFA","ANGPT2","ADM","P4HA1","BNIP3","BNIP3L","PDK1","PFKFB3","ALDOC"],

    # ── Signalling pathways ──────────────────────────────────────────────
    "PI3K_AKT_MTOR":         ["PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","AKT1","AKT2","AKT3","PTEN","MTOR","RPTOR","RPS6KB1","EIF4EBP1","PDPK1","TSC1","TSC2","RHEB","GSK3B","FOXO1","FOXO3","MDM2"],
    "RAS_MAPK":              ["KRAS","NRAS","HRAS","BRAF","RAF1","MAP2K1","MAP2K2","MAPK1","MAPK3","SOS1","GRB2","NF1","DUSP1","DUSP6","SPRY1","ELK1","RSK1","RSK2","ETS1"],
    "WNT_BETA_CATENIN":      ["CTNNB1","APC","AXIN1","AXIN2","GSK3B","CK1A","LRP5","LRP6","FZD1","DVL1","NOTUM","TCF7","LEF1","TCF7L2","CCND1","MYC","AXIN2","DKK1","DKK3","SFRP1"],
    "NOTCH_SIGNALING":       ["NOTCH1","NOTCH2","NOTCH3","NOTCH4","DLL1","DLL3","DLL4","JAG1","JAG2","RBPJ","HES1","HES5","HEY1","HEY2","MAML1","MAML2","NRARP","PSEN1","ADAM10","FURIN"],
    "TGF_BETA":              ["TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD6","SMAD7","ACVR1","BMPR1A","BMP2","BMP4","GDF5","INHBA","LTBP1"],
    "HEDGEHOG":              ["SHH","IHH","DHH","PTCH1","PTCH2","SMO","GLI1","GLI2","GLI3","SUFU","CDON","BOC","GAS1","HHIP","DISP1"],
    "VEGF_SIGNALING":        ["VEGFA","VEGFB","VEGFC","VEGFD","PGF","FLT1","KDR","FLT4","NRP1","NRP2","PLCG1","MAPK1","PIK3CA","AKT1","STAT3","HIF1A"],
    "JAK_STAT":              ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","SOCS1","SOCS2","SOCS3","IL6","IFNG","IL2","IL4","IL10","EPO","GH"],
    "EGFR_ERBB":             ["EGFR","ERBB2","ERBB3","ERBB4","KRAS","BRAF","MAP2K1","MAPK1","MAPK3","PIK3CA","AKT1","STAT3","SRC","GRB2","SOS1","HRAS","MYC"],

    # ── Immune / inflammation ────────────────────────────────────────────
    "INTERFERON_ALPHA_RESP": ["IFIT1","IFIT2","IFIT3","IFIT5","IFITM1","IFITM2","IFITM3","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","OASL","RSAD2","HERC5","XAF1","TRIM22","SAMD9","BST2","CXCL10"],
    "INTERFERON_GAMMA_RESP": ["CIITA","HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-E","B2M","TAP1","TAP2","TAPBP","IRF1","IRF2","CXCL9","CXCL10","CXCL11","GBP1","GBP2","GBP5","IRGM","PSMB8","PSMB9","PSMB10"],
    "TNFA_NF_KB_SIGNALING":  ["NFKB1","RELA","TNFAIP3","TNFAIP6","TNFRSF9","BIRC3","IL8","IL6","CXCL1","CXCL2","CXCL3","CSF2","VCAM1","ICAM1","MMP9","CCL2","CCL20","PTGS2","IL1B","IKBKE","RELB"],
    "COMPLEMENT":            ["C1QA","C1QB","C1QC","C1R","C1S","C2","C3","C4A","C4B","C5","C6","C7","C8A","C8B","C8G","C9","CFB","CFD","CFH","CFHR1","CFI","CFP","CD55","CD59","CR1","C3AR1","C5AR1"],
    "IL2_STAT5_SIGNALING":   ["IL2","IL2RA","IL2RB","IL2RG","IL15","JAK1","JAK3","STAT5A","STAT5B","PI3KR1","PIK3CA","AKT1","MAPK1","MAPK3","MYC","BCL2","MCL1","CCND1","CCND2","CDK4"],
    "IL6_JAK_STAT3":         ["IL6","IL6R","IL6ST","JAK1","JAK2","TYK2","STAT3","SOCS1","SOCS3","CEBPB","BCL3","ICAM1","VEGFA","MMP1","FOS","JUN","TWIST1","SNAI1","CDH2"],
    "COAGULATION":           ["F2","F5","F7","F8","F9","F10","F11","F12","F13A1","SERPINC1","PROC","PROS1","TFPI","VWF","FGB","FGG","FGA","THBD","PLAT","PLAU","PLAUR","SERPINE1"],
    "INFLAMMATORY_RESPONSE": ["IL1A","IL1B","IL6","IL8","TNF","CXCL1","CXCL2","CXCL3","CXCL5","CXCL10","CCL2","CCL5","PTGS2","MMP9","ICAM1","VCAM1","SELE","SELP","TLR2","TLR4","MYD88","NFKB1"],

    # ── T / NK / B cell biology ───────────────────────────────────────────
    "T_CELL_ACTIVATION":     ["CD3D","CD3E","CD3G","CD247","ZAP70","LCK","LAT","SLP76","ITK","PLCG1","NFATC1","NFATC2","NFATC3","IL2","IFNG","TNF","CD28","CD2","CD4","CD8A","CD8B","ICOS","CD40LG"],
    "T_CELL_EXHAUSTION":     ["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","CD244","BTLA","CD160","TOX","TOX2","NR4A1","NR4A2","EOMES","TBX21","PRDM1","BATF","IRF4","ENTPD1","CX3CR1"],
    "NK_CELL_CYTOTOXICITY":  ["NCAM1","NCR1","NCR3","KLRC1","KLRD1","KLRK1","NKG2D","GZMB","GZMK","PRF1","FASLG","IFNG","TNF","KIR2DL1","KIR2DL3","EOMES","TBX21","FCGR3A","CD16A"],
    "B_CELL_RECEPTOR":       ["CD19","CD22","CD79A","CD79B","BLNK","BTK","LYN","SYK","PLCG2","CARD11","BCL10","MALT1","NFKB1","RELA","PI3KD","AKT1","FOXO1","PAX5","EBF1","IKZF1"],
    "CYTOKINE_CYTOKINE_RECV":["IL2","IL4","IL6","IL10","IL12A","IL12B","IL15","IL21","IFNG","TNF","TGFB1","CXCL8","CCL5","CCL19","CCL21","CXCL10","CXCL13","IL7","IL7R","IL2RA","IL6R","IFNGR1","TNFRSF1A"],

    # ── Cell identity / differentiation ──────────────────────────────────
    "EPITHELIAL_MESENCHYMAL": ["CDH1","CDH2","VIM","FN1","TWIST1","SNAI1","SNAI2","ZEB1","ZEB2","FOXC2","TGFB1","MMP2","MMP9","ITGAV","ITGB6","FGFR1","AXL","WNT5A"],
    "ANGIOGENESIS":          ["VEGFA","VEGFB","VEGFC","ANGPT1","ANGPT2","TEK","FLT1","KDR","PDGFB","PDGFRB","HIF1A","NOTCH1","DLL4","CXCL12","CXCR4","MMP2","MMP9","NRP1","NRP2","PECAM1","CDH5"],
    "STEM_CELL":             ["SOX2","OCT4","NANOG","KLF4","LIN28A","TDGF1","CD34","CD44","ALDH1A1","PROM1","ABCG2","NOTCH1","WNT3","LGR5","BMI1","EZH2","SALL4","UTF1","ZFP42","L1TD1"],

    # ── RNA & genome stability ────────────────────────────────────────────
    "TELOMERE_MAINTENANCE":  ["TERT","TERC","DKC1","TINF2","POT1","RAP1","TIN2","TPP1","TNKS","TNKS2","RTEL1","ATRX","DAXX","SMG6"],
    "RNA_DEGRADATION":       ["EXOSC1","EXOSC2","EXOSC3","EXOSC4","EXOSC5","EXOSC6","EXOSC7","EXOSC8","EXOSC9","DIS3","DIS3L","PAPD5","ZCCHC11","PARN","DCP1A","DCP2","XRN1","UPF1","UPF2","UPF3B","SMG1","SMG5","SMG7"],

    # ── Autophagy & lysosome ──────────────────────────────────────────────
    "AUTOPHAGY":             ["ATG1","ATG3","ATG4A","ATG5","ATG6","ATG7","ATG8A","ATG9A","ATG10","ATG12","ATG13","ATG14","ATG16L1","BECN1","PIK3C3","WIPI1","WIPI2","ULK1","ULK2","MAP1LC3A","SQSTM1","NBR1","BNIP3"],
    "LYSOSOME":              ["LAMP1","LAMP2","CTSD","CTSL","CTSB","CTSA","CTSZ","ACP2","HEXB","HEXA","GLA","GBA","SMPD1","NAGLU","GALNS","GNS","ARSB","IDUA","FUCA1","NEU1","MAN2B1"],

    # ── TF-Target gene sets (MSigDB TFT / ENCODE ChIP-seq consensus) ─────
    # Targets of major TFs commonly perturbed in CRISPRi screens
    "TFT_MYC_TARGETS_V1":    ["NPM1","NCL","ODC1","CCND2","CAD","CDK4","LDHA","ENO1","PTMA","HNRNPA1","RPS6","NME1","EIF4E","HSPA4","APEX1","TFAM","MCM5","MCM7","RRM1","RRM2","PCNA","TYMS","TK1","POLE2","POLD1","CDC6","CCNA2","CCNE1","E2F1","CDK2"],
    "TFT_MYC_TARGETS_V2":    ["GART","DHODH","PPAT","IMPDH2","ASNS","CTPS1","SLC7A5","SLC3A2","SLC1A5","SLC43A1","GLUD1","GOT1","PKM","G6PD","FASN","ACACA","MTHFD2","SHMT2","CAD","PAICS","ADSL","ATIC"],
    "TFT_E2F_TARGETS":       ["CCNE1","CCNE2","CCNA2","CDK2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC6","CDC45","ORC1","TYMS","TK1","RRM1","RRM2","DHFR","GINS1","GINS2","POLA1","POLE","POLE2","POLD1","RFC1","FEN1","LIG1","E2F1","E2F2","E2F3"],
    "TFT_STAT3_TARGETS":     ["BCL2","BCL2L1","MCL1","BIRC5","CCND1","MYC","VEGFA","MMP2","MMP9","CDK4","CDKN1A","FOS","JUN","IL6","IL10","IL2RA","SOCS1","SOCS3","PIM1","PIM2","TWIST1","HIF1A","STAT3","IRF1","CXCL10","GBP1","IFIT1","IFIT3"],
    "TFT_STAT5A_TARGETS":    ["CCND1","CCND2","PIM1","PIM2","BCL2","BCL2L1","IRF1","SOCS1","SOCS2","SOCS3","OSM","IL2RA","IL10","IGF1","IRS2","EGR1","CISH","LY6E","PRLR"],
    "TFT_NFKB_TARGETS":      ["TNFAIP3","BIRC2","BIRC3","BCL2A1","TRAF1","TRAF2","IL6","CXCL8","IL1A","IL1B","TNF","LTA","NFKBIA","NFKBIB","ICAM1","VCAM1","CCL2","CCL5","CXCL10","MMP9","PTGS2","NOS2","CD80","CD86","IRF2BP2"],
    "TFT_TP53_TARGETS":      ["CDKN1A","GADD45A","GADD45B","GADD45G","MDM2","BBC3","PMAIP1","BAX","FAS","TNFRSF10B","SESN1","SESN2","DDB2","POLK","RRM2B","BTG2","NOTCH1","SFN","PERP","TIGAR","GDF15","FDXR","TRIAP1","ZMAT3","AEN"],
    "TFT_IRF_TARGETS":       ["IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","OASL","RSAD2","IFI44","IFI44L","IFI6","HERC5","HERC6","XAF1","IRF7","IRF9","STAT1","STAT2","GBP1","GBP2","TRIM5","TRIM22","BST2","CXCL10","CXCL9","CXCL11"],
    "TFT_AP1_TARGETS":       ["FOS","FOSB","FOSL1","FOSL2","JUN","JUNB","JUND","ATF3","ATF4","CREB1","EGR1","MMP1","MMP3","MMP9","VEGFA","CCND1","CDKN1A","BCL2","MCL1","IL2","IL8","IL10","HMGA1","ALDH1A1","TP53","NOS2"],
    "TFT_E2F4_TARGETS":      ["RB1","RBL1","RBL2","HDAC1","HDAC2","HDAC3","SIN3A","RBBP7","RBBP4","TFDP1","TFDP2","MBD3","LIN9","MYBL2","BMYB","TGIF1","RBBP9","RBAP46","DP1","DP2"],
    "TFT_RUNX_TARGETS":      ["CSF1R","SPI1","CEBPA","CEBPB","GATA2","KLF4","MYC","BCL2","MCL1","BIRC5","CDKN1A","CDKN2A","HMGA2","HIF1A","VEGFA","MMP2","MMP9","ITGA4","CD44","IL3","CSF2","IL6","FLT3","KIT"],
    "TFT_GATA_TARGETS":      ["HBB","HBA1","ALAS2","EPOR","GATA1","GATA2","KLF1","FOG1","TAL1","LMO2","LYL1","ETO2","HMGA2","BCL2L1","MCM","CCND3","CDK6","GYPA","GYPB","TFRC","SLC4A1","ANK1","EBF1","PAX5"],
    "TFT_ETS_TARGETS":       ["SPI1","ETS1","ETS2","FLI1","ERG","ETV6","ETV1","ETV4","ETV5","GABPA","ELF1","ELF2","ELF4","ELK1","ELK3","ELK4","SPIB","SPIC","SPDEF","ERF","ETV3","ETV7","FEV","PTTG1","MMP9","VEGFA","BCL2","IL2","IL3","CDKN1A","ICAM1","FAS"],

    # ── KEGG / Reactome representative pathways ─────────────────────────
    "KEGG_T_CELL_RECEPTOR":  ["CD3D","CD3E","CD3G","CD247","ZAP70","LCK","FYN","CD4","CD8A","LAT","SLP76","ITK","PLCG1","RASGRP1","MAP2K1","MAPK3","NFATC1","NFATC2","NFKB1","RELA","PIK3CD","AKT1","PDK1","PTEN","GRB2","SOS1","HRAS","KRAS","RAF1"],
    "KEGG_B_CELL_RECEPTOR":  ["CD79A","CD79B","BTK","SYK","BLNK","PLCG2","CD19","CD21","CD81","PIK3R1","PIK3CD","AKT1","NFKB1","RELA","IKK1","IKK2","IKBKG","VAV1","RAC1","CDC42","MAP2K1","MAPK3","MAPK8","JUN","FOS","BRAF","HRAS"],
    "KEGG_NATURAL_KILLER_CELL":["KLRB1","KLRD1","KIR2DL1","KIR2DL3","KIR3DL1","NCR1","NCR2","NCR3","NKG2D","DNAM1","LFA1","VAV1","SYK","ZAP70","SLP76","LAT","PLCG2","PI3K","AKT1","MAPK1","MAPK3","JUN","FOS","GZMB","PRF1","FASLG","IFNG","TNF"],
    "KEGG_JAK_STAT":         ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","SOCS1","SOCS2","SOCS3","PIAS1","PIAS3","IL2","IL4","IL6","IL7","IL10","IL12A","IL15","IL21","IFNA1","IFNB1","IFNG","LEPR","EPOR","GHR","PRLR"],
    "KEGG_MAPK":             ["MAP2K1","MAP2K2","MAPK1","MAPK3","MAP2K3","MAP2K6","MAPK8","MAPK9","MAPK10","MAP2K4","MAP2K7","MAPK11","MAPK12","MAPK13","MAPK14","BRAF","RAF1","ARAF","KRAS","HRAS","NRAS","SOS1","GRB2","EGFR","FGFR1","PDGFRA","MET","ATF2","ELK1","JUN","FOS","TP53","CDKN1A"],
    "KEGG_PI3K_AKT":         ["PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","AKT1","AKT2","AKT3","PTEN","PDK1","MTOR","TSC1","TSC2","RPTOR","RPS6KB1","EIF4EBP1","GSK3B","FOXO1","FOXO3","MDM2","BCL2","BCL2L1","CCND1","CDK4","CDKN1B"],
    "KEGG_CYTOKINE_RECEPTOR": ["IL1R1","IL2RA","IL4R","IL6R","IL10RA","IL12RB1","IL15RA","IL17RA","IL21R","IFNGR1","IFNAR1","TNFRSF1A","TNFRSF1B","CD40","IL18R1","IL23R","LEPR","GHR","EPOR","CSFR1","FGFR1","EGFR","MET","ALK","RET","NTRK1"],
    "REACTOME_INNATE_IMMUNE": ["TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","MYD88","TRIF","TRAM","TIRAP","IRAK1","IRAK4","TRAF6","NFKB1","RELA","IRF3","IRF7","IFNB1","TNF","IL6","IL12A","IL12B","NLRP3","PYCARD","CASP1","IL1B","IL18","STING1","cGAS"],
    "REACTOME_MHC_CLASS_I":  ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP","CALR","CANX","ERp57","PDIA3","ERAP1","ERAP2","PSME1","PSME2","PSMB8","PSMB9","PSMB10","UBB","UBC","NEDD8","SEC61A1"],
    "REACTOME_TCR_SIGNALING": ["CD3D","CD3E","CD3G","CD247","ZAP70","LCK","LAT","SLP76","ITK","VAV1","GADS","GRB2","SOS1","RASGRP1","PLCG1","IP3K","DAG","PKC","NFATC1","NFKB1","AP1","CD28","CTLA4","PD1","PDL1","PI3K","AKT1","FOXP3"],
    "REACTOME_CYTOKINE_SIGNALING":["IFNG","IFNGR1","IFNGR2","JAK1","JAK2","STAT1","IRF1","GBP1","CXCL9","CXCL10","CXCL11","IL6","IL6R","JAK1","JAK2","STAT3","SOCS1","SOCS3","IL10","IL10RA","IL10RB","IL2","IL2RA","JAK3","STAT5A","STAT5B"],
    "REACTOME_DNA_REPAIR":   ["BRCA1","BRCA2","RAD51","RAD52","PALB2","FANCD2","ATM","ATR","CHEK1","CHEK2","RPA1","RFC1","PCNA","POLD1","POLE","LIG1","LIG3","XRCC1","XRCC4","LIG4","NHEJ1","DCLRE1C","PRKDC","ARTEMIS","H2AFX","53BP1","MDC1","RNF8","RNF168","BRCC3"],
}

# Flat reverse mapping: gene -> list of pathway names
_GENE_TO_PATHWAYS: Dict[str, List[str]] = {}
for _pw, _genes in GENE_SETS.items():
    for _g in _genes:
        _GENE_TO_PATHWAYS.setdefault(_g, []).append(_pw)


# ---------------------------------------------------------------------------
# Hypergeometric / Fisher exact test (pure Python, no scipy)
# ---------------------------------------------------------------------------

def _log_comb(n: int, k: int) -> float:
    """log(C(n,k)) computed via lgamma for large n."""
    from math import lgamma
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def _hypergeom_sf(k: int, M: int, n: int, N: int) -> float:
    """Survival function P(X >= k) for the hypergeometric distribution.

    Args:
        k  -- observed overlap (k-1 so we get P >= k)
        M  -- total background gene count
        n  -- gene set size
        N  -- query gene set size
    """
    # Sum P(X=i) for i in [k, min(n,N)]
    lo = max(0, N + n - M)
    hi = min(n, N)
    if k > hi:
        return 0.0
    total = 0.0
    for i in range(max(k, lo), hi + 1):
        lp = _log_comb(n, i) + _log_comb(M - n, N - i) - _log_comb(M, N)
        total += pow(2.718281828, lp)
    return min(total, 1.0)


def _bh(pvals: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    result = [1.0] * n
    cummin = 1.0
    for rank, idx in enumerate(reversed(order)):
        raw = pvals[idx] * n / (n - rank)
        cummin = min(cummin, raw)
        result[idx] = cummin
    return result


def run_enrichment(
    query_genes: List[str],
    background_genes: List[str],
    top_n: int = 50,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Test all gene sets for enrichment in query_genes vs background.

    Args:
        query_genes     -- DEG gene list (e.g. significant genes).
        background_genes-- All genes tested (full gene list in adata).
        top_n           -- Return only the top N results by p-value.
        fdr_threshold   -- FDR cutoff used to flag significant enrichments.

    Returns:
        DataFrame with columns:
          pathway, n_pathway, n_overlap, n_query, n_background,
          pval, fdr, odds_ratio, significant, overlap_genes
        Sorted by pval ascending.
    """
    if not query_genes or not background_genes:
        return pd.DataFrame()

    bg_set = set(background_genes)
    q_set = set(query_genes) & bg_set
    M = len(bg_set)
    N = len(q_set)
    if N == 0 or M == 0:
        return pd.DataFrame()

    rows = []
    for pw, pw_genes in GENE_SETS.items():
        pw_bg = [g for g in pw_genes if g in bg_set]
        n = len(pw_bg)
        if n == 0:
            continue
        overlap = [g for g in pw_bg if g in q_set]
        k = len(overlap)
        if k == 0:
            continue
        pval = _hypergeom_sf(k, M, n, N)
        # Odds ratio (with pseudocount of 0.5 for stability)
        a = k + 0.5
        b = N - k + 0.5
        c = n - k + 0.5
        d = M - n - N + k + 0.5
        or_ = (a * d) / (b * c)
        rows.append({
            "pathway": pw,
            "n_pathway": n,
            "n_overlap": k,
            "n_query": N,
            "n_background": M,
            "pval": pval,
            "odds_ratio": round(or_, 3),
            "overlap_genes": "|".join(sorted(overlap)),
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("pval")
    pvals = df["pval"].tolist()
    df["fdr"] = _bh(pvals)
    df["significant"] = df["fdr"] < fdr_threshold
    return df.head(top_n).reset_index(drop=True)


def annotate_deg(
    deg_df: pd.DataFrame,
    background_genes: List[str],
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Add per-gene pathway annotations to a DEG DataFrame.

    For each gene in deg_df, looks up which pathway(s) it belongs to in the
    catalogue and records them in a 'pathway_membership' column.

    Also runs set-level enrichment on the significant DEGs and adds two
    columns to deg_df:
        pathway_membership  -- pipe-joined pathways the gene belongs to
        top_pathway         -- the most enriched significant pathway
                               (shown for every gene in the DEG table)

    Args:
        deg_df          -- DataFrame from _compute_deg() (must have 'gene' col).
        background_genes-- All genes tested; used as enrichment background.
        fdr_threshold   -- FDR cutoff for enrichment significance.

    Returns:
        deg_df with two new columns appended (in place copy).
    """
    if deg_df.empty:
        return deg_df

    # Per-gene membership.
    deg_df = deg_df.copy()
    deg_df["pathway_membership"] = deg_df["gene"].map(
        lambda g: "|".join(_GENE_TO_PATHWAYS.get(g, []))
    )

    # Enrichment on significant genes only.
    sig_genes = deg_df.loc[deg_df["significant"], "gene"].tolist()
    enrich_df = run_enrichment(sig_genes, background_genes, top_n=5, fdr_threshold=fdr_threshold)

    if enrich_df.empty or not enrich_df["significant"].any():
        deg_df["top_pathway"] = ""
    else:
        top_pw = enrich_df[enrich_df["significant"]].iloc[0]["pathway"]
        deg_df["top_pathway"] = top_pw

    return deg_df
