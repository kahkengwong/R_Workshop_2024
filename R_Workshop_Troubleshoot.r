#####
# Troubleshooting 1
# If there is a problem with installing the org.Hs.eg.db package, download the file here: https://bioconductor.org/packages/3.18/data/annotation/src/contrib/org.Hs.eg.db_3.18.0.tar.gz  

# # Then put the file into your folder of choice and execute install.packages with the `...` location below that contains the tar.gz file install.packages("C:/.../org.Hs.eg.db_3.18.0.tar.gz", repos = NULL, type = "source")

# After installing, library the package
library(org.Hs.eg.db)

# If there is still problem, update and library the fastmap package and see if it solves the issue
install.packages("fastmap")
library(fastmap)

# Perform the same steps above for any other Bioconductor packages that fail to install using the install.packages command.
#####

#####
# Troubleshooting 2
# If there is a problem with downloading and extracting the list of TGFb genes from Hallmark's JSON file, manually define the genes
# Alternative: Manually define the TGFb genes (extracted manually from Hallmark's database)
hallmark_genes <- c("ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1","CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1","HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2","LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA","PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1","TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3","WWTR1","XIAP") 

# Skip the following steps as you have now defined hallmark_genes as above for the TGFb signalling pathway:
# URL for the HALLMARK_TGF_BETA_SIGNALING gene set JSON file
json_url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_TGF_BETA_SIGNALING&fileType=json"

# Download the JSON file and read it directly into R
json_data <- fromJSON(json_url)

# Extract the list of gene symbols
hallmark_genes <- json_data$HALLMARK_TGF_BETA_SIGNALING$geneSymbols


# If there is also a problem with downloading and extracting the list of TNFa genes from Hallmark's JSON file, manually define the genes
# Alternative: Manually define the TNFa genes (extracted manually from Hallmark's database)
hallmark_genes_TNFa <- c("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIGI","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")

# Skip the following steps as you have now defined hallmark_genes as above for the TNFa signalling pathway
# URL for the HALLMARK_TNFA_SIGNALING_VIA_NFKB gene set JSON file
json_url_TNFa <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_TNFA_SIGNALING_VIA_NFKB&fileType=json"

# Download the JSON file and read it directly into R
json_data_TNFa <- fromJSON(json_url_TNFa)

# Extract the list of gene symbols
hallmark_genes_TNFa <- json_data_TNFa$HALLMARK_TNFA_SIGNALING_VIA_NFKB$geneSymbols