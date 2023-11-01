

#----------------------------------------------------------------------------------
#  Step 0: singleR check
#----------------------------------------------------------------------------------

library(SingleR)
library(ggthemes)

if (sp=='human') {
  
  ref <- celldex::HumanPrimaryCellAtlasData()
  pbmc_counts <- GetAssayData(sce.all.int,
                              slot = "data")
  pred <- SingleR(test=pbmc_counts,
                  ref=ref,
                  labels = ref$label.main)
  pred
  sce.all.int$singleR.labels <- pred$labels[match(rownames(sce.all.int@meta.data),rownames(pred))]
  DimPlot(sce.all.int, label = F, pt.size = 0.5,cols = my36colors,group.by = 'singleR.labels',reduction = method)+
    #NoLegend()+
    labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  ggsave('singleR.pdf',height = 8,width = 10)
  
}else if(sp=='mouse'){
  
  ref <- celldex::MouseRNAseqData()
  pbmc_counts <- GetAssayData(sce.all.int,
                              slot = "data")
  pred <- SingleR(test=pbmc_counts,
                  ref=ref,
                  labels = ref$label.main)
  pred
  sce.all.int$singleR.labels <- pred$labels[match(rownames(sce.all.int@meta.data),rownames(pred))]
  DimPlot(sce.all.int, label = F, pt.size = 0.5,cols = my36colors,group.by = 'singleR.labels',reduction = method)+
    #NoLegend()+
    labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  ggsave('singleR.pdf',height = 8,width = 10)
  
}else {
  print('we only accept human or mouse')
}


#----------------------------------------------------------------------------------
#  Step 0.1: cellcycle check
#----------------------------------------------------------------------------------

library(ggthemes)
sce.all.int <- CellCycleScoring(sce.all.int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

p <- DimPlot(sce.all.int, reduction = method, pt.size = 1, label = FALSE,
             group.by = "Phase", cols = my36colors) + theme_few()
ggsave('02_CC.phase.pdf', p, width = 8.5, height = 7)

p <- VlnPlot(sce.all.int, features = "S.Score", pt.size = 0,
             group.by = "seurat_clusters", cols = my36colors) + theme_few()
ggsave("02_CC.score.S.Score.pdf", p, width = 10, height = 5)

p <- VlnPlot(sce.all.int, features = "G2M.Score", pt.size = 0,
             group.by = "seurat_clusters", cols = my36colors) + theme_few()
ggsave("02_CC.score.G2M.Score.pdf", p, width = 10, height = 5)




#----------------------------------------------------------------------------------
#  Step 0.2: marker check
#----------------------------------------------------------------------------------

library(cowplot)
library(ggpubr)
sce.all.int <- BuildClusterTree(object = sce.all.int)
all.markers <- FindAllMarkers(object = sce.all.int, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1)
all.markers <- all.markers[which(all.markers$p_val_adj < 0.05 & all.markers$avg_log2FC > 0), ]
write.xlsx(all.markers,  "03_top_markers.xlsx", overwrite = T)


all.markers <- all.markers[which(all.markers$pct.1 > 0.25), ]
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unique(top10$gene)
p <- DotPlot(sce.all.int, features = gene.list, dot.scale = 8, cols = c("#DDDDDD", "#003366" ), col.min = -2) + RotatedAxis()
p <- p + theme_few() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))
p <- p + theme(axis.text.y = element_text(size = 20))
p <- p + scale_size(range = c(1, 7))
p <- p + gradient_color(c("#EEEEEE","#ffb459","#e8613c","#b70909"))
ggsave("03_top_markers_dotplot.pdf", p, width = 20, height = 7)




#----------------------------------------------------------------------------------
#  Step 1: marker list
#----------------------------------------------------------------------------------


#doi: 10.1016/j.canlet.2022.215834; 10.1158/2159-8290.CD-19-0094; 10.1038/s41422-019-0195-y; 10.1038/s41467-023-40314-w; 10.1016/j.celrep.2023.112620.; 10.1186/s12967-023-04051-4
# 10.1038/s41420-021-00663-1
pdac_all_markers_list <- list(
  Tcells = c('CD2', 'CD3D','CD3E','CD3D33','CD4','CD8'),
  fibroblasts = c('LUM','COL1A1','COL3A1','DCN'),
  macrophages = c('AIF1', 'CD68','CD14','CD64','LYZ2', 'APOE', 'ARG1', 'CTSS'),
  ductal_all = c('MMP7','TSPAN8','SOX9','LCN2'),
  ductal1 = c('AMBP', 'CFTR', 'MMP7'),
  ductal2 = c('KRT19', 'KRT7', 'TSPAN8', 'SLPI'), #malignant
  endothelial = c('CDH5','CDH5','PLVAP', 'VWF','CLDN5'),
  stellate = c('RGS5','ACTA2', 'PDGFRB', 'ADIRF'),
  acinar = c('PRSS1','CTRB1','CTRB2','REG1B'),
  lymphocytes = c('CD3D','IL7R','CD3G'),
  Bcells = c('MS4A1','CD79A', 'CD79B', 'CD52'),
  plasma = c('MZB1','IGLC1','IGHG1','IGHG1','IGKC','IGHG3','IGHA1','CD27','CD38','SLAMF7','TNFRSF17'),
  mast = c('CPA3','TPASB1','TPSB2','CPA3','TPSAB1','MS4A2','HPGDS','CTSG','AREG'),
  neutrophils = c('G0S2', 'S100A8'),
  granulocytes = c('S100A8', 'S100A9', 'RETNLG','NGP'),
  schwann = c('GAP43','CDH19','PLP1','SOX10','CRYAB','PMP22'),
  endocrine = c('CHGA','CHGB','INS','IAPP'),
  epithelial = c('EPCAM','CDH1','CK19'),  # cellmarker
  NK = c('CD161','NCR1','KLRC1','TMIGD2','FCER1A'),  # cellmarker
  DC = c('CD1C','CD1','CD1A','CD21','CD34'),
  CSC = c('ALDH1A1','NES')
)


cancer_letter_marker <-  list(
  Tcells = c('CD3D','CD3E','CD4','CD8'),
  fibroblasts = c('LUM','COL1A1','DCN'),
  macrophages = c('AIF1', 'CD68','CD14','CD64'),
  ductal1 = c('AMBP', 'CFTR', 'MMP7'),
  ductal2 = c('KRT19', 'KRT7', 'TSPAN8', 'SLPI'), #malignant
  endothelial = c('CDH5','PLVAP', 'VWF','CLDN5'),
  stellate = c('RGS5','ACTA2', 'PDGFRB', 'ADIRF'),
  acinar = c('PRSS1','CTRB1','CTRB2','REG1B'),
  Bcells = c('MS4A1','CD79A', 'CD79B', 'CD52'),
  endocrine = c('CHGA','CHGB','INS','IAPP')
)



# doi:10.1016/j.canlet.2022.215834.
pdac_T_cell_sub <- list(
  All = c('CD2', 'CD3D','CD3E'),
  naive_T =c('CCR7', 'TCF7','SELL'),
  Tregs = c('CTLA4', 'FOXP3', 'IL2RA', 'LAYN', 'TIGIT','TNFRSF4','TNFRSF18','CD127','CD25','CD3','CD4'),
  CD4_T = c('CD4','SELL'),
  CD8_T = c('CD103','CD8A','CD8B'),
  Exhausted_CD8_T = c('LAG3','PDCD1','TIGIT'),
  NKT = c('CD56'),
  TRM = c('CD103','CD49A','PDCD1'),#Non-recirculating tissue-resident memory T cell
  Th = c('TBX21','CD3','CD4','CD45')
  )


pdac_B_cell_sub <- list(
  All = c('CD19', 'CD79A','MS4A1','BLK','CD20','CD22','CD52'),
  Activated_B = c('CD19', 'CD25', 'CD30'),
  Plasma_cell = c('CD27', 'CD38', 'CD78', 'CD138', 'CD319','CD138','IL6'),
  Memory_cell = c('CD20', 'CD27', 'CD40', 'CD80','CXCR3', 'CXCR4', 'CXCR5', 'CXCR6'),
  Marginal_zone_B_cells = c('CD1', 'CD21', 'CD27'),
  Follicular_B_cells = c('CD21', 'CD22', 'CD23'),
  Regulatory_B_cells = c('CD1', 'CD5', 'CD21', 'CD24', 'TLR4','IL10')
)



pdac_macrophage_cell_sub <- list(
  All = c('CD68', 'CD14','AIF1','CD86','CD163','CD206','CD64'),
  M1 = c('CD86', 'IRF1'),
  M2 = c('CD163', 'MRC1', 'CD206')
)


# Tcell marker
Tcells_markers <-  c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
                   'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
                   'IFNG', 'CCL4', 'CCL3' ,
                   'PRF1' , 'NKG7')
### CD4T
CD4_markers_list <- list(
  Tc = c("CD3D","CD3E"),
  CD4 = c("CD4" ),
  Treg = c("TNFRSF4","BATF","TNFRSF18","FOXP3","IL2RA","IKZF2"),
  naive = c("CCR7","SELL","CD5"),
  Tfh = c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST"),#滤泡辅助性T细胞
  ILC = c("TNFRSF25","KRT81","LST1","AREG","LTB","CD69")
)

### CD8T
CD8_markers_list1 <- list(
  CD8 = c("CD8A","CD8B"),
  TN_TCM = c("CCR7","SELL","TCF7","LEF1"),
  TEM = c("GZMK"  ),
  TEFF = c("TBX21","FCGR3A","FGFBP2"),
  TRM = c("XCL1","XCL2","ITGAE","CD69"),
  IEL_T = c("TMIGD2"),
  yT1c = c("GNLY","PTGDS","GZMB","TRDC"),
  yT2c = c("TMN1","HMGB2","TYMS"),
  MAIT_T = c("SLC4A10")
)

CD8_markers_list2 <- list(
  CD8T = c("CD8A","CD8B"),
  MAIT = c("ZBTB16","NCR3","RORA"),
  ExhaustedCD8T = c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4"),
  MemoryCD8 = c("EOMES","ITM2C"),
  Resting_NK = c("XCL1","XCL2","KLRC1"),
  Cytotoxic_NK = c("CX3CR1","FGFBP2","FCGR3A","KLRD1"),
  Pre_exhausted = c("IFNG","PRF1","GNLY","GZMA","NKG7","GZMK")
)

cd4_and_cd8T_markers_list <- list(
  naive = c("CCR7","SELL","TCF7","IL7R","CD27","CD28","LEF1","S1PR1"),
  CD8Trm = c("XCL1","XCL2","MYADM"),
  NKTc = c("GNLY","GZMA"),
  Tfh = c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST"),
  th17 = c("IL17A","KLRB1","CCL20","ANKRD28","IL23R","RORC","FURIN","CCR6","CAPG","IL22"),
  CD8Tem = c("CXCR4","GZMH","CD44","GZMK"),
  Treg = c("FOXP3","IL2RA","TNFRSF18","IKZF2"),
  naive = c("CCR7","SELL","TCF7","IL7R","CD27","CD28"),
  CD8Trm = c("XCL1","XCL2","MYADM"),
  MAIT = c("KLRB1","ZBTB16","NCR3","RORA"),
  yT1c = c("GNLY","PTGDS","GZMB","TRDC"),
  yT2c = c("TMN1","HMGB2","TYMS"),
  yt = c("TRGV9","TRDV2")
)

# CD20 (MS4A1)表达于除plasma B 之外的所有B，很关键的区分naive 和plasma的marker
# SDC1 = CD138 plasma B （接受抗原，可表达抗体）
Bcels_markers_list <- list(
  All = c('MS4A1','SDC1','CD27','CD38','CD19', 'CD79A'),
  GC_B = c('IL4R','TCL1A','LRMP','SUGCT'),
  IGA_plasm_B= c ( 'IGHA1'),
  IGG_plasm_B= c ( 'IGHG1')
)

myeloids_markers_list1 <- list(
  CM=c("TTN","MYH7","MYH6","TNNT2") ,
  EC=c("VWF", "IFI27", "PECAM1","MGP"),
  FB=c("DCN", "C7" ,"LUM","FBLN1","COL1A2"),
  MP=c("CD163", "CCL4", "CXCL8","PTPRC"),
  SMC=c("ACTA2", "CALD1", "MYH11"),
  Tc=c("CD3D","CD3E"),
  DC1 = c( 'Clec9a', 'Xcr1',   'Wdfy4'),
  DC2 = c('Itgax', 'Sirpa',   'Cd209a'),
  mregDCs= c('Ccr7', 'Cd80', 'Cd200',   'Cd247') ,
  hypoxia=c('Hif1a', 'Slc2a1', 'Vegfa', 'Hmox1',
            'Bnip3', 'Nos2', 'Mmp2', 'Sod3',
            'Cited2', 'Ldha'),
  peric=c("ABCC9","PDGFRB","RGS5")
)

myeloids_markers_list2 <- list(pDC = c("CLEC4C","IRF7","TCF4","GZMB"),
                              cDC1 = c("XCR1","CLNK","CLEC9A"),
                              cDC2 = c("FCER1A","HLA-DPB1","HLA-DQB1","CD1E","CD1C","CLEC10A","HLA-DQA2"),
                              DC3 = c("CCL19","LAMP3","IDO1","IDO2","LAD1","FSCN1","CCR7","LY75","CCL22","CD40","BIRC3","NFKB2"),
                              Macrophages = c("APOC1","HLA-DRB5","C1QA","C1QB"),
                              RTMs = c("THBS1"),#Resident tissue macrophages
                              Lam = c("APOE"),#Lipid associated macrophages
                              Monocytes = c("LYZ","HLA-DRB1","TIMP1","S100A11","CXCL8","IL1B","PTGS2","S100A9","S100A8","MMP19"),
                              Mono_C = c('CD14'),#Mono_CD14
                              Mono_F = c('FCGR3A'),#Mono_FCGR3A
                              Mast = c('TPSAB1' , 'TPSB2'))


CAF_marker <- list(
  hsm_iCAF = c('IL6', 'PDGFRA', 'CFD', 'PLA2G2A', 'HAS1', 'CXCL2', 'CCL2', 'CLU', 'EMP1', 'LMNA'),
  hsm_myCAF = c("ACTA2","TAGLN",'MMP11','MYL9','HOPX','POSTN','TPM1','TPM2'),
  mus_iCAF = c("CLEC3B","COL14A1",'HAS1','IL6'),
  mus_myCAF = c("TAGIN","THY1",'COL12A1','THBS2'),
  mus_apCAF = c("H2-AB1",'CD74','SAA3','SLPI'),
  panCAF = c('COL1A1','COL1A2','PDPN','DCN'),
  cCAF = c('COL1A1', 'LUM','MMP11','FAP','SFRP2'),
  csCAF = c('C3', 'C7', 'CFB', 'CFD', 'CFH', 'CFI'),
  PSCs = c('RGS5', 'ADIRF', 'CRIP1', 'NDUFA4L2','NOTCH3' , 'PDGFA')
)


######## MARKER_LIST ########


#----------------------------------------------------------------------------------
#  Step 2: check gene
#----------------------------------------------------------------------------------



#p_umap <- DimPlot(sce.all.int, reduction = method,label = T,repel = T)  #method = tsne or umap
#p_umap

if(sp=='human'){

  if(F){
  pbapply::pblapply(markers, function(x){
    #x=markers[1]
    genes_to_check=str_to_upper(get(x))
    DotPlot(sce.all.int , features = genes_to_check )  +
      coord_flip() +
      theme(axis.text.x=element_text(angle=45,hjust = 1))

    h=length( genes_to_check )/6+3;h
    ggsave(paste('check_for_',x,'.pdf'),height = h)
  })
  }
  pbapply::pblapply(markers_list, function(x){
    # x=markers_list[1]
    genes_to_check = pbapply::pblapply(get(x), str_to_upper)
    dup=names(table(unlist(genes_to_check)))[table(unlist(genes_to_check))>1]
    genes_to_check = pbapply::pblapply(genes_to_check, function(x) x[!x %in% dup])

    DotPlot(sce.all.int , features = genes_to_check )  +
      # coord_flip() +
      theme(axis.text.x=element_text(angle=45,hjust = 1))

    w=length( unique(unlist(genes_to_check)) )/5+6;w
    ggsave(paste('check_for_',x,'.pdf'),width  = w)
  })

 # last_markers_to_check <<- str_to_upper(last_markers )

}else if(sp=='mouse'){
  if(F){
  pbapply::pblapply(markers, function(x){
    #x=markers[1]
    genes_to_check=str_to_title(get(x))
    DotPlot(sce.all.int , features = genes_to_check )  +
      coord_flip() +
      theme(axis.text.x=element_text(angle=45,hjust = 1))

    h=length( genes_to_check )/6+3;h
    ggsave(paste('check_for_',x,'.pdf'),height = h)
  })
  }
  pbapply::pblapply(markers_list, function(x){
    # x=markers_list[1]
    genes_to_check = pbapply::pblapply(get(x), str_to_title)
    dup=names(table(unlist(genes_to_check)))[table(unlist(genes_to_check))>1]
    genes_to_check = pbapply::pblapply(genes_to_check, function(x) x[!x %in% dup])

    DotPlot(sce.all.int , features = genes_to_check )  +
      # coord_flip() +
      theme(axis.text.x=element_text(angle=45,hjust = 1))

    w=length( unique(unlist(genes_to_check)) )/5+6;w
    ggsave(paste('check_for_',x,'.pdf'),width  = w)
  })

 # last_markers_to_check <<- str_to_title(last_markers )
}else {
  print('we only accept human or mouse')
}

if(F){
p_all_markers <- DotPlot(sce.all.int , features = last_markers_to_check )  +
  coord_flip() +
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers+p_umap
h=length( last_markers_to_check )/6+3;h
w=length( unique( Idents(sce.all.int)) )/5+10;w
ggsave(paste('last_markers_and_umap.pdf'),width  = w,height = h)
}









