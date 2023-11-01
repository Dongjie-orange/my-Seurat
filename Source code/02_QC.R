
cdj_basic_qc <- function(input_sce){
  
  input_sce <- sce.all
    
  if (sp=='human') {
    
    # 01 - mito genes
    mito_genes <- rownames(input_sce)[grep("^MT-", rownames(input_sce),ignore.case = T)]
    input_sce <- PercentageFeatureSet(input_sce, "^MT-", col.name = "percent.mt")
    # 02  - hb genes
    hb_genes <- rownames(input_sce)[grep("^HB[^(P)]", rownames(input_sce),ignore.case = T)]
    input_sce <- PercentageFeatureSet(input_sce, "^HB[^(P)]", col.name = "percent.hb")
    # 03 - ribo genes
    ribo_genes <- rownames(input_sce)[grep("^RP[SL]", rownames(input_sce),ignore.case = T)]
    input_sce <- PercentageFeatureSet(input_sce, "^RP[SL]",col.name = "percent.ribo")
    
  }else if(sp=='mouse'){
    
    # 01 - mito genes
    mito_genes <- rownames(input_sce)[grep("^mt-", rownames(input_sce),ignore.case = T)]
    input_sce <- PercentageFeatureSet(input_sce, "^mt-", col.name = "percent.mt")
    # 02  - hb genes
    hb_genes <- rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce),ignore.case = T)]
    input_sce <- PercentageFeatureSet(input_sce, "^Hb[^(p)]", col.name = "percent.hb")
    # 03 - ribo genes
    ribo_genes <- rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce),ignore.case = T)]
    input_sce <- PercentageFeatureSet(input_sce, "^Rp[sl]",col.name = "percent.ribo")
    
  }else {
    print('we only accept human or mouse')
  }
  
  # 04 - log10GenesPerUMI
  input_sce$log10GenesPerUMI <- log10(input_sce$nFeature_RNA)/log10(input_sce$nCount_RNA)
  
  p <- VlnPlot(input_sce,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
               ncol = 5,
               pt.size = 0)
  ggsave('Vlnplot_before_filter.pdf',p,width = 15,height = 5)
  
  # filter
  # [ REMOVE - Cells expressing <100 or >4,000 genes, >10% mitochondrial reads, >2% hemoglobin reads  log10(UMI per gene) < 0.6]
  
  # input_sce.filt <- subset(input_sce, 
  #                          subset = (nFeature_RNA > 200) & (nFeature_RNA < 5000) & (nCount_RNA > 500) & (percent.mt < 10) & (percent.ribo < 25) & (percent.hb < 3) & (log10GenesPerUMI > 0.6))
  # 
  
  input_sce.filt <- subset(input_sce, 
                           subset = (nFeature_RNA > 200) & (percent.mt < 10))

  
  p <- VlnPlot(input_sce.filt,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
               ncol = 5,
               pt.size = 0)
  ggsave('Vlnplot_after_filter.pdf',p,width = 15,height = 5)
  
  
  return(input_sce.filt)

}
