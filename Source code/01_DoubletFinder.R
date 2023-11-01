

cdj_doublet_find <- function(input_list){
  
  datalist.doublet.fileter <- pbapply::pblapply(1:length(input_list),FUN = function(x){

  input_sce <- input_list[[x]] %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst",
                         nfeatures = 2000,
                         mean.cutoff=c(0.0125,3),
                         dispersion.cutoff =c(1.5,Inf)) %>% ScaleData()
  input_sce <- input_sce %>% RunPCA(features = VariableFeatures(object = input_sce)) %>% RunUMAP(dims= 1:15, verbose=TRUE)

  sweep.res.list_sce <- paramSweep_v3(input_sce, PCs = 1:30, sct = FALSE)
  head(sweep.res.list_sce)
  sweep.stats_sce<- summarizeSweep(sweep.res.list_sce, GT = FALSE)
  bcmvn_sce <- find.pK(sweep.stats_sce)
  mpK <- as.numeric(as.vector(bcmvn_sce$pK[which.max(bcmvn_sce$BCmetric)]))

  annotations <- input_sce@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  #DoubletRate <- ncol(input_sce)*8*1e-6  # calculate doublet propotion
  nExp_poi <- round(0.05*length(input_sce@active.ident))  ## ！！！！ tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  input_sce <- doubletFinder_v3(input_sce, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  input_sce <- doubletFinder_v3(input_sce, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  # table(input_sce@meta.data[,5])

  # plot
  input_sce@meta.data[,"DF_hi.lo"] <- input_sce@meta.data[,5]
  input_sce@meta.data$DF_hi.lo[which(input_sce@meta.data$DF_hi.lo == "Doublet" & input_sce@meta.data[,5] == "Singlet")] <- "Doublet-Low Confidience"
  input_sce@meta.data$DF_hi.lo[which(input_sce@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
  table(input_sce@meta.data$DF_hi.lo)
  
  p <- DimPlot(input_sce, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))+
    ggtitle(paste0('Doublet: ',table(input_sce@meta.data[,5])[1],'; ' ,'Singlet: ', table(input_sce@meta.data[,5])[2]))
  
  ggsave(filename=paste0('sample',x,'doubletFinder.pdf'),plot=p,height = 5.5,width = 6)

  sce.doublet.filter <- subset(input_sce, DF_hi.lo=='Singlet')
  
  return(sce.doublet.filter)

  })
  
}


