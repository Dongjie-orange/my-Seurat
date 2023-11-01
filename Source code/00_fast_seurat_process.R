seurat_process <- function(pbmc,is_harmony,var=NULL){
  
  if(is_harmony == 'T'){
    
    
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 4000)
    pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
    pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))
    pbmc <- RunHarmony(pbmc, group.by.vars = var)
    pbmc <- FindNeighbors(pbmc, dims = 1:15, reduction = "harmony")
    pbmc <- FindClusters(pbmc,resolution = 0.8)
    pbmc <- RunUMAP(pbmc, dims = 1:15,reduction = 'harmony')
    pbmc <- RunTSNE(pbmc,dims = 1:20,reduction = 'harmony',check_duplicates = F)
    
  } else {
    
    
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 4000)
    pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
    pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))
    pbmc <- FindNeighbors(pbmc, dims = 1:15)
    pbmc <- FindClusters(pbmc, resolution = 0.8)
    pbmc <- RunUMAP(pbmc,dims = 1:15) 
    pbmc <- RunTSNE(pbmc,dims = 1:15,check_duplicates = F)
    
    
  }
}