
cdj_run_harmony <- function(input_sce,batch){
  
  # run harmony ----
  input_sce <- NormalizeData(input_sce, 
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4) 
  input_sce <- FindVariableFeatures(input_sce)
  input_sce <- ScaleData(input_sce)
  input_sce <- RunPCA(input_sce, features = VariableFeatures(object = input_sce))
  seuratObj <- RunHarmony(input_sce, group.by.vars = batch)
  seuratObj <- RunTSNE(seuratObj,  dims = 1:15, 
                       reduction = "harmony")
  seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                       reduction = "harmony")
  input_sce <- seuratObj
  input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                             dims = 1:15) 
  
  # check resolution ----
  input_sce.all <- input_sce
  for (res in seq(0,1,0.1)) {
    input_sce.all <- FindClusters(input_sce.all,
                               resolution = res, algorithm = 1)
  }
  colnames(input_sce.all@meta.data)
  apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)
  
  # umap plot ----- 
  p1_dim <- plot_grid(ncol = 5, nrow = 2,
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T) + ggtitle("louvain_0.1"), 
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2",label = T) + ggtitle("louvain_0.2"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3",label = T) + ggtitle("louvain_0.3"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.4",label = T) + ggtitle("louvain_0.4"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.5",label = T) + ggtitle("louvain_0.5"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.6",label = T) + ggtitle("louvain_0.6"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.7",label = T) + ggtitle("louvain_0.7"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8",label = T) + ggtitle("louvain_0.8"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.9",label = T) + ggtitle("louvain_0.9"),
                      DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.1",label = T) + ggtitle("louvain_1"))
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_umap.pdf",width = 30,height = 10)
  
  # tsne plot ----- 
  p1_dim <- plot_grid(ncol = 5, nrow = 2,
                   DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.1",label = T) + ggtitle("louvain_0.1"), 
                   DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.2",label = T) + ggtitle("louvain_0.2"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.3",label = T) + ggtitle("louvain_0.3"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.4",label = T) + ggtitle("louvain_0.4"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.5",label = T) + ggtitle("louvain_0.5"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.6",label = T) + ggtitle("louvain_0.6"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.7",label = T) + ggtitle("louvain_0.7"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.8",label = T) + ggtitle("louvain_0.8"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.9",label = T) + ggtitle("louvain_0.9"),
    DimPlot(input_sce.all, reduction = "tsne", group.by = "RNA_snn_res.1",label = T) + ggtitle("louvain_1"))
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_tsne.pdf",width = 30,height = 10)
  
  # tree plot ----- 
  p2_tree <- clustree(input_sce.all@meta.data, prefix = "RNA_snn_res.")
  ggsave(plot=p2_tree, filename="Tree_diff_resolution_tsne.pdf",height = 20,width = 10)

  saveRDS(input_sce.all, "sce.all_int.rds")
  
  return(input_sce.all)
 
  
  
  
}