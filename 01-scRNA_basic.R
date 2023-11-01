rm(list=ls());gc()
options(stringsAsFactors = F) 
source('Source code/00_lib.R')


library(future)
# check the current active plan
plan()
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()

# set the mem
# options(future.globals.maxSize = 10000 * 1024^2)


#sp <- 'human'
sp <- 'mouse'




my36colors <- c("#E31A1C", "#55c2fc", "#A6761D", "#F1E404", "#33A02C", "#1F78B4",
               "#FB9A99", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#F4B3BE",
               "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
               "#F4A11D", "#8DC8ED", "#4C6CB0", "#8A1C1B", "#CBCC2B", "#EA644C",
               "#634795", "#005B1D", "#26418A", "#CB8A93", "#B2DF8A", "#E22826",
               "#A6CEE3", "#F4D31D", "#F4A11D", "#82C800", "#8B5900", "#858ED1",
               "#FF72E1", "#CB50B2", "#007D9B", "#26418A", "#8B495F", "#FF394B")


#----------------------------------------------------------------------------------
#  Step 1: Load the Data
#----------------------------------------------------------------------------------


if(F){
## single sample ----

rt <- read.delim("count-matrix.txt",sep = '',header = T) %>% as.matrix()
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
sce.all <- CreateSeuratObject(counts = data,
                              project = "seurat", 
                              min.cells = 3,
                              min.features = 200)
}


## multiple samples ----

dir_name <- c('Ctrol','D166')
datalist <- list()
for (i in 1:length(dir_name)){
  dir.10x <- paste0("./data/Raw_data/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]] <- CreateSeuratObject(counts = my.data, 
                                   project = dir_name[i], 
                                   min.cells = 3, 
                                   min.features = 200)
}
names(datalist) <- dir_name

# avoid CheckDuplicateCellNames() warning
datalist <- pbapply::pblapply(1:length(datalist),FUN = function(x){
  data <- RenameCells(datalist[[x]],add.cell.id = dir_name[x])
})

sce.all.raw <- merge(datalist[[1]],y=datalist[2:length(datalist)])
as.data.frame(sce.all.raw@assays$RNA@counts[1:10, 1:2])
head(sce.all.raw@meta.data, 10) 
table(sce.all.raw@meta.data$orig.ident)
dim(sce.all.raw)




#----------------------------------------------------------------------------------
#  Step 1.1: doubletFinder
#----------------------------------------------------------------------------------

dir.create("./1-QC")
setwd("./1-QC")

source('Source code/01_DoubletFinder.R')

cdj_doublet_find(input_list = datalist)

names(datalist.doublet.fileter) <- names(datalist)

sce.all <- merge(datalist.doublet.fileter[[1]],y=datalist.doublet.fileter[2:length(datalist.doublet.fileter)])
as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10) 
table(sce.all@meta.data$orig.ident)
dim(sce.all.raw)
dim(sce.all)
# [1] 16787 14028
# [1] 16787 13234




#----------------------------------------------------------------------------------
#  Step 2: QC
#----------------------------------------------------------------------------------

source('Source code/02_QC.R')
sce.all.filt <- cdj_basic_qc(sce.all)
print(dim(sce.all.filt))
# [1] 16787 12792


setwd('../')




#----------------------------------------------------------------------------------
#  Step 3: Harmony
#----------------------------------------------------------------------------------

dir.create("./2-Harmony")
setwd("./2-Harmony")


source('Source code/03_Harmony.R')
sce.all.int <- cdj_run_harmony(sce.all.filt,batch = 'orig.ident')


setwd('../')




#----------------------------------------------------------------------------------
#  Step 4: Dimension reduction & Clustering
#----------------------------------------------------------------------------------


dir.create("./3-Annotation_checking")
setwd("./3-Annotation_checking")

sel.clust <- "RNA_snn_res.0.2" #checking the num of resolution
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) # n=12 





#----------------------------------------------------------------------------------
#  Step 5: Marker annotation & checking
#----------------------------------------------------------------------------------

#method <- 'umap'
method <- 'tsne'

######## MARKER_LIST ########


markers_list <- c(
  'pdac_all_markers_list',
  'pdac_macrophage_cell_sub'
  #'pdac_B_cell_sub',
  #'Bcels_markers_list'
)


source('Source code/04_Cell_type_annotation.R')


setwd('../')



# Define the celltype -----


dir.create("./4-Final_celltype")
setwd("./4-Final_celltype")


if(T){
sce.all.int$celltype <- sce.all.int$RNA_snn_res.0.2  # debug
sce.all.int$celltype <- recode(sce.all.int$celltype,
                               '0' = 'Ductal cell 3',
                               '1' = 'yT cell',
                               '2' = 'T cell',
                               '3' = 'Ductal cell 4',
                               '4' = 'Granulocytes',
                               '5' = 'Macrophages',
                               '6' = 'Acinar cell',
                               '7' = 'Fibroblasts',
                               '8' = 'B cell',
                               '9' = 'Ductal cell 2',
                               '10' = 'Ductal cell 1',
                               '11' = 'Endothelial',
                               '12' = 'DC')
Idents(sce.all.int) <- sce.all.int$celltype
sel.clust <-  "celltype"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
sce.all.int@meta.data$celltype <- factor(sce.all.int@meta.data$celltype,
                                         levels = c('Ductal cell 1','Ductal cell 2','Ductal cell 3','Ductal cell 4',
                                                    'T cell','yT cell','B cell','Granulocytes','Macrophages','Acinar cell','Fibroblasts','Endothelial','DC'))
}






#----------------------------------------------------------------------------------
#  Step 6: Cell fraction
#----------------------------------------------------------------------------------

dir.create("./5-Cell_fraction")
setwd("./5-Cell_fraction")

# calculate the propotion
phe <- sce.all.int@meta.data
head(phe)
table(phe$celltype,phe$orig.ident)
cal_table(phe$orig.ident,phe$celltype,prefix = 'celltype-vs-orig.ident')
#cal_table(phe$group,phe$celltype,prefix = 'celltype-vs-group')





setwd('../')








