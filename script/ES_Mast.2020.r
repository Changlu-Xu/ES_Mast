
library(Seurat)
express.data <- read.csv("Desktop/BianGuohui_ES-Mast/ES_Mast.expression.csv",header = T,row.names = 1)
meta.data <- read.csv("Desktop/BianGuohui_ES-Mast/ES_Mast.meta.csv",header = T, row.names = 1)
embedding.data <- read.csv("Desktop/BianGuohui_ES-Mast/ES_Mast.cell.embeddings.csv",header = T, row.names = 1)

ES_Mast <- CreateSeuratObject(counts = express.data, project = "ES_Mast" )
VlnPlot(ES_Mast, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

FeatureScatter(ES_Mast, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ES_Mast <- NormalizeData(ES_Mast, normalization.method = "LogNormalize", scale.factor = 10000)
ES_Mast <- FindVariableFeatures(ES_Mast, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(ES_Mast)
ES_Mast <- ScaleData(ES_Mast, features = all.genes)

ES_Mast <- RunPCA(ES_Mast, features = VariableFeatures(object = ES_Mast))

DimPlot(ES_Mast, reduction = "pca")

ES_Mast <- RunUMAP(ES_Mast, dims = 1:10)

DimPlot(ES_Mast, reduction = "umap")

head(ES_Mast@reductions$umap@cell.embeddings)
ES_Mast@reductions$umap@cell.embeddings <- as.matrix(embedding.data)
head(ES_Mast@reductions$umap@cell.embeddings)

DimPlot(ES_Mast, reduction = "umap")

Idents(ES_Mast) <- factor(meta.data$seurat_clusters,ordered = T)

FeaturePlot(ES_Mast, features = c("CD34","KDR","KIT"),cols = c("lightgrey", "red"), split.by = "samples")

ES_Mast@meta.data <- meta.data

DimPlot(ES_Mast, reduction = "umap", group.by = "samples")

saveRDS(ES_Mast,"Desktop/BianGuohui_ES-Mast/ES_Mast_2020.rds")

ES_Mast <- readRDS("Desktop/BianGuohui_ES-Mast/ES_Mast_2020.rds")

DefaultAssay(ES_Mast) <- "RNA"

all.markers <- FindAllMarkers(object = ES_Mast, only.pos = T)
write.csv(as.matrix(all.markers),"Desktop/BianGuohui_ES-Mast/All_DEGs.csv")

new.cluster.ids <- c("others", "others", "others", "others", "others", "others", "others",
                     "others", "others", "others", "others", "others","others")
names(new.cluster.ids) <- levels(ES_Mast)
ES_Mast <- RenameIdents(ES_Mast, new.cluster.ids)
DimPlot(ES_Mast, reduction = "umap")

KDR <- subset(ES_Mast, subset =  KDR > 0 & CD34 == 0 &  KIT == 0)
KDR_cells <- row.names(KDR@meta.data)
ES_Mast <- SetIdent(object = ES_Mast, cells = KDR_cells, value = 'KDR+')


CD34 <- subset(ES_Mast, subset =  KDR == 0 & CD34 > 0 &  KIT == 0)
CD34_cells <- row.names(CD34@meta.data)
ES_Mast <- SetIdent(object = ES_Mast, cells = CD34_cells, value = 'CD34+')

KIT <- subset(ES_Mast, subset =  KDR == 0 & CD34 == 0 &  KIT > 0)
KIT_cells <- row.names(KIT@meta.data)
ES_Mast <- SetIdent(object = ES_Mast, cells = KIT_cells, value = 'KIT+')

KDRCD34 <- subset(ES_Mast, subset =  KDR > 0 & CD34 > 0 &  KIT == 0)
KDRCD34_cells <- row.names(KDRCD34@meta.data)
ES_Mast <- SetIdent(object = ES_Mast, cells = KDRCD34_cells, value = 'KDR+CD34+')

CD34KIT <- subset(ES_Mast, subset =  KDR == 0 & CD34 > 0 &  KIT > 0)
CD34KIT_cells <- row.names(CD34KIT@meta.data)
ES_Mast <- SetIdent(object = ES_Mast, cells = CD34KIT_cells, value = 'CD34+KIT+')

CD34KITKDR <- subset(ES_Mast, subset =  KDR > 0 & CD34 > 0 &  KIT > 0)
CD34KITKDR_cells <- row.names(CD34KITKDR@meta.data)
ES_Mast <- SetIdent(object = ES_Mast, cells = CD34KITKDR_cells, value = 'KDR+CD34+KIT+')

DimPlot(ES_Mast, reduction = "umap")
DimPlot(ES_Mast, reduction = "umap", cols= c("black","yellow","purple","red","blue","green","gray"))


target_cells <- subset(x = ES_Mast, idents = c(0,1,4,5,7))
DimPlot(target_cells, reduction = "umap")

target_cells <- NormalizeData(target_cells, normalization.method = "LogNormalize", scale.factor = 10000)
target_cells <- FindVariableFeatures(target_cells, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(target_cells)
target_cells <- ScaleData(target_cells, features = all.genes)

target_cells <- RunPCA(target_cells, features = VariableFeatures(object = ES_Mast))

DimPlot(target_cells, reduction = "pca", pt.size = 1)
DimPlot(target_cells, reduction = "pca", pt.size = 1, group.by = "samples")

Idents(target_cells) <- paste(target_cells@meta.data$seurat_clusters,
                              target_cells@meta.data$samples, sep = "_")
DimPlot(target_cells, reduction = "pca", pt.size = 1)

new.cluster.ids <- c("IM","MPP-pre-MC","EMP","EMP-pre-MC","MPP-pre-MC",
                     "Myeloid","HE")
names(new.cluster.ids) <- levels(target_cells)
target_cells <- RenameIdents(target_cells, new.cluster.ids)
DimPlot(target_cells, reduction = "pca", pt.size = 1)

target_cells@meta.data$defined <- Idents(target_cells)
head(target_cells@meta.data)
levels(x = target_cells) <- c("IM","HE","EMP","EMP-pre-MC","MPP-pre-MC",
                              "Myeloid")
DimPlot(target_cells, reduction = "pca", pt.size = 1)

saveRDS(target_cells,"Desktop/BianGuohui_ES-Mast/ES_Mast_Target_cells.rds")
target_cells <- readRDS("Desktop/BianGuohui_ES-Mast/ES_Mast_Target_cells_DP.rds")

vln_genes <- read.csv("Desktop/BianGuohui_ES-Mast/genes_vln.csv", header = T)

VlnPlot(object = target_cells, features = as.character(vln_genes$genes[1:6]),
         pt.size = 0,ncol = 3)

VlnPlot(object = target_cells, features = as.character(vln_genes$genes[7:12]),
        pt.size = 0,ncol = 3)

target_cells_dp <- subset( x= target_cells, idents = c("IM_DP","HE_DP","EMP_DP",
                                                       "EMP-pre-MC_DP",
                                                       "MPP-pre-MC_DP",
                                                       "Myeloid_DP"))
target_cells_ng <- subset( x= target_cells, idents = c("IM","HE","EMP",
                                                       "EMP-pre-MC",
                                                       "MPP-pre-MC",
                                                       "Myeloid"))


VlnPlot(object = target_cells_dp, features = c("MAPK14","FCER1G","STAT5B"),
        pt.size = 0,ncol = 3)


####
IM <- subset( x = target_cells, idents = "Myeloid")
IM_DP <- subset(IM, subset = CD34 > 0 &  KIT > 0)
IM_DP_cells <- row.names(IM_DP@meta.data)
target_cells <- SetIdent(object = target_cells, cells = IM_DP_cells, 
                         value = 'Myeloid_DP')

levels(x = target_cells) <- c("IM","IM_DP","HE","HE_DP","EMP","EMP_DP",
                              "EMP-pre-MC","EMP-pre-MC_DP",
                              "MPP-pre-MC","MPP-pre-MC_DP",
                              "Myeloid","Myeloid_DP")
DimPlot(target_cells, reduction = "pca", pt.size = 1)

saveRDS(target_cells,"Desktop/BianGuohui_ES-Mast/ES_Mast_Target_cells_DP.rds")
target_cells <- readRDS("Desktop/BianGuohui_ES-Mast/ES_Mast_Target_cells_DP.rds")

target_cells2 <- subset(x=target_cells,idents = c("EMP-pre-MC","EMP-pre-MC_DP",
                                                  "MPP-pre-MC","MPP-pre-MC_DP",
                                                  "Myeloid","Myeloid_DP"))

Average_Expression <- AverageExpression(object = target_cells)
Average_Expression <- Average_Expression$RNA

Myeloid_genes <- read.csv("Desktop/BianGuohui_ES-Mast/Myeloid.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(Myeloid_genes$genes),
          disp.min = -2, disp.max = 2)

Myeloid_genes_average <- Average_Expression[as.character(Myeloid_genes$genes),]

pheatmap(Myeloid_genes_average,scale = "row",cluster_cols  = F)


EMP_pre_MC_genes <- read.csv("Desktop/BianGuohui_ES-Mast/EMP_pre_MC.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(EMP_pre_MC_genes$genes),
          disp.min = -2, disp.max = 2)

EMP_pre_MC_genes_average <- Average_Expression[as.character(EMP_pre_MC_genes$genes),]

pheatmap(EMP_pre_MC_genes_average,scale = "row",cluster_cols  = F)


MPP_pre_MC_genes <- read.csv("Desktop/BianGuohui_ES-Mast/MPP_pre_MC.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(MPP_pre_MC_genes$genes),
          disp.min = -2, disp.max = 2)

MPP_pre_MC_genes_average <- Average_Expression[as.character(MPP_pre_MC_genes$genes),]

pheatmap(MPP_pre_MC_genes_average,scale = "row",cluster_cols  = F)

DE_MC_genes <- read.csv("Desktop/BianGuohui_ES-Mast/MC-GENE_13.csv", header = T)
DoHeatmap(object = target_cells2, features = as.character(DE_MC_genes$DE_MCs),
          disp.min = -1, disp.max = 1)


VlnPlot(object = target_cells2, features = as.character(DE_MC_genes$DE_MCs),
          pt.size = 0.1,ncol = 4)


target_cells_dp <- subset(target_cells, subset = CD34 > 0 &  KIT > 0)
DimPlot(target_cells_dp, reduction = "pca", pt.size = 1)


IM_EG_genes <- read.csv("Desktop/BianGuohui_ES-Mast/IM_EG.csv", header = T)
head(IM_EG_genes)
DoHeatmap(object = target_cells, features = as.character(IM_EG_genes$IM_12),
          disp.min = -1.5, disp.max = 1.5)

EG_genes <- read.csv("Desktop/BianGuohui_ES-Mast/EG.csv", header = T)
head(EG_genes)
DoHeatmap(object = target_cells, features = as.character(EG_genes$EG),
          disp.min = -1, disp.max = 1)


HES_genes <- read.csv("Desktop/BianGuohui_ES-Mast/HES.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(HES_genes$HES),disp.min = -2, disp.max = 2)

DE_MCs_genes <- read.csv("Desktop/BianGuohui_ES-Mast/MC-GENE.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(DE_MCs_genes$DE_MCs),
          disp.min = -1.5, disp.max = 1.5)

Fun_MCs_genes <- read.csv("Desktop/BianGuohui_ES-Mast/Fun_MCs.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(Fun_MCs_genes$Fun_MCs),
          disp.min = -1, disp.max = 1)

Active_MCs_genes <- read.csv("Desktop/BianGuohui_ES-Mast/Activation_MCs.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(Active_MCs_genes$Activation_MCs),
          disp.min = -1.5, disp.max = 1.5)

MCs_merge_genes <- read.csv("Desktop/BianGuohui_ES-Mast/MC_Merge.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(MCs_merge_genes$Activation_MCs),
          disp.min = -1.5, disp.max = 1.5)



My_EMP_genes <- read.csv("Desktop/BianGuohui_ES-Mast/My_EMP.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(My_EMP_genes$Myeloid.progenitor),
          disp.min = -1.5, disp.max = 1.5)

Myeloid_genes <- read.csv("Desktop/BianGuohui_ES-Mast/Myeloid_progenitor.csv", header = T)
DoHeatmap(object = target_cells, features = as.character(Myeloid_genes$Myeloid.progenitor),
          disp.min = -1.5, disp.max = 1.5)




target_cells <- RunUMAP(target_cells, dims = 1:10)

DimPlot(target_cells, reduction = "umap")


target_cells <- FindNeighbors(target_cells, dims = 1:10)
target_cells <- FindClusters(target_cells, resolution = 0.2)
DimPlot(target_cells, reduction = "umap")

cluster0 <- FindNeighbors(cluster0, dims = 1:10)
cluster0 <- FindClusters(cluster0, resolution = 0.1)
DimPlot(cluster0, reduction = "umap",label = T)

all_markers <- FindAllMarkers( object = cluster0,only.pos = T)


########choose specific cluster#####
cluster0_0 <- subset(x = cluster0, idents = 0)
new.cluster.ids <- c("others")
names(new.cluster.ids) <- levels(cluster0_0)
cluster0_0 <- RenameIdents(cluster0_0, new.cluster.ids)
DimPlot(cluster0_0, reduction = "umap")



CD34KIT <- subset(cluster0_0, subset =  CD34 > 0 &  KIT > 0)
CD34KIT_cells <- row.names(CD34KIT@meta.data)
cluster0 <- SetIdent(object = cluster0, cells = CD34KIT_cells, value = '0_CD34+KIT+')

cluster0_dp <- subset( x = cluster0, idents = c("0_CD34+KIT+","1_CD34+KIT+"))

DimPlot(cluster0, reduction = "umap", pt.size = 1, cols = c("red","gray"))

DoHeatmap(object = cluster0, features = as.character(IM_genes$IM_12),
          disp.min = -1, disp.max = 1)

DoHeatmap(object = cluster0, features = as.character(EG_genes$EG),
          disp.min = -1, disp.max = 1)

DoHeatmap(object = cluster0, features = as.character(HES_genes$HES),disp.min = -2, disp.max = 2)

DoHeatmap(object = cluster0, features = as.character(DE_MCs_genes$DE_MCs),
          disp.min = -1, disp.max = 1)

DoHeatmap(object = cluster0, features = as.character(Fun_MCs_genes$Fun_MCs),
          disp.min = -1, disp.max = 1)

DoHeatmap(object = cluster0, features = as.character(Active_MCs_genes$Activation_MCs),
          disp.min = -1.5, disp.max = 1.5)

DoHeatmap(object = cluster0, features = as.character(EMP_genes$EMP),
          disp.min = -1.5, disp.max = 1.5)

DoHeatmap(object = cluster0, features = as.character(Myeloid_genes$Myeloid.progenitor),
          disp.min = -1.5, disp.max = 1.5)





#####cluster1/4/5
cluster4 <- subset(x = ES_Mast, idents = 2)
new.cluster.ids <- c("others")
names(new.cluster.ids) <- levels(cluster4)
cluster4 <- RenameIdents(cluster4, new.cluster.ids)
DimPlot(cluster4, reduction = "umap")

CD34KITCD45 <- subset(cluster4, subset =  PTPRC > 0 & CD34 > 0 &  KIT > 0)
CD34KITCD45_cells <- row.names(CD34KITCD45@meta.data)
cluster4 <- SetIdent(object = cluster4, cells = CD34KITCD45_cells, value = 'CD34+KIT+CD45+')

CD34KITCD45neg <- subset(cluster4, subset =  PTPRC == 0 & CD34 > 0 &  KIT > 0)
CD34KITCD45neg_cells <- row.names(CD34KITCD45neg@meta.data)
cluster4 <- SetIdent(object = cluster4, cells = CD34KITCD45neg_cells, value = 'CD34+KIT+CD45-')

DimPlot(cluster4, reduction = "umap", pt.size = 1, cols = c("blue","red","gray"))

plot_cluster <- function(cluster){
  cluster_id <- subset(x = ES_Mast, idents = cluster)
  new.cluster.ids <- c("others")
  names(new.cluster.ids) <- levels(cluster_id)
  cluster_id <- RenameIdents(cluster_id, new.cluster.ids)
  
  CD34KITCD45 <- subset(cluster_id, subset =  PTPRC > 0 & CD34 > 0 &  KIT > 0)
  CD34KITCD45_cells <- row.names(CD34KITCD45@meta.data)
  cluster_id <- SetIdent(object = cluster_id, cells = CD34KITCD45_cells, value = 'CD34+KIT+CD45+')
  
  CD34KITCD45neg <- subset(cluster_id, subset =  PTPRC == 0 & CD34 > 0 &  KIT > 0)
  CD34KITCD45neg_cells <- row.names(CD34KITCD45neg@meta.data)
  cluster_id <- SetIdent(object = cluster_id, cells = CD34KITCD45neg_cells, value = 'CD34+KIT+CD45-')

  
}

plot_cluster(0)
DimPlot(cluster_id, reduction = "umap", pt.size = 1, cols = c("blue","red","gray"))

plot_cluster(1)
DimPlot(cluster_id, reduction = "umap", pt.size = 1, cols = c("blue","red","gray"))

plot_cluster(0)
DimPlot(cluster_id, reduction = "umap", pt.size = 1, cols = c("blue","red","gray"))

plot_cluster(0)
DimPlot(cluster_id, reduction = "umap", pt.size = 1, cols = c("blue","red","gray"))


######
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

cluster_genes <- read.csv("Desktop/BianGuohui_ES-Mast/Cluster_DEGs.csv",
                          header = T)
head(cluster_genes)

go_result <- compareCluster(ENTREZID~cluster,data =cluster_genes ,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db", 
                     ont= "BP")
dotplot(go_result, showCategory=5, includeAll=FALSE)

write.csv(as.data.frame(go_result),"Desktop/BianGuohui_ES-Mast/Cluster_DEGs_GO-result.csv")

kegg_result <- compareCluster(ENTREZID~cluster,data =cluster_genes ,
                            fun="enrichKEGG",
                            organism="hsa")
dotplot(kegg_result, showCategory=5, includeAll=FALSE)
write.csv(as.data.frame(kegg_result),"Desktop/BianGuohui_ES-Mast/Cluster_DEGs_KEGG-result.csv")
