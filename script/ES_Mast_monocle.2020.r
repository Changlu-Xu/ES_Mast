library(monocle)



expr_matrix <- as.matrix(target_cells@assays$RNA@counts)
sample_sheet <- as.data.frame(target_cells@meta.data)
sample_sheet$seurat_clusters <- paste("Cluster",sample_sheet$seurat_clusters,sep = "")
head(sample_sheet)

gene_annotation <- as.data.frame(row.names(expr_matrix))
colnames(gene_annotation) <- "gene_short_name"
row.names(gene_annotation) <- gene_annotation$gene_short_name
head(gene_annotation)


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(expr_matrix,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 1))
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~samples + num_genes_expressed + S.Score +
                        G2M.Score + percent.mt",
                        verbose = T)

####get ording genes
HSMM_myo <- HSMM[expressed_genes,]
exprs_filtered <- t(t(exprs(HSMM_myo)/pData(HSMM_myo)$Size_Factor))

# Calculate the variance across genes without converting to a dense
# matrix:
expression_means <- Matrix::rowMeans(exprs_filtered)
expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
# Filter out genes that are constant across all cells:
genes_to_keep <- expression_vars > 0
exprs_filtered <- exprs_filtered[genes_to_keep,]
expression_means <- expression_means[genes_to_keep]
expression_vars <- expression_vars[genes_to_keep]
# Here's how to take the top PCA loading genes, but using
# sparseMatrix operations the whole time, using irlba. Note # that the v matrix from irlba is the loading matrix set.seed(0)
irlba_pca_res <- irlba(t(exprs_filtered),
                       nu=0,
                       center=expression_means,
                       scale=sqrt(expression_vars),
                       right_only=TRUE)$v
row.names(irlba_pca_res) <- row.names(exprs_filtered)
# Here, we will just
# take the top 200 genes from components 2 and 3.
# Component 1 usually is driven by technical noise.
# We could also use a more principled approach,
# similar to what dpFeature does below
PC2_genes <- names(sort(abs(irlba_pca_res[, 2]), decreasing = T))[1:200] 
PC3_genes <- names(sort(abs(irlba_pca_res[, 3]), decreasing = T))[1:200]
ordering_genes <- union(PC2_genes, PC3_genes)

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, method = 'DDRTree') 
HSMM_myo <- orderCells(HSMM_myo)

plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime",cell_size = 2)
plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM_myo, color_by = "define",cell_size = 2)


plot_cell_trajectory(HSMM_myo, color_by = "define", cell_size = 2) +
  facet_wrap(~define, nrow = 2)

plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 2)

saveRDS(HSMM_myo,"Desktop/BianGuohui_ES-Mast/ES_Mast_monocle.rds")

HSMM_myo <- readRDS("Desktop/BianGuohui_ES-Mast/ES_Mast_monocle.rds")


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$define)[,"IM"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])) }else { return (1) }
}
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))

plot_cell_trajectory(HSMM_myo, color_by="Pseudotime")

####MITF MAPK14 SRGN FCER1G STAT5B HRH2 CTSG LYZ MPO GYPA GYPB HBA1 HBA2
###HBD HBE1 HBG1 HBG2 EPOR
HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >= 1)) 
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("MITF", "MAPK14", "SRGN", 
                                                    "FCER1G", "STAT5B", "HRH2", 
                                                    "CTSG", "LYZ", "MPO")))
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("GYPA", "GYPB", "HBA1",
                                                    "HBA2","HBD", "HBE1",
                                                    "HBG1", "HBG2", "EPOR")))

cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="define")

plot_genes_branched_pseudotime(cds_subset,
                               branch_point=1,
                               color_by="define",
                               cell_size = 1,
                               ncol=3)

#####
CD34_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "CD34"))
KIT_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "KIT"))
SPN_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "SPN"))
PTPRC_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "PTPRC"))
KDR_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "KDR"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "CD34+KIT+", classify_func=function(x) {x[CD34_id,] > 0 & x[KIT_id,] > 0}) 
cth <- addCellType(cth, "CD34+CD45+", classify_func=function(x) {x[CD34_id,] > 0  & x[PTPRC_id,] > 0 }) 
cth <- addCellType(cth, "CD34+CD43+", classify_func=function(x) {x[CD34_id,] > 0  & x[SPN_id,] > 0 }) 
cth <- addCellType(cth, "CD34+KDR+", classify_func=function(x) {x[CD34_id,] > 0 & x[KDR_id,] > 0 }) 

cth <- addCellType(cth, "CD34+CD45+KIT+", classify_func=function(x) {x[CD34_id,] > 0 & x[KIT_id,] > 0 & x[PTPRC_id,] > 0 }) 
cth <- addCellType(cth, "CD34+CD45+KIT-", classify_func=function(x) {x[CD34_id,] > 0 & x[KIT_id,] == 0 & x[PTPRC_id,] > 0 }) 

HSMM_myo <- classifyCells(HSMM_myo, cth)

table(pData(HSMM_myo)$CellType)

p <- plot_cell_trajectory(HSMM_myo, color_by = "CellType",cell_size = 2)

p + scale_colour_manual(values = c("red", "gray"))
p + scale_colour_manual(values = c("blue", "gray"))
p + scale_colour_manual(values = c("green", "gray"))
p + scale_colour_manual(values = c("purple", "gray"))

p + scale_colour_manual(values = c("blue","red", "gray"))


BEAM_res <- BEAM(HSMM_myo, branch_point=1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(as.matrix(BEAM_res),"Desktop/BianGuohui_ES-Mast/Monocle/Branch_DEGs.csv")

BEAM_res <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/Pseudo_module.csv",
                     header = T,row.names = 1)

head(BEAM_res)
class(BEAM_res)

annotation_row <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/Pseudo_module_annotation.csv",
                                    header = T,row.names = 1)
trace(plot_genes_branched_heatmap,edit = T)
plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res, qval < 1e-7)),], 
                            branch_point = 1,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            add_annotation_row = annotation_row
                            )

heatmap_res <- plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res, qval < 1e-7)),], branch_point = 1,
                            num_clusters = 4,
                            cores = detectCores(),
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = T)
heatmap_res_list <- cutree( heatmap_res$ph_res$tree_row, k=4)
head(as.matrix(heatmap_res_list))
write.csv(as.matrix(heatmap_res_list),"Desktop/BianGuohui_ES-Mast/Monocle/monocle_cluster.csv")

heatmap_res_list <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/monocle_cluster.csv",header = T)
head(heatmap_res_list)
Branch_DEGs_list <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/Branch_DEGs.csv")
head(Branch_DEGs_list)

heatmap_res_list <- merge(heatmap_res_list,Branch_DEGs_list, by = "gene_name")
head(heatmap_res_list)
write.csv(heatmap_res_list,"Desktop/BianGuohui_ES-Mast/Monocle/monocle_cluster.csv")

######
heatmap_res_list <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/monocle_cluster.csv",header = T)
head(heatmap_res_list)

surface_list <- read.csv("Desktop/2018PNAS_2886_surface_marker.csv",header = T)
tf_list <- read.csv("Desktop/Homo_sapiens_TF.csv",header = T)

heatmap_res_sf_list <- merge(heatmap_res_list,surface_list,by = "gene_name",sort = F)
heatmap_res_tf_list <- merge(heatmap_res_list,tf_list,by = "gene_name",sort = F)

write.csv(heatmap_res_sf_list,"Desktop/BianGuohui_ES-Mast/Monocle/monocle_cluster_sf_list.csv")
write.csv(heatmap_res_tf_list,"Desktop/BianGuohui_ES-Mast/Monocle/monocle_cluster_tf_list.csv")

#######
top10_tf <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/cluster_tf_list_top10.csv",
                     header = T,row.names = 1)
head(top10_tf)

annotation_row <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/cluster_tf_list_top10_annotation.csv",
                           header = T,row.names = 1)
head(annotation_row)
class(annotation_row)

##
top10_tf <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/choosed_tf.csv",
                     header = T,row.names = 1)
head(top10_tf)

annotation_row <- read.csv("Desktop/BianGuohui_ES-Mast/Monocle/choosed_tf_annotation.csv",
                           header = T,row.names = 1)
head(annotation_row)




trace(plot_genes_branched_heatmap,edit = T) # false/delete
plot_genes_branched_heatmap(HSMM_myo[as.character(top10_tf$gene_short_name),], branch_point = 1,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            add_annotation_row = annotation_row
                            )
