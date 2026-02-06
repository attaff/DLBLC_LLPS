library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(harmony)
library(ggridges)
library(ggplot2)
library(ggstyle)
library(ggsci)
library(ggpubr) 
library(devtools)
library(scRNAtoolVis)
library(clusterProfiler)
plan("multicore", workers = 3)
options(future.globals.maxSize = 10000 * 1024^2)


DLBCL_1 <- Read10X(data.dir = "DLBCL1")
DLBCL_2 <- Read10X(data.dir = "DLBCL2")
DLBCL_3 <- Read10X(data.dir = "DLBCL3")
rLN_1 <- Read10X(data.dir = "rLN1")
rLN_2 <- Read10X(data.dir = "rLN2")
rLN_3 <- Read10X(data.dir = "rLN3")
DLBCL_1_object <- CreateSeuratObject(DLBCL_1, project = "DLBCL1", min.cells = 0, min.features = 200)
DLBCL_2_object <- CreateSeuratObject(DLBCL_2, project = "DLBCL2", min.cells = 0, min.features = 200)
DLBCL_3_object <- CreateSeuratObject(DLBCL_3, project = "DLBCL3", min.cells = 0, min.features = 200)
rLN_1_object <- CreateSeuratObject(rLN_1, project = "rLN1", min.cells = 0, min.features = 200)
rLN_2_object <- CreateSeuratObject(rLN_2, project = "rLN2", min.cells = 0, min.features = 200)
rLN_3_object <- CreateSeuratObject(rLN_3, project = "rLN3", min.cells = 0, min.features = 200)

DLBCL_1_object[["percent.mt"]] <- PercentageFeatureSet(DLBCL_1_object, pattern = "^MT-")
DLBCL_2_object[["percent.mt"]] <- PercentageFeatureSet(DLBCL_2_object, pattern = "^MT-")
DLBCL_3_object[["percent.mt"]] <- PercentageFeatureSet(DLBCL_3_object, pattern = "^MT-")
rLN_1_object[["percent.mt"]] <- PercentageFeatureSet(rLN_1_object, pattern = "^MT-")
rLN_2_object[["percent.mt"]] <- PercentageFeatureSet(rLN_2_object, pattern = "^MT-")
rLN_3_object[["percent.mt"]] <- PercentageFeatureSet(rLN_3_object, pattern = "^MT-")

DLBCL_1_object <- subset(DLBCL_1_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
DLBCL_2_object <- subset(DLBCL_2_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
DLBCL_3_object <- subset(DLBCL_3_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
rLN_1_object <- subset(rLN_1_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
rLN_2_object <- subset(rLN_2_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
rLN_3_object <- subset(rLN_3_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)

DLBCL_1_object <- RenameCells(DLBCL_1_object, add.cell.id = "DLBCL1")
DLBCL_2_object <- RenameCells(DLBCL_2_object, add.cell.id = "DLBCL2")
DLBCL_3_object <- RenameCells(DLBCL_3_object, add.cell.id = "DLBCL3")
rLN_1_object   <- RenameCells(rLN_1_object, add.cell.id = "rLN1")
rLN_2_object   <- RenameCells(rLN_2_object, add.cell.id = "rLN2")
rLN_3_object   <- RenameCells(rLN_3_object, add.cell.id = "rLN3")

object.list <- list(DLBCL_1_object, DLBCL_2_object, DLBCL_3_object,rLN_1_object,rLN_2_object,rLN_3_object)
object.list <- lapply(object.list, function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)       
  x <- RunPCA(x, npcs = 30)
  return(x)
})
anchors <- FindIntegrationAnchors(object.list = object.list, 
                                  dims = 1:30, 
                                  reduction = "rpca")
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)
DimPlot(integrated, reduction = "umap", group.by = "orig.ident", label = TRUE)

DefaultAssay(integrated) <- "RNA"
integrated <- JoinLayers(integrated)
markers <- FindAllMarkers(integrated, 
                          only.pos = TRUE, 
                          logfc.threshold = 1)
cells <- Cells(integrated)
B_annotation <- read.delim("AnnotationBcells.csv",
                           sep = ";")
T_annotation <- read.delim("Annotation_Tcells.csv",
                           sep = ";")
B_annotation <- B_annotation %>%
  mutate(cell_id = paste0(Sample, "_", Barcode, "-1"))
T_annotation <- T_annotation %>%
  mutate(cell_id = paste0(Sample, "_", Barcode, "-1"))

ann <- bind_rows(B_annotation, T_annotation) %>%
  select(cell_id, Population)

integrated$anno_manual <- NA
integrated$anno_manual[ann$cell_id] <- ann$Population
table(integrated$anno_manual, useNA = "ifany")

DimPlot(integrated, reduction = "umap", group.by = "anno_manual", label = TRUE, repel = TRUE)
DimPlot(integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)

meta <- integrated@meta.data
# 找出 NA 的 cell
na_cells <- rownames(meta)[is.na(meta$anno_manual)]

sce <- as.SingleCellExperiment(integrated)
# 参考数据库
ref <- HumanPrimaryCellAtlasData()
# 只对 NA 细胞跑 SingleR
pred <- SingleR(test = sce[, na_cells], ref = ref, labels = ref$label.main)
# 写回 meta.data
integrated$anno_manual[na_cells] <- pred$labels
DimPlot(integrated, reduction = "umap", group.by = "anno_manual", label = TRUE, repel = TRUE)

unique(integrated$anno_manual) 

integrated$anno_manual <- ifelse(
  integrated$anno_manual %in% c("B_cell", "Pro-B_cell_CD34+","Pre-B_cell_CD34-","B_cell"),
  "MalignantB",
  integrated$anno_manual
)

integrated$anno_manual <- ifelse(
  integrated$anno_manual %in% c("T_cells"),
  "Uncluster",
  integrated$anno_manual
)

integrated$anno_manual <- ifelse(
  integrated$anno_manual %in% c("Macrophage","Myeloid","DC","GMP","Monocyte"),
  "Myeloid_cell",
  integrated$anno_manual
)

unique(integrated$anno_manual) 
DimPlot(integrated, reduction = "umap", group.by = "anno_manual", label = TRUE, repel = TRUE)

llps_model <- readRDS("llps_model.RData")
coefficients <- llps_model$coefficients
coefficients
featurePlot(object = integrated, genes = names(coefficients), nrow = 2, ncol = 3)
# 获取RNA表达矩阵（log-normalized）
gene_expression <- GetAssayData(integrated, assay = "RNA", slot = "data")
# 取交集
common_genes <- intersect(rownames(gene_expression), names(coefficients))
gene_expression <- gene_expression[common_genes, ]
coefficients <- coefficients[common_genes]
coefficients

gene_expression_t <- as.data.frame(t(gene_expression)) 

cell_sums <- rowSums(gene_expression_t)
valid_cells <- cell_sums > 0
gene_expression_t <- gene_expression_t[valid_cells, ]

for (i in seq_along(names(coefficients))) {
  gene_expression_t[, i] <- gene_expression_t[, i] * unname(coefficients)[i]
}
gene_expression_t$risk_score <- rowSums(gene_expression_t[,1:6])

###########################################################################
plot_df <- data.frame(
  celltype = integrated$anno_manual[rownames(gene_expression_t)],
  risk_score = gene_expression_t$risk_score
)
celltype_order <- c("MalignantB", "HealthyB", 
                    "TREG", "TFH","TTOX", "TH","NK_cell", 
                    "Myeloid_cell","Uncluster")
plot_df$celltype <- factor(plot_df$celltype, levels = rev(celltype_order))

malignant_median <- median(plot_df$risk_score[plot_df$celltype == "MalignantB"], na.rm = TRUE)
malignant_median
#-0.2482899

ggplot(plot_df, aes(x = risk_score, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 2, alpha = 0.6, rel_min_height = 0.01) +
  geom_vline(xintercept = malignant_median, linetype = "dashed", color = "black") +
  annotate("text", x = malignant_median, y = length(levels(plot_df$celltype)) + 0.3,
           label = paste0("Median: ", round(malignant_median, 2)),
           color = "black", angle = 0, vjust = -0.5, hjust = 0) +
  scale_fill_sci(
    palette = "observable10",
    type = "discrete",
    modeColor = 1
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10,color = "black"),
    axis.title.x = element_text(size = 12,color = "black"),
    axis.text.y = element_text(size = 10,color = "black")
  ) +
  labs(x = "Risk score", y = "Cell type", title = "Risk score distribution across cell types")



p <- ggplot(plot_df, aes(x = risk_score, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 1, alpha = 0.6, rel_min_height = 0.01, color = NA) +

  geom_boxplot(aes(group = celltype),
               width = 0.15, 
               outlier.shape = NA, 
               alpha = 0.8,
               position = position_nudge(y = -0.25)) +
  geom_vline(xintercept = malignant_median, linetype = "dashed", color = "black") +
  annotate("text", x = malignant_median, y = length(levels(plot_df$celltype)) + 0.3,
           label = paste0("Median: ", round(malignant_median, 2)),
           color = "black", vjust = -0.5, hjust = 0) +
  scale_fill_sci(
    palette = "observable10",
    type = "discrete",
    modeColor = 1
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10,color = "black"),
    axis.title.x = element_text(size = 12,color = "black"),
    axis.text.y = element_text(size = 10,color = "black")
  ) +
  labs(x = "Risk score", y = "Cell type", title = "Risk score distribution across cell types")
p

#############################################################################
#######################################################################################
valid_integrated <- integrated
table(valid_integrated$risk_group)
DimPlot(valid_integrated, reduction = "umap", group.by = "risk_group", label = TRUE, repel = TRUE,cols = c("High" = "firebrick", "Low" = "steelblue"))
DimPlot(valid_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(valid_integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
DimPlot(valid_integrated, reduction = "umap", group.by = "anno_manual", label = TRUE, repel = TRUE)

library(viridis)
FeaturePlot(
  valid_integrated,
  features = "risk_score",
  reduction = "umap"
) + scale_color_viridis(option = "viridis")

cell_types <- setdiff(unique(valid_integrated@meta.data$anno_manual), "Uncluster")
diff_results_list <- list()
for (ct in cell_types) {
  cells_ct <- WhichCells(valid_integrated, expression = anno_manual == ct)
  # 子集对象
  sub_obj <- subset(valid_integrated, cells = cells_ct)
  # FindMarkers 可以做高低危差异
  diff_res <- FindMarkers(
    object = sub_obj,
    ident.1 = "High",
    ident.2 = "Low",
    group.by = "risk_group",
    logfc.threshold = 1,
    min.pct = 0.1,
    test.use = "MAST"
  )
  diff_results_list[[ct]] <- diff_res
}

all_markers <- bind_rows(
  lapply(names(diff_results_list), function(ct) {
    df <- diff_results_list[[ct]]
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df$gene <- rownames(df)  
    df$cluster <- ct

    df <- df %>%
      dplyr::select(
        p_val,
        avg_log2FC,
        pct.1,
        pct.2,
        p_val_adj,
        cluster,
        gene
      )
    return(df)
  })
)

all_markers_test <- all_markers[!all_markers$gene %in% common_genes,]
jjVolcano(
  diffData = all_markers_test,
  log2FC.cutoff = 1,
  topGeneN = 5,
  pSize = 0.7,
  size  = 2.5,
  tile.col = corrplot::COL2('RdBu', 15)[4:12],
  fontface = 'italic',
  aesCol = c("#0070b8","#fbb01a"),

)

MalignB_markers <- all_markers[all_markers$cluster == "MalignantB",]
MalignB_markers <- MalignB_markers[order(MalignB_markers$avg_log2FC, decreasing =TRUE),]

MalignB_markers$regulation <- "normal"
MalignB_markers$regulation[MalignB_markers$avg_log2FC > 0.585 & MalignB_markers$p_val_adj < 0.05] <-"up"
MalignB_markers$regulation[MalignB_markers$avg_log2FC < -0.585 & MalignB_markers$p_val_adj < 0.05] <-"down"

genelist <- MalignB_markers$avg_log2FC
names(genelist) <- rownames(MalignB_markers)

library(msigdbr)
library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(GseaVis)
gmt_file <- "h.all.v2025.1.Hs.symbols.gmt"
hallmark_gene_sets <- read.gmt(gmt_file)
hallmark_gene_sets$term <- gsub(pattern ="HALLMARK_","", hallmark_gene_sets$term)
BiocParallel::register(BiocParallel::SerialParam())
egmt <- GSEA(genelist,
             TERM2GENE = hallmark_gene_sets,
             pvalueCutoff = 0.05,
             minGSSize = 1,
             maxGSSize = 1000,
             pAdjustMethod = "fdr") 

egmt@result$ID
gseaRes <- egmt
term <- c("MTORC1_SIGNALING","MYC_TARGETS_V1","CHOLESTEROL_HOMEOSTASIS",
          "MITOTIC_SPINDLE","INTERFERON_GAMMA_RESPONSE","PI3K_AKT_MTOR_SIGNALING")

gseaNb(object = egmt,
       geneSetID = term,
       addPval = T,
       pvalX = 0.8,
       pvalY = 0.5,
       addGene = T,
       markTopgene = T,
       geneCol = "black",
       subPlot = 2,
       curveCol = jjAnno::useMyCol('stallion2',6),
       pDigit = 3,
       
)
########################################################
library(CellChat)
library(patchwork)
library(future)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

data.input <- GetAssayData(valid_integrated, slot = "data")
meta <- valid_integrated@meta.data[, c("orig.ident", "anno_manual","risk_group")]
colnames(meta) <- c("group", "labels","risk_group")

celltype_order <- c("MalignantB","HealthyB","TTOX","NK_cell",
                    "TFH","TH","TREG","Myeloid_cell","Uncluster")

meta$labels <- factor(meta$labels, levels = celltype_order)
ordered_indices <- order(meta$labels)
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
stopifnot(identical(rownames(meta), colnames(data.input)))

meta$labels_grouped <- with(meta,
                            ifelse(labels == "MalignantB" & risk_group == "High", "MalignantB_High",
                                   ifelse(labels == "MalignantB" & risk_group == "Low",  "MalignantB_Low",
                                          ifelse(labels == "Uncluster", NA, as.character(labels))))
                            )



keep_cells <- !is.na(meta$labels_grouped)
meta <- meta[keep_cells, ]
data.input <- data.input[, keep_cells]
meta$labels_grouped <- factor(meta$labels_grouped, levels = c(celltype_order, "MalignantB_High", "MalignantB_Low"))
stopifnot(identical(rownames(meta), colnames(data.input)))
meta$labels_grouped <- droplevels(meta$labels_grouped)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels_grouped")
cellchat@DB <- subsetDB(CellChatDB.human, search = c("Secreted Signaling","Cell-Cell Contact"))

cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
dim(cellchat@net$count)
rownames(cellchat@net$count)

groupSize <- as.numeric(table(cellchat@idents)) 

netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction weights/strength")

cellchat@netP$pathways

pheatmap::pheatmap(cellchat@net$counts, border_color = "black",
                   cluster_cols = FALSE, fontsize = 10, cluster_rows = FALSE,
                   display_numbers = TRUE, number_color="black", number_format = "%.0f",
                   scale = "none")
pheatmap::pheatmap(cellchat@net$weight, border_color = "black",
                   cluster_cols = FALSE, fontsize = 10, cluster_rows = FALSE,
                   display_numbers = TRUE, number_color="black", number_format = "%.2f",
                   scale = "none")
netVisual_heatmap(cellchat, signaling = "CD22", color.heatmap ="Reds")
plotGeneExpression(cellchat, 
                   signaling = "CD22",
                   enriched.only = FALSE)
netVisual_bubble(cellchat, remove.isolate = FALSE)
