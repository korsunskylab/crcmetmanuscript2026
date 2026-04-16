# 加载 Seurat 和 Matrix 库
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("/projects/b1198/epifluidlab/zhengyue/spring")

# 读取数据
counts <- readMM("GSE243013_NSCLC_immune_scRNA_counts.mtx")
barcodes <- read.csv("GSE243013_barcodes.csv", header = TRUE)
genes <- read.csv("GSE243013_genes.csv", header = TRUE)
metadata <- read.csv("GSE243013_NSCLC_immune_scRNA_metadata.csv", header = TRUE)
counts <- as(counts, "dgCMatrix")
counts <- t(counts)  # Transpose to switch genes and cells
rownames(counts) <- genes$geneSymbol  # 设置基因名
colnames(counts) <- barcodes$barcode  # 设置细胞条形码
NSCLC.object <- CreateSeuratObject(counts = counts, 
                                 meta.data = metadata, # Add metadata
                                 project = "NSCLC.Zhang.2025")
NSCLC.object[["percent.mt"]] <- PercentageFeatureSet(NSCLC.object, pattern = "^MT-")  # Calculate percentage of mitochondrial genes
NSCLC.object[["percent.rb"]] <- PercentageFeatureSet(NSCLC.object, pattern = "^RPL|^RPS")  # Calculate ribosomal genes
VlnPlot(NSCLC.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)  # Violin plots for basic metrics
NSCLC.object <- subset(NSCLC.object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
NSCLC.object <- NormalizeData(NSCLC.object, normalization.method = "LogNormalize", scale.factor = 10000)
NSCLC.object <- FindVariableFeatures(NSCLC.object, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(NSCLC.object), 10)
VariableFeaturePlot(NSCLC.object)
NSCLC.object <- ScaleData(NSCLC.object, features = rownames(NSCLC.object))
NSCLC.object <- RunPCA(NSCLC.object, features = VariableFeatures(object = NSCLC.object))
ElbowPlot(NSCLC.object)  # Look for the "elbow" to determine the number of PCs to use
NSCLC.object <- FindNeighbors(NSCLC.object, dims = 1:10)  # Using the first 10 principal components
NSCLC.object <- FindClusters(NSCLC.object, resolution = 0.5)  # Adjust resolution to change number of clusters
NSCLC.object <- RunUMAP(NSCLC.object, dims = 1:10)
DimPlot(NSCLC.object, reduction = "umap", group.by = "seurat_clusters")

# Visualize cells by major cell type using UMAP
DimPlot(NSCLC.object, reduction = "umap", group.by = "major_cell_type")


# subset T/NK cells
T_NK_cells <- subset(NSCLC.object, subset = major_cell_type == "T/NK cell")
T_NK_cells <- NormalizeData(T_NK_cells, normalization.method = "LogNormalize", scale.factor = 10000)
T_NK_cells <- FindVariableFeatures(T_NK_cells, selection.method = "vst", nfeatures = 2000)
T_NK_cells <- ScaleData(T_NK_cells, features = rownames(T_NK_cells))
T_NK_cells <- RunPCA(T_NK_cells, features = VariableFeatures(object = T_NK_cells))
ElbowPlot(T_NK_cells)
T_NK_cells <- FindNeighbors(T_NK_cells, dims = 1:20)
T_NK_cells <- FindClusters(T_NK_cells, resolution = 0.5)
T_NK_cells <- RunUMAP(T_NK_cells, dims = 1:20)
DimPlot(T_NK_cells, reduction = "umap", label = TRUE)
DimPlot(T_NK_cells, reduction = "umap", group.by = "sub_cell_type", label = TRUE) 
saveRDS(T_NK_cells, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.T-NK.seurat.rds")
#saveRDS(NSCLC.object, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.all.seurat.rds")



# subset
# CD8
CD8_Tcells <- subset(NSCLC.object, subset = sub_cell_type %in% c(
    "CD8T_Tem_GZMK+GZMH+", 
    "CD8T_Tex_CXCL13", 
    "CD8T_terminal_Tex_LAYN", 
    "CD8T_Tm_IL7R", 
    "CD8T_Trm_ZNF683", 
    "CD8T_ISG15", 
    "CD8T_prf_MKI67", 
    "CD8T_Tem_GZMK+NR4A1+", 
    "CD8T_NK-like_FGFBP2", 
    "CD8T_MAIT_KLRB1"
  ))
CD8_Tcells <- NormalizeData(CD8_Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_Tcells <- FindVariableFeatures(CD8_Tcells, selection.method = "vst", nfeatures = 2000)
CD8_Tcells <- ScaleData(CD8_Tcells, features = rownames(CD8_Tcells))
CD8_Tcells <- RunPCA(CD8_Tcells, features = VariableFeatures(object = CD8_Tcells))
ElbowPlot(CD8_Tcells)
CD8_Tcells <- FindNeighbors(CD8_Tcells, dims = 1:20)
CD8_Tcells <- FindClusters(CD8_Tcells, resolution = 0.5)
CD8_Tcells <- RunUMAP(CD8_Tcells, dims = 1:20)
DimPlot(CD8_Tcells, reduction = "umap", label = TRUE) + 
  ggtitle("CD8+ T cell sub cluster") 
DimPlot(CD8_Tcells, reduction = "umap", group.by = "sub_cell_type", label = TRUE) + 
  ggtitle("CD8+ T cell sub cluster") 

FeaturePlot(T_NK_cells, features = c("TREG"), reduction = "umap") 

# CD4
CD4_Tcells <- subset(NSCLC.object, subset = sub_cell_type %in% c(
  "CD4T_Tm_ANXA1",
  "CD4T_Tfh_CXCL13",
  "CD4T_Tn_CCR7" ,
  "CD4T_Treg_CCR8",        
  "CD4T_Treg_FOXP3", 
  "CD4T_Th1-like_CXCL13",  
  "CD4T_Tem_GZMA",
  "CD4T_Tm_XCL1", 
  "CD4T_Treg_MKI67"
))
CD4_Tcells <- NormalizeData(CD4_Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_Tcells <- FindVariableFeatures(CD4_Tcells, selection.method = "vst", nfeatures = 2000)
CD4_Tcells <- ScaleData(CD4_Tcells, features = rownames(CD4_Tcells))
CD4_Tcells <- RunPCA(CD4_Tcells, features = VariableFeatures(object = CD4_Tcells))
ElbowPlot(CD4_Tcells)
CD4_Tcells <- FindNeighbors(CD4_Tcells, dims = 1:20)
CD4_Tcells <- FindClusters(CD4_Tcells, resolution = 0.5)
CD4_Tcells <- RunUMAP(CD4_Tcells, dims = 1:20)
DimPlot(CD4_Tcells, reduction = "umap", label = TRUE) + 
  ggtitle("CD4+ T cell sub cluster") 
DimPlot(CD4_Tcells, reduction = "umap", group.by = "sub_cell_type", label = TRUE) + 
  ggtitle("CD4+ T cell sub cluster") 




NSCLC.object
T_NK_cells
CD8_Tcells

FeaturePlot(CD8_Tcells, features = c("FOXP3", "CXCL13"), reduction = "umap", order = TRUE) 
VlnPlot(NSCLC.object, features = c("FOXP3", "CXCL13"), pt.size = 0)



#saveRDS(CD8_Tcells, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.CD8T.seurat.rds")
#saveRDS(NSCLC.object, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.all.seurat.rds")




# how many cells are there 
# Count cells with FOXP3 expression > 0 in CD8_Tcells
foxp3_expr_counts <- GetAssayData(CD8_Tcells, assay = "RNA", slot = "counts")["FOXP3", ]
num_foxp3_positive <- sum(foxp3_expr_counts > 0)
total_cd8_cells <- ncol(CD8_Tcells)
cat("FOXP3+ CD8 T cells:", num_foxp3_positive, "out of", total_cd8_cells, "\n")

cxcl13_expr_counts <- GetAssayData(CD8_Tcells, assay = "RNA", slot = "counts")["CXCL13", ]
num_cxcl13_positive <- sum(cxcl13_expr_counts > 0)
total_cd8_cells <- ncol(CD8_Tcells)
cat("CXCL13+ CD8 T cells:", num_cxcl13_positive, "out of", total_cd8_cells, "\n")

combined_expr_counts <- GetAssayData(CD8_Tcells, assay = "RNA", slot = "counts")[c("FOXP3", "CXCL13"), ]
num_combined_positive <- sum(combined_expr_counts > 0)
total_cd8_cells <- ncol(CD8_Tcells)
cat("FOXP3+ CXCL13+ CD8 T cells:", num_foxp3_positive, "out of", total_cd8_cells, "\n")

str(combined_expr_counts)
str(foxp3_expr_counts)


# Violin
violin_data <- data.frame(
  FOXP3 = as.numeric(foxp3_expr_counts),
  CXCL13 = as.numeric(cxcl13_expr_counts)
)
library(tidyr)
violin_data_long <- pivot_longer(violin_data, cols = everything(),
                                 names_to = "Gene", values_to = "Expression")
library(ggplot2)
ggplot(violin_data_long, aes(x = Gene, y = Expression)) +
  geom_violin(fill = "skyblue", color = "black", scale = "width") +
  geom_jitter(height = 0, width = 0.2, alpha = 0.2, shape = 16) +  # shape 16 = filled circle, no border
  theme_minimal() +
  ylab("Expression Level") +
  ggtitle("FOXP3 vs CXCL13 Expression in CD8 T cells")



# Venn
library(VennDiagram)
foxp3_only <- 2382
cxcl13_only <- 56949
both <- 598
# Draw Venn diagram with smaller text and adjusted positioning
venn.plot <- draw.pairwise.venn(
  area1 = foxp3_only,             # FOXP3+
  area2 = cxcl13_only,            # CXCL13+
  cross.area = both,        # overlap
  category = c("", ""),      # Empty strings for categories (hide them)
  fill = c("#1F78B4", "#33A02C"),
  lty = "blank",
  cex = 1.2,                   # Count size
  cat.cex = 0,                 # Hide labels (extra safe)
  cat.dist = c(0, 0),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5))
)



# Treg Subset
Treg <- subset(NSCLC.object, subset = sub_cell_type %in% c(
  "CD4T_Treg_CCR8", 
  "CD4T_Treg_FOXP3", 
  "CD4T_Treg_MKI67"))
Treg <- NormalizeData(Treg, normalization.method = "LogNormalize", scale.factor = 10000)
Treg <- FindVariableFeatures(Treg, selection.method = "vst", nfeatures = 2000)
Treg <- ScaleData(Treg, features = rownames(Treg))
Treg <- RunPCA(Treg, features = VariableFeatures(object = Treg))
ElbowPlot(Treg)
Treg <- FindNeighbors(Treg, dims = 1:20)
Treg <- FindClusters(Treg, resolution = 0.5)
Treg <- RunUMAP(Treg, dims = 1:20)
DimPlot(Treg, reduction = "umap", label = TRUE) + 
  ggtitle("Regulatory T cell sub cluster") 
DimPlot(Treg, reduction = "umap", group.by = "sub_cell_type", label = TRUE) + 
  ggtitle("Regulatory T cell sub cluster") 
FeaturePlot(Treg, features = c("CD8A", "CD8B"), reduction = "umap", order = TRUE) 
FeaturePlot(Treg, features = c("CD4"), reduction = "umap", order = TRUE) 
saveRDS(CD8_Treg, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.CD8_Treg.seurat.rds")





4/5
#导入之前存好的rds
NSCLC.object <- readRDS("NSCLC.zemin.all.seurat.rds")
CD8_Treg <- readRDS("/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.CD8_Treg.seurat.rds")


# 看眼metadata在不在 不在的话
# CD8_Treg Clonotype
CD8_Treg <- subset(Treg, subset = CD8A > 0 | CD8B > 0)
rownames(TCR.data) <- TCR.data$cellID
clonotype.meta <- TCR.data["clonotype"]
clonotype_number.meta <- TCR.data["clonotype_number"]
CD8_Treg <- AddMetaData(CD8_Treg, metadata = clonotype.meta)
CD8_Treg <- AddMetaData(CD8_Treg, metadata = clonotype_number.meta)
CD8_Treg$clonotype <- sub("^[^_]+_", "", CD8_Treg$clonotype)
head(CD8_Treg@meta.data)


# FOXP3+ /CXCL13+ CD8+ Treg cells (based on expression)
FOXP3_Cells <- WhichCells(CD8_Tcells, expression = FOXP3 > 0)
CXCL13_Cells <- WhichCells(CD8_Tcells, expression = CXCL13 > 0)
FOXP3_clones <- CD8_Tcells$clonotype[FOXP3_Cells]
CXCL13_clones <- CD8_Tcells$clonotype[CXCL13_Cells]

shared_clones <- intersect(FOXP3_clones, CXCL13_clones)
length(CXCL13_clones)
length(FOXP3_clones)
length(shared_clones)
length(unique(CD8_Tcells$clonotype))

draw.pairwise.venn(
  area1 = length(FOXP3_clones),
  area2 = length(CXCL13_clones),
  cross.area = length(shared_clones),
  category = c("", ""),
  fill = c("#1F78B4", "#33A02C"),
  alpha = 0.5,
  cat.cex = 1.2
)

# CXCL13+ CD8+ Treg in Treg
DimPlot(CD8_Tcells, cells.highlight = CXCL13_Cells, 
        cols.highlight = "red", 
        cols = "grey80", 
        pt.size = 1) + 
  ggtitle("CXCL13+ CD8_Tcells Cells")


# TCR - heatmap 但value太低了不方便看
foxp3_expr <- FetchData(CD8_Tcells, vars = "FOXP3")
cxcl13_expr <- FetchData(CD8_Tcells, vars = "CXCL13")
treg_marker <- ifelse(foxp3_expr$FOXP3 > 0 & cxcl13_expr$CXCL13 == 0, "FOXP3+",
                      ifelse(foxp3_expr$FOXP3 == 0 & cxcl13_expr$CXCL13 > 0, "CXCL13+",
                             ifelse(foxp3_expr$FOXP3 > 0 & cxcl13_expr$CXCL13 > 0, "FOXP3+.CXCL13+",
                                    "neither")))
CD8_Tcells$treg.marker <- treg_marker
CD8_Treg_TCR.matrix <- CD8_Tcells@meta.data[, c("cellID", "clonotype", "clonotype_number", "treg.marker")]
CD8_Treg_TCR.matrix <- CD8_Treg_TCR.matrix[!is.na(CD8_Treg_TCR.matrix$clonotype), ]
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(tibble)
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("/projects/b1198/epifluidlab/zhengyue/spring")

pivot_matrix <- CD8_Treg_TCR.matrix %>%
  select(cellID, clonotype, clonotype_number) %>%
  pivot_wider(names_from = clonotype, values_from = clonotype_number, values_fill = list(clonotype_number = 0))
annotation_row <- CD8_Treg_TCR.matrix %>%
  select(cellID, treg.marker) %>%
  distinct()
annotation_row <- as.data.frame(annotation_row)
rownames(annotation_row) <- annotation_row$cellID
annotation_row$cellID <- NULL
ordered_cells <- annotation_row %>%
  arrange(treg.marker) %>%
  rownames()
library(tibble)
pivot_matrix <- pivot_matrix %>%
  column_to_rownames(var = "cellID")

pivot_matrix <- pivot_matrix[ordered_cells, ]
annotation_row <- annotation_row[ordered_cells, , drop=FALSE]
annotation_colors = list(
  treg.marker = c(
    "FOXP3+" = "#E31A1C",
    "FOXP3+.CXCL13+" = "#FDBF6F",
    "CXCL13+" = "#33A02C",
    "neither" = "gray70"
  )
)
pheatmap(
  pivot_matrix,
  cluster_rows = FALSE,    # group by treg.marker
  cluster_cols = TRUE,     # cluster clonotypes
  color = viridis(100),
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,   # cleaner plot
  show_colnames = FALSE,   # optionally hide clonotype names if too many
  fontsize = 10,
  border_color = NA,
  main = "Clonal Origin Tracking of CD8+ Treg"  
)



# Dotplot
library(dplyr)
library(tidyr)
library(ggplot2)  # For dot plot
pivot_matrix <- CD8_Treg_TCR.matrix %>%
  select(cellID, clonotype, clonotype_number) %>%
  pivot_wider(names_from = clonotype, values_from = clonotype_number, values_fill = list(clonotype_number = 0))
annotation_row <- CD8_Treg_TCR.matrix %>%
  select(cellID, treg.marker) %>%
  distinct()
annotation_row <- as.data.frame(annotation_row)
rownames(annotation_row) <- annotation_row$cellID
annotation_row$cellID <- NULL
ordered_cells <- annotation_row %>%
  arrange(treg.marker) %>%
  rownames()
pivot_matrix <- pivot_matrix %>%
  column_to_rownames(var = "cellID")
pivot_matrix <- pivot_matrix[ordered_cells, ]
annotation_row <- annotation_row[ordered_cells, , drop=FALSE]
pivot_matrix_long <- pivot_matrix %>%
  rownames_to_column(var = "cellID") %>%
  pivot_longer(cols = -cellID, names_to = "clonotype", values_to = "clonotype_number") %>%
  filter(clonotype_number > 0)  # Keep only rows with clonotype_number > 0
annotation_row <- annotation_row %>%
  tibble::rownames_to_column(var = "cellID")
pivot_matrix_long <- pivot_matrix_long %>%
  left_join(annotation_row, by = "cellID")
str(pivot_matrix_long)
pivot_matrix_long$cellID <- factor(pivot_matrix_long$cellID, levels = ordered_cells)
ggplot(pivot_matrix_long, aes(x = clonotype, y = cellID, color = treg.marker)) +
  geom_point(aes(size = clonotype_number), alpha = 0.7) +
  scale_size_continuous(range = c(1, 6)) +
  scale_color_manual(values = c("FOXP3+" = "#33A02C" , 
                                "FOXP3+.CXCL13+" = "#FDBF6F",
                                "CXCL13+" = "#E31A1C", 
                                "neither" = "gray70")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Hide clonotype names if too many
    axis.text.y = element_blank(), # Hide cellID names if too many
    axis.ticks = element_blank(),
    panel.grid = element_blank() # Optional, remove grid lines
  ) +
  labs(
    x = "Clonotype",
    y = "Cell",
    size = "Clonotype\nNumber",
    title = "Clonal Origin Tracking of CD8+ Treg"
  )





4/7
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(tibble)
library(Seurat)
library(Matrix)
library(ggplot2)
setwd("/projects/b1198/epifluidlab/zhengyue/spring")
library(data.table)
#导入之前存好的rds
NSCLC.object <- readRDS("NSCLC.zemin.all.seurat.rds") #需要补上clonotype metadata
CD8_Tcells <- readRDS("NSCLC.zemin.CD8T.seurat.rds") #需要补上clonotype metadata

# 给NSCLC.object加上metadata
TCR.data <- fread("GSE243013_T_with_TCR_annotation.csv")
rownames(TCR.data) <- TCR.data$cellID
clonotype.meta <- TCR.data[["clonotype"]]
clonotype_number.meta <- TCR.data[["clonotype_number"]]
names(clonotype.meta) <- TCR.data$cellID  # Important: set names = cell IDs!!
NSCLC.object <- AddMetaData(NSCLC.object, metadata = clonotype.meta, col.name = "clonotype")
names(clonotype_number.meta) <- TCR.data$cellID  # Important: set names = cell IDs!!
NSCLC.object <- AddMetaData(NSCLC.object, metadata = clonotype_number.meta, col.name = "clonotype_number")
NSCLC.object$clonotype <- sub("^[^_]+_", "", NSCLC.object$clonotype)
head(NSCLC.object@meta.data)

# 给CD8_Tcells加上metadata
TCR.data <- fread("GSE243013_T_with_TCR_annotation.csv")
rownames(TCR.data) <- TCR.data$cellID
clonotype.meta <- TCR.data[["clonotype"]]
clonotype_number.meta <- TCR.data[["clonotype_number"]]
names(clonotype.meta) <- TCR.data$cellID  # Important: set names = cell IDs!!
CD8_Tcells <- AddMetaData(CD8_Tcells, metadata = clonotype.meta, col.name = "clonotype")
names(clonotype_number.meta) <- TCR.data$cellID  # Important: set names = cell IDs!!
CD8_Tcells <- AddMetaData(CD8_Tcells, metadata = clonotype_number.meta, col.name = "clonotype_number")
CD8_Tcells$clonotype <- sub("^[^_]+_", "", CD8_Tcells$clonotype)
head(CD8_Tcells@meta.data)

# 重新跑Treg
# Treg Subset
Treg <- subset(NSCLC.object, subset = sub_cell_type %in% c(
  "CD4T_Treg_CCR8", 
  "CD4T_Treg_FOXP3", 
  "CD4T_Treg_MKI67"))
Treg <- NormalizeData(Treg, normalization.method = "LogNormalize", scale.factor = 10000)
Treg <- FindVariableFeatures(Treg, selection.method = "vst", nfeatures = 2000)
Treg <- ScaleData(Treg, features = rownames(Treg))
Treg <- RunPCA(Treg, features = VariableFeatures(object = Treg))
ElbowPlot(Treg)
Treg <- FindNeighbors(Treg, dims = 1:20)
Treg <- FindClusters(Treg, resolution = 0.5)
Treg <- RunUMAP(Treg, dims = 1:20)
DimPlot(Treg, reduction = "umap", label = TRUE) + 
  ggtitle("Regulatory T cell sub cluster") 
DimPlot(Treg, reduction = "umap", group.by = "sub_cell_type", label = TRUE) + 
  ggtitle("Regulatory T cell sub cluster") 
FeaturePlot(Treg, features = c("CD8A", "CD8B"), reduction = "umap", order = TRUE) 
FeaturePlot(Treg, features = c("CD4"), reduction = "umap", order = TRUE) 
#saveRDS(Treg, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.Treg.seurat.rds")

# Alex request:
# subcluster CD8+ FOXP3+ Treg cells. Find the TCRs for the CD8+ FOXP3+ cells. Save that as a list.
Treg_CD8 <- subset(Treg, subset = CD8A > 0 | CD8B > 0)
Treg_CD8_FOXP3 <- subset(Treg_CD8, subset = FOXP3 > 0)
Treg_CD8_FOXP3.TCR <- Treg_CD8_FOXP3@meta.data$clonotype
Treg_CD8_FOXP3.TCR <- Treg_CD8_FOXP3.TCR[!is.na(Treg_CD8_FOXP3.TCR)]
length(Treg_CD8_FOXP3.TCR)

Treg_CD8_FOXP3.TCR.unique <- unique(Treg_CD8_FOXP3.TCR)
length(Treg_CD8_FOXP3.TCR.unique)

# Now take the CXCL13+ CD8s (from the original CD8 clustering) and get their TCRs.
CD8_CXCL13 <- subset(CD8_Tcells, subset = CXCL13 > 0)
CD8_CXCL13.TCR <- CD8_CXCL13@meta.data$clonotype
CD8_CXCL13.TCR <- CD8_CXCL13.TCR[!is.na(CD8_CXCL13.TCR)]
length(CD8_CXCL13.TCR)
CD8_CXCL13.TCR.unique <- unique(CD8_CXCL13.TCR)
length(CD8_CXCL13.TCR.unique)

# see all CD8 TCR
CD8.TCR <- CD8_Tcells@meta.data$clonotype
CD8.TCR <- CD8.TCR[!is.na(CD8.TCR)]
length(CD8.TCR)
CD8.TCR.unique <- unique(CD8.TCR)
length(CD8.TCR.unique)



# is there overlap? 
overlap1 <- intersect(CD8_CXCL13.TCR, Treg_CD8_FOXP3.TCR)
length(overlap1) 

matching_cells <- CD8_Tcells@meta.data$clonotype %in% Treg_CD8_FOXP3.TCR
sum(matching_cells) # how many cells in CD8 has the matching clonotype from

overlap3 <- intersect(CD8.TCR.unique, Treg_CD8_FOXP3.TCR.unique)
length(overlap3) 



# is there overlap? unique
overlap.unique <- intersect(CD8_CXCL13.TCR.unique, Treg_CD8_FOXP3.TCR.unique)
length(overlap.unique) 

library(VennDiagram)
grid::grid.newpage()
draw.pairwise.venn(
  area1 = length(CD8_CXCL13.TCR.unique),
  area2 = length(Treg_CD8_FOXP3.TCR.unique),
  cross.area = length(overlap.unique),
  category = c("", ""),
  fill = c("#1F78B4", "#33A02C"),
  alpha = 0.5,
  cat.cex = 1.2
)



# Umap 最终版
# Object: CD8; Orange: CXCL13+;  Black: TCR matches with FOXP3+ CD8
# TCR matching
# CD8_FOXP3 <- subset(CD8_Tcells, subset = FOXP3 > 0)

# FOXP3+ CD8
CD8_FOXP3 <- subset(CD8_Tcells, subset = FOXP3 > 0)
CD8_FOXP3.TCR <- CD8_FOXP3@meta.data$clonotype
CD8_FOXP3.TCR <- CD8_FOXP3.TCR[!is.na(CD8_FOXP3.TCR)]
length(CD8_FOXP3.TCR)
CD8_FOXP3.TCR.unique <- unique(CD8_FOXP3.TCR)
length(CD8_FOXP3.TCR.unique)

CD8_FOXP3.TCR <- data.frame(
  CellID = rownames(CD8_FOXP3@meta.data),
  Clonotype = CD8_FOXP3@meta.data$clonotype
)
CD8_FOXP3.TCR <- CD8_FOXP3.TCR[!is.na(CD8_FOXP3.TCR$Clonotype), ]
str(CD8_FOXP3.TCR)

CD8_FOXP3.TCR.ClonoList <- CD8_FOXP3.TCR$Clonotype
length(CD8_FOXP3.TCR.ClonoList)
head(CD8_FOXP3.TCR.ClonoList)

FeaturePlot(
  object = CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
)  + 
  ggtitle("CXCL13 Expression in CD8 T Cells")

# Create a data frame for the cell IDs in CD8_FOXP3.TCR
marked_cells <- CD8_FOXP3.TCR$CellID

# Feature plot for CXCL13 expression in CD8 T cells
p <- FeaturePlot(
  object = CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
) + 
  ggtitle("CXCL13 Expression in CD8 T Cells")
p

# Add the little black triangles for the matching cells using the correct UMAP coordinates
p + 
  geom_point(data = subset(CD8_Tcells@meta.data, rownames(CD8_Tcells@meta.data) %in% marked_cells), 
             aes(x = CD8_Tcells@reductions$umap@cell.embeddings[rownames(CD8_Tcells@meta.data) %in% marked_cells, 1], 
                 y = CD8_Tcells@reductions$umap@cell.embeddings[rownames(CD8_Tcells@meta.data) %in% marked_cells, 2]), 
             shape = 17, color = "black", size = 0.3)

pdf("FOXP3_CD8_and_CD8_overlap_umap.pdf", width = 6, height = 5)  # 可以调整宽高

p + 
  geom_point(
    data = subset(CD8_Tcells@meta.data, rownames(CD8_Tcells@meta.data) %in% marked_cells),
    aes(
      x = CD8_Tcells@reductions$umap@cell.embeddings[rownames(CD8_Tcells@meta.data) %in% marked_cells, 1],
      y = CD8_Tcells@reductions$umap@cell.embeddings[rownames(CD8_Tcells@meta.data) %in% marked_cells, 2]
    ),
    shape = 17,   # 三角形
    color = "black", 
    size = 0.3    # 大小
  )

# 关闭PDF文件，保存图
dev.off()


#Jon: now, do the umap, merge CD8+ Treg to the CD8 T
CD8_Tcells_raw <- subset(NSCLC.object, subset = sub_cell_type %in% c(
  "CD8T_Tem_GZMK+GZMH+", 
  "CD8T_Tex_CXCL13", 
  "CD8T_terminal_Tex_LAYN", 
  "CD8T_Tm_IL7R", 
  "CD8T_Trm_ZNF683", 
  "CD8T_ISG15", 
  "CD8T_prf_MKI67", 
  "CD8T_Tem_GZMK+NR4A1+", 
  "CD8T_NK-like_FGFBP2"))
Treg_CD8_raw <- subset(NSCLC.object, subset = sub_cell_type %in% c(
  "CD4T_Treg_CCR8", 
  "CD4T_Treg_FOXP3", 
  "CD4T_Treg_MKI67"))
Treg_CD8_raw <- subset(Treg_CD8_raw, subset = CD8A > 0 | CD8B > 0)
combined_CD8_Tcells <- merge(CD8_Tcells_raw, y = Treg_CD8_raw)

combined_CD8_Tcells <- NormalizeData(combined_CD8_Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
combined_CD8_Tcells <- FindVariableFeatures(combined_CD8_Tcells, selection.method = "vst", nfeatures = 2000)
combined_CD8_Tcells <- ScaleData(combined_CD8_Tcells, features = rownames(combined_CD8_Tcells))
combined_CD8_Tcells <- RunPCA(combined_CD8_Tcells, features = VariableFeatures(object = combined_CD8_Tcells))
ElbowPlot(combined_CD8_Tcells)
combined_CD8_Tcells <- FindNeighbors(combined_CD8_Tcells, dims = 1:20)
combined_CD8_Tcells <- FindClusters(combined_CD8_Tcells, resolution = 0.5)
combined_CD8_Tcells <- RunUMAP(combined_CD8_Tcells, dims = 1:20)
DimPlot(combined_CD8_Tcells, reduction = "umap", label = TRUE) + 
  ggtitle("CD8+ Treg and CD8 cell sub cluster") 
DimPlot(combined_CD8_Tcells, reduction = "umap", group.by = "sub_cell_type", label = TRUE) + 
  ggtitle("CD8+ Treg and CD8 cell sub cluster") 

FeaturePlot(
  object = combined_CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
)  + 
  ggtitle("CXCL13 Expression in CD8 T Cells")
FeaturePlot(
  object = combined_CD8_Tcells,
  features = "FOXP3",
  cols = c("grey", "orange")
)  + 
  ggtitle("FOXP3 Expression in CD8 T Cells")

# Create a data frame for the cell IDs in CD8_FOXP3.TCR
marked_cells <- CD8_FOXP3.TCR$CellID

# Feature plot for CXCL13 expression in CD8 T cells
p <- FeaturePlot(
  object = combined_CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
) + 
  ggtitle("CXCL13 Expression in Combined CD8 T Cells")

# Add the little black triangles for the matching cells using the correct UMAP coordinates
p + 
  geom_point(data = subset(combined_CD8_Tcells@meta.data, rownames(combined_CD8_Tcells@meta.data) %in% marked_cells), 
             aes(x = combined_CD8_Tcells@reductions$umap@cell.embeddings[rownames(combined_CD8_Tcells@meta.data) %in% marked_cells, 1], 
                 y = combined_CD8_Tcells@reductions$umap@cell.embeddings[rownames(combined_CD8_Tcells@meta.data) %in% marked_cells, 2]), 
             shape = 17, color = "black", size = 0.3)

saveRDS(Treg, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.Treg.seurat.rds")
saveRDS(combined_CD8_Tcells, file = "/projects/b1198/epifluidlab/zhengyue/spring/NSCLC.zemin.combined_CD8.seurat.rds")

FeaturePlot(
  object = combined_CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
) 
gc()


4/9
# redo-umap!
# subcluster CD8+ FOXP3+ Treg cells. Find the TCRs for the CD8+ FOXP3+ cells. Save that as a list.
Treg_CD8 <- subset(Treg, subset = CD8A > 0 | CD8B > 0)
Treg_CD8_FOXP3 <- subset(Treg_CD8, subset = FOXP3 > 0)
Treg_CD8_FOXP3.TCR <- Treg_CD8_FOXP3@meta.data$clonotype
Treg_CD8_FOXP3.TCR <- Treg_CD8_FOXP3.TCR[!is.na(Treg_CD8_FOXP3.TCR)]
length(Treg_CD8_FOXP3.TCR)

Treg_CD8_FOXP3.TCR <- data.frame(
  CellID = rownames(Treg_CD8_FOXP3@meta.data),
  Clonotype = Treg_CD8_FOXP3@meta.data$clonotype
)
Treg_CD8_FOXP3.TCR <- Treg_CD8_FOXP3.TCR[!is.na(Treg_CD8_FOXP3.TCR$Clonotype), ]
str(Treg_CD8_FOXP3.TCR)

FeaturePlot(
  object = CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
)  + 
  ggtitle("CXCL13 Expression in CD8 T Cells")

# List of TCR in FOXP3+ CD8+ Treg
marked_clonotypes <- Treg_CD8_FOXP3.TCR$Clonotype
marked_clonotypes.unique <- unique(marked_clonotypes)
length(marked_clonotypes.unique)

# Feature plot for CXCL13 expression in CD8 T cells
p <- FeaturePlot(
  object = CD8_Tcells,
  features = "CXCL13",
  cols = c("grey", "orange")
) + 
  ggtitle("CXCL13 Expression in CD8 T Cells")
p

# Step 1: Create a logical vector to find which clonotypes are in marked_clonotypes
marked_cells <- rownames(CD8_Tcells@meta.data)[CD8_Tcells@meta.data$clonotype %in% marked_clonotypes.unique]

# Step 2: Plot the cells
p + 
  geom_point(data = subset(CD8_Tcells@meta.data, rownames(CD8_Tcells@meta.data) %in% marked_cells), 
             aes(x = CD8_Tcells@reductions$umap@cell.embeddings[rownames(CD8_Tcells@meta.data) %in% marked_cells, 1], 
                 y = CD8_Tcells@reductions$umap@cell.embeddings[rownames(CD8_Tcells@meta.data) %in% marked_cells, 2]), 
             shape = 17, color = "black", size = 0.05)




# null hypothesis: CXCL13+ CD8 from Zemin Zhang lung cancer dataset are clonally related to FOXP3+ CD8+ Treg.
# clonal overlap enrichment analysis 

# 1. collect clonal distribution of FOXP3+ CD8+ Treg TCR
TCR.data <- fread("GSE243013_T_with_TCR_annotation.csv")
rownames(TCR.data) <- TCR.data$cellID
clonotype.meta <- TCR.data[["clonotype"]]
clonotype_number.meta <- TCR.data[["clonotype_number"]]
TRB_cdr3.meta <- TCR.data[["TRB_cdr3"]]
names(clonotype.meta) <- TCR.data$cellID
names(clonotype_number.meta) <- TCR.data$cellID
names(TRB_cdr3.meta) <- TCR.data$cellID

# FOXP3+ CD8+ Treg 
Treg_CD8_FOXP3 <- AddMetaData(Treg_CD8_FOXP3, metadata = clonotype.meta, col.name = "clonotype")
Treg_CD8_FOXP3 <- AddMetaData(Treg_CD8_FOXP3, metadata = clonotype_number.meta, col.name = "clonotype_number")
Treg_CD8_FOXP3 <- AddMetaData(Treg_CD8_FOXP3, metadata = TRB_cdr3.meta, col.name = "TRB_cdr3")
Treg_CD8_FOXP3$clonotype <- sub("^[^_]+_", "", Treg_CD8_FOXP3$clonotype)
Treg_CD8_FOXP3.TCR <- data.frame(
  CellID = rownames(Treg_CD8_FOXP3@meta.data),
  Clonotype = Treg_CD8_FOXP3@meta.data$clonotype,
  TRB_cdr3 = Treg_CD8_FOXP3@meta.data$TRB_cdr3
)
Treg_CD8_FOXP3.TCR <- Treg_CD8_FOXP3.TCR[!is.na(Treg_CD8_FOXP3.TCR$Clonotype), ]
str(Treg_CD8_FOXP3.TCR)

# CD8 all clusters
CD8_Tcells <- AddMetaData(CD8_Tcells, metadata = clonotype.meta, col.name = "clonotype")
CD8_Tcells <- AddMetaData(CD8_Tcells, metadata = clonotype_number.meta, col.name = "clonotype_number")
CD8_Tcells <- AddMetaData(CD8_Tcells, metadata = TRB_cdr3.meta, col.name = "TRB_cdr3")
CD8_Tcells$clonotype <- sub("^[^_]+_", "", CD8_Tcells$clonotype)
CD8_Tcells.TCR <- data.frame(
  CellID = rownames(CD8_Tcells@meta.data),
  Cluster = CD8_Tcells@meta.data$sub_cell_type,
  Clonotype = CD8_Tcells@meta.data$clonotype,
  TRB_cdr3 = CD8_Tcells@meta.data$TRB_cdr3
)
CD8_Tcells.TCR <- CD8_Tcells.TCR[!is.na(CD8_Tcells.TCR$Clonotype), ]
str(CD8_Tcells.TCR)

set.seed(1234)

#statistical analysis
CD8_data <- CD8_Tcells.TCR
FOXP3_data <- Treg_CD8_FOXP3.TCR

clusters <- unique(CD8_data$Cluster)
n_permutations <- 1000

# Step 1: Observed overlap
observed_overlap <- sapply(clusters, function(cluster) {
  cluster_cdr3 <- CD8_data$TRB_cdr3[CD8_data$Cluster == cluster]
  sum(cluster_cdr3 %in% FOXP3_data$TRB_cdr3)
})

# Step 2: Permutations
null_distributions <- matrix(NA, nrow = n_permutations, ncol = length(clusters))
colnames(null_distributions) <- clusters

for (i in 1:n_permutations) {
  shuffled_cdr3 <- sample(CD8_data$TRB_cdr3)
  CD8_data$Shuffled_cdr3 <- shuffled_cdr3
  
  null_distributions[i, ] <- sapply(clusters, function(cluster) {
    cluster_cdr3 <- CD8_data$Shuffled_cdr3[CD8_data$Cluster == cluster]
    sum(cluster_cdr3 %in% FOXP3_data$TRB_cdr3)
  })
}

# Step 3: Calculate p-values and Z-scores
mean_null <- colMeans(null_distributions)
sd_null <- apply(null_distributions, 2, sd)

Z_scores <- (observed_overlap - mean_null) / sd_null

p_values <- sapply(1:length(observed_overlap), function(j) {
  mean(null_distributions[, j] >= observed_overlap[j])
})

FDR_values <- p.adjust(p_values, method = "fdr")

results <- data.frame(
  Cluster = clusters,
  Observed = observed_overlap,
  Z_score = Z_scores,
  P_value = p_values,
  FDR = FDR_values
)

# Step 4: Plotting
ggplot(results, aes(x = Cluster, y = Z_score)) +
  geom_point(color = "black") +   # 先画所有点（黑色）
  geom_hline(yintercept = 0, linetype = "dotted") +  # 中间水平线
  geom_point(data = subset(results, FDR < 0.2),
             aes(x = Cluster, y = Z_score),
             color = "red", size = 3) +  # 再覆盖画FDR<0.2的红点
  geom_text_repel(data = subset(results, FDR < 0.05),
                  aes(label = paste0(Cluster, "\nFDR=", signif(FDR, 2))),
                  color = "red", box.padding = 0.5, max.overlaps = Inf) +  # FDR<0.05才加label
  theme_classic() +
  labs(title = "TCR Overlap between FOXP3+ CD8+ Treg and CD8 Clusters",
       y = "Z-score", x = "Clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#画null distribution
library(ggplot2)

# 假设你要画的 cluster
target_cluster <- "CD8T_Tex_CXCL13"
target_cluster <- "CD8T_ISG15"

# 提取这个 cluster 的 null distribution 和 observed value
null_values <- null_distributions[, target_cluster]
observed_value <- observed_overlap[target_cluster]
# 画 null distribution
df <- data.frame(overlap = null_values)

ggplot(df, aes(x = overlap)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", size = 1) +
  theme_classic() +
  labs(title = paste0("Null Distribution for ", target_cluster),
       x = "Overlap after permutation",
       y = "Count") +
  annotate("text", 
           x = observed_value, 
           y = max(table(cut(null_values, 30))), 
           label = paste0("Observed = ", observed_value),
           vjust = -1, color = "red")











