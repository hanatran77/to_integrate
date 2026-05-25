# ── Q3.7: DATA INTEGRATION MULTI-METHOD WORKFLOW ENGINE ──────────────────────
library(Seurat)
library(harmony)
library(patchwork)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(dplyr)

# 1. Load the multi-sample dataset
seurat_list <- readRDS("C:/Users/Hana/Desktop/Screening/data_scRNA/to_integrate.rds")
print(seurat_list)

# 2. Add an explicit donor/batch flag to each object metadata before joining them
# (This guarantees we have a clean batch variable for group processing)
#for (i in 1:length(seurat_list)) {
#  seurat_list[[i]]$donor_batch <- paste0("Donor_", i)
#}

# 3. Merge the list into a single unified Seurat v5 multi-layer object
# Using join="outer" creates the standard Seurat v5 multi-layer infrastructure
seurat_integrated <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[2:length(seurat_list)], 
  add.cell.ids = c("C1", "C2", "AD1", "AD2"), 
  project = "AD_GSE287652"
)

#Set the cell types column as our active classification system
Idents(seurat_integrated) <- seurat_integrated$paper_annot
# Join the molecule layers to establish a clean base tracking profile
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])

#==============================================================================
# ── APPROACH 1: SIMPLE MERGE (WITHOUT BATCH CORRECTION) ──────────────────────
print("Processing Approach 1: Simple Merge...")
seurat_merge <- NormalizeData(seurat_integrated, verbose = FALSE)
seurat_merge <- FindVariableFeatures(seurat_merge, nfeatures = 2000, verbose = FALSE)
seurat_merge <- ScaleData(seurat_merge, verbose = FALSE)
seurat_merge <- RunPCA(seurat_merge, reduction.name = "pca.merge", verbose = FALSE)
seurat_merge <- RunUMAP(seurat_merge, reduction = "pca.merge", dims = 1:15, reduction.name = "umap.merge", verbose = FALSE)

#==============================================================================
# ── APPROACH 2: RECIPROCAL PCA (RPCA) INTEGRATION ────────────────────────────
print("Processing Approach 2: Reciprocal PCA...")
seurat_rpca <- seurat_integrated

# Split layers based on our sample_id variable (Mandatory for Seurat v5)
seurat_rpca[["RNA"]] <- split(seurat_rpca[["RNA"]], f = seurat_rpca$sample_id)

# Process pre-requisite layers for RPCA
seurat_rpca <- NormalizeData(seurat_rpca, verbose = FALSE)
seurat_rpca <- FindVariableFeatures(seurat_rpca, nfeatures = 2000, verbose = FALSE)
seurat_rpca <- ScaleData(seurat_rpca, verbose = FALSE)
seurat_rpca <- RunPCA(seurat_rpca, verbose = FALSE)

# Adjust memory limits due to the large matrix load
#options(future.globals.maxSize = Inf)

# Discover anchors and compute layer alignment
seurat_rpca <- IntegrateLayers(
  object = seurat_rpca, 
  assay = "RNA", 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca", 
  verbose = FALSE
)
seurat_rpca <- RunUMAP(seurat_rpca, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca", verbose = FALSE)

#==============================================================================
# ── APPROACH 3: HARMONY INTEGRATION ──────────────────────────────────────────
print("Processing Approach 3: Harmony...")
seurat_harmony <- seurat_integrated

# Re-verify layout tracks are un-split for Harmony execution loops
DefaultAssay(seurat_harmony)<- "RNA"
seurat_harmony[["RNA"]] <- JoinLayers(seurat_harmony[["RNA"]])

seurat_harmony <- NormalizeData(seurat_harmony, verbose = FALSE)
seurat_harmony <- FindVariableFeatures(seurat_harmony, nfeatures = 2000, verbose = FALSE)
seurat_harmony <- ScaleData(seurat_harmony, verbose = FALSE)
seurat_harmony <- RunPCA(seurat_harmony, reduction.name = "pca", verbose = FALSE)

# Group specifically by our synchronized batch metadata parameter
seurat_harmony <- RunHarmony(
  object = seurat_harmony, 
  group.by.vars = "sample_id", 
  reduction.use = "pca",
  dims.use = 1:15,
  reduction.save = "harmony", 
  verbose = FALSE
)
seurat_harmony <- RunUMAP(
  object = seurat_harmony, 
  reduction = "harmony", 
  dims = 1:15, 
  reduction.name = "umap.harmony", 
  verbose = FALSE
)

# ── VISUALIZE CORRECTION STYLES SIDE-BY-SIDE ─────────────────────────────────
p1 <- DimPlot(seurat_merge, reduction = "umap.merge", group.by = "sample_id") + ggtitle("1. Uncorrected Merge") + theme_minimal()
p2 <- DimPlot(seurat_rpca, reduction = "umap.rpca", group.by = "sample_id") + ggtitle("2. Integrated RPCA") + theme_minimal()
p3 <- DimPlot(seurat_harmony, reduction = "umap.harmony", group.by = "sample_id") + ggtitle("3. Harmonized Space") + theme_minimal()

# Combined visualization dashboard matrix panel
print((p1 | p2 | p3) + plot_layout(guides = "collect"))

#===============================================================================
# ── DGE DISCOVERY AND MULTI-TAB EXCEL EXPORT (Q3.10) ──────────────────
DefaultAssay(seurat_harmony) <- "RNA"
seurat_harmony[["RNA"]] <- JoinLayers(seurat_harmony[["RNA"]])

# Lock cell identities to the true biological labels stored by the author
Idents(seurat_harmony) <- seurat_harmony$paper_annot
cell_types <- unique(as.character(Idents(seurat_harmony)))

# Initialize empty workbook
wb <- createWorkbook()

# EXPLICIT COHORT VALUES DEFINED FROM CONSOLE MATRIX
AD_label   <- "AD"
ctrl_label <- "C."

for (cell_type in cell_types) {
  print(paste("Running DEGs for:", cell_type))
  
  # Subset to current cell type
  cell_subset <- subset(seurat_harmony, idents = cell_type)
  
  # Switch group identities to the disease condition field ('group')
  Idents(cell_subset) <- cell_subset$group
  
  # Ensure both comparison groups are actively present inside this specific subset
  if(all(c(AD_label, ctrl_label) %in% Idents(cell_subset))) {
    
    deg_results <- FindMarkers(
      object = cell_subset,
      ident.1 = AD_label,
      ident.2 = ctrl_label,
      test.use = "MAST",
      min.pct = 0.10,
      logfc.threshold = 0.25,
      verbose = FALSE
    )
    
    # Format table data
    deg_results$Gene <- rownames(deg_results)
    deg_results <- deg_results %>% select(Gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2)
    
    # Excel sheet naming cleanup (max 31 characters, clean formatting)
    clean_tab_name <- substr(gsub("[^A-Za-z0-9_]", "_", cell_type), 1, 30)
    
    addWorksheet(wb, sheetName = clean_tab_name)
    writeData(wb, sheet = clean_tab_name, x = deg_results)
  }
}

# Write workbook to disk
excel_out <- "C:/Users/Hana/Desktop/Screening/data_scRNA/AD_Control_CellType_DEGs.xlsx"
saveWorkbook(wb, file = excel_out, overwrite = TRUE)

#========================================================

