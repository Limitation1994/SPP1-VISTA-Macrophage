rm(list = ls())

library(Seurat)
library(monocle)
library(tidyverse)
library(data.table)
library(ggplot2)
library(plyr)
library(scRNAtoolVis)
library(viridis)
library(cols4all)
library(ggpubr)
library(readxl)
library(patchwork)
library(harmony)

setwd("/home/data/t020437/ascite")

sce <- readRDS("data/ascites.rds")

table(sce$orig.ident)
table(sce$group_old)
table(sce$sample_old)
table(sce$group)

sce <- CreateSeuratObject(
  counts = sce@assays$RNA@counts,
  meta.data = sce@meta.data[, c(1:4, 11:13, 21, 22)]
)

sce <- NormalizeData(
  sce,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

sce <- FindVariableFeatures(sce)

sce <- ScaleData(sce)

sce <- RunPCA(
  sce,
  features = VariableFeatures(sce)
)

sce <- RunHarmony(
  sce,
  group.by.vars = "orig.ident"
)

sce <- RunUMAP(
  sce,
  reduction = "harmony",
  dims = 1:15
)

sce <- RunTSNE(
  sce,
  reduction = "harmony",
  dims = 1:15
)

sce <- FindNeighbors(
  sce,
  reduction = "harmony",
  dims = 1:15
)

for (res in c(0.01, 0.05, 0.1, 0.5, 1, 1.5, 2)) {
  sce <- FindClusters(
    sce,
    resolution = res,
    algorithm = 1
  )
}

dir.create(
  "1cell",
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  sce,
  file = "1cell/sce.RDS"
)

p_cluster <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "RNA_snn_res.1.5",
  dims = c(1, 2),
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  label.color = "black",
  label.box = TRUE,
  sizes.highlight = 1
)

p_cluster

mainmarkers <- c(
  "PTPRC",
  "CD3D",
  "CD3E",
  "CD3G",
  "TRAC",
  "CD4",
  "CD8A",
  "NKG7",
  "KLRF1",
  "KLRD1",
  "CLEC9A",
  "XCR1",
  "CADM1",
  "DNASE1L3",
  "CLEC10A",
  "CD1C",
  "FCER1A",
  "CD1E",
  "STMN1",
  "MKI67",
  "TOP2A",
  "TUBA1B",
  "TUBB",
  "S100B",
  "GSN",
  "TMEM97",
  "PAK1",
  "PRDM16",
  "LAMP3",
  "CCL22",
  "FSCN1",
  "MARCKSL1",
  "BASP1",
  "CSF1R",
  "CD68",
  "CD163",
  "CD14",
  "S100A8",
  "S100A9",
  "LYZ",
  "FCN1",
  "FCGR3A",
  "LST1",
  "LILRB2",
  "CD79B",
  "MS4A1",
  "CD79A",
  "JCHAIN",
  "XBP1",
  "MZB1",
  "COL3A1",
  "NNMT",
  "PDGFRA",
  "COL1A2"
)

p_markers <- DotPlot(
  sce,
  features = mainmarkers,
  group.by = "RNA_snn_res.1.5"
) +
  scale_color_continuous_c4a_seq(
    "linear_yl_mg_bu",
    reverse = FALSE
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

p_markers

ggsave(
  filename = "1cell/cluster_type_1.png",
  plot = p_markers,
  width = 15,
  height = 8
)

ggsave(
  filename = "1cell/cluster_type_1.pdf",
  plot = p_markers,
  width = 15,
  height = 8
)

Idents(sce) <- sce$RNA_snn_res.1.5

pbmc <- sce

pbmc.markers <- FindAllMarkers(
  pbmc,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers_2 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = abs(avg_log2FC),
    n = 2,
    with_ties = FALSE
  )

write.csv(
  pbmc.markers,
  file = "1cell/pbmc.markers.csv",
  row.names = FALSE
)

top_markers_5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  )

Idents(sce) <- sce$RNA_snn_res.1.5

annotation_curated_main <- read_excel(
  "1cell/curated_annotation_main1.xlsx"
)

new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(sce)

sce <- RenameIdents(
  sce,
  new_ids_main
)

sce$main_cell_type <- Idents(sce)

table(sce$main_cell_type)

saveRDS(
  sce,
  file = "1cell/sce_cell_type.Rds"
)

cols <- c(
  "#C1E6F3",
  "#6778AE",
  "#91D0BE",
  "#DCC1DD",
  "#58A4C3",
  "#CCE0F5",
  "#CCC9E6",
  "#F1BB72",
  "#F3B1A0",
  "#D6E7A3",
  "#E63863",
  "#AB3282",
  "#68A180",
  "#23452F",
  "#BD956A",
  "#8C549C",
  "#585658",
  "#9FA3A8",
  "#E0D4CA",
  "#C5DEBA",
  "#E4C755",
  "#F7F398",
  "#AA9A59",
  "#E63863",
  "#E39A35",
  "#B53E2B",
  "#712820",
  "#476D87",
  "#625D9E",
  "#68A180",
  "#3A6963",
  "#968175"
)

p_celltype <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "main_cell_type",
  dims = c(1, 2),
  label.size = 4,
  pt.size = 0.1,
  label.color = "black",
  cols = cols,
  label.box = TRUE,
  sizes.highlight = 1
)

p_celltype

ggsave(
  filename = "1cell/cell_type_umap.png",
  plot = p_celltype,
  width = 6,
  height = 5
)

ggsave(
  filename = "1cell/cell_type_umap.pdf",
  plot = p_celltype,
  width = 6,
  height = 5
)

p_class <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "main_cell_type",
  split.by = "class",
  dims = c(1, 2),
  label.size = 4,
  pt.size = 0.1,
  label.color = "black",
  cols = cols,
  label.box = TRUE,
  sizes.highlight = 1
)

p_class

ggsave(
  filename = "1cell/cell_type_umap_organ.png",
  plot = p_class,
  width = 11,
  height = 5
)

ggsave(
  filename = "1cell/cell_type_umap_organ.pdf",
  plot = p_class,
  width = 11,
  height = 5
)

p_sample <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "main_cell_type",
  split.by = "sample_old",
  dims = c(1, 2),
  label.size = 4,
  pt.size = 0.1,
  label.color = "black",
  cols = cols,
  label.box = TRUE,
  sizes.highlight = 1
)

p_sample

ggsave(
  filename = "1cell/cell_type_umap_sample.png",
  plot = p_sample,
  width = 30,
  height = 5
)

ggsave(
  filename = "1cell/cell_type_umap_sample.pdf",
  plot = p_sample,
  width = 30,
  height = 5
)

p_group <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "main_cell_type",
  split.by = "group_old",
  dims = c(1, 2),
  label.size = 4,
  pt.size = 0.1,
  label.color = "black",
  cols = cols,
  label.box = TRUE,
  sizes.highlight = 1
)

p_group

ggsave(
  filename = "1cell/cell_type_umap_cemip.png",
  plot = p_group,
  width = 11,
  height = 5
)

ggsave(
  filename = "1cell/cell_type_umap_cemip.pdf",
  plot = p_group,
  width = 11,
  height = 5
)

saveRDS(
  sce,
  file = "1cell/sce_cell_type.Rds"
)