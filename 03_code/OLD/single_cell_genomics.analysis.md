



```bash
cd /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/analysis/skb2_hcontortus/all-sample/DGE_filtered


cat cell_metadata.csv | sed -e 's/hcontortus_adult_female/AdultFemale/g' -e 's/hcontortus_adult_male/AdultMale/g' -e 's/hcontortus_L3/L3/g' > cell_metadata_fixed.csv
```

# haemonchus
```R
library(pals)
library(Seurat)
library(tidyverse)
library(patchwork)


mat_path <- "/nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/analysis/skb2_hcontortus/all-sample/DGE_filtered"
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata_fixed.csv"), row.names = 1)

# Create object
data <- CreateSeuratObject(mat, min.cells = 100, names.field = 0, meta.data = cell_meta)


data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(data), 10)
# plot variable features with and without labels
#plot <- VariableFeaturePlot(data)
#plot <- LabelPoints(plot = plot, points = top10, repel = TRUE)

data <- ScaleData(data)

data <- RunPCA(data)


#plot <- DimPlot(data, reduction = "pca", group.by = "orig.ident")


data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)

data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 1.2)

data <- BuildClusterTree(data, reorder = TRUE, reorder.numeric = TRUE)


data <- RunUMAP(data, dims = 1:30)
plot <- DimPlot(data, reduction = "umap")


data <- SetIdent(data, value = "sample")



plot1 <- DimPlot(data, reduction = "umap") 
plot1 <- plot1 + scale_colour_viridis_d(begin=0.8, end=0.2, option="plasma")


data <- SetIdent(data, value = "seurat_clusters")


plot2 <- DimPlot(data, reduction = "umap")

plot1 + plot2
```


FeaturePlot(data, features = c("HCON-00085310", "HCON-00171910", "HCON-00099020")) + 
    patchwork::plot_layout(ncol = 3, nrow = 1)

ggsave("figure_hc_gene-expression_examples.pdf", height=100, width=330, units="mm")
ggsave("figure_hc_gene-expression_examples.png")






# trichuris
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/analysis/skb1_tmuris/all-sample/DGE_filtered/

cat cell_metadata.csv | sed -e 's/tmuris_adult_female/AdultFemale/g' -e 's/tmuris_adult_male/AdultMale/g' > cell_metadata_fixed.csv


/nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/analysis/skb1_tmuris/all-sample/DGE_filtered
```


```R
library(Seurat)
library(tidyverse)
library(patchwork)


mat_path <- "/nfs/users/nfs_s/sd21/lustre_link/parse_single_cell/analysis/skb1_tmuris/all-sample/DGE_filtered"
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata_fixed.csv"), row.names = 1)

# Create object
data <- CreateSeuratObject(mat, min.cells = 100, names.field = 0, meta.data = cell_meta)


data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)
# plot variable features with and without labels
plot <- VariableFeaturePlot(data)
plot <- LabelPoints(plot = plot, points = top10, repel = TRUE)

data <- ScaleData(data)

data <- RunPCA(data)


plot <- DimPlot(data, reduction = "pca", group.by = "orig.ident")


data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)

data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 1.2)

data <- BuildClusterTree(data, reorder = TRUE, reorder.numeric = TRUE)


data <- RunUMAP(data, dims = 1:30)
plot <- DimPlot(data, reduction = "umap")
plot




data <- SetIdent(data, value = "sample")
plot3 <- DimPlot(data, reduction = "umap")
plot3 <- plot3 + scale_colour_viridis_d(begin=0.8, end=0.2)


data <- SetIdent(data, value = "seurat_clusters")
plot4 <- DimPlot(data, reduction = "umap")


plot3 + plot4

```



```R

plot1 + plot2 + plot3 + plot4 + plot_layout(ncol=2)

ggsave("figure_haemonchus_trichuris_lifestage_nuclei-clusters.pdf", height=250, width=330, units="mm")
ggsave("figure_haemonchus_trichuris_lifestage_nuclei-clusters.png")
```

```R
FeaturePlot(data, features = c("Gene:WBGene00290339","Gene:WBGene00295760", "Gene:WBGene00299858" )) +
patchwork::plot_layout(ncol = 3, nrow = 1)


ggsave("figure_tm_gene-expression_examples.pdf", height=100, width=330, units="mm")
ggsave("figure_tm_gene-expression_examples.png")
```




```R
markers <- FindAllMarkers(data, only.pos = TRUE)

markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10



markers <- FindAllMarkers(data, only.pos = TRUE)

markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20



    DoHeatmap(data, features = top10$gene) + NoLegend()


FeaturePlot(data, features = "HCON_00176790")

FeaturePlot(data, features = c(
"HCON-00149900", 
"HCON-00171910", 
"HCON-00188550", 
"HCON-00184500",
"HCON-00085310",
"HCON-00164760",
"HCON-00184500",
"HCON-00057880",
"HCON-00175410",
"HCON-00032930",
"HCON-00099380",
"HCON-00113700",
"HCON-00160940",
"HCON-00137670",
"HCON-00084180",
"HCON-00024010",
"HCON-00099020",
"HCON-00114720",
"HCON-00182130",
"HCON-00168040"))


FeaturePlot(object = data,
    features = c("Gene:WBGene00300018", "Gene:WBGene00299858"))

```