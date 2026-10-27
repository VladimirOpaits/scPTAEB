library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(dplyr)
library(viridis)
library(RANN)

harmony_coords = as.matrix(read.csv("data_bio/harmony_coords.csv"))
harmony_reduced = harmony_coords[, 1:10] # most of info in first dimensions
cluster_labels = read.csv("data_bio/cluster_labels.csv")$full_clustering

#downsampling
set.seed(67) 
idx_sub <- sample(nrow(harmony_reduced), 30000)
coords_sub <- harmony_reduced[idx_sub, ]
labels_sub <- cluster_labels[idx_sub]


sce_sub <- SingleCellExperiment(assays = list(counts = matrix(0, nrow = 1, ncol = nrow(coords_sub))))
reducedDim(sce_sub, "HARMONY") <- coords_sub
colData(sce_sub)$full_clustering <- labels_sub


sce_sub <- slingshot(sce_sub, clusterLabels = 'full_clustering', reducedDim = "HARMONY", start.clus = "CD14_mono", approx_points = 150) 





library(ggplot2)
library(viridis)
library(cowplot)

umap_coords <- as.matrix(read.csv("data_bio/umap_coords.csv"))

pt <- slingPseudotime(sce_sub)
crv <- slingCurves(sce_sub)

pseudotime_vec <- pt[, 1]

df_plot <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Cluster = cluster_labels,
  Pseudotime = NA
)

df_plot$Pseudotime[idx_sub] <- pseudotime_vec

curve_coords <- crv[[1]]$s[crv[[1]]$ord, ]
harmony_sub <- harmony_reduced[idx_sub, , drop = FALSE]

nn_result <- nn2(data = curve_coords, query = harmony_sub, k = 1)
nearest_cells <- nn_result$nn.idx[, 1]

umap_curve <- umap_coords[idx_sub[nearest_cells], , drop = FALSE]

curve_umap <- data.frame(
  UMAP1 = umap_curve[, 1],
  UMAP2 = umap_curve[, 2],
  PT = pseudotime_vec[nearest_cells] 
)

curve_umap <- curve_umap[order(curve_umap$PT), ]
curve_umap$UMAP1 <- zoo::rollmean(curve_umap$UMAP1, k = 5, fill = "extend")
curve_umap$UMAP2 <- zoo::rollmean(curve_umap$UMAP2, k = 5, fill = "extend")

p1 <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "Кластеры моноцитов", x = "UMAP1", y = "UMAP2")

p2 <- ggplot() +
  geom_point(data = df_plot, aes(x = UMAP1, y = UMAP2, color = Pseudotime), 
             size = 0.5, alpha = 0.6) +
  geom_path(data = curve_umap, aes(x = UMAP1, y = UMAP2), 
            color = "black", linewidth = 0.5, alpha = 0.3) +
  scale_color_viridis_c(option = "plasma", na.value = "lightgrey") +
  theme_classic() +
  labs(title = "Slingshot Pseudotime", x = "UMAP1", y = "UMAP2")

plot_grid(p1, p2, ncol = 2)