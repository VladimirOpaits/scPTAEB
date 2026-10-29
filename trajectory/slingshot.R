library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(dplyr)
library(viridis)
library(RANN)
library(cowplot)
library(zoo)

harmony_df <- read.csv("data_bio/harmony_coords_mono_authors.csv", row.names = 1)
harmony_reduced <- as.matrix(harmony_df[, 1:20]) 
labels_df <- read.csv("data_bio/cluster_labels_mono.csv")
umap_df <- read.csv("data_bio/covid_mono_umap.csv")

stopifnot(identical(rownames(harmony_df), labels_df$cell_id))
stopifnot(identical(rownames(harmony_df), umap_df[[1]]))

cluster_labels <- labels_df$full_clustering
umap_coords <- as.matrix(umap_df[, c("UMAP1", "UMAP2")])

set.seed(67) 
idx_sub <- sample(nrow(harmony_reduced), 30000)
coords_sub <- harmony_reduced[idx_sub, ]
labels_sub <- cluster_labels[idx_sub]

sce_sub <- SingleCellExperiment(assays = list(counts = matrix(0, nrow = 1, ncol = nrow(coords_sub))))
reducedDim(sce_sub, "HARMONY") <- coords_sub
colData(sce_sub)$full_clustering <- labels_sub

sce_sub <- slingshot(sce_sub, clusterLabels = 'full_clustering', 
                     reducedDim = "HARMONY", start.clus = "CD14_mono", approx_points = 150) 

umap_sub <- umap_coords[idx_sub, ]
colnames(umap_sub) <- c("UMAP1", "UMAP2")
reducedDim(sce_sub, "UMAP") <- umap_sub

emb <- embedCurves(sce_sub, "UMAP")
curves <- do.call(rbind, lapply(seq_along(slingCurves(emb)), function(i) {
  cv <- slingCurves(emb)[[i]]
  d <- as.data.frame(cv$s[cv$ord, ])
  colnames(d) <- c("UMAP1", "UMAP2")
  d$Lineage <- paste0("Lineage", i)
  d
}))

colnames(sce_sub) <- rownames(harmony_df)[idx_sub]

pt <- slingPseudotime(sce_sub)

df <- data.frame(
  UMAP1    = umap_sub[, 1],
  UMAP2    = umap_sub[, 2],
  Cluster  = labels_sub,
  Lineage1 = pt[, "Lineage1"],
  Lineage2 = pt[, "Lineage2"]
)

write.csv(curves, "data_bio/slingshot_curves_author_umap.csv", row.names = FALSE)

cent <- aggregate(cbind(UMAP1, UMAP2) ~ Cluster, df, median)
arr  <- arrow(length = unit(0.25, "cm"), type = "closed")

p_cl <- ggplot(df, aes(UMAP1, UMAP2, color = Cluster)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_path(data = curves, aes(UMAP1, UMAP2, group = Lineage),
            color = "black", linewidth = 1, arrow = arr, inherit.aes = FALSE) +
  geom_text(data = cent, aes(UMAP1, UMAP2, label = Cluster),
            color = "black", size = 3, fontface = "bold", inherit.aes = FALSE) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Clusters + lineages")

plot_pt <- function(lname) {
  d <- df[order(!is.na(df[[lname]])), ]
  ggplot(d, aes(UMAP1, UMAP2, color = .data[[lname]])) +
    geom_point(size = 0.4, alpha = 0.6) +
    geom_path(data = subset(curves, Lineage == lname), aes(UMAP1, UMAP2),
              color = "black", linewidth = 1, arrow = arr, inherit.aes = FALSE) +
    scale_color_viridis_c(option = "plasma", na.value = "grey88") +
    theme_classic() + labs(title = paste(lname, "pseudotime"), color = "PT")
}

p_all <- plot_grid(p_cl, plot_pt("Lineage1"), plot_pt("Lineage2"), ncol = 3)
print(p_all)

ggsave("data_bio/slingshot_author_umap.png", p_all, width = 18, height = 7, dpi = 200)

pt_out <- as.data.frame(slingPseudotime(sce_sub))
pt_out$cell_id <- colnames(sce_sub)
write.csv(pt_out, "data_bio/slingshot_pt_authors.csv", row.names = FALSE)

lin <- slingLineages(sce_sub)
writeLines(sapply(names(lin), function(n) paste(n, paste(lin[[n]], collapse = " > "))),
           "data_bio/slingshot_lineages_authors.txt")
print(lin)

curves20 <- do.call(rbind, lapply(seq_along(slingCurves(sce_sub)), function(i) {
  cv <- slingCurves(sce_sub)[[i]]
  d <- as.data.frame(cv$s[cv$ord, ])
  d$Lineage <- paste0("Lineage", i)
  d$order <- seq_len(nrow(d))
  d
}))
write.csv(curves20, "data_bio/slingshot_curves_author_20D.csv", row.names = FALSE)