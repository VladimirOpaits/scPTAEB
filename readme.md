Single cell RNA seq project to identify differential states of cells with AutoEncoder latent space


*TODO:*

-- Recalculate Slingshot with 20 Harmony components   #DONE  , much better



-- Protein ground truth

-- exclude top 10 genes and repeat experiment  # DONE, but requires to check:


> slingLineages(sce_sub)
$Lineage1
[1] "CD14_mono"      "Mono_prolif"    "CD83_CD14_mono"

$Lineage2
[1] "CD14_mono"    "C1_CD16_mono" "CD16_mono"   
> table(labels_sub, is.na(pt[, 1]))
                
labels_sub       FALSE  TRUE
  C1_CD16_mono     101   773
  CD14_mono      12927   227
  CD16_mono        183  3287
  CD83_CD14_mono 12321    35
  Mono_prolif      146     0
> slingLineages(sce_sub)
$Lineage1
[1] "CD14_mono"      "Mono_prolif"    "CD83_CD14_mono"

$Lineage2
[1] "CD14_mono"    "C1_CD16_mono" "CD16_mono"   
> colnames(pt)
[1] "Lineage1" "Lineage2"
> pt_all <- as.data.frame(slingPseudotime(sce_sub))
> pt_all$cell_id <- rownames(harmony_df)[idx_sub]
> write.csv(pt_all, "data_bio/slingshot_pseudotime_all_lineages.csv", row.names = FALSE)

