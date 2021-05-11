source("utils.R")

# Figure 1c

ggplot(meta_raw, aes(x = percent_mito, y = nGene)) +
    theme_classic() +
    geom_point(color = "lightgrey",
               data = meta_raw[meta_raw$percent_mito >= .2 | meta_raw$nGene <= 500,], size = .3, aes(percent_mito, nGene)) +
    geom_hline(yintercept = 500, linetype = "dashed", color = "darkred") +
    geom_vline(xintercept = .2, linetype = "dashed", color = "darkred") +
    xlab("Percent mitochondrial UMIs per cell") + ylab("Number of genes per cell") +
    theme(legend.position = "NA")

# Extended Data Figure 1a

ggplot(data.frame(CD3=adt_exprs_norm["CD3",], CD45RO=adt_exprs_norm["CD45RO",]), aes(CD3, CD45RO)) +
    geom_point(shape = ".") +
    theme_classic() +
    geom_hline(yintercept = 0.75, color = "darkred", linetype = "dashed") +
    theme(legend.position = "NA") +
    geom_vline(xintercept = 1.75, color = "darkred", linetype = "dashed") +
    scale_color_viridis()

# Step 5 and 6 from Figure 1b

meta_raw$correct <- !meta_raw$donor %in% c("M0013046-9", "M0013353-9", "M0010790-5", "M0014748-9", "M0008352-8")
meta_raw$memTgate <- meta_raw$correct & adt_exprs_norm["CD3",] > 1.75 & adt_exprs_norm["CD45RO",] > 0.75
meta_raw$memTgate[meta_raw$correct] <- meta_raw$memTgate[meta_raw$correct] & umap_res[,1] < 8 & umap_res[,2] < 10
