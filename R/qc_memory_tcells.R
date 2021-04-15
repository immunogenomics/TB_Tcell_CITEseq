

ggplot(data.frame(CD3=adt_exprs_norm["CD3",], CD45RO=adt_exprs_norm["CD45RO",]), aes(CD3, CD45RO, color = plot_dens)) +
    geom_point(shape = ".") +
    theme_classic() +
    geom_hline(yintercept = 0.75, color = "darkred", linetype = "dashed") +
    theme(legend.position = "NA") +
    geom_vline(xintercept = 1.75, color = "darkred", linetype = "dashed") +
    scale_color_viridis()

hist(adt_exprs_norm["CD45RO",], breaks = 80) + abline(v = 0.75)
hist(adt_exprs_norm["CD3",], breaks = 80) + abline(v = 1.75)

meta_data$correct <- !meta_data$donor %in% c("M0013046-9", "M0013353-9", "M0010790-5", "M0014748-9", "M0008352-8")
meta_data$memTgate <- meta_data$correct & adt_exprs_norm["CD3",] > 1.75 & adt_exprs_norm["CD45RO",] > 0.75
meta_data$memTgate[meta_data$correct] <- meta_data$memTgate[meta_data$correct] & umap_res[,1] < 8 & umap_res[,2] < 10
