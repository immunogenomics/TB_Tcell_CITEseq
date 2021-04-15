
adt_exprs_norm <- readRDS("adt_exprs_norm_01.12.rds")
meta_data <- readRDS("meta_data_01.12.rds")

plot_dens <- get_density(adt_exprs_norm["CD3",], adt_exprs_norm["TCRab",], n = 100)
fig.size(3,4)
ggplot(data = data.frame(CD3 = adt_exprs_norm["CD3",], 
                         TCRab = adt_exprs_norm["TCRab",]),
       aes(x = CD3, y = TCRab, color = plot_dens)) +
    geom_point(shape = ".") +
    scale_color_viridis() +
    labs(color = "density")

ab_idx <- rep(1, length(meta_data$cell_id))
ab_idx[adt_exprs_norm["TCRab",] < (1.5*adt_exprs_norm["CD3",]-5.6)] <- 2
meta_data$ab <- ab_idx

cd4_idx <- rep(1, sum(ab_idx == 1))
cd4_idx[adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] > 1.5 & adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id] < .4*adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] + .75] <- 2
cd4_idx[adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] < 1.5 & adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id] > .4*adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] + .75] <- 3
cd4_idx[adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] < 1.5 & adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id] < .4*adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] + .75] <- 4
meta_data$cd4 <- 0
meta_data$cd4[meta_data$ab == 1] <- cd4_idx

model <- rpart(cl ~ ., data = data.frame(t(adt_exprs_norm), cl = factor(ids_ref$res_2.00 == 12)), cp = .000001)

cd26_idx <- rep(1, sum(cd4_idx == 2))
cd26_idx[adt_exprs_norm["CD26",subset(meta_data, cd4 == 2)$cell_id] > 2] <- 2
meta_data$cd26 <- 0
meta_data$cd26[meta_data$cd4 == 2] <- cd26_idx

cd161_idx <- rep(1, sum(cd4_idx == 2))
cd161_idx[adt_exprs_norm["CD161",subset(meta_data, cd4 == 2)$cell_id] > .125] <- 2
meta_data$cd161 <- 0
meta_data$cd161[meta_data$cd4 == 2] <- cd161_idx

ccr6_idx <- rep(1, sum(cd4_idx == 2))
ccr6_idx[adt_exprs_norm["CD196/CCR6",subset(meta_data, cd4 == 2)$cell_id] > .5] <- 2
meta_data$ccr6 <- 0
meta_data$ccr6[meta_data$cd4 == 2] <- ccr6_idx
