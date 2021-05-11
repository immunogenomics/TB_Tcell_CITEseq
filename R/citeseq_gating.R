source("utils.R")

adt_exprs_norm <- readRDS("../data/adt_exprs_norm.rds")

# Gating TCRab+CD3+
plot_dens <- get_density(adt_exprs_norm["CD3",], adt_exprs_norm["TCRab",], n = 100)
ggplot(data = data.frame(CD3 = adt_exprs_norm["CD3",], 
                         TCRab = adt_exprs_norm["TCRab",]),
       aes(x = CD3, y = TCRab, color = plot_dens)) +
    geom_point(shape = ".") +
    scale_color_viridis() +
    labs(color = "density")

ab_idx <- rep(1, length(meta_data$cell_id))
ab_idx[adt_exprs_norm["TCRab",] < (1.5*adt_exprs_norm["CD3",]-5.6)] <- 2
meta_data$ab <- ab_idx

# Gating CD4 and CD8
tmp <- data.frame(CD4 = adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id],
                  CD8 = adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id])
tmp <- subset(tmp, CD4 != 0 & CD8 != 0)
plot_dens <- get_density(tmp[,1], tmp[,2], n = 100)
ggplot(data = tmp,
       aes(x = CD4, y = CD8, color = plot_dens)) +
    geom_point(shape = ".") +
    scale_color_viridis() +
    labs(color = "density") +
    geom_vline(xintercept = 1.5, color = "darkred") +
    geom_hline(yintercept = 1.75, color = "darkred") +
    geom_hline(yintercept = 1, color = "darkred")

cd4_idx <- rep(1, sum(ab_idx == 1))
cd4_idx[adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] > 1.5 & adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id] < 1.75] <- 2
cd4_idx[adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] < 1.5 & adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id] > 1] <- 3
cd4_idx[adt_exprs_norm["CD4.1",subset(meta_data, ab == 1)$cell_id] < 1.5 & adt_exprs_norm["CD8a",subset(meta_data, ab == 1)$cell_id] < 1] <- 4
meta_data$cd4 <- 0
meta_data$cd4[meta_data$ab == 1] <- cd4_idx

# Gating central memory
cd4_cd62l_idx <- rep(1, sum(meta_data$cd4 == 2))
cd4_cd62l_idx[adt_exprs_norm["CD45RO",subset(meta_data, cd4 == 2)$cell_id] > 0.75 & adt_exprs_norm["CD62L",subset(meta_data, cd4 == 2)$cell_id] < .6] <- 2
cd4_cd62l_idx[adt_exprs_norm["CD45RO",subset(meta_data, cd4 == 2)$cell_id] < 0.75 & adt_exprs_norm["CD62L",subset(meta_data, cd4 == 2)$cell_id] > .6] <- 3
cd4_cd62l_idx[adt_exprs_norm["CD45RO",subset(meta_data, cd4 == 2)$cell_id] < 0.75 & adt_exprs_norm["CD62L",subset(meta_data, cd4 == 2)$cell_id] < .6] <- 4
meta_data$cd4_cd62l <- 0
meta_data$cd4_cd62l[meta_data$cd4 == 2] <- cd4_cd62l_idx

cd8_cd62l_idx <- rep(1, sum(meta_data$cd4 == 3))
cd8_cd62l_idx[adt_exprs_norm["CD45RO",subset(meta_data, cd4 == 3)$cell_id] > 1.25 & adt_exprs_norm["CD62L",subset(meta_data, cd4 == 3)$cell_id] < .6] <- 2
cd8_cd62l_idx[adt_exprs_norm["CD45RO",subset(meta_data, cd4 == 3)$cell_id] < 1.25 & adt_exprs_norm["CD62L",subset(meta_data, cd4 == 3)$cell_id] > .6] <- 3
cd8_cd62l_idx[adt_exprs_norm["CD45RO",subset(meta_data, cd4 == 3)$cell_id] < 1.25 & adt_exprs_norm["CD62L",subset(meta_data, cd4 == 3)$cell_id] < .6] <- 4
meta_data$cd8_cd62l <- 0
meta_data$cd8_cd62l[meta_data$cd4 == 3] <- cd8_cd62l_idx

# Define populations
totalseq_prop <- data.frame(T_memory_gd_memory = table(factor(subset(meta_data, ab == 2)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_memory = table(factor(subset(meta_data, ab == 1)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_CD4_memory = table(factor(subset(meta_data, cd4 == 2)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_CD8_memory = table(factor(subset(meta_data, cd4 == 3)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_CD4_CM_memory = table(factor(subset(meta_data, cd4_cd62l == 1 | cd4_cd62l == 3)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_CD4_EM_memory = table(factor(subset(meta_data, cd4_cd62l == 2 | cd4_cd62l == 4)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_CD8_CM_memory = table(factor(subset(meta_data, cd8_cd62l == 1 | cd4_cd62l == 3)$id, levels = names(table(meta_data$id))))/table(meta_data$id),
			T_memory_ab_CD8_EM_memory = table(factor(subset(meta_data, cd8_cd62l == 2 | cd4_cd62l == 4)$id, levels = names(table(meta_data$id))))/table(meta_data$id))

totalseq_prop[is.na(totalseq_prop)] <- 0
colnames(totalseq_prop) <- gsub("\\.Freq", "", colnames(totalseq_prop))
row.names(totalseq_prop) <- names(table(meta_data$id))

flow_prop <- read.csv("../data/10.29_UpdatedFlowProportionsWithAbundances.csv")
flow_prop$batch <- Reduce(c, lapply(1:46, function(x) {rep(x, 6)}))
row.names(flow_prop) <- paste(flow_prop$sample_ID, flow_prop$batch, sep = "-")

flow_prop_sub <- flow_prop[row.names(totalseq_prop),]
totalseq_prop <- totalseq_prop[!is.na(flow_prop_sub$sample_ID),]
flow_prop_sub <- flow_prop_sub[!is.na(flow_prop_sub$sample_ID),]

# Extended Data Figure 2a
tmp = c("T_memory_ab_memory", "T_memory_gd_memory", "T_memory_ab_CD4_memory", "T_memory_ab_CD8_memory", "T_memory_ab_CD4_CM_memory", "T_memory_ab_CD4_EM_memory", "T_memory_ab_CD8_CM_memory", "T_memory_ab_CD8_EM_memory")
tmp_label = c("ab memory", "gd memory", "ab CD4+ memory", "ab CD8+ memory", "ab CD4+ central memory", "ab CD4+ effector memory", "ab CD8+ central memory", "ab CD8+ effector memory")
ggplot(data = data.frame(totalseq = 100*colMeans(totalseq_prop)[tmp],
           flow = 100*colMeans(flow_prop_sub[,-c(1:80)])[tmp]),
       aes(x = flow, y = totalseq, label = tmp_label)) +
    geom_point() +
    scale_x_log10(limits = c(0.02,100)) + scale_y_log10(limits = c(0.02,100)) +
    geom_text_repel(size = 2.5) +
    geom_abline(slope = 1, color = "darkred", linetype = "dashed") +
    xlab("Mean proportion in flow (log10 scaled)") + ylab("Mean proportion in Total-seq (log10 scaled)")

# Extended Data Figure 2b
tmp = c("T_memory_ab_memory", "T_memory_gd_memory", "T_memory_ab_CD4_memory", "T_memory_ab_CD4_CM_memory", "T_memory_ab_CD4_EM_memory", "T_memory_ab_CD8_memory", "T_memory_ab_CD8_CM_memory", "T_memory_ab_CD8_EM_memory")
Reduce(`+`, lapply(tmp, function(x) {
    ggplot(data = data.frame(totalseq = totalseq_prop[,x],
                         flow = flow_prop_sub[,x]),
       aes(flow, totalseq)) +
    geom_point(size = .5) +
    geom_abline(slope = 1, linetype = "dashed") +
    xlim(0,1) + ylim(0,1) +
    theme(axis.title = element_blank()) +
    labs(title = paste0("r = ", round(cor(log10(flow_prop_sub[,x] + .0001), log10(totalseq_prop[,x] + .0001), method = "pearson"), digits = 2)))})) + 
    plot_layout(ncol = 3)
