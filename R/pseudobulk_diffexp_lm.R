library(presto)
library(DESeq2)
library(edgeR)

meta_data <- readRDS("meta_data_01.12.rds")
ids_ref <- readRDS("ids_ref_cca_01.12.rds")

meta_data$clus <- ids_ref$res_2.00
meta_data <- meta_data[meta_data$clus %in% names(table(meta_data$clus))[table(meta_data$clus) > 5],]

exprs_raw <- readRDS("exprs_raw_01.12.rds")
all_collapse <- collapse_counts(exprs_raw, meta_data, c("clus", "donor", "batch"))
rm(exprs_raw)
gc()

cutoff <- 30
drop <- which(apply(all_collapse$counts_mat, 1, max) < cutoff)
d <- all_collapse$counts_mat[-drop,]
dim(d)
tmp <- log2(cpm(d)+1)

group <- interaction(all_collapse$meta_data$clus)
cell_umi <- apply(all_collapse$counts_mat, 2, function(x) {log2(sum(x) + 1)})
donor <- interaction(all_collapse$meta_data$donor)
batch <- interaction(all_collapse$meta_data$batch)

mm <- model.matrix(~0 + group + batch + donor + cell_umi)
colnames(mm) <- gsub("-", ".", colnames(mm))

diffexp_lm <- function(gexp, gname, min_clus, max_clus, mm) {
    Reduce(rbind, lapply(min_clus:max_clus, function(x) {
        m <- lm(gexp~., data = data.frame(mm[,colnames(mm) == paste0("group", x) | !grepl("group", colnames(mm))]))
        m_null <- lm(gexp~., data = data.frame(mm[,!grepl("group", colnames(mm))]))

        m_chisq <- anova(m_null, m, test = "LRT")
        m_f <- anova(m_null, m, test = "F")
        m_f_var <- anova(m, test = "F")
        c(summary(m)$coefficients[2,], m_chisq$`Sum of Sq`[2], m_chisq$P[2], m_f$F[2], m_f$P[2], m_f_var$F[1], m_f_var$P[1], gname, x)
    }))
}

results <- Reduce(rbind, lapply(start:end, function(x) {
    diffexp_lm(tmp[x,], row.names(tmp)[x], 0, 30, mm)}
))
