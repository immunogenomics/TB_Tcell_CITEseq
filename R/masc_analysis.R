source("utils.R")

meta_data <- read.table("../data/meta_data.txt", sep = "\t", header = T)
ids_ref <- readRDS("../data/ids_ref_cca.rds")

meta_data$nGene <- scale(log2(meta_data$nGene))
meta_data$nUMI <- scale(log2(meta_data$nUMI))
meta_data$TB_STATUS <- factor(meta_data$TB_STATUS, levels = c("CONTROL", "CASE"))

# Set variable of interest
var_name = "TB_STATUS"

masc_df_status <- data.frame(meta_data, clus = ids_ref$res_2.00)
masc_df_status <- masc_df_status[!masc_df_status$clus %in% names(table(ids_ref$res_2.00))[table(ids_ref$res_2.00) < 5],]
masc_df_status <- masc_df_status[!is.na(masc_df_status[,var_name]),]

# Baseline model (first column in Figure 3f)

masc_res <- MASC.me(masc_df_status, factor(masc_df_status$clus),
                contrast = var_name,
                random_effects = c("donor", "batch"),
                fixed_effects = c("nUMI", "percent_mito"),
                verbose = TRUE,
                save_models = F) %>%
    dplyr::mutate(bonferroni = p.adjust(model.pvalue, method = "bonferroni")) %>%
    dplyr::arrange(model.pvalue)

# Full model (last column in Figure 3f)

masc_res <- MASC.me(masc_df_status, factor(masc_df_status$clus),
                contrast = var_name,
                random_effects = c("donor", "batch"),
                fixed_effects = c("nUMI", "percent_mito", "tbru_age", "I(tbru_age^2)", "Sex", "season_winter", "EURad4KR"),
                verbose = TRUE,
                save_models = F) %>%
    dplyr::mutate(bonferroni = p.adjust(model.pvalue, method = "bonferroni")) %>%
    dplyr::arrange(model.pvalue)

# Figure 4a (and similar volcano plots)

ggplot(data = masc_res, aes(x = statusCASE.OR, y = -log10(model.pvalue))) +
    theme_classic() + ylab("-log(p-value)") + xlab("Odds ratio case vs. control") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkblue") +
    geom_errorbarh(aes(xmin=statusCASE.OR.95pct.ci.lower, xmax=statusCASE.OR.95pct.ci.upper), col = "darkgrey") +
    geom_point() + 
    geom_hline(yintercept = -log10(.05/31), linetype = "dashed", color = "darkred") + 
    scale_x_log10(breaks = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)) + geom_label_repel(size = 3)

# Gamma p value (Figure 3f)

pgamma(-sum(log(masc_res$model.pvalue)), 31, 1, lower.tail = F))
