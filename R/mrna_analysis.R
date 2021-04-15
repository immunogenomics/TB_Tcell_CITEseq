library(singlecellmethods)
library(Matrix)
library(irlba)
library(harmony)
library(uwot)
library(dplyr)
source("utils.R")

meta_data <- readRDS("meta_data_01.12.rds")
exprs_raw <- readRDS("exprs_raw_01.12.rds")

exprs_norm <- exprs_raw[,meta_data$cell_id[meta_data$qc]] %>% singlecellmethods::normalizeData(method = "log")
var_genes <- FindVariableGenesBatch(exprs_norm, meta_data, "cell_id", "batch", 1000, 0.1) %>% with(unique(gene))
exprs_cosine <- exprs_norm[var_genes, ] %>% singlecellmethods::scaleData() %>% cosine_normalize(2)
pca_res <- irlba::prcomp_irlba(t(exprs_cosine), 20)

harmony_res <- HarmonyMatrix(pca_res, meta_data, c("donor","batch"), theta = c(2,2), lambda = c(1,1),
                               plot_convergence = TRUE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)
umap_res <- umap(harmony_res, n_neighbors = 30L, metric = "euclidean", min_dist = .1)
