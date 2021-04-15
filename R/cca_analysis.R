library(singlecellmethods)
library(Matrix)
library(irlba)
library(harmony)
library(uwot)
library(dplyr)
library(parallel)
source("utils.R")

meta_data <- readRDS("meta_data_01.12.rds")
exprs_raw <- readRDS("exprs_raw_01.12.rds")

exprs_norm <- exprs_raw[,meta_data$cell_id[meta_data$qc]] %>% singlecellmethods::normalizeData(method = "log")
rm(exprs_raw);gc()
var_genes <- FindVariableGenesBatch(exprs_norm, meta_data, "cell_id", "batch", 1000, 0.1) %>% with(unique(gene))
exprs_scaled <- exprs_norm[var_genes, ] %>% singlecellmethods::scaleData()

adt_exprs_raw <- readRDS("adt_exprs_raw_01.12.rds")
adt_exprs_norm <- singlecellmethods::normalizeData(adt_exprs_raw[,meta_data$cell_id[meta_data$memTgate]], method = "cellCLR")
rm(adt_exprs_raw);gc()
var_prots <- row.names(adt_exprs_norm)[!row.names(adt_exprs_norm) %in% c("MouseIgG")]
adt_exprs_scaled <- adt_exprs_norm[var_prots,] %>% singlecellmethods::scaleData()

tcr_genes <- read.csv("../data/TCRgenes.csv")
matched.genes <- tcr_genes$IMGT.GENE.DB[grepl("T cell receptor", tcr_genes$IMGT.and.HGNC.gene.definition)]

res_cca = cc(t(exprs_scaled[!row.names(exprs_scaled) %in% matched.genes,]), t(adt_exprs_scaled))
res_cca = res_cca$scores$xscores[,1:20]

harmony_res <- HarmonyMatrix(res_cca, meta_data, c("donor","batch"), theta = c(2,2), lambda = c(1,1),
                               plot_convergence = TRUE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)
umap_res <- umap(harmony_res, n_neighbors = 30L, metric = "euclidean", min_dist = .1)

snn_ref <- BuildSNNSeurat(harmony_res, nn.eps = 0)
resolution_list <- c(2.0, 2.4)
ids_ref <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
    Seurat:::RunModularityClustering(SNN = snn_ref, modularity = 1,
        resolution = res_use, algorithm = 1, n.start = 20,
        n.iter = 20, random.seed = 100, print.output = FALSE,
        temp.file.location = NULL, edge.file.name = NULL)
}, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))
ids_ref <- data.frame(ids_ref)
rm(snn_ref)
gc()
colnames(ids_ref) <- sprintf("res_%.2f", resolution_list)
