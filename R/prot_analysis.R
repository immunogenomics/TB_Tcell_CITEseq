source("utils.R")

# Load data

meta_data <- read.table("../data/meta_data.txt", sep = "\t", header = T)
adt_exprs_raw <- readRDS("../data/adt_exprs_raw.rds")

# Normalize and scale
adt_exprs_norm <- singlecellmethods::normalizeData(adt_exprs_raw, method = "cellCLR")
var_prots <- row.names(adt_exprs_norm)[!row.names(adt_exprs_norm) %in% c("MouseIgG")]
adt_exprs_scaled <- adt_exprs_norm[var_prots,] %>% singlecellmethods::scaleData()

# Dimensionality reduction

pca_res <- irlba::prcomp_irlba(t(adt_exprs_scaled), 20)

# Batch correction

harmony_res <- HarmonyMatrix(pca_res, meta_data, c("donor","batch"), theta = c(2,2), lambda = c(1,1),
                               plot_convergence = TRUE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)
umap_res <- umap(harmony_res, n_neighbors = 30L, metric = "euclidean", min_dist = .1)

# Clustering

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

# Mutual information comparison with CCA clusters

ids_ref_cca <- readRDS("../data/ids_ref_cca.rds")
mutinformation(ids_ref_cca$res_2.00, ids_ref$res_2.00)
