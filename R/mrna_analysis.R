source("utils.R")

# Load data

meta_data <- read.table("../data/meta_data.txt", sep = "\t", header = T)
exprs_raw <- readRDS("../data/exprs_raw.rds")

# Normalize and scale

exprs_norm <- exprs_raw %>% singlecellmethods::normalizeData(method = "log")
var_genes <- FindVariableGenesBatch(exprs_norm, meta_data, "cell_id", "batch", 1000, 0.1) %>% with(unique(gene))
exprs_cosine <- exprs_norm[var_genes, ] %>% singlecellmethods::scaleData() %>% cosine_normalize(2)

# Dimensionality reduction

pca_res <- irlba::prcomp_irlba(t(exprs_cosine), 20)

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
