source("utils.R")

# Load data

meta_data <- read.table("../data/meta_data.txt", sep = "\t", header = T) # From GEO
exprs_raw <- readRDS("../data/exprs_raw.rds")
adt_exprs_raw <- readRDS("../data/adt_exprs_raw.rds")

# Normalize and scale mRNA

exprs_norm <- exprs_raw %>% singlecellmethods::normalizeData(method = "log")
rm(exprs_raw);gc()
var_genes <- FindVariableGenesBatch(exprs_norm, meta_data, "cell_id", "batch", 1000, 0.1) %>% with(unique(gene))
exprs_scaled <- exprs_norm[var_genes, ] %>% singlecellmethods::scaleData()

# Normalize and scale protein

adt_exprs_norm <- singlecellmethods::normalizeData(adt_exprs_raw, method = "cellCLR")
rm(adt_exprs_raw);gc()
var_prots <- row.names(adt_exprs_norm)[!row.names(adt_exprs_norm) %in% c("MouseIgG")]
adt_exprs_scaled <- adt_exprs_norm[var_prots,] %>% singlecellmethods::scaleData()

# Dimensionality reduction with CCA

tcr_genes <- read.csv("../data/TCRgenes.csv") # From IMGT
matched.genes <- tcr_genes$IMGT.GENE.DB[grepl("T cell receptor", tcr_genes$IMGT.and.HGNC.gene.definition)]

res_cca = cc(t(exprs_scaled[!row.names(exprs_scaled) %in% matched.genes,]), t(adt_exprs_scaled))
res_cca = res_cca$scores$xscores[,1:20]

# Batch correction

harmony_res <- HarmonyMatrix(res_cca, meta_data, c("donor","batch"), theta = c(2,2), lambda = c(1,1),
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
