source("utils.R")

# Load data

meta_data <- read.table("../data/meta_data.txt", sep = "\t", header = T)
ids_ref <- readRDS("../data/ids_ref.rds")

# Bulk RNA-seq processing

bulk_exprs <- read.csv("../data/bulk_exprs.csv", check.names = F)
bulk_exprs <- bulk_exprs %>% group_by(hgnc_symbol) %>% summarise('M0022318-1' = sum(`M0022318-1`),
                                              'M0007191-1' = sum(`M0007191-1`),
                                              'M0017830-2' = sum(`M0017830-2`),
                                              'M0017902-9' = sum(`M0017902-9`),
                                              'M0007229-9' = sum(`M0007229-9`),
                                              'M0023596-1' = sum(`M0023596-1`),
                                              'M0010677-4' = sum(`M0010677-4`),
                                              'M0021395-0' = sum(`M0021395-0`),
                                              'M0011249-1' = sum(`M0011249-1`),
                                              'M0018532-3' = sum(`M0018532-3`),
                                              'M0014787-7' = sum(`M0014787-7`),
                                              'M0012464-5' = sum(`M0012464-5`),
                                              'M0014094-8' = sum(`M0014094-8`),
                                              'M0014768-7' = sum(`M0014768-7`),
                                              'M0012131-0' = sum(`M0012131-0`))

bulk_exprs <- data.frame(bulk_exprs[-which(bulk_exprs$hgnc_symbol == "" | is.na(bulk_exprs$hgnc_symbol)),], check.names = F)
row.names(bulk_exprs) <- unlist(bulk_exprs[,1])
bulk_exprs <- log2(bulk_exprs[,-1]+1)

gene_summary <- data.frame(
  mean    = rowMeans(bulk_exprs),
  sd      = apply(bulk_exprs,1,sd),
  samples = rowSums(bulk_exprs > 0),
  gene = rownames(bulk_exprs)
)
gene_summary$mean_quantile <- rank(gene_summary$mean) / nrow(gene_summary)
gene_summary$sd_quantile <- rank(gene_summary$sd) / nrow(gene_summary)

ind <- gene_summary$mean == 0 
gene_summary_no_zero <- gene_summary[!ind,]
genes <- rownames(gene_summary_no_zero[which(gene_summary_no_zero$mean > 2),])
gene_summary_final <- gene_summary_no_zero[genes,]

bulk_exprs_final <- bulk_exprs[row.names(bulk_exprs) %in% rownames(gene_summary_final),]

# CITE-seq data processing

exprs_raw <- readRDS("../data/exprs_raw.rds")
all_collapse <- collapse_counts(exprs_raw, meta_data, c("id"))
cutoff <- 30
d <- all_collapse$counts_mat[-which(apply(all_collapse$counts_mat, 1, max) < cutoff),]
colnames(d) <- all_collapse$meta_data$id
d <- d[,order(colnames(d))]

# Berry data processing (for Extended Data Figures 8b and c)

berry1_meta <- fread("GSE19439_series_matrix.txt.gz", fill = T) # From GEO
berry1_data <- berry1_meta[74:(nrow(berry1_meta)-1),]
colnames(berry1_data) <- unlist(berry1_meta[73,])
berry1_data <- data.frame(berry1_data)
berry1_data[,2:ncol(berry1_data)] <- apply(berry1_data[,2:ncol(berry1_data)], 2, as.numeric)
berry1_data$gene <- data.frame(Gene=unlist(mget(x = unlist(berry1_meta[74:(nrow(berry1_meta)-1),1]),envir = illuminaHumanv3SYMBOL)))$Gene
berry1_data <- berry1_data[!is.na(berry1_data$gene),]
berry1_data <- berry1_data[!duplicated(berry1_data$gene),]
row.names(berry1_data) <- berry1_data$gene
berry1_data <- berry1_data[,-c(1,ncol(berry1_data))]

# Berry data processing (for Extended Data Figure 8d)

berry2_meta <- fread("GSE19435_series_matrix.txt.gz", fill = T) # From GEO
berry2_data <- berry2_meta[65:(nrow(berry2_meta)-1),]
colnames(berry2_data) <- unlist(berry2_meta[64,])
berry2_data <- data.frame(berry2_data)
berry2_data[,2:ncol(berry2_data)] <- apply(berry2_data[,2:ncol(berry2_data)], 2, as.numeric)
berry2_data$gene <- data.frame(Gene=unlist(mget(x = unlist(berry2_meta[65:(nrow(berry2_meta)-1),1]),envir = illuminaHumanv3SYMBOL)))$Gene
berry2_data <- berry2_data[!is.na(berry2_data$gene),]
berry2_data <- berry2_data[!duplicated(berry2_data$gene),]
row.names(berry2_data) <- berry2_data$gene
berry2_data <- berry2_data[,-c(1,ncol(berry2_data))]

# Scriba data processing

scriba_ea <- data.frame(fread("GSE103147_Tcells-EA-rawCounts_GEO.txt.gz"), check.names = F) # From GEO
scriba_ea <- scriba_ea[!is.na(scriba_ea$symbol),]
row.names(scriba_ea) <- scriba_ea$symbol
scriba_ea <- scriba_ea[,-c(1:2)]
scriba_ea <- log2(edgeR::cpm(scriba_ea)+1)

scriba_bgi <- data.frame(fread("GSE103147_TcellsMono-BGI-rawCounts_GEO.txt.gz"), check.names = F) # From GEO
scriba_bgi <- scriba_bgi[!is.na(scriba_bgi$symbol),]
scriba_bgi <- scriba_bgi[!duplicated(scriba_bgi$symbol),]
row.names(scriba_bgi) <- scriba_bgi$symbol
scriba_bgi <- scriba_bgi[,-c(1:2)]
scriba_bgi <- log2(edgeR::cpm(scriba_bgi)+1)

scriba <- data.frame(fread("GSE103147_series_matrix.txt.gz", fill = T), check.names = F) # From GEO
scriba <- t(scriba[-c(1:32,80:82),])
colnames(scriba) <- scriba[1,]
scriba <- scriba[-1,]
row.names(scriba) <- scriba[,1]
scriba <- scriba[,-1]

scriba_meta <- scriba[c(colnames(scriba_ea), colnames(scriba_bgi)),]
keep_scriba <- scriba_meta[,14] != "group: Not a PP case" & scriba_meta[,10] == "stimulation: unstim" & scriba_meta[,9] == "cell type: Tcells"
scriba_meta <- scriba_meta[keep_scriba,]

# Cross validation within CITE-seq T cell and bulk PBMC data

i = 14
lambdas <- 10^seq(5, -5, by = -.1)
var_genes <- unlist(read.csv("../data/var_genes.csv"))
c_prop <- table(factor(meta_data$id[ids_ref$res_2.00 == i], levels = names(table(meta_data$id))))/table(meta_data$id)
row.names(c_prop) <- substring(row.names(c_prop), 1, 10)
select_genes = intersect(intersect(var_genes, row.names(scriba_ea)), row.names(bulk_exprs_final))

opt_lambdas <- sapply(1:500, function(x){
	cv.glmnet(rbind(scale(t(d[select_genes,])),
                        scale(t(bulk_exprs_final[select_genes,]))), 
                  c(as.vector(c_prop), 
		  as.vector(c_prop[colnames(bulk_exprs_final)])), alpha = 0, lambda = lambdas)$lambda.min
	}
)

opt_lambda <- names(sort(table(opt_lambdas), decreasing = T))[1]
#opt_lambda <- 0.00199526231496888

#select_genes = intersect(intersect(var_genes, row.names(scriba_ea)), row.names(bulk_exprs_final))
select_genes = intersect(select_genes, row.names(scriba_bgi))

idx <- sample(c(rep(1:10, 28), 1, 2, 3, 4, 5, 6), 286)
y_predicted_test_cv <- c()
c_prop_test_cv <- c()
names_cv <- c()
for(x in 1:10) {
    fit_test <- glmnet(rbind(scale(t(d[select_genes,])),
                             scale(t(bulk_exprs_final[select_genes,])))[idx != x,], 
                     c(as.vector(c_prop), as.vector(c_prop[colnames(bulk_exprs_final)]))[idx != x], alpha = 0, lambda = lambdas)
    y_predicted_test <- predict(fit_test, s = opt_lambda, newx = rbind(scale(t(d[select_genes,])),
                             scale(t(bulk_exprs_final[select_genes,])))[idx == x,])
    y_predicted_test_cv <- c(y_predicted_test_cv, y_predicted_test)
    c_prop_test_cv <- c(c_prop_test_cv, c(as.vector(c_prop), as.vector(c_prop[colnames(bulk_exprs_final)]))[idx == x])
    names_cv <- c(names_cv, colnames(cbind(tmp[select_genes,], exprs_final[select_genes,]))[idx == x])
}

names(y_predicted_test_cv) <- names_cv
cor(y_predicted_test_cv, c_prop_test_cv)

# bulk PBMC correlation

cor(y_predicted_test_cv[nchar(names(y_predicted_test_cv)) != 10], c_prop_test_cv[nchar(names(y_predicted_test_cv)) != 10])

# single-cell T cell correlation

cor(y_predicted_test_cv[nchar(names(y_predicted_test_cv)) == 10], c_prop_test_cv[nchar(names(y_predicted_test_cv)) == 10])

# Extended Data Figure 8a

options(repr.plot.height = 8, repr.plot.width = 8)
ggplot(data.frame(y = y_predicted_test_cv,
                  x = c_prop_test_cv,
                  tech = factor(ifelse(nchar(names(y_predicted_test_cv)) == 10, 2, 1))) %>% arrange(tech), aes(x, y, color = tech)) +
    geom_point(size = 3) + geom_abline(slope = 1, intercept = 0) + xlim(-.005, 0.09) + ylim(-.005, 0.09) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), legend.position = "NA") + scale_color_manual(values = c("black", "deepskyblue3"))+
    xlab("C-12 proportion") + ylab("Predicted C-12 proportion")

opt_lambdas <- sapply(1:100, function(x){cv.glmnet(scale(t(tmp[intersect(var_genes, row.names(berry1_data)),])), as.vector(c_prop), alpha = 0, lambda = lambdas)$lambda.min})
#opt_lambda_berry1 <- 0.00398107170553497

# Berry analysis (1): Train model on CITE-seq T cell and bulk PBMC data

select_genes = intersect(intersect(var_genes, row.names(berry1_data)), row.names(bulk_exprs_final))
fit_berry1 <- glmnet(rbind(scale(t(d[select_genes,])),
                             scale(t(bulk_exprs_final[select_genes,]))),
                     c(as.vector(c_prop), as.vector(c_prop[colnames(bulk_exprs_final)])), alpha = 0, lambda = lambdas)

# Berry analysis (1): Apply model to Berry data

y_predicted_berry1 <- predict(fit_berry1, s = opt_lambda, newx = scale(t(berry1_data[select_genes, ])))

# Berry analysis (2): Train model on CITE-seq T cell and bulk PBMC data

fit_berry2 <- glmnet(rbind(scale(t(d[select_genes,])),
                             scale(t(bulk_exprs_final[select_genes,]))),
                     c(as.vector(c_prop), as.vector(c_prop[colnames(bulk_exprs_final)])), alpha = 0, lambda = lambdas)
y_predicted_berry2 <- predict(fit_berry2, s = opt_lambda, newx = scale(t(berry2_data[select_genes, ])))

# Scriba analysis: Train model on CITE-seq T cell and bulk PBMC data

fit_ea_bgi <- glmnet(rbind(scale(t(d[select_genes,])),
                             scale(t(bulk_exprs_final[select_genes,]))),
                     c(as.vector(c_prop), as.vector(c_prop[colnames(bulk_exprs_final)])), alpha = 0, lambda = lambdas)

# Scriba analysis: Apply model to Scriba data

y_predicted_ea_bgi <- predict(fit_ea_bgi, s = opt_lambda,
                              newx = scale(t(cbind(scriba_ea[select_genes,], scriba_bgi[select_genes,])[, keep_scriba])))

tmp_lmer_ea_bgi <- data.frame(X1 = y_predicted_ea_bgi,
                       status = scriba_meta[,14],
                       timept = as.numeric(gsub(".*: ", "", scriba_meta[,12])),
                       age = as.numeric(gsub(".*: ", "", scriba_meta[,15])),
                       sex = gsub(".*: ", "", scriba_meta[,16]),
                       donor = scriba_meta[,11],
                       dataset = factor(c(rep("ea", 256), rep("bgi", 24))))

tmp_tbru <- data.frame(meta_data[!duplicated(meta_data$id),])
row.names(tmp_tbru) <- tmp_tbru$id

# Downsampling for power estimation (Extended Data Figure 8e)

tbru_ds_true <- t(sapply(1:1000, function(x){
    x1 = sample(as.vector(c_prop)[tmp_tbru[names(c_prop), "TB_STATUS"] == "CASE"], 7)
    x2 = sample(as.vector(c_prop)[tmp_tbru[names(c_prop), "TB_STATUS"] == "CONTROL"], 12)
    c(t.test(x1, x2)$p.value, mean(x1)/mean(x2))}
))

options(repr.plot.height = 6, repr.plot.width =  10)
ggplot(data.frame(tbru_ds_true), aes(X1)) + 
    geom_histogram(bins = 100, color = "black", fill = "darkgrey") +
    xlab("P-value") + ylab("Count") + geom_vline(xintercept = 0.05, linetype = "dashed", color = "darkred") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size =  16))

# Permutation for p value

donor_status <- data.frame(donor = tmp_lmer_ea_bgi[["donor"]][!duplicated(tmp_lmer_ea_bgi[["donor"]])],
                           status = tmp_lmer_ea_bgi[["status"]][!duplicated(tmp_lmer_ea_bgi[["donor"]])])

perm_beta <- c()
perm_t <- c()
for(x in 1:10000) {
    perm_lmer <- lmer(X1~perm_status+age+sex+dataset+(1|donor), donor_status %>% mutate(perm_status = sample(status)) %>% inner_join(., (tmp_lmer_ea_bgi  %>% filter(timept %in% c(360,540))), by = c("donor", "status")))
    perm_beta <- c(perm_beta, fixef(perm_lmer)[2])
    perm_t <- c(perm_t, t.stat(perm_lmer)[2])
}
sum(perm_beta > 0.0031365)/10000

# Extended Data Figure 8f

ggplot(tmp_lmer_ea_bgi %>% filter(timept %in% c(540, 360)), aes(x = factor(-timept), y = X1, fill = status)) + 
    geom_boxplot() +
    xlab("") + ylab("Predicted C-12 proportion") +
    theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 16), 
          axis.title = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_blank()) +
    scale_x_discrete(labels = c("540 days\n(n = 47)", "360 days\n(n = 51)")) +
    scale_fill_manual(values = c("firebrick2", "dodgerblue3")) + ylim(0, .055) 


# Data for Figure 4e

summary(lmerTest::lmer(X1~status+age+sex+dataset+(1|donor), tmp_lmer_ea_bgi[tmp_lmer_ea_bgi$timept %in% c(360,540),]))
summary(lmerTest::lmer(as.vector(y_predicted_test_cv[row.names(tmp_tbru)])~TB_STATUS+tbru_age+Sex+(1|donor), tmp_tbru))

# Extended Data Figure 8g

tmp_all_early <- data.frame(prop = c(tmp_lmer_ea_bgi[tmp_lmer_ea_bgi$timept %in% c(360,540),][["X1"]], as.vector(y_predicted_test_cv[row.names(tmp_tbru)])), 
                      dataset = c(rep("Scriba", nrow(tmp_lmer_ea_bgi[tmp_lmer_ea_bgi$timept %in% c(360,540),])), rep("TBRU", nrow(tmp_tbru))),
                      status = c(gsub(".*: ", "", tmp_lmer_ea_bgi[tmp_lmer_ea_bgi$timept %in% c(360,540),][["status"]]), tolower(tmp_tbru[names(y_predicted_test_cv[row.names(tmp_tbru)]),"TB_STATUS"])))

ggplot(tmp_all_early, aes(x = dataset, y = prop, fill = status)) + 
    geom_boxplot() +
    xlab("") + ylab("Predicted C-12 proportion") +
    theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 22), 
          axis.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_blank()) +
    scale_x_discrete(labels = c("Pre-disease\n(Scriba 2017)", "Post-disease\nsteady state\n(TBRU LIMAA)")) +
    scale_fill_manual(values = c("firebrick2", "dodgerblue3"))
