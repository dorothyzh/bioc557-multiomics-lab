/* R equivalents for each chapter. The full script is r/walkthrough.R; data are in r/data/. */
window.RCODE = {
  study: `# Packages: Bioconductor + CRAN
# BiocManager::install(c("MOFA2", "limma", "fgsea", "msigdbr"))
library(data.table); library(MOFA2); library(limma)

ds <- "human_luad"   # or "honeybee"
samples <- fread(sprintf("r/data/%s_samples.csv", ds))
read_view <- function(file) { x <- fread(file); m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m }
views <- if (ds == "human_luad") list(
  RNA     = read_view("r/data/human_luad_rna.csv"),
  Protein = read_view("r/data/human_luad_protein.csv"),
  Phospho = read_view("r/data/human_luad_phospho.csv")
) else list(
  Head    = read_view("r/data/honeybee_head.csv"),
  Abdomen = read_view("r/data/honeybee_abdomen.csv"),
  Gut     = read_view("r/data/honeybee_gut.csv")
)
# align every view to the same sample order BY ID, never by position
views <- lapply(views, function(m) m[, samples$id])
sapply(views, dim)`,
  quality: `# Missingness per sample and per feature
miss_sample  <- sapply(views, function(m) colMeans(is.na(m)))
miss_feature <- lapply(views, function(m) rowMeans(is.na(m)))

# Is missingness related to abundance? (left-censoring / MNAR in proteomics)
m <- views[[2]]
plot(rowMeans(m, na.rm = TRUE), rowMeans(is.na(m)),
     xlab = "Mean log2 abundance", ylab = "Fraction missing", pch = 16, cex = .4)

# Label-free filtering: <= 30% missing, then top-variance features
prep <- function(m, max_missing = 0.3, n_top = 200) {
  m <- m[rowMeans(is.na(m)) <= max_missing, ]
  v <- apply(m, 1, var, na.rm = TRUE)
  m[order(-v)[seq_len(min(n_top, nrow(m)))], ]
}
model_views <- lapply(views, prep)

# Sample PCA needs complete data: impute with the feature median (for PCA only)
imp <- function(m) { m[is.na(m)] <- apply(m, 1, median, na.rm = TRUE)[row(m)[is.na(m)]]; m }
pc <- prcomp(t(imp(model_views[[2]])), scale. = TRUE)
summary(pc)$importance[2, 1:3]`,
  layers: `# Moderated t-tests (limma). Paired design: block on the pairing variable.
group <- factor(samples[[if (ds == "human_luad") "tissue" else "group"]])
pair  <- factor(samples[[if (ds == "human_luad") "patient" else "pair"]])
design_paired   <- model.matrix(~ pair + group)
design_unpaired <- model.matrix(~ group)

fit <- eBayes(lmFit(views[[2]], design_paired))
tt  <- topTable(fit, coef = ncol(design_paired), number = Inf)
head(tt)
sum(tt$adj.P.Val < 0.05)

# Compare with the unpaired test: how many fewer discoveries?
tt_u <- topTable(eBayes(lmFit(views[[2]], design_unpaired)), coef = 2, number = Inf)
sum(tt_u$adj.P.Val < 0.05)`,
  concordance: `# Per-gene Spearman correlation between RNA and protein across samples
shared <- intersect(rownames(views$RNA), rownames(views$Protein))
rho <- sapply(shared, function(g) cor(views$RNA[g, ], views$Protein[g, ],
                                      method = "spearman", use = "pairwise.complete.obs"))
median(rho, na.rm = TRUE)
hist(rho, breaks = 40, main = "", xlab = "Spearman rho (RNA vs protein)")

# Which gene sets are coupled or decoupled? Rank genes by rho and test sets
library(fgsea); library(msigdbr)
h <- msigdbr(species = "Homo sapiens", collection = "H")
sets <- split(h$gene_symbol, h$gs_name)
fgsea(sets, stats = sort(rho, decreasing = TRUE), minSize = 5)[order(padj)][1:10]`,
  integration: `# MOFA2: multi-view factor analysis
mofa <- create_mofa(lapply(model_views, function(m) t(scale(t(m)))))   # features centred + unit variance
mo   <- get_default_model_options(mofa); mo$num_factors <- 6
to   <- get_default_training_options(mofa); to$seed <- 557
mofa <- prepare_mofa(mofa, model_options = mo, training_options = to)
mofa <- run_mofa(mofa, outfile = tempfile(fileext = ".hdf5"), use_basilisk = TRUE)

samples_metadata(mofa) <- as.data.frame(samples)[, setdiff(names(samples), "")] |> transform(sample = id)
plot_variance_explained(mofa)
plot_factor(mofa, factors = 1:2, color_by = names(samples)[3])
correlate_factors_with_covariates(mofa, covariates = c("sex", "age"), plot = "log_pval")
plot_top_weights(mofa, view = 1, factor = 1, nfeatures = 10)`,
  network: `# Cross-view co-variation among the top-weighted features of one factor
W <- get_weights(mofa, factors = 1)
top <- unlist(lapply(names(W), function(v) {
  w <- W[[v]][, 1]; paste0(v, "::", names(sort(abs(w), decreasing = TRUE))[1:8]) }))
X <- do.call(rbind, lapply(names(model_views), function(v) {
  m <- model_views[[v]]; rownames(m) <- paste0(v, "::", rownames(m)); m }))
C <- cor(t(X[intersect(top, rownames(X)), ]), method = "spearman", use = "pairwise.complete.obs")
library(igraph)
g <- graph_from_adjacency_matrix(abs(C) >= 0.6 & upper.tri(C), mode = "undirected")
plot(g, vertex.size = 6, vertex.label.cex = .6)`,
  pathways: `# Rank-based enrichment of factor weights (or moderated t)
w <- get_weights(mofa, views = 2, factors = 1)[[1]][, 1]
names(w) <- sub("_Protein$", "", names(w))
library(fgsea)
fgsea(sets, stats = sort(w), minSize = 5, maxSize = 400)[order(padj)][1:15]

# Over-representation of the top 40 positive weights against the modelled background
top <- names(sort(w, decreasing = TRUE))[1:40]
res <- clusterProfiler::enricher(top, universe = names(w),
         TERM2GENE = data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets)))`,
  predict: `# Honest cross-validation: feature selection INSIDE each training fold
library(glmnet)
cv_auc <- function(X, y, k = 5, top = 20, leaky = FALSE) {
  folds <- sample(rep(1:k, length.out = length(y)))
  pick <- function(idx) { t <- apply(X[, idx], 1, function(r) abs(t.test(r ~ y[idx])$statistic)); names(sort(t, decreasing = TRUE))[1:top] }
  if (leaky) sel_all <- pick(seq_along(y))            # <- the mistake
  score <- numeric(length(y))
  for (f in 1:k) {
    tr <- folds != f; sel <- if (leaky) sel_all else pick(which(tr))
    fit <- glmnet(t(imp(X[sel, tr])), y[tr], family = "binomial", alpha = 0, lambda = 0.1)
    score[!tr] <- predict(fit, t(imp(X[sel, !tr, drop = FALSE])))
  }
  pROC::auc(y, score)
}`,
  challenge: `# Break the pairing on purpose: shuffle the sample labels of one view
broken <- model_views
perm <- sample(ncol(broken[[3]]))
colnames(broken[[3]]) <- colnames(broken[[3]])[perm]   # values move, IDs no longer match
broken <- lapply(broken, function(m) m[, samples$id])  # re-align by ID

# Agreement between views: correlation of sample-distance matrices (Mantel statistic)
d <- lapply(broken, function(m) dist(t(imp(m))))
cor(c(d[[1]]), c(d[[3]]), method = "spearman")`,
  export: `# Everything the browser computed is in the results bundle (CSV + JSON + SVG figures).
# The complete R walkthrough lives in r/walkthrough.R and reproduces every chapter with MOFA2 and limma.
source("r/walkthrough.R")`
};
