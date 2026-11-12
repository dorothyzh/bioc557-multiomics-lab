## BIOC 557 · Multi-omics Lab — companion R walkthrough
## Reproduces the browser lab with limma, MOFA2 and fgsea.
## Run from the folder that contains r/ (setwd to the site folder), or edit DATA_DIR.
## install.packages(c("data.table", "ggplot2", "patchwork"))
## BiocManager::install(c("limma", "MOFA2", "fgsea", "msigdbr"))

suppressPackageStartupMessages({ library(data.table); library(limma); library(MOFA2); library(ggplot2) })
DATA_DIR <- "r/data"
DATASET  <- "honeybee"          # "honeybee" or "human_luad"
set.seed(557)

## ---- 1. Load and align by sample ID ------------------------------------------------------------
read_view <- function(file) { x <- fread(file.path(DATA_DIR, file)); m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m }
if (DATASET == "human_luad") {
  samples  <- fread(file.path(DATA_DIR, "human_luad_samples.csv"))
  views    <- list(RNA = read_view("human_luad_rna.csv"), Protein = read_view("human_luad_protein.csv"), Phospho = read_view("human_luad_phospho.csv"))
  group    <- factor(samples$tissue, levels = c("Normal", "Tumor")); pair <- factor(samples$patient)
} else {
  samples  <- fread(file.path(DATA_DIR, "honeybee_samples.csv"))
  views    <- list(Head = read_view("honeybee_head.csv"), Abdomen = read_view("honeybee_abdomen.csv"), Gut = read_view("honeybee_gut.csv"))
  group    <- factor(samples$timepoint, levels = c("Timepoint t2", "Timepoint t4")); pair <- factor(samples$pair)
}
views <- lapply(views, function(m) m[, samples$id])     # never align by position
stopifnot(all(sapply(views, function(m) identical(colnames(m), samples$id))))

## Optional: centre each feature within batch (honey bee: year = MS batch)
BATCH_ADJUST <- FALSE
if (BATCH_ADJUST && DATASET == "honeybee") views <- lapply(views, function(m) {
  for (b in unique(samples$year)) { j <- samples$year == b; m[, j] <- m[, j] - rowMeans(m[, j], na.rm = TRUE) + rowMeans(m, na.rm = TRUE) }; m })

## ---- 2. Quality control ---------------------------------------------------------------------------
sapply(views, function(m) mean(is.na(m)))
mnar <- rbindlist(lapply(names(views), function(v) data.table(layer = v, mean = rowMeans(views[[v]], na.rm = TRUE), missing = rowMeans(is.na(views[[v]])))))
ggplot(mnar, aes(mean, missing)) + geom_point(size = .3, alpha = .3) + geom_smooth(se = FALSE, colour = "black", linewidth = .5) +
  facet_wrap(~layer, scales = "free_x") + labs(x = "Mean log2 abundance", y = "Fraction missing") + theme_classic(base_size = 7)

prep <- function(m, max_missing = 0.3, n_top = 200) {
  m <- m[rowMeans(is.na(m)) <= max_missing, ]
  v <- apply(m, 1, var, na.rm = TRUE); m[order(-v)[seq_len(min(n_top, nrow(m)))], ]
}
model_views <- lapply(views, prep)
impute_median <- function(m) { med <- apply(m, 1, median, na.rm = TRUE); idx <- which(is.na(m), arr.ind = TRUE); m[idx] <- med[idx[, 1]]; m }

## ---- 3. Each layer: PCA and moderated tests (paired) ------------------------------------------------
pcs <- lapply(model_views, function(m) prcomp(t(impute_median(m)), scale. = TRUE))
sapply(pcs, function(p) round(summary(p)$importance[2, 1:3], 3))

design  <- model.matrix(~ pair + group)
de <- lapply(views, function(m) { m <- m[rowMeans(is.na(m)) <= 0.3, ]
  topTable(eBayes(lmFit(m, design)), coef = ncol(design), number = Inf) })
sapply(de, function(t) sum(t$adj.P.Val < 0.05))
de_unpaired <- lapply(views, function(m) { m <- m[rowMeans(is.na(m)) <= 0.3, ]
  topTable(eBayes(lmFit(m, model.matrix(~ group))), coef = 2, number = Inf) })
sapply(de_unpaired, function(t) sum(t$adj.P.Val < 0.05))   # compare with the paired analysis

## ---- 4. Across layers (feature-level concordance) ----------------------------------------------------
a <- views[[1]]; b <- views[[2]]; shared <- intersect(rownames(a), rownames(b))
rho <- sapply(shared, function(g) suppressWarnings(cor(a[g, ], b[g, ], method = "spearman", use = "pairwise.complete.obs")))
summary(rho)

## ---- 5. Integration with MOFA2 -------------------------------------------------------------------------
mofa <- create_mofa(lapply(model_views, function(m) t(scale(t(m)))))
mo <- get_default_model_options(mofa); mo$num_factors <- ifelse(DATASET == "honeybee", 4, 6)
to <- get_default_training_options(mofa); to$seed <- 557; to$verbose <- FALSE
mofa <- prepare_mofa(mofa, model_options = mo, training_options = to)
mofa <- run_mofa(mofa, outfile = tempfile(fileext = ".hdf5"), use_basilisk = TRUE)
md <- as.data.frame(samples); md$sample <- md$id
if ("group" %in% names(md)) { md$colony_group <- md$group }; md$group <- "group1"   # MOFA2 reserves the column "group"
samples_metadata(mofa) <- md
plot_variance_explained(mofa)
covs <- intersect(c("tissue", "sex", "age", "smoking", "timepoint", "colony_group", "year", "nosema", "virus"), names(md))
correlate_factors_with_covariates(mofa, covariates = covs, plot = "log_pval")
plot_factor(mofa, factors = 1:2, color_by = covs[1])

## ---- 6. Functional interpretation (rank-based) -----------------------------------------------------------
if (DATASET == "human_luad") {
  library(fgsea); library(msigdbr)
  h <- msigdbr(species = "Homo sapiens", collection = "H"); sets <- split(h$gene_symbol, h$gs_name)
  w <- get_weights(mofa, views = "Protein", factors = 1)[[1]][, 1]; names(w) <- sub("_Protein$", "", names(w))
  print(head(fgsea(sets, stats = sort(w), minSize = 5)[order(padj)], 10))
  t_stat <- setNames(de$Protein$t, rownames(de$Protein))
  print(head(fgsea(sets, stats = sort(t_stat), minSize = 5)[order(padj)], 10))
}

## ---- 7. Honest cross-validation (feature selection inside folds) ------------------------------------------
cv_auc <- function(m, y, k = 5, top = 20, leaky = FALSE) {
  m <- impute_median(m[rowMeans(is.na(m)) <= 0.3, ]); folds <- sample(rep(seq_len(k), length.out = length(y)))
  pick <- function(idx) { t <- apply(m[, idx], 1, function(r) abs(t.test(r[y[idx] == 1], r[y[idx] == 0])$statistic)); names(sort(t, decreasing = TRUE))[1:top] }
  sel_all <- if (leaky) pick(seq_along(y)) else NULL; score <- numeric(length(y))
  for (f in seq_len(k)) { tr <- folds != f; sel <- if (leaky) sel_all else pick(which(tr))
    x <- t(m[sel[1:5], ]); mu <- colMeans(x[tr, ]); sdv <- apply(x[tr, ], 2, sd); x <- sweep(sweep(x, 2, mu), 2, sdv, "/")
    b <- coef(suppressWarnings(glm.fit(cbind(1, x[tr, ]), y[tr], family = binomial()))); b[is.na(b)] <- 0
    score[!tr] <- cbind(1, x[!tr, , drop = FALSE]) %*% b }   # training-fold statistics only
  r <- rank(score); (sum(r[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) / (sum(y == 1) * sum(y == 0))
}
y <- as.integer(group == levels(group)[2])
sapply(views, cv_auc, y = y)
sapply(views, cv_auc, y = sample(y), leaky = TRUE)   # shuffled labels: should be ~0.5, is it?
