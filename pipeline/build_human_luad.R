## Build the human teaching dataset: CPTAC lung adenocarcinoma (harmonized pan-cancer release v1.2)
## Li Y, Dou Y, et al. Proteogenomic data and resources for pan-cancer analysis. Cancer Cell 41:1397 (2023)
## Gillette MA, et al. Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. Cell 182:200 (2020)
suppressMessages({ library(data.table); library(jsonlite); library(org.Hs.eg.db); library(msigdbr) })
raw <- "raw/cptac_luad"; out <- "../site/data"; dir.create(out, showWarnings = FALSE)
POOL <- 1500; MAXMISS <- 0.5

rd <- function(f) { x <- fread(file.path(raw, f), na.strings = c("NA", "")); m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m }
pair_views <- function(prefix) {
  tu <- rd(sprintf("%s_Tumor.txt", prefix)); no <- rd(sprintf("%s_Normal.txt", prefix))
  colnames(tu) <- paste0(colnames(tu), ".T"); colnames(no) <- paste0(colnames(no), ".N")
  ids <- intersect(rownames(tu), rownames(no)); cbind(tu[ids, ], no[ids, ])
}
rna  <- pair_views("LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2")
prot <- pair_views("LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized")
phos <- pair_views("LUAD_phospho_site_abundance_log2_reference_intensity_normalized")
cat("raw dims RNA", dim(rna), " protein", dim(prot), " phospho", dim(phos), "\n")

## ---- samples: patients with tumour AND normal in all three layers
meta <- fread(file.path(raw, "LUAD_meta.txt"))[-1]
pts <- Reduce(intersect, lapply(list(rna, prot, phos), function(m) sub("\\.[TN]$", "", colnames(m))))
both <- pts[sapply(pts, function(p) all(sapply(list(rna, prot, phos), function(m) all(paste0(p, c(".T", ".N")) %in% colnames(m)))))]
both <- sort(intersect(both, meta$case_id))
sids <- c(paste0(both, ".T"), paste0(both, ".N"))
cat("patients with T+N in all layers:", length(both), "\n")
rna <- rna[, sids]; prot <- prot[, sids]; phos <- phos[, sids]

## ---- feature annotation
ens <- sub("\\..*$", "", rownames(rna)); sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
gtype <- mapIds(org.Hs.eg.db, ens, "GENETYPE", "ENSEMBL", multiVals = "first"); sym[is.na(gtype) | gtype != "protein-coding"] <- NA
ensP <- sub("\\..*$", "", rownames(prot)); symP <- mapIds(org.Hs.eg.db, ensP, "SYMBOL", "ENSEMBL", multiVals = "first")
collapse <- function(m, s) { keep <- !is.na(s); m <- m[keep, , drop = FALSE]; s <- s[keep]
  mu <- rowMeans(m, na.rm = TRUE); o <- order(-mu); m <- m[o, ]; s <- s[o]; d <- !duplicated(s); m <- m[d, ]; rownames(m) <- s[d]; m }
rnaS <- collapse(rna, sym); protS <- collapse(prot, symP)
cat("mapped symbols RNA", nrow(rnaS), " protein", nrow(protS), "\n")
## phosphosites: ENSG|ENSP|site|peptide|n -> SYMBOL site
pp <- tstrsplit(rownames(phos), "|", fixed = TRUE)
symPh <- mapIds(org.Hs.eg.db, sub("\\..*$", "", pp[[1]]), "SYMBOL", "ENSEMBL", multiVals = "first")
phLab <- ifelse(is.na(symPh), NA, paste0(symPh, " ", pp[[3]]))
keep <- !is.na(phLab) & !duplicated(phLab); phos <- phos[keep, ]; rownames(phos) <- phLab[keep]; phGene <- symPh[keep]; names(phGene) <- rownames(phos)
cat("phosphosites with symbol", nrow(phos), "\n")

## ---- label-free teaching pools: expressed / <=50% missing, then top-variance
missf <- function(m) rowMeans(is.na(m))
rnaE <- rnaS[apply(rnaS, 1, median) >= 4 & rowMeans(rnaS > 0) >= 0.9, ]
pool <- function(m, n) { m <- m[missf(m) <= MAXMISS, ]; v <- apply(m, 1, var, na.rm = TRUE); m[order(-v)[seq_len(min(n, nrow(m)))], ] }
rnaP <- pool(rnaE, POOL); protP <- pool(protS, POOL); phosP <- pool(phos, POOL)

## ---- full-matrix RNA-protein concordance (offline, all shared genes)
shared <- intersect(rownames(rnaS), rownames(protS))
conc <- rbindlist(lapply(shared, function(g) {
  r <- rnaS[g, ]; p <- protS[g, ]; ok <- !is.na(p) & !is.na(r)
  okT <- ok & grepl("\\.T$", sids); okN <- ok & grepl("\\.N$", sids)
  data.table(gene = g, n = sum(ok),
    rho_all = if (sum(ok) >= 20) cor(r[ok], p[ok], method = "spearman") else NA_real_,
    rho_tumor = if (sum(okT) >= 20) cor(r[okT], p[okT], method = "spearman") else NA_real_,
    rho_normal = if (sum(okN) >= 20) cor(r[okN], p[okN], method = "spearman") else NA_real_,
    protein_mean = mean(p, na.rm = TRUE), rna_mean = mean(r, na.rm = TRUE), protein_missing = mean(is.na(p)))
}))
conc <- conc[!is.na(rho_tumor)]
cat("RNA-protein concordance genes:", nrow(conc), " median rho (tumour):", round(median(conc$rho_tumor), 3), "\n")

## ---- metadata
mt <- meta[match(sub("\\.[TN]$", "", sids), case_id)]
smoke <- fifelse(mt$Tobacco_smoking_history %in% c("current smoker", "past smoker"), "Ever", fifelse(mt$Tobacco_smoking_history == "non-smoker", "Never", NA_character_))
samples <- data.table(id = sids, patient = mt$case_id, tissue = fifelse(grepl("\\.T$", sids), "Tumor", "Normal"),
  sex = mt$Sex, age = suppressWarnings(as.numeric(mt$Age)), smoking = smoke, stage = mt$Stage,
  EGFR = fifelse(mt$EGFR_mutation == "1", "Mutant", "Wild type"), KRAS = fifelse(mt$KRAS_mutation == "1", "Mutant", "Wild type"),
  TP53 = fifelse(mt$TP53_mutation == "1", "Mutant", "Wild type"), STK11 = fifelse(mt$STK11_mutation == "1", "Mutant", "Wild type"))
# tumour-only genotype labels: normals carry the same germline, so mark them "Normal tissue" to avoid misreading
for (g in c("EGFR", "KRAS", "TP53", "STK11", "stage")) samples[tissue == "Normal", (g) := "Normal tissue"]
print(samples[, .N, by = .(tissue, sex)])

## ---- gene sets restricted to the universe of genes present in any pool
univ <- unique(c(rownames(rnaP), rownames(protP), unname(phGene[rownames(phosP)]), conc$gene))
gs <- function(coll, sub = NULL, prefix) {
  x <- if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll) else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  x <- as.data.table(x)[gene_symbol %in% univ]
  sets <- split(x$gene_symbol, x$gs_name); sets <- lapply(sets, unique); sets <- sets[lengths(sets) >= 5 & lengths(sets) <= 400]
  nm <- names(sets); pretty <- gsub("_", " ", sub(paste0("^", prefix, "_"), "", nm))
  pretty <- paste0(substr(pretty, 1, 1), tolower(substring(pretty, 2)))
  list(names = nm, pretty = pretty, sets = sets)
}
catalogs <- list(hallmark = gs("H", NULL, "HALLMARK"), gobp = gs("C5", "GO:BP", "GOBP"), reactome = gs("C2", "CP:REACTOME", "REACTOME"))
gidx <- setNames(seq_along(univ) - 1L, univ)
catJSON <- lapply(catalogs, function(c) list(ids = c$names, names = c$pretty, members = unname(lapply(c$sets, function(s) unname(gidx[s])))))
cat("gene sets:", sapply(catalogs, function(c) length(c$sets)), "\n")

## ---- write JSON
r3 <- function(m) { m <- round(m, 3); lapply(seq_len(ncol(m)), function(j) { v <- m[, j]; v[is.na(v)] <- NA; v }) }
view <- function(key, name, short, assay, unit, m, symbols, labels, note) list(
  key = key, name = name, short = short, assay = assay, unit = unit, note = note,
  sample_ids = colnames(m), features = rownames(m), symbols = unname(symbols), labels = unname(labels), values = r3(m))
ds <- list(
  id = "human_luad", organism = "Homo sapiens", species_short = "Human",
  title = "Lung adenocarcinoma tumours and matched normal lung",
  subtitle = sprintf("%d patients · tumour and adjacent normal tissue · RNA, proteins and phosphosites", length(both)),
  contrast = list(key = "tissue", levels = c("Normal", "Tumor"), labels = c("Normal lung", "Tumour")),
  pairing = list(key = "patient", label = "Patient", note = "Each patient contributes one tumour and one normal-adjacent sample."),
  covariates = list(
    list(key = "tissue", label = "Tissue", type = "categorical"), list(key = "patient", label = "Patient", type = "id"),
    list(key = "sex", label = "Sex", type = "categorical"), list(key = "age", label = "Age", type = "numeric"),
    list(key = "smoking", label = "Smoking", type = "categorical"), list(key = "stage", label = "Stage", type = "categorical", tumour_only = TRUE),
    list(key = "EGFR", label = "EGFR", type = "categorical", tumour_only = TRUE), list(key = "KRAS", label = "KRAS", type = "categorical", tumour_only = TRUE),
    list(key = "TP53", label = "TP53", type = "categorical", tumour_only = TRUE), list(key = "STK11", label = "STK11", type = "categorical", tumour_only = TRUE)),
  samples = samples,
  views = list(
    view("rna", "Transcriptome", "RNA", "RNA-seq, RSEM gene level", "log₂ RSEM (upper-quartile)", rnaP, rownames(rnaP), rownames(rnaP), "Protein-coding genes; expressed (median ≥ 4, detected in ≥ 90% of samples)."),
    view("protein", "Proteome", "Protein", "TMT-11 LC–MS/MS", "log₂ reference-normalised intensity", protP, rownames(protP), rownames(protP), "Gene-level protein abundance; missing values are not imputed."),
    view("phospho", "Phosphoproteome", "Phospho", "IMAC-enriched TMT-11 LC–MS/MS", "log₂ reference-normalised intensity", phosP, unname(phGene[rownames(phosP)]), rownames(phosP), "Single phosphosites, not normalised to protein abundance.")),
  concordance = list(gene = conc$gene, n = conc$n, rho_all = round(conc$rho_all, 3), rho_tumor = round(conc$rho_tumor, 3), rho_normal = round(conc$rho_normal, 3),
    protein_mean = round(conc$protein_mean, 2), rna_mean = round(conc$rna_mean, 2), protein_missing = round(conc$protein_missing, 3)),
  source = list(
    papers = list(
      list(cite = "Li Y, Dou Y, et al. Proteogenomic data and resources for pan-cancer analysis. Cancer Cell 41, 1397–1406 (2023).", doi = "10.1016/j.ccell.2023.06.009"),
      list(cite = "Gillette MA, Satpathy S, Cao S, et al. Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. Cell 182, 200–225 (2020).", doi = "10.1016/j.cell.2020.06.013")),
    files = list(
      list(name = "CPTAC pan-cancer harmonized data freeze v1.2 (LUAD)", url = "https://kb.linkedomics.org/download"),
      list(name = "Public S3 bucket", url = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_meta.txt")),
    processing = c(
      sprintf("Kept the %d patients with a tumour and a normal-adjacent sample in all three layers (%d samples).", length(both), length(sids)),
      "Mapped Ensembl gene IDs to HGNC symbols with org.Hs.eg.db; duplicated symbols keep the most abundant row. Phosphosites are labelled 'SYMBOL site'.",
      sprintf("Teaching pools are label-free: features with ≤ %d%% missing values (RNA: expressed genes), ranked by variance across all samples; the top %d per layer are shipped.", MAXMISS * 100, POOL),
      sprintf("RNA–protein concordance was computed offline on all %d genes measured in both layers (Spearman, ≥ 20 paired observations).", nrow(conc)),
      "Values are rounded to three decimals. No imputation is applied before shipping."),
    limitations = c(
      "Normal-adjacent tissue is not healthy lung: it can carry field effects of the tumour and of smoking.",
      "Tumour purity and immune/stromal content vary between patients and drive part of the variation.",
      "Phosphosite abundance partly reflects protein abundance; the shipped values are not protein-normalised.",
      "Genotype labels (EGFR, KRAS, TP53, STK11) and stage describe the tumour and are shown as 'Normal tissue' for normal samples."),
    license = "CPTAC data are public; cite the CPTAC publications above when reusing them."),
  genesets_file = "human_luad_genesets.json")
writeLines(toJSON(ds, auto_unbox = TRUE, na = "null", digits = NA, dataframe = "rows"), file.path(out, "human_luad.json"))
writeLines(toJSON(list(genes = univ, catalogs = catJSON, labels = list(hallmark = "MSigDB Hallmark", gobp = "GO biological process", reactome = "Reactome"), version = "MSigDB via msigdbr 26.1.0"), auto_unbox = TRUE), file.path(out, "human_luad_genesets.json"))
## also keep tidy matrices for the R walkthrough
dir.create("../r/data", showWarnings = FALSE, recursive = TRUE)
fwrite(as.data.table(rnaP, keep.rownames = "feature"), "../r/data/human_luad_rna.csv")
fwrite(as.data.table(protP, keep.rownames = "feature"), "../r/data/human_luad_protein.csv")
fwrite(as.data.table(phosP, keep.rownames = "feature"), "../r/data/human_luad_phospho.csv")
fwrite(samples, "../r/data/human_luad_samples.csv")
cat("done:", file.size(file.path(out, "human_luad.json")) / 1e6, "MB\n")
