## Build the honey bee teaching dataset from the public DIA proteomics deposit of
## Zhong H, et al. Proteomics (2025) doi:10.1002/pmic.70033 — MassIVE MSV000097590 / PXD062819 (CC BY 4.0 article)
## Highbush blueberry (HBB) apiaries, 2020 (HBB01–05) and 2021 (HBB06–10); colony groups e/u; timepoints t2/t4;
## head (H), abdomen (A) and midgut (G) dissected from the same pool of 12 nurse bees.
suppressMessages({ library(data.table); library(jsonlite); library(GO.db) })
raw <- "raw/zhong2025"; out <- "../site/data"; POOL <- 1500; MAXMISS <- 0.5

## ---- protein annotation from the search database
fa <- readLines(file.path(raw, "bee_db.fasta")); hdr <- sub("^>", "", fa[startsWith(fa, ">")])
acc <- ifelse(grepl("^(sp|tr)\\|", hdr), sub("^(sp|tr)\\|([^|]+)\\|.*", "\\2", hdr), sub("\\.[0-9]+$", "", sub(" .*", "", hdr)))
org <- ifelse(grepl("OS=", hdr), sub(".*OS=(.*?) OX=.*", "\\1", hdr), sub(".*\\[(.*)\\]$", "\\1", hdr))
gn <- ifelse(grepl(" GN=", hdr), sub(".* GN=([^ ]+).*", "\\1", hdr), NA)
desc <- ifelse(grepl("^(sp|tr)\\|", hdr), sub("^(sp|tr)\\|[^ ]+ (.*?) OS=.*", "\\2", hdr), sub("^[^ ]+ (.*) \\[.*$", "\\1", hdr))
ann <- data.table(acc, org, gn, desc)[!duplicated(acc)]
ann[, kind := fifelse(org == "Apis mellifera", "bee", fifelse(grepl("Nosema|Vairimorpha", org), "nosema", "virus"))]

## ---- read the two highbush-blueberry searches (each year was searched and normalised separately)
rd <- function(f, year) {
  x <- fread(file.path(raw, f)); s <- grep("BeeCSI", names(x), value = TRUE)
  code <- regmatches(s, regexpr("HBB[0-9]+[eu]_t[0-9]_[AGH]", s))
  stopifnot(length(code) == length(s))
  m <- as.matrix(x[, ..s]); colnames(m) <- code; rownames(m) <- x$Protein.Group
  plate <- regmatches(s, regexpr("_P[0-9]+_", s)); list(m = m, plate = setNames(gsub("_", "", plate), code), year = year)
}
a <- rd("hbb2020_pg.tsv", 2020); b <- rd("hbb2021_pg.tsv", 2021)
clean <- function(ids) vapply(strsplit(ids, ";"), function(v) paste(sub(" .*", "", trimws(v)), collapse = ";"), "")
rownames(a$m) <- clean(rownames(a$m)); rownames(b$m) <- clean(rownames(b$m))
pg <- union(rownames(a$m), rownames(b$m))
M <- cbind(a$m[match(pg, rownames(a$m)), ], b$m[match(pg, rownames(b$m)), ]); rownames(M) <- pg
M[M <= 0] <- NA; M <- log2(M)
plate <- c(a$plate, b$plate)
cat("protein groups", nrow(M), " runs", ncol(M), "\n")

## first accession of each group decides species & annotation
first <- sub("\\.[0-9]+$", "", sub(";.*", "", pg)); A <- ann[match(first, acc)]
kind <- fifelse(is.na(A$kind), "unknown", A$kind)
cat("groups by origin:"); print(table(kind))

## ---- sample design
cd <- colnames(M); pool <- sub("_[AGH]$", "", cd); tissue <- c(A = "Abdomen", G = "Gut", H = "Head")[sub(".*_", "", cd)]
pools <- sort(unique(pool)); stopifnot(all(table(pool) == 3))
site <- sub("[eu]_t[0-9]$", "", pools); grp <- sub("^HBB[0-9]+([eu])_.*", "\\1", pools); tp <- sub(".*_", "", pools)
yr <- ifelse(as.integer(sub("HBB", "", site)) <= 5, 2020L, 2021L)

## pathogen protein signal per pool (summed linear intensity of Nosema / virus protein groups, log2(1 + x), gut and whole pool)
lin <- 2^M; lin[is.na(lin)] <- 0
path_sig <- function(k, tis) sapply(pools, function(p) { cols <- which(pool == p & (is.null(tis) | tissue %in% tis)); log2(1 + sum(lin[kind == k, cols])) })
samples <- data.table(id = pools, site = site, group = paste("Group", grp), timepoint = paste("Timepoint", tp), year = as.character(yr), pair = paste(site, grp, sep = "_"),
  plate = unname(plate[paste0(pools, "_H")]), nosema = round(path_sig("nosema", "Gut"), 2), virus = round(path_sig("virus", NULL), 2))
print(samples[, .N, by = .(group, timepoint, year)])

## ---- one view per tissue: Apis mellifera protein groups, label-free pools
bee <- kind == "bee"
xr_all <- fread(file.path(raw, "uniprot_xref.tsv"), quote = "")
gn_up <- xr_all[, .(gn = { g <- unlist(strsplit(gene_names, " ")); g <- g[!grepl("^LOC|^GB[0-9]", g) & nzchar(g)]; if (length(g)) g[1] else NA_character_ }), by = refseq][!is.na(gn)][!duplicated(refseq)]
label_of <- function(ids) { f <- sub("\\.[0-9]+$", "", sub(";.*", "", ids)); r <- ann[match(f, acc)]
  up <- gn_up$gn[match(f, gn_up$refseq)]
  d <- sub("^LOW QUALITY PROTEIN: ", "", r$desc); d <- sub(" isoform.*$", "", d); d <- sub(", partial$", "", d)
  d <- sub("^uncharacterized protein (LOC[0-9]+)$", "\\1", d, ignore.case = TRUE)
  d <- sub("-like$", "-like", d); d <- ifelse(nchar(d) > 26, paste0(substr(d, 1, 25), "…"), d)
  lab <- ifelse(!is.na(up), up, d); lab[is.na(lab)] <- f[is.na(lab)]
  make.unique(lab, sep = " #") }
views <- list(); used <- character()
for (t in c("Head", "Abdomen", "Gut")) {
  m <- M[bee, tissue == t]; colnames(m) <- pool[tissue == t]; m <- m[, pools]
  m <- m[rowMeans(is.na(m)) <= MAXMISS, ]; v <- apply(m, 1, var, na.rm = TRUE); m <- m[order(-v)[seq_len(min(POOL, nrow(m)))], ]
  views[[t]] <- m; used <- union(used, sub(";.*", "", rownames(m)))
  cat(t, "proteins in pool", nrow(m), " missing", round(mean(is.na(m)), 3), "\n")
}

## ---- GO biological process from UniProt-GOA (A. mellifera), propagated to ancestors with GO.db
xr <- fread(file.path(raw, "uniprot_xref.tsv"), quote = "")[refseq %in% used & go_p != ""]
goa <- xr[, .(go = unlist(regmatches(go_p, gregexpr("GO:[0-9]{7}", go_p)))), by = .(acc = refseq)]
## better labels: UniProt protein name where the RefSeq description is generic
anc <- as.list(GOBPANCESTOR); prop <- goa[, .(go = unique(c(go, unlist(anc[go])))), by = acc][go != "all"]
sets <- split(prop$acc, prop$go); sets <- lapply(sets, unique); sets <- sets[lengths(sets) >= 5 & lengths(sets) <= 300]
nm <- Term(names(sets)); keep <- !is.na(nm); sets <- sets[keep]; nm <- nm[keep]
## drop near-duplicate parent/child sets with identical membership
sig <- vapply(sets, function(s) paste(sort(s), collapse = ","), ""); d <- !duplicated(sig); sets <- sets[d]; nm <- nm[d]
univ <- sort(used); gidx <- setNames(seq_along(univ) - 1L, univ)
cat("GO BP sets", length(sets), " proteins with any BP annotation (UniProt via RefSeq cross-reference)", length(unique(goa$acc)), "of", length(used), "\n")

r3 <- function(m) lapply(seq_len(ncol(m)), function(j) { v <- round(m[, j], 3); v[is.na(v)] <- NA; v })
mkview <- function(key, t, m) list(key = key, name = paste(t, "proteome"), short = t, assay = "DIA LC–MS/MS (DIA-NN, MaxLFQ)", unit = "log₂ MaxLFQ intensity",
  note = sprintf("Apis mellifera protein groups quantified in ≤ %d%% missing of %s samples; top %d by variance.", MAXMISS * 100, tolower(t), POOL),
  sample_ids = colnames(m), features = rownames(m), symbols = sub(";.*", "", rownames(m)), labels = label_of(rownames(m)), values = r3(m))
ds <- list(
  id = "honeybee_hbb", organism = "Apis mellifera", species_short = "Honey bee",
  title = "Nurse bees from highbush blueberry apiaries: three tissues from one pool",
  subtitle = sprintf("%d pools of 12 nurse bees · head, abdomen and midgut proteomes · 10 apiaries, 2020–2021", length(pools)),
  contrast = list(key = "timepoint", levels = c("Timepoint t2", "Timepoint t4"), labels = c("Timepoint t2", "Timepoint t4")),
  pairing = list(key = "pair", label = "Colony set", note = "Each apiary × colony group was sampled at both timepoints."),
  batch = "year",
  covariates = list(
    list(key = "timepoint", label = "Timepoint", type = "categorical"), list(key = "pair", label = "Colony set", type = "id"),
    list(key = "group", label = "Colony group (e/u)", type = "categorical"), list(key = "year", label = "Year / MS batch", type = "categorical"),
    list(key = "site", label = "Apiary", type = "categorical"), list(key = "nosema", label = "Nosema protein signal (gut)", type = "numeric"),
    list(key = "virus", label = "Virus protein signal", type = "numeric")),
  vp_covariates = c("group", "year"),
  samples = samples,
  views = list(mkview("head", "Head", views$Head), mkview("abdomen", "Abdomen", views$Abdomen), mkview("gut", "Gut", views$Gut)),
  source = list(
    papers = list(list(cite = "Zhong H, Shi Y, Rogalski JC, et al., Zayed A, Foster LJ. Omics insights into the effects of highbush blueberry and cranberry crop agroecosystems on honey bee health and physiology. Proteomics (2025).", doi = "10.1002/pmic.70033")),
    files = list(list(name = "MassIVE MSV000097590 (ProteomeXchange PXD062819): search/SEARCH/Results_HBB_2020 and Results_HBB_2021, DIA-NN pg_matrix", url = "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000097590"),
      list(name = "UniProtKB GO biological process terms, matched through RefSeq cross-references", url = "https://rest.uniprot.org/uniprotkb/search?query=organism_id:7460")),
    processing = c(
      "Run names were parsed as apiary (HBB01–10), colony group (e/u), timepoint (t2/t4) and tissue (H/A/G); every pool has all three tissues.",
      "The 2020 and 2021 matrices come from separate DIA-NN searches; they were joined by protein group and log₂-transformed without further normalisation.",
      "Protein groups were assigned to honey bee, Nosema or virus using the search FASTA. Pathogen protein signal = log₂(1 + summed intensity) per pool.",
      sprintf("Each tissue view keeps honey bee protein groups with ≤ %d%% missing values in that tissue, ranked by variance without labels; top %d shipped.", MAXMISS * 100, POOL),
      "GO biological process terms were taken from UniProtKB entries cross-referenced to each RefSeq protein (mostly electronic annotation), propagated to ancestor terms with GO.db; sets of 5–300 proteins."),
    limitations = c(
      "The deposit contains more runs than the 60 analysed in the paper; colony-group codes e/u and timepoints t2/t4 are used as recorded in the file names.",
      "Year is fully confounded with apiary set, MS plate and DIA-NN search: treat 'Year / MS batch' as a batch variable.",
      "Each sample is a pool of 12 bees, so between-bee variation is averaged out and n = 40 pools.",
      "This module currently has one molecular layer (protein) in three tissues; transcript counts from the same pools can be added as further views."),
    license = "Article CC BY 4.0; MassIVE/ProteomeXchange public dataset. Cite Zhong et al. (2025) when reusing."),
  genesets_file = "honeybee_hbb_genesets.json")
writeLines(toJSON(ds, auto_unbox = TRUE, na = "null", digits = NA, dataframe = "rows"), file.path(out, "honeybee_hbb.json"))
writeLines(toJSON(list(genes = univ, catalogs = list(gobp = list(ids = names(sets), names = unname(nm), members = unname(lapply(sets, function(s) unname(gidx[s]))))),
  labels = list(gobp = "GO biological process (UniProtKB)"), version = "UniProtKB (Nov 2026) via RefSeq cross-references; GO.db ancestors"), auto_unbox = TRUE), file.path(out, "honeybee_hbb_genesets.json"))
dir.create("../r/data", showWarnings = FALSE, recursive = TRUE)
for (t in names(views)) fwrite(as.data.table(views[[t]], keep.rownames = "protein_group"), sprintf("../r/data/honeybee_%s.csv", tolower(t)))
fwrite(samples, "../r/data/honeybee_samples.csv")
cat("done:", round(file.size(file.path(out, "honeybee_hbb.json")) / 1e6, 2), "MB\n")
