# BIOC 557 · Multi-omics Lab

An interactive, browser-only multi-omics teaching lab with two public datasets:

| | Human | Honey bee |
|---|---|---|
| Study | CPTAC lung adenocarcinoma, harmonized release (Li et al., *Cancer Cell* 2023) | Zhong et al., *Proteomics* 2025 (MassIVE MSV000097590) |
| Samples | 101 patients × tumour/normal = 202 | 40 pools of 12 nurse bees (10 apiaries × 2 colony sets × 2 timepoints) |
| Layers | RNA-seq · TMT proteome · phosphoproteome | Head · abdomen · midgut DIA proteomes |
| Contrast | Tumour vs normal, paired by patient | Timepoint t4 vs t2, paired by colony set; year = MS batch |

## Folder layout

```
site/                 ← the website (static; deploy this folder)
  index.html          ten chapters, all analysis runs in a Web Worker
  sources.html        data provenance and methods
  js/engine.js        statistics: moderated t (limma-style), PCA, Bayesian group factor analysis (MOFA-type), enrichment, CV
  js/worker.js        pipeline per run; js/charts.js publication-style SVG figures; js/app.js UI
  data/*.json         teaching matrices + gene-set catalogues
  r/walkthrough.R     companion R script (limma, MOFA2, fgsea) + r/data/*.csv
pipeline/             build scripts (R) that download → process → write site/data
```

## Run locally

```bash
cd site && python3 -m http.server 8557
```

Then open http://127.0.0.1:8557 (a server is needed because the analysis runs in a Web Worker).

## Rebuild the data

Raw files go in `pipeline/raw/` (download URLs are at the top of each script and on `sources.html`). Run under a UTF-8 locale:

```bash
cd pipeline
LC_ALL=en_US.UTF-8 Rscript build_human_luad.R && node prettify_genesets.js
LC_ALL=en_US.UTF-8 Rscript build_honeybee.R
```

After rebuilding, bump `DATA_VERSION` in `site/js/app.js` so browsers fetch the new files.

## Deployment

The site is served by GitHub Pages. Every push to `main` runs `.github/workflows/pages.yml`, which publishes the `site/` folder as is (there is no build step). One-time setup: Settings → Pages → Source → **GitHub Actions**.

### Custom domain (optional)

1. Register the domain with any registrar.
2. Settings → Pages → Custom domain → enter it and save (with Actions deployment a `CNAME` file is not needed).
3. At the registrar, add DNS records:
   - apex domain (`example.org`): four `A` records → `185.199.108.153`, `185.199.109.153`, `185.199.110.153`, `185.199.111.153`
   - subdomain (`www.example.org` or `lab.example.org`): one `CNAME` record → `<username>.github.io`
4. When the certificate is issued (minutes to an hour), tick **Enforce HTTPS**.

## Notes

- Honey bee colony-set codes (`e`/`u`) are not defined in the source article; the lab uses them only as labels.
- More layers can be added per dataset (for example transcript counts from the same honey bee pools); the pipeline and the factor model handle any number of views.
