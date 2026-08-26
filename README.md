# GBM-Surprisal
This is the repository for the "Spatial Omics Decodes Human Brain Tumor Interactypes via Thermodynamic Surprisal Analysis"

## Data Preprocessing Part
Please refer to https://github.com/sdparekh/zUMIs for the raw fastq to the gene count matrix processes and zUMIs setup. Our parameters are stated in the template yaml file.

## Downstream Processing
1. We run each of the sample 
2. For the NICHES analysis, we run NICHES based on log normalized assay "Spatial" with k=9 as neighboring spots, but put the label as "NeighborhoodToCell_SCT_9"
3. The surprisal analysis is based on raw expression matrix and the scripts can be found here @https://github.com/ShouryoGhose/Surprisal-Analysis-of-Human-GBM-samples

### Important: the file numbers are not the run order

The `01.`–`06.` prefixes are naming labels. Following the actual input/output paths
inside each script, the real order is:

01 → 02 → (a NICHES pass that is not in the repo) → 05 → 03 → 06 → 04

Why:

- `05.scType.R` reads `*_gen5_05_noMT_on_tissue_NICHES.rds` and writes `gen5_rdslist_06_scType.rds`
- `03.*.Rmd` reads `gen5_rdslist_06_scType.rds`
- `06.Boundary_label_*.R` reads `07.maunally_labeled_GLSB_v2_20250218.rds`, which is the saved output of step 03
- `04.NICHES_P1_Core.R` reads `P1_Core_20250409.rds`, which is the output of `06.Boundary_label_P1_Core.R`

There are two NICHES runs. The first one (before scType) has no committed script; it only
supplies the `seurat_clusters` column that scType scores against, and step 03 throws that
clustering away anyway. The second one is the `04.NICHES_*.R` set and produces the final objects.

---

## The files

### zUMIs_Patho-DBiT_template.yaml.zip

The zUMIs config template. It is password-encrypted and there is no plain-text copy, so
nobody cloning the repo can actually read it. Same for the other five zips (those are
redundant copies of scripts that are already present unencrypted).

### 01.gen5_rds_generation_20241206.R

Turns zUMIs output into on-tissue Seurat objects. Run once per sample.

- Expects a variable `sample` to already exist holding an ExperimentID like `"R1110C"`.
  The line that would define it is commented out at lines 37–39, so it has to come from
  the caller.
- Builds a sample lookup table (SampleID, TumorType, SampleType) and matches the current sample.
- Reads the DBiT barcode file and builds `x.y.cord`. Not used later in this script.
- Reads the 10x-format image folder with `Read10X_Image(filter.matrix = TRUE)` to get the
  on-tissue barcodes.
- Reads the zUMIs `.rds` and takes `umicount$inex$all` (intron+exon UMI counts).
- Creates a Seurat object with assay `Spatial`, adds metadata (`nucleotide_name`,
  `expriment_id` — the typo is in the actual column name, `sample_id`, `TumorType`,
  `SampleType`, `MT_del`), attaches the image, and subsets to on-tissue spots.
- Writes a QC png and `*_gen5_01_original_on_tissue.rds`.
- Applies a gene filter, then repeats the whole thing with `MT_del = "Yes"` and writes
  `*_gen5_02_noMT_on_tissue.rds`.

Note: the filter regex is `^(MT|AL|H3|EEF|RPS|RPL|FP|AC|ENSM|ENSG|RNA|RN5|RN7|RNU|RNV)`,
so it removes ribosomal proteins, histone H3, elongation factors, small nuclear RNAs and
unnamed Ensembl IDs as well as mitochondrial genes. The printed "MT gene percentage"
covers all of that.

Bug: line 148 subsets the filtered object with `common_barcodes` (the *unfiltered*
intersection) and line 149 attaches `img_origin` instead of the `img_noMT` computed at 143.
`common_barcodes_filtered` and `img_noMT` are never used. Change 148 to
`cells = common_barcodes_filtered` and 149 to `img_noMT`.

Minor: lines 60–62 overwrite `message` twice, so only the SampleType line prints.

### 02.gen5_rds_normalization_20241207.Rmd

Parameterised Rmd, one render per sample, selected by `params$sample_name`.

Runs two normalisations on each of the two objects from step 01:

- LogNormalize → FindVariableFeatures → ScaleData → RunPCA → FindNeighbors(1:30) → FindClusters → RunUMAP(1:30)
- SCTransform(assay="Spatial") → RunPCA → FindNeighbors(1:30) → FindClusters(res=0.8) → RunUMAP(1:30)

Writes four `.rds` and four pngs per sample:
`*_gen5_03_original_on_tissue_log_res0.8.rds`, `*_gen5_03_original_on_tissue_SCT_res0.8.rds`,
`*_gen5_04_noMT_on_tissue_log_res0.8.rds`, `*_gen5_04_noMT_on_tissue_SCT_res0.8.rds`.

The SCT chunk runs on the same object the log chunk just modified, so the SCT file has both
reductions and both cluster columns; `seurat_clusters` reflects whichever ran last.

Note: `output_dir:` and `output_file:` sit at the top level of the YAML instead of under
`html_document:`, so they do nothing. Pass them to `rmarkdown::render()` instead.

### 05.scType.R

Automated cell typing across all ten samples. Runs *before* the manual labelling step.

- Sources `gene_sets_prepare.R` and `sctype_score_.R` from the sc-type GitHub raw URLs, so
  it needs network access at runtime.
- Globs every `*_gen5_05_noMT_on_tissue_NICHES.rds` recursively, reads them into `rds_list`,
  names each element after its parent folder.
- Builds a `stats_table` of mean nCount/nFeature. Diagnostic only, not used again.
- Positive gene sets come from a local `mega_data_set.xlsx` (not in the repo). Negatives are
  empty except `Macrophage`, which gets `NCAM1`. `"Endothel"` is renamed `"Endothelial"`.
  `gs_list_DB` is built from the public ScTypeDB but never used.
- Runs `sctype_score` on `SCT@scale.data`, sums per `seurat_clusters`, takes the top type per
  cluster, and marks it `"Unknown"` if the score is below `0.1 * ncells` (`cutoff_proportion`,
  line 71). Writes the result to the `scType_SCT` metadata column.
- Saves `gen5_rdslist_06_scType.rds`.

Important: `list.files(recursive = TRUE)` returns paths alphabetically, and steps 03 and 06
index this list positionally assuming 1 = P1_Core through 10 = P5_Peripheral. That only works
if the sample subfolders are named by SampleID (`P1_Core`, ...). Steps 01 and 02 write into
ExperimentID-named folders, which would sort differently. Either name the folders by SampleID
or reorder `rds_list` before saving.

### 03.manual_celltype_GLSB_20250212_*.Rmd (10 files)

Interactive notebooks, not batch scripts. The header says "run the below one by one, find the
best PCA dim and clusters." Knitting straight through will not reproduce the published labels,
because the PC sets and resolutions in each file are the *result* of looking at the plots.

Each file does roughly:

1. Load `samp_list[[i]]`, set DefaultAssay to SCT, re-run ScaleData and FindVariableFeatures.
2. RunPCA (`pca.ph.SCT`), then ElbowPlot and two DimHeatmap panels — you read these to pick PCs.
3. RunUMAP / FindNeighbors on the chosen PCs, FindClusters at a chosen resolution into
   `test_SCT_snn_res.<res>`, FindAllMarkers with a `ratio = pct.1/pct.2` column added.
4. Find and drop the low-count technical artifact population (either a specific cluster or an
   `nCount_SCT` threshold).
5. RunPCA again (`pca.ph.SCT.v2`), new UMAP, new clusters at a lower resolution, new markers.
6. Look at markers with FeaturePlot / SpatialFeaturePlot / SpatialDimPlot(cells.highlight).
7. `RenameIdents` to biological names, assign `CellType_GLSB_manual`.
8. `new_samp_list[[i]] <- samp`.

Per-sample settings (round 1 PCs / res, artifact removed, round 2 PCs / res):

- P1_Core (1): 2:8,10,16 / 0.15, cluster 2, 1:9 / 0.2
- P1_Peripheral (2): 2:7,9,10,18 / 0.25, cluster 2, 1:11 / 0.2
- P2_Core (3): 2:9,11 / 1, cluster 8, 1:10 / 0.2
- P2_Peripheral (4): 2:7,9,13 / 1, cluster 8, 1:10 / 0.3
- P3_Core (5): 2:8 / 1, nCount_SCT <= 3200, 1:11 / 0.3
- P3_Peripheral (6): 2:8,14 / 1, cluster 8, 1:8 / 0.2
- P4_Core (7): 2:8,10,11,13 / 1, nCount_SCT <= 7000, 1:9 / 0.5
- P4_Peripheral (8): 2:8,11 / 0.6, cluster 4, 1:9,11 / 0.1
- P5_Core (9): 2:10 / 1, cluster 8, 1:11 / 0.5
- P5_Peripheral (10): 2:10 / 1, cluster 11 then cluster 4, 1:9,11 then 1:8 / 0.5 then 0.4

P5_Peripheral is the only one needing a third round (`pca.ph.SCT.v3`, dims 1:8, res 0.4).

Final `CellType_GLSB_manual` labels:

- P1_Core: Oligodendrocytes, Invasive-front, MES-like
- P1_Peripheral: Invasive-front, Neuron, Astrocytes
- P2_Core: Invasive-front, Neuron, Astrocytes
- P2_Peripheral: RBFOX3+, Oligodendrocytes (x2), NPC-like, Perivascular, Macrophages/Microglia
- P3_Core: Invasive-front, NPC-like (x2), MES-like, Endothelial
- P3_Peripheral: Neuron (x2), Unidentified, NPC-like
- P4_Core: AC-like, Invasive-front (x3), MES-like_hypoxia, Endothelial
- P4_Peripheral: NPC-like, Invasive-front
- P5_Core: NPC-like, Astrocytes, Macrophage/Microglia, MES-like, Inflammation, Endothelial
- P5_Peripheral: AC-like, Inflammation, Neuron, NPC-like, Perivascular

Gotchas:

- Only P1_Core has an active `samp_list <- readRDS(...)` and an active `new_samp_list <- list()`.
  The other nine assume both already exist, so run all ten in one session, in index order,
  P1_Core first.
- The first chunk calls `rm(samp)` / `rm(mark)`, which errors in a clean session.
- None of the ten saves anything. After the tenth you have to
  `saveRDS(new_samp_list, "07.maunally_labeled_GLSB_v2_20250218.rds")` yourself — that exact
  filename, typo included, is what step 06 reads.

### 06.Boundary_label_*.R (10 files)

Two jobs. All ten do the first; four also do the second.

**Job A — removing failed microfluidic channels (all ten).**

Blank rows/columns in the 50x50 grid create fake gaps that corrupt spatial neighbour
calculations. The fix:

1. Take centroid coords from `samp@images$<ExperimentID>@boundaries$centroids@coords`.
2. Bin x and y into 50 intervals each.
3. Look at `SpatialDimPlot(samp)`, note which bins are empty, hard-code them into
   `bins_to_remove_row` / `bins_to_remove_col`.
4. Filter those spots out.
5. Close the gaps: for each remaining spot, count removed bins below it (rows) or left of it
   (columns) and shift its coordinate by that count times the bin size. This slides the tissue
   together rather than just deleting spots.
6. Subset the object, write the new coords back into the centroids slot, saveRDS.

Per file (image slot, rows removed, cols removed, output):

- P1_Core: GBM1027C, rows 7 and 9, col 4 → `P1_Core_20250409.rds`
- P1_Peripheral: GBM1027N, no rows, col 4 → no saveRDS in the file
- P2_Core: GBM1031C, no rows, col 4 → saveRDS commented out at line 283
- P2_Peripheral: GBM1031N, no rows, cols 4, 44, 45 → `P2_Peripheral_20250414.rds`
- P3_Core: R1102C, no rows, cols 4, 44 → `P3_Core_20250414.rds`
- P3_Peripheral: R1102N, no rows, cols 4, 44 → `P3_Peripheral_20250414.rds`
- P4_Core: R1110C, row 43, col 4 → no saveRDS in the file
- P4_Peripheral: R1110N, no rows, cols 4, 17, 44, 45 → `P4_Peripheral_20250414.rds`
- P5_Core: R0103N, no rows, cols 4, 47 → `P5_Core_20250414.rds`
- P5_Peripheral: R0103C, rows 47–50, cols 4, 44, 45 → `P5_Peripheral_20250414.rds`

Column 4 fails on every single section.

**Job B — hand-painting regions (P1_Core, P1_Peripheral, P2_Core, P4_Core).**

These embed a Shiny app for selecting regions by hand. P1_Core uses a rectangular brush but
never calls `shinyApp()`, so it's a scaffold only. The other three use plotly with
`dragmode = "lasso"` and do launch.

Workflow: launch the app, lasso a region on the `CellType_GLSB_manual` scatter, copy the
printed x/y coordinates from the console, paste them back into the script as
`selected_x_coords_*` / `selected_y_coords_*`. The script loops over those pairs, matches
them against `samp@meta.data$x` / `$y`, and writes a new metadata column:

- P1_Peripheral → `neuron_gbm`: Other, Neuron-Healthy, Neuron-1A, Neuron-1B, GBM-1A, GBM-1B, GBM-1C
- P2_Core → `neuron_new`: Other, Neuron-Healthy, Neuron-1A, Neuron-1B, Neuron-1C, GBM
- P4_Core → `vascular_tumor`: Other, Invasive-front bulk 1–2, Invasive-front near vasculature 1–3

Each of the three also writes the column out with `write.csv`. These columns are what let
step 04 compare interface niches against bulk tumour.

Bugs:

- `06.Boundary_label_P1_Core.R:75` does `samp <- samp_list[[1]]`, but `samp_list` is never
  defined in that file (line 23 loaded `samples`). Should be `samples[[1]]`.
- `06.Boundary_label_P3_Peripheral.R` never assigns `samp` at all. Add `samp <- samples[[6]]`
  after line 19.
- All ten have `library(layer)`, which isn't on CRAN or Bioconductor and isn't used. Delete it.
- P1_Peripheral, P2_Core and P4_Core have no active saveRDS, yet step 04 reads exactly those
  outputs. Add or uncomment them.

Because the coordinates are pasted in by hand, these files are a *record* of what was selected,
not a repeatable selection procedure. Treat the pasted vectors as data.

### 04.NICHES_*.R (10 files)

The final analysis step.

```r
NICHES_output <- NICHES::RunNICHES(samp,
    assay = 'Spatial', LR.database = 'fantom5', species = 'human',
    position.x = 'x', position.y = 'y',
    k = 9, rad.set = NULL,
    cell_types = "CellType_GLSB_manual",
    CellToCellSpatial = T, CellToNeighborhood = T, NeighborhoodToCell = T,
    CellToCell = F, SystemToCell = F, CellToSystem = F)
```

Then:

1. Take the `data` layer of `NICHES_output$NeighborhoodToCell`, rename its columns to the
   `ReceivingCell` barcodes so they line up with the parent object, and add it back as an
   assay called `NeighborhoodToCell_SCT_9`.
2. Set that assay default, then ScaleData → FindVariableFeatures(selection.method="disp") →
   RunPCA → RunUMAP(1:20) → FindNeighbors(1:15) → FindClusters at a per-sample resolution.
3. FindAllMarkers(logfc.threshold=0.25, min.pct=0.1), plus `ratio = pct.1/pct.2` and in some
   files `power = ratio * avg_log2FC`.
4. saveRDS.

Of the three NICHES outputs requested, only `NeighborhoodToCell` is grafted back as an assay,
because it's per-receiving-cell and therefore alignable with the spatial object.

Per file (input, resolution, output):

- P1_Core: `P1_Core_20250409.rds`, 0.3 → `P1_Core_20250414_NICHES.rds`
- P1_Peripheral: `P1_Peripheral_20250409.rds`, 0.3 → `P1_Peripheral_20250414_NICHES.rds`
- P2_Core: `P2_Core_20250409.rds`, 0.3 → `P2_Core_20250414_NICHES.rds`
- P2_Peripheral: `P2_Peripheral_20250414.rds`, 0.2 → `P2_Peripheral_20250414_NICHES.rds`
- P3_Core: `P3_Core_20250414.rds`, 0.15 → `P3_Core_20250414_NICHES.rds`
- P3_Peripheral: `P3_Peripheral_20250414.rds`, 0.1 → `P3_Peripheral_20250414_NICHES.rds`
- P4_Core: `P4_Core_20250409.rds`, 0.2 → `P4_Core_20250414_NICHES.rds`
- P4_Peripheral: `P4_Peripheral_20250414.rds`, 0.1 → `P4_Peripheral_20250414_NICHES.rds`
- P5_Core: `P5_Core_20250414.rds`, 0.1 → `P5_Core_20250414_NICHES.rds`
- P5_Peripheral: `P5_Peripheral_20250414.rds`, 0.1 → `P5_Peripheral_20250414_NICHES.rds`

Three files have extra code beyond the common template:

- `04.NICHES_P4_Core.R` also analyses the standalone `NeighborhoodToCell` object on its own
  (scale, variable features, PCA, UMAP) to make the "Cellular Microenvironment" embedding.
- `04.NICHES_P2_Core.R` builds an `Interface` column by folding the Shiny `Neuron-GBM`
  selection into `CellType_GLSB_manual`, re-runs FindAllMarkers grouped by it, and correlates
  `MAP2` against the `CALM2—GRM7` NICHES signal on non-zero spots.
- `04.NICHES_P1_Peripheral.R` and `04.NICHES_P2_Core.R` both have an EnhancedVolcano block
  reading `03.25.findmarkers_neuron_v_other.txt`, a ranked-ratio ggplot, and SpatialFeaturePlot
  panels for named pairs (NLGN3—NRXN1/2/3, SEMA5A—MET, VEGFA—ITGA9, THBS2—CD47, HGF—MET).

Notes:

- NICHES pair names use an em dash (—, U+2014), not a hyphen. Copy them, don't retype them.
- The comment on the `k = 9` line says "10 nearest neighbors" — that's 9 neighbours plus the
  cell itself. The parameter is right, the comment counts inclusively.
- `BiocManager::install(...)` is inline in two files. Install once separately and comment out.

### Duplicate files

`01.gen5_rds_generation_20241206 2.R`, `... 3.R`, `05.scType 2.R`, `05.scType 3.R`,
`02.gen5_rds_normalization_20241207 2.Rmd`, `... 3.Rmd` are byte-identical copies from a file
sync client. Use the un-suffixed ones.

---

## How to run it

### Packages

```r
install.packages(c("Seurat","SeuratObject","tidyverse","dplyr","magrittr","patchwork",
                   "cowplot","data.table","ggplot2","ggpubr","colorspace","stringr",
                   "purrr","janitor","R.utils","shiny","plotly","openxlsx","HGNChelper",
                   "scales","knitr","rmarkdown"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

remotes::install_github("msraredon/NICHES", ref = "master")   # not on CRAN
```

R >= 4.2 with Seurat v5. Don't install `layer`. If your compute nodes have no outbound
network, download `gene_sets_prepare.R` and `sctype_score_.R` from the sc-type repo and
source them locally instead of over HTTPS.

Every path in every script points at a Yale HPC (`/gpfs/gibbs/pi/fan/sb2723/...`,
`/vast/palmer/...`) and needs rewriting.

### Files you have to supply

- `barcodes-AB.xls` — barcode, y.cord, x.cord for the DBiT grid (used by 01)
- a 10x-format image folder per sample at `02.10x_fmt_img/<ID>/spatial/` — tissue positions
  csv, scalefactors json, hires/lowres pngs (used by 01)
- `mega_data_set.xlsx` — the scType marker workbook (used by 05)
- `03.25.findmarkers_neuron_v_other.txt` — for the volcano plots in 04
- `x` and `y` metadata columns — see below

### The missing x / y columns

`04.NICHES_*.R` passes `position.x = 'x'`, `position.y = 'y'`, and the step-06 Shiny apps plot
and match against `samp@meta.data$x` / `$y`. No committed script creates these columns. Add
them from the image centroids:

```r
coords <- samp@images[[<ExperimentID>]]@boundaries$centroids@coords
samp$x <- coords[, 1]
samp$y <- coords[, 2]
```

Add them **after** step 06's gap-closing shift, otherwise NICHES computes neighbours on the
uncorrected grid, which defeats the point of step 06.

### Step 0 — zUMIs

Run zUMIs per sample using the YAML template (you'll need an unencrypted copy). The pipeline
consumes `<zUMIs_out>/<ID>/zUMIs_output/expression/<ID>.dgecounts.rds`, specifically
`umicount$inex$all`.

### Step 1 — run 01

`sample` has to come from outside the script:

```r
samples <- c("GBM1027C","GBM1027N","GBM1031C","GBM1031N",
             "R1102C","R1102N","R1110C","R1110N","R0103C","R0103N")
for (sample in samples) source("01.gen5_rds_generation_20241206.R", local = new.env())
```

or per sample on a scheduler:

```bash
Rscript -e 'sample <- commandArgs(TRUE)[1]; source("01.gen5_rds_generation_20241206.R")' GBM1027C
```

Check the two QC pngs per sample: nCount_Spatial should follow tissue structure, not the
channel grid.

### Step 2 — run 02

```r
for (s in samples) {
  rmarkdown::render("02.gen5_rds_normalization_20241207.Rmd",
    params      = list(sample_name = s),
    output_file = paste0(s, "_gen5_normalization.html"),
    output_dir  = "<your_root>/04.rds_gen5_normalized_20241207/html_report")
}
```

Pass `output_file` / `output_dir` to `render()`, not via the YAML.

### Step 3 — the missing NICHES pass

`05.scType.R` needs files named `*_gen5_05_noMT_on_tissue_NICHES.rds`. Reconstruct them by
running the `RunNICHES` block from any `04.NICHES_*.R` against
`*_gen5_04_noMT_on_tissue_SCT_res0.8.rds`.

If you just want the published result, this pass barely matters — its only contribution
downstream is the `seurat_clusters` column scType scores against, and step 03 re-clusters
from scratch. Renaming the step-02 SCT object to match the glob works.

Name the sample subfolders by SampleID (`P1_Core`, `P1_Peripheral`, ...), not ExperimentID,
so the alphabetical ordering matches what steps 03 and 06 assume.

### Step 4 — run 05

Put `mega_data_set.xlsx` in the working directory, fix the paths, run once. Then check:

```r
l <- readRDS("gen5_rdslist_06_scType.rds")
length(l)                 # 10
names(l)                  # P1_Core ... P5_Peripheral in that order
table(l[[1]]$scType_SCT)  # should not be all "Unknown"
```

If everything comes back "Unknown", lower `cutoff_proportion` on line 71. If the order is
wrong, fix it now — everything after this indexes positionally.

### Step 5 — run the ten 03 notebooks

One R session, in index order, P1_Core first. Work chunk by chunk:

1. Run the PCA chunk, read the ElbowPlot and DimHeatmaps, pick PCs. The published choices
   mostly start at PC 2 — PC 1 is dominated by sequencing depth here.
2. UMAP / neighbours / clusters at your resolution.
3. Check `VlnPlot(nCount_SCT, group.by = ...)`. A cluster with low counts and no coherent
   spatial footprint is a technical artifact — drop it.
4. Re-run PCA / UMAP / clusters on the cleaned object at a lower resolution.
5. Inspect markers (FindAllMarkers ranked by `ratio`, SpatialFeaturePlot,
   SpatialDimPlot(cells.highlight)).
6. RenameIdents, assign `CellType_GLSB_manual`.
7. `new_samp_list[[i]] <- samp`.

To reproduce the paper, use the per-sample settings listed above — they're already in the
files. Then save it yourself:

```r
saveRDS(new_samp_list, file = "07.maunally_labeled_GLSB_v2_20250218.rds")
```

### Step 6 — run the ten 06 scripts

1. `samp <- samples[[i]]`. Fix P1_Core and P3_Peripheral first (see bugs above).
2. Run the binning block, look at `SpatialDimPlot(samp)`, identify blank rows/cols as 1-based
   bin indices in a 50x50 grid. Expect column 4 every time.
3. Set `bins_to_remove_row` / `bins_to_remove_col`. Use `NULL`, not `c()`, when there are none.
4. Run the filter/shift/write-back block and re-plot — the tissue should be gap-free with no
   visible discontinuity.
5. saveRDS. Add this to P1_Peripheral, P2_Core and P4_Core.

For the four with region painting: launch `shinyApp(ui, server)`, lasso the regions, copy the
coordinates into the `selected_*_coords_*` vectors, run the metadata loops, check with
SpatialDimPlot, then write.csv.

Compare `dim(samp)` before and after to see how many spots the channel removal cost.

### Step 7 — run the ten 04 scripts

Check that `x`/`y` exist and `CellType_GLSB_manual` is populated first. `RunNICHES` is the
expensive call — run it on a compute node.

```r
samp <- readRDS("P1_Core_20250414_NICHES.rds")
Assays(samp)                                        # includes NeighborhoodToCell_SCT_9
dim(samp[["NeighborhoodToCell_SCT_9"]])             # LR pairs x spots
head(rownames(samp[["NeighborhoodToCell_SCT_9"]]))  # e.g. "CALM2—GRM7"
SpatialFeaturePlot(samp, "NLGN3—NRXN1")
```

### Step 8 — surprisal analysis

Runs in https://github.com/ShouryoGhose/Surprisal-Analysis-of-Human-GBM-samples on the raw
count matrix — the `counts` layer of the `Spatial` assay from step 01, not log-normalised or
SCT values.

```r
samp <- readRDS("<ID>_gen5_01_original_on_tissue.rds")
mat  <- GetAssayData(samp, assay = "Spatial", layer = "counts")
write.csv(as.matrix(mat), "P1_Core_raw_counts.csv")
```

Decide deliberately whether to use the `_01_original_` (all genes) or `_02_noMT_` (filtered)
object. The step-01 regex strips ribosomal and unnamed Ensembl genes too, which changes the
count distribution. The README says raw expression but doesn't say which of the two.

---

## Known issues, collected

- `01.R:148–149` — filtered object subset with unfiltered barcodes, `img_noMT` unused
- `01.R:60–62` — `message` overwritten twice, only SampleType prints
- `01.R:117` — regex removes far more than MT genes but the printout says "MT"
- `02.Rmd:6–7` — `output_dir`/`output_file` at top level of YAML, inert
- `05.scType.R:45` — `gs_list_DB` computed and never used
- `05.scType.R:17–20` — list order must be P1_Core...P5_Peripheral; depends on folder naming
- `03.*.Rmd:22–23` — `rm(samp)` errors in a clean session
- `03.*.Rmd` — none of them saves `new_samp_list`
- `06.P1_Core.R:75` — `samp_list` undefined, should be `samples`
- `06.P3_Peripheral.R` — `samp` never assigned
- `06.*.R` — `library(layer)` in all ten, not a real package here
- `06.P1_Peripheral / P2_Core / P4_Core` — no active saveRDS, but step 04 reads their outputs
- `04.NICHES_P1_Peripheral / P2_Core` — inline `BiocManager::install`
- all `.zip` files are password-encrypted; the zUMIs yaml has no plain-text copy
- `* 2.*` and `* 3.*` files are byte-identical duplicates
