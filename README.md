# Small particle aerosol infection of nonhuman primates with Lassa virus identifies previremic biomarkers and correlates of survival

[![DOI](https://zenodo.org/badge/1033915430.svg)](https://doi.org/10.5281/zenodo.19359338)

> Manuscript under review. Citation TBD.

## Analysis methods

All analyses were performed in R v4.5.3. Analyses other than Nanostring and LEGENDplex were performed with [`rstatix v0.7.3`](https://rpkgs.datanovia.com/rstatix/) and can be found in [`analysis.r`](analysis.r). Nonparametric tests were used where possible.

### Transcriptomics (Nanostring)

Targeted transcriptomic analysis of the expression of 770 host mRNAs was quantified via the Nanostring NHP Immunology v2 panel (Bruker, #531-115000276). The concentration and purity of whole blood RNA extracted from AVL following RT-qPCR were assessed using a NanoDrop spectrophotometer (Thermo Fisher Scientific). Samples were then diluted per the manufacturer’s instructions. Samples with <20 ng/µL RNA were used as-is. All samples were then prepared using the nCounter Prep Station and assayed via the nCounter Pro Analysis System (Bruker) according to the manufacturer’s instructions. 
Output files were loaded into [`nSolver v4.0`](https://brukerspatialbiology.com/products/ncounter-analysis-system/ncounter-analysis-solutions/), and background thresholding using several housekeeping gene abundances was performed using the default parameters. Thresholded normalized count matrices were exported from nSolver as a CSV file and analyzed with [`limma v3.66.0`](https://doi.org/10.1093/nar/gkv007). Principal component analysis figures, volcano plots, and gene expression figures were made using [`ggplot2 v.4.0.2`](https://ggplot2.tidyverse.org). Normalized data (log fold change values and FDR-adjusted p-values) for early and late gene expression were exported as CSV files for [Ingenuity Pathway Analysis (IPA)](https://doi.org/10.1093/bioinformatics/btt703)-based (QIAGEN) pathway enrichment analysis. Significantly upregulated or downregulated mRNAs were filtered by p<sub>adj</sub> < 0.05 and log<sub>2</sub> fold change > 1 or < -1. Z-scores for enriched pathways were imported into R and visualized with a radar plot using [`fmsb v0.7.6`](https://doi.org/10.32614/CRAN.package.fmsb).  Code for Nanostring analyses can be found in [`nanostring.r`](nanostring.r).

### Proteomics (LEGENDplex)

Circulating cytokines, thrombosis markers, and fibrinolysis analytes were measured by LEGENDplex bead-based immunoassays (BioLegend). Gamma-irradiated plasma samples were assessed in duplicate for each NHP Inflammation (#741491, 1:4 dilution), Human Thrombosis (#740892, 1:100 dilution), NHP chemokine/cytokine (#740388, 1:4 dilution), Human Immune Checkpoint (#740962, 1:2 dilution), Human Vascular Inflammation (#740590, 1:10,000) and Human Fibrinolysis (#740761, 1:40,000 dilution) panel according to the manufacturer instructions. Assay standards were mixed in batches and aliquoted across all plates to ensure consistency. Optional wash steps were incorporated to reduce the background signal. Assay samples were analyzed on FACS Canto-II or Accuri C6 Plus flow cytometers (BD Biosciences). Raw `.fcs` files from each assay were imported into LEGENDplex Qognit cloud-based Data Analysis Software suite (BioLegend), which automatically determined analyte concentrations in experimental samples following 5-parameter logistic regression curve fitting of assay standards. Analyte concentration data from assayed plasma samples were exported from Qognit and analyzed in R using [`limma v3.66.0`](https://doi.org/10.1093/nar/gkv007). Principal-component analysis figures, volcano plots, and protein abundance figures were made using [`ggplot2 v.4.0.2`](https://ggplot2.tidyverse.org). Significantly upregulated or downregulated proteins were filtered by p<sub>adj</sub> < 0.05 and log<sub>2</sub> fold change > 1 or < -1. Code for LEGENDplex analyses can be found in [`legendplex.r`](legendplex.r).

## Data availability

The LV0043 genome consensus sequence is deposited in GenBank ([PX279630](https://www.ncbi.nlm.nih.gov/nuccore/PX279630) and [PX279631](https://www.ncbi.nlm.nih.gov/nuccore/PX279631)). Nanostring transcriptomic data can be found in NCBI GEO (ACCESSION). Circulating protein data can be found HERE. 
