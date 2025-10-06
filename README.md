# Assessment 4 – R Bioinformatics Project

**Unit:** Applied Bioinformatics

**Report title (R Markdown):** *Assessment_4_R_Project* (HTML output)

**Author:** **N. Koppert**

**Repository purpose:** Public, reproducible workflow demonstrating data wrangling, statistical analysis, and sequence analysis in R for Assessment 4 (Parts 1 & 2).

---
  
## Purpose of this Project
  
  This repository provides a fully reproducible R Markdown workflow that implements **all required tasks for Assessment 4**:
  
**Part 1 (Q1–Q10):** Tabular and graphical analysis of two datasets:
  
* RNA-seq counts (`gene_expression.tsv`): import, summary statistics, ranking by mean expression, thresholding, and histogramming (log-transform included).
* Tree growth (`growth_data.csv`): import, descriptive statistics per site and year, three comparative boxplots, 10-year growth computation, and a **Welch two-sample t-test** after variance diagnostics.

**Part 2 (Q1–Q6):** Comparative **genome coding sequence** (CDS) analysis for *Escherichia coli* K-12 MG1655 vs **Corynebacterium bovis**:
  
* Programmatic retrieval from **Ensembl Genomes** (FTP), decompression, and parsing of CDS FASTA files.
* Genome‐scale summaries (CDS counts; total coding DNA).
* CDS length statistics (mean/median) and distribution visualisation.
* **Nucleotide** and **amino-acid** frequency profiling with comparative bar charts.
* **Codon usage** quantification via **RSCU** (Relative Synonymous Codon Usage) and heatmap/bar visualisations of codon-bias differences.
* Data-driven identification and visualisation of **over- and under-represented protein k-mers** (3–5 aa) per organism, with biological interpretation.

---

## Inputs / Outputs (by code chunk)

### **Setup and Package Loading**

**Chunk: `setup`**

* **Inputs:** None (initialization).
* **Outputs:** Global `knitr` options set (`echo = TRUE`).

**Chunk: `load_packages`**

* **Inputs:** None.
* **Outputs:** Required CRAN packages installed (if missing) and loaded: `seqinr`, `R.utils`, `ggplot2`, and `knitr`.

---

### **Part 1 — Gene Expression Analysis**

**Chunk: `import_gene_expression`**

* **Inputs:** URL `"https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv"`.
* **Outputs:** Data frame `raw_gene_expression`; Markdown table of the first six genes rendered via `kable()`.

**Chunk: `calculate_mean_expression`**

* **Inputs:** Object `raw_gene_expression`.
* **Outputs:** Data frame `gene_expression` (with new column `mean_expression`); Markdown table of first six rows.

**Chunk: `top10_mean_expression`**

* **Inputs:** Object `gene_expression`.
* **Outputs:** Markdown table listing the 10 genes with the highest mean expression values.

**Chunk: `count_low_expression_genes`**

* **Inputs:** `genes_sorted_by_mean$mean_expression`.
* **Outputs:** Console scalar output: number of genes with mean expression < 10.

**Chunk: `plot_histogram_raw_means`**

* **Inputs:** `gene_expression$mean_expression`.
* **Outputs:** Figure — histogram of untransformed mean expression values.

**Chunk: `plot_histogram_log_means`**

* **Inputs:** `gene_expression$mean_expression`.
* **Outputs:** Figure — histogram of log-transformed mean expression values.

---

### **Part 1 — Tree Growth Analysis**

**Chunk: `import_growth_data`**

* **Inputs:** URL `"https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv"`.
* **Outputs:** Data frame `raw_growth_data`; column names printed to console.

**Chunk: `growth_summary_stats`**

* **Inputs:** `raw_growth_data`.
* **Outputs:** Data frame `summary_table` containing mean and SD of circumference by site and year; rendered Markdown table.

**Chunk: `growth_boxplot_2005`**

* **Inputs:** `growth_boxplot_data` subset for `Year == "2005"`.
* **Outputs:** Figure — box plot comparing circumference between sites (2005).

**Chunk: `growth_boxplot_2020`**

* **Inputs:** `growth_boxplot_data` subset for `Year == "2020"`.
* **Outputs:** Figure — box plot comparing circumference between sites (2020).

**Chunk: `growth_boxplot_combined`**

* **Inputs:** `growth_boxplot_data` for both years.
* **Outputs:** Figure — combined box plot comparing 2005 vs 2020 for both sites.

**Chunk: `growth_10yr_mean`**

* **Inputs:** Subsets `northeast` and `southwest`.
* **Outputs:** Markdown table of mean 10-year growth (2010–2020) per site.

**Chunk: `variance_ratio_test`**

* **Inputs:** Growth vectors `northeast_growth` and `southwest_growth`.
* **Outputs:** Console output of variance ratio (F) test.

**Chunk: `welch_t_test`**

* **Inputs:** Growth vectors `northeast_growth` and `southwest_growth`.
* **Outputs:** Console output of Welch’s two-sample t-test (mean growth comparison).

---

### **Part 2 — Biological Sequence Analysis**

**Chunk: `download_ecoli_cds`**

* **Inputs:** Ensembl FTP URL for *E. coli* CDS FASTA file.
* **Outputs:** Downloaded file `ecoli_cds.fa.gz`; decompressed file `ecoli_cds.fa`; console confirmation via `list.files()`.

**Chunk: `download_cbovis_cds`**

* **Inputs:** Ensembl FTP URL for *C. bovis* CDS FASTA file.
* **Outputs:** Downloaded file `cbovis_cds.fa.gz`; decompressed file `cbovis_cds.fa`; console confirmation via `list.files()`.

**Chunk: `read_fasta_sequences`**

* **Inputs:** Files `ecoli_cds.fa` and `cbovis_cds.fa`.
* **Outputs:** FASTA objects `cds_ecoli` and `cds_cbovis`.

**Chunk: `count_cds_entries`**

* **Inputs:** `cds_ecoli`, `cds_cbovis`.
* **Outputs:** Markdown table of CDS counts per organism (`kable()`).

**Chunk: `compare_total_coding_dna`**

* **Inputs:** `cds_ecoli`, `cds_cbovis`.
* **Outputs:** Markdown table of total coding DNA (bp) per organism (`kable()`).

**Chunk: `analyze_cds_length`**

* **Inputs:** `cds_ecoli`, `cds_cbovis`.
* **Outputs:** Markdown table of mean and median CDS lengths (`kable()`).

**Chunk: `plot_cds_length_boxplot`**

* **Inputs:** Vectors `len_ecoli` and `len_cbovis`.
* **Outputs:** Figure — box plot of CDS length distributions.

**Chunk: `nucleotide_frequency_plot`**

* **Inputs:** `cds_ecoli`, `cds_cbovis`.
* **Outputs:** Figure — ggplot bar chart comparing nucleotide composition (A, T, G, C).

**Chunk: `amino_acid_frequency_plot`**

* **Inputs:** Translated protein sequences derived from `cds_ecoli` and `cds_cbovis`.
* **Outputs:** Figure — ggplot bar chart comparing amino acid composition.

---

### **Part 2 — Codon and k-mer Analyses**

**Chunk: `codon_usage_compute`**

* **Inputs:** `cds_ecoli`, `cds_cbovis`.
* **Outputs:**

  * Two RSCU tables (`kable()` for *E. coli* and *C. bovis*).
  * Combined `codon_usage_all` data frame for plotting.

**Chunk: `codon_usage_heatmap`**

* **Inputs:** `codon_usage_all`.
* **Outputs:** ggplot heatmap showing codon usage bias across organisms.

**Chunk: `codon_usage_barplot`**

* **Inputs:** `codon_usage_all`.
* **Outputs:** ggplot bar chart of RSCU values per codon for both species.

**Chunk: `define_kmer_functions`**

* **Inputs:** None.
* **Outputs:** Three user-defined functions: `get_kmer_extremes()`, `create_kmer_df()`, and `plot_kmer_bar()`.

**Chunk: `compute_kmer_frequencies`**

* **Inputs:** Protein vectors `prot_ecoli_all` and `prot_cbovis_all`.
* **Outputs:** Lists of k-mer frequency data (`freq_ecoli`, `freq_cbovis`), extracted top/bottom sets, and merged data frames (e.g., `ecoli_kmer_3`, `cbovis_kmer_4`, etc.).

**Chunks: `plot_ecoli_kmer3`, `plot_cbovis_kmer3`, `plot_ecoli_kmer4`, `plot_cbovis_kmer4`, `plot_ecoli_kmer5`, `plot_cbovis_kmer5`**

* **Inputs:** Corresponding k-mer data frames (`ecoli_kmer_n`, `cbovis_kmer_n`).
* **Outputs:** Six ggplot bar charts showing over- and under-represented k-mers (3–5 aa) for both species.

**Chunk: `combine_kmer_tables`**

* **Inputs:** Individual k-mer data frames.
* **Outputs:** Two Markdown tables summarising top and bottom 10 k-mers for *E. coli* and *C. bovis* (via `kable()`).

**Chunk: `session_info`**

* **Inputs:** None.
* **Outputs:** Console listing of R environment details and loaded package versions.

---

## Files included
  
The analysis downloads external data at runtime.

* `Assessment_4_R_Project.Rmd`
Self-contained R Markdown source implementing all steps for Parts 1 & 2, including figure/table generation and explanation.

* `README.md` (this file)
Comprehensive documentation (purpose, contents, how to run, contributors, reproducibility, and references).

 **Generated during execution (not committed):**
  
* `gene_expression.tsv` and `growth_data.csv` are **fetched from GitHub** at run time (Part 1).
 * `ecoli_cds.fa.gz` → `ecoli_cds.fa` and `cbovis_cds.fa.gz` → `cbovis_cds.fa` are **downloaded from Ensembl Genomes FTP** and decompressed locally (Part 2).
 * HTML report, figures, and tables are created by knitting the Rmd.

---
  
## How to run the code
  
### RStudio (recommended)
  
1. Open the project folder in **RStudio**.
2. Open `Assessment_4_R_Project.Rmd`.
3. Ensure internet access (files are downloaded programmatically).
4. Click **Knit** → **Knit to HTML**.
5. Can take a minute or two to knit to HTML due to table formatting with the knitr package

The Rmd:
  
* Sets a CRAN mirror and installs/loads required packages as needed.
* Downloads data and CDS FASTA files to the **working directory**.
* Produces an HTML report with all figures, tables, and statistical outputs.


### System & package requirements

* **R ≥ 4.0** (tested with the version shown in `sessionInfo()` at the end of the Rmd).
* Packages (auto-installed by the Rmd):
  `seqinr`, `R.utils`, `ggplot2`, `knitr`

---
  
## Contributors
  
  * **N. Koppert**
  
  ---
  
## Additional information
  
### Methodological notes
  
  * **Part 1.
* ** RNA-seq counts are treated as **raw integers**; histograms of mean counts therefore exhibit **right-skewness** typical of count data. A log(·+1) transform is used solely for **distributional visualisation**.
* ** Ten-year growth comparison employs a **variance ratio test** to justify **Welch’s t-test** (heteroscedasticity-robust). This preserves nominal Type-I error under unequal variances and sample sizes.

* **Part 2.**
  
  * CDS lengths are derived from FASTA sequence records (nucleotides), with summary statistics reported alongside a **box-plot** to depict central tendency and dispersion.
* **RSCU** is computed to quantify codon-usage bias and visualised via a **heatmap** and grouped bars; values >1 indicate codon over-representation relative to uniform synonymous usage.
* **k-mer** analyses (3–5 aa) use frequency‐based ranking; bar plots display the most over-/under-represented motifs per organism to reveal proteomic compositional signatures.

### Limitations

* **Scope:** The workflow focuses on descriptive genomics and basic inferential statistics; it does **not** perform genome annotation, orthology mapping, or phylogenetics.

---
  
### References (Software & Tools, Datasets & Primary Literature)
  
Bengtsson, H 2025, *R.utils: Various Programming Utilities* [R package], version 2.13.0, CRAN, viewed 6 October 2025, [https://CRAN.R-project.org/package=R.utils](https://CRAN.R-project.org/package=R.utils).

Bohlin, J, Brynildsrud, O, Vesth, T, Skjerve, E & Ussery, DW 2013, ‘Amino acid usage is asymmetrically biased in AT- and GC-rich microbial genomes’, *PLoS ONE*, vol. 8, no. 7, e69878, doi:10.1371/journal.pone.0069878.

Charif, D & Lobry, JR 2023, *seqinr: Biological Sequences Retrieval and Analysis* [R package], version 4.2-36, CRAN, viewed 6 October 2025, [https://CRAN.R-project.org/package=seqinr](https://CRAN.R-project.org/package=seqinr).

Ensembl Genomes 2025, *Corynebacterium bovis (GCA_003932475) coding DNA sequences (CDS), Release 62* [data set], viewed 6 October 2025, [https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_34_collection/corynebacterium_bovis_gca_003932475/cds/Corynebacterium_bovis_gca_003932475.ASM393247v1.cds.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_34_collection/corynebacterium_bovis_gca_003932475/cds/Corynebacterium_bovis_gca_003932475.ASM393247v1.cds.all.fa.gz).

Ensembl Genomes 2025, *Escherichia coli K-12 (GCA_000005845) coding DNA sequences (CDS), Release 62* [data set], viewed 6 October 2025, [https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz](https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz).

ghazkha n.d., *gene_expression.tsv* [data set], GitHub, viewed 6 October 2025, [https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv](https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv).

ghazkha n.d., *growth_data.csv* [data set], GitHub, viewed 6 October 2025, [https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv](https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv).

Hershberg, R & Petrov, DA 2008, ‘Selection on codon bias’, *Annual Review of Genetics*, vol. 42, pp. 287–299, doi:10.1146/annurev.genet.42.110807.091442.

Knight, RD, Freeland, SJ & Landweber, LF 2001, ‘A simple model based on mutation and selection explains trends in codon and amino-acid usage and GC composition within and across genomes’, *Genome Biology*, vol. 2, no. 4, p. RESEARCH0010, doi:10.1186/gb-2001-2-4-research0010.

R Core Team 2025, *R: A language and environment for statistical computing* [computer software], version 4.5.1, R Foundation for Statistical Computing, Vienna, viewed 6 October 2025, [https://www.r-project.org/](https://www.r-project.org/).

Rocha, EPC & Danchin, A 2002, ‘Base composition bias might result from competition for metabolic resources’, *Trends in Genetics*, vol. 18, no. 6, pp. 291–294, doi:10.1016/S0168-9525(02)02690-2.

Tadeo, X, López-Méndez, B, Trigueros, T, Laín, A, Castaño, D & Millet, O 2009, ‘Structural basis for the amino acid composition of proteins from halophilic archaea’, *PLoS Biology*, vol. 7, no. 12, e1000257, doi:10.1371/journal.pbio.1000257.

Wickham, H, Chang, W, Henry, L, Pedersen, TL, Takahashi, K, Wilke, C, Woo, K, Yutani, H & Dunnington, D 2025, *ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics* [R package], version 4.0.0, CRAN, viewed 6 October 2025, [https://CRAN.R-project.org/package=ggplot2](https://CRAN.R-project.org/package=ggplot2).

Xie, Y 2025, *knitr: A General-Purpose Package for Dynamic Report Generation in R* [R package], version 1.50, CRAN, viewed 6 October 2025, [https://CRAN.R-project.org/package=knitr](https://CRAN.R-project.org/package=knitr).