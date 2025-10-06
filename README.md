# Assessment 4 – R Bioinformatics Project

> **Unit:** Applied Bioinformatics
> **Report title (R Markdown):** *Assessment_4_R_Project* (HTML output)
> **Author:** **N. Koppert**
  > **Repository purpose:** Public, reproducible workflow demonstrating data wrangling, statistical analysis, and sequence analysis in R for Assessment 4 (Parts 1 & 2).

---
  
## Purpose of this Project
  
  This repository provides a fully reproducible R Markdown workflow that implements **all required tasks for Assessment 4**:
  
  * **Part 1 (Q1–Q10):** Tabular and graphical analysis of two datasets:
  
  * RNA-seq counts (`gene_expression.tsv`): import, summary statistics, ranking by mean expression, thresholding, and histogramming (log-transform included).
* Tree growth (`growth_data.csv`): import, descriptive statistics per site and year, three comparative boxplots, 10-year growth computation, and a **Welch two-sample t-test** after variance diagnostics.

* **Part 2 (Q1–Q6):** Comparative **genome coding sequence** (CDS) analysis for *Escherichia coli* K-12 MG1655 vs **Corynebacterium bovis**:
  
  * Programmatic retrieval from **Ensembl Genomes** (FTP), decompression, and parsing of CDS FASTA files.
* Genome‐scale summaries (CDS counts; total coding DNA).
* CDS length statistics (mean/median) and distribution visualisation.
* **Nucleotide** and **amino-acid** frequency profiling with comparative bar charts.
* **Codon usage** quantification via **RSCU** (Relative Synonymous Codon Usage) and heatmap/bar visualisations of codon-bias differences.
* Data-driven identification and visualisation of **over- and under-represented protein k-mers** (3–5 aa) per organism, with biological interpretation.

---

## Inputs / Outputs

**P1 Q1**

* **Inputs:** URL to `gene_expression.tsv`.
* **Outputs:** Markdown table showing the first six genes (via `kable`).

**P1 Q2**

* **Inputs:** In-memory object `raw_gene_expression`.
* **Outputs:** Markdown table of `head()` including the appended mean column (via `kable`).

**P1 Q3**

* **Inputs:** In-memory object from Q2.
* **Outputs:** Markdown table of the top 10 genes ranked by mean expression (via `kable`).

**P1 Q4**

* **Inputs:** In-memory mean expression vector.
* **Outputs:** Scalar count printed to console.

**P1 Q5**

* **Inputs:** In-memory mean expression vector.
* **Outputs:** Figure `p1_q5_hist_mean.png` (histogram of mean expression).

**P1 Q6**

* **Inputs:** URL to `growth_data.csv`.
* **Outputs:** Column names printed (structure check).

**P1 Q7**

* **Inputs:** In-memory `growth_data`.
* **Outputs:** Markdown table of mean ± SD by site × year.

**P1 Q8**

* **Inputs:** In-memory `growth_data`.
* **Outputs:** Figure `p1_q8_dbh_boxplot.png` (comparative box plots).

**P1 Q9**

* **Inputs:** In-memory `growth_data`.
* **Outputs:** Markdown table of 10-year growth by site.

**P1 Q10**

* **Inputs:** In-memory growth summaries per site.
* **Outputs:** Welch two-sample t-test output (including p-value).

**P2 Q1**

* **Inputs:** Ensembl FTP URLs for CDS FASTA files.
* **Outputs:** Markdown table of CDS counts for the two organisms.

**P2 Q2**

* **Inputs:** Decompressed `.fa` files.
* **Outputs:** Markdown table of total coding DNA (bp) per organism.

**P2 Q3**

* **Inputs:** `.fa` files parsed to CDS records.
* **Outputs:** Figure `p2_q3_cds_len_boxplot.png` and reported mean/median CDS lengths.

**P2 Q4**

* **Inputs:** `.fa` files (nucleotide and translated amino-acid sequences).
* **Outputs:** Bar-chart figures for nucleotide and amino-acid frequencies.

**P2 Q5**

* **Inputs:** `.fa` files (coding sequences).
* **Outputs:** RSCU table(s) and associated plot(s) of codon-usage bias.

**P2 Q6**

* **Inputs:** Protein sequences derived from CDS.
* **Outputs:** Plot(s) of most over- and under-represented 3–5 aa k-mers per organism.


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