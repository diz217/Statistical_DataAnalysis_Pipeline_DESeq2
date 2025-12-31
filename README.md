# Data Analysis Pipeline for GEO Accession Data 

## Overview
This project provides an end-to-end data analysis pipeline for GEO (Gene Expression Omnibus) accession datasets.

GEO data are often heterogeneous, inconsistently annotated, and frequently unsuitable for standard analysis workflows without extensive manual cleaning. 
This pipeline is designed to handle such real-world conditions and to analyze any GEO dataset that can be reasonably interpreted.

This pipeline automatically cleans and aligns sample metadata, preprocesses count matrices, constructs appropriate experimental designs, and selects the correct statistical method based on data type and replicate availability. 
Depending on the input, the pipeline generates either full DESeq2 differential expression results or log2 fold-change summaries using limma-based methods.

## Methodology & Design
### Metadata Cleaning 

### Count Matrix Preprocessing Pipeline

### DESeq2 and Limma analysis

## Requirements
### R Environment
- R (>=4.2)
- RStudio 
Required R packages:
- DESeq2
- limma
- edgeR
- GEOquery
Install R dependencies via Bioconductor:
```r
install.packages("BiocManager")
BiocManager::install(c(
      "DESeq2",
      "limma",
      "edgeR",
      "GEOquery"))
```
### Python Environment
- python(>=3.9)
  - **Add python executable to system path**
```
python -m pip install -r requirements.txt
```
## Usage
### Working Directory Assumption
This pipeline assumes the existence of an `/R space/` directory (see in `.R` scripts), where all the source codes need to be placed into. 
This working diretory is expected to be manually created by the user, under the path: `C:/Users/{username}/Documents`, which is the default home directory set by RStudio. 

### Running the Pipeline
The pipeline is executed by the python file `main_v2.py`. Run the pipeline by providing a GEO accession ID:
```bash
python main_v2.py GSE269485
```
Alternatively, if no GSE accession ID is provided, the pipeline will prompt for user input. The input can be pure numbers or numbers mixed with characters.

### Interactive Prompts
While the pipeline is designed for full automation, there are certain (rare) cases where assumptions cannot be made safely and user interaction might be required. 
1. **Metadata cleaning**
   - If more than 20 unique experimental conditions are detected after all the cleaning procedures, the pipeline will prompt user to either:
     - continue the analysis, or
     - stop and inspect sample titles for potential un-cleanable complications
2. **Multiple series matrix files**
   - If more than one series matrix file is detected for a GEO accession, it is because the samples are interpreted by different GPLs (GEO Platform), and they cannot be compared or combined in the analysis. Therefore, the pipeline will prompt the user to select one matrix at a time (0-indexed).
3. **Count matrix column selection in single GSM files**
   - When aggregating GSM-level count matrices, we expect two columns, 1st column being Genes and 2nd column being counts. If there are more columns, the pipeline will prompt the user to check and specify which column represents the column data.

## Outputs
The pipeline generates a set of intermediate and final result files under the working directory for each GEO accession. 
### Intermediate files
- `{GSE_ID}_meta_gsm_condition.tsv` Metadata matrix with GSM IDs, pre-cleaned titles, descriptions if applicable, pose-cleaned conditions. Note, only the GSMs detected in the count matrix files are kept.
- `{GSE_ID}_conditions.tsv` Unique conditions from post-cleaned titles.
- `{GSE_ID}_control_experiments.tsv` Autoamtically inferred case-control pairs from unique conditions. This file is used in pipeline case-control design for statistical analysis.
- `{GSE_ID}_(count|var)_gsm_matrix_{i}.tsv` Preprocessed aggregated count matrix (i-th downloaded file) with Gene IDs as the first column and Sample IDs (GSM) as the first row. Note, the GSMs follow the same order as in the 'geo_accession' column in the metadata tsv.
  - 'count' stands for integer count data and 'val' stands for normalized count data.

### Final analysis results
Depending on the data type and replicate availability, the pipeline generates two or more statistical result tables:
- DESeq2 results:
  - Generated for integer count data with biological replicated experiments.
  - p-values
  - p-adj
  - baseMean
  - logFC
  - lfcSE
- limma/ limma-voom results:
  - Generated when replicates are unavailable. limma-voom is used for integer count data, and limma is used for normalized count data.
  - logFC
  - logMean (AMean)
### Example output directory
![Example output directory](docs/Screenshot_outputs_directory.png)
## Directory 
```
├── src/
|   ├── main_v2.py
│   ├── save_gse_meta.py
│   ├── Cleantitle_v3.py
│   ├── DownloadCounts_v2.py
│   ├── run_deseq2_v2.py
│   ├── run_limmaVoom_v2.py
│   ├── run_limma_v2.py                             
│
├── script/
│   ├── Downloadcount_v2_script.py                   
│
├── docs/
|   ├── Screenshot_outputs_directory.png
|
├── LICENSE
├── README.md                     
└── requirements.txt
```
## Notes
- This pipeline operation relies on count matrices of Gene IDs. It does not handle FASTQ files, read alignment, or quantification from raw sequencing data.
- As long as a GEO dataset provides Gene-level count data that can be parsed and aggregated at the GSM level, the pipeline is expected to work.
- Due to the heterogeneous and often noisy nature of public GEO data, users are encouraged to inspect intermediate results and pay attention to the console outputs for warning or error messages. 
