# README: DNA Microarray Data Analysis with Bioconductor and R

This repository contains the R code and data files for analyzing DNA microarray data using Bioconductor and R. The project is based on the example provided in Chapter 16 of the *Limma User Guide* and focuses on the "Swirl Zebrafish" experiment, which investigates differential gene expression between mutant and wild-type zebrafish.

## Project Overview

The primary goal of this project is to:
- Perform genome-wide analyses of differentially expressed genes.
- Use permutation distributions and False Discovery Rate (FDR) corrections.
- Generate visualizations to interpret the data and analysis results.

The repository includes:
- R Markdown file (`DNA_Microarray_Analysis.Rmd`) with detailed analysis steps and visualizations.
- Data files necessary to replicate the analysis.

## Key Features

1. **Data Preprocessing**:
   - Reads and inspects raw data files (e.g., targets and GAL files).
   - Creates image plots of background intensities and raw M-values.

2. **Normalization**:
   - Implements within-array normalization using "printtiploess."
   - Applies scale normalization across arrays as recommended by the *Limma User Guide*.

3. **Linear Modeling and Statistical Analysis**:
   - Fits a linear model to the normalized data.
   - Applies empirical Bayes moderation to improve statistical reliability.
   - Identifies differentially expressed genes.

4. **FDR Estimation**:
   - Compares results from Limma's FDR estimation with other methods using the `FDRestimation` package.

5. **Visualizations**:
   - MA-plots, boxplots, volcano plots, and image plots for data interpretation.
   - Effect of empirical Bayes moderation visualized through scatter plots.

## Getting Started

### Prerequisites

Ensure you have the following R packages installed:
- `limma`
- `knitr`
- `ggokabeito`
- `tidyverse`
- `FDRestimation`


### Folder Structure
```
repo/
├── data/                  # Contains raw data files (e.g., SwirlSample.txt, fish.gal)
├── DNA_Microarray_Analysis.Rmd  # R Markdown file for analysis
└── README.md              # Project documentation
```

### Running the Analysis
1. Open terminal.
2. Clone this repository: `git clone https://github.com/yourusername/your-repo-name.git`
3. `cd your-repo-name`
4. Open the DNA_Microarray_Analysis.Rmd file in RStudio.
5. Run the analysis by knitting the R Markdown file to the desired output format (HTML, PDF, Word).

## Results

The results include:
- Normalized gene expression data.
- Lists of significant genes with adjusted p-values.
- Visualizations to assess data quality and normalization effects.

## References

1. Smyth, G. K. (2005). Limma: Linear models for microarray data. In Bioinformatics and Computational Biology Solutions Using R and Bioconductor (pp. 397-420). Springer.
2. Limma User Guide.

## Acknowledgments

Special thanks to the Bioconductor community and the authors of the Limma User Guide for providing the foundation for this analysis.
