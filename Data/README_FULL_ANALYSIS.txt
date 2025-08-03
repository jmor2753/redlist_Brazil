
README – Data, Methods, and Outputs for Conservation Red List Analysis

This project accompanies the manuscript "Insects are underrepresented across red lists of threatened biodiversity in Brazil" and includes the complete data and R scripts necessary to replicate all analyses and visualizations.

---

## Overview

The study compares Brazil’s national (ICMBIO and MMA) and international (IUCN) red lists of threatened species with expectations based on the country's known biodiversity. It focuses on identifying biases in taxonomic representation, particularly for insect orders, by applying:

- Jaccard similarity analysis (species overlap across lists)
- Proportional comparison of expected vs. observed species
- Concordance Correlation Coefficient (CCC) for assessing alignment
- Faceted diagnostic scatterplots with custom labeling

---

## Input Data Files and Column Descriptions

### 1. `spps_dataset copy.csv`

Used to build presence/absence matrices and calculate Jaccard similarity.

| Column         | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| scientificName | Species binomial name                                                      |
| class          | Taxonomic class (e.g., Mammalia, Insecta)                                   |
| order          | Taxonomic order                                                             |
| family         | Taxonomic family                                                            |
| status         | Conservation status on the list                                             |
| list           | Source list (ICMBIO, MMA, or IUCN)                                          |

---

### 2. `scatter_plot_data copy.csv`

Used to compute expected vs. observed proportions and generate bar plots and scatterplots.

| Column         | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| class          | Taxonomic class of the species                                              |
| order          | Taxonomic order                                                             |
| list           | Conservation list source (ICMBIO, MMA, IUCN)                                |
| iucn_category  | IUCN threat level (CR, EN, VU, NT, LC, DD)                                   |
| exp            | Expected count of species (based on risk proportion)                        |
| obt            | Observed count of species on the list                                       |

---

### 3. `ccc_dataset copy.csv`

Used to plot the Concordance Correlation Coefficient (CCC) per class and list.

| Column            | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| class             | Taxonomic class                                                             |
| order             | Taxonomic order                                                             |
| list              | Conservation list source                                                    |
| regbr_class       | National registry classification (class level)                              |
| regbr_order       | National registry classification (order level)                              |
| iucn_category     | IUCN category under analysis                                                 |
| type              | Estimate type (e.g., CCC estimate, lower, or upper bound)                   |
| n_sp1             | Number of species used in the calculation                                   |
| list_prop_class   | Proportion of listed species within the class                               |

---

## Analytical Workflow

1. **Species Cleaning**: Duplicate rows removed to create a clean binary matrix.
2. **Similarity Index**: Jaccard distance calculated and transformed to similarity; bootstrapped CIs included.
3. **Proportions**: Computed expected vs. observed proportions of each order within classes per list.
4. **Diagnostic Scatterplots**:
   - Points colored and shaped by class and list.
   - Labels applied to orders where expected ≈ observed, and where expected ≫ observed.
   - Shapes: ICMBIO = circle, MMA = triangle, IUCN = square.
5. **Concordance Testing**:
   - Lin’s CCC computed to measure agreement between expected and observed species counts.
   - Bar charts show CCC values and 95% CIs per class and list.

---

## Outputs

- `similarity_plot`: Heatmap of Jaccard similarity with CI labels
- `prop_plot`: Bar plot of expected vs. observed order proportions
- `diagplot_labeled`: Faceted scatterplot with selected label colors
- `cccplot`: CCC bar chart with 95% confidence intervals

---

## Reproducibility

All scripts use R version 4.3.1 and the following packages:

- tidyverse
- vegan
- janitor
- ggh4x
- DescTools
- ggrepel
- RColorBrewer
- colorspace

---

## Contact

For questions about the data or analysis, contact:
Juliano Morimoto – juliano.morimoto@abdn.ac.uk

---

## Citation

If using this data or code, please cite the main manuscript:
Justino, G. L., Zanco, B., Guzman, L. M., et al. (2025). *Insects are underrepresented across red lists of threatened biodiversity in Brazil*.

