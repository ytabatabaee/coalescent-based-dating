## Coalescent-based dating datasets

This repository contains the datasets and scripts used in the following paper:

- Y. Tabatabaee, S. Claramunt, S. Mirarab (2024). Coalescent-based branch length estimation improves dating of species trees. https://www.biorxiv.org/content/10.1101/2025.02.25.640207v1.abstract

For experiments in this study, we generated three sets of simulated datasets with gene tree discordance due to incomplete lineage sorting (ILS) and analyzed two avian biological datasets from [Harvey et. al. (2020)](https://www.science.org/doi/10.1126/science.aaz6970) and [Stiller et. al. (2024)](https://www.nature.com/articles/s41586-024-07323-1). The simulated datasets have model species trees with substitution-unit, generation-unit, and time-unit branch lengths. All datasets are available [here](https://drive.google.com/drive/folders/1fwU1Nc6BtzvXqQ1KYWPmkZqpC5w3c6_s?usp=sharing).

### Simulated datasets

**30-taxon dataset**

- Original dataset is from [Mai at al. (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182238) and available at [https://uym2.github.io/MinVar-Rooting/](https://uym2.github.io/MinVar-Rooting/).
- Results and intermediate data from the experiments in the paper are available [here](https://drive.google.com/drive/folders/1fwU1Nc6BtzvXqQ1KYWPmkZqpC5w3c6_s?usp=sharing).

**101-taxon dataset**
- Original dataset is from [Zhang et. al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) and available at [https://gitlab.com/esayyari/ASTRALIII/](https://gitlab.com/esayyari/ASTRALIII/).
- Results and intermediate data from the experiments in the paper are available [here](https://drive.google.com/drive/folders/1fwU1Nc6BtzvXqQ1KYWPmkZqpC5w3c6_s?usp=sharing).

**Large dataset**
- This dataset has 8 model conditions with 50, 100, 200, 500, 1K, 2K, 5K, and 10K-taxon trees.
- Raw dataset and results and intermediate data from the experiments in the paper are available [here](https://drive.google.com/drive/folders/1fwU1Nc6BtzvXqQ1KYWPmkZqpC5w3c6_s?usp=sharing).

### Biological datasets

- **Neoavian**: 363-taxon neoavian dataset from [Stiller et al. (2024)](https://www.nature.com/articles/s41586-024-07323-1) with 63,430 single-copy genes. The original data is available [here](https://sid.erda.dk/cgi-sid/ls.py?share_id=ENhZODU9YE). Results from the analysis in this study is available at [/biological/avian-stiller](https://github.com/ytabatabaee/coalescent-based-dating/tree/main/biological/avian-stiller).
- **Suboscines**: 1683-taxon suboscines dataset from [Harvey et. al. (2020)](https://www.science.org/doi/10.1126/science.aaz6970) with 2,389 single-copy genes. The original data is available at https://github.com/mgharvey/tyranni. Results from the analysis in this study is available at [/biological/suboscines-harvey](https://github.com/ytabatabaee/coalescent-based-dating/tree/main/biological/suboscines-harvey).
