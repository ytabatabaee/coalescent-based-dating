## Coalescent-based dating datasets

This repository contains the datasets and scripts used in the following paper:

- Y. Tabatabaee, S. Claramunt, S. Mirarab (2025). Coalescent-based branch length estimation improves dating of species trees. https://www.biorxiv.org/content/10.1101/2025.02.25.640207v1.abstract

For experiments in this study, we generated three sets of simulated datasets with gene tree discordance due to incomplete lineage sorting (ILS) and analyzed two avian biological datasets from [Harvey et. al. (2020)](https://www.science.org/doi/10.1126/science.aaz6970) and [Stiller et. al. (2024)](https://www.nature.com/articles/s41586-024-07323-1). The simulated datasets have model species trees with substitution-unit, generation-unit, and time-unit branch lengths. All datasets are available [here](https://drive.google.com/drive/folders/1fwU1Nc6BtzvXqQ1KYWPmkZqpC5w3c6_s?usp=sharing).

### Simulated datasets

**30-taxon dataset.**
This dataset has six model conditions with varying deviation from the molecular clock and inclusion of an outgroup, each with 100 replicates. The model conditions are specified as `outgroup.[has-OG].species.[DEV].genes.[DEV]` where `[has-OG]` is 1 when the dataset has an outgroup and 0 otherwise, and `[DEV]` shows the level of deviation from the clock (parameter α of the gamma distribution) that is set to 5 (low), 1.5 (medium), or 0.15 (high). Original dataset is from [Mai at al. (2017)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182238) and available at [https://uym2.github.io/MinVar-Rooting/](https://uym2.github.io/MinVar-Rooting/). Below is a description of files in each directory.

- `truegenetrees`: true gene trees
- `estimatedgenetre.gtr`: gene trees estimated under GTR evolution model
- `s_tree.trees`: true species tree in substitution units
- `s_tree_gu.trees`: true species tree in generation units
- `s_tree_tu_5.trees`: true species tree in time units (million years), assuming an average generation time of 5 years
- `s_tree_unit_ultrametric.trees`: unit ultrametric true species tree
- `s_tree_tu_5_calibrations_n[num-calib]_[root_unfixed].txt`: calibration information (given in the format of (node, time) pairs)
- `s_tree_tu_5_calib_n[num-calib]_[root_unfixed]_mcmctree.ctl`: control file for MCMCTree
- `s_tree_tu_5_calib_n[num-calib]_[root_unfixed]_mcmctree.date.nwk[.normalized]`: tree dated with MCMCTree with [num-calib] calibration points. The .[normalized] flag specifies the unit-ultrametric version of the dated tree.
- `RAxML_result.concat_align_s_tree.trees.rooted.labeled`: true species tree furnished with RAxML SU branch lengths
- `castlespro_estimatedgenetre.gtr_s_tree.trees.rooted.labeled`: true species tree furnished with CASTLES-Pro SU branch lengths
- `[dating-method]_n[num-calib]_[root_unfixed]_RAxML_result.concat_align_s_tree.trees.rooted.labeled.[normalized]`: RAxML SU tree dated with [dating-method] (can be treepl, wlogdate, mdcat, and lsd2) with [num-calib] calibration points. The .[normalized] flag specifies the unit-ultrametric version of the dated tree. Trees dated with lsd2 have a `.date.nwk` extension.
- `[dating-method]_n[num-calib]_[root_unfixed]_castlespro_estimatedgenetre.gtr_s_tree.trees.rooted.labeled.[normalized]`: CASTLES-Pro SU tree dated with [dating-method] (can be treepl, wlogdate, mdcat, and lsd2) with [num-calib] calibration points. The .[normalized] flag specifies the unit-ultrametric version of the dated tree. Trees dated with lsd2 have a `.date.nwk` extension.
- `ad.txt` : average RF distance between the model species tree and true gene trees
- `gtee_gtr.txt`: average RF distance between true and estimated gene trees


**101-taxon dataset.**
This dataset has four model conditions with varying sequence lengths (1600bp, 800bp, 400bp, 200bp) corresponding to different levels of gene tree estimation error (23%, 31%, 42%, and 55%). The original dataset is from [Zhang et. al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) and available at [https://gitlab.com/esayyari/ASTRALIII/](https://gitlab.com/esayyari/ASTRALIII/). Below is a description of files in each directory.

- `truegenetrees`: true gene trees
- `fasttree_genetrees_[seq-len]_non.[num-genes]`: [num-genes] gene trees estimated from sequence alignments with length [seq-len]bp
- `s_tree.trees`: true species tree in substitution units
- `s_tree_gu.trees`: true species tree in generation units
- `s_tree_tu_5.trees`: true species tree in time units (million years), assuming an average generation time of 5 years
- `s_tree_unit_ultrametric.trees`: unit ultrametric true species tree
- `s_tree_tu_5_calibrations_n[num-calib]_[root_unfixed].txt`: calibration information (given in the format of (node, time) pairs)
- `RAxML_result.concat_for_fasttree_[seq-len].[num-genes]_s_tree.trees.rooted.labeled`: true species tree furnished with RAxML SU branch lengths given a concatenation of [seq-len]bp sequences for the first [num-genes] genes
- `castlespro_fasttree_genetrees_[seq-len]_non.[num-genes]_s_tree.trees.rooted.labeled`: true species tree furnished with CASTLES-Pro SU branch lengths given [num-genes] FastTree gene trees estimated from [seq-len]bp sequences
- `ad.txt` : average RF distance between the model species tree and true gene trees
- `[dating-method]_n[num-calib]_[root_unfixed]_castlespro_fasttree_genetrees_[seq-len]_non.[num-genes]_s_tree.trees.rooted.labeled.[normalized]`: CASTLES-Pro SU tree dated with [dating-method] (can be treepl, wlogdate, mdcat, and lsd2) with [num-calib] calibration points for the model condition corresponding to [num-genes] genes and [seq-len]bp sequences. The .[normalized] flag specifies the unit-ultrametric version of the dated tree. Trees dated with lsd2 have a `.date.nwk` extension. Files with additional extensions include config files and log files for each method.
- `[dating-method]_n[num-calib]_[root_unfixed]_RAxML_result.concat_for_fasttree_[seq-len].[num-genes]_s_tree.trees.rooted.labeled.[normalized]`: RAxML SU tree dated with [dating-method] (can be treepl, wlogdate, mdcat, and lsd2) with [num-calib] calibration points for the model condition corresponding to [num-genes] genes and [seq-len]bp sequences. The .[normalized] flag specifies the unit-ultrametric version of the dated tree. Trees dated with lsd2 have a `.date.nwk` extension. Files with additional extensions include config files and log files for each method.

**Large dataset.**
This dataset has 8 model conditions with 50, 100, 200, 500, 1K, 2K, 5K, and 10K-taxon trees with 20 replicates in each condition. Below is a description of files in each directory `large/[num-taxa]/[rep-num]/` where `[num-taxa]` is the model condition (number of taxa) and `[rep-num]` is the replicate index.

- `truegenetrees`: true gene trees
- `estimatedgenetre`: estimated gene trees
- `s_tree.trees`: true species tree in substitution units
- `s_tree_gu.trees`: true species tree in generation units
- `s_tree_tu_5.trees`: true species tree in time units (million years), assuming an average generation time of 5 years
- `s_tree_unit_ultrametric.trees`: unit ultrametric true species tree
- `s_tree_tu_5_calibrations_n[num-calib].txt`: calibration information
- `RAxML_result.concat_s_tree.trees.rooted.labeled`: true species tree furnished with RAxML SU branch lengths
- `castlespro_estimatedgenetre_s_tree.trees.rooted.labeled`: true species tree furnished with CASTLES-Pro SU branch lengths
- `treepl_n[num-calib]_RAxML_result.concat_s_tree.trees.rooted.labeled`: RAxML SU tree dated with TreePL with [num-calib] calibration points
- `treepl_n[num-calib]_RAxML_result.concat_s_tree.trees.rooted.labeled.config`: config file for TreePL for the RAxML dated tree
- `treepl_n[num-calib]_castlespro_estimatedgenetre_s_tree.trees.rooted.labeled`: CASTLES-Pro SU tree dated with TreePL with [num-calib] calibration points
- `treepl_n[num-calib]_castlespro_estimatedgenetre_s_tree.trees.rooted.labeled.config`: config file for TreePL for the CASTLES-Pro dated tree
- `ad.txt` : average RF distance between the model species tree and true gene trees
- `gtee.txt` : average RF distance between true and estimated gene trees


### Biological datasets

- **Neoavian**: 363-taxon neoavian dataset from [Stiller et al. (2024)](https://www.nature.com/articles/s41586-024-07323-1) with 63,430 single-copy genes. The original data is available [here](https://sid.erda.dk/cgi-sid/ls.py?share_id=ENhZODU9YE). Results from the analysis in this study is available at [/biological/avian-stiller](https://github.com/ytabatabaee/coalescent-based-dating/tree/main/biological/avian-stiller). Below is a description of files in this directory.

- `mdcat_median_40Kl_caml_stiller.rooted.tre`: ASTRAL tree furnished with ConBL branch lengths dated with MD-Cat
- `mdcat_median_40Kl_astral4_stiller_nolabel.rooted.tre`: ASTRAL tree furnished with CASTLES-Pro branch lengths dated with MD-Cat
- `wlogdate_median_astral4_stiller.rooted.no_label.tre`: ASTRAL tree furnished with CASTLES-Pro branch lengths dated with wLogDate
- `wlogdate_median_astral_63K_concat_bl.rooted.no_label.tre`: ASTRAL tree furnished with ConBL branch lengths dated with wLogDate
- `treepl_median_astral_63K_concat_bl.rooted.tre`: ASTRAL tree furnished with ConBL branch lengths dated with TreePL
- `treepl_median_castlespro_stiller.rooted.tre`: ASTRAL tree furnished with CASTLES-Pro branch lengths dated with TreePL
- `genera_castlespro_concat.csv`: Age of genera estimated using different dating methods with CASTLES-Pro or ConBL branch lengths on concatenation or ASTRAL topology
- `families_castlespro_concat.csv`: Age of families estimated using different dating methods with CASTLES-Pro or ConBL branch lengths on concatenation or ASTRAL topology
- `orders_castlespro_concat.csv`: Age of orders estimated using different dating methods with CASTLES-Pro or ConBL branch lengths on concatenation or ASTRAL topology
- `ltt_stiller.csv`: Lineage-through-time information for different dated trees
- `fossil_list_median.txt`: Fossil calibration information for median (50\% quantiles)


- **Suboscines**: 1683-taxon suboscines dataset from [Harvey et. al. (2020)](https://www.science.org/doi/10.1126/science.aaz6970) with 2,389 single-copy genes. The original data is available at https://github.com/mgharvey/tyranni. Results from the analysis in this study is available at [/biological/suboscines-harvey](https://github.com/ytabatabaee/coalescent-based-dating/tree/main/biological/suboscines-harvey). Below is a description of files in this directory.

- `concat_T400F.examl.rooted.tre`: Concatenation tree furnished with CAML branch lengths
- `castlespro_T400F.examl.rooted.tre`: Concatenation tree furnished with CASTLES-Pro branch lengths
- `treepl_castlespro_T400F.astral.rooted.tre`: ASTRAL tree furnished with CASTLES-Pro branch lengths dated with TreePL
- `treepl_examl_T400F.astral.rooted.tre`: ASTRAL tree furnished with CAML branch lengths dated with TreePL
- `treepl_castlespro_T400F.examl.rooted.tre`: Concatenation tree furnished with CASTLES-Pro branch lengths dated with TreePL
- `treepl_concat_T400F.examl.rooted.tre`: Concatenation tree furnished with CAML branch lengths dated with TreePL
- `genera_treepl_castlespro_concat.csv`: Age of different genera estimated using TreePL+CASLTES-Pro andTreePL+ConBL on the concatenation topology
- `families_treepl_castlespro_concat.csv`: Age of different families estimated using TreePL+CASLTES-Pro andTreePL+ConBL on the concatenation topology
- `suboscines_ltt.csv`: Lineage-through-time information for the four different dated trees (ASTRAL or CAML furnished with CASTLES-Pro or ConBL)
