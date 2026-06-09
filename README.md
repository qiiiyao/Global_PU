Project for the article entitled 
"Global loss of phylogenetic uniqueness by naturalized species across taxonomic groups"
by
Qi Yao, Lirong Cai, Shao-peng Li, and Patrick Weigelt

Authorship of code: R Code was written by Qi Yao

Details of R script and their function stored in the code file:

1.data preparation/data preparation_ages1_35_equal_interval.R - Prepares the BSS dataset for modeling. Includes filtering for the 50 focal species, formatting species cover data into population growth data, determining the invasion stage of focal species and restricting the analysis to early successional years (ages 1–35).

2.data preparation/Phylo_Func/Calculate_funcd_pd.R - Calculates pairwise phylogenetic and functional distances among all species in the BSS community using phylogenetic trees and trait data.

3.fit/stancode - Contains Stan scripts for fitting six commonly used density-dependent population models, including variants of the Beverton–Holt model and Ricker models. These models were used to estimate interaction coefficients required for deriving ND and RFD.

4.fit/fit_plot_top50_ages1_35_equal_interval_all_models - Includes R scripts for model fitting and comparison across six candidate population models and different numbers of focal species using long-term species cover time series. Scripts were executed on a server due to computational demands (thousands of posterior samples). Resulting posteriors are not included in this archive due to file size.

5.fit/fit_plot_ages1_35_top50_mod_equal_interval_mod_comparison.R - Compares model performance for each focal species across all plots using Expected log pointwise predictive density (ELPD) from the loo R package, and additional filtering of best performing models based on Rhat convergence diagnostics.

6.results_analysing/results_analysing_ages1_35_top50_equal_interval_bh_partialb.R - Primary analysis pipeline. Evaluates how abundance-weighted phylogenetic/functional (PD/FD) and niche/relative fitness (ND/RFD) differences influence establishment and dominance of exotic species. Generates Figures 1, 2, S5, and S9 from the main text and supplement. Based on the best-performing model: the partial exponential Beverton–Holt (BH_partialb) model.

7.results_analysing/figures_code - Scripts to generate all remaining figures (both main and supplementary), except Figures 1, 2, S5, and S9. Includes visualization of SEM results, model comparisons, trait/phylogenetic distances, and robustness checks.
