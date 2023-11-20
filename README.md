 # Phylogenetic conservatism and complex pattern of trait evolution in sharks: repository 

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Testing for phylogenetic signal](#1-Testing-for-phylogenetic-signal)
	- [1.1 Prerequisite](#11-Prerequisite)
	- [1.2 Discretization](#111-Discretization)
   	- [1.3 Discretization](#111-Discretization)
- [2 Diversification analysis](#2-Diversification-analysis)
	- [2.1 Selecting the likeliest model](#21-Selecting-the-likeliest-model)
	- [2.2 Bayesian analysis](#22-Bayesian-analysis)
	- [2.3 Ancestral state estimation](#23-Ancestral-state-estimation)
- [3 Sensitivity analyses](#3-Sensitivity-analyses)
	- [3.1 Revisiting MuSSE analysis](#31-Revisiting-MuSSE-analysis)
	- [3.2 Testing a single multi-state trait](#32-Testing-a-single-multi-state-trait)
	- [3.3 Testing multiple binary traits](#33-Testing-multiple-binary-traits)
- [Reference](#Reference)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "Phylogenetic conservatism and complex pattern of trait evolution in sharks" but can be generalized to other similar data as the scripts are not specific. All of the presented scripts are written in R language (R Core Team, 2022). If you are planning to use any of these scripts, please cite "Marion et al., 2024". </p>

## Overview

<p align="justify"> This repository contains scripts and markdown for performing the following comparative analyses:

**1**: Testing for phylogenetic signal

**2**: Modeling trait evolution

**3**: Conducting diversification analyses

<p align="justify"> Example data used for our article are available at "figshare.com".</p>

## 1 Testing for phylogenetic signal

`package requirement (mclust, geiger, stringr)`

`used script (Pagel's_Lambda.r)`

<p align="justify"> One of the most simple comparative analyses. Pagel's lambda is a metric whose purpose is to account for phylogenetic signals for a trait, namely testing whether traits are phylogenetically structured across the tree of life. This can be considered as a preliminary test when accounting for phylogenetic conservatism</p>

### 1.1 Prerequisite

<p align="justify"> The first step in any comparative analysis is to clean the data. Here, we want to test whether maximum body size is phylogenetically conserved across the phylogeny. To do so, we need to extract such values and isolate them in a vector. However maximum body size is a continuous data, thus we may need to discretize it. In this script, we use a Gaussian mixture model with the mclust function of the mclust package to cluster continuous data. This package selects automatically the most likely number of groups. Then,  we name each trait value according to its species name </p> 

### 1.2 Creating the white noise and Lambda differ models

<p align="justify"> The second step is to perform a test for phylogenetic signals using Pagel's lambda. To perform this test, we will first indicate which transition structure we want for our model. Here, we rely on the most exhaustive framework: all rates differ (ARD), which allow all possible transition. Secondly, we build two models, one with no phylogenetic signal (White noise) and one including phylogenetic signals (Lambda model). We then compare the relative fit of each model with AICC. </p> 

### 1.3 Making simulation

<p align="justify"> We also performed simulations using randomized datasets to test the significance of the Lambda value. The results is then plotted, with the barplot representing 1000 randomized replicates, and the redline the empirical Lambda value. </p> 

<p align="center">
    <img src="Lambda_bds.png" \>
</p>

## 2 Modeling trait evolution

`package requirement (corHMM, mclust, string, phytools, qpcR, ggtree, secsse)`

`used script (Diversity_analysis_all.r)`

<p align="justify"> The models used for trait-dependant diversification are known as "SSE" and originated from the original BiSSE model (Maddison et al., 2007). These are complex models using both trait data and branch length from a calibrated tree to correlate trait and diversification. Here we used the MuSSE model which was more fitted for our analysis. The second script is known as "Diversity_analysis.r" and allows the user to conduct a diversification analysis and ancestral state estimation for a multi-state trait under MUSSE. </p>

### 2.1 Selecting the likeliest model

<p align="justify"> The first step in conducting this analysis is to select the likeliest model. To do so, we build several variations of a complete model where all the rates (i.e. speciation, extinction, and transition between states) can vary. The best model is then selected based on statistical metrics such as AICc, delta AIC, and AIC weight. In our case, with the shark dataset, the best model is the full model with no variation in extinction rates.  </p>

### 2.2 Bayesian analysis

<p align="justify"> The second step consists in running the best maximum likelihood model into a bayesian framework. To do so, we need to reference priors, as they are mandatory for bayesian analyses. We decided to make an exponential prior, and a run of 10000 generations with 10 % burn-in. To see if convergence was achieved, we visualized the likelihood plot minus the burn-in period. Since the likelihood seems to have stabilized shortly after the burn-in, we considered that this run has converged. </p>

<p align="justify"> Secondly, we plotted the net diversification rates for each pair of trait and the pairwise net diversification rate difference between each trait. The first plot allows the user to visualize easily the range of net diversification rate for each trait. This plot is very qualitative and cannot account for real statistical differences, and hence, it is needed to characterize quantitatively the difference between rates.</p>


<p align="center">
    <img src="Figure_MuSSE_git.png" \>
</p>

<p align="justify"> The second plot accounts for such differences, as they are significant if the plot does not overlap 0 (represented by a red line). Here we can see two main groups, with respectively low and high net diversification rates. </p>

<p align="center">
    <img src="Differences-MuSSE-diversification-rates-MCMC-git.png" \>
</p>

### 2.3 Ancestral state estimation

<p align="justify"> The last step of this analysis is to compute ancestral states considering the phylogeny. Using bayesian data generated in part 2.2, we estimated each category's ancestral state with the help of the asr.marginal function. Computing A.S.E for each node of the phylogeny results in a probability table. Since the table is large, we decided to plot it directly over the phylogeny. Several ways to represent probability for each node in a phylogeny are available, we chose pie charts, as they are readable. </p>

## 3 Sensitivity analyses 

<p align="justify"> SSE models are sensitive to several mathematical biases (such as rejection of the null hypothesis). We carried out two additional analyses, to minimize such biases. </p>

#### 3.1 Revisiting MuSSE analysis

`package requirement (ape, secsse, DDD, tidyverse, parallel, qgraph, optional(stringr))`

`used script (SecSSE_Full.r)`

Accounting for trait effect is subject to numerous methodological biases (Beaulieu and Donoghue, 2013). Indeed, SSE models can falsely indicate an effect of the focal trait on diversification. Models with hidden traits, such as SeCSSE (Herrera-Alsina et al., 2019) or HiSSE (Beaulieu and O'Meara, 2016) can account for hidden variables in trait-dependant diversification. Here we use SeCSSE to detect : 

1 - an effect of the trait on diversification, 

2 - the possible existence of other hidden variables in our dataset (other variables than the trait-association variable).

3 - the coexistence of focal and hidden variables in our dataset (trait-association and hidden variables). 

Following Herrera-Alsina et al. (2019) and most importantly Liedtke et al. (2022), we revisited our MuSSE original results. The present script, as well as the "SecSSE_Reproduction.r" and "SecSSE_Body-size.r" scripts are adapted from Liedtke et al. (2022). As such, if you you want to reuse these scripts, please cite the original authors beforehand. Morevoer, some scripts from Liedtke et al. (2022) are required to run the analysis (aux_scripts).

#### 3.2 Testing a single multi-state trait

`package requirement (ape, secsse, DDD, tidyverse, parallel, qgraph, optional(stringr))`

`used script (SecSSE_Reproduction.r, SecSSE_Body-size.r)`

Individual trait analyses were conducted reuse the same archetype as for the trait-association variable.

<p align="justify">  </p>

#### 3.3 Testing multiple binary traits

`package requirement (diversitree, qpcR, optional(stringr))`

`used script (Multi_State_MUSSE.r)`

<p align="justify"> As habitat is hardly discriminable, the most efficient way to account for, this trait was to subdivide it into several binary traits, with one state being the presence of a species in a certain habitat, and the other its absence. Because of its nature, such a trait could not be analyzed with the previous method. Consequently, we used the function musse.multitrait presented in the package diversitree. Using the same dataset as before and with multitrait_binary_analysis.r script, you will be able to conduct this analysis.  </p>

### Reference

Beaulieu, J. M. Donoghue, M. J. (2013). Fruit evolution and diversification in campanulid angiosperms: Campanulid fruit evolution. Evolution. 67(11): 3132-3144.

Herrera-Alsina, L. van Els, P. Etienne, R. S. (2019). Detecting the dependence of diversification on multiple traits from phylogenetic trees and trait data. Systematic Biology. 68(2): 317-328.

Liedtke, H. C. Wiens, J. J. & Gomez-Mestre, I. (2022). The evolution of reproductive modes and life cycles in amphibians. Nature Communications, 13(1): 1-15.

Maddison, W.P., Midford, P.E. & Otto, S.P. (2007) Estimating a binary character's effect on speciation and extinction. Systematic Biology, 56, 701–710

Menardi, G. (2011). Density-based Silhouette diagnostics for clustering methods. Statistics and Computing. 21(3): 295-308.

Mérigot, B. Durbec, J. P. Gaertner, J. C. (2010). On goodness-of-fit measure for dendrogram-based analyses. Ecology. 91(6): 1850-1859.

R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.

Zambelli, A. E. (2016). A data-driven approach to estimating the number of clusters in hierarchical clustering. F1000 Research 5: 2809.
