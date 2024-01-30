 # Tempo and mode of trait evolution across the shark Tree of Life: repository 

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Modeling trait evolution](#1-Modeling-trait-evolution)
	- [1.1 Data cleaning and model construction](#11-Data-cleaning-and-model-construction)
	- [1.2 Selecting the likeliest model](#12-Selecting-the-likeliest-model)
	- [1.3 Ancestral state estimation](#13-Ancestral-state-estimation)
- [2 Diversification analyses](#2-Diversification-analyses)
	- [2.1 Data cleaning and model construction](#21-Data-cleaning-and-model-construction)
	- [2.2 Running and selecting the likeliest model](#22-Running-and-selecting-the-likeliest-model)
	- [2.3 Ancestral state estimation](#23-Ancestral-state-estimation)
- [3 Sensitivity analyses](#3-Sensitivity-analyses)
	- [3.1 Revisiting MuSSE analysis](#31-Revisiting-MuSSE-analysis)
	- [3.2 Testing a single multi-state trait](#32-Testing-a-single-multi-state-trait)
	- [3.3 Testing multiple binary traits](#33-Testing-multiple-binary-traits)
- [Reference](#Reference)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "Tempo and mode of trait evolution across the shark Tree of Life" but can be generalized to other similar data as the scripts are not specific. All of the presented scripts are written in R language (R Core Team, 2022). If you plan to use any of these scripts, please cite "Marion et al., 2024". </p>

## Overview

<p align="justify"> This repository contains scripts and html files for performing the following comparative analyses:

**1**: Modeling trait evolution

**2**: Analyses of diversification

**3**: Sensitivity analyses

<p align="justify"> Example data used for our article are available at "figshare.com".</p>

## 1 Modeling trait evolution

`package requirement (corHMM, mclust, string, phytools, qpcR, ggtree, secsse)`

`used script (corHMM_script.r/corHMM_ASE_script.r)`

<p align="justify">  The purpose of this analysis is twofold. First, we want to test how traits evolved and what evolutionary scenario they followed. Second, using corHMM output we will perform which can be later compared with the ancestral state estimation from SecSSE</p>

### 1.1 Data cleaning and model construction

<p align="justify"> The first step in any comparative analysis is to clean and prepare the trait data. To do so, we must extract and isolate such values in a vector. As maximum body size is continuous data, we may want to discretize it. In this script, we use a Gaussian mixture model with the mclust function of the mclust package to cluster continuous data. This package automatically selects the most likely number of groups. Then, we name each trait value according to its species name.

The second step is to build the actual evolutionary model we want to test. Usually, three models are tested when performing trait evolution analyses: an equal rate, a symmetric rate (transition from A to B is the same that B to A), and all rates differ, which is self-explanatory. However, such models may not be adequate when testing specific evolutionary scenarios. For example, we may want to test whether the direct transition from small size to large size is impossible. To do so we create a custom transition matrix by indicating that direct transitions between opposed states are impossible. In our case of Shark evolution, all these models will be referred to from now on as "sequential models", namely if one wants to go from A to C, there are two transitions, A to B, and B to C. Each model parameter specification was then entered in a standardised corHMM fashion. To properly model trait evolution, one should always keep in mind that evolutionary rates may vary across the phylogeny, thus it is best practice to introduce another parameter, the "rate" parameter. This parameter allows for heterogeneity of transition rates across the phylogeny, that is if one part of the phylogeny has short branching patterns as opposed to the rest of the tree, corHMM will take into account such heterogeneity, regardless of the trait state. Introducing such parameters is a great way to avoid false positives but is also very resource-consuming as these models are parameters rich (often more than twice the number of parameters for a one-rate model). To keep the likelihood space to allow and avoid over-parametrization, we did not go as far as 3 rates.

</p> 

### 1.2 Selecting the likeliest model

<p align="justify"> Then, we may perform each model and save relevant metrics for comparison. This script saves each model output into a data frame filled with the log-likelihood, the number of parameters, the AICc, the $\Delta$ AICc, and the $\omega$ AICc.
We then automatically save the output from the best-fitting model for each trait.</p> 

### 1.3 Ancestral state estimation

Lastly, one may be interested in estimating the ancestral condition of his clade of interest. Fortunately, corHMM jointly estimates transition rates and ancestral states, thus it is pretty easy to extract these values and plot directly on the phylogeny. To do this, we must use *corHMM_ASE_script.r*. While it is pretty easy to estimate them, ancestral state estimations are and remain **estimations**, thus one should always interpret them with utmost care.

<p align="center">
    <img src="ASE_BDS.png" \>
</p>

## 2 Diversification analyses

`package requirement (ape, secsse, DDD, tidyverse, parallel, qgraph, optional(stringr))`

`used script (SecSSE_bds.r/SecSSE_rp.r/SecSSE_ht.r)`

<p align="justify"> Now that we pictured trait evolution dynamics, we may want to assess whether our examined traits are responsible for extant diversity patterns or not. To do so we rely on SSE models (State-dependent speciation and extinction models), which are well-known for accounting for the impact that trait evolution has on patterns of lineage diversification. However, accounting for trait-dependent diversification is subject to numerous methodological biases (Beaulieu and Donoghue, 2013). Indeed, SSE models can falsely indicate an effect of the focal trait on diversification. Models with hidden traits (aka concealed traits), such as SeCSSE (Herrera-Alsina et al., 2019) or HiSSE (Beaulieu and O'Meara, 2016) can account for hidden variables in trait-dependant diversification. Thus, we will be using SecSSE, an SSE implementation for detecting trait-dependent diversification using phylogenies on multi-state characters. </p>

### 2.1 Data cleaning and model construction

<p align="justify"> Highly similar to corHMM, we will first need to discretise our continuous traits and clean our data. However, unlike corHMM, SecSSE estimates jointly transition rates and diversification rates (speciation and extinction). This is both a good thing, as it has been shown that accounting for diversification allows for better estimation of trait evolution (Maddison, 2006), and a bad thing as it dramatically increases the number of parameters in our analyses. Over-parametrization is a common issue in statistics, but it has not often been addressed in macroevolution. Here, we want to avoid as much as possible this issue by limiting the upper number of estimated parameters. To do so, we kept the structure of the best-fitting corHMM model for each trait, while still estimating its transition rates. Namely, if the best-fitting model for reproduction is a two-rate sequential symmetrical model, we will forbid direct transition from A to C, while constraining transitions from A to B to be equal to B to A. To minimize the effect of structuration on the output of the model, we kept the same structure of the transition matrix across all variants for each trait.

For each trait, we constructed seven models, one for constant rates (CR), three for examined-trait diversification (ETD) and three for concealed-trait diversification (CTD). The three variants for CTD and ETD models include a pure speciation model (different speciation for each state, one shared extinction), a pure extinction model (different extinction for each state, one shared speciation) and a net diversification model (different speciation and extinction for each state). For each of the seven models to o avoid finding a local optimum, following Herrera-Alsina et al. (2019), we used three sets of initial parameters. One set was estimated using the standard birth-death model (bd_ML function in the R package “DDD” 5.2.2; Etienne et al., 2012), one was halved, and the other was doubled.

</p>

### 2.2 Running and selecting the likeliest model

<p align="justify"> Then, we may perform each model and save relevant metrics for comparison. Here, we save the output of each of the 21 models for each trait in a separate directory. SecSSE does not compute the AICc directly, thus we used the following formula to estimate the AICc "((2*(-(ll)+2*k+(2*k*(k+1))/(n-k-1)),3)" where ll is the log-likelihood of the model, k its number of parameters and n the number of taxa included in the analysis. Then directly in the script, we construct a data frame filled with the log-likelihood, the number of parameters, the AICc, the $\Delta$ AICc, the $\omega$ AICc and the estimated rates. </p>


### 2.3 Ancestral state estimation

<p align="justify"> After the selection of the best-fitting model, one may want to estimate ancestral states across the phylogeny. However, unlike corHMM, SecSSE does not directly infer ancestral state, rather we may require an other procedure to properly estimate ancestral state. At the end of each script, there is the code for filling proper parameters. This procedure is a bit tricky, but you will find in the code that everything is annotated. Finally, after running our model (which is quite fast), we may plot our ancestral state directly on the phylogeny, similarly to what was done with corHMM  </p>

## 3 Sensitivity analyses 

<p align="justify"> SSE models are sensitive to several mathematical biases (such as rejection of the null hypothesis). We carried out two additional analyses, to minimize such biases. </p>

#### 3.1 Revisiting MuSSE analysis

`package requirement (ape, secsse, DDD, tidyverse, parallel, qgraph, optional(stringr))`

`used script (SecSSE_Full.r)`

Accounting for trait effect is subject to numerous methodological biases (Beaulieu and Donoghue, 2013). Indeed, SSE models can falsely indicate an effect of the focal trait on diversification. Models with hidden traits, such as SeCSSE (Herrera-Alsina et al., 2019) or HiSSE (Beaulieu and O'Meara, 2016) can account for hidden variables in trait-dependant diversification. Here we use SeCSSE to detect : 

1 - an effect of the trait on diversification, 

2 - the possible existence of other hidden variables in our dataset (other variables than the trait-association variable).

3 - the coexistence of focal and hidden variables in our dataset (trait association and hidden variables). 



#### 3.2 Comparing the 

`package requirement (ape, secsse, DDD, tidyverse, parallel, qgraph, optional(stringr))`

`used script (SecSSE_Reproduction.r, SecSSE_Body-size.r)`

Individual trait analyses were conducted reuse the same archetype as for the trait-association variable.

<p align="justify">  </p>

#### 3.3 Testing multiple binary traits

`package requirement (diversitree, qpcR, optional(stringr))`

`used script (Multi_State_MUSSE.r)`

<p align="justify"> As habitat is hardly discriminable, the most efficient way to account for, this trait was to subdivide it into several binary traits, with one state being the presence of a species in a certain habitat, and the other its absence. Because of its nature, such a trait could not be analyzed with the previous method. Consequently, we used the function musse.multitrait presented in the package diversitree. Using the same dataset as before and with multitrait_binary_analysis.r script, you will be able to conduct this analysis.  </p>


<p align="justify"> The second step consists in running the best maximum likelihood model into a bayesian framework. To do so, we need to reference priors, as they are mandatory for bayesian analyses. We decided to make an exponential prior, and a run of 10000 generations with 10 % burn-in. To see if convergence was achieved, we visualized the likelihood plot minus the burn-in period. Since the likelihood seems to have stabilized shortly after the burn-in, we considered that this run has converged. </p>

<p align="justify"> Secondly, we plotted the net diversification rates for each pair of trait and the pairwise net diversification rate difference between each trait. The first plot allows the user to visualize easily the range of net diversification rate for each trait. This plot is very qualitative and cannot account for real statistical differences, and hence, it is needed to characterize quantitatively the difference between rates.</p>

### Reference

Beaulieu, J. M. Donoghue, M. J. (2013). Fruit evolution and diversification in campanulid angiosperms: Campanulid fruit evolution. Evolution. 67(11): 3132-3144.

Herrera-Alsina, L. van Els, P. Etienne, R. S. (2019). Detecting the dependence of diversification on multiple traits from phylogenetic trees and trait data. Systematic Biology. 68(2): 317-328.

Maddison, W.P., Midford, P.E. & Otto, S.P. (2007) Estimating a binary character's effect on speciation and extinction. Systematic Biology, 56, 701–710

R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.
