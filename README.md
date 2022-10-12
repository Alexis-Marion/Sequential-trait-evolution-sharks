 # A dense time-calibrated phylogeny of sharks provides insights into the role of traits on their deep-time diversification : repository 

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Creating trait-syndrome data](#1-Creating-trait-syndrome-data)
	- [1.1 Discretization](#1.1-Discretization)
	- [1.2 Hierarchal clustering](#1.2-Hierarchal-clustering)
		- [1.2.1 Choosing the best algorithm](#1.2.1-Choosing-the-best-algorithm)
		- [1.2.2 Post-analysis group determination](#1.2.2-Post-analysis-group-determination)
- [2 Diversification analysis](#2-Diversification-analysis)
	- [2.1 Selecting the likeliest model](#2.1-Selecting-the-likeliest-model)
	- [2.2 Bayesian analysis](#2.2-Bayesian-analysis)
	- [2.3 Ancestral state estimation](#2.3-Ancestral-state-estimation)
- [3 Robustness analysis](#3-Robustness-analysis)
	- [3.1 Permutation analysis](#3.1-Permutation-analysis)
	- [3.2 Testing each trait alone](#3.2-Testing-each-trait-alone)
		- [3.2.1 Testing a single multi-state trait](#3.2.1-Testing-a-single-multi-state-trait)
		- [3.2.2 Testing multiple binary traits](#3.2.2-Testing-multiple-binary-traits)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "A dense time-calibrated phylogeny of sharks provides insights into the role of traits on their deep-time diversification" but can be generalized to other similar data as the scripts are not specific. All of the presented scripts are written in R language (R Core Team, 2022). You will also gain access to rdata files and notebook. If you are planning to use any of these scripts, please cite "XXX". </p>



## Overview

<p align="justify"> To analyze the role of traits in shark diversification, we had to conduct several trait-related diversification models. Those models can be cumbersome and, can account mostly for only a multi-state trait: testing several traits at once can be impossible. As such post-phylogeny analysis for diversification purposes is divided into three steps :

**1**: Merging trait data into a single trait-syndrome data

**2**: Conducting the actual diversification analysis

**3**: Conducting robustness analysis

<p align="justify"> Example data used for our article are available at "figshare.com".</p>

## 1 Creating trait-syndrome data

`package requirement (mclust, Rtsne, ggplot2, reshape2, dplyr, dendextend, cluster, fpc)`

`used script (Multitrait_analysis)`

<p align="justify"> As explained earlier, today's models are not able to account for multiple traits, as such, testing their effect on diversification may need a workaround. To do so and using statistical tools, we have created composite data summarizing all traits presented in the article with the first dedicated script called: "Multitrait_analysis.r".</p>

The following script and explanation are inspired by the excellent  ["Hierarchical Clustering on Categorical Data in R"](https://towardsdatascience.com/hierarchical-clustering-on-categorical-data-in-r-a27e578f2995), written by Anastasia Reusova.

### 1.1 Discretization

<p align="justify"> The first, and optional, step is to discretize continuous data. If you are using traits, you will probably handle continuous data, as such, you need to discretize them. In this script, we use a Gaussian mixture model with the mclust function of the mclust package. This package selects automatically the most likely number of groups. </p> 

### 1.2 Hierarchal clustering

<p align="justify"> Multivariate and clustering analysis can be quite complicated, especially for mixed and discrete data. Only a few options are available and here we choose to explore one of them: hierarchal clustering using Gower's distance. Gower's distance is appropriate for this kind of data and can be used through hierarchal clustering afterward. The problem is, hierarchal clustering can create groups, but the user may not know which number of groups is optimal. To avoid this and be as rigorous as we could select the best distance-based clustering algorithm and determine the optimal groups through different criteria. </p>

#### 1.2.1 Choosing the best algorithm

<p align="justify"> Distance based-clustering, as the phenetic method in phylogeny, uses several different reconstruction algorithms. Here we tried to choose the best algorithm based on two different criteria: the 2-norm criterion (Mérigot et al., 2010) and the least-square criterion. For these two criteria, the lower the better. The script allows the user to calculate these metrics, and here the best method is the UPGMA </p>

#### 1.2.2 Post-analysis group determination 

<p align="justify"> Now that we possess the best algorithm for our dataset, we should find the best number of groups. We present two complementary methods to asses which is the optimal number of groups for our dataset.
The elbow method (Zambelli, 2016) computes a score for each cluster determined with the algorithm, the higher the difference between successive clusters the better. As such, the user will probably select the optimal number of clusters when he sees a break in the curve. Similarly, the silhouette method (Menardi, 2011) is a measure of how similar a data point is within-cluster compared to other clusters. Here the higher the curve, the better. Most of the time the two methods will come to similar if not identical results. Again, the composition of each group should be carefully examined, as they should represent any biological reality. For our dataset, the best number of groups is five. </p>

## 2 Diversification analysis

`package requirement (diversitree, qpcR, ggtree, ggplot2, ggpmisc, optional(stringr))`

`used script (Diversity_analysis.r)`

<p align="justify"> The models used for trait-dependant diversification are known as "SSE" and originated from the original BISSE model (Maddison et al., 2007). These are complex models using both trait data and branch length from a calibrated tree to correlate trait and diversification. Here we used the MUSSE model which was more fitted for our analysis. The second script is known as "Diversity_analysis.r" and allows the user to conduct a diversification analysis and ancestral state estimation for a multi-state trait under MUSSE. </p>

### 2.1 Selecting the likeliest model

<p align="justify"> The first step in conducting this analysis is to select the likeliest model. To do so we build several variations of a complete model where all the rates (ie speciation, extinction, and transition between states) can vary. The best model is then selected based on statistical metrics such as AICc, AIC delta, and AIC weight. In our case, with the shark dataset, the best model is the full model with no variation in extinction rates.  </p>

### 2.2 Bayesian analysis

<p align="justify"> The second step consists in running the best maximum likelihood model into a bayesian framework. To do so we need to reference priors, as they are mandatory for bayesian analysis. We decided to make an exponential prior, and a run of 1000 generations with 10 % burn-in. To see if convergence was achieved, we visualized the likelihood plot minus the burn-in period. Since the likelihood seems to have stabilized shortly after the burn-in, we decided to consider that this run has converged. </p>

<p align="justify"> Secondly, we plotted the speciation rates for each trait and the speciation difference between each trait. The first plot allows the user to visualize easily the range of speciation rate for each trait. Here, the grey group possesses the smallest speciation rate among all the clusters. Conversely, the yellow group possesses the largest diversification rate, almost 3 times more than the grey group. This plot is very qualitative, but cannot account for real statistical differences, therefore, there is a need to characterize quantitatively the difference between rates. The second plot accounts for such differences, as they are significant if the plot is not adjacent to 0 (represented by a red line). Here we can see two main groups, with respectively low and large speciation rates. </p>

### 2.3 Ancestral state estimation

<p align="justify"> The last step of this analysis is to compute the ancestral state for the phylogeny. Using bayesian data generated in part 2.2, we estimated each category's ancestral state with the help of the asr.marginal function. Computing A.S.E for each node of the phylogeny results in a probability table. Since the table is large, we decided to plot it directly in the phylogeny. Several ways to represent probability for each node in a phylogeny are available, we chose pie charts, as they are readable. </p>

## 3 Robustness analysis 

<p align="justify"> SSE models are sensitive to several mathematical biases (such as rejection of the null hypothesis). To minimize such biases we carried out two additional analysis. </p>

### 3.1 Permutation analysis

`package requirement ()`

`used script (permutation.r)`

<p align="justify"> As a first way to examine for overestimating the rejection of the Null hypothesis, we performed a permutation analysis on the same dataset as previously mentioned.  </p>

### 3.2 Testing each trait alone

<p align="justify"> Since we tested three traits at once, we wanted to measure the effect of each trait on diversification to account for their impact. Here we provide two ways to deal with individual trait data. </p>

#### 3.2.1 Testing a single multi-state trait

`package requirement (secsse, qpcR, optional(stringr))`

`used script (SecSSE_Reproduction, SecSSE_Body-size)`

<p align="justify"> Accounting for trait effect is subject to numerous methodological biases (Beaulieu and Donoghue, 2013). Indeed, SSE models can falsely indicate an effect of the focal trait on diversification. Models with hidden traits, such as SECSSE (Herrera-Alsina et al., 2019) or HISSE (Beaulieu and O'Meara, 2016) can account for hidden variables in trait-dependant diversification. Here we use SECSSE to detect : 
1 - an effect of the trait on diversification, 
2 - the possible existence of other hidden variables in our dataset (the other variable of the trait-syndrome variable). </p>

#### 3.2.2 Testing multiple binary traits

*package requirement* (diversitree, qpcR, optional(stringr))

*used script* (Multi_State_MUSSE)

<p align="justify"> As the habitat is hardly discriminable, the most efficient way to account for it was to subdivide it into several binary traits, with one state being the presence of a species in a certain habitat, and the other its absence. Because of its nature, such traits could not be analyzed with the previous method. As such, we used the function musse.multitrait presented in the package diversitree. Using the same dataset as before and with multitrait_binary_analysis.r script, you will be able to conduct this analysis.  </p>

### Reference

R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.



Maddison, W.P., Midford, P.E. & Otto, S.P. (2007) Estimating a binary character's effect on speciation and extinction. Systematic Biology, 56, 701–710
