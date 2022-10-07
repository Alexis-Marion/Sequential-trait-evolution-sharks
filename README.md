# Phylogeny and evolution of Sharks : deciphering the role of traits in past diversification 

Phylogeny and evolution of Sharks: deciphering the role of traits in past diversification

This repository's purpose is to give a means of replicability to the article "placeholder" but can be generalized to other similar data as the scripts are not specific. All of the presented scripts are written in R language (source). You will also gain access to rdata files and notebook.

Overview

To analyze the role of traits in shark diversification, we had to conduct several trait-related diversification models. Those models can be cumbersome and, can account mostly for only a multi-state trait: testing several traits at once can be impossible. As such post-phylogeny analysis for diversification purposes is divided into three steps :

1: Merging trait data into a single trait-syndrome data
2: Conducting the actual diversification analysis
3: Conducting robustness analysis

Example data used for our article are available at "figshare.com".
For the first part, you can use the trait_data_0.tsv file.

1: Merging trait data 

package requirement (mclust, Rtsne, ggplot2, reshape2, dplyr, dendextend, cluster, fpc)

As explained earlier, today's models are not able to account for multiple traits, as such, testing their effect on diversification may need a workaround. To do so and by using statistical tools, we have created composite data summarizing all traits presented in the article with the first dedicated script called: "Multitrait_analysis".  

1.1 Discretization

The first, and optional, step in this analysis is to discretize continuous data. If you are using traits, you will probably handle continuous data, as such, you will need to discretize them. In this script, we use a Gaussian mixture model with the mclust function of the mclust package. This package selects automatically the most likely number of groups but beware as biologists we must consider only groups reflecting any biological reality. In this case, the body size trait is continuous and should be discretized.

1.2 Hierarchal clustering

Multivariate and clustering analysis can be quite complicated, especially for mixed and discrete data. Only a few options are available and here we choose to explore one of them: hierarchal clustering using Gower's distance. Gower's distance is appropriate for this kind of data and can be used through hierarchal clustering afterward. The problem is, hierarchal clustering can create groups, but the user may not know which number of groups is optimal. To avoid this and be as rigorous as we could select the best distance-based clustering algorithm and determine the optimal groups through different criteria.

1.2.1 Choosing the best algorithm

Distance based-clustering, as the phenetic method in phylogeny, uses several different reconstruction algorithms. Here we tried to choose the best algorithm based on two different criteria: the 2-norm criteria (sources) and the least-square criteria (sources). For these two criteria, the lower the better. The script allows the user to calculate these metrics, and here the best method is the UPGMA:

1.2.2 Post-analysis group determination 

Now that we possess the best algorithm for our dataset, we should find the best number of groups. Here we present two complementary methods to asses which is the optimal number of groups for our dataset.
The elbow method (sources) computes a score for each cluster determined with the algorithm, the higher the difference between successive clusters the better. As such, the user will probably select the optimal number of clusters when he sees a break in the curve. Similarly, the silhouette method is a measure of how similar a data point is within-cluster compared to other clusters. Here the higher the curve, the better.
Most of the time the two methods will come to similar if not identical results. Again the composition of each group should be carefully examined, as they normally should represent any biological reality. For our dataset, the best number of groups is five. 
As you can see here: 

2: Diversification analysis

package requirement ()

The models used for trait-dependant diversification are known as "SSE" models and originated all from the original BISSE model (source). These are complex models using both trait data and branch length from a calibrated tree to correlate trait and diversification. Here we used the MUSSE model which was more fitted for our analysis. The second script is known as "XXX" and allows the user to conduct a diversification analysis and ancestral state estimation for a multi-state trait under MUSSE.

2.1 Selecting the likeliest model



2.2 Bayesian analysis

2.3 Ancestral state estimation

3: Robustness analysis 

package requirement ()

SSE models are rather sensitive to several mathematical biases (such as rejection of the null hypothesis). To counteract this we carried out, for lack of better words, two robustness analyses.

3.1: True robustness

3.2: Testing each trait alone

Since we tested three traits at once, we wanted to measure the effect of each trait on diversification to account for their impact. Here we provide two ways to deal with individual trait data.




