# Multivariate analysis on traits

# I Data preparation

# Packages loading
library("mclust")
library("Rtsne")
library("ggplot2")
library("reshape2")
library("dplyr")
library("dendextend")
library("cluster")
library("fpc")

# Data loading

df1<-read.csv("Path/to/your/file.tsv", sep ="\t") # omit sep ="\t" for .csv files

# Clustering continous data

mb1 = Mclust(as.numeric(df1$"Your_variable"))
summary(mb1, parameters = TRUE)

plot(mb1, what=c("classification")) # plot the discretization

df1$"Your_variable"<-mb1$classification

df1<-na.omit(df1)

# Hierarchical clustering

df1<-(as.data.frame(unclass(df1),stringsAsFactors=TRUE))

df1$"Your_variable"<-as.factor(df1$"Your_variable")

rownames(df1)<-df1[,1]
df1<-df1[,-1]
df1.1<-df1

# Dissimilarity matrix creation

gower.dist <- daisy(as.data.frame(df1, metric = c("gower")))

divisive.clust <- diana(as.matrix(gower.dist),diss = TRUE, keep.diss = TRUE)
aggl.clust.c <- hclust(gower.dist, method = "complete")
aggl.clust.m <- hclust(gower.dist, method = "average")
aggl.clust.w <- hclust(gower.dist, method = "ward.D2")

# 2-norm and least square criterion

# Complete

sigma <- var(gower.dist)+var(cophenetic(aggl.clust.c))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(eigen(gower.dist-cophenetic(aggl.clust.c))$values)))
col.c<-c("Complete", thres, sign, sum((eigen(gower.dist-cophenetic(aggl.clust.c))$values)**2))

# UPGMA

sigma <- var(gower.dist)+var(cophenetic(aggl.clust.m))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(eigen(gower.dist-cophenetic(aggl.clust.m))$values)))
col.m<-c("UPGMA", thres, sign, sum((eigen(gower.dist-cophenetic(aggl.clust.m))$values)**2))

# Ward

sigma <- var(gower.dist)+var(cophenetic(aggl.clust.w))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(eigen(gower.dist-cophenetic(aggl.clust.w))$values)))
col.w<-c("Ward", thres, sign, sum((eigen(gower.dist-cophenetic(aggl.clust.w))$values)**2))

# Divisive

sigma <- var(gower.dist)+var(cophenetic(divisive.clust))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(eigen(gower.dist-cophenetic(divisive.clust))$values)))
col.div<-c("Divisive", thres, sign, sum((eigen(gower.dist-cophenetic(divisive.clust))$values)**2))

# Comparaison of algorithm dataframe

algo_sel <- rbind(col.c, col.m, col.w, col.div)
colnames(algo_sel)<-c("Name", "Threshold value", "Significance", "Least square")

algo_sel

# Select the best algorithm, here : UPGMA

# Following the method presented by Anastasia Reusova in her blogpost in Towardsdatascience

cstats.table <- function(dist, tree, k) {
clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                  "wb.ratio","dunn2","avg.silwidth")
clust.size <- c("cluster.size")
stats.names <- c()
row.clust <- c()
output.stats <- matrix(ncol = k, nrow = length(clust.assess))
cluster.sizes <- matrix(ncol = k, nrow = k)
    for(i in c(1:k)){
  row.clust[i] <- paste("Cluster-", i, " size")
}
    for(i in c(2:k)){
  stats.names[i] <- paste("Test", i-1)
  
  for(j in seq_along(clust.assess)){
    output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])[j]
    
  }
  
  for(d in 1:k) {
    cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
    dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
    cluster.sizes[d, i]
    
  }
}
    output.stats.df <- data.frame(output.stats)
    cluster.sizes <- data.frame(cluster.sizes)
cluster.sizes[is.na(cluster.sizes)] <- 0
    rows.all <- c(clust.assess, row.clust)
output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
colnames(output) <- stats.names[2:k]
rownames(output) <- rows.all
is.num <- sapply(output, is.numeric)
output[is.num] <- lapply(output[is.num], round, 2)
output
}

stats.df.agglm <-cstats.table(gower.dist, aggl.clust.m, 8) 
stats.df.agglm

# Elbow visualization

ggplot(data = data.frame(t(cstats.table(gower.dist, aggl.clust.m, 15))), 
  aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("UPGMA") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))

# Silhouette visualization

ggplot(data = data.frame(t(cstats.table(gower.dist, aggl.clust.m, 15))), 
  aes(x=cluster.number, y=avg.silwidth)) + 
  geom_point()+
  geom_line()+
  ggtitle("UPGMA") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))

# Dendrogramme visualization
dendro <- as.dendrogram(aggl.clust.m)
dendro.col <- dendro %>%
  set("branches_k_color", k = 5, value =   c("darkslategray", "darkslategray4", "darkslategray3", "gold3", "darkcyan")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram, k = 4")

# Traits & repartition for optimal group number

clust.num <- cutree(aggl.clust.m, k = 5)
df<-cbind(df1.1, clust.num)
df<-as.data.frame(unclass(df),stringsAsFactors=TRUE)
summary(df[df[,length(df[1,])]== 1,])
summary(df[df[,length(df[1,])]== 2,])
summary(df[df[,length(df[1,])]== 3,])
summary(df[df[,length(df[1,])]== 4,])
summary(df[df[,length(df[1,])]== 5,])

rownames(df)<-rownames(df1.1)

# Saving data
write.table(df, "Path/to/your/exit.file", sep ="\t")
