library("mclust")
library("Rtsne")
library("ggplot2")
library("reshape2")
library("dplyr")
library("dendextend")
library("cluster")
library("fpc")
library("FactoMineR")
library("factoextra")

df1<-read.csv("trait.tsv", sep ="\t") # omit sep ="\t" for .csv files

df1<-na.omit(df1)

df1<-df1 %>% arrange(bino)

mb1 = Mclust(as.numeric(df1$Consensus))
summary(mb1, parameters = TRUE)

plot(mb1, what=c("classification")) # plot the discretization

df1<-(as.data.frame(unclass(df1),stringsAsFactors=TRUE))

df1$"body_size"<-log(df1$Consensus)

rownames(df1)<-df1[,1]
df1<-df1[-c(1,2)]

df1.1.1<-df1

df1.1.1$"body_size"<-as.factor(mb1$classification)

res.mca<-MCA(df1.1.1, ncp = 5, graph = TRUE)

eig.val <- get_eigenvalue(res.mca)

fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))

fviz_mca_var(res.mca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # avoid text overlapping (slow)
             ggtheme = theme_minimal()
             )

fviz_mca_var(res.mca, axes = c(1,3), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # avoid text overlapping (slow)
             ggtheme = theme_minimal()
             )

fviz_mca_var(res.mca, axes = c(2,3), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # avoid text overlapping (slow)
             ggtheme = theme_minimal()
             )

gower.dist <- daisy(as.data.frame(df1), metric = c("gower"))

divisive.clust <- diana(as.matrix(gower.dist),diss = TRUE, keep.diss = TRUE)
aggl.clust.c <- hclust(gower.dist, method = "complete")
aggl.clust.m <- hclust(gower.dist, method = "average")
aggl.clust.w <- hclust(gower.dist, method = "ward.D2")

sigma <- var(gower.dist)+var(cophenetic(aggl.clust.c))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(svd(gower.dist-cophenetic(aggl.clust.c))$d)))
col.c<-c("Complete", thres, sign, sum((gower.dist-cophenetic(aggl.clust.c))**2))

sigma <- var(gower.dist)+var(cophenetic(aggl.clust.m))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(svd(gower.dist-cophenetic(aggl.clust.m))$d)))
col.m<-c("UPGMA", thres, sign, sum((gower.dist-cophenetic(aggl.clust.m))**2))

sigma <- var(gower.dist)+var(cophenetic(aggl.clust.w))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(svd(gower.dist-cophenetic(aggl.clust.w))$d)))
col.w<-c("Ward", thres, sign, sum((gower.dist-cophenetic(aggl.clust.w))**2))

sigma <- var(gower.dist)+var(cophenetic(divisive.clust))
thres <- 2*sqrt(nrow(as.matrix(gower.dist))*sigma)
sign<-(thres > max(abs(svd(gower.dist-cophenetic(divisive.clust))$d)))
col.div<-c("Divisive", thres, sign, sum((gower.dist-cophenetic(divisive.clust))**2))

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

ggplot(data = data.frame(t(cstats.table(gower.dist, aggl.clust.m, 15))), 
  aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("UPGMA") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = data.frame(t(cstats.table(gower.dist, aggl.clust.m, 15))), 
  aes(x=cluster.number, y=avg.silwidth)) + 
  geom_point()+
  geom_line()+
  ggtitle("UPGMA") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))

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

### Clearly indicates 5 groups

clust.num <- cutree(aggl.clust.m, k = 5)
df<-cbind(df1, clust.num)
df<-as.data.frame(unclass(df),stringsAsFactors=TRUE)
summary(df[df[,length(df[1,])]== 1,])
summary(df[df[,length(df[1,])]== 2,])
summary(df[df[,length(df[1,])]== 3,])
summary(df[df[,length(df[1,])]== 4,])
summary(df[df[,length(df[1,])]== 5,])

write.table(df, "Trait_syndrom_tab.tsv", sep ="\t")
