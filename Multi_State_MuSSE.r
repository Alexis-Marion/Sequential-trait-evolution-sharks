# Multi binary state : habitat

## I Data preparation

### Loading packages

library("diversitree")
library("stringr")
library("qpcR")

### Loading data

my_df<-read.csv("Trait_syndrom_tab.tsv", sep="\t")

aicc<-function(ll,k){round((2*(-round(ll,4))+2*k+(2*k*(k+1))/(Ntip(phy2.0)-k-1)),3)}

rownames(my_df)<-str_replace(rownames(my_df), " ", "_")

phy<-read.nexus("382sp_16C_20FC.tree")

my_df

states<-data.frame(my_df$Lower.Slope.Deep, my_df$Upper.Slope, my_df$Shelf, my_df$Reef, my_df$Coastal, my_df$Oceano.Pelagic)
states[states == "No"] <- 0
states[states == "Yes"] <- 1
rownames(states)<-rownames(my_df)
names(states)<-rownames(df)
states2<-states[!rownames(states) %in% setdiff(rownames(states), phy$tip.label),]

### Computing sampling fraction

samp_frac<-c(length(states2[,1])/length(states[,1]))

### Computing states percentages

states2[,1]<-as.numeric(states2[,1])
states2[,2]<-as.numeric(states2[,2])
states2[,3]<-as.numeric(states2[,3])
states2[,4]<-as.numeric(states2[,4])
states2[,5]<-as.numeric(states2[,5])
states2[,6]<-as.numeric(states2[,6])

phy2.0<-drop.tip(phy, setdiff(phy$tip.label, (rownames(states))))

colnames(states2)<-c("A", "B", "C", "D", "E", "F")

## II Running the analysis

### Null model (no variation in any rates)

samp_frac

lik0 <- make.musse.multitrait(phy2.0, states2, sampling.f = samp_frac, depth = 0)
argnames(lik0)
p <- starting.point.musse.multitrait(phy2.0, lik0)
fit0 <- find.mle(lik0, p)

### Speciation 1 model (no combinaison)

lik_sp1 <- make.musse.multitrait(phy2.0, states2, sampling.f = samp_frac, depth = c(1, 0, 0))
argnames(lik_sp1)
p <- starting.point.musse.multitrait(phy2.0, lik_sp1)
fit_sp1 <- find.mle(lik_sp1, p)

### Speciation 2 model (combinaison of two traits allowed)

lik_sp2 <- make.musse.multitrait(phy2.0, states2, sampling.f = samp_frac, depth = c(2, 0, 0))
argnames(lik_sp2)
p <- starting.point.musse.multitrait(phy2.0, lik_sp2)
fit_sp2 <- find.mle(lik_sp2, p)

sort(fit_sp2$par[c(2:22)])

### Speciation 2 model (combinaison of three traits allowed)

lik_sp3 <- make.musse.multitrait(phy2.0, states2, sampling.f = samp_frac, depth = c(3, 0, 0))
argnames(lik_sp3)
p <- starting.point.musse.multitrait(phy2.0, lik_sp3)
fit_sp3 <- find.mle(lik_sp3, p)

### Comparing models

# Aicc function
aicc<-function(ll,k){round((2*(-round(ll,4))+2*k+(2*k*(k+1))/(Ntip(phy2.0)-k-1)),3)}

aicc_vec<-c(aicc(fit0$lnLik, length(argnames(lik0))), aicc(fit_sp1$lnLik, length(argnames(lik_sp1))), aicc(fit_sp2$lnLik, length(argnames(lik_sp2))),aicc(fit_sp3$lnLik, length(argnames(lik_sp3))))

habitat_binary_model_comp<-as.data.frame(cbind(c(fit0$lnLik, fit_sp1$lnLik, fit_sp2$lnLik, fit_sp3$lnLik), c(length(argnames(lik0)), length(argnames(lik_sp1)), length(argnames(lik_sp2)), length(argnames(lik_sp3))), aicc_vec, round(akaike.weights(aicc_vec)$deltaAIC, 3), round(akaike.weights(aicc_vec)$weights, 3)))

colnames(habitat_binary_model_comp)<-c("log likelihood", "number of parameters", "AICc", "Delta AIC", "AIC weights")
rownames(habitat_binary_model_comp)<-c("Null model", "Speciation 1", "Speciation 2", "Speciation 3")

habitat_binary_model_comp

write.table(habitat_binary_model_comp, "MuSSE_Multi_State.tsv", sep ="\t")
