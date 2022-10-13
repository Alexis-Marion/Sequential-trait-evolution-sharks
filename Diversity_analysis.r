# Diversity analysis

# I Data preparation

# Loading Packages
library("diversitree")
library("stringr")
library("qpcR")
library("ggtree")
library("ggplot2")
library("ggpmisc")

# Loading Packages
df<-read.csv("Path/to/your/file.tsv", sep="\t")
nex<-read.nexus("Path/to/your/tree.tree")

# Aicc function
aicc<-function(ll,k){round((2*(-round(ll,4))+2*k+(2*k*(k+1))/(Ntip(nex)-k-1)),3)}

# Optional step (cleaning data)
rownames(df)<-str_replace(rownames(df), " ", "_")
states<-as.numeric(df$clust.num)
names(states)<-rownames(df)
states2<-states[!names(states) %in% setdiff(names(states), nex$tip.label)]
nex3<-drop.tip(nex, setdiff(nex$tip.label, (names(states))))

# Computing sampling fractions
samp_frac<-c(
length(states2[states2==1])/length(df$clust.num[df$clust.num==1]),
length(states2[states2==2])/length(df$clust.num[df$clust.num==2]),
length(states2[states2==3])/length(df$clust.num[df$clust.num==3]),
length(states2[states2==4])/length(df$clust.num[df$clust.num==4]),
length(states2[states2==5])/length(df$clust.num[df$clust.num==5]))

# initialization
lik<-make.musse(nex3, states, 5, sampling.f = samp_frac, strict=TRUE)
p <- starting.point.musse(nex3, 5)

# Complete model
lik.all.free<-constrain(lik)
fit.all.free<-find.mle(lik.all.free, p[argnames(lik.all.free)])

p<-fit.all.free$par

# Null model
lik.base <- constrain(lik,
lambda5 ~ lambda1, lambda4 ~ lambda1, lambda3 ~ lambda1, lambda2 ~ lambda1,
mu5 ~ mu1, mu4 ~ mu1, mu3 ~ mu1, mu2 ~ mu1,
q15 ~ q12, q14 ~ q12, q13 ~ q12,
q25 ~ q12, q24 ~ q12, q23 ~ q12, q21 ~ q12,
q35 ~ q12, q34 ~ q12, q32 ~ q12, q31 ~ q12,
q45 ~ q12, q43 ~ q12, q42 ~ q12, q41 ~ q12,
q54 ~ q12, q53 ~ q12, q52 ~ q12, q51 ~ q12)

fit.base <- find.mle(lik.base, p[argnames(lik.base)])

# Free speciation model
lik.sp <- constrain(lik,
mu5 ~ mu1, mu4 ~ mu1, mu3 ~ mu1, mu2 ~ mu1,
q15 ~ q12, q14 ~ q12, q13 ~ q12,
q25 ~ q12, q24 ~ q12, q23 ~ q12, q21 ~ q12,
q35 ~ q12, q34 ~ q12, q32 ~ q12, q31 ~ q12,
q45 ~ q12, q43 ~ q12, q42 ~ q12, q41 ~ q12,
q54 ~ q12, q53 ~ q12, q52 ~ q12, q51 ~ q12)
fit.sp <- find.mle(lik.sp, p[argnames(lik.sp)])

# Free extinction model
lik.ext <- constrain(lik,
lambda5 ~ lambda1, lambda4 ~ lambda1, lambda3 ~ lambda1, lambda2 ~ lambda1,
q15 ~ q12, q14 ~ q12, q13 ~ q12,
q25 ~ q12, q24 ~ q12, q23 ~ q12, q21 ~ q12,
q35 ~ q12, q34 ~ q12, q32 ~ q12, q31 ~ q12,
q45 ~ q12, q43 ~ q12, q42 ~ q12, q41 ~ q12,
q54 ~ q12, q53 ~ q12, q52 ~ q12, q51 ~ q12)           
fit.ext <- find.mle(lik.ext, p[argnames(lik.ext)])

# Free transition model
lik.qt<-constrain(lik,                  
lambda5 ~ lambda1, lambda4 ~ lambda1, lambda3 ~ lambda1, lambda2 ~ lambda1,
mu5 ~ mu1, mu4 ~ mu1, mu3 ~ mu1, mu2 ~ mu1)
fit.qt <- find.mle(lik.qt, p[argnames(lik.qt)])

# Free speciation-extinction model
lik.sp_ext<-constrain(lik,
q15 ~ q12, q14 ~ q12, q13 ~ q12,
q25 ~ q12, q24 ~ q12, q23 ~ q12, q21 ~ q12,
q35 ~ q12, q34 ~ q12, q32 ~ q12, q31 ~ q12,
q45 ~ q12, q43 ~ q12, q42 ~ q12, q41 ~ q12,
q54 ~ q12, q53 ~ q12, q52 ~ q12, q51 ~ q12)
fit.sp_ext<- find.mle(lik.sp_ext, p[argnames(lik.sp_ext)])

# Free speciation-transition model
lik.sp_qt<-constrain(lik, 
mu5 ~ mu1, mu4 ~ mu1, mu3 ~ mu1, mu2 ~ mu1)
fit.sp_qt<- find.mle(lik.sp_qt, p[argnames(lik.sp_qt)])

# Free transition-extinction model
lik.qt_ext<-constrain(lik,
lambda5 ~ lambda1, lambda4 ~ lambda1, lambda3 ~ lambda1, lambda2 ~ lambda1)
fit.qt_ext<- find.mle(lik.qt_ext, p[argnames(lik.qt_ext)])

# Comparaison table
results<-matrix(NA,8,35) 


colnames(results)<-c("Nbparam","logL","AICc","DeltaAIC", "WAIC",argnames(lik))

rownames(results)<-c("base","sp","ext","qt","sp_qt","sp_ext","qt_ext","all.free")

for(j in 1:length(rownames(results)))
{
eval(parse(text=paste("results[j,names(fit.",rownames(results)[j],"$par)]<-fit.",rownames(results)[j],"$par",sep="")))
eval(parse(text=paste("results[j,'Nbparam']<-length(fit.",rownames(results)[j],"$par)",sep="")))
eval(parse(text=paste("results[j,'logL']<-fit.",rownames(results)[j],"$lnLik",sep="")))
eval(parse(text=paste("results[j,'AICc']<-aicc(fit.",rownames(results)[j],"$lnLik,length(fit.",rownames(results)[j],"$par) )",sep="")))
}

results[,4] <- round(akaike.weights(results[,3])$deltaAIC, 3)
results[,5] <- round(akaike.weights(results[,3])$weights, 3)
print(results)

# Saving data
write.csv(results, "Path/to/Result_tab.csv", sep ="\t")

# III Bayesian analysis

# Setting priors
p<-fit.sp_qt$par
prior <- make.prior.exponential(1/(2*(p[3]-p[8])))

# Running the analysis for 1000 generation and screening every 10 generation
set.seed(1)
tmp <- mcmc(lik.sp_qt, p, nsteps=100, prior=prior, w=1, print.every=10)

w <- diff(sapply(tmp[2:(length(tmp)-1)], quantile, c(0.025, 0.975)))

mcmc_MuSSE <- mcmc(lik.sp_qt, p, nsteps=1000, prior=prior, w=w)
save(mcmc_MuSSE, file="Path/to/Your_MCMC_DATA")

#Burn-in period
mcmc_MuSSE_new <- subset(mcmc_MuSSE, i > 100)

pdf("Path/to/Appendix_burnin.pdf")

par(mfrow=c(1,1), mar=c(3,4,1,1))

plot(mcmc_MuSSE_new$i, mcmc_MuSSE_new$p, xlim=c(0,10000), ty="l", xlab="Generations", ylab="log-likelihood", bty="n", main="After a 10% burn-in", cex.main="0.8")

dev.off()

# Setting the diversification difference
mcmc_MuSSEdiff <- with(mcmc_MuSSE_new, data.frame(S1 = lambda1, S2 = lambda2,  S3 = lambda3,   S4 = lambda4,   S5 = lambda5))

mcmc_MuSSEdiffS1 <- with(mcmc_MuSSE_new, data.frame(S1_2 = lambda1 - lambda2))
mcmc_MuSSEdiffS2 <- with(mcmc_MuSSE_new, data.frame(S1_3 = lambda1 - lambda3))
mcmc_MuSSEdiffS3 <- with(mcmc_MuSSE_new, data.frame(S2_3 = lambda2 - lambda3))
mcmc_MuSSEdiffS4 <- with(mcmc_MuSSE_new, data.frame(S1_4 = lambda1 - lambda4))
mcmc_MuSSEdiffS5 <- with(mcmc_MuSSE_new, data.frame(S1_5 = lambda1 - lambda5))
mcmc_MuSSEdiffS7 <- with(mcmc_MuSSE_new, data.frame(S2_4 = lambda2 - lambda4))
mcmc_MuSSEdiffS8 <- with(mcmc_MuSSE_new, data.frame(S1_5 = lambda2 - lambda5))
mcmc_MuSSEdiffS10 <- with(mcmc_MuSSE_new, data.frame(S3_4 = lambda3 - lambda4))
mcmc_MuSSEdiffS11 <- with(mcmc_MuSSE_new, data.frame(S3_5 = lambda3 - lambda5))
mcmc_MuSSEdiffS13 <- with(mcmc_MuSSE_new, data.frame(S4_5 = lambda4 - lambda5))

# Plot speciation 

pdf("Path/to/Figure-MuSSE.pdf")

par( mar=c(3,4,1,1))

colors=c("deepskyblue","grey", "yellow", "red", "green")

profiles.plot(mcmc_MuSSEdiff[1:5], col.line=colors, xlim=c(0,0.07), xlab="", ylab="", las=1, bty="n", main="a) Speciation rates", cex.main=0.8, n.br = 100)

legend("topright", bty="n", c("Shelf shark", "Oceanic shark", "Small deep shark", "Small reef shark", "Big reef shark"), col= colors, lty=1,lwd="4",cex=0.6)

dev.off()

# Plot speciation difference
pdf("Path/to/Appendix-Differences-MuSSE-speciation-rates-MCMC.pdf")

par(mfrow=c(5,3), mar=c(3,4,1,1))

colors=c("deepskyblue","grey", "yellow", "red", "green")

#diff r1 and r2
profiles.plot(mcmc_MuSSEdiffS1[1], col.line=c("grey"), xlim=c(-0.05,0.05), xlab="", ylab="", las=1, bty="n", main="a) Difference between blue and grey diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r1 and r3
profiles.plot(mcmc_MuSSEdiffS2[1], col.line=c("grey"), xlim=c(-0.15,0.05), xlab="", ylab="", las=1, bty="n", main="b) Difference between blue and yellow diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r1 and r4
profiles.plot(mcmc_MuSSEdiffS4[1], col.line=c("grey"), xlim=c(-0.05,0.05), xlab="", ylab="", las=1, bty="n", main="d) Difference between blue and red diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r1 and r5
profiles.plot(mcmc_MuSSEdiffS5[1], col.line=c("grey"), xlim=c(-0.15,0.05), xlab="", ylab="", las=1, bty="n", main="e) Difference between blue and green diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r2 and r3
profiles.plot(mcmc_MuSSEdiffS3[1], col.line=c("grey"), xlim=c(-0.15,0.05), xlab="", ylab="", las=1, bty="n", main="c) Difference between grey and yellow diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r2 and r4
profiles.plot(mcmc_MuSSEdiffS7[1], col.line=c("grey"), xlim=c(-0.05,0.05), xlab="", ylab="", las=1, bty="n", main="f) Difference between grey and red diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r2 and r5
profiles.plot(mcmc_MuSSEdiffS8[1], col.line=c("grey"), xlim=c(-0.15,0.05), xlab="", ylab="", las=1, bty="n", main="g) Difference between grey and green diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r3 and r4
profiles.plot(mcmc_MuSSEdiffS10[1], col.line=c("grey"), xlim=c(-0.05,0.05), xlab="", ylab="", las=1, bty="n", main="h) Difference between yellow and red diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r3 and r5
profiles.plot(mcmc_MuSSEdiffS11[1], col.line=c("grey"), xlim=c(-0.15,0.05), xlab="", ylab="", las=1, bty="n", main="i) Difference between yellow and green diversification", cex.main=0.8)
abline(v=c(0),col="red")

#diff r4 and r5
profiles.plot(mcmc_MuSSEdiffS13[1], col.line=c("grey"), xlim=c(-0.05,0.05), xlab="", ylab="", las=1, bty="n", main="j) Difference between red and green diversification", cex.main=0.8)
abline(v=c(0),col="red")

dev.off()

# IV Ancestral state estimmation

attach("Shark_MCMC_DATA")

# Applying the ancestral state estimmation for each trait
st1<-apply(mcmc_MuSSE_new[2:27], 1, function(x) asr.marginal(lik.sp_qt, x)[1,])
st2<-apply(mcmc_MuSSE_new[2:27], 1, function(x) asr.marginal(lik.sp_qt, x)[2,])
st3<-apply(mcmc_MuSSE_new[2:27], 1, function(x) asr.marginal(lik.sp_qt, x)[3,])
st4<-apply(mcmc_MuSSE_new[2:27], 1, function(x) asr.marginal(lik.sp_qt, x)[4,])
st5<-apply(mcmc_MuSSE_new[2:27], 1, function(x) asr.marginal(lik.sp_qt, x)[5,])

st.m.avg1<-rowMeans(st1)
st.m.avg2<-rowMeans(st2)
st.m.avg3<-rowMeans(st3)
st.m.avg4<-rowMeans(st4)
st.m.avg5<-rowMeans(st5)

table_ASE<-as.data.frame(t(rbind(st.m.avg1, st.m.avg2, st.m.avg3, st.m.avg4, st.m.avg5)))

# Saving the probability table
write.csv(as.data.frame(table_ASE), "Path/to/Table_MCMC_ASE.tsv", sep ="\t")

# V Plotting the Ancestral State Estimmation

# Generating the tree (semicircular)
phy_cil<-ggtree(nex3, layout="fan", open.angle=15)  + 
theme_bw() +
      theme(panel.border = element_blank(),
            legend.key = element_blank(),
           axis.ticks = element_blank(),
           axis.text.y = element_blank(),
           axis.text.x = element_blank(),
           panel.grid = element_blank(),
           panel.grid.minor = element_blank(), 
           panel.grid.major = element_blank(),
                   panel.background = element_blank(),
               plot.background = element_rect(fill = "transparent",colour = NA))

# Generating pie chart for each node
table_ASE$node<-c(376:750) # selecting internal nodes
pies <- nodepie(table_ASE, cols=1:5, color=c('#9E0142', '#F46D43', '#FEE08B', '#ABDDA4', '#5E4FA2'), alpha=1)

# Integrating pie chart into the phylogeny
df_ASE<-tibble::tibble(node=as.numeric(table_ASE$node), pies=pies)
phy_cil2 <- phy %<+% df_ASE
ASE_plot<-phy_cil2 + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.03, vp.height=0.03, hjust=0.5, vjust=0.5)

# Saving the tree
ggsave(ASE_plot, filename = "Path/to/output.pdf",  bg = "transparent", width = 10, height = 10)
