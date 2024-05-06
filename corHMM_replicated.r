library("corHMM")
library("mclust")
library("stringr")
library("phytools")
library("qpcR")
library("secsse")

# Loading data

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("Species_data.tsv", sep="\t") # omit sep ="\t" for .csv files

phy<-read.nexus(args[1])

phy<-phy[[as.numeric(args[2])]]

mb1 = Mclust(as.numeric(df$Body.size))
summary(mb, parameters = TRUE)
df$Body.size <- as.factor(mb1$classification)

df$Species<-str_replace(df$Species, " ", "_")

setdiff(phy$tip.label, df$Species)

states<-cbind(df$Species, df$Body.size, df$Reproduction, df$Habitat)

states_traits<-states[!states[,1] %in% setdiff(states[,1], phy$tip.label),]

states_traits[is.na(states_traits)]<-"?"

states_traits<-states_traits[match(phy$tip.label,states_traits[,1]),]

states_size<-states_traits[,c(1,2)]
LegendAndRateMat <- getStateMat4Dat(states_size)
RateMat <- LegendAndRateMat$rate.mat
RateMat_trans <- dropStateMatPars(RateMat, c(2,5))
pars2equal <- list(c(1,2), c(3,4))
RateMat_trans_sym <- equateStateMatPars(RateMat_trans, pars2equal)
pars2equal <- list(c(1:4))
RateMat_trans_eq <- equateStateMatPars(RateMat_trans, pars2equal)

bds_1_eq<-corHMM(phy, states_size, rate.cat = 1, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_1_sym<-corHMM(phy, states_size, rate.cat = 1, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_1_ard<-corHMM(phy, states_size, rate.cat = 1, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_1_tran_s_l<-corHMM(phy, states_size, rate.cat = 1, rate.mat=RateMat_trans, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_1_tran_s_l_sym<-corHMM(phy, states_size, rate.cat = 1, rate.mat=RateMat_trans_sym, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_1_tran_s_l_eq<-corHMM(phy, states_size, rate.cat = 1, rate.mat=RateMat_trans_eq, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

RateMat_trans_2_rates<-list(RateMat_trans, RateMat_trans)
RateMat_trans_sym_2_rates<-list(RateMat_trans_sym, RateMat_trans_sym)
RateMat_trans_eq_2_rates<-list(RateMat_trans_eq, RateMat_trans_eq)
RateClassMat <- getRateCatMat(2) 
RateMat_trans_2_rates <- getFullMat(RateMat_trans_2_rates, RateClassMat)
RateMat_trans_sym_2_rates <- getFullMat(RateMat_trans_sym_2_rates, RateClassMat)
RateMat_trans_eq_2_rates <- getFullMat(RateMat_trans_eq_2_rates, RateClassMat)

bds_2_eq<-corHMM(phy, states_size, rate.cat = 2, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_2_sym<-corHMM(phy, states_size, rate.cat = 2, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_2_ard<-corHMM(phy, states_size, rate.cat = 2, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_2_tran_s_l<-corHMM(phy, states_size, rate.cat = 2, rate.mat= RateMat_trans_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_2_tran_s_l_sym<-corHMM(phy, states_size, rate.cat = 2, rate.mat=RateMat_trans_sym_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_2_tran_s_l_eq<-corHMM(phy, states_size, rate.cat = 2, rate.mat=RateMat_trans_eq_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

RateMat_trans_3_rates<-list(RateMat_trans, RateMat_trans, RateMat_trans)
RateMat_trans_sym_3_rates<-list(RateMat_trans_sym, RateMat_trans_sym, RateMat_trans_sym)
RateMat_trans_eq_3_rates<-list(RateMat_trans_eq, RateMat_trans_eq, RateMat_trans_eq)
RateClassMat <- getRateCatMat(3) 
RateMat_trans_3_rates <- getFullMat(RateMat_trans_3_rates, RateClassMat)
RateMat_trans_sym_3_rates <- getFullMat(RateMat_trans_sym_3_rates, RateClassMat)
RateMat_trans_eq_3_rates <- getFullMat(RateMat_trans_eq_3_rates, RateClassMat)

bds_3_eq<-corHMM(phy, states_size, rate.cat = 3, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_3_sym<-corHMM(phy, states_size, rate.cat = 3, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_3_ard<-corHMM(phy, states_size, rate.cat = 3, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_3_tran_s_l<-corHMM(phy, states_size, rate.cat = 3, rate.mat=RateMat_trans_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_3_tran_s_l_sym<-corHMM(phy, states_size, rate.cat = 3, rate.mat=RateMat_trans_sym_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bds_3_tran_s_l_eq<-corHMM(phy, states_size, rate.cat = 3, rate.mat=RateMat_trans_eq_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

df_size<-data.frame(cbind(c(bds_1_eq$loglik, bds_1_sym$loglik, bds_1_ard$loglik, bds_1_tran_s_l$loglik, bds_1_tran_s_l_sym$loglik, bds_1_tran_s_l_eq$loglik,
                          bds_2_eq$loglik, bds_2_sym$loglik, bds_2_ard$loglik, bds_2_tran_s_l$loglik, bds_2_tran_s_l_sym$loglik, bds_2_tran_s_l_eq$loglik,
                          bds_3_eq$loglik, bds_3_sym$loglik, bds_3_ard$loglik, bds_3_tran_s_l$loglik, bds_3_tran_s_l_sym$loglik, bds_3_tran_s_l_eq$loglik),
                        c(bds_1_eq$AICc, bds_1_sym$AICc, bds_1_ard$AICc, bds_1_tran_s_l$AICc, bds_1_tran_s_l_sym$AICc, bds_1_tran_s_l_eq$AICc,
                          bds_2_eq$AICc, bds_2_sym$AICc, bds_2_ard$AICc, bds_2_tran_s_l$AICc, bds_2_tran_s_l_sym$AICc, bds_2_tran_s_l_eq$AICc,
                          bds_3_eq$AICc, bds_3_sym$AICc, bds_3_ard$AICc, bds_3_tran_s_l$AICc, bds_3_tran_s_l_sym$AICc, bds_3_tran_s_l_eq$AICc),
                akaike.weights(c(bds_1_eq$AICc, bds_1_sym$AICc, bds_1_ard$AICc, bds_1_tran_s_l$AICc, bds_1_tran_s_l_sym$AICc, bds_1_tran_s_l_eq$AICc,
                          bds_2_eq$AICc, bds_2_sym$AICc, bds_2_ard$AICc, bds_2_tran_s_l$AICc, bds_2_tran_s_l_sym$AICc, bds_2_tran_s_l_eq$AICc,
                          bds_3_eq$AICc, bds_3_sym$AICc, bds_3_ard$AICc, bds_3_tran_s_l$AICc, bds_3_tran_s_l_sym$AICc, bds_3_tran_s_l_eq$AICc))$deltaAIC,
                akaike.weights(c(bds_1_eq$AICc, bds_1_sym$AICc, bds_1_ard$AICc, bds_1_tran_s_l$AICc, bds_1_tran_s_l_sym$AICc, bds_1_tran_s_l_eq$AICc,
                          bds_2_eq$AICc, bds_2_sym$AICc, bds_2_ard$AICc, bds_2_tran_s_l$AICc, bds_2_tran_s_l_sym$AICc, bds_2_tran_s_l_eq$AICc,
                          bds_3_eq$AICc, bds_3_sym$AICc, bds_3_ard$AICc, bds_3_tran_s_l$AICc, bds_3_tran_s_l_sym$AICc, bds_3_tran_s_l_eq$AICc))$weights,
                c((max(as.vector(bds_1_eq$index.mat)[!is.na(as.vector(bds_1_eq$index.mat))])), (max(as.vector(bds_1_sym$index.mat)[!is.na(as.vector(bds_1_sym$index.mat))])), (max(as.vector(bds_1_ard$index.mat)[!is.na(as.vector(bds_1_ard$index.mat))])), (max(as.vector(bds_1_tran_s_l$index.mat)[!is.na(as.vector(bds_1_tran_s_l$index.mat))])), (max(as.vector(bds_1_tran_s_l_sym$index.mat)[!is.na(as.vector(bds_1_tran_s_l_sym$index.mat))])), (max(as.vector(bds_1_tran_s_l_eq$index.mat)[!is.na(as.vector(bds_1_tran_s_l_eq$index.mat))])),
(max(as.vector(bds_2_eq$index.mat)[!is.na(as.vector(bds_2_eq$index.mat))])), (max(as.vector(bds_2_sym$index.mat)[!is.na(as.vector(bds_2_sym$index.mat))])), (max(as.vector(bds_2_ard$index.mat)[!is.na(as.vector(bds_2_ard$index.mat))])), (max(as.vector(bds_2_tran_s_l$index.mat)[!is.na(as.vector(bds_2_tran_s_l$index.mat))])), (max(as.vector(bds_2_tran_s_l_sym$index.mat)[!is.na(as.vector(bds_2_tran_s_l_sym$index.mat))])), (max(as.vector(bds_2_tran_s_l_eq$index.mat)[!is.na(as.vector(bds_2_tran_s_l_eq$index.mat))])),
(max(as.vector(bds_3_eq$index.mat)[!is.na(as.vector(bds_3_eq$index.mat))])), (max(as.vector(bds_3_sym$index.mat)[!is.na(as.vector(bds_3_sym$index.mat))])), (max(as.vector(bds_3_ard$index.mat)[!is.na(as.vector(bds_3_ard$index.mat))])), (max(as.vector(bds_3_tran_s_l$index.mat)[!is.na(as.vector(bds_3_tran_s_l$index.mat))])), (max(as.vector(bds_3_tran_s_l_sym$index.mat)[!is.na(as.vector(bds_3_tran_s_l_sym$index.mat))])), (max(as.vector(bds_3_tran_s_l_eq$index.mat)[!is.na(as.vector(bds_3_tran_s_l_eq$index.mat))])))
                ))
rownames(df_size)<-c("bds_1_eq", "bds_1_sym", "bds_1_ard", "bds_1_tran_s_l", "bds_1_tran_s_l_sym", "bds_1_tran_s_l_eq",
                          "bds_2_eq", "bds_2_sym", "bds_2_ard", "bds_2_tran_s_l", "bds_2_tran_s_l_sym", "bds_2_tran_s_l_eq",
                          "bds_3_eq", "bds_3_sym", "bds_3_ard", "bds_3_tran_s_l", "bds_3_tran_s_l_sym", "bds_3_tran_s_l_eq")
colnames(df_size)<-c("loglik", "AICc", "Delta_AICc", "AICcWt", "K_rates")

write.table(df_size, paste("Results/", "df_size", args[2], ".tsv", sep =""), sep ="\t")

states_reproduction<-states_traits[,c(1,3)]
LegendAndRateMat <- getStateMat4Dat(states_reproduction)
RateMat <- LegendAndRateMat$rate.mat
RateMat_trans_o <- dropStateMatPars(RateMat, c(1,2,4,7))
pars2equal <- list(c(1,6), c(2,4), c(3,7), c(5,8))
RateMat_trans_o_sym <- equateStateMatPars(RateMat_trans_o, pars2equal)
pars2equal <- list(c(1:8))
RateMat_trans_o_eq <- equateStateMatPars(RateMat_trans_o, pars2equal)

rp_1_eq<-corHMM(phy, states_reproduction, rate.cat = 1, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_1_sym<-corHMM(phy, states_reproduction, rate.cat = 1, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_1_ard<-corHMM(phy, states_reproduction, rate.cat = 1, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_1_trans_o<-corHMM(phy, states_reproduction, rate.cat = 1, rate.mat= RateMat_trans_o, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_1_trans_o_sym<-corHMM(phy, states_reproduction, rate.cat = 1, rate.mat= RateMat_trans_o_sym, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_1_trans_o_eq<-corHMM(phy, states_reproduction, rate.cat = 1, rate.mat= RateMat_trans_o_eq, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

RateMat_trans_o_2_rates<-list(RateMat_trans_o, RateMat_trans_o)
RateMat_trans_o_sym_2_rates<-list(RateMat_trans_o_sym, RateMat_trans_o_sym)
RateMat_trans_o_eq_2_rates<-list(RateMat_trans_o_eq, RateMat_trans_o_eq)
RateClassMat <- getRateCatMat(2) 
RateMat_trans_o_2_rates <- getFullMat(RateMat_trans_o_2_rates, RateClassMat)
RateMat_trans_o_sym_2_rates <- getFullMat(RateMat_trans_o_sym_2_rates, RateClassMat)
RateMat_trans_o_eq_2_rates <- getFullMat(RateMat_trans_o_eq_2_rates, RateClassMat)

rp_2_eq<-corHMM(phy, states_reproduction, rate.cat = 2, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_2_sym<-corHMM(phy, states_reproduction, rate.cat = 2, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_2_ard<-corHMM(phy, states_reproduction, rate.cat = 2, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_2_trans_o<-corHMM(phy, states_reproduction, rate.cat = 2, rate.mat= RateMat_trans_o_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_2_trans_o_sym<-corHMM(phy, states_reproduction, rate.cat = 2, rate.mat= RateMat_trans_o_sym_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_2_trans_o_eq<-corHMM(phy, states_reproduction, rate.cat = 2, rate.mat= RateMat_trans_o_eq_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

RateMat_trans_o_3_rates<-list(RateMat_trans_o, RateMat_trans_o, RateMat_trans_o)
RateMat_trans_o_sym_3_rates<-list(RateMat_trans_o_sym, RateMat_trans_o_sym, RateMat_trans_o_sym)
RateMat_trans_o_eq_3_rates<-list(RateMat_trans_o_eq, RateMat_trans_o_eq, RateMat_trans_o_eq)
RateClassMat <- getRateCatMat(3) 
RateMat_trans_o_3_rates <- getFullMat(RateMat_trans_o_3_rates, RateClassMat)
RateMat_trans_o_sym_3_rates <- getFullMat(RateMat_trans_o_sym_3_rates, RateClassMat)
RateMat_trans_o_eq_3_rates <- getFullMat(RateMat_trans_o_eq_3_rates, RateClassMat)

rp_3_eq<-corHMM(phy, states_reproduction, rate.cat = 3, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_3_sym<-corHMM(phy, states_reproduction, rate.cat = 3, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_3_ard<-corHMM(phy, states_reproduction, rate.cat = 3, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_3_trans_o<-corHMM(phy, states_reproduction, rate.cat = 3, rate.mat= RateMat_trans_o_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_3_trans_o_sym<-corHMM(phy, states_reproduction, rate.cat = 3, rate.mat= RateMat_trans_o_sym_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

rp_3_trans_o_eq<-corHMM(phy, states_reproduction, rate.cat = 3, rate.mat= RateMat_trans_o_eq_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

df_reproduction<-data.frame(cbind(c(rp_1_eq$loglik, rp_1_sym$loglik, rp_1_ard$loglik, rp_1_trans_o$loglik, rp_1_trans_o_sym$loglik, rp_1_trans_o_eq$loglik,
                          rp_2_eq$loglik, rp_2_sym$loglik, rp_2_ard$loglik, rp_2_trans_o$loglik, rp_2_trans_o_sym$loglik, rp_2_trans_o_eq$loglik,
                          rp_3_eq$loglik, rp_3_sym$loglik, rp_3_ard$loglik, rp_3_trans_o$loglik, rp_3_trans_o_sym$loglik, rp_3_trans_o_eq$loglik),
                        c(rp_1_eq$AICc, rp_1_sym$AICc, rp_1_ard$AICc, rp_1_trans_o$AICc, rp_1_trans_o_sym$AICc, rp_1_trans_o_eq$AICc,
                          rp_2_eq$AICc, rp_2_sym$AICc, rp_2_ard$AICc, rp_2_trans_o$AICc, rp_2_trans_o_sym$AICc, rp_2_trans_o_eq$AICc,
                          rp_3_eq$AICc, rp_3_sym$AICc, rp_3_ard$AICc, rp_3_trans_o$AICc, rp_3_trans_o_sym$AICc, rp_3_trans_o_eq$AICc),
                akaike.weights(c(rp_1_eq$AICc, rp_1_sym$AICc, rp_1_ard$AICc, rp_1_trans_o$AICc, rp_1_trans_o_sym$AICc, rp_1_trans_o_eq$AICc,
                          rp_2_eq$AICc, rp_2_sym$AICc, rp_2_ard$AICc, rp_2_trans_o$AICc, rp_2_trans_o_sym$AICc, rp_2_trans_o_eq$AICc,
                          rp_3_eq$AICc, rp_3_sym$AICc, rp_3_ard$AICc, rp_3_trans_o$AICc, rp_3_trans_o_sym$AICc, rp_3_trans_o_eq$AICc))$deltaAIC,
                akaike.weights(c(rp_1_eq$AICc, rp_1_sym$AICc, rp_1_ard$AICc, rp_1_trans_o$AICc, rp_1_trans_o_sym$AICc, rp_1_trans_o_eq$AICc,
                          rp_2_eq$AICc, rp_2_sym$AICc, rp_2_ard$AICc, rp_2_trans_o$AICc, rp_2_trans_o_sym$AICc, rp_2_trans_o_eq$AICc,
                          rp_3_eq$AICc, rp_3_sym$AICc, rp_3_ard$AICc, rp_3_trans_o$AICc, rp_3_trans_o_sym$AICc, rp_3_trans_o_eq$AICc))$weights,
                c((max(as.vector(rp_1_eq$index.mat)[!is.na(as.vector(rp_1_eq$index.mat))])), (max(as.vector(rp_1_sym$index.mat)[!is.na(as.vector(rp_1_sym$index.mat))])), (max(as.vector(rp_1_ard$index.mat)[!is.na(as.vector(rp_1_ard$index.mat))])), (max(as.vector(rp_1_trans_o$index.mat)[!is.na(as.vector(rp_1_trans_o$index.mat))])), (max(as.vector(rp_1_trans_o_sym$index.mat)[!is.na(as.vector(rp_1_trans_o_sym$index.mat))])), (max(as.vector(rp_1_trans_o_eq$index.mat)[!is.na(as.vector(rp_1_trans_o_eq$index.mat))])),
(max(as.vector(rp_2_eq$index.mat)[!is.na(as.vector(rp_2_eq$index.mat))])), (max(as.vector(rp_2_sym$index.mat)[!is.na(as.vector(rp_2_sym$index.mat))])), (max(as.vector(rp_2_ard$index.mat)[!is.na(as.vector(rp_2_ard$index.mat))])), (max(as.vector(rp_2_trans_o$index.mat)[!is.na(as.vector(rp_2_trans_o$index.mat))])), (max(as.vector(rp_2_trans_o_sym$index.mat)[!is.na(as.vector(rp_2_trans_o_sym$index.mat))])), (max(as.vector(rp_2_trans_o_eq$index.mat)[!is.na(as.vector(rp_2_trans_o_eq$index.mat))])),
(max(as.vector(rp_3_eq$index.mat)[!is.na(as.vector(rp_3_eq$index.mat))])), (max(as.vector(rp_3_sym$index.mat)[!is.na(as.vector(rp_3_sym$index.mat))])), (max(as.vector(rp_3_ard$index.mat)[!is.na(as.vector(rp_3_ard$index.mat))])), (max(as.vector(rp_3_trans_o$index.mat)[!is.na(as.vector(rp_3_trans_o$index.mat))])), (max(as.vector(rp_3_trans_o_sym$index.mat)[!is.na(as.vector(rp_3_trans_o_sym$index.mat))])), (max(as.vector(rp_3_trans_o_eq$index.mat)[!is.na(as.vector(rp_3_trans_o_eq$index.mat))])))
                ))
rownames(df_reproduction)<-c("rp_1_eq", "rp_1_sym", "rp_1_ard", "rp_1_trans_o", "rp_1_trans_o_sym", "rp_1_trans_o_eq",
                          "rp_2_eq", "rp_2_sym", "rp_2_ard", "rp_2_trans_o", "rp_2_trans_o_sym", "rp_2_trans_o_eq",
                          "rp_3_eq", "rp_3_sym", "rp_3_ard", "rp_3_trans_o", "rp_3_trans_o_sym", "rp_3_trans_o_eq")
colnames(df_reproduction)<-c("logLik","AICc", "Delta_AICc", "AICcWt", "K_rates")

write.table(df_size, paste("Results/", "df_reproduction", args[2], ".tsv", sep =""), sep ="\t")

states_habitat<-states_traits[,c(1,4)]
LegendAndRateMat <- getStateMat4Dat(states_habitat)
RateMat <- LegendAndRateMat$rate.mat
RateMat_trans_d <- dropStateMatPars(RateMat, c(5,9,17,1,2,4))
pars2equal <- list(c(1,8), c(2,5), c(3,9), c(4,12), c(6,10), c(13,7), c(11,14))
RateMat_trans_d_sym <- equateStateMatPars(RateMat_trans_d, pars2equal)
pars2equal <- list(c(1:14))
RateMat_trans_d_eq <- equateStateMatPars(RateMat_trans_d, pars2equal)

ht_1_eq<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_1_sym<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_1_ard<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_1_d<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_1_d_sym<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_sym, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_1_d_eq<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_eq, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

RateMat_trans_d_2_rates<-list(RateMat_trans_d, RateMat_trans_d)
RateMat_trans_d_sym_2_rates<-list(RateMat_trans_d_sym, RateMat_trans_d_sym)
RateMat_trans_d_eq_2_rates<-list(RateMat_trans_d_eq, RateMat_trans_d_eq)
RateClassMat <- getRateCatMat(2) 
RateMat_trans_d_2_rates <- getFullMat(RateMat_trans_d_2_rates, RateClassMat)
RateMat_trans_d_sym_2_rates <- getFullMat(RateMat_trans_d_sym_2_rates, RateClassMat)
RateMat_trans_d_eq_2_rates <- getFullMat(RateMat_trans_d_eq_2_rates, RateClassMat)

ht_2_eq<-corHMM(phy, states_habitat, rate.cat = 2, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_2_sym<-corHMM(phy, states_habitat, rate.cat = 2, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_2_ard<-corHMM(phy, states_habitat, rate.cat = 2, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_2_d<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_2_d_sym<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_sym_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_2_d_eq<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_eq_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

RateMat_trans_d_3_rates<-list(RateMat_trans_d, RateMat_trans_d, RateMat_trans_d)
RateMat_trans_d_sym_3_rates<-list(RateMat_trans_d_sym, RateMat_trans_d_sym, RateMat_trans_d_sym)
RateMat_trans_d_eq_3_rates<-list(RateMat_trans_d_eq, RateMat_trans_d_eq, RateMat_trans_d_eq)
RateClassMat <- getRateCatMat(3) 
RateMat_trans_d_3_rates <- getFullMat(RateMat_trans_d_3_rates, RateClassMat)
RateMat_trans_d_sym_3_rates <- getFullMat(RateMat_trans_d_sym_3_rates, RateClassMat)
RateMat_trans_d_eq_3_rates <- getFullMat(RateMat_trans_d_eq_3_rates, RateClassMat)

ht_3_eq<-corHMM(phy, states_habitat, rate.cat = 3, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_3_sym<-corHMM(phy, states_habitat, rate.cat = 3, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_3_ard<-corHMM(phy, states_habitat, rate.cat = 3, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_3_d<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_3_d_sym<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_sym_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

ht_3_d_eq<-corHMM(phy, states_habitat, rate.cat = 1, rate.mat=RateMat_trans_d_eq_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

df_habitat<-data.frame(cbind(c(ht_1_eq$loglik, ht_1_sym$loglik, ht_1_ard$loglik, ht_1_d$loglik, ht_1_d_sym$loglik, ht_1_d_eq$loglik,
                          ht_2_eq$loglik, ht_2_sym$loglik, ht_2_ard$loglik, ht_2_d$loglik, ht_2_d_sym$loglik, ht_2_d_eq$loglik,
                          ht_3_eq$loglik, ht_3_sym$loglik, ht_3_ard$loglik, ht_3_d$loglik, ht_3_d_sym$loglik, ht_3_d_eq$loglik),
                        c(ht_1_eq$AICc, ht_1_sym$AICc, ht_1_ard$AICc, ht_1_d$AICc, ht_1_d_sym$AICc, ht_1_d_eq$AICc,
                          ht_2_eq$AICc, ht_2_sym$AICc, ht_2_ard$AICc, ht_2_d$AICc, ht_2_d_sym$AICc, ht_2_d_eq$AICc,
                          ht_3_eq$AICc, ht_3_sym$AICc, ht_3_ard$AICc, ht_3_d$AICc, ht_3_d_sym$AICc, ht_3_d_eq$AICc),
                akaike.weights(c(ht_1_eq$AICc, ht_1_sym$AICc, ht_1_ard$AICc, ht_1_d$AICc, ht_1_d_sym$AICc, ht_1_d_eq$AICc,
                          ht_2_eq$AICc, ht_2_sym$AICc, ht_2_ard$AICc, ht_2_d$AICc, ht_2_d_sym$AICc, ht_2_d_eq$AICc,
                          ht_3_eq$AICc, ht_3_sym$AICc, ht_3_ard$AICc, ht_3_d$AICc, ht_3_d_sym$AICc, ht_3_d_eq$AICc))$deltaAIC,
                akaike.weights(c(ht_1_eq$AICc, ht_1_sym$AICc, ht_1_ard$AICc, ht_1_d$AICc, ht_1_d_sym$AICc, ht_1_d_eq$AICc,
                          ht_2_eq$AICc, ht_2_sym$AICc, ht_2_ard$AICc, ht_2_d$AICc, ht_2_d_sym$AICc, ht_2_d_eq$AICc,
                          ht_3_eq$AICc, ht_3_sym$AICc, ht_3_ard$AICc, ht_3_d$AICc, ht_3_d_sym$AICc, ht_3_d_eq$AICc))$weights,
                c((max(as.vector(ht_1_eq$index.mat)[!is.na(as.vector(ht_1_eq$index.mat))])), (max(as.vector(ht_1_sym$index.mat)[!is.na(as.vector(ht_1_sym$index.mat))])), (max(as.vector(ht_1_ard$index.mat)[!is.na(as.vector(ht_1_ard$index.mat))])), (max(as.vector(ht_1_d$index.mat)[!is.na(as.vector(ht_1_d$index.mat))])), (max(as.vector(ht_1_d_sym$index.mat)[!is.na(as.vector(ht_1_d_sym$index.mat))])), (max(as.vector(ht_1_d_eq$index.mat)[!is.na(as.vector(ht_1_d_eq$index.mat))])),
(max(as.vector(ht_2_eq$index.mat)[!is.na(as.vector(ht_2_eq$index.mat))])), (max(as.vector(ht_2_sym$index.mat)[!is.na(as.vector(ht_2_sym$index.mat))])), (max(as.vector(ht_2_ard$index.mat)[!is.na(as.vector(ht_2_ard$index.mat))])), (max(as.vector(ht_2_d$index.mat)[!is.na(as.vector(ht_2_d$index.mat))])), (max(as.vector(ht_2_d_sym$index.mat)[!is.na(as.vector(ht_2_d_sym$index.mat))])), (max(as.vector(ht_2_d_eq$index.mat)[!is.na(as.vector(ht_2_d_eq$index.mat))])),
(max(as.vector(ht_3_eq$index.mat)[!is.na(as.vector(ht_3_eq$index.mat))])), (max(as.vector(ht_3_sym$index.mat)[!is.na(as.vector(ht_3_sym$index.mat))])), (max(as.vector(ht_3_ard$index.mat)[!is.na(as.vector(ht_3_ard$index.mat))])), (max(as.vector(ht_3_d$index.mat)[!is.na(as.vector(ht_3_d$index.mat))])), (max(as.vector(ht_3_d_sym$index.mat)[!is.na(as.vector(ht_3_d_sym$index.mat))])), (max(as.vector(ht_3_d_eq$index.mat)[!is.na(as.vector(ht_3_d_eq$index.mat))])))
                ))
rownames(df_habitat)<-c("ht_1_eq", "ht_1_sym", "ht_1_ard", "ht_1_d", "ht_1_d_sym", "ht_1_d_eq",
                          "ht_2_eq", "ht_2_sym", "ht_2_ard", "ht_2_d", "ht_2_d_sym", "ht_2_d_eq",
                          "ht_3_eq", "ht_3_sym", "ht_3_ard", "ht_3_d", "ht_3_d_sym", "ht_3_d_eq")
colnames(df_habitat)<-c("logLik","AICc", "Delta_AICc", "AICcWt", "K_rates")

write.table(df_size, paste("Results/", "df_habitat", args[2], ".tsv", sep =""), sep ="\t")
