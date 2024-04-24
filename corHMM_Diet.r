library("corHMM")
library("mclust")
library("stringr")
library("qpcR")

# Loading data
phy<-read.nexus("16FC_16C_374_sp.tree")
df<-read.table("table_data_diet.tsv", header = TRUE)

df$Species<-str_replace((df$Species), " ", "_")

states<-cbind(c(str_replace((df$Species), " ", "_"), setdiff(phy$tip.label, str_replace((df$Species), " ", "_"))), c(df$Diet, rep("?", length(setdiff(phy$tip.label, str_replace((df$Species), " ", "_"))))))

states_traits<-as.data.frame(states[!states[,1] %in% setdiff(states[,1], phy$tip.label),])

colnames(states_traits)<-c("Species", "Diet")

states_traits<-states_traits[match(phy$tip.label,states_traits[,1]),]

states_traits<-states_traits[,c(1,2)]
LegendAndRateMat <- getStateMat4Dat(states_traits)
RateMat <- LegendAndRateMat$rate.mat
RateMat_trans <- dropStateMatPars(RateMat, c(1,3))
pars2equal <- list(c(1,3), c(2,4))
RateMat_trans_sym <- equateStateMatPars(RateMat_trans, pars2equal)
pars2equal <- list(c(1:4))
RateMat_trans_eq <- equateStateMatPars(RateMat_trans, pars2equal)

diet_1_eq<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_1_sym<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_1_ard<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_1_tran_I_M<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=RateMat_trans, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_1_tran_I_M_sym<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=RateMat_trans_sym, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_1_tran_I_M_eq<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=RateMat_trans_eq, node.states = "marginal",
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

diet_2_eq<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_2_sym<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_2_ard<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_2_tran_I_M<-corHMM(phy, states_traits, rate.cat = 2, rate.mat= RateMat_trans_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_2_tran_I_M_sym<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=RateMat_trans_sym_2_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_2_tran_I_M_eq<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=RateMat_trans_eq_2_rates, node.states = "marginal",
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

diet_3_eq<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_3_sym<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_3_ard<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_3_tran_I_M<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=RateMat_trans_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_3_tran_I_M_sym<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=RateMat_trans_sym_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

diet_3_tran_I_M_eq<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=RateMat_trans_eq_3_rates, node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

df_diet<-data.frame(cbind(c(diet_1_eq$loglik, diet_1_sym$loglik, diet_1_ard$loglik, diet_1_tran_I_M$loglik, diet_1_tran_I_M_sym$loglik, diet_1_tran_I_M_eq$loglik,
                          diet_2_eq$loglik, diet_2_sym$loglik, diet_2_ard$loglik, diet_2_tran_I_M$loglik, diet_2_tran_I_M_sym$loglik, diet_2_tran_I_M_eq$loglik,
                          diet_3_eq$loglik, diet_3_sym$loglik, diet_3_ard$loglik, diet_3_tran_I_M$loglik, diet_3_tran_I_M_sym$loglik, diet_3_tran_I_M_eq$loglik),
                        c(diet_1_eq$AICc, diet_1_sym$AICc, diet_1_ard$AICc, diet_1_tran_I_M$AICc, diet_1_tran_I_M_sym$AICc, diet_1_tran_I_M_eq$AICc,
                          diet_2_eq$AICc, diet_2_sym$AICc, diet_2_ard$AICc, diet_2_tran_I_M$AICc, diet_2_tran_I_M_sym$AICc, diet_2_tran_I_M_eq$AICc,
                          diet_3_eq$AICc, diet_3_sym$AICc, diet_3_ard$AICc, diet_3_tran_I_M$AICc, diet_3_tran_I_M_sym$AICc, diet_3_tran_I_M_eq$AICc),
                akaike.weights(c(diet_1_eq$AICc, diet_1_sym$AICc, diet_1_ard$AICc, diet_1_tran_I_M$AICc, diet_1_tran_I_M_sym$AICc, diet_1_tran_I_M_eq$AICc,
                          diet_2_eq$AICc, diet_2_sym$AICc, diet_2_ard$AICc, diet_2_tran_I_M$AICc, diet_2_tran_I_M_sym$AICc, diet_2_tran_I_M_eq$AICc,
                          diet_3_eq$AICc, diet_3_sym$AICc, diet_3_ard$AICc, diet_3_tran_I_M$AICc, diet_3_tran_I_M_sym$AICc, diet_3_tran_I_M_eq$AICc))$deltaAIC,
                akaike.weights(c(diet_1_eq$AICc, diet_1_sym$AICc, diet_1_ard$AICc, diet_1_tran_I_M$AICc, diet_1_tran_I_M_sym$AICc, diet_1_tran_I_M_eq$AICc,
                          diet_2_eq$AICc, diet_2_sym$AICc, diet_2_ard$AICc, diet_2_tran_I_M$AICc, diet_2_tran_I_M_sym$AICc, diet_2_tran_I_M_eq$AICc,
                          diet_3_eq$AICc, diet_3_sym$AICc, diet_3_ard$AICc, diet_3_tran_I_M$AICc, diet_3_tran_I_M_sym$AICc, diet_3_tran_I_M_eq$AICc))$weights,
                c((max(as.vector(diet_1_eq$index.mat)[!is.na(as.vector(diet_1_eq$index.mat))])), (max(as.vector(diet_1_sym$index.mat)[!is.na(as.vector(diet_1_sym$index.mat))])), (max(as.vector(diet_1_ard$index.mat)[!is.na(as.vector(diet_1_ard$index.mat))])), (max(as.vector(diet_1_tran_I_M$index.mat)[!is.na(as.vector(diet_1_tran_I_M$index.mat))])), (max(as.vector(diet_1_tran_I_M_sym$index.mat)[!is.na(as.vector(diet_1_tran_I_M_sym$index.mat))])), (max(as.vector(diet_1_tran_I_M_eq$index.mat)[!is.na(as.vector(diet_1_tran_I_M_eq$index.mat))])),
(max(as.vector(diet_2_eq$index.mat)[!is.na(as.vector(diet_2_eq$index.mat))])), (max(as.vector(diet_2_sym$index.mat)[!is.na(as.vector(diet_2_sym$index.mat))])), (max(as.vector(diet_2_ard$index.mat)[!is.na(as.vector(diet_2_ard$index.mat))])), (max(as.vector(diet_2_tran_I_M$index.mat)[!is.na(as.vector(diet_2_tran_I_M$index.mat))])), (max(as.vector(diet_2_tran_I_M_sym$index.mat)[!is.na(as.vector(diet_2_tran_I_M_sym$index.mat))])), (max(as.vector(diet_2_tran_I_M_eq$index.mat)[!is.na(as.vector(diet_2_tran_I_M_eq$index.mat))])),
(max(as.vector(diet_3_eq$index.mat)[!is.na(as.vector(diet_3_eq$index.mat))])), (max(as.vector(diet_3_sym$index.mat)[!is.na(as.vector(diet_3_sym$index.mat))])), (max(as.vector(diet_3_ard$index.mat)[!is.na(as.vector(diet_3_ard$index.mat))])), (max(as.vector(diet_3_tran_I_M$index.mat)[!is.na(as.vector(diet_3_tran_I_M$index.mat))])), (max(as.vector(diet_3_tran_I_M_sym$index.mat)[!is.na(as.vector(diet_3_tran_I_M_sym$index.mat))])), (max(as.vector(diet_3_tran_I_M_eq$index.mat)[!is.na(as.vector(diet_3_tran_I_M_eq$index.mat))])))
                ))
rownames(df_diet)<-c("diet_1_eq", "diet_1_sym", "diet_1_ard", "diet_1_tran_I_M", "diet_1_tran_I_M_sym", "diet_1_tran_I_M_eq",
                          "diet_2_eq", "diet_2_sym", "diet_2_ard", "diet_2_tran_I_M", "diet_2_tran_I_M_sym", "diet_2_tran_I_M_eq",
                          "diet_3_eq", "diet_3_sym", "diet_3_ard", "diet_3_tran_I_M", "diet_3_tran_I_M_sym", "diet_3_tran_I_M_eq")
colnames(df_diet)<-c("loglik", "AICc", "Delta_AICc", "AICcWt", "K_rates")

write.table(df_diet, "Results/df_diet.tsv", sep ="\t")
saveRDS(eval(parse(text = rownames(df_diet)[which.min(df_diet$AICc)])), "ASE/diet_corHMM.rds")
