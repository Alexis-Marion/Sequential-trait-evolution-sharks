# Loading packages
library("Rtsne")
library("ape")
library("secsse")
library("DDD")
library("doMC")
library("tidyverse")
library("qgraph")
library("stringr")

# Loading data

phy<-read.nexus("16FC_16C_374_sp.tree")
df<-read.csv("table_diet.tsv", sep="\t")

# Optional step (cleaning data)
df$Species<-str_replace((df$Species), " ", "_")
states<-cbind(c(str_replace((df$Species), " ", "_"), setdiff(phy$tip.label, str_replace((df$Species), " ", "_"))), c(df$Diet, rep(NA, length(setdiff(phy$tip.label, str_replace((df$Species), " ", "_"))))))
states_traits<-as.data.frame(states[!states[,1] %in% setdiff(states[,1], phy$tip.label),])
colnames(states_traits)<-c("Species", "Diet")
states<- states_traits$Diet
names(states)= states_traits$Species

states2<-states

# Computing sampling fractions
f<-c(
(length(states2[states2=="Invertebrate-feeders"])-length(states2[is.na(states2)]))/length(df$Diet[df$Diet=="Invertebrate-feeders"]),
(length(states2[states2=="Mesopredator"])-length(states2[is.na(states2)]))/length(df$Diet[df$Diet=="Mesopredator"]),
(length(states2[states2=="Macropredator"])-length(states2[is.na(states2)]))/length(df$Diet[df$Diet=="Macropredator"])
)

states3<-as.data.frame(cbind(names(states2), as.factor(states2)))
colnames(states3)<-c("species","states")
rownames(states3)<-NULL
traits <- sortingtraits(states3, phy)
traits<-traits-1

write.table(cbind(phy$tip.label,traits), "table_data_diet.tsv")

BD_model <- bd_ML(brts = ape::branching.times(phy))
BD_lambda <- BD_model$lambda0
BD_mu <- BD_model$mu0

strt_lambda<-c(BD_lambda, BD_lambda/2, BD_lambda*2)
strt_mu <- c(BD_mu, BD_mu/2, BD_mu*2)
strt_q <- strt_lambda/3

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 1))
spec_matrix <- rbind(spec_matrix, c(2, 2, 2, 1))

lambda_list_cr <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "CR")
lambda_list_cr

mu_vec_cr <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "CR",
                                   lambda_list = lambda_list_cr)
mu_vec_cr

shift_matrix_cr  <- c()
shift_matrix_cr  <- rbind(shift_matrix_cr , c(0, 2, 3))
shift_matrix_cr  <- rbind(shift_matrix_cr , c(2, 0, 4))
shift_matrix_cr  <- rbind(shift_matrix_cr , c(1, 2, 5))
shift_matrix_cr  <- rbind(shift_matrix_cr , c(2, 1, 6))

q_matrix_cr <- secsse::create_q_matrix(state_names = c(0, 1, 2),
                                    num_concealed_states = 3,
                                    shift_matrix = shift_matrix_cr,
                                    diff.conceal = TRUE)
q_matrix_cr[q_matrix_cr>=7]<-7
q_matrix_cr

### Try 1

idparsopt <- c(1:7) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(strt_lambda[1], strt_mu[1], rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_cr
idparslist[[3]] <- q_matrix_cr

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 7
saveRDS(model, "SecSSE_Results/Diet_results/diet_model_try1.rds")

### Try 2

idparsopt <- c(1:7) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(strt_lambda[2], strt_mu[2], rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_cr
idparslist[[3]] <- q_matrix_cr

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 7
saveRDS(model, "SecSSE_Results/Diet_results/diet_model_try2.rds")

### Try 3

idparsopt <- c(1:7) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(strt_lambda[3], strt_mu[3], rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_cr
idparslist[[3]] <- q_matrix_cr

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 7
saveRDS(model, "SecSSE_Results/Diet_results/diet_model_try3.rds")

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
spec_matrix <- rbind(spec_matrix, c(2, 2, 2, 3))

lambda_list_etd <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "ETD")
lambda_list_etd

mu_vec_etd <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "ETD",
                                   lambda_list = lambda_list_etd)
mu_vec_etd

shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 2, 7))
shift_matrix <- rbind(shift_matrix, c(2, 0, 8))
shift_matrix <- rbind(shift_matrix, c(1, 2, 9))
shift_matrix <- rbind(shift_matrix, c(2, 1, 10))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1, 2),
                                    num_concealed_states = 3,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix[q_matrix>=11]<-11
q_matrix

mu_vec_etd_sp<-rep(4, length(mu_vec_etd))

### Try 1
idparsopt <- c(1:4, 7:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1],3), rep(strt_mu[1],1), rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_sp_try1.rds")

### Try 2
idparsopt <- c(1:4, 7:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2],3), rep(strt_mu[2],1), rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_sp_try2.rds")

### Try 3
idparsopt <- c(1:4, 7:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3],3), rep(strt_mu[3],1), rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_sp_try3.rds")

### Try 1
idparsopt <- c(1, 4:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1],1), rep(strt_mu[1],3), rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_mu_try1.rds")

### Try 2
idparsopt <- c(1, 4:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2],1), rep(strt_mu[2],3), rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_mu_try2.rds")

### Try 3
idparsopt <- c(1, 4:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3],1), rep(strt_mu[3],3), rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_mu_try3.rds")

### Try 1
idparsopt <- c(1:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1],3), rep(strt_mu[1],3), rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 11
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_net_diversification_try1.rds")

### Try 2
idparsopt <- c(1:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2],3), rep(strt_mu[2],3), rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 11
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_net_diversification_try2.rds")

### Try 3
idparsopt <- c(1:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3],3), rep(strt_mu[3],3), rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 11
saveRDS(model, "SecSSE_Results/Diet_results/diet_examined_trait_net_diversification_try3.rds")

lambda_list_ctd <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")
lambda_list_ctd

mu_vec_ctd <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "CTD",
                                   lambda_list = lambda_list_ctd)
mu_vec_ctd

mu_vec_ctd_sp<-rep(4, length(mu_vec_ctd))

### Try 1
idparsopt <- c(1:4, 7:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1],3), rep(strt_mu[1],1), rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_sp_try1.rds")

### Try 2
idparsopt <- c(1:4, 7:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2],3), rep(strt_mu[2],1), rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_sp_try2.rds")

### Try 3
idparsopt <- c(1:4, 7:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3],3), rep(strt_mu[3],1), rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_sp_try3.rds")

### Try 1
idparsopt <- c(1, 4:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1],1), rep(strt_mu[1],3), rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_mu_try1.rds")

### Try 1
idparsopt <- c(1, 4:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2],1), rep(strt_mu[2],3), rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_mu_try2.rds")

### Try 3
idparsopt <- c(1, 4:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3],1), rep(strt_mu[3],3), rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 9
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_mu_try3.rds")

### Try 1
idparsopt <- c(1:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1],3), rep(strt_mu[1],3), rep(strt_q[1],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 11
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_net_diversification_try1.rds")

### Try 2
idparsopt <- c(1:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2],3), rep(strt_mu[2],3), rep(strt_q[2],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 11
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_net_diversification_try2.rds")

### Try 3
idparsopt <- c(1:11) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3],3), rep(strt_mu[3],3), rep(strt_q[3],5))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 3,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 11
saveRDS(model, "SecSSE_Results/Diet_results/diet_concealed_trait_net_diversification_try3.rds")
