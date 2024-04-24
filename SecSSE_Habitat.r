# Loading packages
library("mclust")
library("Rtsne")
library("ape")
library("secsse")
library("DDD")
library("doMC")
library("tidyverse")
library("parallel")
library("qgraph")
library("stringr")

# Loading data

df<-read.csv("Species_data.tsv", sep="\t") # omit sep ="\t" for .csv files

phy<-read.nexus("Data2/16FC_16C_374_sp.tree")

# Optional step (cleaning data)
states<-df$Habitat
names(states)=str_replace((df$Species), " ", "_")
states2<-states[!names(states) %in% setdiff(names(states), phy$tip.label)]

setdiff(phy$tip.label,names(states2))

states3<-as.data.frame(cbind(names(states2), as.factor(states2)))
colnames(states3)<-c("species","states")
rownames(states3)<-NULL
traits <- sortingtraits(states3, phy)
traits<-traits-1

# Computing sampling fractions
f<-c(
length(na.omit(states2[states2=="D"]))/length(na.omit(df$Habitat[df$Habitat=="D"])),
length(na.omit(states2[states2=="IS"]))/length(na.omit(df$Habitat[df$Habitat=="IS"])),
length(na.omit(states2[states2=="O"]))/length(na.omit(df$Habitat[df$Habitat=="O"])),
length(na.omit(states2[states2=="OS"]))/length(na.omit(df$Habitat[df$Habitat=="OS"])),
length(na.omit(states2[states2=="R"]))/length(na.omit(df$Habitat[df$Habitat=="R"])))

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
spec_matrix <- rbind(spec_matrix, c(3, 3, 3, 1))
spec_matrix <- rbind(spec_matrix, c(4, 4, 4, 1))

lambda_list_cr <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3, 4),
                                          num_concealed_states = 5,
                                          transition_matrix = spec_matrix,
                                          model = "CR")
lambda_list_cr

mu_vec_cr <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3, 4),
                                   num_concealed_states = 5,
                                   model = "CR",
                                   lambda_list = lambda_list_cr)
mu_vec_cr

shift_matrix_cr <- c()
shift_matrix_cr <- rbind(shift_matrix_cr, c(0, 3, 3))
shift_matrix_cr <- rbind(shift_matrix_cr, c(1, 2, 4))
shift_matrix_cr <- rbind(shift_matrix_cr, c(1, 3, 5))
shift_matrix_cr <- rbind(shift_matrix_cr, c(1, 4, 6))
shift_matrix_cr <- rbind(shift_matrix_cr, c(2, 3, 7))
shift_matrix_cr <- rbind(shift_matrix_cr, c(2, 4, 8))
shift_matrix_cr <- rbind(shift_matrix_cr, c(3, 0, 9))
shift_matrix_cr <- rbind(shift_matrix_cr, c(3, 1, 10))
shift_matrix_cr <- rbind(shift_matrix_cr, c(4, 1, 11))
shift_matrix_cr <- rbind(shift_matrix_cr, c(4, 2, 12))
shift_matrix_cr <- rbind(shift_matrix_cr, c(4, 3, 13))

q_matrix_cr <- secsse::create_q_matrix(state_names = c(0, 1, 2, 3, 4),
                                    num_concealed_states = 5,
                                    shift_matrix = shift_matrix_cr,
                                    diff.conceal = TRUE)
q_matrix_cr[q_matrix_cr>=14]<-14
q_matrix_cr

### Try 1

idparsopt <- c(1:14) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(strt_lambda[1], strt_mu[1], rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_cr
idparslist[[3]] <- q_matrix_cr

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 14
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_model_try1.rds")

### Try 2

idparsopt <- c(1:14) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(strt_lambda[2], strt_mu[2], rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_cr
idparslist[[3]] <- q_matrix_cr

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 14
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_model_try2.rds")

### Try 3

idparsopt <- c(1:14) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(strt_lambda[3], strt_mu[3], rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_cr
idparslist[[3]] <- q_matrix_cr

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 14
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_model_try3.rds")

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
spec_matrix <- rbind(spec_matrix, c(2, 2, 2, 3))
spec_matrix <- rbind(spec_matrix, c(3, 3, 3, 4))
spec_matrix <- rbind(spec_matrix, c(4, 4, 4, 5))

lambda_list_etd <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3, 4),
                                          num_concealed_states = 5,
                                          transition_matrix = spec_matrix,
                                          model = "ETD")
lambda_list_etd

mu_vec_etd <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3, 4),
                                   num_concealed_states = 5,
                                   model = "ETD",
                                   lambda_list = lambda_list_etd)
mu_vec_etd

shift_matrix <- c()
shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 3, 11))
shift_matrix <- rbind(shift_matrix, c(1, 2, 12))
shift_matrix <- rbind(shift_matrix, c(1, 3, 13))
shift_matrix <- rbind(shift_matrix, c(1, 4, 14))
shift_matrix <- rbind(shift_matrix, c(2, 3, 15))
shift_matrix <- rbind(shift_matrix, c(2, 4, 16))
shift_matrix <- rbind(shift_matrix, c(3, 0, 17))
shift_matrix <- rbind(shift_matrix, c(3, 1, 18))
shift_matrix <- rbind(shift_matrix, c(4, 1, 19))
shift_matrix <- rbind(shift_matrix, c(4, 2, 20))
shift_matrix <- rbind(shift_matrix, c(4, 3, 21))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1, 2, 3, 4),
                                    num_concealed_states = 5,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix[q_matrix>=22]<-22
q_matrix

mu_vec_etd_sp<-rep(6, length(mu_vec_etd))

### Try 1
idparsopt <- c(1:6, 11:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1], 5), rep(strt_mu[1], 1), rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_sp_try1.rds")

### Try 2
idparsopt <- c(1:6, 11:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2], 5), rep(strt_mu[2], 1), rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_sp_try2.rds")

### Try 3
idparsopt <- c(1:6, 11:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3], 5), rep(strt_mu[3], 1), rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_sp_try3.rds")

### Try 1
idparsopt <- c(1, 6:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1], 1), rep(strt_mu[1], 5), rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_mu_try1.rds")

### Try 2
idparsopt <- c(1, 6:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2], 1), rep(strt_mu[2], 5), rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_mu_try2.rds")

### Try 3
idparsopt <- c(1, 6:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3], 1), rep(strt_mu[3], 5), rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_mu_try3.rds")

### Try 1
idparsopt <- c(1:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1], 5), rep(strt_mu[1], 5), rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 22
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_net_diversification_try1.rds")

### Try 2
idparsopt <- c(1:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2], 5), rep(strt_mu[2], 5), rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 22
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_net_diversification_try2.rds")

### Try 3
idparsopt <- c(1:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3], 5), rep(strt_mu[3], 5), rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_etd
idparslist[[2]] <- mu_vec_etd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 22
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_examined_trait_net_diversification_try3.rds")

lambda_list_ctd <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3, 4),
                                          num_concealed_states = 5,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")
lambda_list_ctd

mu_vec_ctd <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3, 4),
                                   num_concealed_states = 5,
                                   model = "CTD",
                                   lambda_list = lambda_list_ctd)
mu_vec_ctd

mu_vec_ctd_sp<-rep(6, length(mu_vec_ctd))

### Try 1
idparsopt <- c(1:6, 11:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1], 5), rep(strt_mu[1], 1), rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_sp_try1.rds")

### Try 2
idparsopt <- c(1:6, 11:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2], 5), rep(strt_mu[2], 1), rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_sp_try2.rds")

### Try 2
idparsopt <- c(1:6, 11:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3], 5), rep(strt_mu[3], 1), rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd_sp
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_sp_try3.rds")

### Try 1
idparsopt <- c(1, 6:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1], 1), rep(strt_mu[1], 5), rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_mu_try1.rds")

### Try 2
idparsopt <- c(1, 6:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2], 1), rep(strt_mu[2], 5), rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_mu_try2.rds")

### Try 3
idparsopt <- c(1, 6:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3], 1), rep(strt_mu[3], 5), rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_cr
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 18
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_mu_try3.rds")

### Try 1
idparsopt <- c(1:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[1], 5), rep(strt_mu[1], 5), rep(strt_q[1], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 22
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_net_diversification_try1.rds")

### Try 2
idparsopt <- c(1:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[2], 5), rep(strt_mu[2], 5), rep(strt_q[2], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 22
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_net_diversification_try2.rds")

### Try 3
idparsopt <- c(1:22) # our maximum rate parameter was 3
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- c(rep(strt_lambda[3], 5), rep(strt_mu[3], 5), rep(strt_q[3], 12))
initparsfix <- c(0.0) # all zeros remain at zero.
idparslist <- list()
idparslist[[1]] <- lambda_list_ctd
idparslist[[2]] <- mu_vec_ctd
idparslist[[3]] <- q_matrix

model <- secsse::cla_secsse_ml(phy = phy,
                              traits = traits,
                              num_concealed_states = 5,
                              idparslist = idparslist,
                              idparsopt = idparsopt,
                              initparsopt = initparsopt,
                              idparsfix = idparsfix,
                              parsfix = initparsfix,
                              sampling_fraction = f,
                              verbose = FALSE,
                              num_threads =  8)
model$k <- 22
saveRDS(model, "SecSSE_Results/Habitat_results/habitat_concealed_trait_net_diversification_try3.rds")
