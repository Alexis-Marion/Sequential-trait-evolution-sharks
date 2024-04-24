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
library("qpcR")

# Loading data

df<-read.csv("Species_data.tsv", sep="\t") # omit sep ="\t" for .csv files

phy<-read.nexus("16FC_16C_374_sp.tree")

aicc<-function(ll,k){round((2*(-round(ll,4))+2*k+(2*k*(k+1))/(374-k-1)),3)}

clnames<-c("Log Likelihood", "Parameters", "AICc", "Delta AICc", "AICcWt")

mb = Mclust(as.numeric(df$Body.size))
summary(mb, parameters = TRUE)
df$Body.size <- as.factor(mb$classification)

# Optional step (cleaning data)
states<-df$Body.size
names(states)=str_replace((df$Species), " ", "_")
states2<-states[!names(states) %in% setdiff(names(states), phy$tip.label)]

states2<-states[!names(states) %in% setdiff(names(states), phy$tip.label)]

states3<-as.data.frame(cbind(names(states2), as.factor(states2)))
colnames(states3)<-c("species","states")
rownames(states3)<-NULL
traits <- sortingtraits(states3, phy)
traits<-traits-1

# Computing sampling fractions
f<-c(
length(states2[states2=="1"])/length(df$Body.size[df$Body.size=="1"]),
length(states2[states2=="2"])/length(df$Body.size[df$Body.size=="2"]),
length(states2[states2=="3"])/length(df$Body.size[df$Body.size=="3"]))

spec_matrix_cr <- c()
spec_matrix_cr <- rbind(spec_matrix_cr, c(0, 0, 0, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(1, 1, 1, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(2, 2, 2, 1))

lambda_list_cr <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix_cr,
                                          model = "CR")
mu_vec_cr <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "CR",
                                   lambda_list = lambda_list_cr)

shift_matrix_cr  <- c()
shift_matrix_cr  <- rbind(shift_matrix_cr , c(0, 1, 3))
shift_matrix_cr  <- rbind(shift_matrix_cr , c(1, 0, 3))
shift_matrix_cr  <- rbind(shift_matrix_cr , c(1, 2, 4))
shift_matrix_cr  <- rbind(shift_matrix_cr , c(2, 1, 4))

q_matrix_cr <- secsse::create_q_matrix(state_names = c(0, 1, 2),
                                    num_concealed_states = 3,
                                    shift_matrix = shift_matrix_cr,
                                    diff.conceal = TRUE)
q_matrix_cr[q_matrix_cr>=5]<-5

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
spec_matrix <- rbind(spec_matrix, c(2, 2, 2, 3))

lambda_list_etd <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "ETD")

mu_vec_etd <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "ETD",
                                   lambda_list = lambda_list_etd)

mu_vec_etd_sp<-rep(4, length(mu_vec_etd))

lambda_list_ctd <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")

mu_vec_ctd <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "CTD",
                                   lambda_list = lambda_list_ctd)

### Speciation rates

mu_vec_ctd_sp<-rep(4, length(mu_vec_ctd))

shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 1, 7))
shift_matrix <- rbind(shift_matrix, c(1, 0, 7))
shift_matrix <- rbind(shift_matrix, c(1, 2, 8))
shift_matrix <- rbind(shift_matrix, c(2, 1, 8))


q_matrix <- secsse::create_q_matrix(state_names = c(0, 1, 2),
                                    num_concealed_states = 3,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix[q_matrix>=9]<-9

vec_AICc<-c()
vec_k<-c()
vec_L<-c()
vec_name<-c()
for (file in list.files("SecSSE_Results/Size_results/", full.names = TRUE, pattern = "size")){
    vec_AICc<-c(vec_AICc, aicc(readRDS(file)$ML,readRDS(file)$k))
    vec_k<-c(vec_k, readRDS(file)$k)
    vec_L<-c(vec_L, readRDS(file)$ML)
    vec_name<-c(vec_name, file)
}

data_max_Body_size<-cbind(vec_L,vec_k,vec_AICc, vec_name)

vec_AICc2<-c()
vec_k2<-c()
vec_L2<-c()
vec_name2<-c()

for (i in 1:nrow(data_max_Body_size)){
    x<-as.numeric(data_max_Body_size[i,3])
    if (i%%3 == 0){
        minval<-which.min(c(data_max_Body_size[i-2,3], data_max_Body_size[i-1,3], data_max_Body_size[i,3]))       
        minval<-i+minval-3
        vec_AICc2<-c(vec_AICc2, vec_AICc[minval])
        vec_k2<-c(vec_k2, vec_k[minval])
        vec_L2<-c(vec_L2, vec_L[minval])
        vec_name2<-c(vec_name2, vec_name[minval])
            
        if (grepl("_model_",vec_name[minval])){
            idparslist <- list()
            idparslist[[1]] <- lambda_list_cr
            idparslist[[2]] <- mu_vec_cr
            idparslist[[3]] <- q_matrix_cr
        }

        else{
            idparslist <- list()
            if (grepl("_examined_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_etd
                idparslist[[2]] <- mu_vec_etd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_etd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }
            if (grepl("_concealed_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_ctd
                idparslist[[2]] <- mu_vec_ctd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_ctd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }

            idparslist[[3]] <- q_matrix
        }
        if((min(vec_AICc2)==vec_AICc[minval])){
            print(vec_name[minval])
            temp_list<-idparslist
        }
        print(round(secsse::extract_par_vals(idparslist, readRDS(vec_name[minval])$MLpars), 3))
    }
}



data_max_Body_size<-data.frame(cbind(round(vec_L2,2),vec_k2,round(vec_AICc2,2),round(akaike.weights(vec_AICc2)$deltaAIC, 2), round(akaike.weights(vec_AICc2)$weights,2)))

colnames(data_max_Body_size)<-clnames

write.table(data_max_Body_size[order(data_max_Body_size[,3]),], "SecSSE_Results/data_max_Body_size.tsv", sep ="\t")

model<-readRDS("SecSSE_Results/Size_results/size_concealed_trait_sp_try3.rds")

CTD_par <- secsse::extract_par_vals(temp_list, model$MLpars)

CTD_par

spec_rates <- CTD_par[1:3]
ext_rates <- CTD_par[4:6]
Q_Examined <- CTD_par[7:8]
Q_Concealed <- CTD_par[9]

model$MLpars[[3]]

params <- secsse::id_paramPos(c(0, 1, 2), 3)

params[[1]][1:3]<-spec_rates[1]
params[[1]][4:6]<-spec_rates[2]
params[[1]][7:9]<-spec_rates[3]
params[[2]][1:3]<-ext_rates[1]
params[[2]][4:6]<-ext_rates[1]
params[[2]][7:9]<-ext_rates[1]
params[[3]]<-model$MLpars[[3]]

ll <- secsse::secsse_loglik(parameter = params,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = 3,
                             see_ancestral_states = TRUE,
                             sampling_fraction = f,
                             num_threads =  8)

saveRDS(ll$states, "ASE/bds_SecSSE.rds")

# Optional step (cleaning data)
states<-df$Reproduction
names(states)=str_replace((df$Species), " ", "_")
states2<-states[!names(states) %in% setdiff(names(states), phy$tip.label)]

states2<-states[!names(states) %in% setdiff(names(states), phy$tip.label)]

states3<-as.data.frame(cbind(names(states2), as.factor(states2)))
colnames(states3)<-c("species","states")
rownames(states3)<-NULL
traits <- sortingtraits(states3, phy)
traits<-traits-1

# Computing sampling fractions
f<-c(
length(na.omit(states2[states2=="O"]))/length(na.omit(df$Reproduction[df$Reproduction=="O"])),
length(na.omit(states2[states2=="Oo"]))/length(na.omit(df$Reproduction[df$Reproduction=="Oo"])),
length(na.omit(states2[states2=="PV"]))/length(na.omit(df$Reproduction[df$Reproduction=="PV"])),
length(na.omit(states2[states2=="YV"]))/length(na.omit(df$Reproduction[df$Reproduction=="YV"])))

spec_matrix_cr <- c()
spec_matrix_cr <- rbind(spec_matrix_cr, c(0, 0, 0, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(1, 1, 1, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(2, 2, 2, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(3, 3, 3, 1))

lambda_list_cr <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3),
                                          num_concealed_states = 4,
                                          transition_matrix = spec_matrix_cr,
                                          model = "CR")

mu_vec_cr <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3),
                                   num_concealed_states = 4,
                                   model = "CR",
                                   lambda_list = lambda_list_cr)

shift_matrix_cr <- c()
shift_matrix_cr <- rbind(shift_matrix_cr, c(0, 3, 3))
shift_matrix_cr <- rbind(shift_matrix_cr, c(3, 0, 3))
shift_matrix_cr <- rbind(shift_matrix_cr, c(1, 3, 3))
shift_matrix_cr <- rbind(shift_matrix_cr, c(3, 1, 3))
shift_matrix_cr <- rbind(shift_matrix_cr, c(2, 3, 3))
shift_matrix_cr <- rbind(shift_matrix_cr, c(3, 2, 3))

q_matrix_cr <- secsse::create_q_matrix(state_names = c(0, 1, 2, 3),
                                    num_concealed_states = 4,
                                    shift_matrix = shift_matrix_cr,
                                    diff.conceal = TRUE)
q_matrix_cr[q_matrix_cr>=4]<-4

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
spec_matrix <- rbind(spec_matrix, c(2, 2, 2, 3))
spec_matrix <- rbind(spec_matrix, c(3, 3, 3, 4))

lambda_list_etd <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3),
                                          num_concealed_states = 4,
                                          transition_matrix = spec_matrix,
                                          model = "ETD")

mu_vec_etd <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3),
                                   num_concealed_states = 4,
                                   model = "ETD",
                                   lambda_list = lambda_list_etd)

shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 3, 9))
shift_matrix <- rbind(shift_matrix, c(3, 0, 9))
shift_matrix <- rbind(shift_matrix, c(1, 3, 9))
shift_matrix <- rbind(shift_matrix, c(3, 1, 9))
shift_matrix <- rbind(shift_matrix, c(2, 3, 9))
shift_matrix <- rbind(shift_matrix, c(3, 2, 9))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1, 2, 3),
                                    num_concealed_states = 4,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = TRUE)
q_matrix[q_matrix>=10]<-10

mu_vec_etd_sp<-rep(5, length(mu_vec_etd))

lambda_list_ctd <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3),
                                          num_concealed_states = 4,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")

mu_vec_ctd <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3),
                                   num_concealed_states = 4,
                                   model = "CTD",
                                   lambda_list = lambda_list_ctd)

mu_vec_ctd_sp<-rep(5, length(mu_vec_ctd))

vec_AICc<-c()
vec_k<-c()
vec_L<-c()
vec_name<-c()
for (file in list.files("SecSSE_Results/Reproduction_results/", full.names = TRUE, pattern = "reproduction")){
    vec_AICc<-c(vec_AICc, aicc(readRDS(file)$ML,readRDS(file)$k))
    vec_k<-c(vec_k, readRDS(file)$k)
    vec_L<-c(vec_L, readRDS(file)$ML)
    vec_name<-c(vec_name, file)
}

data_max_Reproduction<-cbind(vec_L,vec_k,vec_AICc, vec_name)

vec_AICc2<-c()
vec_k2<-c()
vec_L2<-c()
vec_name2<-c()

for (i in 1:nrow(data_max_Reproduction)){
    x<-as.numeric(data_max_Reproduction[i,3])
    if (i%%3 == 0){
        minval<-which.min(c(data_max_Reproduction[i-2,3], data_max_Reproduction[i-1,3], data_max_Reproduction[i,3]))       
        minval<-i+minval-3
        vec_AICc2<-c(vec_AICc2, vec_AICc[minval])
        vec_k2<-c(vec_k2, vec_k[minval])
        vec_L2<-c(vec_L2, vec_L[minval])
        vec_name2<-c(vec_name2, vec_name[minval])
            
        if (grepl("_model_",vec_name[minval])){
            idparslist <- list()
            idparslist[[1]] <- lambda_list_cr
            idparslist[[2]] <- mu_vec_cr
            idparslist[[3]] <- q_matrix_cr
        }

        else{
            idparslist <- list()
            if (grepl("_examined_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_etd
                idparslist[[2]] <- mu_vec_etd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_etd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }
            if (grepl("_concealed_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_ctd
                idparslist[[2]] <- mu_vec_ctd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_ctd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }

            idparslist[[3]] <- q_matrix
        }
        if((min(vec_AICc2)==vec_AICc[minval])){
            print(vec_name[minval])
            temp_list<-idparslist
        }
        print(round(secsse::extract_par_vals(idparslist, readRDS(vec_name[minval])$MLpars),3))
    }
}



data_max_Reproduction<-data.frame(cbind(round(vec_L2,2),vec_k2,round(vec_AICc2,2),round(akaike.weights(vec_AICc2)$deltaAIC, 2), round(akaike.weights(vec_AICc2)$weights,2)))

colnames(data_max_Reproduction)<-clnames

write.table(data_max_Reproduction[order(data_max_Reproduction[,3]),], "SecSSE_Results/data_max_Reproduction.tsv", sep ="\t")

model<-readRDS("SecSSE_Results/Reproduction_results/reproduction_concealed_trait_sp_try1.rds")

CTD_par <- secsse::extract_par_vals(temp_list, model$MLpars)

CTD_par

spec_rates <- CTD_par[1:4]
ext_rates <- CTD_par[5:8]

model$MLpars[[3]]

params <- secsse::id_paramPos(c(0, 1, 2, 3), 4)

params[[1]][1:4]<-spec_rates[1]
params[[1]][5:8]<-spec_rates[2]
params[[1]][9:12]<-spec_rates[3]
params[[1]][13:16]<-spec_rates[4]
params[[2]][1:4]<-ext_rates[1]
params[[2]][5:8]<-ext_rates[1]
params[[2]][9:12]<-ext_rates[1]
params[[2]][13:16]<-ext_rates[1]
params[[3]]<-model$MLpars[[3]]

ll <- secsse::secsse_loglik(parameter = params,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = 4,
                             see_ancestral_states = TRUE,
                             sampling_fraction = f,
                             num_threads =  8)
ll

saveRDS(ll$states, "ASE/rp_SecSSE.rds")

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

spec_matrix_cr <- c()
spec_matrix_cr <- rbind(spec_matrix_cr, c(0, 0, 0, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(1, 1, 1, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(2, 2, 2, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(3, 3, 3, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(4, 4, 4, 1))
lambda_list_cr <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3, 4),
                                          num_concealed_states = 5,
                                          transition_matrix = spec_matrix_cr,
                                          model = "CR")
mu_vec_cr <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3, 4),
                                   num_concealed_states = 5,
                                   model = "CR",
                                   lambda_list = lambda_list_cr)
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
mu_vec_etd <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3, 4),
                                   num_concealed_states = 5,
                                   model = "ETD",
                                   lambda_list = lambda_list_etd)
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
mu_vec_etd_sp<-rep(6, length(mu_vec_etd))
lambda_list_ctd <- secsse::create_lambda_list(state_names = c(0, 1, 2, 3, 4),
                                          num_concealed_states = 5,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")
mu_vec_ctd <- secsse::create_mu_vector(state_names = c(0, 1, 2, 3, 4),
                                   num_concealed_states = 5,
                                   model = "CTD",
                                   lambda_list = lambda_list_ctd)
mu_vec_ctd_sp<-rep(6, length(mu_vec_ctd))

vec_AICc<-c()
vec_k<-c()
vec_L<-c()
vec_name<-c()
for (file in list.files("SecSSE_Results/Habitat_results/", full.names = TRUE, pattern = "habitat")){
    vec_AICc<-c(vec_AICc, aicc(readRDS(file)$ML,readRDS(file)$k))
    vec_k<-c(vec_k, readRDS(file)$k)
    vec_L<-c(vec_L, readRDS(file)$ML)
    vec_name<-c(vec_name, file)
}

data_max_Habitat<-cbind(vec_L,vec_k,vec_AICc, vec_name)

vec_AICc2<-c()
vec_k2<-c()
vec_L2<-c()
vec_name2<-c()

for (i in 1:nrow(data_max_Habitat)){
    x<-as.numeric(data_max_Habitat[i,3])
    if (i%%3 == 0){
        minval<-which.min(c(data_max_Habitat[i-2,3], data_max_Habitat[i-1,3], data_max_Habitat[i,3]))       
        minval<-i+minval-3
        vec_AICc2<-c(vec_AICc2, vec_AICc[minval])
        vec_k2<-c(vec_k2, vec_k[minval])
        vec_L2<-c(vec_L2, vec_L[minval])
        vec_name2<-c(vec_name2, vec_name[minval])
            
        if (grepl("_model_",vec_name[minval])){
            idparslist <- list()
            idparslist[[1]] <- lambda_list_cr
            idparslist[[2]] <- mu_vec_cr
            idparslist[[3]] <- q_matrix_cr
        }

        else{
            idparslist <- list()
            if (grepl("_examined_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_etd
                idparslist[[2]] <- mu_vec_etd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_etd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }
            if (grepl("_concealed_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_ctd
                idparslist[[2]] <- mu_vec_ctd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_ctd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }

            idparslist[[3]] <- q_matrix
        }
        if((min(vec_AICc2)==vec_AICc[minval])){
            print(vec_name[minval])
            temp_list<-idparslist
        }
        print(round(secsse::extract_par_vals(idparslist, readRDS(vec_name[minval])$MLpars), 3))
    }
}



data_max_Habitat<-data.frame(cbind(round(vec_L2,2),vec_k2,round(vec_AICc2,2),round(akaike.weights(vec_AICc2)$deltaAIC, 2), round(akaike.weights(vec_AICc2)$weights,2)))

colnames(data_max_Habitat)<-clnames

write.table(data_max_Habitat[order(data_max_Habitat[,3]),], "SecSSE_Results/data_max_Habitat.tsv", sep ="\t")

model<-readRDS("SecSSE_Results/Habitat_results/habitat_concealed_trait_sp_try3.rds")

CTD_par <- secsse::extract_par_vals(temp_list, model$MLpars)

CTD_par

spec_rates <- CTD_par[1:5]
ext_rates <- CTD_par[6:10]

model$MLpars[[3]]

params <- secsse::id_paramPos(c(0, 1, 2, 3, 4), 5)

params[[1]][1:5]<-spec_rates[1]
params[[1]][6:10]<-spec_rates[2]
params[[1]][11:15]<-spec_rates[3]
params[[1]][16:20]<-spec_rates[4]
params[[1]][21:25]<-spec_rates[5]
params[[2]][1:5]<-ext_rates[1]
params[[2]][6:10]<-ext_rates[1]
params[[2]][11:16]<-ext_rates[1]
params[[2]][16:20]<-ext_rates[1]
params[[1]][21:25]<-spec_rates[1]
params[[3]]<-model$MLpars[[3]]

ll <- secsse::secsse_loglik(parameter = params,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = 5,
                             see_ancestral_states = TRUE,
                             sampling_fraction = f,
                             num_threads =  8)

saveRDS(ll$states, "ASE/ht_SecSSE.rds")

df<-read.csv("table_diet.tsv", sep="\t")

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

spec_matrix_cr <- c()
spec_matrix_cr <- rbind(spec_matrix_cr, c(0, 0, 0, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(1, 1, 1, 1))
spec_matrix_cr <- rbind(spec_matrix_cr, c(2, 2, 2, 1))

lambda_list_cr <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix_cr,
                                          model = "CR")
mu_vec_cr <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "CR",
                                   lambda_list = lambda_list_cr)

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

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 2))
spec_matrix <- rbind(spec_matrix, c(2, 2, 2, 3))

lambda_list_etd <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "ETD")

mu_vec_etd <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "ETD",
                                   lambda_list = lambda_list_etd)

mu_vec_etd_sp<-rep(4, length(mu_vec_etd))

lambda_list_ctd <- secsse::create_lambda_list(state_names = c(0, 1, 2),
                                          num_concealed_states = 3,
                                          transition_matrix = spec_matrix,
                                          model = "CTD")

mu_vec_ctd <- secsse::create_mu_vector(state_names = c(0, 1, 2),
                                   num_concealed_states = 3,
                                   model = "CTD",
                                   lambda_list = lambda_list_ctd)

### Speciation rates

mu_vec_ctd_sp<-rep(4, length(mu_vec_ctd))

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

vec_AICc<-c()
vec_k<-c()
vec_L<-c()
vec_name<-c()
for (file in list.files("SecSSE_Results/Diet_results/", full.names = TRUE, pattern = "diet")){
    vec_AICc<-c(vec_AICc, aicc(readRDS(file)$ML,readRDS(file)$k))
    vec_k<-c(vec_k, readRDS(file)$k)
    vec_L<-c(vec_L, readRDS(file)$ML)
    vec_name<-c(vec_name, file)
}

data_Diet<-cbind(vec_L,vec_k,vec_AICc, vec_name)

vec_AICc2<-c()
vec_k2<-c()
vec_L2<-c()
vec_name2<-c()

for (i in 1:nrow(data_Diet)){
    x<-as.numeric(data_Diet[i,3])
    if (i%%3 == 0){
        minval<-which.min(c(data_Diet[i-2,3], data_Diet[i-1,3], data_Diet[i,3]))       
        minval<-i+minval-3
        vec_AICc2<-c(vec_AICc2, vec_AICc[minval])
        vec_k2<-c(vec_k2, vec_k[minval])
        vec_L2<-c(vec_L2, vec_L[minval])
        vec_name2<-c(vec_name2, vec_name[minval])
            
        if (grepl("_model_",vec_name[minval])){
            idparslist <- list()
            idparslist[[1]] <- lambda_list_cr
            idparslist[[2]] <- mu_vec_cr
            idparslist[[3]] <- q_matrix_cr
        }

        else{
            idparslist <- list()
            if (grepl("_examined_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_etd
                idparslist[[2]] <- mu_vec_etd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_etd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }
            if (grepl("_concealed_",vec_name[minval])){
                idparslist <- list()
                idparslist[[1]] <- lambda_list_ctd
                idparslist[[2]] <- mu_vec_ctd
                if (grepl("_sp_",vec_name[minval])){
                    idparslist[[2]] <- mu_vec_ctd_sp
                }
                if (grepl("_ex_",vec_name[minval])){
                    idparslist[[1]] <- lambda_list_cr
                }
            }

            idparslist[[3]] <- q_matrix
        }
        if((min(vec_AICc2)==vec_AICc[minval])){
            print(vec_name[minval])
            temp_list<-idparslist
        }
        print(round(secsse::extract_par_vals(idparslist, readRDS(vec_name[minval])$MLpars), 3))
    }
}



data_Diet<-data.frame(cbind(round(vec_L2,2),vec_k2,round(vec_AICc2,2),round(akaike.weights(vec_AICc2)$deltaAIC, 2), round(akaike.weights(vec_AICc2)$weights,2)))

colnames(data_Diet)<-clnames

write.table(data_Diet[order(data_Diet[,3]),], "SecSSE_Results/data_Diet.tsv", sep ="\t")

model<-readRDS("SecSSE_Results/Diet_results/diet_concealed_trait_sp_try1.rds")

CTD_par <- secsse::extract_par_vals(temp_list, model$MLpars)

CTD_par

spec_rates <- CTD_par[1:3]
ext_rates <- CTD_par[4:6]
Q_Examined <- CTD_par[7:10]
Q_Concealed <- CTD_par[11]

model$MLpars[[3]]

params <- secsse::id_paramPos(c(0, 1, 2), 3)

params[[1]][1:3]<-spec_rates[1]
params[[1]][4:6]<-spec_rates[2]
params[[1]][7:9]<-spec_rates[3]
params[[2]][1:3]<-ext_rates[1]
params[[2]][4:6]<-ext_rates[1]
params[[2]][7:9]<-ext_rates[1]
params[[3]]<-model$MLpars[[3]]

ll <- secsse::secsse_loglik(parameter = params,
                             phy = phy,
                             traits = traits,
                             num_concealed_states = 3,
                             see_ancestral_states = TRUE,
                             sampling_fraction = f,
                             num_threads =  8)

saveRDS(ll$states, "ASE_diet.rds")
