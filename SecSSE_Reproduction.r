## SecSSE : Reproduction 

## I Data preparation

### Loading packages

library("secsse")
library("stringr")
library("qpcR")
library("ape")
library("DDD")

### Loading data

df<-read.csv("Path/to/your/file.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.nexus("Path/to/your/tree.tree")

states<-df$Reproduction
names(states)=str_replace((rownames(df)), " ", "_")
states2<-states[!names(states) %in% setdiff(names(states), phy$tip.label)]

phy2.0<-drop.tip(phy, setdiff(phy$tip.label, (names(states))))

### Computing sampling fraction

samp_frac<-c(
length(states2[states2=="O"])/length(df$Reproduction[df$Reproduction=="O"]),
length(states2[states2=="Oo"])/length(df$Reproduction[df$Reproduction=="Oo"]),
length(states2[states2=="PV"])/length(df$Reproduction[df$Reproduction=="PV"]),
length(states2[states2=="YV"])/length(df$Reproduction[df$Reproduction=="YV"]))

### Naming states according to Secsse

states3<-as.data.frame(cbind(names(states2), as.factor(states2)))
colnames(states3)<-c("species","states")
rownames(states3)<-NULL
traits <- sortingtraits(states3, phy2.0)

traits
states3 0.0358367783113309 0.0237273957947573 0.0423775017413011 0.040433164010602

## II Running the analysis

### Initialisation

startingpoint <- bd_ML(brts = ape::branching.times(phy2.0))
intGuessLambda <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0
#COND=conditioning on the state of the root
cond<-"proper_cond"
#root_state_weight
root_state_weight<-"proper_weights"
#SAMPLING_FRACTION
#TOL
tol = c(1e-04, 1e-05, 1e-07)

### Constant rate model

  idparslist <- id_paramPos(traits, num_concealed_states = 4)
  idparslist[[1]][] <- 1 #same speciation rate
  idparslist[[2]][] <- 2 #same extinction rate
  
  #same transition rate
  masterBlock<-matrix(99,ncol=4,nrow=4,byrow=T) 
  diag(masterBlock) <- NA
  masterBlock[1,2] <- 3
  masterBlock[1,3] <- 3
  masterBlock[1,4] <- 3
  masterBlock[2,1] <- 3
  masterBlock[2,3] <- 3
  masterBlock[2,4] <- 3
  masterBlock[3,1] <- 3
  masterBlock[3,2] <- 3
  masterBlock[3,4] <- 3
  masterBlock[4,1] <- 3
  masterBlock[4,2] <- 3
  masterBlock[4,3] <- 3
  diff.conceal <- TRUE
  myQ<-q_doubletrans(traits,masterBlock,diff.conceal)
  myQ<-replace(myQ, myQ == 4, 3)
  idparslist[[3]] <- myQ
  myQ

#INITPARSOPT = initiale values of each parameters
initparsopt <- c(intGuessLambda, intGuessMu, 0.25) # speciation rate, extinction rate, transition rate
#IDPARSOPT = parameters we want to optimise
idparsopt <- c(1:3)
#IDPARSFIX = parameters we want to fix
idparsfix <- 0
#PARSFIX = values of the fix parameters
parsfix<- 0

CR_sp<-secsse_ml(phy2.0, traits, num_concealed_states=4, idparslist,#here the matrix
          idparsopt, initparsopt, idparsfix, parsfix,
          cond=cond,root_state_weight = root_state_weight, tol = tol,
          sampling_fraction=samp_frac, 
          maxiter = 1000 * round((1.25)^length(idparsopt)), 
          use_fortran=TRUE,methode="ode45", optimmethod = "simplex", 
          num_cycles = 1, run_parallel=FALSE)

### CTD2 model (Two conceal states)

  idparslist <- id_paramPos(traits, num_concealed_states = 2)
  idparslist[[1]][1:4] <- 1 # Speciation varies only between concealed states
  idparslist[[1]][5:8] <- 2 
  idparslist[[2]][] <- 3 # same extinction rate

  diag(idparslist[[3]]) <- NA # remove same state transition rate (diagonal)
  # same transition rates
  idparslist[[3]][1,c(2,3,4,5)] <- 4
  idparslist[[3]][1,c(6,7,8)] <- 0
  idparslist[[3]][2,c(1,3,4,6)] <- 4
  idparslist[[3]][2,c(5,7,8)] <- 0
  idparslist[[3]][3,c(1,2,4,7)] <- 4
  idparslist[[3]][3,c(5,6,8)] <- 0
  idparslist[[3]][4,c(1,2,3,8)] <- 4
  idparslist[[3]][4,c(5,6,7)] <- 0
  idparslist[[3]][5,c(1,6,7,8)] <- 4
  idparslist[[3]][5,c(2,3,4)] <- 0
  idparslist[[3]][6,c(2,5,7,8)] <- 4
  idparslist[[3]][6,c(1,3,4)] <- 0
  idparslist[[3]][7,c(3,5,6,8)] <- 4
  idparslist[[3]][7,c(1,2,4)] <- 0
  idparslist[[3]][8,c(4,5,6,7)] <- 4
  idparslist[[3]][8,c(1,2,3,8)] <- 0
  idparslist[[3]]

#INITPARSOPT = initiale values of each parameters
initparsopt <- c(rep(intGuessLambda,2), intGuessMu, 0.25)# speciation rate, extinction rate, transition rate
#IDPARSOPT = parameters we want to optimise
idparsopt <- c(1:4)
#IDPARSFIX = parameters we want to fix
idparsfix <- 0
#PARSFIX = values of the fix parameters
parsfix<- 0

CTD2_sp<-secsse_ml(phy2.0,traits, num_concealed_states=2, idparslist,#here the matrix
          idparsopt, initparsopt, idparsfix, parsfix,
          cond=cond,root_state_weight = root_state_weight, tol = tol,
          sampling_fraction=samp_frac, 
          maxiter = 1000 * round((1.25)^length(idparsopt)), 
          use_fortran=TRUE,methode="ode45", optimmethod = "simplex", 
          num_cycles = 1, run_parallel=FALSE)

### CTD3 model (Three conceal states)

  idparslist <- id_paramPos(traits, num_concealed_states = 3)
  idparslist[[1]][1:4] <- 1 #Speciation varies only between concealed states
  idparslist[[1]][5:8] <- 2 
  idparslist[[1]][9:12] <- 3
  idparslist[[2]][] <- 4 #same extinction rate
  
  diag(idparslist[[3]]) <- NA#remove same state transition rate (diagonal)
  #all transition rates equal
  idparslist[[3]][1,c(2,3,4,5,9)] <- 5
  idparslist[[3]][1,c(6,7,8,10,11,12)] <- 0
  idparslist[[3]][2,c(1,3,4,6,10)] <- 5
  idparslist[[3]][2,c(5,7,8,9,11,12)] <- 0
  idparslist[[3]][3,c(1,2,4,7,11)] <- 5
  idparslist[[3]][3,c(5,6,8,9,10,12)] <- 0
  idparslist[[3]][4,c(1,2,3,8,12)] <- 5
  idparslist[[3]][4,c(5,6,7,9,10,11)] <- 0
  idparslist[[3]][5,c(1,6,7,8,9)] <- 5
  idparslist[[3]][5,c(2,3,4,10,11,12)] <- 0
  idparslist[[3]][6,c(2,5,7,8,10)] <- 5
  idparslist[[3]][6,c(1,3,4,9,11,12)] <- 0
  idparslist[[3]][7,c(3,5,6,8,11)] <- 5
  idparslist[[3]][7,c(1,2,4,9,10,12)] <- 0
  idparslist[[3]][8,c(4,5,6,7,12)] <- 5
  idparslist[[3]][8,c(1,2,3,9,10,11)] <- 0
  idparslist[[3]][9,c(1,5,10,11,12)] <- 5
  idparslist[[3]][9,c(2,3,4,6,7,8)] <- 0
  idparslist[[3]][10,c(2,6,9,11,12)] <- 5
  idparslist[[3]][10,c(1,3,4,5,7,8)] <- 0
  idparslist[[3]][11,c(3,7,9,10,12)] <- 5
  idparslist[[3]][11,c(1,2,4,5,6,8)] <- 0
  idparslist[[3]][12,c(4,8,9,10,11)] <- 5
  idparslist[[3]][12,c(1,2,3,5,6,7)] <- 0
  idparslist[[3]]

#INITPARSOPT=initiale values of each parameters
initparsopt <- c(rep(intGuessLambda,3), intGuessMu, 0.25)#speciation rate, transition rate
#IDPARSOPT = parameters we want to optimise
idparsopt <- c(1:5)
#IDPARSFIX = parameters we want to fix
idparsfix <- 0
#PARSFIX = values of the fix parameters
parsfix<- 0

CTD3_sp<-secsse_ml(phy2.0,traits, num_concealed_states=3, idparslist,#here the matrix
          idparsopt, initparsopt, idparsfix, parsfix,
          cond=cond,root_state_weight = root_state_weight, tol = tol,
          sampling_fraction=samp_frac, 
          maxiter = 1000 * round((1.25)^length(idparsopt)), 
          use_fortran=TRUE,methode="ode45", optimmethod = "simplex", 
          num_cycles = 1, run_parallel=FALSE)

### CTD4 model (Four conceal states)

  idparslist <- id_paramPos(traits, num_concealed_states = 4)
  idparslist[[1]][1:4] <- 1 #Speciation varies only between concealed states
  idparslist[[1]][5:8] <- 2 
  idparslist[[1]][9:12] <- 3  
  idparslist[[1]][13:16] <- 4 
  idparslist[[2]][] <- 5 #same extinction rate
  
  masterBlock<-matrix(99,ncol=4,nrow=4,byrow=T) 
  diag(masterBlock) <- NA
  masterBlock[1,2] <- 6
  masterBlock[1,3] <- 6
  masterBlock[1,4] <- 6
  masterBlock[2,1] <- 6
  masterBlock[2,3] <- 6
  masterBlock[2,4] <- 6
  masterBlock[3,1] <- 6
  masterBlock[3,2] <- 6
  masterBlock[3,4] <- 6
  masterBlock[4,1] <- 6
  masterBlock[4,2] <- 6
  masterBlock[4,3] <- 6
  diff.conceal <- TRUE
  myQ<-q_doubletrans(traits,masterBlock,diff.conceal)
  myQ<-replace(myQ, myQ == 7, 6)
  idparslist[[3]] <- myQ
  myQ

#INITPARSOPT=initiale values of each parameters
initparsopt <- c(rep(intGuessLambda,4), intGuessMu, 0.25)#speciation rate, transition rate
#IDPARSOPT = parameters we want to optimise
idparsopt <- c(1:6)
#IDPARSFIX = parameters we want to fix
idparsfix <- 0
#PARSFIX = values of the fix parameters
parsfix<- 0

CTD4_sp<-secsse_ml(phy2.0,traits, num_concealed_states=4, idparslist,#here the matrix
          idparsopt, initparsopt, idparsfix, parsfix,
          cond=cond,root_state_weight = root_state_weight, tol = tol,
          sampling_fraction=samp_frac, 
          maxiter = 1000 * round((1.25)^length(idparsopt)), 
          use_fortran=TRUE,methode="ode45", optimmethod = "simplex", 
          num_cycles = 1, run_parallel=FALSE)

### ETD model (Focal trait)

 idparslist <- id_paramPos(traits, num_concealed_states = 4)
  idparslist[[1]][c(1,5,9,13)] <- 1 #Speciation varies only between the 4 states
  idparslist[[1]][c(2,6,10,14)] <- 2 
  idparslist[[1]][c(3,7,11,15)] <- 3
  idparslist[[1]][c(4,8,12,16)] <- 4
  idparslist[[2]][] <- 5 #same extinction rate
  
  #all transition rates equal
  masterBlock<-matrix(99,ncol=4,nrow=4,byrow=T) 
  diag(masterBlock) <- NA
  masterBlock[1,2] <- 6
  masterBlock[1,3] <- 6
  masterBlock[1,4] <- 6
  masterBlock[2,1] <- 6
  masterBlock[2,3] <- 6
  masterBlock[2,4] <- 6
  masterBlock[3,1] <- 6
  masterBlock[3,2] <- 6
  masterBlock[3,4] <- 6
  masterBlock[4,1] <- 6
  masterBlock[4,2] <- 6
  masterBlock[4,3] <- 6
  diff.conceal <- TRUE
  myQ<-q_doubletrans(traits,masterBlock,diff.conceal)
  myQ<-replace(myQ, myQ == 7, 6)
  idparslist[[3]] <- myQ
  myQ

#INITPARSOPT=initiale values of each parameters
initparsopt <- c(rep(intGuessLambda,4), intGuessMu, 0.25)#speciation rate, transition rate
#IDPARSOPT = parameters we want to optimise
idparsopt <- c(1:6)
#IDPARSFIX = parameters we want to fix
idparsfix <- 0
#PARSFIX = values of the fix parameters
parsfix<- 0

ETD_sp<-secsse_ml(phy2.0,traits, num_concealed_states=4, idparslist,#here the matrix
          idparsopt, initparsopt, idparsfix, parsfix,
          cond=cond,root_state_weight = root_state_weight, tol = tol,
          sampling_fraction=samp_frac, 
          maxiter = 1000 * round((1.25)^length(idparsopt)), 
          use_fortran=TRUE,methode="ode45", optimmethod = "simplex", 
          num_cycles = 1, run_parallel=FALSE)

### Saving Data

save(CR_sp, CTD2_sp, CTD3_sp, CTD4_sp, ETD_sp, file = "reprod_data.rds")

### Comparing models

# Aicc function
aicc<-function(ll,k){round((2*(-round(ll,4))+2*k+(2*k*(k+1))/(Ntip(phy2.0)-k-1)),3)}

aicc_vec<-c(aicc(CR_sp$ML, 3), aicc(CTD2_sp$ML, 4), aicc(CTD3_sp$ML, 5), aicc(CTD4_sp$ML, 6), aicc(ETD_sp$ML, 6))

reprod_model_comp<-as.data.frame(cbind(c(CR_sp$ML, CTD2_sp$ML, CTD3_sp$ML, CTD4_sp$ML, ETD_sp$ML), c(3, 4, 5, 6, 6), aicc_vec, round(akaike.weights(aicc_vec)$deltaAIC, 3), round(akaike.weights(aicc_vec)$weights, 3)))

colnames(reprod_model_comp)<-c("log likelihood", "number of parameters", "AICc", "Delta AIC", "AIC weights")
rownames(reprod_model_comp)<-c("Constant rate model", "Concealt trait 2", "Conceal trait 3", "Conceal trait 4", "Focal trait")

### Saving Data

write.table(reprod_model_comp, "Path/to/your/exit.file", sep ="\t")
