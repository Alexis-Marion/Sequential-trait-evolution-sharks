# SecSSE Reproduction

# I Data preparation

### Loading packages
library("ape")
library("secsse")
library("DDD")
library("doMC")
library("tidyverse")
library("parallel")
library("qgraph")
library("stringr")

# Loading data
df<-read.csv("Trait_syndrom_tab.tsv", sep="\t") # omit sep ="\t" for .csv files
phy1.0<-read.nexus("382sp_16C_20FC.tree")

# Optional step (cleaning data)
states<-df$Reproduction
names(states)=str_replace((rownames(df)), " ", "_")
states2<-states[!names(states) %in% setdiff(names(states), phy1.0$tip.label)]

# Computing sampling fractions
f<-c(
length(states2[states2=="O"])/length(df$Reproduction[df$Reproduction=="O"]),
length(states2[states2=="Oo"])/length(df$Reproduction[df$Reproduction=="Oo"]),
length(states2[states2=="PV"])/length(df$Reproduction[df$Reproduction=="PV"]),
length(states2[states2=="YV"])/length(df$Reproduction[df$Reproduction=="YV"]))

states3<-as.data.frame(cbind(names(states2), as.factor(states2)))
colnames(states3)<-c("species","states")
rownames(states3)<-NULL
traits <- sortingtraits(states3, phy)

n_states<-length(unique(na.omit(traits)))

states<-traits

# II SecSSE analyses

# Loads the models to be tested
source("aux_scripts/secsse_functions.R")
source("aux_scripts/secsse_base_models.R")

# Constraining some transitions to be null
mask<-matrix(1, nrow=4, ncol=4)
diag(mask)<-NA
mask[1,3]<-0
mask[3,1]<-0
mask[1,4]<-0
mask[4,1]<-0

# Applying the previous matrix
for(i in 1:length(models)){
  models[[i]]$idparslist$Q<-mask_q(q=models[[i]]$idparslist$Q, mask=mask, n_states = n_states)
}

for(i in 1:length(models)){
  
  idparsopt = sort(unique(na.omit(unname(unlist(models[[i]]$idparslist)))))
  if(any(idparsopt %in% 0)){
    idparsopt = idparsopt[idparsopt!=0]
    idparsfix = c(0)
    parsfix = c(0)
  } else {
    idparsfix = c()
    parsfix = c()
  }
  
  models[[i]]<-append(models[[i]],
                      list(
                        idparsopt=idparsopt,
                        idparsfix=idparsfix,
                        parsfix=parsfix
                      )
  )
}

# Assigning starting values and initialization
startingpoint <- bd_ML(brts = ape::branching.times(phy))
intGuessLambda <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0

intGuessLambdas<-c(intGuessLambda, intGuessLambda/2, intGuessLambda*2)
intGuessMus <- c(intGuessMu, intGuessMu/2, intGuessMu*2)
initTrans <- intGuessLambdas/n_states

nrep=length(intGuessLambdas)

replicate_models<-rep(models,each=nrep)
names(replicate_models)<-paste0(rep(names(models), each=nrep), rep(paste0("_try", 1:nrep),nrep))
names(replicate_models)
models<-replicate_models
iterator=rep(1:nrep, length(models)/nrep)
for(i in 1:length(models)){
  initparsopt = c(
    rep(intGuessLambdas[iterator[i]],length(unique(models[[i]]$idparslist$lambdas))),
    ifelse(unique(models[[i]]$idparslist$mus)!=0,
      rep(intGuessMus[iterator[i]],length(unique(models[[i]]$idparslist$mus))),
      0
    ),
    rep(initTrans[iterator[i]],length(which(unique(c(models[[i]]$idparslist$Q))>0)))
  )
  initparsopt<-initparsopt[initparsopt!=0]

  models[[i]]$initparsopt<-initparsopt

}

all(sapply(models, FUN = function(x) length(x$initparsopt)==length(x$idparsopt)))

lapply(models, `[[`, "idparslist")

models<-lapply(models, FUN=append,
               list(phy = phy,
                    traits = states,
                    cond="proper_cond",
                    root_state_weight = c(1,0,0,0),
                    sampling_fraction=f,
                    num_cycles = 3,
                    run_parallel=T)
               )

# Assigning starting value and initialization
run_secsse<-function(i) {

  if(!dir.exists("secsse_out_Reproduction")) dir.create("secsse_out_Reproduction")

  try(
    saveRDS(do.call(secsse_ml, args=models[[i]]),
            file = paste0("secsse_out_Reproduction/",names(models)[i],"_out.rds")))
}

# Running the models in parallel
n_cores=ifelse(length(models)>30, 30, length(models))
n_cores
mclapply(FUN=run_secsse,
         X=1:length(models),
         mc.cores = n_cores)
