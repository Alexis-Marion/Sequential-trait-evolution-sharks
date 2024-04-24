library("ape")
library("ggtree")
library("stringr")
library("phytools")
library("mclust")
library("treedataverse")
library("RColorBrewer")
library("ggplot2")
library("ggpmisc")

df<-read.csv("Species_data.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.nexus("16FC_16C_374_sp.tree")
df_diet<-read.table("table_data_diet.tsv", header = TRUE)
colnames(df_diet)<-c("Species", "Diet")
df_diet[,2]<-as.character(df_diet[,2])

mb1 = Mclust(as.numeric(df$Body.size))
df$Body.size <- as.factor(mb1$classification)

df$Species<-str_replace(df$Species, " ", "_")

states<-cbind(df$Species, df$Body.size, df$Reproduction, df$Habitat)

states_traits<-states[!states[,1] %in% setdiff(states[,1], phy$tip.label),]

states_traits[is.na(states_traits)]<-"?"

states_traits<-states_traits[match(phy$tip.label,states_traits[,1]),]

phylo<-ggtree(phy, layout="fan", open.angle=180)  + 
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

# For corHMM

node_states_size<-readRDS("ASE/bds_corHMM.rds")$states 

node_states<-c(375:747)

col1 <- (node_states_size[,1]+node_states_size[,4])
col2 <- (node_states_size[,2]+node_states_size[,5])
col3 <- (node_states_size[,3]+node_states_size[,6])

# For SecSSE

#node_states_size<-readRDS("ASE/bds_SecSSE.rds")

#node_states_size<-node_states_size[375:747,c(10:18)]

#col1 <- (node_states_size[,1]+node_states_size[,4]+node_states_size[,7])
#col2 <- (node_states_size[,2]+node_states_size[,5]+node_states_size[,8])
#col3 <- (node_states_size[,3]+node_states_size[,6]+node_states_size[,9])

node_states_size<-cbind(col1, col2, col3, node_states)

node_states_size<-as.data.frame(node_states_size)

colnames(node_states_size)<-c("Small", "Medium", "Large", "node")

pies <- nodepie(node_states_size, cols=1:3, color=c("#A2EBD4", "#D5E3DB", "#74BDCB"), alpha=1)

df<-tibble::tibble(node=as.numeric(node_states_size$node), pies=pies)

phylo_node_bds <- phylo %<+% df

states_traits<-as.data.frame(states_traits)

colnames(states_traits)<-c("Species", "Size", "Reproduction", "Habitat")

phylo_bds_complete<-phylo_node_bds + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.0125, vp.height=0.0125, hjust=0.5, vjust=0.5)

phylo_bds_complete<-phylo_bds_complete %<+% states_traits[,c(1,2)]
ASE_plot<-phylo_bds_complete + geom_tippoint(data=td_filter(isTip),aes(color=Size), size=0.4) + scale_color_manual(values=c("#74BDCB", "#D5E3DB", "#A2EBD4"))

ggsave(ASE_plot, filename = "Ase_body_size.pdf",  bg = "transparent", width = 10, height = 10)

# For corHMM

node_states_reproduction<-readRDS("ASE/rp_corHMM.rds")

node_states_reproduction<-cbind(node_states_reproduction$states, node_states)

col1 <- (node_states_reproduction[,1]+node_states_reproduction[,5])
col2 <- (node_states_reproduction[,2]+node_states_reproduction[,6])
col3 <- (node_states_reproduction[,3]+node_states_reproduction[,7])
col4 <- (node_states_reproduction[,4]+node_states_reproduction[,8])

# For SecSSE

#node_states_reproduction<-readRDS("ASE/rp_SecSSE.rds")

#node_states_reproduction<-node_states_reproduction[375:747,c(17:32)]

#col1 <- (node_states_reproduction[,1]+node_states_reproduction[,5]+node_states_reproduction[,9]+node_states_reproduction[,13])
#col2 <- (node_states_reproduction[,2]+node_states_reproduction[,6]+node_states_reproduction[,10]+node_states_reproduction[,14])
#col3 <- (node_states_reproduction[,3]+node_states_reproduction[,7]+node_states_reproduction[,11]+node_states_reproduction[,15])
#col4 <- (node_states_reproduction[,4]+node_states_reproduction[,8]+node_states_reproduction[,12]+node_states_reproduction[,16])

node_states_reproduction<-cbind(col1, col2, col3, col4, node_states)

colnames(node_states_reproduction)<-c("Oviparous", "Oophageous", "Placental", "Yolk_sack_dependant", "node")

node_states_reproduction<-as.data.frame(node_states_reproduction)

pies <- nodepie(node_states_reproduction, cols=1:4, color=c("#C91036", "#F4EBBC", "#FF918B", "#FF775E"), alpha=1)

df<-tibble::tibble(node=as.numeric(node_states_reproduction$node), pies=pies)

phylo_node_rp <- phylo %<+% df

phylo_rp_complete<-phylo_node_rp + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.0125, vp.height=0.0125, hjust=0.5, vjust=0.5) 

phylo_rp_complete<-phylo_rp_complete %<+% states_traits[,c(1,3)]
ASE_plot<-phylo_rp_complete + geom_tippoint(data=td_filter(isTip),aes(color=Reproduction), size=0.4) + scale_color_manual(values=c("#E5E5E5", "#F4EBBC", "#C91036", "#FF918B", "#FF775E"))

ggsave(ASE_plot, filename = "Ase_reproduction.pdf",  bg = "transparent", width = 10, height = 10)

# For corHMM

node_states_habitat<-readRDS("ASE/habitat_corHMM.rds")

node_states_habitat<-cbind(node_states_habitat$states, node_states)

# For SecSSE

#node_states_habitat<-readRDS("ASE/habitat_SecSSE.rds")

#node_states_habitat<-node_states_habitat[375:747,c(26:50)]

#col1 <- (node_states_habitat[,1]+node_states_habitat[,6]+node_states_habitat[,11]+node_states_habitat[,16]+node_states_habitat[,21])
#col2 <- (node_states_habitat[,2]+node_states_habitat[,7]+node_states_habitat[,12]+node_states_habitat[,17]+node_states_habitat[,22])
#col3 <- (node_states_habitat[,3]+node_states_habitat[,8]+node_states_habitat[,13]+node_states_habitat[,18]+node_states_habitat[,23])
#col4 <- (node_states_habitat[,4]+node_states_habitat[,9]+node_states_habitat[,14]+node_states_habitat[,19]+node_states_habitat[,24])
#col5 <- (node_states_habitat[,5]+node_states_habitat[,10]+node_states_habitat[,15]+node_states_habitat[,20]+node_states_habitat[,25])

colnames(node_states_habitat)<-c("Deepwater", "Inner_Shelf", "Oceanic", "Outer_shelf", "Reef", "node")

node_states_habitat<-as.data.frame(node_states_habitat)

pies <- nodepie(node_states_habitat, cols=1:5, color=c("#5E4FA2", "#FEE08B", "#BEDCEB", "#74BDCB", "#FFA384"), alpha=1)

df<-tibble::tibble(node=as.numeric(node_states_habitat$node), pies=pies)

phylo_node_ht <- phylo %<+% df

phylo_ht_complete<-phylo_node_ht + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.0125, vp.height=0.0125, hjust=0.5, vjust=0.5)

phylo_ht_complete<-phylo_ht_complete %<+% as.data.frame(states_traits[,c(1,4)])
ASE_plot<-phylo_ht_complete + geom_tippoint(data=td_filter(isTip),aes(color=Habitat), size=0.4) + scale_color_manual(values=c("#5E4FA2", "#FEE08B", "#BEDCEB", "#74BDCB", "#FFA384"))

ggsave(ASE_plot, filename = "Ase_habitat.pdf",  bg = "transparent", width = 10, height = 10)

# For corHMM

node_states_diet<-readRDS("ASE/diet_corHMM.rds")

node_states_diet<-cbind(node_states_diet$states, node_states)

df_diet[,2]<-as.character(df_diet[,2])

col1 <- (node_states_diet[,1]+node_states_diet[,4])
col2 <- (node_states_diet[,2]+node_states_diet[,5])
col3 <- (node_states_diet[,3]+node_states_diet[,6])

# For SecSSE

node_states_diet<-readRDS("ASE/diet_SecSSE.rds")

node_states_diet<-node_states_diet[375:747,c(10:18)]

col1 <- (node_states_diet[,1]+node_states_diet[,4]+node_states_diet[,7])
col2 <- (node_states_diet[,2]+node_states_diet[,5]+node_states_diet[,8])
col3 <- (node_states_diet[,3]+node_states_diet[,6]+node_states_diet[,9])

node_states_diet<-cbind(col1, col2, col3, node_states)

node_states_diet<-as.data.frame(node_states_diet)

colnames(node_states_diet)<-c("Invertebrate-feeder", "Macropredator", "Mesopredator", "node")

pies <- nodepie(node_states_diet, cols=1:3, color=c("#52B2CF", "#D4AFB9", "#9CADCE"), alpha=1)

df<-tibble::tibble(node=as.numeric(node_states_diet$node), pies=pies)

phylo_node_diet <- phylo %<+% df

phylo_diet_complete<-phylo_node_diet + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.0125, vp.height=0.0125, hjust=0.5, vjust=0.5) 

phylo_diet_complete<-phylo_diet_complete %<+% df_diet[,c(1,2)]
ASE_plot<-phylo_diet_complete + geom_tippoint(data=td_filter(isTip),aes(color=Diet), size=0.4) + scale_color_manual(values=c("#52B2CF", "#D4AFB9", "#9CADCE"))

ggsave(ASE_plot, filename = "Ase_diet.pdf",  bg = "transparent", width = 10, height = 10)
