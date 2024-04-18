#### phylogentic diveristy and tree stuff for the mass. pneumos

library(ggtree)
library(treeio)
library(ggplot2)

#harmonic numbers for scaling -- see wakeley eqn 3.23
harm <- function(x){
  thing <- rep(NA,(length(x)-1))
 for (i in 1:(length(x)-1)) {
   thing[i] <- sum(1/i)
 }
  return(sum(thing)*2)
}

#small weights 
# read in the iqtree files
tree.files <- list.files('smcog/', pattern = ".fasta.treefile")
tree.list.sm <- list()
for (i in 1:length(tree.files)) {
  file <- paste0("smcog/",tree.files[i])
  tree.list.sm[[i]] <- treeio::read.newick(file)
}
names(tree.list.sm) <- gsub(".fasta.treefile","", tree.files)

pd.cog.sm <-rep(NA,length(tree.list.sm))
# calculate phylo diversity
for (i in 1:length(tree.list.sm)){
 pd.cog.sm[i] <- sum(tree.list.sm[[i]]$edge.length)
 pd.cog.sm[i] <- pd.cog.sm[i]*harm(pneumo.seq.df.sm[[i]]$sequence)
}

#large weights 
# read in the iqtree files
tree.files <- list.files('lgcog/', pattern = ".fasta.treefile")
tree.list.lg <- list()
for (i in 1:length(tree.files)) {
  file <- paste0("lgcog/",tree.files[i])
  tree.list.lg[[i]] <- treeio::read.newick(file)
}
names(tree.list.lg) <- gsub(".fasta.treefile","", tree.files)

pd.cog.lg <- rep(NA,length(tree.list.lg))
# calculate phylo diversity
for (i in 1:length(tree.list.lg)){
  pd.cog.lg[i] <- sum(tree.list.lg[[i]]$edge.length)
  pd.cog.lg[i] <- pd.cog.lg[i]*harm(pneumo.seq.df.lg[[i]]$sequence)
}



### make figs #####
co <- c("#46ACC2","#4A1942","#935116","#F5B041","#5C8001","#D44D5C")

ggplot(rbind(data.frame(PD=pd.cog.lg,Weight=rep("Large")),
             data.frame(PD=pd.cog.sm,Weight=rep("Small"))),
       aes(x=PD, fill=Weight, color=Weight)) +
   #geom_density(alpha = 0.2) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.01) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Phylogenetic diversity", y = "Density", fill="Weight", color="Weight") +
  theme_minimal(base_size = 12) +
theme(legend.position = "none")+
xlim(0,10)

ggsave("figs/PD-den-MA.png",dpi=300,bg = "white",height=3,width = 3)

# there was one in the large group with really high div. ill look at the tree for those 
######## TREES #######
which.max(pd.cog.lg) #the super highly diverse one 
names(tree.list.lg)[which.max(pd.cog.lg)]
ggtree(tree.list.lg[["CLS01665"]])
#low div one for comparison 
names(tree.list.lg)[which.min(pd.cog.lg)]
ggtree(tree.list.lg[["CLS00602"]])


####### KEEP ONLY TIPS FOR PRE VAX SAMPLES ######

s.data <- read.table("s_study_Croucher.txt", header = T,sep ="\t" )

prevax.isolates <- dplyr::filter(s.data,Factor.Value.isolation.year. ==2001)
prevax.isolate.names <- prevax.isolates[,1]

prevax.tree.list.lg <- list()
for (i in 1:221){
  prevax.tree.list.lg[[i]] <- ape::keep.tip(tree.list.lg[[i]], which(str_detect(tree.list.lg[[i]]$tip.label, str_c(prevax.isolate.names, collapse ="|"))))
}
prevax.tree.list.sm <- list()
for (i in 1:221){
  prevax.tree.list.sm[[i]] <- ape::keep.tip(tree.list.sm[[i]], which(str_detect(tree.list.sm[[i]]$tip.label, str_c(prevax.isolate.names, collapse ="|"))))
}



#### pre vax PD #####
#large 
pvax.pd.cog.lg <- rep(NA,length(prevax.tree.list.lg))
for (i in 1:length(prevax.tree.list.lg)){
  pvax.pd.cog.lg[i] <- sum(prevax.tree.list.lg[[i]]$edge.length)
}
#small
pvax.pd.cog.sm <- rep(NA,length(prevax.tree.list.sm))
for (i in 1:length(prevax.tree.list.sm)){
  pvax.pd.cog.sm[i] <- sum(prevax.tree.list.sm[[i]]$edge.length)
}

#plot 
ggplot(rbind(data.frame(PD=pvax.pd.cog.lg,Weight=rep("Large")),
             data.frame(PD=pvax.pd.cog.sm,Weight=rep("Small"))),
       aes(x=PD, fill=Weight, color=Weight)) +
  #geom_density(alpha = 0.2) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.01) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Phylogenetic diversity", y = "Density", fill="Weight", color="Weight") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")#+
#xlim(0,1)

ggsave("figs/PD-den-MA.png",dpi=300,bg = "white",height=3,width = 3)



###### POST VAX ##########
# check how many isolates in each year 
library(dplyr)
s.data %>% count(Factor.Value.isolation.year.)
#names of the 2007 ones 
postvax.isolates <- dplyr::filter(s.data,Factor.Value.isolation.year. ==2007)
postvax.isolate.names <- postvax.isolates[,1]
#2004 names 
inter.isolates <- dplyr::filter(s.data,Factor.Value.isolation.year. ==2004)
inter.isolate.names <- inter.isolates[,1]

postvax.isolates <- dplyr::filter(s.data,Factor.Value.isolation.year. ==2007)
postvax.isolate.names <- postvax.isolates[,1]

postvax.tree.list.sm <- list()
for (i in 1:221){
  postvax.tree.list.sm[[i]] <- ape::keep.tip(tree.list.sm[[i]], which(str_detect(tree.list.sm[[i]]$tip.label, str_c(inter.isolate.names, collapse ="|"))))
}

postvax.tree.list.lg <- list()
for (i in 1:221){
  postvax.tree.list.lg[[i]] <- ape::keep.tip(tree.list.lg[[i]], which(str_detect(tree.list.lg[[i]]$tip.label, str_c(inter.isolate.names, collapse ="|"))))
}



#### post vax PD #####
#large 
postvax.pd.cog.lg <- rep(NA,length(postvax.tree.list.lg))
for (i in 1:length(postvax.tree.list.lg)){
  postvax.pd.cog.lg[i] <- sum(postvax.tree.list.lg[[i]]$edge.length)
}
#small
postvax.pd.cog.sm <- rep(NA,length(postvax.tree.list.sm))
for (i in 1:length(postvax.tree.list.sm)){
  postvax.pd.cog.sm[i] <- sum(postvax.tree.list.sm[[i]]$edge.length)
}

#plot 
ggplot(rbind(data.frame(PD=postvax.pd.cog.lg,Weight=rep("Large")),
             data.frame(PD=postvax.pd.cog.sm,Weight=rep("Small"))),
       aes(x=PD, fill=Weight, color=Weight)) +
  #geom_density(alpha = 0.2) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.01) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Phylogenetic diversity", y = "Density", fill="Weight", color="Weight") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")#+
#xlim(0,1)

ggsave("figs/PD-den-MA.png",dpi=300,bg = "white",height=3,width = 3)



