#### phylogentic diveristy and tree stuff for the mass. pneumos

library(ggtree)
library(treeio)
library(ggplot2)

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
theme(legend.position = "none")#+
#xlim(0,1)

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

