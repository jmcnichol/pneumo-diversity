
require(stringr)
require(seqinr)
require(ape)
library(dplyr)
library(reshape2)
# Extracting the Mass. isolates from the cogs based on the 616 in the nick paper
# supp data from Nick
s.data <- read.table("s_study_Croucher.txt", header = T,sep ="\t" )

sum(str_detect(pneumo.seq.df.lg[["CLS02117"]][["name"]], str_c(s.data$V1, collapse ="|")))

# extract the Mass isolates
for (i in 1:length(pneumo.seq.df.lg)){
  pneumo.seq.df.lg[[i]] <- pneumo.seq.df.lg[[i]][which(str_detect(pneumo.seq.df.lg[[i]][["name"]], str_c(s.data$V1, collapse ="|"))),]
}
gsub("\\_.*","",pneumo.seq.df.lg[[1]]$name) #they are in order omg
for (i in 1:length(pneumo.seq.df.sm)){
  pneumo.seq.df.sm[[i]] <- pneumo.seq.df.sm[[i]][which(str_detect(pneumo.seq.df.sm[[i]][["name"]], str_c(s.data$V1, collapse ="|"))),]
}

# 
weight.small <- read.table("gCOG_sequences/sm_weight_file_names.txt")
weight.small <- apply(weight.small,2,function(x)gsub(".aligned.fasta","",x))

weight.large <- read.table("gCOG_sequences/lg_weight_file_names.txt")
weight.large <- apply(weight.large,2,function(x)gsub(".aligned.fasta","",x))

# read in the aligned sequences with small weights 
raw.files <- list.files('gCOG_sequences/SmallWeights/')
filenames <- paste0(weight.small,".aligned.fasta")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df.sm <- list()
for (i in 1:length(files.all)) {
  file <- paste0("gCOG_sequences/SmallWeights/",files.all[i])
  file.list[[i]] <- seqinr::read.fasta(file,seqtype= "AA", as.string = T, forceDNAtolower = F)
  
  pneumo.seq.df.sm[[i]] <- data.frame(name=paste(getAnnot(file.list[[i]])), sequence=paste0(file.list[[i]]))
  pneumo.seq.df.sm[[i]]$name <- gsub(">", "", pneumo.seq.df.sm[[i]]$name)
}
names(pneumo.seq.df.sm) <- gsub(".aligned.fasta","", files.all)

# read in large
raw.files <- list.files('gCOG_sequences/LargeWeights/')
filenames <- paste0(weight.large,".aligned.fasta")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df.lg <- list()
for (i in 1:length(files.all)) {
  file <- paste0("gCOG_sequences/LargeWeights/",files.all[i])
  file.list[[i]] <- seqinr::read.fasta(file,seqtype= "AA", as.string = T, forceDNAtolower = F)
  
  pneumo.seq.df.lg[[i]] <- data.frame(name=paste(getAnnot(file.list[[i]])), sequence=paste0(file.list[[i]]))
  pneumo.seq.df.lg[[i]]$name <- gsub(">", "", pneumo.seq.df.lg[[i]]$name)
}
names(pneumo.seq.df.lg) <- gsub(".aligned.fasta","", files.all)


#### colours for plots ####

co <- c("#46ACC2","#4A1942","#935116","#F5B041","#5C8001","#D44D5C")


######## TAJIMA'S D ##########

# create a list of matrices instead of char strings
lg.mat <- list()
for (c in 1:221) {
  mat1 <- matrix(NA, nrow=length(pneumo.seq.df.lg[[c]][["sequence"]]),ncol=str_length(pneumo.seq.df.lg[[c]][["sequence"]][1]))
  for (m in 1:length(pneumo.seq.df.lg[[c]][["sequence"]])) {
    mat1[m,] <- s2c(pneumo.seq.df.lg[[c]][["sequence"]][m])
  }
  lg.mat[[c]] <- mat1
}
sm.mat <- list()
for (c in 1:684) {
  mat1 <- matrix(NA, nrow=length(pneumo.seq.df.sm[[c]][["sequence"]]),ncol=str_length(pneumo.seq.df.sm[[c]][["sequence"]][1]))
  for (m in 1:length(pneumo.seq.df.sm[[c]][["sequence"]])) {
    mat1[m,] <- s2c(pneumo.seq.df.sm[[c]][["sequence"]][m])
  }
  sm.mat[[c]] <- mat1
}

tp.lg = rep(NA,length(pneumo.seq.df.lg))
tp.sm = rep(NA,length(pneumo.seq.df.sm))
for (l in 1:length(lg.mat)) {
  #convert character strings to binary DNA
  #returns indices of segregating sites -- i just want to know how many there are, so add em up
  # td.lg[l] <- pegas::tajima.test(as.DNAbin(lg.mat[[l]]))$D
  #  td.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$D
  tp.lg[l] <- pegas::tajima.test(as.DNAbin(lg.mat[[l]]))$Pval.beta
  # tp.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$Pval.normal
}

for (l in 1:length(sm.mat)) {
  #convert character strings to binary DNA
  # td.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$D
  tp.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$Pval.beta
}

# look at pvalues
pdata <- c(psmall=tp.sm, plarge=tp.lg)
pdata <- data.frame(Weight=c(rep("Small",length(td.sm)),rep("Large",length(td.lg))),pval=pdata)

pdata <- data.frame(psmall=tp.sm, plarge=tp.lg)
sigps <- melt(pdata) %>% group_by(variable) %>% count(sig=value<=0.05) %>% group_by(variable) 
ggplot(filter(sigps, is.na(sig)==F), aes(x=n, y=sig, fill=variable)) +
  geom_bar(stat="identity",position = position_dodge(),alpha = 0.4) +
  #geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(labels = c("Small","Large"), values=c(co[2],co[1]))+
  scale_color_manual(labels = c("Small","Large"), values=c(co[2],co[1]))+
  labs(x="Number of p-values", y = " ", fill="Weight") +
  theme_minimal(base_size = 12) +
  scale_y_discrete(labels=c("TRUE" = "p < 0.05", "FALSE" = "p > 0.05"))

# D values 
library(reshape2)
ddata <- c(psmall=td.sm, plarge=td.lg)
ddata <- data.frame(Weight=c(rep("Small",length(td.sm)),rep("Large",length(td.lg))),D=ddata)

Ddata <- data.frame(Large=td.lg,Small=td.sm)
ggplot(filter(melt(Ddata), is.na(value)==F), aes(x=value, fill=variable, color=variable)) +
 # geom_density(alpha = 0.2) +
  geom_histogram( alpha=0.4, position = "identity", binwidth = 0.1) +
  scale_fill_manual( values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Tajima's D", y = "Count", fill="Weight",color="Weight") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
ggsave("figs/D-hist.png",dpi=300,bg = "white",height=4,width = 6)

#overlay the density
ggplot(filter(melt(Ddata), is.na(value)==F), aes(x=value, fill=variable, color=variable)) +
    geom_density(alpha = 0.2) +
 # geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  scale_fill_manual( values=c(co[1],co[2]))+
  scale_color_manual( values=c(co[1],co[2]))+
  labs(x="Tajima's D", y = "Density", fill="Weight",color="Weight") +
  theme_minimal(base_size = 12)+
  theme(legend.position = "none")
ggsave("figs/D-den.png",dpi=300,bg = "white",height=3,width = 3)

##### SHANNON #####
#calculates Shannon entropy 
shan.div <- function(list.of.dfs){
  div <- rep(0,length(list.of.dfs))
  for (i in 1:length(list.of.dfs)){
    # unique sequences
    pneumo.freq <- list.of.dfs[[i]] %>% dplyr::count(sequence)
    pneumo.freq$n <- pneumo.freq$n/nrow(list.of.dfs[[i]])
    div[i] <- vegan::diversity(pneumo.freq$n)
  }
  div <- data.frame(locus = names(list.of.dfs),diversity = div)
  return(div)
}


pneumo.shan.lg <- shan.div(pneumo.seq.df.lg)
pneumo.shan.sm <- shan.div(pneumo.seq.df.sm)

#combine large and small for one plot
shan.comb <- rbind(data.frame(pneumo.shan.sm, weight= rep("Small")),
                   data.frame(pneumo.shan.lg, weight=rep("Large")))

ggplot(shan.comb, aes(x=diversity, fill=weight, color=weight)) +
 # geom_density(alpha = 0.2) +
  geom_histogram( alpha=0.4, position = "identity", binwidth = 0.1) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Shannon diversity", y = "Count", fill="Weight", color="Weight") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("figs/shannon-hist.png",dpi=300,bg = "white",height=4,width = 6)


ggplot(shan.comb, aes(x=diversity, fill=weight, color=weight)) +
   geom_density(alpha = 0.2) +
  #geom_histogram( alpha=0.4, position = "identity", binwidth = 0.1) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Shannon diversity", y = "Diversity", fill="Weight", color="Weight") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
ggsave("figs/shannon-div.png",dpi=300,bg = "white",height=3,width = 3)

### sort by VT/NVT ####

#im kind of confused. i can group the sequences by VT or time but if i calc diversity im calculating diversity of the cogs
#should i calc diversity of the taxa instead? 

## if i want to calculate PD should i export the fasta files? 
# are there loci that dont have weights at all? 


for (i in 1:length(pneumo.seq.df.lg)) {
  seqinr::write.fasta(names=pneumo.seq.df.lg[[i]]$name,sequences=as.list(pneumo.seq.df.lg[[i]]$sequence),file.out=paste0("lgcog/",names(pneumo.seq.df.lg)[i],".fasta"))
  seqinr::write.fasta(names=pneumo.seq.df.sm[[i]]$name,sequences=as.list(pneumo.seq.df.sm[[i]]$sequence),file.out=paste0("smcog/",names(pneumo.seq.df.sm)[i],".fasta"))
  
}


######## TEST ########
lg.mean <- mean((td.lg),na.rm = T)
sm.mean <- mean((td.sm),na.rm = T)

sd((td.lg),na.rm = T)

abs(lg.mean-sm.mean)/sqrt(1/221+1/684)
dnorm(abs(lg.mean-sm.mean)/sqrt(1/221+1/684))

tp.lg.p <- p.adjust(tp.lg)
tp.sm.p <- p.adjust(tp.sm)

sum(tp.lg.p < 0.05, na.rm = T)/221*100
sum(tp.sm.p < 0.05, na.rm = T)/684*100
#use these ones and see if they do as well. they cant use all the weights. make model of NFDS using this approach. 


#### summary stats and tests ########
shan.comb %>% group_by(weight) %>% summarise(sd(diversity))
ks.test(x=pneumo.shan.lg$diversity,y=pneumo.shan.sm$diversity, exact = T)

#from tree-diversity.R 
mean(pd.cog.lg)
mean(pd.cog.sm)

sd(pd.cog.lg)
sd(pd.cog.sm)

ks.test(x=pd.cog.lg[-62],y=pd.cog.sm)

which.max(pd.cog.lg)
pd.cog.lg[62]

sd(pd.cog.lg[-62])
mean(pd.cog.lg[-62])
