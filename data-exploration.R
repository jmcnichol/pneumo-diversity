library(dplyr)
library(seqinr)
library(parallel)
library(vegan)
library(ggplot2)
library(ggpubr)
library(insect) #for char2dna function 
library(ape)
# gene has to do something so could be less diverse. 
# pushes against the notion that this is driven by nfds 

### Example/test case ####
# have a look at an example file
pneumo.test <- seqinr::read.fasta("gCOG_sequences/LargeWeights/CLS01244.aligned.fasta", seqtype= "AA", as.string = T, forceDNAtolower = F)
# a lot of these are apparently going to be the same 
# some like 40% of one string and 60 of another
# compute plog where there a term for each of the unique sequences. 
# frequency of each unique one then compute the entropy for them. 
# but that may not be a great idea if there is low diversity -- so first i need to 
# see how unique the rows are then see what diversity to use. 
# -- surprising if they're all unique -- don't use Shannon in that case. 
# go through all them and compare unique sequences to taxa. 

# make a ML phylogeny instead? 

pneumo.test <- data.frame(name=paste(getAnnot(pneumo.test)), sequence=paste0(pneumo.test))
pneumo.test$name <- gsub(">", "", pneumo.test$name)
seq.unique <- unique(pneumo.test$sequence)
nrow(pneumo.test) #4243
length(seq.unique) #454
#freq of each unique thing 
pneumo.freq <- pneumo.test %>% dplyr::count(sequence)
pneumo.freq$n <- pneumo.freq$n/nrow(pneumo.test)
#calculate Shannon entropy 
div <- vegan::diversity(pneumo.freq$n)

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



##### ------------------------------------------------------------#
### Explore the data ####
# read in data 
massdata <- read.csv("massdata.csv")
# create df for all the stuff that aren't loci id cols 
mass.info <- massdata[,2:6]
select.strength <- read.table("intLoci_selection_strengths.csv")
colnames(select.strength) <- c("locus","weight")
hist(select.strength$weight) #two distinct groups  
unique(select.strength$weight) # there are only two weights 

# read in the ones that we have weights for
raw.files <- list.files('gCOG_sequences2/')
filenames <- paste0(select.strength$locus ,".out")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df <- list()
for (i in 1:length(files.all)) {
  file <- paste0("gCOG_sequences2/",files.all[i])
  file.list[[i]] <- seqinr::read.fasta(file,seqtype= "AA", as.string = T, forceDNAtolower = F)
  
  pneumo.seq.df[[i]] <- data.frame(name=paste(getAnnot(file.list[[i]])), sequence=paste0(file.list[[i]]))
  pneumo.seq.df[[i]]$name <- gsub(">", "", pneumo.seq.df[[i]]$name)
}
#get the list of loci that were read in with out the file extension 
loci.weigted <- gsub(".out", "", files.all)
names(pneumo.seq.df) <- loci.weigted #905 of 1090 files matched with weights
#now we have a df for each one 

### data manipulation ####
# choose the cols from massdata that match the weights and were read in 
weight.large <- select.strength |> filter(weight > 0.06) |> filter(locus %in% loci.weigted)
weight.small <- select.strength |> filter(weight < 0.06) |> filter(locus %in% loci.weigted)
##### New direction -- smaller sample size #####
# turns out the above sequences weren't aligned. we will sample so that each 
# high/low weight have the same number then align those files and delete some 
# of the loci (?) so that the number of tips is reduced 
# (takes a long time to generate a tree)
dim(weight.large) #there are 221 large ones, lets sample 221 small ones
#set.seed(12345)
#weight.sample <- sample.int(221)
#weight.small <- weight.small[-weight.sample,]

nrow(weight.large) #270 --> 221 read in 
nrow(weight.small) #820 --> 684 read in 

# extract the corresponding cols from massdata 
# big weights (keep only the ones that were read in)
data.large.weight <- massdata |> select(weight.large$locus)
data.large.weight <- cbind(mass.info,data.large.weight)
#small weights
data.small.weight <- massdata |> select(weight.small$locus)
data.small.weight <- cbind(mass.info,data.small.weight)

### Apply shannon entropy ####
# all data
pneumo.shan <- shan.div(pneumo.seq.df)
# split into small and large weight groups 
large.weight.div <- pneumo.shan |> filter(locus %in% weight.large$locus)
small.weight.div <- pneumo.shan |> filter(locus %in% weight.small$locus)
# make a nice df 
div.data <- rbind(data.frame(large.weight.div, weight = rep(as.factor(0.1363),nrow(large.weight.div))),
      data.frame(small.weight.div, weight = rep(as.factor(0.0023),nrow(small.weight.div))))

### Plots ####
ggplot(div.data, aes(x=diversity, fill=weight)) +
 geom_density(alpha = 0.4) +
  #geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
 # geom_density(alpha=0,aes(color=vax)) +
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

ggplot(div.data, aes(x=diversity, fill=weight)) +
  #geom_density(alpha = 0.4) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

# they look similar lets do a KS test
ks.test(large.weight.div$diversity,small.weight.div$diversity)
#p = 0.4981 --- in the sample the pvalue is 0.118 (but these are not aligned)

mean(large.weight.div$diversity)
mean(small.weight.div$diversity)

var(large.weight.div$diversity)
var(small.weight.div$diversity)

##### save selected file names so i can align only those ones. 
# get rid of taxa with few characters 

# max number of chars for the file, or um delete everything less than 2/3 of that
# most or any of taxa are full length sequences (?)
# remove seq that are "too short" <3/4 of predom length, for some gene that's 1500 bases.

#### find out how many bases each taxon one has ####
seq.length <- list()
for (k in 1:length(pneumo.seq.df)) {
  seq.length[[k]] <- nchar(pneumo.seq.df[[k]]$sequence)
  
}
lapply(seq.length, range)

median(seq.length[[900]])

#### Align ###### 
# start with the high weight ones, get their file names for alignment 
lg_weight_file_names <- paste0(weight.large$locus,".fasta")
write.table(lg_weight_file_names, "lg_weight_file_names.txt", sep="\t", col.names = F, row.names = F)

# shell script to align files--- but putting it here to remeber 
#for i in *.fasta
#do    
#mafft $i > ${i%.fasta}.aligned.fasta
#done

# start with the high weight ones, get their file names for alignment 
sm_weight_file_names <- paste0(weight.small$locus,".fasta")
write.table(sm_weight_file_names, "sm_weight_file_names.txt", sep="\t", col.names = F, row.names = F)

# do this again for the rest of the small weight ones -- may.7 2024
sm_weight_file_names2 <- paste0(weight.small$locus,".fasta")
write.table(sm_weight_file_names2, "sm_weight_file_names2.txt", sep="\t", col.names = F, row.names = F)

########### Aligned Sequences #########
#### replicating the original script for the aligned sequences 

# read in the aligned sequences with large weights 
raw.files <- list.files('gCOG_sequences/LargeWeights/')
filenames <- paste0(weight.large$locus,".aligned.fasta")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df.lg <- list()
for (i in 1:length(files.all)) {
  file <- paste0("gCOG_sequences/LargeWeights/",files.all[i])
  file.list[[i]] <- seqinr::read.fasta(file,seqtype= "AA", as.string = T, forceDNAtolower = F)
  
  pneumo.seq.df.lg[[i]] <- data.frame(name=paste(getAnnot(file.list[[i]])), sequence=paste0(file.list[[i]]))
  pneumo.seq.df.lg[[i]]$name <- gsub(">", "", pneumo.seq.df.lg[[i]]$name)
}
names(pneumo.seq.df.lg) <- gsub(".aligned.fasta","", files.all)

pneumo.shan.lg <- shan.div(pneumo.seq.df.lg)

ggplot(pneumo.shan.lg, aes(x=diversity)) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

# read in the aligned sequences with small weights 
raw.files <- list.files('SmallCOGS/')
filenames <- paste0(weight.small$locus,".aligned.fasta")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df.sm <- list()
for (i in 1:length(files.all)) {
  file <- paste0("SmallCOGS/",files.all[i])
  file.list[[i]] <- seqinr::read.fasta(file,seqtype= "AA", as.string = T, forceDNAtolower = F)
  
  pneumo.seq.df.sm[[i]] <- data.frame(name=paste(getAnnot(file.list[[i]])), sequence=paste0(file.list[[i]]))
  pneumo.seq.df.sm[[i]]$name <- gsub(">", "", pneumo.seq.df.sm[[i]]$name)
}
names(pneumo.seq.df.sm) <- gsub(".aligned.fasta","", files.all)

pneumo.shan.sm <- shan.div(pneumo.seq.df.sm)

ggplot(pneumo.shan.sm, aes(x=diversity)) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

#combine large and small for one plot
shan.comb <- rbind(data.frame(pneumo.shan.sm, weight= rep("Small")),
      data.frame(pneumo.shan.lg, weight=rep("Large")))

colos <- c("#754668","#935116","#5DADE2","#F5B041","#5C8001","#D44D5C")

co <- c("#46ACC2","#4A1942")


ggplot(shan.comb, aes(x=diversity, fill=weight)) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  labs(x="Shannon diversity", y = "Density") +
  scale_fill_manual(values=c(colos[1],colos[3]))+
  theme_minimal() 

shan.den <- ggplot(shan.comb, aes(x=diversity, fill=weight, color=weight)) +
  geom_density(alpha = 0.4) +
  labs(x="Shannon diversity", y = "Density") +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  theme_minimal() +
  theme(legend.position = "bottom")

ks.test(pneumo.shan.lg$diversity,pneumo.shan.sm$diversity)
#p-value = 0.1182

median(pneumo.shan.lg$diversity)


#### Segregating Sites ####
# number of segregating sites -- number of col with at least two different letters 
num.seg.lg = num.seg.sm = rep(NA,length(pneumo.seq.df.lg))
for (l in 1:length(pneumo.seq.df.lg)) {
  #convert character strings to binary DNA
  #returns indices of segregating sites -- i just want to know how many there are, so add em up
  num.seg.lg[l] <- length(ape::seg.sites(char2dna(pneumo.seq.df.lg[[l]]$sequence)))/stringr::str_length((pneumo.seq.df.lg[[l]]$sequence[1]))
  num.seg.sm[l] <- length(ape::seg.sites(char2dna(pneumo.seq.df.sm[[l]]$sequence)))/stringr::str_length((pneumo.seq.df.sm[[l]]$sequence[1]))
}

#combine large and small for one plot
segsites.comb <- rbind(data.frame(seg = num.seg.sm, weight= rep("Small")),
                   data.frame(seg = num.seg.lg, weight=rep("Large")))

seg.hist <- ggplot(segsites.comb, aes(x=seg, fill=weight)) +
  geom_histogram( alpha=0.5, position = "identity") +
  labs(x="Proportion of segregating sites", y = "Count") +
  scale_fill_manual(values=c(co[1],co[2])) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom")
# what may be more useful is something relative to the length of the sequences. 

library(ggpubr)

ggarrange(shan.den,seg.hist,legend = "bottom",common.legend = T)


#### Simpler version ####
# calc all the same stuff with 0/1 data instead 

binary.lg <- dplyr::select(massdata, weight.large$locus)
binary.sm <- dplyr::select(massdata, weight.small$locus)

#shannon diversity for each COG  
div.lg <- vegan::diversity(t(binary.lg)) 
div.sm <- vegan::diversity(t(binary.sm)) 

library(reshape2)
div2.lg <- data.frame(melt(div.lg),weight=rep("Large"))
div2.sm <- data.frame(melt(div.sm),weight=rep("Small"))
div2 <- rbind(div2.lg,div2.sm)


ggplot(div2, aes(x=value, fill=weight, color=weight)) +
  geom_density(alpha = 0.2) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="Shannon diversity", y = "Density", fill="Weight", color="Weight") +
  theme_minimal() 

ggplot(div2, aes(x=value, fill=weight)) +
  #geom_density(alpha = 0.4) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
   geom_density(alpha=0,aes(color=weight)) +
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

ks.test(div.lg,div.sm) # 0.07

#only include rows that were at time 0 
# then the last time 
# not vax type 

#### Group by time ####

#pre vax 
prevax <- select(dplyr::filter(massdata, Time=="Y0"), c(weight.small$locus,weight.large$locus))
postvax <- select(dplyr::filter(massdata, Time=="Y6"), c(weight.small$locus,weight.large$locus))

#shannon diversity for each COG  
pre.div <- vegan::diversity(t(prevax)) 
post.div <- vegan::diversity(t(postvax)) 

ks.test(pre.div[,1],post.div[,1])

pre.div <- data.frame(melt(pre.div),Time=rep("Pre-vax"))
post.div <- data.frame(melt(post.div),Time=rep("Post-vax"))
div3 <- rbind(pre.div,post.div)

ggplot(div3, aes(x=value, fill=Time)) +
  #geom_density(alpha = 0.4) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  geom_density(alpha=0,aes(color=Time)) +
  scale_fill_manual(values=c(colos[2],colos[4]))+
  scale_color_manual(values=c(colos[2],colos[4]))+
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

ggplot(div3, aes(x=value, fill=Time, color=Time)) +
  geom_density(alpha = 0.2) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(values=c(colos[2],colos[5]))+
  scale_color_manual(values=c(colos[2],colos[5]))+
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 

### Group by VT vs NVT ###
VT <- select(dplyr::filter(massdata, VT=="VT"), c(weight.small$locus,weight.large$locus))
NVT <- select(dplyr::filter(massdata, VT=="NVT"), c(weight.small$locus,weight.large$locus))

#shannon diversity for each COG  
VT.div <- vegan::diversity(t(VT)) 
NVT.div <- vegan::diversity(t(NVT)) 

ks.test(VT.div[,1],NVT.div[,1])


VT.div <- data.frame(melt(VT.div),Type=rep("VT"))
NVT.div <- data.frame(melt(NVT.div),Type=rep("NVT"))
div4 <- rbind(VT.div,NVT.div)

ggplot(div4, aes(x=value, fill=Type, color=Type)) +
  geom_density(alpha = 0.2) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(values=c(colos[4],colos[6]))+
  scale_color_manual(values=c(colos[4],colos[6]))+
  labs(x="Shannon diversity", y = "Density") +
  theme_minimal() 


######## TAJIMA'S D ##########

# create a list of matrices instead of char strings

require(stringr)
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
  tp.lg[l] <- pegas::tajima.test(as.DNAbin(lg.mat[[l]]))$Pval.normal
 # tp.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$Pval.normal
}

for (l in 1:length(sm.mat)) {
  #convert character strings to binary DNA
 # td.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$D
  tp.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$Pval.normal
}


library(reshape2)
ddata <- c(psmall=td.sm, plarge=td.lg)
ddata <- data.frame(Weight=c(rep("Small",length(td.sm)),rep("Large",length(td.lg))),D=ddata)
ggplot(ddata, aes(x=D, y=Weight, fill=Weight)) +
  geom_boxplot(alpha = 0.4) +
  #geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="D value", y = "Weight") +
  theme_minimal() 

ggplot(pdata, aes(x=pval, fill=Weight)) +
  geom_histogram(alpha = 0.4, bins=200) +
  #geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(values=c(co[1],co[2]))+
  scale_color_manual(values=c(co[1],co[2]))+
  labs(x="p-value", y = "Count") +
  theme_minimal() 


sigps <- pdata %>% group_by(Weight) %>% dplyr::count(sig=pval<=0.05) %>% group_by(Weight) 
ggplot(filter(sigps, is.na(sig)==F), aes(x=n, y=sig, fill=Weight)) +
  geom_bar(stat="identity",position = position_dodge(),alpha = 0.4) +
  #geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  # geom_density(alpha=0,aes(color=vax)) +
  scale_fill_manual(labels = c("Small","Large"), values=c(co[2],co[1]))+
  scale_color_manual(labels = c("Small","Large"), values=c(co[2],co[1]))+
  labs(x="Number of p-values", y = " ", fill="Weight") +
  theme_minimal() +
  scale_y_discrete(labels=c("TRUE" = "p < 0.05", "FALSE" = "p > 0.05"))


pdata <- c(psmall=tp.sm, plarge=tp.lg)
pdata <- data.frame(Weight=c(rep("Small",length(td.sm)),rep("Large",length(td.lg))),pval=pdata)


## look at all the ones that were significant by tajima
sigD <- data.frame(COG = c(names(pneumo.seq.df.sm),names(pneumo.seq.df.lg)),melt(pdata)) %>% group_by(variable)  %>% filter(value<=0.05) 

# now look at other properties of them --- hmm properties are for each isolate
sigD.isoinfo <- data.frame(mass.info,select(massdata, sigD$COG))


### hmm should i look at diversity of taxa instead of diversity of COGs ??? 
### then i could compare across groups and stuff??? hmmmmmmmm --- i can just group first -- lets do that 
            

#### Tajima's D for VT/NVT #####


x =pneumo.seq.df.lg[["CLS02665"]][["name"]]
y=as.character(mass.info$Taxon) 
sum(mapply(function(x, y) any(x %in% y), 
       strsplit(pneumo.seq.df.lg[["CLS02665"]][["name"]], "\\s+"), strsplit(as.character(mass.info$Taxon) , "\\s+")))







### ignor this bit for now
seq.lg.VT <-  pneumo.seq.df.lg[[1]]$name

seq.lg.VT <- pneumo.seq.df.lg[names(VT)]

seq.lg.VT <- list()
for (i in 1:length(pneumo.seq.df.lg)){
  if (names(pneumo.seq.df.lg[i]) %in% names(VT)){
    seq.lg.VT[[i]] <- pneumo.seq.df.lg[[i]]
  }
}

lg.mat <- list()
for (c in 1:221) {
  mat1 <- matrix(NA, nrow=length(pneumo.seq.df.lg[[c]][["sequence"]]),ncol=str_length(pneumo.seq.df.lg[[c]][["sequence"]][1]))
  for (m in 1:length(pneumo.seq.df.lg[[c]][["sequence"]])) {
    mat1[m,] <- s2c(pneumo.seq.df.lg[[c]][["sequence"]][m])
  }
  lg.mat[[c]] <- mat1
}
sm.mat <- list()
for (c in 1:221) {
  mat1 <- matrix(NA, nrow=length(pneumo.seq.df.sm[[c]][["sequence"]]),ncol=str_length(pneumo.seq.df.sm[[c]][["sequence"]][1]))
  for (m in 1:length(pneumo.seq.df.sm[[c]][["sequence"]])) {
    mat1[m,] <- s2c(pneumo.seq.df.sm[[c]][["sequence"]][m])
  }
  sm.mat[[c]] <- mat1
}

            
            
            
            
            
