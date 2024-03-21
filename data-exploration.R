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
pneumo.freq <- pneumo.test |> count(sequence)
pneumo.freq$n <- pneumo.freq$n/nrow(pneumo.test)
#calculate Shannon entropy 
div <- vegan::diversity(pneumo.freq$n)

#calculates Shannon entropy 
shan.div <- function(list.of.dfs){
  div <- rep(0,length(list.of.dfs))
  for (i in 1:length(list.of.dfs)){
    # unique sequences
    pneumo.freq <- list.of.dfs[[i]] |> count(sequence)
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
raw.files <- list.files('gCOG_sequences/')
filenames <- paste0(select.strength$locus ,".out")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df <- list()
for (i in 1:length(files.all)) {
  file <- paste0("gCOG_sequences/",files.all[i])
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
set.seed(12345)
weight.sample <- sample.int(221)
weight.small <- weight.small[weight.sample,]

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
raw.files <- list.files('gCOG_sequences/SmallWeights/')
filenames <- paste0(weight.small$locus,".aligned.fasta")
files.all <- na.omit(raw.files[match(filenames,raw.files)])
file.list <- pneumo.seq.df.sm <- list()
for (i in 1:length(files.all)) {
  file <- paste0("gCOG_sequences/SmallWeights/",files.all[i])
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

colos <- c("#754668","#935116","#5DADE2","#F5B041")

ggplot(shan.comb, aes(x=diversity, fill=weight)) +
  geom_histogram(aes(y=..density..), alpha=0.4, position = "identity", binwidth = 0.1) +
  labs(x="Shannon diversity", y = "Density") +
  scale_fill_manual(values=c(colos[1],colos[3]))+
  theme_minimal() 

ggplot(shan.comb, aes(x=diversity, fill=weight)) +
  geom_density(alpha = 0.4) +
  labs(x="Shannon diversity", y = "Density") +
  scale_fill_manual(values=c(colos[4],colos[3]))+
  theme_minimal() 

ks.test(pneumo.shan.lg$diversity,pneumo.shan.sm$diversity)
#p-value = 0.1182

median(pneumo.shan.lg$diversity)

# number of segregating sites -- number of col with at least two different letters 
num.seg.lg = num.seg.sm = rep(NA,length(pneumo.seq.df.lg))
for (l in 1:length(pneumo.seq.df.lg)) {
  #convert character strings to binary DNA
  #returns indices of segregating sites -- i just want to know how many there are, so add em up
  num.seg.lg[l] <- sum(ape::seg.sites(char2dna(pneumo.seq.df.lg[[l]]$sequence)))
  num.seg.sm[l] <- sum(ape::seg.sites(char2dna(pneumo.seq.df.sm[[l]]$sequence)))
}




