library(dplyr)
library(seqinr)
library(parallel)
library(vegan)
library(ggplot2)
library(ggpubr)


### Example/test case ####
# have a look at an example file
pneumo.test <- seqinr::read.fasta("gCOG_sequences/CLS00005.fasta", seqtype= "AA", as.string = T, forceDNAtolower = F)
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
#p = 0.4981 

mean(large.weight.div$diversity)
mean(small.weight.div$diversity)

var(large.weight.div$diversity)
var(small.weight.div$diversity)


