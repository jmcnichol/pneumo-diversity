
library(dplyr)
# read in data 
massdata <- read.csv("massdata.csv")
# create df for all the stuff that arent loci id cols 
mass.info <- massdata[,2:14]
select.strength <- read.table("intLoci_selection_strengths.csv")
colnames(select.strength) <- c("locus","weight")
hist(select.strength$weight) #two distinct groups  
unique(select.strength$weight) # there are only two weights 

### data manipulation ####
# choose the cols from massdata that match the weights=
weight.large <- select.strength %>% filter(weight > 0.06)
weight.small <- select.strength %>% filter(weight < 0.06)

nrow(weight.large) #270
nrow(weight.small) #820

# extract the corresponding cols from massdata 
# big weights
data.large.weight <- massdata %>% select(weight.large$locus)
data.large.weight <- cbind(mass.info,data.large.weight)
#small weights
data.small.weight <- massdata %>% select(weight.small$locus)
data.small.weight <- cbind(mass.info,data.small.weight)

# now we need to read in the fasta files for each set. 
