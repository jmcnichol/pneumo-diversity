#functions for pneumo diversity calculations
# packages 
library(entropy) #for shannon entorpy 

pneumo.test <- seqinr::read.fasta("gCOG_sequences/CLS00005.fasta", forceDNAtolower = F)




# Calculate shannon entropy