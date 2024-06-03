#' June 3, 2024
#' This script performs Tajim's D with benjamini-hochberg correction on the "small"
#' and "large" weight NFDS COGs in Corander et al. 2017 using the Mass. data only.
#' At the end we save a list of all the COGs that were significant and had D>=2. 
#' We use these COGs in the MATLAB model. 

#locus weights
select.strength <- read.table("intLoci_selection_strengths.csv")
colnames(select.strength) <- c("locus","weight")
unique(select.strength$weight) # there are only two weights 

### Small weight COGs ####
#read in small weight COGs
raw.files <- list.files('SmallCOGS/')
file.list <- pneumo.seq.df.sm <- list()
for (i in 1:length(raw.files)) {
  file <- paste0("SmallCOGS/",raw.files[i])
  file.list[[i]] <- seqinr::read.fasta(file,seqtype= "AA", as.string = T, forceDNAtolower = F)
  
  pneumo.seq.df.sm[[i]] <- data.frame(name=paste(getAnnot(file.list[[i]])), sequence=paste0(file.list[[i]]))
  pneumo.seq.df.sm[[i]]$name <- gsub(">", "", pneumo.seq.df.sm[[i]]$name)
}
names(pneumo.seq.df.sm) <- gsub(".aligned.fasta","", raw.files)

#### Tajima's D #####
sm.mat <- list()
for (c in 1:684) {
  mat1 <- matrix(NA, nrow=length(pneumo.seq.df.sm[[c]][["sequence"]]),ncol=str_length(pneumo.seq.df.sm[[c]][["sequence"]][1]))
  for (m in 1:length(pneumo.seq.df.sm[[c]][["sequence"]])) {
    mat1[m,] <- s2c(pneumo.seq.df.sm[[c]][["sequence"]][m])
  }
  sm.mat[[c]] <- mat1
}

#compute tajimas D and p value
tp.sm = rep(NA,length(pneumo.seq.df.sm))
td.sm = rep(NA,length(pneumo.seq.df.sm))
for (l in 1:length(sm.mat)) {
  #convert character strings to binary DNA
  td.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$D
  tp.sm[l] <- pegas::tajima.test(as.DNAbin(sm.mat[[l]]))$Pval.beta
}

# significant p values #####
#multiple test correction
tp.sm.p <- p.adjust(tp.sm,method = "BH")

hist(td.sm[which(tp.sm.p < 0.05)]) # most of the sig p vlaues correspond to D values below -2. 
hist(td.sm[which(tp.sm.p < 0.05 & td.sm >= 2)])
# we will keep D values above 2
seq.sm.D <- pneumo.seq.df.sm[which(tp.sm.p < 0.05 & td.sm >= 2)]
cog.sm.D <- names(pneumo.seq.df.sm)[which(tp.sm.p < 0.05 & td.sm >= 2)]

### Large weight COGs ####
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
lg.mat <- list()
for (c in 1:221) {
  mat1 <- matrix(NA, nrow=length(pneumo.seq.df.lg[[c]][["sequence"]]),ncol=str_length(pneumo.seq.df.lg[[c]][["sequence"]][1]))
  for (m in 1:length(pneumo.seq.df.lg[[c]][["sequence"]])) {
    mat1[m,] <- s2c(pneumo.seq.df.lg[[c]][["sequence"]][m])
  }
  lg.mat[[c]] <- mat1
}
tp.lg = rep(NA,length(pneumo.seq.df.lg))
td.lg = rep(NA,length(pneumo.seq.df.lg))
for (l in 1:length(lg.mat)) {
  #convert character strings to binary DNA
  #returns indices of segregating sites -- i just want to know how many there are, so add em up
  td.lg[l] <- pegas::tajima.test(as.DNAbin(lg.mat[[l]]))$D
  tp.lg[l] <- pegas::tajima.test(as.DNAbin(lg.mat[[l]]))$Pval.beta
}
# significant p values #####
tp.lg.p <- p.adjust(tp.lg,method = "BH")
hist(td.lg[which(tp.lg.p < 0.05)]) # most of the sig p vlaues correspond to D values below -2. 
hist(td.lg[which(tp.lg.p < 0.05 & td.lg >= 2)])
# we will keep D values above 2
seq.lg.D <- pneumo.seq.df.lg[which(tp.lg.p < 0.05 & td.lg >= 2)]
cog.lg.D <- names(pneumo.seq.df.lg)[which(tp.lg.p < 0.05 & td.lg >= 2)]
### COGs with significant D ####

write.csv(c(cog.lg.D,cog.sm.D),"sigDcognames.csv", col.names = F, row.names = F,quote=F)







