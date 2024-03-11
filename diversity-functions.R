#functions for pneumo diversity calculations
# packages 
#library(entropy) #for shannon entorpy 
library(seqinr)
library(vegan)

positionInfo <- function(data, verbose=TRUE){
 # datadf <-  seqinr::read.fasta(file, seqtype= "AA", as.string = T)
#  datadf <- data.frame(name=paste(getAnnot(datadf)), sequence=paste0(datadf))
#  datadf$name <- gsub(">", "", datadf$name)
#  data = datadf
  info <- list()
  for (i in 1:(nchar(as.character(data$sequence[1])))){
    ss <- substr(data[,2],i,i)
    ss <- toupper(ss)
    df <- data.frame(base=unique(ss))
    tmp <- c()
    for ( j in 1:length(unique(ss))){
      tmp[j] <- length(grep(df[j,1], ss))
    }
    df$rep <- tmp
    rownames(df) <- as.character(df[,1])
    df[,1] <- NULL
    
    dff <- subset(df, rownames(df)!="-")
    
    info$shan[i] <- vegan::diversity(t(df))                                    # Shannon entropy of position
    info$shanc[i] <- ifelse(nrow(dff)==0, NA, vegan::diversity(t(dff)))        # Shannon entropy of position removing gaps ("-")
    info$rich[i] <- length(unique(ss))                                  # Position richness
    info$richc[i] <- ifelse(nrow(dff)==0, NA, length(rownames(dff)))    # Position richness removing gaps ("-")
    info$uniq[i] <- list(unique(ss))                                    # Unique bases in position
    info$repe[i] <- list(tmp)                                           # Repetitions of the unique bases in position
    
  }
  
  df <- data.frame(posi=c(1:length(info$shan)),
                   shan=info$shan,
                   shanc=info$shanc,
                   rich=info$rich,
                   richc=info$richc,
                   uniq=unlist(lapply(info$uniq, function(x) paste(x, collapse="|"))),
                   repe=unlist(lapply(info$repe, function(x) paste(x, collapse="|"))))
  return(df)
  
}




# Calculate shannon entropy