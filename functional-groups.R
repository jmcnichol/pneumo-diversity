# pneumo functional groups 
# looking at the cogs by functional annotation for the 616 mass isolates 

cog.fun <- read.csv("cog-functional.csv")

#names for the cogs i used 
cog.lg <- cog.fun %>% filter(COG %in% weight.large)
cog.sm <- cog.fun %>% filter(COG %in% weight.small)


cog.lg <- cog.lg[,-c(1,2,3,4)]
clg <- matrix(0,ncol=ncol(cog.lg),nrow = nrow(cog.lg))
for (i in 1:nrow(cog.lg)) {
  for (j in 1:ncol(cog.lg)) {
    if(cog.lg[i,j] =="N"){
      clg[i,j] = 0
    }else if (cog.lg[i,j]=="Y"){
      clg[i,j]= 1
    }
  }
}
cog.sm <- cog.sm[,-c(1,2,3,4)]
csm <- matrix(0,ncol=ncol(cog.sm),nrow = nrow(cog.sm))
for (i in 1:nrow(cog.sm)) {
  for (j in 1:ncol(cog.sm)) {
    if(cog.sm[i,j] =="N"){
      csm[i,j] = 0
    }else if (cog.sm[i,j]=="Y"){
      csm[i,j]= 1
    }
  }
}


#%>% group_by(COG) %>% tidyr::gather(value)
  ggplot(data.frame(freq=c(apply(clg, 2, sum), apply(csm, 2, sum)),
                    Weight=c(rep("Large",25),rep("Small",25)), 
                    Function = rep(c(1:25),2)), aes(Function, fill=Weight)) + 
  geom_bar(aes(freq))+
    #scale_alpha_identity(guide = "none") +
  theme( axis.text.x = element_text(angle = 45, hjust = 1))
