###Hill climbing algorithim to fit S-locus genotypes to plants crossed using a connected small diallel design
###Developed by John Hill Price


library(dplyr)
library(foreach)
library(parallel)
library(doParallel)
####read in data and rule set
dir <- "C:\\home\\"


data <- read.csv(paste(dir,"CrossData.csv", sep = ""), header = T)
rules_all <- read.csv(paste(dir,"all_rules.csv", sep = ""))

#Number of runs per rule
num <- 4000
#Number of cores to use for parallel processing
cor <- 24 

#Penalty for a mismatch between predicted and observed results
cost <- .75
#Threshold seed set value. Observations below this number will be considered "incompatible"
Threshold <- 0.2

#Each hill-climbing chain will end if it goes through this many chnages without improving the score
End <- 1000
##Anchor individual (number is the place of the selected individual in a list of all unique individuals) . 
##Select one individual in the population to be set to the same S-locus genotype across all runs. 
##We recommend selecting an individual used in a relatively high number of crosses
anchor <- 2
###output info (include directory and file prefix)

output <- "C:\\home\\results\\test"

###Here begins the code, should be able to just hit Go
for(h in 1:(ncol(rules_all) - 2)){

prefix <- colnames(rules_all)[h+2]
rules <- rules_all[c(1,2,h+2)]
colnames(rules) <- c("MomGeno", "DadGeno", "Odds")




####recode seed set percentage to binary. 
data2 <- data
for( i in 1:length(data2$Plant)){
  if(data2[i,3] > Threshold){
    data2[i,4] <- 1
  }else( data2[i,4] <- 0)
}

w <- max(rules[,3])

for(i in 1:length(rules[,1])){
 if(rules[i,3] == w){
   rules[i,3] <- cost} else{
     rules[i,3] <- (1-cost)
   }
}

###Set up parent list
scorelist <- vector(length = num)
rulelist <- matrix(nrow = length(parents[,1]), ncol= (num+1))
rulelist[,1] <- parents[,1]

parents <- as.vector(unique(data$Plant))
dad <- as.vector(unique(data$FatherPlant))
parents <- as.data.frame(unique(c(parents, dad)))                 
colnames(parents) <- "Plant"
###assign genotypes at random, then assign an "anchor" individual a particular genotype (arbitrarily chosen as AC)


registerDoParallel(cor)
rule_list <- foreach(i=1:num, .combine = 'cbind') %dopar% {
  
  parents$MomGeno <- sample(c("AC","AD","BC","BD"), length(parents$Plant), replace = T)
  
  
  parents[anchor,2] <- "AC"
  set <- c(1:length(parents$Plant))[!c(1:length(parents$Plant)) %in% anchor ]
  
  ####set up the other files, giving each cross in the dataset genotypes for the two parents
  dads <- parents
  colnames(dads) <- c("FatherPlant", "DadGeno")
  
  rules$MomGeno <- as.character(rules$MomGeno)
  rules$DadGeno <- as.character(rules$DadGeno)
  
  
  data3 <- dplyr::left_join(data2, parents, by = "Plant")
  data3 <- dplyr::left_join(data3, dads, by="FatherPlant")
  
  
  data3 <- dplyr::left_join(data3, rules, by=c("MomGeno", "DadGeno"))
  data3$scores <- abs(data3[,4]-data3[,7])
  score <- sum(data3$scores)
  
  
  p <- 0
  
  ##Be sure to remove score2 before a new run
  rm(score2)
  
  
  #############
  ##Randomly change one individual to a different genotype. If that reduces the score, change it and proceed. If not, trya  different change.
  ##If the cycle repeats "End" number of times without the score decreasing, end
  repeat{
    parents2 <- parents
    q <- sample(set, 1)
    parents2[q,2] <- sample(c("AC","AD","BC","BD")[!c("AC","AD","BC","BD") %in% parents2[q,2]], 1)
    dads2 <- parents2
    colnames(dads2) <- c("FatherPlant", "DadGeno")
    data4 <- dplyr::left_join(data2, parents2, by = "Plant")
    data4 <- dplyr::left_join(data4, dads2, by="FatherPlant")
    data4 <- dplyr::left_join(data4, rules, by=c("MomGeno", "DadGeno"))
    data4$scores <- abs(data4[,4]-data4[,7])
    score2 <- sum(data4$scores)
    
    if(score2 < score){
      score <- score2
      data3 <- data4
      parents <- parents2
      p <-0}else(p <- p+1)
    
    if(p == End){
      break
    }}
  
  
  #scorelist[i] <- score
  rulelist[,i+1] <- parents[,2]
}


rule_list <-as.data.frame(cbind(parents[,1], rule_list))

#####score

score_list <- foreach(i=2:num, .combine = 'cbind') %dopar% {
parentscore <- as.data.frame(rule_list[,c(1,i)])
colnames(parentscore) <- c("Plant", "MomGeno")
parentscore$Plant <- as.integer(as.character(parentscore$Plant))
parentscore$MomGeno <- as.character(parentscore$MomGeno)

dadscore <- parentscore
colnames(dadscore) <- c("FatherPlant", "DadGeno")



datascore <- dplyr::left_join(data2, parentscore, by = "Plant")
datascore <- dplyr::left_join(datascore, dadscore, by="FatherPlant")


datascore <- dplyr::left_join(datascore, rules, by=c("MomGeno", "DadGeno"))
datascore$scores <- abs(datascore[,4]-datascore[,7])
sum(datascore$scores)}


######Export
write.csv(score_list, paste(output,"_scores.csv", sep = ""))
write.csv(rule_list, paste(output,"_genotypes.csv", sep = ""))
}
