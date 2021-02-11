###Single marker regression approach to mapping the S-locus in plants crossed using a connected small diallel design
###Developed by John Hill Price

library(tidyverse)
library(car)

genodata <- read.csv("three_alleles_recode.csv")
crossdata <- read.csv("CrossData.csv")

#Recode seed set values to compatible or incompatible, based on a threshold value
Threshold <- .2
for( i in 1:length(crossdata$Plant)){
  if(crossdata[i,3] > Threshold){
    crossdata[i,4] <- 1
  }else(crossdata[i,4] <- 0)
}
#Remove any observations without marker data
filtercrossdata <- crossdata %>% filter(Plant %in% genodata$Individual) %>% filter(FatherPlant %in% genodata$Individual) 



#Create tables of maternal and paternal genotypes
Mgeno <- left_join(filtercrossdata, genodata, b = c("Plant" = "Individual"))
Mgeno <- as.data.frame(Mgeno[,-c(1:6)])


Dgeno <- left_join(filtercrossdata, genodata, b = c("FatherPlant" = "Individual"))
Dgeno <- as.data.frame(Dgeno[,-c(1:6)])

#Impute the mean for missing marker observations
Mgenoi <- Mgeno %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) 
Dgenoi <- Dgeno %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) 

effects1 <- vector(length = ncol(Mgenoi))

success <- filtercrossdata[,4]

#####Regression
for( i in 1:length(colnames(Mgenoi))){
  model12 <- glm(success ~ Mgenoi[,i] + Dgenoi[,i] + Mgenoi[,i]:Dgenoi[,i], family = binomial(link = "logit"))
  if(nrow(summary(model12)$coefficients) == 4){
    effects1[i] <- summary(model12)$coefficients[4,4]
  }else{effects1[i] <- NA}
}

View(effects1)