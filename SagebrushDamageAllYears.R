library("rethinking")
setwd("/Users/mmcmunn/Desktop/Yang_Lab/Sagebrush/")
data2012<-read.csv(file = "Sagebrush.Damage.Data.2012.csv", header = T)
data2013<-read.csv(file = "2013 Sagehen/Day.Night.Induction.2013.ALL.DATA.csv", header = T)
data2014<-read.csv(file = "2014 Sagehen/2014.ALL.DAMAGE.DATA.csv", header = T)


#a general strategy for splitting observations from treatments into vectors of proportion damage
#encode the following variables in the same way in each of the 3 datasets: 

      #Year
      #Damge time of day
      #water/no water
      #month of damage

#############2012 Data
#remove "branch removal" treatments and all damage, these don't make any sense.

head(data2012,20)
data2012<- data2012[which(c(data2012$Treatment== 2 | 
                  data2012$Treatment== 4 | 
                  data2012$Treatment== 6 | 
                  data2012$Treatment== 8)==TRUE) , ]

#split the text treatment label into elements and grab the first one, this is time of day
spltrt2012 <- strsplit(as.character(data2012$Treatment.Name),split = "\\.")

for(i in 1:length(spltrt2012)){
data2012$damageTime[i] <-   spltrt2012[[ i ]][1] 
}

#change the treatment names to be universal across years
data2012[data2012$damageTime=="No",]$damageTime <- "control"
data2012[data2012$damageTime=="Night",]$damageTime <- "night"
data2012[data2012$damageTime=="Morning",]$damageTime <- "morning"
data2012[data2012$damageTime=="Afternoon",]$damageTime <- "afternoon"

data2012$year <- 2012
data2012$water <- "noWater"
data2012$monthDamage <- "july"
data2012[data2012$damageTime=="control",]$monthDamage <-"control"

data2012$damLeaves <- data2012$Damaged.Leaves
data2012$totLeaves <- data2012$Total.Leaves

#############2013 Data
data2013$damageTime <- NA

data2013[data2013$Treatment=="Control",]$damageTime <- "control"
data2013[data2013$Treatment=="Night",]$damageTime <- "night"
data2013[data2013$Treatment=="Morning",]$damageTime <- "morning"
data2013[data2013$Treatment=="Day",]$damageTime <- "afternoon"



data2013$year <- 2013
data2013$water <- "noWater"
data2013$monthDamage <- "july"
data2013[data2013$Treatment=="Control",]$monthDamage <-"control"

data2013$damLeaves <- data2013$Damaged.Leaves
data2013$totLeaves <- data2013$Total.Leaves

#############2014 Data

data2014$damageTime <- NA
data2014$water <- NA
data2014$monthDamage <- NA

data2014[data2014$time=="Control",]$damageTime <- "control"
data2014[data2014$time=="Night",]$damageTime <- "night"
data2014[data2014$time=="Afternoon",]$damageTime <- "afternoon"
data2014[data2014$time=="Morning",]$damageTime <- "morning"

data2014[data2014$w.treatment=="no water",]$water <- "noWater"
data2014[data2014$w.treatment=="WATER",]$water <- "water"

data2014[data2014$month=="May",]$monthDamage <- "may"
data2014[data2014$month=="July",]$monthDamage <- "july"
data2014[data2014$month=="Control",]$monthDamage <- "control"

data2014$year <- 2014

data2014$damLeaves <- data2014$dam.leaves
data2014$totLeaves <- data2014$total.leaves

#combine into one data.frame, creating unique columns and collapsing synonyms across years manually
data.temp <- merge(data2012, data2013, all=TRUE)
allData <- merge(data.temp, data2014, all=TRUE)

#get rid of NA's, there are 4
allData <- allData[-which(is.na(allData$damLeaves)==TRUE) , ]

#use aggregate to create a list of unique combinations of the four key factors
treatments <- aggregate(allData$totLeaves, 
                        by = list(allData$year, allData$monthDamage, allData$damageTime, allData$water ), 
                        FUN=length)

treatments <-  treatments[ order(treatments[,1], treatments[,2], treatments[,3], treatments[,4]) , ]

colnames(treatments) <- c( "year", "monthDamage","damageTime", "water","reps" )


head(allData)
forLeafCount <- allData[allData$damageTime!="allDamage",]
aggregate(forLeafCount$totLeaves, by=list(forLeafCount$year),FUN=mean)
aggregate(forLeafCount$totLeaves, by=list(forLeafCount$year),FUN=sum)
length(forLeafCount$totLeaves)
length(forLeafCount$totLeaves)


#create a global treatment number for the purposes of resampling within these data groupings 
#(all the means we want to investigate are within these groups)

treatments$treatmentGlobal <- 1:nrow(treatments)

#match this number back to original data using paste and those 4 columns
treatments$toMatch <- paste(treatments$year, treatments$monthDamage, treatments$damageTime, treatments$water, sep=".")
write.csv(treatments[,-ncol(treatments)], file="printtreatments.csv")

allData$toMatch <- paste(allData$year, allData$monthDamage, allData$damageTime, allData$water, sep=".")

allData$treatmentGlobal <- treatments[match(allData$toMatch, treatments$toMatch),"treatmentGlobal"]

damTreatData <- allData[ ,c("treatmentGlobal","damageTime", "year", "water", "monthDamage", "damLeaves", "totLeaves" )]
damTreatData <- damTreatData[ order(damTreatData$treatmentGlobal) , ]
damTreatData$percentDam <- damTreatData$damLeaves / damTreatData$totLeaves

aggregate(damTreatData$percentDam, by=list(damTreatData$treatmentGlobal),mean)


allD <- damTreatData
head(allD)

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##################################################
#Stan model to estimate year damage probability parameters
nowaterD <- allD[ allD$water=="noWater" , ]
nowaterDJul <- nowaterD[nowaterD$monthDamage!="may" , ]

modelString<-
"
data{
  int<lower=1> nPlants; //

  int damLeaves[nPlants];
  int totLeaves[nPlants];
  int year[nPlants];

  int<lower=1> nYears;
}
parameters{
  real<lower=0> theta;
  vector[nYears] b;
}
model{
  vector[nPlants] pbar;
  theta ~ exponential( 1 );
  b ~ normal( 0 , 10 );

  for ( i in 1:nPlants ) {
    pbar[i] <- b[year[i]];
    pbar[i] <- inv_logit(pbar[i]);

  }
  damLeaves ~ beta_binomial( totLeaves , pbar*theta , (1-pbar)*theta );
}
generated quantities{
  vector[nPlants] pbar;
  real dev;
  dev <- 0;
  for ( i in 1:nPlants ) {
    pbar[i] <- b[year[i]];
    pbar[i] <- inv_logit(pbar[i]);
  }
  dev <- dev + (-2)*beta_binomial_log( damLeaves , totLeaves , pbar*theta , (1-pbar)*theta );
}
"
head(allD)
year <- as.numeric(as.factor(nowaterDJul$year))
totLeaves <- nowaterDJul$totLeaves
damLeaves <- nowaterDJul$damLeaves
nYears <- length(unique(year))
nPlants <- length(damLeaves)


writeLines(modelString, con = "sagebrushYears.stan")
fit <- stan(file = "sagebrushYears.stan", data = list(damLeaves, totLeaves,year, nYears, nPlants), 
            iter = 2000, chains = 4, cores = 4)

print(fit)
quantile(logistic(fit @sim$samples$b[1]), c(0.025, 0.5, 0.975))

lines(fit @sim$samples[[2]]$a, col = "red")


plot(fit, pars = "b")
pairs(fit)
print(fit, digits = 3)


library(rethinking)
data(UCBadmit)
d <- UCBadmit

m11.5 <- map2stan(
  alist(
    admit ~ dbetabinom(applications,pbar,theta),
    logit(pbar) <- a + beta[dept],
    a ~ dnorm(0,2),
    theta ~ dexp(1)
  ),
  data=d,
  constraints=list(theta="lower=0"),
  start=list(theta=3),
  iter=4000 , warmup=1000 , chains=2 , cores=2 )


library(rethinking)
data(reedfrogs)
d <- reedfrogs
# make the tank cluster variable
d$tank <- 1:nrow(d)
# fit
m12.1 <- map2stan(
  alist(
    surv ~ dbinom( density , p ) ,
    logit(p) <- a_tank[tank] ,
    a_tank[tank] ~ dnorm( 0 , 5 )
  ), data=d )

stancode(m12.1)



precis(m11.5)

stancode(m11.5)
