
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
#Bayesian beta-binomial model - just estimate p across all samples
nowaterD <- allD[ allD$water=="noWater" , ]
nowaterDJul <- nowaterD[nowaterD$monthDamage!="may" , ]

modelString <- "
data {
  int<lower=0> I; // number of plants 
  int<lower=0> n[I]; // number of leaves measured
  int<lower=0> y[I]; // number of leaves damaged
}
parameters {
  real<lower=0> p;
  real<lower=0> a;
  real<lower=0> b;
}
transformed parameters {
}
model {
  y ~ binomial(n , p);
  p ~ beta(a, b);
  a ~ gamma(1, 0.1);
  b ~ gamma(1, 0.1);
}
"
writeLines(modelString, con = "allTogetherSage.stan")

#specify data
n <- nowaterDJul$totLeaves
y <- nowaterDJul$damLeaves
sage_dat <- list(I = dim(nowaterDJul)[1], y = y, n = n)

fit <- stan(file = "allTogetherSage.stan", data = sage_dat, 
            iter = 2000, chains = 4)

plot(fit, pars = "p")
pairs(fit)
print(fit, digits = 3)
###################################################
#Bayesian beta-binomial model - estimate p across all years

nowaterD <- allD[ allD$water=="noWater" , ]
nowaterDJul <- nowaterD[nowaterD$monthDamage!="may" , ]

modelString <- "
data {
int<lower=0> I; // number of plants 
int<lower=0> Y; // number of years 
int<lower=0> n[I]; // number of leaves measured
int<lower=0> dam[I]; // number of leaves damaged
int<lower=0> year[I]; // years of experiment

}
parameters {
real<lower=0> p[Y];
real<lower=0> a[Y];
real<lower=0> b[Y];
}
transformed parameters {
}
model {
dam ~ binomial(n , p);
p ~ beta(a, b);
a ~ gamma(1, 0.1);
b ~ gamma(1, 0.1);
}
"
writeLines(modelString, con = "byYearSage.stan")

#specify data
n <- nowaterDJul$totLeaves
dam <- nowaterDJul$damLeaves
year <- as.numeric(as.factor(nowaterDJul$year))

sage_dat <- list(I = dim(nowaterDJul)[1], dam = dam, n = n,
                 year = year, Y = length(unique(year)))

fit <- stan(file = "byYearSage.stan", data = sage_dat, 
            iter = 2000, chains = 4)

plot(fit, pars = "p")
pairs(fit)
print(fit, digits = 3)
















nowaterD <- allD[ allD$water=="noWater" , ]
nowaterDJul <- nowaterD[nowaterD$monthDamage!="may" , ]


n <- nowaterDJul$totLeaves
y <- nowaterDJul$damLeaves
year <- nowaterDJul$year
year <- year-2011

dataList <- list(y = y, n = n, year = year)
modelString = "
model{
        for (i in 1:length(y)){

        y[i] ~ dbin(p[year[i]], n[i])
        }

                for(year in 1:3){

  p[year] ~ dbeta(a[year],b[year])
  a[year] ~ dgamma(1, 0.1)
  b[year] ~ dgamma(1, 0.1)
}}
"
writeLines( modelString , con="betabinomial.jags.txt" )

jagsModel = jags.model( file="betabinomial.jags.txt" , data=dataList , 
                        n.chains=3 , n.adapt=500 )

update( jagsModel , n.iter=1000 )
codaSamples = coda.samples( jagsModel , variable.names=c("p","a","b") ,
                            n.iter=10000 )
save( codaSamples , file=paste("Year.Mcmc.Rdata") )

parameter.names <- gsub("\\[[0-9]*\\]", "", colnames(codaSamples[[1]]))

#polt for convergence
plot(codaSamples[,"p[1]"])
plot(codaSamples[,"p[2]"])
plot(codaSamples[,"p[3]"])

#plot summary statistics
plotPost(codaSamples[,"p[1]" ], cenTend="median")
plotPost(codaSamples[,"p[2]"] , cenTend="median")
plotPost(codaSamples[,"p[3]"] , cenTend="median")

####################################################
##Bayesian beta-binomial model - estimate p across all years, times, season, waterings
#a non-hierarchical option, where each of the 17 p's are estimated completly seperately. 
#(water and no water 2014 + 12 timeXyear reps + 3 more seasonal May inductions)

n <- allD$totLeaves
y <- allD$damLeaves
treatment <- allD$treatmentGlobal

dataList <- list(y = y, n = n, treatment = treatment)

modelString = "
model{
for (i in 1:length(y)){

y[i] ~ dbin(p[treatment[i]], n[i])
}

for(treatment in 1:17){

p[treatment] ~ dbeta(a[treatment],b[treatment])
a[treatment] ~ dgamma(1, 0.1)
b[treatment] ~ dgamma(1, 0.1)

}}
"
writeLines( modelString , con="betabinomial.jags.txt" )

jagsModel = jags.model( file="betabinomial.jags.txt" , data=dataList , 
                        n.chains=3 , n.adapt=500 )

update( jagsModel , n.iter=10000 )
codaSamples = coda.samples( jagsModel , variable.names=c("p","a","b") ,
                            n.iter=10000 )

save( codaSamples , file=paste("Trt.Estimates.Mcmc.Rdata") )


####################
#plot of p[i] lower and upper 95% HDI limits and medians averaged over 3 chains for all 17 treatments

lowCI <- vector()
highCI <- vector()
medianPost <- vector()
for(h in 1:(17*3)){
  
colname <- colnames(codaSamples[[1]])[h]

CI <- list()
median <- vector()
for(i in 1:3){
  
CI[[i]] <- HDIofMCMC(codaSamples[[i]][,colname])

median[i] <- as.numeric(median(codaSamples[[i]][,colname]))

}

lowCI[h] <- mean(c(CI[[1]][2], CI[[2]][2], CI[[3]][2]))
highCI[h] <- mean(c(CI[[1]][1], CI[[2]][1], CI[[3]][1]))
medianPost[h] <- mean(c(median[1], median[2], median[3]))

}

CIavgs <- cbind(lowCI, highCI, medianPost)
rownames(CIavgs) <-  colnames(codaSamples[[1]])

CIavgs <- as.data.frame(CIavgs)
par(oma = c(10,0,0,0))
plot(CIavgs$median[35:51], pch=16,cex=.5,xaxt="n",xlab="", ylab = "mean percent herbivory", ylim = c(0, .18))
arrows(x0=1:17,,x1=1:17, y0=CIavgs$lowCI[35:51],y1=CIavgs$highCI[35:51],angle=90,length=.1,code=3)
axis(1, at = 1:17, labels = treatments$toMatch, las = 2)
?axis
#go get the descriptors for each of the 17 treatments and make a treatment list
variable.list <- allD[1,]
for (i in 1:17){
  temp <- allD[allD$treatmentGlobal==i,]
  variable.list[i,] <- temp[i,]
  variable.list <- variable.list[ , 1:5]
}
CIavgs <- cbind(variable.list,CIavgs)

####################################################
##Bayesian beta-binomial model - use only july induction treatments in each of 3 years
#just time of day
#use only july data

n <- nowaterDJul$totLeaves
y <- nowaterDJul$damLeaves

unique(nowaterDJul$damageTime)
timeOfDay <- as.numeric(as.factor(nowaterDJul$damageTime))


dataList <- list(y = y, n = n, timeOfDay = timeOfDay)
modelString = "
model{
for (i in 1:length(y)){

y[i] ~ dbin(p[timeOfDay[i]], n[i])
}

for(timeOfDay in 1:4){

p[timeOfDay] ~ dbeta(a[timeOfDay],b[timeOfDay])
a[timeOfDay] ~ dgamma(1, 0.1)
b[timeOfDay] ~ dgamma(1, 0.1)
}}
"
writeLines( modelString , con="betabinomial.jags.txt" )

jagsModel = jags.model( file="betabinomial.jags.txt" , data=dataList , 
                        n.chains=3 , n.adapt=1000 )

update( jagsModel , n.iter=10000 )
codaSamples = coda.samples( jagsModel , variable.names=c("p","a","b") ,
                            n.iter=10000 )
save( codaSamples , file=paste("time.Mcmc.Rdata") )



parameter.names <- gsub("\\[[0-9]*\\]", "", colnames(codaSamples[[1]]))

#polt for convergence
plot(codaSamples[,"p[1]"])
plot(codaSamples[,"p[2]"])
plot(codaSamples[,"p[3]"])
plot(codaSamples[,"p[4]"])

#plot summary statistics
plotPost(codaSamples[,"p[1]" ], cenTend="median")
plotPost(codaSamples[,"p[2]"] , cenTend="median")
plotPost(codaSamples[,"p[3]"] , cenTend="median")
plotPost(codaSamples[,"p[4]"] , cenTend="median")


####################################################
##Bayesian beta-binomial model - use only july induction treatments in each of 3 years
## using year + time of day + year*time of day

nowaterDJul <- nowaterD[nowaterD$monthDamage!="may" , ]

n <- nowaterDJul$totLeaves
y <- nowaterDJul$damLeaves

unique(nowaterDJul$damageTime)
timeOfDay <- as.numeric(as.factor(nowaterDJul$damageTime))
year <- nowaterDJul$year
year <- year-2011

dataList <- list(y = y, n = n, timeOfDay = timeOfDay, year = year)

modelString = "
model{
for (i in 1:length(y)){

y[i] ~ dbin(p[i], n[i])
p[i] <- a1[timeOfDay[i]] + a2[year[i]] + a1a2[timeOfDay[i], year[i]]

}

for(timeOfDay in 1:4){
a1[timeOfDay] ~ dbeta(ag1[timeOfDay],bg1[timeOfDay])
ag1[timeOfDay] ~ dgamma(1, 0.1)
bg1[timeOfDay] ~ dgamma(1, 0.1)
}

for(year in 1:3){
a2[year] ~ dbeta(ag2[year],bg2[year])
ag2[year] ~ dgamma(1, 0.1)
bg2[year] ~ dgamma(1, 0.1)
}

for(timeOfDay in 1:4){
for(year in 1:3){
a1a2[timeOfDay,year] ~ dbeta(ag12[timeOfDay, year],bg12[timeOfDay, year])
ag12[timeOfDay,year] ~ dgamma(1, 0.1)
bg12[timeOfDay,year] ~ dgamma(1, 0.1)
}}
}



"

inits <- list("a1" = c(.01, .01, .01, .01) , "a2" = c(.02, .02, .02)  , 
              "ag1" = c(1,1,1,1)  , "bg1" = c(2,2,2,2) , "ag2" = c(1,1,1), "bg2" = c(2,2,2),
              "a1a2" = matrix(nrow = 4, ncol = 3, data = c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)),
              "ag12" = matrix(nrow = 4, ncol = 3, data = c(1,1,1,1,1,1,1,1,1,1,1,1)), 
              "bg12" = matrix(nrow = 4, ncol = 3, data = c(2,2,2,2,2,2,2,2,2,2,2,2)))
writeLines( modelString , con="betabinomial.jags.txt" )

jagsModel = jags.model( file="betabinomial.jags.txt" , data=dataList , inits = inits,
                        n.chains=3 , n.adapt=100000 )


update( jagsModel , n.iter=20000 )
codaSamples = coda.samples( jagsModel , variable.names=c("a1","a2","ag1","bg1", "ag2", "bg2","a1a2","a0") ,
                            n.iter=20000 )
save( codaSamples , file=paste("time.year.Mcmc.Rdata") )



parameter.names <- gsub("\\[[0-9]*\\]", "", colnames(codaSamples[[1]]))

#polt for convergence
plot(codaSamples[,"a1[2]"])
plot(codaSamples[,"a2[3]"])

head(codaSamples)

#plot summary statistics
plotPost(codaSamples[,"a2[1]" ], cenTend="median")
plotPost(codaSamples[,"a2[2]" ], cenTend="median")
plotPost(codaSamples[,"a2[3]" ], cenTend="median")
plotPost(codaSamples[,"a1[4]" ], cenTend="median")





plotPost(codaSamples[,"a2[1]" ], cenTend="median")



####################################################
##Bayesian beta-binomial model - use only july induction treatments in each of 3 years
## using HIERARCHICAL STRUCTURE, TIME OF DAY WITHIN TIME OF YEAR

nowaterDJul <- nowaterD[nowaterD$monthDamage!="may" , ]

n <- nowaterDJul$totLeaves
y <- nowaterDJul$damLeaves

unique(nowaterDJul$damageTime)
timeOfDay <- as.numeric(as.factor(nowaterDJul$damageTime))
year <- nowaterDJul$year
year <- year-2011

dataList <- list(y = y, n = n, timeOfDay = timeOfDay, year = year)

modelString = "
model{
for (i in 1:length(y)){

y[i] ~ dbin(p[i], n[i])
p[i] ~ dbeta(omega[year[timeOfday[i]]]*(kappa[year[timeOfday[i]]]-2)+1 , 
                        (1-omega[year[timeOfday[i]]])*(kappa[year[itimeOfday[i]]]-2)+1)

}

for(year in 1:3){

omega[year[timeOfDay]] ~ dbeta(omegaY[year]*(kappaY[year]-2)+1 , 
                        (1-omegaY[year])*(kappaY[year]-2)+1)

kappa[year[timeOfDay]] <- kappaMinusTwo[year[timeOfDay]] + 2
kappaMinusTwo[year[timeOfDay]] ~ dgamma( 0.01 , 0.01 )

omegaY[year] ~ dbeta( 1, 1)
kappaY[year] <- kappaMinusTwoY[year] + 2
kappaMinusTwoY[year] ~ dgamma( 0.01 , 0.01 )
}}
"

writeLines( modelString , con="betabinomial.jags.txt" )

jagsModel = jags.model( file="betabinomial.jags.txt" , data=dataList ,
                        n.chains=3 , n.adapt=1000 )


update( jagsModel , n.iter=20000 )
codaSamples = coda.samples( jagsModel , variable.names=c("a1","a2","omega",
                            "omegaT", "omegaY", "kappa","kappaT","kappaY") ,
                            n.iter=20000 )
save( codaSamples , file=paste("time.year.Mcmc.Rdata") )



parameter.names <- gsub("\\[[0-9]*\\]", "", colnames(codaSamples[[1]]))

#polt for convergence
plot(codaSamples[,"a1[2]"])
plot(codaSamples[,"a2[3]"])

head(codaSamples)

#plot summary statistics
plotPost(codaSamples[,"a2[1]" ], cenTend="median")
plotPost(codaSamples[,"a2[2]" ], cenTend="median")
plotPost(codaSamples[,"a2[3]" ], cenTend="median")
plotPost(codaSamples[,"a1[4]" ], cenTend="median")





plotPost(codaSamples[,"a2[1]" ], cenTend="median")




