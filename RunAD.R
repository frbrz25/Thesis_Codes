############################################ RUNNING THE AGE-DEPENDENT BACK-CALCULATION MODELS ###############################

### NOTE: please modify the directories
### (for loading the files) as suitable

### Models that can be run
### the number in bracket denotes the respective model index
## 1) Tps - yr, no age diag
## 2) Tps - yr, age diag 1
## 3) Tps - yr, age diag 2
## 4) Ptens - yr, no age diag
## 5) Ptens - yr, age diag 1
## 6) Ptens - yr, age diag 2
## 7) Ptens - qt, no age diag
## 8) Ptens - qt, age diag 1
## 9) Ptens - qt, age diag 2

rm(list=ls())

### R packages necessary
library(rstan)
library(mgcv)

### STAN setup
### Specifying number of cores
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

### picking one of the incidence models
## model.ind denotes which model is considered
## this must be between 1 and 9
## (see description at the top of the file)
model.ind <- 1
if(! model.ind %in% c(1:9) )stop("model.ind incorrectly specified")

if(model.ind==1) model.txt <- "tps.stan"
if(model.ind==2) model.txt <- "tps_agediag.stan"
if(model.ind==3) model.txt <- "tps_agediag1.stan"
if(model.ind==4) model.txt <- "ptens.stan"
if(model.ind==5) model.txt <- "ptens_agediag.stan"
if(model.ind==6) model.txt <- "ptens_agediag1.stan"
if(model.ind==7) model.txt <- "ptens_quar.stan"
if(model.ind==8) model.txt <- "ptens_quar_agediag.stan"
if(model.ind==9) model.txt <- "ptens_quar_agediag1.stan"

### loading a simulated data-set
## (real data cannot be shared, as not allowed by PHE)
### NOTE: different datasets for quarterly and yearly models
### data are available in the Github repo
if(model.ind %in% c(1:6)){
  load("/home/fran/Thesis/RealData/Github/TempDataAD.RData") 
}else{
  load("/home/fran/Thesis/RealData/Github/TempDataADQt.RData") 
}

### Data-set (data) contains
ls(data)
## nyr - number of years considered
## nquar - number of quarters considered (only for qt models)
## nage - number of ages considered
## AIDS - matrix (nyr x nage) time-and-age matrix of AIDS diagnoses
## HIV - matrix (nyr x nage) time-and-age matrix of HIV diagnoses
## CD4 - array (nyr x nage x 4) time-and-age CD4 diagnoses, further stratified by CD4 count group
## init.prev - matrix (nage x 4), initial number of undiangosed individuals in the 4 CD4 states, stratified by age, at time 1
## q - the age-dependent progression probabilities for each CD4 state
## NOTE: nyr replaced by nquar for quarterly models

### Creating design matrix X and penalty matrix S1
### using jagam function from the mgcv package  

### To do so need to create some
### "fake" data in mgcv format
if(model.ind %in% c(1:6)){
  ### Yearly data
  
  ### random data
  tmp <- runif(data$nyr*data$nage)
  
  yrs <- 1:data$nyr
  ages <- 1:data$nage
  
  ### Creating data
  jags.data <- list(
    age = rep(ages, each=data$nyr),
    yrs = rep(yrs, times=data$nage),
    D = tmp
  )
  
}else{
  ### Quarterly data
  ### Creating design matrix X and penalty matrix S1
  
  ### random data
  tmp <- runif(data$nquar*data$nage)
  
  yrs <- 1:data$nquar ## yrs is actually quarters, but convenient to call it like his
  ages <- 1:data$nage
  
  ### Creating data
  jags.data <- list(
    age = rep(ages, each=data$nquar),
    yrs = rep(yrs, times=data$nage),
    D = tmp
  )
}

### This creates the spline object
### and a bug file with some stan code
### diagonalize = TRUE makes a reparameterization
### so that the priors iid normal
### This is ignored for tensor product splines
### (this reparameterization can not be done)

### If model.ind in 1:3 deal with thin plate spline
### otherwhise tensor product with marginal cubic B-spline
##' with first order penalty
if(model.ind %in% c(1:3)){
  jagam.out <- jagam(D ~ s(yrs, age, bs="ts", k=80),
                     family=gaussian, data=jags.data, 
                     file="/home/fran/HIV backcalculation/Age dependency/MSMdata14/Bayesian Stuff/STAN/ts.bug", diagonalize=TRUE)
}else{
  m.list <- list(c(2,1),c(2,1)) 
  jagam.out <- jagam(D ~ te(yrs, age, bs=c("ps","ps"), k=c(10,8), m=m.list),
                     family=gaussian, data=jags.data, 
                     file="/home/fran/HIV backcalculation/Age dependency/MSMdata14/Bayesian Stuff/STAN/ptens.bug", diagonalize=TRUE)
}



x <- jagam.out$jags.data$X
s1 <- jagam.out$jags.data$S1
rm(jagam.out, jags.data, tmp, yrs, ages)


### Only use 1000 niter for quarterly model (long running times)
### 2000 for yearly model
if(model.ind %in% c(7:9)){
  n.iter <- 1000
}else{
  n.iter <- 2000
}

n.iter <- 10 ### testing purposes

stan.data <- list("HIV"=data$HIV, "AIDS"=data$AIDS,"CD4"=data$CD4,
                  "nyr"=data$nyr, "nage"=data$nage,"q"=data$q, "X"=x,
                  "init_prev"=data$init.prev, "ninfpars"=ncol(x), "S1"=s1)

p.save <- c("beta","lambda","vardelta","d","under_rep")

if(model.ind %in% c(1:3)) stan.data$S1 <- NULL ## No need for S1 for thin plate splines (iid prior)
if(model.ind %in% c(7:9)){ ## adding quarters to data for quarterly model
  stan.data$nquar <- stan.data$nyr*4
  stan.data$nyr <- NULL 
} 
if(model.ind %in% c(2,3,5,6,8,9)) p.save <- c(p.save,"alpha") ## saving age-specific intercept dx parameters

fit <- stan(file = paste0("/home/fran/Thesis/RealData/Models/", model.txt),
            data = stan.data, iter = n.iter, chains = 3, 
            pars= p.save)
