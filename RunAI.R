############################################ RUNNING THE AGE-INDEPENDENT BACK-CALCULATION MODELS ###############################

### NOTE: please modify the directories as suitable

### Four incidence MODELS
## 1) RW 1 order
## 5) ts 10k - thin plate spl linear shrink
## 6) bspl 10k - first order penaltyS
## 4) Gaussian process

rm(list=ls())

### R packages necessary
library(rstan)
library(mgcv)

### STAN setup
### Specifying number of cores
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

### loading a simulated data-set
## (real data cannot be shared, as not allowed by PHE)
### this is available in the Github repo
load("/home/fran/Thesis/RealData/Github/TempDataAI.RData") 

### Data-set (data) contains
ls(data)
## nquar - number of quarters considered
## AIDS - time series of AIDS diagnoses
## HIV - time series of HIV diagnoses
## CD4 - matrix (nquar x 4) of CD4 diagnoses by state
## init.prev - initial number of undiangosed individuals in the 4 CD4 states at time 1
## q - the progression probabilities

### picking one of the incidence models
## model.ind denotes which model is considered
## this must be between 1 and 4
## (see description at the top of the file)
model.ind <- 1
if(!model.ind %in% 1:4)stop("model.ind must be in 1:4")

## RW 1st order
if(model.ind==1){
  stan.data <- list("HIV"=data$HIV, "AIDS"=data$AIDS, "CD4"=data$CD4,
                    "nquar"=data$nquar, "q"=data$q, "init_prev"=data$init.prev)
  model.txt <- "RW1AI.stan"
  pars.save <- c("gamma","var1","vardelta","d","log_lik","exp_HIV_dx","exp_AIDS_dx","exp_p")
}


if(model.ind==2){
  model <- "ts_10k"
  b <- "ts"
  k <- 10
  m <- 2
  model.txt <- "splAI.stan"
  pars.save <- c("beta","gamma","lambda","sigma","vardelta","d","log_lik","sigma","exp_HIV_dx","exp_AIDS_dx","exp_p")
}


if(model.ind==3){
  model <- "bs_10k_ord1"
  b <- "bs"
  k <- 10
  m <- c(2,1)
  model.txt <- "splAI.stan"
  pars.save <- c("beta","gamma","lambda","sigma","vardelta","d","log_lik","sigma","exp_HIV_dx","exp_AIDS_dx","exp_p")
}

## GP
if(model.ind==4){
  model <- "GP"
  sigma.sq <- 0.001 ## fixed nugget for gp
  range01 <- function(x){(x-min(x))/(max(x)-min(x))} ## scale inputs of GP in 0-1
  inputs <- range01(1:data$nquar)
  stan.data <- list("HIV"=data$HIV, "AIDS"=data$AIDS, "CD4"=data$CD4, "nquar"=data$nquar,
                    "q"=data$q, "init_prev"=data$init.prev, "x"=inputs, "sigma_sq"=sigma.sq)
  model.txt <- "GPAI.stan"
  # model.txt <- "GP1AI.stan"
  pars.save <- c("gamma","inv_rho","eta","vardelta","d","log_lik","exp_HIV_dx","exp_AIDS_dx","exp_p")
}


### Creating design matrix (x) for splines
### using jagam function from the mgcv package  

### To do so need to create some
### "fake" data in mgcv format
if(model.ind %in% c(2:3)){
    ### random data
  tmp <- runif(data$nquar)
  
  ### Creating data
  jags.data <- list(
    qts = 1:data$nquar,
    D = tmp
  )
  
  ### This creates the spline object
  ### and a bug file with some stan code
  ### diagonalize = TRUE makes a reparameterization
  ### so that the priors iid normal
  jagam.out <- jagam(D ~ s(qts, bs=b, k=k, m=m),
                       family=gaussian, data=jags.data, 
                       file="/scratch/fran/FinalSimStudy/splAI.bug", diagonalize=TRUE)
  
  ### Grepping the design matrix x
  x <- jagam.out$jags.data$X
  
  ### Removing mgcv specific stuff
  rm(jagam.out, jags.data, tmp)
  
  ### Augmenting data to include x
  stan.data <- list("HIV"=data$HIV, "AIDS"=data$AIDS, "CD4"=data$CD4, "nquar"=data$nquar,
                    "q"=data$q, "init_prev"=data$init.prev, "X"=x, "ninfpars"=ncol(x))
}

### Running the model
fit <- stan(file = paste("/home/fran/Thesis/RealData/Models/",model.txt,sep=""), data = stan.data,
            iter = 10, chains = 3, seed=1204, control = list("adapt_delta"=0.99),
            pars = pars.save, save_warmup=FALSE)