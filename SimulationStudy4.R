##################################################################################
##        Semiparametric Mixed-Effects Models with Censored Responses           ##
## ---------------------------------------------------------------------------- ## 
## Article: A semiparametric mixed-effects model for censored longitudinal data ##
## Authors: Thalita B. Mattos, Larissa A. Matos and Victor H. Lachos            ##
## ---------------------------------------------------------------------------- ##
## SIMULATION STUDY 4                                                           ##
##################################################################################

rm(list=ls(all=TRUE))

# loading functions 

source("smec.R")
source("utils_smec.R")
source("em_nlmec.R")

# loading packages

library(nlme)

## --------------------------- ## 
##  Function to generate data  ##
## --------------------------- ## 

genNLMix <- function(betas,sigmae,D1,ttc,nj,semente,percCens)
{
  p <- length(betas)
  q1 <- dim(D1)[2]
  m <- length(nj)
  n <- nj[1]
    
  y <- matrix(0,m,n)
  
  set.seed(semente)
  for(i in 1:m)
  {
    tt1 <- as.vector(ttc[(sum(nj[(1:i)-1])+1) : (sum(nj[1:i])),  ])
    Ei <-  diag(1,nrow=length(tt1),ncol=length(tt1))
    Omega <- sigmae*Ei
    b <- matrix(rmvnorm(1,c(0,0),D1),q1,1)
    e <- matrix(rmvnorm(1,rep(0,nj[i]),Omega),nj[i],1)
    y[i,] <- exp(betas[1] + b[1,1]) + exp(betas[2])/(1+exp((tt1-exp(betas[3]))/exp(betas[4] + b[2,1]))) + e
  } 
  
  y <- matrix(t(y),m*n,1)
  
  aa <- sort(y,decreasing = FALSE)
  cutoff <- aa[ceiling(percCens*length(y))]
  cc <- matrix(1,length(y),1)*(y <= cutoff)
  y[cc==1] <- cutoff
  LL <- rep(-Inf,length(y))
  LU <- as.vector(y)
  
  x <- ttc
  
  return(list(y=y,cc=cc,x=x,nj=nj,ttc=ttc,cutoff=cutoff,LL=LL,LU=LU))
}

## ---------------- ## 
##  Create folders  ##
## ---------------- ##

CreateDir <- function(mainDir,subDir){
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir),showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))
  }
}

## ------------------- ## 
##  SIMULATION STUDY 4 ##
## ------------------- ##

mainDir <- "~./Simulation/Simulation4"

nSim <- 2

## ----------------------------------------------------

subDir <- paste("Cens",0.10,"m",100,"NL",sep = "_")
CreateDir(mainDir,subDir)

betas <- c(1.6094,0.6931,3.8067,2.3026)
sigmae <- 0.55
D1 <- matrix(c(0.0025,-0.001,-0.001,0.01),2,2)
tt <- seq(0,90,10)
m <- 100
ttc <- matrix(rep(tt,m),ncol=1)
n <- length(tt)
nj <- n*rep(1,m)
percCens <- 0.10

for(count in 1:nSim){

  semente <- 145 + count
  
  dados <- genNLMix(betas,sigmae,D1,ttc,nj,semente,percCens)
  y <-dados$y
  x <-dados$x
  cc <- dados$cc 
  nj <- dados$nj
  ttc <- dados$ttc
  cluster <- kronecker(as.matrix(seq(1:m)),as.matrix(rep(1,n)))
  LL <- dados$LL
  LU <- dados$LU 
  
  Ns <- matrix(0,nrow=length(y),ncol=length(tt)) 
  for(j in 1:length(tt))
  {
    for(i in 1:length(ttc))
    {
      if(ttc[i] == tt[j]){Ns[i,j] <- 1} else {Ns[i,j] <- 0}
    }
  }
  
  ### ------- SAVE DATA
  
  yG <- y
  ccG <- cc
  
  write.table(t(yG),file="yG.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
  write.table(t(ccG),file="ccG.txt",sep=" ",row.names=F,col.names=F, append = TRUE)

  ## -------------- ##
  ## INITIAL VALUES ##
  ## -------------- ##
  
  dataf <- data.frame(y=y,x=x,cluster=cluster)
  betastart <- c(betas)
  
  fit.initial <- try(nlme(y ~ exp(A) + exp(B)/(1 + exp((x-exp(C))/exp(D))), data=dataf,
                       fixed = A + B + C + D ~ 1,
                       random = pdSymm(A + D ~ 1),
                       groups = ~ cluster,
                       start = list(fixed=betastart)))
  
  if(class(fit.initial) != "try-error"){
    
    betasI <- fixed.effects(fit.initial)
    sigma2I <- (fit.initial$sigma)^2
    DI <- round(as.matrix(fit.initial$modelStruct$reStruct$cluster),6)
    
    initial <- list(beta=betasI,sigma2=sigma2I,alphas=DI,lambda=100)
    
    write.table(count,file="countG.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
    
    ## -------------------------------------------------------
    ## Nonlinear - normal : yij = nonlinearfunction + error
    ## -------------------------------------------------------
    
    fit.normal <- try(EM_nl(cc=cc,tt=tt,y=y,x=x,nj=nj,initial=initial))
    
    if(class(fit.normal) != "try-error"){
      
      theta_nl <- c(fit.normal$beta1, fit.normal$sigmae, fit.normal$dd)
      ep_nl <- sqrt(diag(fit.normal$varbeta))
      loglik_nl <- fit.normal$loglik
      AIC_nl <- fit.normal$AIC
      BIC_nl <- fit.normal$BIC
      bi_nl <- apply(fit.normal$ubi,1,sum)
      yi_nl <- apply(fit.normal$uyi,1,sum)

      write.table(t(theta_nl),file="theta_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      write.table(t(ep_nl),file="ep_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      write.table(loglik_nl,file="loglik_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      write.table(AIC_nl,file="AIC_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      write.table(BIC_nl,file="BIC_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      write.table(t(bi_nl),file="bi_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      write.table(t(yi_nl),file="yi_nl.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
      
      write.table(count,file="countNL.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
      
      ## -------------------------------------------------------------
      ## SEMIPARAMETRIC - normal: yij = f(tij) + b0 + b1*tij + error
      ## -------------------------------------------------------------
      
      z <- cbind(1,ttc)
      
      fit.semi <- try(nsmec(y=y,cc=cc,x=NULL,z=z,tt=tt,ttc=ttc,nj=nj,LL=LL,LU=LU,Ns=Ns,initial=initial,struc="unc",lambda.fixed=FALSE,iter.max=200,precision=1e-6))
      
      if(class(fit.semi) != "try-error")
      {
        
        theta_semi <- c(fit.semi$ff, fit.semi$sigmae, fit.semi$dd,fit.semi$lambda)
        ep_semi <- sqrt(diag(solve(fit.semi$Infbetasff)))
        loglik_semi <- c(fit.semi$loglik,fit.semi$loglikp)
        AIC_semi <- fit.semi$AIC
        BIC_semi <- fit.semi$BIC
        bi_semi <- apply(fit.semi$ubi,1,sum)
        yi_semi <- apply(fit.semi$uyi,1,sum)
        
        write.table(t(theta_semi),file="theta_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        write.table(t(ep_semi),file="ep_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        write.table(loglik_semi,file="loglik_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        write.table(AIC_semi,file="AIC_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        write.table(BIC_semi,file="BIC_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        write.table(t(bi_semi),file="bi_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        write.table(t(yi_semi),file="yi_semi.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        
        write.table(count,file="countS.txt",sep=" ",row.names=F,col.names=F,append = TRUE)
        
      } #end try semi
      
    } #end try Normal
    
  } #end try chute
  
}    #end count




