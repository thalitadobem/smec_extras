##################################################################################
##        Semiparametric Mixed-Effects Models with Censored Responses           ##
## ---------------------------------------------------------------------------- ## 
## Article: A semiparametric mixed-effects model for censored longitudinal data ##
## Authors: Thalita B. Mattos, Larissa A. Matos and Victor H. Lachos            ##
## ---------------------------------------------------------------------------- ##
## SIMULATION STUDY 3                                                           ##
##################################################################################

rm(list=ls(all=TRUE))

# loading SMEC functions 

source("smec.R")
source("utils_smec.R")

# loading packages

library(nlme)

## --------------------------- ## 
##  Function to generate data  ##
## --------------------------- ## 

genSemiMix <- function(m,x1,x2,tt,beta,sigma,D1,funcSemi,semente,percCens,struc,phi1,phi2)
{
  set.seed(semente)
  
  p <- length(beta)
  q <- dim(D1)[2]
  times <- funcSemi(tt)
  times <- matrix(times,ncol=1)
  
  ttc <- matrix(rep(tt,m),ncol=1)
  
  n <- length(tt)
  nj <- n*rep(1,m)
  
  z <- cbind(matrix(1,n*m,1),rep(tt,m))
  x <- cbind(matrix(t(x1),m*n,1),matrix(t(x2),m*n,1))
  
  y <-matrix(0,m,n)
  
  for(i in 1:m)
  {
    tt1 <- ttc[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ]
    Gama <- MatDec(tt1,phi1,phi2,struc)
    Sigma <- sigmae*Gama
    b <- matrix(rmvnorm(1,c(0,0),D1),q,1)
    e <- matrix(rmvnorm(1,rep(0,nj[i]),Sigma),nj[i],1)
    x1 <- matrix(x[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ],ncol=p)
    z1 <- matrix(z[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ],ncol=q)
    y[i,] <- x1%*%beta + z1%*%b + times + e
  }
  
  y <- matrix(t(y),m*n,1)
  
  aa <- sort(y,decreasing = FALSE)
  cutoff <- aa[ceiling(percCens*m*n)]
  cc <- matrix(1,m*n,1)*(y <= cutoff)
  y[cc==1] <- cutoff
  LL <- rep(-Inf,length(y))
  LU <- as.vector(y)
  
  cluster <- kronecker(as.matrix(seq(1:m)),as.matrix(rep(1,n)))
  
  Ns <- matrix(0,nrow=length(y),ncol=length(tt)) 
  for(j in 1:length(tt))
  {
    for(i in 1:length(ttc))
    {
      if(ttc[i] == tt[j]){Ns[i,j] <- 1} else {Ns[i,j] <- 0}
    }
  }
  
  return(list(y=y,cc=cc,x=x,z=z,nj=nj,tt=tt,cluster=cluster,cutoff=cutoff,Ns=Ns, LL=LL,LU=LU))
}

## ---------------- ## 
##  Initial values  ##
## ---------------- ## 

inits <- function(y,x,ttc,cluster,init,q) {
  
  dataf <- data.frame(y=y,x=x,ttc = ttc, cluster=cluster)
  
  fitN <- lme(y ~ -1 + x, data=dataf, random = list(cluster=pdSymm(~ ttc)), correlation=corCAR1())
  
  if (missing(init) || is.na(match("beta", names(init) ) ) )
  {
    beta <- as.vector(fixed.effects(fitN));   names(beta) <- NULL
  }
  else {
    beta <- as.vector(init$beta); names(beta) <- NULL
  }
  
  if (missing(init) || is.na(match("sigma", names(init) ) ) )
    sigma2 <- (sigma(fitN))^2
  else
    sigma2 <- init$sigma^2
  
  if (missing(init) || is.na(match("D1", names(init) ) ) )
  {
    dd <- as.matrix(getVarCov(fitN, type = "random.effects")[1:q,1:q])
    D1 <- matrix(dd,q,q,dimnames = NULL)
  }
  else
    D1 <- init$D1
     
  if (missing(init) || is.na(match("phi1", names(init) ) ) )
    phi1 <- abs(as.numeric(coef(fitN$modelStruct$corStruct, unconstrained=FALSE)))
  else
    phi1 <- init$phi1 
  
  return(list(betas=beta,sigma2=sigma2,alphas=D1, phi1=phi1)) 
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
##  SIMULATION STUDY 3 ##
## ------------------- ##

mainDir <- "~./Simulation/Simulation3"

nSim <- 500

tij <- list(t1 = c(2,4,6,8,10,12), t2 = c(2,3,6,9,10,12), t3=c(1,5,9,13,17,21,25,29))

for(k in 1:length(tij)){
  
  subDir <- paste("Cens",0.15,"m",200,"phi",0.6,"t",k,sep = "_")
  CreateDir(mainDir,subDir)
  
  funcSemi <- function(tt) cos(pi*sqrt(tt))
  beta0 <- 2
  beta1 <- -1.5
  betas <- c(beta0,beta1)
  sigmae <- 0.55
  D1 <- matrix(c(0.25,0.1,0.1,0.2),2,2)
  phi1 <- 0.6
  phi2 <- 1
  tt <- as.vector(tij[[k]])
  n <- length(tt)
  m <- 200
  percCens <- 0.15
  ttc <- matrix(rep(tt,m),ncol=1)
  
  set.seed(25)
  x01 <- matrix(runif(n*m,0,1),m,n)
  x02 <- matrix(runif(n*m,-1,2),m,n)
  
  theta_ar <- matrix(0,nSim,(7+length(tt)))
  lambda_ar <- matrix(0,nSim,1)
  dp_Infbf_ar <-  matrix(0,nSim,(length(tt)+2))
  criterio_ar <- matrix(0,nSim,2)
  AIC_criterio_ar <- matrix(0,nSim,1)
  
  for(i in 1:nSim)
  {
    semente <- i*25
    
    dados <- genSemiMix(m=m,x1=x01,x2=x02,tt=tt,beta=betas,sigma=sigmae,D1=D1,funcSemi=funcSemi,semente=semente,percCens=percCens,struc="ar",phi1=phi1,phi2=phi2)
    
    x <-dados$x
    y <-dados$y
    z <-dados$z
    cc <- dados$cc
    nj <-dados$nj
    tt <- dados$tt
    cluster <- dados$cluster
    Ns <- dados$Ns
    LL <- dados$LL
    LU <- dados$LU
    
    initial1 <- NA
    initial1 <- inits(y,x,ttc,cluster,initial1,ncol(as.matrix(z)))
    
    initial <- list(betas=initial1$betas, sigma2=initial1$sigma2, alphas=initial1$alphas, phi1=initial1$phi1, phi2=1, lambda=50)
    
    ## FIT - AR
    
    EST_ar <- nsmec(y=y, cc=cc, x=x, z=z, tt=tt, ttc=ttc, nj=nj, LL=LL, LU=LU, Ns=Ns, initial=initial, struc="ar", lambda.fixed=FALSE, iter.max=300,precision=1e-6)

    theta_ar[i,] <- c(EST_ar$beta1, EST_ar$ff, EST_ar$sigmae,EST_ar$dd,EST_ar$phi1)
    lambda_ar[i,] <-  EST_ar$lambda
    criterio_ar[i,] <- c(EST_ar$loglikp,EST_ar$loglik)
    AIC_criterio_ar[i,] <- c(EST_ar$AICp)
    dp_Infbf_ar[i,] <- sqrt(diag(solve(EST_ar$Infbetasff)))
    
    write.table(theta_ar,file="theta_cens_ar.txt",sep=" ",row.names=F,col.names=F)
    write.table(lambda_ar,file="lambda_cens_ar.txt",sep=" ",row.names=F,col.names=F)
    write.table(criterio_ar,file="criterio_cens_ar.txt",sep=" ",row.names=F,col.names=F)
    write.table(AIC_criterio_ar,file="AICcriterio_cens_ar.txt",sep=" ",row.names=F,col.names=F)
    write.table(dp_Infbf_ar,file="dp_Infbf_cens_ar.txt",sep=" ",row.names=F,col.names=F)
    write.table(i,file="count.txt",sep=" ",row.names=F,col.names=F, append = TRUE)
    
  } #end for nSim
}    #end k




