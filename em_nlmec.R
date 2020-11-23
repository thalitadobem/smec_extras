### ---------------------------------------------------------------- ###
### Censored mixed-effects models for irregularly observed repeated  ###
### measures with applications to HIV viral loads - TEST (2016)      ### 
### ---------------------------------------------------------------- ###                                                   
### Authors: Matos LA, Castro LM and Lachos VH.                      ###
### ---------------------------------------------------------------- ###  

library(mnormt)

## Momentos da Normal Truncada ##

MomemNT <- function(u=c(0,0),S=diag(2),qc=c(1,2)) {
  
  nic=length(u)
  
  if (nic==1) {
    
    qq <- (1/sqrt(S))*(-qc+u)
    R<-1
    alpha <- pnorm(-qq)
    
    dd <- dnorm(-qq)
    H <- qq*dd
    EX <- (1/alpha)*dd   # a vector with a length of nic
    EXX <- 1+1/alpha*H
    varX <- EXX-EX^2
    Eycens <- -sqrt(S)*EX+u
    varyic<- varX*S
    E2yy<-varyic+Eycens^2
    
  }
  
  else {
    
    qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
    R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
    alpha <- pmvnorm(upper=as.vector(-qq), corr=R)
    dd <- rep(0, nic)   #derivative vector
    
    for (j in 1:nic){
      V <- R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu <- -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V)
    }
    
    H <- matrix(rep(0, nic*nic), nrow=nic)
    RH <- matrix(rep(0, nic*nic), nrow=nic)
    
    if(nic==2)     {
      H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }
    
    else {
      for( s in 1:(nic-1)){
        for (t in (s+1):nic){
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V <-  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
          RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
        }
      }
    }
    
    h <- qq*dd-apply(RH, 1, sum)
    diag(H) <- h
    EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
    EXX <- R+1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(S)))%*%EX+u
    varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
    E2yy <- varyic+Eycens%*%t(Eycens)
    
  }
  
  return(list(Ey=Eycens,Eyy=E2yy,Vary=varyic))
  
}


################################################################################
## Agoritmo EM Dados Censurados - Erros não correlacionados 
################################################################################


##  Função do modelo NL usado ##

nlf<-function(t,u,beta1){
  resp<-(exp(beta1[1]+u[1]))+(exp(beta1[2])/(1 + exp((t - exp(beta1[3]))/exp(beta1[4]+u[2]))))
  return(resp) ## transformando os dados em censuras a direita
}

funcH<-function(t,u,beta1){
  
  ni<-length(t)
  resp<-matrix(0,ni,2)
  resp[,1]<- rep(exp(beta1[1]+u[1]),ni)
  resp[,2]<-((t-exp(beta1[3]))*exp(beta1[2]+exp(-beta1[4]-u[2])*(t-exp(beta1[3]))-beta1[4]-u[2]))/((1 + exp((t - exp(beta1[3]))/exp(beta1[4]+u[2])))^2)
  
  return(resp) ## transformando os dados em censuras a direita
}

funcW<-function(t,u,beta1){
  
  ni<-length(t)
  resp<-matrix(0,ni,length(beta1))
  resp[,1]<- rep(exp(beta1[1]+u[1]),ni)
  resp[,2]<- (exp(beta1[2])/(1 + exp((t - exp(beta1[3]))/exp(beta1[4]+u[2]))))
  resp[,3]<- (exp(beta1[2]+exp(-beta1[4]-u[2])*(t-exp(beta1[3]))-beta1[4]-u[2]+beta1[3]))/((1 + exp((t - exp(beta1[3]))/exp(beta1[4]+u[2])))^2)
  resp[,4]<- ((t-exp(beta1[3]))*exp(beta1[2]+exp(-beta1[4]-u[2])*(t-exp(beta1[3]))-beta1[4]-u[2]))/((1 + exp((t - exp(beta1[3]))/exp(beta1[4]+u[2])))^2)
  
  
  return(resp) ## transformando os dados em censuras a direita
  
}


EM_nl <- function(cc,tt,y,x,nj,initial)
{
  
  m<-length(nj)
  N<-sum(nj)
    
    #valores iniciais
    beta1 <- initial$beta 
    sigmae <- initial$sigma2
    D1 <- initial$alphas
    iD1<-solve(D1)
    
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)])
    tetacrit <- c(beta1,sigmae)
    
    q1<-dim(D1)[2]
    p<-length(beta1)  
    m1<-m*p
    m2<-m*q1
    
    ubi=matrix(0,m2,m)  
    
    criterio<-1
    count<-0
    
    while(criterio > 0.0000001){
      
      count <- count + 1
      #print(count)
      soma1<-matrix(0,q1,q1)
      soma2<-0
      soma3<-matrix(0,p,p)
      soma4<-matrix(0,p,1)
      soma5<-matrix(0,p,p) 
      
      ub1<-ubi
      
      ubi=matrix(0,m2,m)     #inicial.b
      ubbi=matrix(0,m2,m2)
      uybi=matrix(0,N,m2)
      uyyi=matrix(0,N,N)
      uyi=matrix(0,N,m)
      xi=matrix(0,N,m1)
      zi=matrix(0,N,m2)
      
      ver<-matrix(0,m,1)
      
      
      for (j in 1:m ){
        
        cc1=cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
        # cc1=c(0,0,0,0,0)
        
        y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
        x1=x[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
        
        W=funcW(x1,(ub1[(((j-1)*q1)+1) : (j*q1), j]),beta1) #X
        H=funcH(x1,(ub1[(((j-1)*q1)+1) : (j*q1), j]),beta1) #Z
        yt=y1-nlf(x1,(ub1[(((j-1)*q1)+1) : (j*q1), j]),beta1)+W%*%beta1+H%*%(ub1[(((j-1)*q1)+1) : (j*q1), j])
        
        x1<-W
        z1<-H
        y1<-yt
        
        gammai=x1%*%beta1
        
        if(sum(cc1)==0){
          
          Psi<-(sigmae*diag(nj[j])+(z1)%*%D1%*%t(z1))
          Psi<-(Psi+t(Psi))/2
          
          delta<- solve(iD1+(t(z1)%*%((z1*(1/sigmae)))))
          uy<- matrix(y1,nj[j],1)
          uyy<- y1%*%t(y1)
          ub<- delta%*%(t(z1)*(1/sigmae))%*%(uy-gammai)
          ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%z1%*%delta)
          uyb<- (uyy-uy%*%t(gammai))%*%(z1*(1/sigmae))%*%delta
          ver[j,]<- dmvnorm(as.vector(y1),gammai,Psi)
          
        }
        
        if(sum(cc1)>=1){
          
          Psi<-(sigmae*diag(nj[j])+(z1)%*%D1%*%t(z1))
          Psi<-(Psi+t(Psi))/2
          
          if(sum(cc1)==nj[j]){
            muc=x1%*%beta1
            Sc<-Psi
            delta<- solve(iD1+(t(z1)%*%((z1*(1/sigmae)))))
            aux<- MomemNT(muc,Sc,y1)
            uy<-aux$Ey
            uyy<- aux$Eyy
            ub<- delta%*%(t(z1)*(1/sigmae))%*%(uy-gammai)
            ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%z1%*%delta)
            uyb<- (uyy-uy%*%t(gammai))%*%(z1*(1/sigmae))%*%delta
            ver[j,]<-pmnorm(y1,as.vector(muc),Sc)
          }
          
          else {
            muc=x1[cc1==1,]%*%beta1+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)
            Sc <-Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]
            delta<- solve(iD1+(t(z1)%*%((z1*(1/sigmae)))))
            aux <-MomemNT(muc,Sc,y1[cc1==1])
            uy <-matrix(y1,nj[j],1)
            uy[cc1==1]<-aux$Ey
            
            uyy<-matrix(0,nj[j],nj[j])
            uyy[cc1==1,cc1==1]<-aux$Vary
            uyy<- uyy+uy%*%t(uy)
            ub<- delta%*%(t(z1)*(1/sigmae))%*%(uy-gammai)
            ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%z1%*%delta)
            uyb<- (uyy-uy%*%t(gammai))%*%(z1*(1/sigmae))%*%delta
            ver[j,]<-dmvnorm(y1[cc1==0],gammai[cc1==0],as.matrix(Psi[cc1==0,cc1==0]))*pmnorm(y1[cc1==1],as.vector(muc),Sc)
          }
          
        }
        
        soma1<- soma1 + ubb
        soma2<- soma2 + (sum(diag(uyy))-t(uy)%*%gammai-t(gammai)%*%uy-sum(diag(t(uyb)%*%z1))-sum(diag(uyb%*%t(z1)))
                         +t(gammai)%*%z1%*%ub+t(ub)%*%t(z1)%*%gammai+t(gammai)%*%gammai+sum(diag(ubb%*%t(z1)%*%z1)))
        soma3<- soma3 + (t(x1)%*%x1)
        soma4<- soma4 + (t(x1)%*%(uy-z1%*%ub))
        soma5<- soma5 + (t(x1)%*%solve(Psi)%*%x1-t(x1)%*%solve(Psi)%*%(uyy-uy%*%t(uy))%*%solve(Psi)%*%x1)
        
        
        uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), (sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<-uyy
        uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), j]<-uy
        uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), (((j-1)*q1)+1) : (j*q1)]<-uyb
        ubi[(((j-1)*q1)+1) : (j*q1), j]<-ub
        ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<-ubb
        zi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), (((j-1)*q1)+1) : (j*q1)]<-z1
        xi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), (((j-1)*p)+1) : (j*p)]<-x1
      }
      
      beta1<- solve(soma3)%*%soma4
      sigmae<- (1/(N))*as.numeric(soma2)
      D1<- (1/(m))*(soma1)
      iD1<-solve(D1)
      
      teta1 <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)])
      teta1crit <- c(beta1,sigmae)
      
      #print(teta1)
      logver <- sum(log(ver))
      #print(logver)
      varbeta<-solve(soma5)
      
      if (count>1){
        #criterio <- sqrt((teta1crit/tetacrit-1)%*%(teta1crit/tetacrit-1))
        criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1)) 
      }
      
      if (count==200){
        criterio <- 0.00000000001
      }
      
      teta<-teta1
      tetacrit <- teta1crit
      logver1<-logver
      
      gamma <- 0
      rho <- 0    
      
  } 
  
  dd<-D1[upper.tri(D1, diag = T)]
  
  npar<-length(c(teta1))
  
  ni<-sum(nj)
  
  loglik<-logver1
  
  AICc<- -2*loglik +2*npar
  AICcorr<- AICc + ((2*npar*(npar+1))/(ni-npar-1))
  BICc <- -2*loglik +log(ni)*npar
  
  obj.out <- list(beta1 = beta1, sigmae= sigmae, rho=rho,gamma=gamma, dd = dd, iter = count, varbeta=varbeta, loglik=loglik, AIC=AICc, BIC=BICc, AICcorr=AICcorr, 
                  xi = xi, zi = zi, ubi = ubi, ubbi = ubbi, uybi = uybi, uyi = uyi, uyyi = uyyi  )
  
  class(obj.out) <- "EM_NCens"
  
  return(obj.out)
  
}

