#########################################################################
######### Option Calibration of InHomogeneous Levy Models #################
#########################################################################
library(cobs)  #Needed for interpolation
library(ggplot2) #plot package
library(zoo)
library(ggExtra)
library(latex2exp)
library(ggthemes)
library(gridExtra)
library(ggmatplot)
library(reshape2)
library(tidyquant) #Used to easily download option prices
library(stringr)
library(pdqr)

#Set seed
set.seed(4)

### Following function provides the continuous logarithm (logc) needed in the function calibration. 
# We want that psiHat(-i)=0 (Martingale condition)
logc <- function(y,k){
  #the continuous logarithm evaluated at [y1,..,yn] trimmed at k
  nargin <- length(as.list(match.call())) -1 #number of input arguments
  if (nargin<2){
    k <- 0
  }
  n <- length(y)
  for (i in 1:n){
    if(abs(y[i]) == 0){
      errorCondition("Don't evaluate the logarithm at 0")
     } else if(abs(y[i])<k^2){
       y[i] <- k*y[i]/abs(y[i])
     } else{
      y[i] <- y[i]
    }
  }
  a <- Re(log(y))
  b <- Im(log(y))
  for(j in 1:(n-1)){
    b[j+1] <- b[j+1] - 2*pi*floor((b[j+1]-b[j])/(2*pi)+ 0.5)
  }
  a+1i*b
}

### Following function provides the Fourier Transform (FT) using the fft (See Cont Tankov p.370) 
# y is the to be transformed vector, and x is the original grid.
FT <- function(x,y,noInverse){
  #Integral bounds
  nargin <- length(as.list(match.call())) -1 #number of input arguments
  u <- rep(0,length(x))
  iv <- rep(0,length(x))
  if (nargin < 3){
    noInverse = TRUE
  }
  N <- length(x)
  A <- x[1]
  B <- x[length(x)]
  Delta <- x[2] - x[1]
  
  #Trapzoidal Rule
  w <- rep(1,N)
  w[1] <- 0.5
  w[length(w)] <- 0.5
    
  y <- w*y
  sig <- 1-2*noInverse # inverse:-1 ; !inverse: 1
  coef <- 1-(1-2*pi)*(!noInverse) # inverse: 1 ; !inverse: 2*pi

  #Negative values
  nn <- seq(N/2,1,-1)
  if (noInverse){ # !inverse to run on negative half
    fyn <- fft(y)} else{
    fyn <- fft(y,inverse=TRUE)}
  un <-  -2*pi*(nn-1)/(N*Delta) # grid points where fourier transform is calculated
  ivn <- (B-A)/N*exp(-1i*un*A*sig)*fyn[nn]
  
  #Positive values
  np <- seq(2,N/2+1,1)
  if (noInverse){
    fyp <- fft(y, inverse=TRUE)} else{
    fyp <- fft(y)}
  up <- 2*pi*(np-1)/(N*Delta)
  ivp <- (B-A)/N*exp(-1i*up*A*sig)*fyp[np]
  
  # concatenate
  u <- c(un,up)
  iv <- c(ivn,ivp)
  # mesh of new grid: 2*pi/(N*Delta)
  list(u =u, fy = iv/coef, x = x, y =y)
}

###Weight Functions for calibration and confidence Intervals

#weights Soehl
wSigma <- function(x){
  cSigma <- -2/(2*(2*s+1)/(2*s+3) - 8*(2*s+3)/(2*s+5) + 12*(2*s+5)/(2*s+7) - 8*(2*s+7)/(2*s+9) + 2*(2*s+9)/(2*s+11))
  cSigma*((2*s+1)*x^(2*s)-4*(2*s+3)*x^(2*s+2)+6*(2*s+5)*x^(2*s+4)-4*(2*s+7)*x^(2*s+6)+(2*s+9)*x^(2*s+8))*(abs(x)<=1)
  }
wGamma <- function(x){
  cGamma <- 1/2*1/(1/(2*s+3)-3/(2*s+5)+3/(2*s+7)-1/(2*s+9))
  cGamma*(x^(2*s+1)-3*x^(2*s+3)+3*x^(2*s+5)-x^(2*s+7))*(abs(x)<=1)
  }
wLambda <- function(x){
  cLambda <- 1/(2*(2*s+3)/(2*s+1) - 8*(2*s+5)/(2*s+3) + 12*(2*s+7)/(2*s+5) - 8*(2*s+9)/(2*s+7) + 2*(2*s+11)/(2*s+9))
  cLambda*((2*s+3)*x^(2*s)-4*(2*s+5)*x^(2*s+2)+6*(2*s+7)*x^(2*s+4)-4*(2*s+9)*x^(2*s+6)+(2*s+11)*x^(2*s+8))*(abs(x)<=1)
  }
flatTopKernel <- function(x){
  ((abs(x)<=0.05)+(abs(x)<1 & abs(x)>0.05)*(exp(-exp(-(abs(x)-0.05)^(-2))/(abs(x)-1)^2)))*(abs(x)<=1)
  }

# #weights Belomestny
# wSigma <- function(x){(s+3)/(1-2^(-2/(s+1)))*abs(x)^s*(1-2*(abs(x)>(2^(-1/(s+1)))))*(abs(x)<=1)}
# wGamma <- function(x){0.5*(s+2)*abs(x)^s*sign(x)*(abs(x)<=1)}
# wLambda <- function(x){(s+1)/(2*(2^(2/(s+3))-1))*abs(x)^s*(1-2*(abs(x)<(2^(-1/(s+3)))))*(abs(x)<=1)}
# flatTopKernel <- function(x){((abs(x)<=0.05)+(abs(x)<1 & abs(x)>0.05)*(exp(-exp(-(abs(x)-0.05)^(-2))/(abs(x)-1)^2)))*(abs(x)<=1)}

# #Plot of Weight functions
# s <- 2
# x <- seq(-1.1,1.1,0.00001)
# ggWSigma <- ggplot(data.frame(x,y=wSigma(x)), aes(x,y)) + geom_line(size=1.5) +theme_bw() + xlab(TeX("$v$"))+ ylab(TeX("$w_{\\sigma_j^2}^{1}(v)$"))+
#   theme(text = element_text(size=20))
# ggWGamma <- ggplot(data.frame(x,y=wGamma(x)), aes(x,y)) + geom_line(size=1.5) +theme_bw() + xlab(TeX("$v$"))+ ylab(TeX("$w_{\\gamma_j}^{1}(v)$"))+
#   theme(text = element_text(size=20))
# ggWLambda <- ggplot(data.frame(x,y=wLambda(x)), aes(x,y)) + geom_line(size=1.5) +theme_bw() + xlab(TeX("$v$"))+ ylab(TeX("$w_{\\lambda_j}^{1}(v)$"))+
#   theme(text = element_text(size=20))
# ggWeights <- arrangeGrob(ggWSigma,ggWGamma,ggWLambda, nrow=1)
# ggsave("weights.jpg",ggWeights, width = 12, height=4)


confIntervals <- function(k,mi){
  if(k== 1){
    if(model == "simulations"){
      TNew <- time[k]
      xNew <- sk1[mi,]
      UNew <- U[mi] #Not really justified....
      #UNew <- 20
      DeltaNew <- median(diff(xNew))
      deltaNew <- stdNoise[mi,]
    } else{
      TNew <- time[k]
      xNew <- sk1
      UNew <- pmin(U[k]) #Not really justified....
      #UNew <- 20
      DeltaNew <- median(diff(xNew))
      deltaNew <- stdNoise1
    }
    fSigmaNew <- 1/(UNew^3)*wSigma(v/UNew)*1i*v*(1+1i*v)/((TNew)*phiNum)
    fGammaNew <- 1/(UNew^2)*wGamma(v/UNew)*1i*v*(1+1i*v)/((TNew)*phiNum)
    fLambdaNew <- 1/UNew*wLambda(v/UNew)*1i*v*(1+1i*v)/((TNew)*phiNum)
    
    FSigmaNew <- FT(v,fSigmaNew,noInverse = FALSE)
    FGammaNew <- FT(v,fGammaNew,noInverse = FALSE)
    FLambdaNew <- FT(v,fLambdaNew,noInverse = FALSE)
    
    valFSigmaNew <- approxfun(FSigmaNew$u,Re(FSigmaNew$fy))(-xNew)+1i*approxfun(FSigmaNew$u,Im(FSigmaNew$fy))(-xNew)
    valFSigmaNew[is.na(valFSigmaNew)] <- 0
    valFGammaNew <- approxfun(FGammaNew$u,Re(FGammaNew$fy))(-xNew)+1i*approxfun(FGammaNew$u,Im(FGammaNew$fy))(-xNew)
    valFGammaNew[is.na(valFGammaNew)] <- 0
    valFLambdaNew <- approxfun(FLambdaNew$u,Re(FLambdaNew$fy))(-xNew)+1i*approxfun(FLambdaNew$u,Im(FLambdaNew$fy))(-xNew)
    valFLambdaNew[is.na(valFLambdaNew)] <- 0
    
    varSigma2 <- 4*pi^2*sum((deltaNew*Re(valFSigmaNew)*DeltaNew)^2)
    varGamma <- 4*pi^2*sum((-Re(valFSigmaNew)+Im(valFGammaNew))^2*deltaNew^2*DeltaNew^2)
    varLambda <- 4*pi^2*sum((-Re(0.5*valFSigmaNew+valFLambdaNew)+Im(valFGammaNew))^2*deltaNew^2*DeltaNew^2)
    
    #Estimate Nu
    if(model == "simulations"){
      UNuNew <- UNu[mi]
      x0 <-x[mi,] 
    } else{
      UNuNew <- UNu[k]
      x0 <-x[k,] 
    }
    fNuNew <- -flatTopKernel(v/UNuNew)*v*(1i+v)/((TNew)*phiNu)
    FNuNew <- FT(v,fNuNew,noInverse = FALSE)
    varNu <- rep(0,N)
    
    g0New <- FT(v,v^0*flatTopKernel(v/UNuNew), noInverse = FALSE)
    valg0New <- approxfun(g0New$u,Re(g0New$fy))(x0)+1i*approxfun(g0New$u,Im(g0New$fy))(x0)
    g1New <- FT(v,v^1*flatTopKernel(v/UNuNew), noInverse = FALSE)
    valg1New <- approxfun(g1New$u,Re(g1New$fy))(x0)+1i*approxfun(g1New$u,Im(g1New$fy))(x0)
    g2New <- FT(v,v^2*flatTopKernel(v/UNuNew), noInverse = FALSE)
    valg2New <- approxfun(g2New$u,Re(g2New$fy))(x0)+1i*approxfun(g2New$u,Im(g2New$fy))(x0)
    
    for(j in 1:N){
    valFNuNew <- approxfun(FNuNew$u,Re(FNuNew$fy))(x0[j]-xNew)+1i*approxfun(FNuNew$u,Im(FNuNew$fy))(x0[j]-xNew) #What is x0?
    valFNuNew[is.na(valFNuNew)] <- 0
    
    #Still an imaginary term, really small though
    varNu[j] <- 4*pi^2*sum((exp(-xNew)/(2*pi)*valFNuNew+Re(valFSigmaNew)*(0.5*valg2New[j]+1i*valg1New[j]-0.5*valg0New[j])+Im(valFSigmaNew)*(-1i*valg1New[j]+valg0New[j])-Re(valFLambdaNew)*valg0New[j])^2*deltaNew^2*DeltaNew^2)
    }
  } else{
    TOld <- time[k-1]
    TNew <- time[k]
    if(model=="simulations"){
      xOld <- sk0[mi,]
      xNew <- sk1[mi,]
      UOld <- pmin(eval(as.name(paste("U",k-1,sep="")))[mi]) #Not really justified...
      UNew <- U[mi] 
      UNuNew <- UNu[mi]
      deltaOld <- eval(as.name(paste("stdNoise",k-1,sep="")))[mi,]
      deltaNew <- stdNoise[mi,]
    } else{
      xOld <- sk0
      xNew <- sk1
      UOld <- pmin(U[k-1]) 
      UNew <- pmin(U[k])
      UNuNew <- UNu[k]
      deltaOld <- stdNoise0
      deltaNew <- stdNoise1
    }
    DeltaOld <- median(diff(xOld))
    DeltaNew <- median(diff(xNew))
    
    fSigmaOld <- 1/(UOld^3)*wSigma(v0/UOld)*1i*v0*(1+1i*v0)/((TNew-TOld)*phiNum0)
    fGammaOld <- 1/(UOld^2)*wGamma(v0/UOld)*1i*v0*(1+1i*v0)/((TNew-TOld)*phiNum0)
    fLambdaOld <- 1/UOld*wLambda(v0/UOld)*1i*v0*(1+1i*v0)/((TNew-TOld)*phiNum0)
    fSigmaNew <- 1/(UNew^3)*wSigma(v/UNew)*1i*v*(1+1i*v)/((TNew-TOld)*phiNum)
    fGammaNew <- 1/(UNew^2)*wGamma(v/UNew)*1i*v*(1+1i*v)/((TNew-TOld)*phiNum)
    fLambdaNew <- 1/UNew*wLambda(v/UNew)*1i*v*(1+1i*v)/((TNew-TOld)*phiNum)
    
    FSigmaOld <- FT(v0,fSigmaOld ,noInverse = FALSE)
    FGammaOld <- FT(v0,fGammaOld,noInverse = FALSE)
    FLambdaOld <- FT(v0,fLambdaOld,noInverse = FALSE)
    FSigmaNew <- FT(v,fSigmaNew,noInverse = FALSE)
    FGammaNew <- FT(v,fGammaNew,noInverse = FALSE)
    FLambdaNew <- FT(v,fLambdaNew,noInverse = FALSE)
    
    valFSigmaOld <- approxfun(FSigmaOld$u,Re(FSigmaOld$fy))(-xOld)+1i*approxfun(FSigmaOld$u,Im(FSigmaOld$fy))(-xOld)
    valFSigmaOld[is.na(valFSigmaOld)] <- 0
    valFGammaOld <- approxfun(FGammaOld$u,Re(FGammaOld$fy))(-xOld)+1i*approxfun(FGammaOld$u,Im(FGammaOld$fy))(-xOld)
    valFGammaOld[is.na(valFGammaOld)] <- 0
    valFLambdaOld <- approxfun(FLambdaOld$u,Re(FLambdaOld$fy))(-xOld)+1i*approxfun(FLambdaOld$u,Im(FLambdaOld$fy))(-xOld)
    valFLambdaOld[is.na(valFLambdaOld)] <- 0
    valFSigmaNew <- approxfun(FSigmaNew$u,Re(FSigmaNew$fy))(-xNew)+1i*approxfun(FSigmaNew$u,Im(FSigmaNew$fy))(-xNew)
    valFSigmaNew[is.na(valFSigmaNew)] <- 0
    valFGammaNew <- approxfun(FGammaNew$u,Re(FGammaNew$fy))(-xNew)+1i*approxfun(FGammaNew$u,Im(FGammaNew$fy))(-xNew)
    valFGammaNew[is.na(valFGammaNew)] <- 0
    valFLambdaNew <- approxfun(FLambdaNew$u,Re(FLambdaNew$fy))(-xNew)+1i*approxfun(FLambdaNew$u,Im(FLambdaNew$fy))(-xNew)
    valFLambdaNew[is.na(valFLambdaNew)] <- 0
    
    varSigma2 <- 4*pi^2*sum((deltaNew*Re(valFSigmaNew)*DeltaNew)^2)+4*pi^2*sum((deltaOld*Re(valFGammaOld)*DeltaOld)^2)
    varGamma <- 4*pi^2*sum((-Re(valFSigmaNew)+Im(valFGammaNew))^2*deltaNew^2*DeltaNew^2) + 4*pi^2*sum((-Re(valFSigmaOld)+Im(valFGammaOld))^2*deltaOld^2*DeltaOld^2)
    varLambda <- 4*pi^2*sum((Re(-0.5*valFSigmaNew+valFLambdaNew)+Im(valFGammaNew))^2*deltaNew^2*DeltaNew^2)+ 4*pi^2*sum((Re(-0.5*valFSigmaOld+valFLambdaOld)+Im(valFGammaOld))^2*deltaOld^2*DeltaOld^2)
    
    #Estimate Nu
    if(model == "simulations"){
      UNuNew <- UNu[mi]
      UNuOld <- eval(as.name(paste("UNu",k-1,sep="")))[mi]
      x0New <-x[mi,] 
      x0Old <- eval(as.name(paste("x",k-1,sep="")))[mi,]
    } else{
      UNuNew <- UNu[k]
      UNuOld <- UNu[k-1]
      x0New <-x[k,]
      x0Old <- x[k-1,]
    }
    
    
    fNuNew <- -flatTopKernel(v/UNuNew)*v*(1i+v)/((TNew-TOld)*phiNum)
    fNuOld <- -flatTopKernel(v0/UNuOld)*v0*(1i+v0)/((TNew-TOld)*phiNum0)
    FNuNew <- FT(v,fNuNew,noInverse = FALSE)
    FNuOld <- FT(v,fNuOld,noInverse = FALSE)
    varNu <- rep(0,N)
    
    g0New <- FT(v,v^0*flatTopKernel(v/UNuNew), noInverse = FALSE)
    g0Old <- FT(v0,v0^0*flatTopKernel(v0/UNuOld), noInverse = FALSE)
    valg0New <- approxfun(g0New$u,Re(g0New$fy))(x0New)+1i*approxfun(g0New$u,Im(g0New$fy))(x0New)
    valg0Old <- approxfun(g0Old$u,Re(g0Old$fy))(x0Old)+1i*approxfun(g0Old$u,Im(g0Old$fy))(x0Old)
    g1New <- FT(v,v^1*flatTopKernel(v/UNuNew), noInverse = FALSE)
    g1Old <- FT(v0,v0^1*flatTopKernel(v/UNuOld), noInverse = FALSE)
    valg1New <- approxfun(g1New$u,Re(g1New$fy))(x0New)+1i*approxfun(g1New$u,Im(g1New$fy))(x0New)
    valg1Old <- approxfun(g1Old$u,Re(g1Old$fy))(x0Old)+1i*approxfun(g1Old$u,Im(g1Old$fy))(x0Old)
    g2New <- FT(v,v^2*flatTopKernel(v/UNuNew), noInverse = FALSE)
    g2Old <- FT(v0,v0^2*flatTopKernel(v0/UNuOld), noInverse = FALSE)
    valg2New <- approxfun(g2New$u,Re(g2New$fy))(x0New)+1i*approxfun(g2New$u,Im(g2New$fy))(x0New)
    valg2Old <- approxfun(g2Old$u,Re(g2Old$fy))(x0Old)+1i*approxfun(g2Old$u,Im(g2Old$fy))(x0Old)
    
    for(j in 1:N){
      valFNuNew <- approxfun(FNuNew$u,Re(FNuNew$fy))(x0New[j]-xNew)+1i*approxfun(FNuNew$u,Im(FNuNew$fy))(x0New[j]-xNew) #What is x0?
      valFNuNew[is.na(valFNuNew)] <- 0
      valFNuOld <- approxfun(FNuOld$u,Re(FNuOld$fy))(x0Old[j]-xOld)+1i*approxfun(FNuOld$u,Im(FNuOld$fy))(x0Old[j]-xOld) #What is x0?
      valFNuOld[is.na(valFNuOld)] <- 0
      
      #Still an imaginary term, really small though
      varNu[j] <- 4*pi^2*sum((exp(-xNew)/(2*pi)*valFNuNew+Re(valFSigmaNew)*(0.5*valg2New[j]+1i*valg1New[j]-0.5*valg0New[j])+Im(valFSigmaNew)*(-1i*valg1New[j]+valg0New[j])-Re(valFLambdaNew)*valg0New[j])^2*deltaNew^2*DeltaNew^2) +
        4*pi^2*sum((exp(-xOld)/(2*pi)*valFNuOld+Re(valFSigmaOld)*(0.5*valg2Old[j]+1i*valg1Old[j]-0.5*valg0Old[j])+Im(valFSigmaOld)*(-1i*valg1Old[j]+valg0Old[j])-Re(valFLambdaOld)*valg0Old[j])^2*deltaOld^2*DeltaOld^2)
        
    }
  }
  list(varSigma2 = varSigma2, varGamma = varGamma, varLambda = varLambda, varNu = Re(varNu))
}

coverageProbability <- function(k,xNu1,xNu2){
  sigma2 <- eval(as.name(paste("sigma",k,sep="")))^2
  sigma2Hat <- eval(as.name(paste("sigmaHat",k,sep="")))^2
  varSigma2Hat <-eval(as.name(paste("varSigma2Hat",k,sep="")))
  gamma <- eval(as.name(paste("gamma",k,sep="")))
  gammaHat <- eval(as.name(paste("gammaHat",k,sep="")))
  varGammaHat <-eval(as.name(paste("varGammaHat",k,sep="")))
  lambda <- eval(as.name(paste("lambda",k,sep="")))
  lambdaHat <- eval(as.name(paste("lambdaHat",k,sep="")))
  varLambdaHat <-eval(as.name(paste("varLambdaHat",k,sep="")))
  
  #Also Levy density at x=0 and x=1/4
  x0 <- eval(as.name(paste("x",k,sep="")))
  Nu0 <- eval(as.name(paste("Nu",k,sep="")))(xNu1)
  nuHat0 <- eval(as.name(paste("nuHat",k,sep="")))[,which.min(x0[1,]< xNu1)]
  varNuHat0 <- eval(as.name(paste("varNuHat",k,sep="")))[,which.min(x0[1,]< xNu1)]
  Nu4 <- eval(as.name(paste("Nu",k,sep="")))(xNu2)
  nuHat4 <- eval(as.name(paste("nuHat",k,sep="")))[,which.min(x0[1,]< xNu2)]
  varNuHat4 <- eval(as.name(paste("varNuHat",k,sep="")))[,which.min(x0[1,]< xNu2)]
  
  counterNu0 <- sum((nuHat0+lowerQ*sqrt(varNuHat0))<=Nu0 & (nuHat0+upperQ*sqrt(varNuHat0))>=Nu0)
  counterNu4 <- sum((nuHat4+lowerQ*sqrt(varNuHat4))<=Nu4 & (nuHat4+upperQ*sqrt(varNuHat4))>=Nu4)
  counterSigma2 <- sum((sigma2Hat+lowerQ*sqrt(varSigma2Hat))<=sigma2 & (sigma2Hat+upperQ*sqrt(varSigma2Hat))>=sigma2)
  counterGamma <- sum((gammaHat+lowerQ*sqrt(varGammaHat))<=gamma & (gammaHat+upperQ*sqrt(varGammaHat))>=gamma)
  counterLambda <- sum((lambdaHat+lowerQ*sqrt(varLambdaHat))<=lambda & (lambdaHat+upperQ*sqrt(varLambdaHat))>=lambda)
  covVector <- c(counterSigma2/length(sigmaHat)*100,counterGamma/length(gammaHat)*100,counterLambda/length(lambdaHat)*100, counterNu0/length(nuHat0)*100, counterNu4/length(nuHat4)*100)
  return(covVector)
}

### Model selection 
modelParameters <- function(r,T,k){
  if(k==1){
    # #Merton Model
    sigma <- 0.3    #diffusion volatility
    lambda <- 5     #jump intensity
    mu <- -0.1        #mean jump size
    delta <- 0.2           #std deviation of jump size
    gamma <- -sigma^2/2 - lambda*(exp(delta^2/2 + mu)-1) #drift -> martingale condition
    Nu <- function(x){lambda/(sqrt(2*pi)*delta)*exp(-(x-mu)^2/(2*delta^2))} #Levy Density
    psi <- function(u){-sigma^2*u^2/2+1i*gamma*u+lambda*(exp(-delta^2*u^2/2+1i*mu*u)-1)} #characteristic exponent
    phi <- function(u){exp(T*psi(u))} #Characteristic function
  } else if(k==2){
    sigma <- 0.1       #diffusion volatility
    lambda <- 2      #jump intensity
    mu <- -0.2        #mean jump size
    delta <- 0.4          #std deviation of jump size
    gamma <- -sigma^2/2 - lambda*(exp(delta^2/2 + mu)-1) #drift -> martingale condition
    Nu <- function(x){lambda/(sqrt(2*pi)*delta)*exp(-(x-mu)^2/(2*delta^2))} #Levy Density
    psi <- function(u){-sigma^2*u^2/2+1i*gamma*u+lambda*(exp(-delta^2*u^2/2+1i*mu*u)-1)} #characteristic exponent
    phi <- function(u){exp(T*psi(u))} #Characteristic function
    # sigma <- 0.25     #diffusion volatility
    # lambda <- 5        #jump intensity
    # lambdaPlus <- 1/0.07 #exponential decay positive x-axis
    # lambdaMin <-  1/0.13  #exponential decay negative x-axis
    # p <- lambdaMin/(lambdaPlus+lambdaMin)      #probability of positive jump, chosen such that we are continuous at x=0
    # gamma <- -sigma^2/2 - lambda*(p/(lambdaPlus-1)-(1-p)/(lambdaMin+1)) #drift -> martingale condition
    # Nu <- function(x){p*lambda*lambdaPlus*exp(-lambdaPlus*x)*(x>=0)+(1-p)*lambda*lambdaMin*exp(-lambdaMin*abs(x))*(x<0)} #Levy Density
    # psi <- function(u){-0.5*sigma^2*u^2+1i*gamma*u+1i*u*lambda*(p/(lambdaPlus-1i*u)-(1-p)/(lambdaMin+1i*u))} #characteristic exponent
    # phi <- function(u){exp(T*psi(u))} #Characteristic function
  } else if(k==3){
    sigma <- 0.2     #diffusion volatility
    lambda <- 3      #jump intensity
    mu <- -0.1      #mean jump size
    delta <- 0.3        #std deviation of jump size
    #delta <- 0.1
    gamma <- -sigma^2/2 - lambda*(exp(delta^2/2 + mu)-1) #drift -> martingale condition
    Nu <- function(x){lambda/(sqrt(2*pi)*delta)*exp(-(x-mu)^2/(2*delta^2))} #Levy Density
    psi <- function(u){-sigma^2*u^2/2+1i*gamma*u+lambda*(exp(-delta^2*u^2/2+1i*mu*u)-1)} #characteristic exponent
    phi <- function(u){exp(T*psi(u))} #Characteristic function
  } else if(k==4){
    sigma <- 0.05     #diffusion volatility
    lambda <- 1     #jump intensity
    mu <- -0.1        #mean jump size
    delta <- 0.2           #std deviation of jump size
    gamma <- -sigma^2/2 - lambda*(exp(delta^2/2 + mu)-1) #drift -> martingale condition
    Nu <- function(x){lambda/(sqrt(2*pi)*delta)*exp(-(x-mu)^2/(2*delta^2))} #Levy Density
    psi <- function(u){-sigma^2*u^2/2+1i*gamma*u+lambda*(exp(-delta^2*u^2/2+1i*mu*u)-1)} #characteristic exponent
    phi <- function(u){exp(T*psi(u))} #Characteristic function
    # sigma <- 0.8    #diffusion volatility
    # lambda <- 8        #jump intensity
    # lambdaPlus <- 1/0.10 #exponential decay positive x-axis
    # lambdaMin <-  1/0.15  #exponential decay negative x-axis
    # p <- lambdaMin/(lambdaPlus+lambdaMin)      #probability of positive jump, chosen such that we are continuous at x=0
    # gamma <- -sigma^2/2 - lambda*(p/(lambdaPlus-1)-(1-p)/(lambdaMin+1)) #drift -> martingale condition
    # Nu <- function(x){p*lambda*lambdaPlus*exp(-lambdaPlus*x)*(x>=0)+(1-p)*lambda*lambdaMin*exp(-lambdaMin*abs(x))*(x<0)} #Levy Density
    # psi <- function(u){-0.5*sigma^2*u^2+1i*gamma*u+1i*u*lambda*(p/(lambdaPlus-1i*u)-(1-p)/(lambdaMin+1i*u))} #characteristic exponent
    # phi <- function(u){exp(T*psi(u))} #Characteristic function
  }else{
    print('Out Of Range') 
  }
  list(sigma = sigma, lambda = lambda, gamma=gamma, phi = phi, Nu = Nu, mu = mu, delta=delta)
}

### Simulation of Option prices with specified model (Merton/Kou)

simulateOptionPrices <- function(r,T,model, linear, sampleSize, noiseLevel, design, Nu){
  #Get initial Option prices to use for simulation
  #Parameters for Fourier Transform
  A <- 2^10 #(See page 370 Cont and Tankov for interpretation of parameters)
  M <- 12
  N <- 2^M
  l <- 0:(N-1)
  Delta <- A/(N-1)
  x <- -A/2 + l*Delta
  Zeta <- function(v,r,phi){
    #Case when v =0:
    eps <- 1e-10
    v <- v + (v==0)*eps
    #value
    exp(1i*v*r*T)*(phi(v-1i)-1)/(1i*v*(1+1i*v))
  }
  y <- Zeta(x,r,phi)
  fxy <- FT(x,y,FALSE) 
  k <- fxy$u 
  stdNoise <- noiseLevel*abs(Re(fxy$fy))
  op <- Re(fxy$fy) #+ pmax(0.,1-exp(k-r*T)) # Call Option prices(equation (11.17) in Cont Tankov)
  nop <- op + rnorm(length(op),0,stdNoise)
  #nop <- op + rnorm(length(op),0,0.005)
  
  M1 <- round(length(op)/2)-which.min(abs(op[1:round(length(op)/2)]-10^(-6)))
  M2 <- which.max(op)+M1
  M1 <- which.max(op)-M1
  
  if(design == "random"){
    prob <- exp(-(k[M2:M1]-r*T)^2)
    prob <- prob/sum(prob)
    ind_sub <- sort(sample(M2:M1,sampleSize,prob=prob))
    sk <- k[ind_sub]
    snop <- nop[ind_sub]
  } else if(design == "deterministic"){
    quantils <- qnorm(c(1:sampleSize)/(sampleSize+1),0,2^(-0.5))
    ind_sub <- sapply(quantils,function(x){which.min(abs(k-r*T-x))})
    sk <- k[ind_sub]
    snop <- nop[ind_sub]
  }
  snop <- snop + pmax(0.,1-exp(sk-r*T)) 
  list(sk = sk, snop = snop, stdNoise = stdNoise[ind_sub])
}


### following function performs the calibration of the model. 
#Input: sk is the log of the strike prices and snop are the corresponding call options prices
#Output: U cut-off value, UNu cut-off value Levy density, sigmaHat estimator of volatility,
#        gammaHat estimator of drift, lambdaHat estimator of intensity, x is grid of Levy density
#        nuHat estimator of Levy density, nuError2 L2 error of nuHat in case of simulation

calibrationHom <- function(sk, snop, mode, T){
  ##Approximation of function O by quadratic spine interpolation
  #Adding points x_0=1 und x_{N+1}=0
  #extrapolation <- 0.005
  extrapolation <- 0.005
  sK <- exp(sk)
  
  extraK <- rep(0,length(sK)+2)
  extraK[1] <- extrapolation * sK[1]
  extraK[2:(length(sK)+1)] <- sK
  extraK[length(sK)+2] <- sK[length(sK)]*1/extrapolation
  
  extraOp <- rep(0,length(sK)+2)
  extraOp[1] <- 1-exp(-r*T)*extraK[1]
  extraOp[2:(length(sK)+1)] <- snop
  extraOp[length(sK)+2] <- 0
  
  BI <- log(extraK[1])   
  BE <- log(extraK[length(extraK)])
  knew <- seq(BI,BE,length=2^12) #only positive values
  
  if(linear){
    OpofK <- approxfun(extraK, extraOp, method = "linear", rule = 2, ties = mean)
    opnew <- OpofK(exp(knew))
  } else{
    css <- cobs(extraK, extraOp, constraint = "convex", nknots = min(100, length(sK)-2), degree =2)
    opnew <- predict(css, exp(knew))[,2]
  }
  
  opnew <- opnew - pmax(0.,1-exp(knew-r*T)) 
  opnew <- pmax(0,opnew)
  
  ##Fourier Transform to estimate psi
  z <- FT(knew-r*T,opnew,TRUE) #x=k-rT, z= F O(v)
  
  v <- z$u
  
  if(model=="simulations"){phiTrue <- testPhi(v-1i)}
  phi <- z$fy*1i*v*(1+1i*v)+1
  
  #temporary plots
  
  #frame()
  #par(mfrow =c(1,2))
  #plot(v[which.max(-65<v):which.max(v>65)],Re(phi[which.max(-65<v):which.max(v>65)]), xlab = TeX('$v$'), ylab = TeX('$Re(\\tilde{\\psi}_1(v))$') ) #These two are not the same
  #grid()
  #if(model=="simulations"){lines(v[which.max(-65<v):which.max(v>65)],Re(phiTrue[which.max(-65<v):which.max(v>65)]), lwd = 2)}
  #plot(v[which.max(-65<v):which.max(v>65)],Im(phi[which.max(-65<v):which.max(v>65)]), xlab = TeX('$v$'), ylab = TeX('$Re(\\tilde{\\psi}_1(v))$') ) #These two are not the same
  #grid()
  #if(model=="simulations"){lines(v[which.max(-65<v):which.max(v>65)],Im(phiTrue[which.max(-65<v):which.max(v>65)]),lwd=2)
  
  psi <- rep(0,length(phi))
  #complex log and use symmetry for the other part
  psi[(length(v)/2):length(v)] <- 1/T*logc(phi[(length(v)/2):length(v)])
  psi[1:(length(v)/2-1)] <- Re(psi[(length(v)-1):(length(v)/2+1)]) -1i*Im(psi[(length(v)-1):(length(v)/2+1)])
  #psiTrue <- 1/T*logc(phiTrue)
  L <- length(v)/2 # place where v has 0 in complex log
  
  # # #temporary plots
  # frame()
  # par(mfrow =c(1,2))
  # plot(v[which.max(-65<v):which.max(v>65)],Re(psi[which.max(-65<v):which.max(v>65)]), xlab = "Fourier Transformed Grid v", ylab = "real part of psi" ) #These two are not the same
  # grid()
  # if(model=="simulations"){lines(v[which.max(-65<v):which.max(v>65)],Re(psiTrue[which.max(-65<v):which.max(v>65)]), lwd = 2)}
  # plot(v[which.max(-65<v):which.max(v>65)],Im(psi[which.max(-65<v):which.max(v>65)]), xlab = "Fourier Transformed Grid v", ylab = "imaginary part of psi" ) #These two are not the same
  # grid()
  # if(model=="simulations"){lines(v[which.max(-65<v):which.max(v>65)],Im(psiTrue[which.max(-65<v):which.max(v>65)]), lwd = 2)}
  
  
  ## Calculation of Estimators for numerous cut-off values U
  if(mode=="oracle" | mode=="flat"){
    cutOffIterations <- 130
  } else if (mode=="PLS"){
    cutOffIterations <- 100 #was first 90, but 30 is better number 
  } else if(mode=="fix"){
    cutOffIterations <- 1}
  
  sigma2HatCur <- rep(0,cutOffIterations)
  gammaHatCur <- rep(0,cutOffIterations)
  lambdaHatCur <- rep(0,cutOffIterations)
  bestError <- -1
  for(i in 1:cutOffIterations){
    if(mode=="fix"){
      UCur <- Ufix
      gridU <- round(Ufix/(v[2]-v[1]))
    } else{
      UCur <- (v[2]-v[1])*(i+10)
      gridU <- i+10
    }
    v_cut <- v[L:(L+gridU)]
    psi_cut <- psi[L:(L+gridU)]
    
    ##Approximation of the estimator integrals with composite trapezoidal rule
    w <- rep(1,length(psi_cut))
    w[1] <- 0.5
    w[length(psi_cut)] <- 0.5 
    
    #weights functions volatility:
    weight <- wSigma(v_cut/UCur) #/U^3 does not need to be done, this will cancel in denominator and numerator.
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight) #Why -2*sum? Normalisation?
    #volatility estimate:
    sigma2HatCur[i] <- sum(w*weight*Re(psi_cut))
    sigma2HatCur[i] <- -2*sigma2HatCur[i]/sum(w*abs(v_cut)^{2}*weight) 
    
    #weight function drift:
    weight <- wGamma(v_cut/UCur) #/U^2 does not need to be done, this will cancel in denominator and numerator.
    #drift estimate:
    gammaHatCur[i] <- sum(w*weight*Im(psi_cut))
    gammaHatCur[i] <- gammaHatCur[i]/sum(w*weight*v_cut)-sigma2HatCur[i]
    
    #weight function intensity:
    weight <- wLambda(v_cut/UCur) #/U does not need to be done, this will cancel in denominator and numerator.
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight*v_cut^2)/v[L+gridU]^2
    lambdaHatCur[i] <- sum(w*Re(psi_cut)*weight)
    lambdaHatCur[i] <- -lambdaHatCur[i]/sum(w*weight)+gammaHatCur[i]+sigma2HatCur[i]/2 
    
    ##Different methods of selecting cut-off U
    #Oracle Method(Real parameter values are known):
    if(mode == "oracle"){
      error <- 1*abs(sigma2HatCur[i]-sigma^2)+1*abs(gammaHatCur[i]-gamma)+1*abs(lambdaHatCur[i]-lambda) #Hij kiest hier de hele tijd i=1
      if((i==1) | error<bestError){
        U <- UCur
        sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
        gammaHat <- gammaHatCur[i]
        lambdaHat <- lambdaHatCur[i]
        bestError <- error
      }
    }
    #Fix Method(U is fixed):
    if(mode=="fix"){
      U <- UCur
      sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
      gammaHat <- gammaHatCur[i]
      lambdaHat <- lambdaHatCur[i]
      #bestError <- error
    }
    #Flat Method(moving average of order 2 of the difference of the discrete curve of the volatility):
    if(mode == "flat"){
      diff <- diff(sigma2HatCur[1:length(sigma2HatCur)])
      ma <- rep(0,length(diff))
      maOrder <- 2
      for(i in (1:length(diff))){
        ma[i]<- sum(abs(diff[max(i-maOrder,1):min(i+maOrder-1,length(diff))]))/(min(i+maOrder,length(diff))-max(i-maOrder,1))
      }
      Us <- (v[2]-v[1])*(10+(1:(cutOffIterations-1)))
      alpha <- 1e-5
      i <- which.min(ma+alpha*Us) # with penalty for high values of U
      U <- Us[i]
      sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
      gammaHat <- gammaHatCur[i]
      lambdaHat <- lambdaHatCur[i]
    }
  }
  
  ## Estimation of Levy Density, al procedures almost analogue as above, i.e., defining complex log, multiple cut-off values etc
  ## only different estimator and usage of flat method
  z <- FT(knew-r*T,exp(-knew+r*T)*opnew,TRUE)
  phiNu <- 1-v*(v+1i)*z$fy
  psiNu <- rep(0,length(phiNu))
  psiNu[(length(v)/2):length(v)] <- 1/T*logc(phiNu[(length(v)/2):length(v)])
  psiNu[1:(length(v)/2-1)] <- Re(psiNu[(length(v)-1):(length(v)/2+1)])-1i*Im(psiNu[(length(v)-1):(length(v)/2+1)])
  
  # Iteration to find best cut-off value for Levy Density
  if(mode == "oracle"| mode == "flat"){
    cutOffIterations <- max(ceiling(U/(v[2]-v[1])),20)-10}
  nuHatCur <- matrix(0,cutOffIterations,length(v))
  nuError2 <- -1 
  if(mode == "oracle"){
    nuError2Cur <- rep(0,cutOffIterations)}
  if(mode == "PLS"){
    objective <- -1}
  for(i in 1:cutOffIterations){
    if(mode == "fix"){
      UCur <- UfixNu
    } else{
      UCur <- (v[2]-v[1])*(i+10)}
    if(mode == "oracle" | mode == "flat" | mode == "fix"){
      lam <- lambdaHat
      fnuHat <- (psiNu-1i*gammaHat*v+lam+sigmaHat^2*v^{2}/2)*flatTopKernel(v/UCur)
    } else if (mode == "PLS"){
      lam <- lambdaHatCur[i]
      fnuHat <- (psiNu-1i*gammaHatCur[i]*v+lam+sigma2HatCur[i]*v^{2}/2)*flatTopKernel(v/UCur)
    }
    fxy_nu <- FT(v,fnuHat,noInverse=FALSE)
    xCur <- fxy_nu$u
    nuHatCur[i,] <- Re(fxy_nu$fy)
    #Assuring that the Density is positive
    wNu <- rep(xCur[2]-xCur[1],length(xCur))
    wNu[1] <- 0.5*wNu[1]
    wNu[length(wNu)] <- 0.5*wNu[length(wNu)]
    
    if(lam>0){
      eq <- function(xi){sum(wNu*(pmax(0,nuHatCur[i,]-xi)))-lam}
      xi0 <- 0
      fxi0 <- eq(xi0)
      xi1 <- 1
      fxi1 <- eq(xi1)
      while(abs(fxi1)>0.001){
        xi2 <- xi1 - (xi1-xi0)/(fxi1-fxi0)*fxi1
        xi0 <- xi1
        xi1 <- max(0,xi2)
        fxi0 <- fxi1
        fxi1 <- eq(xi1)
      }
    } else{
      xi2 <- 0}
    nuHatCur[i,] <- pmax(0,nuHatCur[i,]-xi2)
    if(mode == "oracle"){
      nuTrue <- Nu(xCur) 
      nuError2Cur[i]<-sum(wNu*(nuHatCur[i,]-nuTrue)^2) 
      if(nuError2 < 0 | nuError2Cur[i]<nuError2){
        UNu <- UCur
        x <- xCur
        nuHat <- nuHatCur[i,]
        nuError2 <- nuError2Cur[i]
      }
    }
    if(mode == "fix"){
      UNu <- UCur
      x <- xCur
      nuHat <- nuHatCur[i,]
      if(model!="real_data"){
        nuTrue <- Nu(xCur)
        nuError2 <- sum(wNu*(nuHatCur[i,]-nuTrue)^2)
      }
    }
    if(mode == "PLS"){
      beta <- (0.5*sigma2HatCur[i]+gammaHatCur[i])/(lambdaHatCur[i]-sum(wNu*exp(xCur)*nuHatCur[i,]))
      lambdaHatCur[i] <- beta*lambdaHatCur[i]
      nuHatCur[i,] <- beta*nuHatCur[i,]
      
      vE <- -2^9/2+(0:(2^9-1))*2^9/(2^9-1)
      
      ZetaHat <- function(v,sigma2Hat, gammaHat, lambdaHat, x, nuHat){
        #Case when v =0:
        eps <- 1e-10
        v <- v + (v==0)*eps
        FNu <- FT(x,nuHat*exp(x)) # F(nu(x))(v-i) = F(nu(x)*exp(x))(v)
        valFNu <- splinefun(FNu$u, Re(FNu$fy), method = "fmm")(v)+1i*splinefun(FNu$u, Im(FNu$fy), method = "fmm")(v)
        CFHat <- exp(T*(-sigma2Hat*(v-1i)^2/2+1i*gammaHat*(v-1i)-lambdaHat+valFNu))
        return(exp(1i*v*r*T)*(CFHat-1)/(1i*v*(1+1i*v)))
        #(exp(-r*T)*phi(v-1i)-exp(1i*v*r*T))/(1i*v*(1+1i*v))
        #plot(v,Re(exp(1i*v*r*T)*(CFHat-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
        #plot(v,Re(exp(1i*v*r*T)*(phi1(v-1i)-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
      }
      yE <- sapply(vE, function(v){ZetaHat(v,sigma2HatCur[i],gammaHatCur[i],lambdaHatCur[i],xCur,nuHatCur[i,])})
      #plot(Re(yE))
      fxy <- FT(vE,yE,FALSE)
      kE <- fxy$u
      opE <- Re(fxy$fy)
      #plot(opE)
      kHat <- rep(0,length(sk))
      opHat <- rep(0,length(sk))
      for(j in (1:length(sk))){
        index <- which.min(abs(kE-sk[j]))
        kHat[j] <- kE[index]
        opHat[j] <- opE[index]
      }
      nuDeriv <- diff(nuHatCur[i,],2)
      wE <- rep(xCur[2]-xCur[1],length(nuDeriv))
      wE[1] <- 0.5*wE[1]
      wE[length(wE)] <- 0.5*wE[length(wE)]
      alpha <- 0 #1e-8
      objectiveCur <- sum((opHat-snop+pmax(0.,1-exp(sk-r*T)))^2)+alpha*sum(wE*(nuDeriv/(xCur[2]-xCur[1])^2)^2)
      if(!is.nan(objectiveCur) & (objective < 0 | objectiveCur < objective)){
        U <- UCur
        UNu <- UCur
        sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
        gammaHat <- gammaHatCur[i]
        lambdaHat <- lambdaHatCur[i]
        x <- xCur
        nuHat <- nuHatCur[i,]
        objective <- objectiveCur
        u <- i
        
        
        if(model!="real_data"){
          nuTrue <- Nu(xCur)
          nuError2 <- sum(wNu*(nuHat-nuTrue)^2)
        }
      }
    }
    if(mode == "flat"){
      derivativU <- diff(nuHatCur)/(v[2]-v[1])
      derivativL2 <- rep(0,length(derivativU[,1]))
      for(i in (1:length(derivativU[,1])))
        derivativL2[i] <- sum(wNu*derivativU[i,]^2)
      i <- which.min(derivativL2)
      UNu <- (v[2]-v[1])*(i+10)
      x <- xCur
      nuHat <- nuHatCur[i,]
      if(model!="real_data"){
        nuTrue <- Nu(xCur)
        nuError2 <- sum(wNu*(nuHat-nuTrue)^2)
      }
    }
  }
  # return
  if (model == "simulations"){
    list(U=U, UNu = UNu, sigmaHat = sigmaHat, gammaHat = gammaHat, lambdaHat = lambdaHat, x=x, 
         nuHat = nuHat, nuError2 = nuError2, nuTrue = nuTrue, v=v, phiNum = phi, phiNu = phiNu)
  } else{
    list(U=U, UNu = UNu, sigmaHat = sigmaHat, gammaHat = gammaHat, lambdaHat = lambdaHat, x=x, 
         nuHat = nuHat, nuError2 = nuError2,v=v, phiNum = phi, phiNu = phiNu)
  }
}

calibrationInHom <- function(sk1, snop1,T1, sk0, snop0, T0, mode){
  # 1 are the new option prices and 0 are the previous
  ##Approximation of function O by quadratic spine interpolation
  #Adding points x_0=1 und x_{N+1}=0
  #extrapolation <- 0.005
  extrapolation <- 0.005
  sK1 <- exp(sk1)
  sK0 <- exp(sk0)
  
  extraK1 <- rep(0,length(sK1)+2)
  extraK1[1] <- extrapolation * sK1[1]
  extraK1[2:(length(sK1)+1)] <- sK1
  extraK1[length(sK1)+2] <- sK1[length(sK1)]*1/extrapolation
  
  extraK0 <- rep(0,length(sK0)+2)
  extraK0[1] <- extrapolation * sK0[1]
  extraK0[2:(length(sK0)+1)] <- sK0
  extraK0[length(sK0)+2] <- sK0[length(sK0)]*1/extrapolation
  
  
  extraOp1 <- rep(0,length(sK1)+2)
  extraOp1[1] <- 1-exp(-r*T1)*extraK1[1]
  extraOp1[2:(length(sK1)+1)] <- snop1
  extraOp1[length(sK1)+2] <- 0
  
  extraOp0 <- rep(0,length(sK0)+2)
  extraOp0[1] <- 1-exp(-r*T0)*extraK0[1]
  extraOp0[2:(length(sK0)+1)] <- snop0
  extraOp0[length(sK0)+2] <- 0
  
  BI1 <- log(extraK1[1])   
  BE1 <- log(extraK1[length(extraK1)])
  knew1 <- seq(BI1,BE1,length=2^12) #only positive values
  
  BI0 <- log(extraK0[1])   
  BE0 <- log(extraK0[length(extraK0)])
  knew0 <- seq(BI0,BE0,length=2^12) #only positive values
  
  if(linear){
    OpofK1 <- approxfun(extraK1, extraOp1, method = "linear", rule = 2, ties = mean)
    opnew1 <- OpofK1(exp(knew1))
    OpofK0 <- approxfun(extraK0, extraOp0, method = "linear", rule = 2, ties = mean)
    opnew0 <- OpofK0(exp(knew0))
  } else{
    css1 <- cobs(extraK1, extraOp1, constraint = "convex", nknots = min(100, length(sK1)-2), degree =2)
    opnew1 <- predict(css1, exp(knew1))[,2]
    css0 <- cobs(extraK0, extraOp0, constraint = "convex", nknots = min(100, length(sK0)-2), degree =2)
    opnew0 <- predict(css0, exp(knew0))[,2]
  }
  
  opnew1 <- opnew1 - pmax(0.,1-exp(knew1-r*T1)) 
  opnew1 <- pmax(0,opnew1)
  
  opnew0 <- opnew0 - pmax(0.,1-exp(knew0-r*T0)) 
  opnew0 <- pmax(0,opnew0)
  # plot(knew1,opnew1)
  # plot(knew0, opnew0)
  #Plots of interpolated option prices
  # frame()
  # par(mfrow=c(1,2))
  # plot(sk,snop-pmax(0.,1-exp(sk1-r*T1)), main = 'simulated')
  # grid()
  # plot(knew,opnew, main = 'Interpolated')
  # grid()
  
  ##Fourier Transforms to estimate psi
  z1 <- FT(knew1-r*T1,opnew1,TRUE) #x=k-rT, z= F O(v)
  v <- z1$u
  
  z0 <- FT(knew0-r*T0,opnew0,TRUE) #x=k-rT, z= F O(v)
  v0 <- z0$u
  
  phi1 <- (z1$fy*1i*v*(1+1i*v)+1)
  phi0 <- (z0$fy*1i*v0*(1+1i*v0)+1)
  
  if(model == "simulations"){psiTrue <- 1/(T1-T0)*log(testPhi(v-1i))}
  
  # #temporary plots
  #frame()
  #par(mfrow =c(1,2))
  #plot(v[which.max(-65<v):which.max(v>65)],Re(phiTest[which.max(-65<v):which.max(v>65)]), xlab = "Fourier Transformed Grid v", ylab = "real part of phi" ) #These two are not the same
  #grid()
  # lines(v[which.max(-65<v):which.max(v>65)],Re(phiTrue[which.max(-65<v):which.max(v>65)]), lwd = 2)
  #plot(v[which.max(-65<v):which.max(v>65)],Im(phiTest[which.max(-65<v):which.max(v>65)]), xlab = "Fourier Transformed Grid v", ylab = "imaginary part of phi" ) #These two are not the same
  #grid()
  # lines(v[which.max(-65<v):which.max(v>65)],Im(phiTrue[which.max(-65<v):which.max(v>65)]),lwd=2)
  #dev.off()
  
  psi <- rep(0,length(phi1))
  #complex log and use symmetry for the other part
  psi[(length(v)/2):length(v)] <- 1/(T1-T0)*(logc(phi1[(length(v)/2):length(v)])-logc(phi0[(length(v)/2):length(v)]))
  psi[1:(length(v)/2-1)] <- Re(psi[(length(v)-1):(length(v)/2+1)]) -1i*Im(psi[(length(v)-1):(length(v)/2+1)])
  L <- length(v)/2 # place where v has 0 in complex log
  
  # #temporary plots
  frame()
  par(mfrow =c(1,2))
  plot(v[which.max(-65<v):which.max(v>65)],Re(psi[which.max(-65<v):which.max(v>65)]), xlab = "Fourier Transformed Grid v", ylab = "real part of psi" ) #These two are not the same
  grid()
  if(model=="simulations"){lines(v[which.max(-65<v):which.max(v>65)],Re(psiTrue[which.max(-65<v):which.max(v>65)]), lwd = 2)}
  plot(v[which.max(-65<v):which.max(v>65)],Im(psi[which.max(-65<v):which.max(v>65)]), xlab = "Fourier Transformed Grid v", ylab = "imaginary part of psi" ) #These two are not the same
  grid()
  if(model=="simulations"){lines(v[which.max(-65<v):which.max(v>65)],Im(psiTrue[which.max(-65<v):which.max(v>65)]),lwd=2)}
  
  
  ## Calculation of Estimators for numerous cut-off values U
  if(mode=="oracle" | mode=="flat"){
    cutOffIterations <- 130
  } else if (mode=="PLS"){
    cutOffIterations <- 100 #was 90
  } else if(mode=="fix"){
    cutOffIterations <- 1}
  
  sigma2HatCur <- rep(0,cutOffIterations)
  gammaHatCur <- rep(0,cutOffIterations)
  lambdaHatCur <- rep(0,cutOffIterations)
  bestError <- -1
  for(i in 1:cutOffIterations){
    if(mode=="fix"){
      UCur <- Ufix
      gridU <- round(Ufix/(v[2]-v[1]))
    } else{
      UCur <- (v[2]-v[1])*(i+10)
      gridU <- i+10
    }
    v_cut <- v[L:(L+gridU)]
    psi_cut <- psi[L:(L+gridU)]
    
    ##Approximation of the estimator integrals with composite trapezoidal rule
    w <- rep(1,length(psi_cut))
    w[1] <- 0.5
    w[length(psi_cut)] <- 0.5 
    
    #weights functions volatility:
    weight <- wSigma(v_cut/UCur)
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight) #Why -2*sum? Normalisation?
    #volatility estimate:
    sigma2HatCur[i] <- sum(w*weight*Re(psi_cut))
    sigma2HatCur[i] <- -2*sigma2HatCur[i]/sum(w*abs(v_cut)^{2}*weight) 
    
    #weight function drift:
    weight <- wGamma(v_cut/UCur) 
    #drift estimate:
    gammaHatCur[i] <- sum(w*weight*Im(psi_cut))
    gammaHatCur[i] <- gammaHatCur[i]/sum(w*weight*v_cut)-sigma2HatCur[i]
    
    #weight function intensity:
    weight <- wLambda(v_cut/UCur)
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight*v_cut^2)/v[L+gridU]^2
    lambdaHatCur[i] <- sum(w*Re(psi_cut)*weight)
    lambdaHatCur[i] <- -lambdaHatCur[i]/sum(w*weight)+gammaHatCur[i]+sigma2HatCur[i]/2 
    
    ##Different methods of selecting cut-off U
    #Oracle Method(Real parameter values are known):
    if(mode == "oracle"){
      error <- 1*abs(sigma2HatCur[i]-sigma^2)+1*abs(gammaHatCur[i]-gamma)+1*abs(lambdaHatCur[i]-lambda) #Hij kiest hier de hele tijd i=1
      if((i==1) | error<bestError){
        U <- UCur
        sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
        gammaHat <- gammaHatCur[i]
        lambdaHat <- lambdaHatCur[i]
        bestError <- error
      }
    }
    #Fix Method(U is fixed):
    if(mode=="fix"){
      U <- UCur
      sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
      gammaHat <- gammaHatCur[i]
      lambdaHat <- lambdaHatCur[i]
      #bestError <- error
    }
    #Flat Method(moving average of order 2 of the difference of the discrete curve of the volatility):
    if(mode == "flat"){
      diff <- diff(sigma2HatCur[1:length(sigma2HatCur)])
      ma <- rep(0,length(diff))
      maOrder <- 2
      for(i in (1:length(diff))){
        ma[i]<- sum(abs(diff[max(i-maOrder,1):min(i+maOrder-1,length(diff))]))/(min(i+maOrder,length(diff))-max(i-maOrder,1))
      }
      Us <- (v[2]-v[1])*(10+(1:(cutOffIterations-1)))
      alpha <- 1e-5
      i <- which.min(ma+alpha*Us) # with penalty for high values of U
      U <- Us[i]
      sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
      gammaHat <- gammaHatCur[i]
      lambdaHat <- lambdaHatCur[i]
    }
  }
  
  ## Estimation of Levy Density, al procedures almost analogue as above, i.e., defining complex log, multiple cut-off values etc
  ## only different estimator and usage of flat method
  z1 <- FT(knew1-r*T1,exp(-knew1+r*T1)*opnew1,TRUE)

  z0 <- FT(knew0-r*T0,exp(-knew0+r*T0)*opnew0,TRUE)

  phiNu1 <- 1-v*(v+1i)*z1$fy
  phiNu0 <- 1-v0*(v0+1i)*z0$fy
  psiNu <- rep(0,length(phiNu1))
  psiNu[(length(v)/2):length(v)] <- 1/(T1-T0)*(logc(phiNu1[(length(v)/2):length(v)])-logc(phiNu0[(length(v)/2):length(v)]))
  psiNu[1:(length(v)/2-1)] <- Re(psiNu[(length(v)-1):(length(v)/2+1)])-1i*Im(psiNu[(length(v)-1):(length(v)/2+1)])
  
  # Iteration to find best cut-off value for Levy Density
  if(mode == "oracle"| mode == "flat"){
    cutOffIterations <- max(ceiling(U/(v[2]-v[1])),20)-10}
  nuHatCur <- matrix(0,cutOffIterations,length(v))
  nuError2 <- -1  
  if(mode == "oracle"){
    nuError2Cur <- rep(0,cutOffIterations)}
  if(mode == "PLS"){
    objective <- -1}
  for(i in 1:cutOffIterations){
    if(mode == "fix")
      UCur <- UfixNu
    else
      UCur <- (v[2]-v[1])*(i+10)
    if(mode == "oracle" | mode == "flat" | mode == "fix"){
      lam <- lambdaHat
      fnuHat <- (psiNu-1i*gammaHat*v+lam+sigmaHat^2*v^{2}/2)*flatTopKernel(v/UCur)
    } 
    else if (mode == "PLS"){
      lam <- lambdaHatCur[i]
      fnuHat <- (psiNu-1i*gammaHatCur[i]*v+lam+sigma2HatCur[i]*v^{2}/2)*flatTopKernel(v/UCur)
    }
    fxy_nu <- FT(v,fnuHat,noInverse=FALSE)
    xCur <- fxy_nu$u
    nuHatCur[i,] <- Re(fxy_nu$fy)
    #Assuring that the Density is positive
    wNu <- rep(xCur[2]-xCur[1],length(xCur))
    wNu[1] <- 0.5*wNu[1]
    wNu[length(wNu)] <- 0.5*wNu[length(wNu)]
    
    if(lam>0){
      eq <- function(xi){sum(wNu*(pmax(0,nuHatCur[i,]-xi)))-lam}
      xi0 <- 0
      fxi0 <- eq(xi0)
      xi1 <- 1
      fxi1 <- eq(xi1)
      while(abs(fxi1)>0.001){
        xi2 <- xi1 - (xi1-xi0)/(fxi1-fxi0)*fxi1
        xi0 <- xi1
        xi1 <- max(0,xi2)
        fxi0 <- fxi1
        fxi1 <- eq(xi1)
      }
    }else{
      xi2 <- 0}
    nuHatCur[i,] <- pmax(0,nuHatCur[i,]-xi2)
    if(mode == "oracle"){
      nuTrue <- Nu(xCur) 
      nuError2Cur[i]<-sum(wNu*(nuHatCur[i,]-nuTrue)^2) 
      if(nuError2 < 0 | nuError2Cur[i]<nuError2){
        UNu <- UCur
        x <- xCur
        nuHat <- nuHatCur[i,]
        nuError2 <- nuError2Cur[i]
      }
    }
    if(mode == "fix"){
      UNu <- UCur
      x <- xCur
      nuHat <- nuHatCur[i,]
      if(model!="real_data"){
        nuTrue <- Nu(xCur)
        nuError2 <- sum(wNu*(nuHatCur[i,]-nuTrue)^2)
      }
    }
    if(mode == "PLS"){
      beta <- (0.5*sigma2HatCur[i]+gammaHatCur[i])/(lambdaHatCur[i]-sum(wNu*exp(xCur)*nuHatCur[i,]))
      lambdaHatCur[i] <- beta*lambdaHatCur[i]
      nuHatCur[i,] <- beta*nuHatCur[i,]
      
      vE <- -2^9/2+(0:(2^9-1))*2^9/(2^9-1)
      ZetaHat <- function(v,sigma2Hat, gammaHat, lambdaHat, x, nuHat){
        #Case when v =0:
        eps <- 1e-10
        v <- v + (v==0)*eps
        FNu <- FT(x,nuHat*exp(x)) # F(nu(x))(v-i) = F(nu(x)*exp(x))(v)
        valFNu <- splinefun(FNu$u, Re(FNu$fy), method = "fmm")(v)+1i*splinefun(FNu$u, Im(FNu$fy), method = "fmm")(v)
        valFNu0 <- splinefun(FNu$u, Re(phi0), method = "fmm")(v)+1i*splinefun(FNu$u, Im(phi0), method = "fmm")(v)
        CFHat <- exp((T1-T0)*(-sigma2Hat*(v-1i)^2/2+1i*gammaHat*(v-1i)-lambdaHat+valFNu))*valFNu0 #Not correct yet, we still need the previous phi
        return(exp(1i*v*r*(T1-T0))*(CFHat-1)/(1i*v*(1+1i*v)))
        #(exp(-r*T)*phi(v-1i)-exp(1i*v*r*T))/(1i*v*(1+1i*v))
        #plot(v,Re(exp(1i*v*r*T)*(CFHat-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
        #plot(v,Re(exp(1i*v*r*T)*(phi1(v-1i)-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
      }
      yE <- sapply(vE, function(v){ZetaHat(v,sigma2HatCur[i],gammaHatCur[i],lambdaHatCur[i],xCur,nuHatCur[i,])})
      #plot(Re(yE))
      fxy <- FT(vE,yE,FALSE)
      kE <- fxy$u
      opE <- Re(fxy$fy)
      #plot(opE)
      kHat <- rep(0,length(sk1))
      opHat <- rep(0,length(sk1))
      for(j in (1:length(sk1))){
        index <- which.min(abs(kE-sk1[j]))
        kHat[j] <- kE[index]
        opHat[j] <- opE[index]
      }
      nuDeriv <- diff(nuHatCur[i,],2)
      wE <- rep(xCur[2]-xCur[1],length(nuDeriv))
      wE[1] <- 0.5*wE[1]
      wE[length(wE)] <- 0.5*wE[length(wE)]
      alpha <- 0 #1e-8
      objectiveCur <- sum((opHat-snop1+pmax(0.,1-exp(sk1-r*T1)))^2)+alpha*sum(wE*(nuDeriv/(xCur[2]-xCur[1])^2)^2)
      if(!is.nan(objectiveCur) & (objective < 0 | objectiveCur < objective)){
        U <- UCur
        UNu <- UCur
        sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
        gammaHat <- gammaHatCur[i]
        lambdaHat <- lambdaHatCur[i]
        x <- xCur
        nuHat <- nuHatCur[i,]
        objective <- objectiveCur
        
        O_observed = snop1-pmax(0.,1-exp(sk1-r*T1))
        epsilonGrid = 1/30*max(O_observed)
        lowerB = kHat[which.max(O_observed>epsilonGrid)]
        #lowerB = min(kHat)
        upperB = kHat[length(O_observed)-which.max(rev(O_observed)>epsilonGrid)]
        #upperB = max(kHat)
        gridd <- seq(from = lowerB, to = upperB, by = 0.001)
        O_fun = splinefun(kHat, snop1-pmax(0.,1-exp(sk1-r*T1)))
        opHat_fun = splinefun(kHat, opHat)
        
        #jpeg(paste("O_PLS",k, ".jpeg",sep=""), width = 2, height = 4, units = 'in', res = 600)
        plot(gridd, O_fun(gridd), col = "grey60", xlab = TeX("$x$"), ylab=TeX("Option Function $O(x)$"))
        grid()
        lines(gridd,opHat_fun(gridd), col = "black", lwd = 2)
        #dev.off()
        
        if(model!="real_data"){
          nuTrue <- Nu(xCur)
          nuError2 <- sum(wNu*(nuHat-nuTrue)^2)
        }
      }
    }
    if(mode == "flat"){
      derivativU <- diff(nuHatCur)/(v[2]-v[1])
      derivativL2 <- rep(0,length(derivativU[,1]))
      for(i in (1:length(derivativU[,1])))
        derivativL2[i] <- sum(wNu*derivativU[i,]^2)
      i <- which.min(derivativL2)
      UNu <- (v[2]-v[1])*(i+10)
      x <- xCur
      nuHat <- nuHatCur[i,]
      if(model!="real_data"){
        nuTrue <- Nu(xCur)
        nuError2 <- sum(wNu*(nuHat-nuTrue)^2)
      }
    }
  }
  # return
  if (model == "simulations"){
    list(U=U, UNu = UNu, sigmaHat = sigmaHat, gammaHat = gammaHat, lambdaHat = lambdaHat, x=x, 
         nuHat = nuHat, nuError2 = nuError2, nuTrue = nuTrue, v = v, v0 = v0,phiNum = phi1, phiNum0 = phi0, phiNu = phiNu1, phiNu0 = phiNu0)
  } else{
    list(U=U, UNu = UNu, sigmaHat = sigmaHat, gammaHat = gammaHat, lambdaHat = lambdaHat, x=x, 
         nuHat = nuHat, nuError2 = nuError2, v = v, v0 = v0,phiNum = phi1, phiNum0 = phi0, phiNu = phiNu1, phiNu0 = phiNu0)
  }
}


estNoiseHom <- function(sk, snop, T){
  ##Approximation of function O by quadratic spine interpolation
  #Adding points x_0=1 und x_{N+1}=0
  extrapolation <- 0.005
  #extrapolation <- 0.5
  sK <- exp(sk)
  
  extraK <- rep(0,length(sK)+2)
  extraK[1] <- extrapolation * sK[1]
  extraK[2:(length(sK)+1)] <- sK
  extraK[length(sK)+2] <- sK[length(sK)]*1/extrapolation
  
  extraOp <- rep(0,length(sK)+2)
  extraOp[1] <- 1-exp(-r*T)*extraK[1]
  extraOp[2:(length(sK)+1)] <- snop
  extraOp[length(sK)+2] <- 0
  
  BI <- log(extraK[1])   
  BE <- log(extraK[length(extraK)])
  knew <- seq(BI,BE,length=2^12) #only positive values
  
  if(linear){
    OpofK <- approxfun(extraK, extraOp, method = "linear", rule = 2, ties = mean)
    opnew <- OpofK(exp(knew))
  } else{
    css <- cobs(extraK, extraOp, constraint = "convex", nknots = min(100, length(sK)-2), degree =2)
    opnew <- predict(css, exp(knew))[,2]
  }
  
  opnew <- opnew - pmax(0.,1-exp(knew-r*T))
  opnew <- pmax(0,opnew)
  z <- FT(knew-r*T,opnew,TRUE) #x=k-rT, z= F O(v)
  
  v <- z$u
  phi <- z$fy*1i*v*(1+1i*v)+1
  
  psi <- rep(0,length(phi))
  #complex log and use symmetry for the other part
  psi[(length(v)/2):length(v)] <- 1/T*logc(phi[(length(v)/2):length(v)])
  psi[1:(length(v)/2-1)] <- Re(psi[(length(v)-1):(length(v)/2+1)]) -1i*Im(psi[(length(v)-1):(length(v)/2+1)])
  #psiTrue <- 1/T*logc(phiTrue)
  L <- length(v)/2 # place where v has 0 in complex log
  
  
  cutOffIterations <- 100
  #cutOffIterantions <- 30
  
  sigma2HatCur <- rep(0,cutOffIterations)
  gammaHatCur <- rep(0,cutOffIterations)
  lambdaHatCur <- rep(0,cutOffIterations)
  bestError <- -1
  for(i in 1:cutOffIterations){
    UCur <- (v[2]-v[1])*(i+10)
    gridU <- i+10
    
    v_cut <- v[L:(L+gridU)]
    psi_cut <- psi[L:(L+gridU)]
    
    ##Approximation of the estimator integrals with composite trapezoidal rule
    w <- rep(1,length(psi_cut))
    w[1] <- 0.5
    w[length(psi_cut)] <- 0.5 
    
    #weights functions volatility:
    weight <- wSigma(v_cut/UCur) #/U^3 does not need to be done, this will cancel in denominator and numerator.
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight) 
    #volatility estimate:
    sigma2HatCur[i] <- sum(w*weight*Re(psi_cut))
    sigma2HatCur[i] <- -2*sigma2HatCur[i]/sum(w*abs(v_cut)^{2}*weight) 
    
    #weight function drift:
    weight <- wGamma(v_cut/UCur) #/U^2 does not need to be done, this will cancel in denominator and numerator.
    #drift estimate:
    gammaHatCur[i] <- sum(w*weight*Im(psi_cut))
    gammaHatCur[i] <- gammaHatCur[i]/sum(w*weight*v_cut)-sigma2HatCur[i]
    
    #weight function intensity:
    weight <- wLambda(v_cut/UCur) #/U does not need to be done, this will cancel in denominator and numerator.
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight*v_cut^2)/v[L+gridU]^2
    lambdaHatCur[i] <- sum(w*Re(psi_cut)*weight)
    lambdaHatCur[i] <- -lambdaHatCur[i]/sum(w*weight)+gammaHatCur[i]+sigma2HatCur[i]/2 
  }
  
  ## Estimation of Levy Density, al procedures almost analogue as above, i.e., defining complex log, multiple cut-off values etc
  ## only different estimator and usage of flat method
  z <- FT(knew-r*T,exp(-knew+r*T)*opnew,TRUE)
  phiNu <- 1-v*(v+1i)*z$fy
  psiNu <- rep(0,length(phiNu))
  psiNu[(length(v)/2):length(v)] <- 1/T*logc(phiNu[(length(v)/2):length(v)])
  psiNu[1:(length(v)/2-1)] <- Re(psiNu[(length(v)-1):(length(v)/2+1)])-1i*Im(psiNu[(length(v)-1):(length(v)/2+1)])
  
  nuHatCur <- matrix(0,cutOffIterations,length(v))
  nuError2 <- -1 
  objective <- -1
  for(i in 1:cutOffIterations){
    UCur <- (v[2]-v[1])*(i+10)
    lam <- lambdaHatCur[i]
    fnuHat <- (psiNu-1i*gammaHatCur[i]*v+lam+sigma2HatCur[i]*v^{2}/2)*flatTopKernel(v/UCur)
    
    fxy_nu <- FT(v,fnuHat,noInverse=FALSE)
    xCur <- fxy_nu$u
    nuHatCur[i,] <- Re(fxy_nu$fy)
    #Assuring that the Density is positive
    wNu <- rep(xCur[2]-xCur[1],length(xCur))
    wNu[1] <- 0.5*wNu[1]
    wNu[length(wNu)] <- 0.5*wNu[length(wNu)]
    
    if(lam>0){
      eq <- function(xi){sum(wNu*(pmax(0,nuHatCur[i,]-xi)))-lam}
      xi0 <- 0
      fxi0 <- eq(xi0)
      xi1 <- 1
      fxi1 <- eq(xi1)
      while(abs(fxi1)>0.001){
        xi2 <- xi1 - (xi1-xi0)/(fxi1-fxi0)*fxi1
        xi0 <- xi1
        xi1 <- max(0,xi2)
        fxi0 <- fxi1
        fxi1 <- eq(xi1)
      }
    } else{
      xi2 <- 0}
    nuHatCur[i,] <- pmax(0,nuHatCur[i,]-xi2)
    beta <- (0.5*sigma2HatCur[i]+gammaHatCur[i])/(lambdaHatCur[i]-sum(wNu*exp(xCur)*nuHatCur[i,]))
    lambdaHatCur[i] <- beta*lambdaHatCur[i]
    nuHatCur[i,] <- beta*nuHatCur[i,]
    vE <- -2^9/2+(0:(2^9-1))*2^9/(2^9-1)
      
    ZetaHat <- function(v,sigma2Hat, gammaHat, lambdaHat, x, nuHat){
        #Case when v =0:
      eps <- 1e-10
      v <- v + (v==0)*eps
      FNu <- FT(x,nuHat*exp(x)) # F(nu(x))(v-i) = F(nu(x)*exp(x))(v)
      valFNu <- splinefun(FNu$u, Re(FNu$fy), method = "fmm")(v)+1i*splinefun(FNu$u, Im(FNu$fy), method = "fmm")(v)
      CFHat <- exp(T*(-sigma2Hat*(v-1i)^2/2+1i*gammaHat*(v-1i)-lambdaHat+valFNu))
      return(exp(1i*v*r*T)*(CFHat-1)/(1i*v*(1+1i*v)))
      #(exp(-r*T)*phi(v-1i)-exp(1i*v*r*T))/(1i*v*(1+1i*v))
      #plot(v,Re(exp(1i*v*r*T)*(CFHat-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
      #plot(v,Re(exp(1i*v*r*T)*(phi1(v-1i)-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
    }
    yE <- sapply(vE, function(v){ZetaHat(v,sigma2HatCur[i],gammaHatCur[i],lambdaHatCur[i],xCur,nuHatCur[i,])})
    #plot(Re(yE))
    fxy <- FT(vE,yE,FALSE)
    kE <- fxy$u
    opE <- Re(fxy$fy)
    #plot(opE)
    kHat <- rep(0,length(sk))
    opHat <- rep(0,length(sk))
    for(j in (1:length(sk))){
      index <- which.min(abs(kE-sk[j]))
      kHat[j] <- kE[index]
      opHat[j] <- opE[index]
    }
    
    nuDeriv <- diff(nuHatCur[i,],2)
    wE <- rep(xCur[2]-xCur[1],length(nuDeriv))
    wE[1] <- 0.5*wE[1]
    wE[length(wE)] <- 0.5*wE[length(wE)]
    alpha <- 0 #1e-8
    objectiveCur <- sum((opHat-snop+pmax(0.,1-exp(sk-r*T)))^2)+alpha*sum(wE*(nuDeriv/(xCur[2]-xCur[1])^2)^2)
    #print(objectiveCur)
    if(!is.nan(objectiveCur) & (objective < 0 | objectiveCur < objective)){
      U <- UCur
      UNu <- UCur
      sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
      gammaHat <- gammaHatCur[i]
      lambdaHat <- lambdaHatCur[i]
      x <- xCur
      nuHat <- nuHatCur[i,]
      objective <- objectiveCur
      u <- i
      
      O_observed = snop-pmax(0.,1-exp(sk-r*T))
      #epsilonGrid = 1/30*max(O_observed) #1/50
      epsilonGrid = 1/30*max(O_observed) #1/50
      lowerB = kHat[which.max(O_observed>epsilonGrid)]
      #lowerB = min(kHat)
      upperB = kHat[length(kHat)-which.max(rev(O_observed)>epsilonGrid)]
      #upperB = max(kHat)
      gridd <- seq(from = lowerB, to = upperB, by = 0.001)
      O_fun = splinefun(kHat, snop-pmax(0.,1-exp(sk-r*T)))
      opHat_fun = splinefun(kHat, opHat)
      
      #jpeg("O_PLS_emp1.jpeg", width = 2, height = 4, units = 'in', res = 600)
      plot(gridd, O_fun(gridd), col = "grey60", xlab = TeX("$x$"), ylab=TeX("Option Function $O(x)$"))
      grid()
      lines(gridd,opHat_fun(gridd), col = "black", lwd = 2)
      #dev.off()
      
      residuals2 = function(x){(((O_fun)(x)-opHat_fun(x))/opHat_fun(x))^2}
      plot(gridd,residuals2(gridd))
    
      estNoiseLevel <-sqrt(mean(residuals2(gridd)))
    }
  }
  # return
  print(estNoiseLevel)
  print(c(U,sigmaHat,gammaHat,lambdaHat))
  estNoiseLevel*opHat
}

estNoiseInHom <- function(sk1, snop1,T1, sk0, snop0, T0){
  # 1 are the new option prices and 0 are the previous
  ##Approximation of function O by quadratic spine interpolation
  #Adding points x_0=1 und x_{N+1}=0
  #extrapolation <- 0.005
  extrapolation <- 0.005
  sK1 <- exp(sk1)
  sK0 <- exp(sk0)
  
  extraK1 <- rep(0,length(sK1)+2)
  extraK1[1] <- extrapolation * sK1[1]
  extraK1[2:(length(sK1)+1)] <- sK1
  extraK1[length(sK1)+2] <- sK1[length(sK1)]*1/extrapolation
  
  extraK0 <- rep(0,length(sK0)+2)
  extraK0[1] <- extrapolation * sK0[1]
  extraK0[2:(length(sK0)+1)] <- sK0
  extraK0[length(sK0)+2] <- sK0[length(sK0)]*1/extrapolation
  
  
  extraOp1 <- rep(0,length(sK1)+2)
  extraOp1[1] <- 1-exp(-r*T1)*extraK1[1]
  extraOp1[2:(length(sK1)+1)] <- snop1
  extraOp1[length(sK1)+2] <- 0
  
  extraOp0 <- rep(0,length(sK0)+2)
  extraOp0[1] <- 1-exp(-r*T0)*extraK0[1]
  extraOp0[2:(length(sK0)+1)] <- snop0
  extraOp0[length(sK0)+2] <- 0
  
  BI1 <- log(extraK1[1])   
  BE1 <- log(extraK1[length(extraK1)])
  knew1 <- seq(BI1,BE1,length=2^12) #only positive values
  
  BI0 <- log(extraK0[1])   
  BE0 <- log(extraK0[length(extraK0)])
  knew0 <- seq(BI0,BE0,length=2^12) #only positive values
  
  if(linear){
    OpofK1 <- approxfun(extraK1, extraOp1, method = "linear", rule = 2, ties = mean)
    opnew1 <- OpofK1(exp(knew1))
    OpofK0 <- approxfun(extraK0, extraOp0, method = "linear", rule = 2, ties = mean)
    opnew0 <- OpofK0(exp(knew0))
  } else{
    css1 <- cobs(extraK1, extraOp1, constraint = "convex", nknots = min(100, length(sK1)-2), degree =2)
    opnew1 <- predict(css1, exp(knew1))[,2]
    css0 <- cobs(extraK0, extraOp0, constraint = "convex", nknots = min(100, length(sK0)-2), degree =2)
    opnew0 <- predict(css0, exp(knew0))[,2]
  }
  
  opnew1 <- opnew1 - pmax(0.,1-exp(knew1-r*T1)) 
  opnew1 <- pmax(0,opnew1)
  
  opnew0 <- opnew0 - pmax(0.,1-exp(knew0-r*T0)) 
  opnew0 <- pmax(0,opnew0)
  
  ##Fourier Transforms to estimate psi
  z1 <- FT(knew1-r*T1,opnew1,TRUE) #x=k-rT, z= F O(v)
  v <- z1$u
  
  z0 <- FT(knew0-r*T0,opnew0,TRUE) #x=k-rT, z= F O(v)
  v0 <- z0$u
  
  phi1 <- (z1$fy*1i*v*(1+1i*v)+1)
  phi0 <- (z0$fy*1i*v0*(1+1i*v0)+1)
  psi <- rep(0,length(phi1))
  #complex log and use symmetry for the other part
  psi[(length(v)/2):length(v)] <- 1/(T1-T0)*(logc(phi1[(length(v)/2):length(v)])-logc(phi0[(length(v)/2):length(v)]))
  psi[1:(length(v)/2-1)] <- Re(psi[(length(v)-1):(length(v)/2+1)]) -1i*Im(psi[(length(v)-1):(length(v)/2+1)])
  L <- length(v)/2 # place where v has 0 in complex log
  
  cutOffIterations <- 100
  
  sigma2HatCur <- rep(0,cutOffIterations)
  gammaHatCur <- rep(0,cutOffIterations)
  lambdaHatCur <- rep(0,cutOffIterations)
  bestError <- -1
  for(i in 1:cutOffIterations){
    UCur <- (v[2]-v[1])*(i+10)
    gridU <- i+10
    v_cut <- v[L:(L+gridU)]
    psi_cut <- psi[L:(L+gridU)]
    
    ##Approximation of the estimator integrals with composite trapezoidal rule
    w <- rep(1,length(psi_cut))
    w[1] <- 0.5
    w[length(psi_cut)] <- 0.5 
    
    #weights functions volatility:
    weight <- wSigma(v_cut/UCur)
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight) #Why -2*sum? Normalisation?
    #volatility estimate:
    sigma2HatCur[i] <- sum(w*weight*Re(psi_cut))
    sigma2HatCur[i] <- -2*sigma2HatCur[i]/sum(w*abs(v_cut)^{2}*weight) 
    
    #weight function drift:
    weight <- wGamma(v_cut/UCur) 
    #drift estimate:
    gammaHatCur[i] <- sum(w*weight*Im(psi_cut))
    gammaHatCur[i] <- gammaHatCur[i]/sum(w*weight*v_cut)-sigma2HatCur[i]
    
    #weight function intensity:
    weight <- wLambda(v_cut/UCur)
    weight[length(weight)] <- weight[length(weight)]-2*sum(w*weight*v_cut^2)/v[L+gridU]^2
    lambdaHatCur[i] <- sum(w*Re(psi_cut)*weight)
    lambdaHatCur[i] <- -lambdaHatCur[i]/sum(w*weight)+gammaHatCur[i]+sigma2HatCur[i]/2 
  }
  
  ## Estimation of Levy Density, al procedures almost analogue as above, i.e., defining complex log, multiple cut-off values etc
  ## only different estimator and usage of flat method
  z1 <- FT(knew1-r*T1,exp(-knew1+r*T1)*opnew1,TRUE)
  
  z0 <- FT(knew0-r*T0,exp(-knew0+r*T0)*opnew0,TRUE)
  
  phiNu1 <- 1-v*(v+1i)*z1$fy
  phiNu0 <- 1-v0*(v0+1i)*z0$fy
  psiNu <- rep(0,length(phiNu1))
  psiNu[(length(v)/2):length(v)] <- 1/(T1-T0)*(logc(phiNu1[(length(v)/2):length(v)])-logc(phiNu0[(length(v)/2):length(v)]))
  psiNu[1:(length(v)/2-1)] <- Re(psiNu[(length(v)-1):(length(v)/2+1)])-1i*Im(psiNu[(length(v)-1):(length(v)/2+1)])
  
  # Iteration to find best cut-off value for Levy Density
  nuHatCur <- matrix(0,cutOffIterations,length(v))
  nuError2 <- -1  #nuError2 seems not to work in case of flat
  objective <- -1
  for(i in 1:cutOffIterations){
    UCur <- (v[2]-v[1])*(i+10)
    lam <- lambdaHatCur[i]
    fnuHat <- (psiNu-1i*gammaHatCur[i]*v+lam+sigma2HatCur[i]*v^{2}/2)*flatTopKernel(v/UCur)
    fxy_nu <- FT(v,fnuHat,noInverse=FALSE)
    xCur <- fxy_nu$u
    nuHatCur[i,] <- Re(fxy_nu$fy)
    #Assuring that the Density is positive
    wNu <- rep(xCur[2]-xCur[1],length(xCur))
    wNu[1] <- 0.5*wNu[1]
    wNu[length(wNu)] <- 0.5*wNu[length(wNu)]
    
    if(lam>0){
      eq <- function(xi){sum(wNu*(pmax(0,nuHatCur[i,]-xi)))-lam}
      xi0 <- 0
      fxi0 <- eq(xi0)
      xi1 <- 1
      fxi1 <- eq(xi1)
      while(abs(fxi1)>0.001){
        xi2 <- xi1 - (xi1-xi0)/(fxi1-fxi0)*fxi1
        xi0 <- xi1
        xi1 <- max(0,xi2)
        fxi0 <- fxi1
        fxi1 <- eq(xi1)
      }
    }else{
      xi2 <- 0}
    nuHatCur[i,] <- pmax(0,nuHatCur[i,]-xi2)
    beta <- (0.5*sigma2HatCur[i]+gammaHatCur[i])/(lambdaHatCur[i]-sum(wNu*exp(xCur)*nuHatCur[i,]))
    lambdaHatCur[i] <- beta*lambdaHatCur[i]
    nuHatCur[i,] <- beta*nuHatCur[i,]
    
    vE <- -2^9/2+(0:(2^9-1))*2^9/(2^9-1)
    ZetaHat <- function(v,sigma2Hat, gammaHat, lambdaHat, x, nuHat){
      #Case when v =0:
      eps <- 1e-10
      v <- v + (v==0)*eps
      FNu <- FT(x,nuHat*exp(x)) # F(nu(x))(v-i) = F(nu(x)*exp(x))(v)
      valFNu <- splinefun(FNu$u, Re(FNu$fy), method = "fmm")(v)+1i*splinefun(FNu$u, Im(FNu$fy), method = "fmm")(v)
      valFNu0 <- splinefun(FNu$u, Re(phi0), method = "fmm")(v)+1i*splinefun(FNu$u, Im(phi0), method = "fmm")(v)
      CFHat <- exp((T1-T0)*(-sigma2Hat*(v-1i)^2/2+1i*gammaHat*(v-1i)-lambdaHat+valFNu))*valFNu0 #Not correct yet, we still need the previous phi
      return(exp(1i*v*r*(T1-T0))*(CFHat-1)/(1i*v*(1+1i*v)))
      #(exp(-r*T)*phi(v-1i)-exp(1i*v*r*T))/(1i*v*(1+1i*v))
      #plot(v,Re(exp(1i*v*r*T)*(CFHat-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
      #plot(v,Re(exp(1i*v*r*T)*(phi1(v-1i)-1)/(1i*v*(1+1i*v))), ylim = c(0,0.004))
    }
    yE <- sapply(vE, function(v){ZetaHat(v,sigma2HatCur[i],gammaHatCur[i],lambdaHatCur[i],xCur,nuHatCur[i,])})
    #plot(Re(yE))
    fxy <- FT(vE,yE,FALSE)
    kE <- fxy$u
    opE <- Re(fxy$fy)
    #plot(opE)
    kHat <- rep(0,length(sk1))
    opHat <- rep(0,length(sk1))
    for(j in (1:length(sk1))){
      index <- which.min(abs(kE-sk1[j]))
      kHat[j] <- kE[index]
      opHat[j] <- opE[index]
    }
    nuDeriv <- diff(nuHatCur[i,],2)
    wE <- rep(xCur[2]-xCur[1],length(nuDeriv))
    wE[1] <- 0.5*wE[1]
    wE[length(wE)] <- 0.5*wE[length(wE)]
    alpha <- 0 #1e-8
    objectiveCur <- sum((opHat-snop1+pmax(0.,1-exp(sk1-r*T1)))^2)+alpha*sum(wE*(nuDeriv/(xCur[2]-xCur[1])^2)^2)
    
    if(!is.nan(objectiveCur) & (objective < 0 | objectiveCur < objective)){
      U <- UCur
      UNu <- UCur
      sigmaHat <- sqrt(pmax(0,sigma2HatCur[i]))
      gammaHat <- gammaHatCur[i]
      lambdaHat <- lambdaHatCur[i]
      x <- xCur
      nuHat <- nuHatCur[i,]
      objective <- objectiveCur
      u <- i 
      
      O_observed = snop1-pmax(0.,1-exp(sk1-r*T1))
      epsilonGrid = 1/30*max(O_observed)
      lowerB = kHat[which.max(O_observed>epsilonGrid)]
      #lowerB = min(kHat)
      upperB = kHat[length(O_observed)-which.max(rev(O_observed)>epsilonGrid)]
      #upperB = max(kHat)
      gridd <- seq(from = lowerB, to = upperB, by = 0.001)
      O_fun = splinefun(kHat, snop1-pmax(0.,1-exp(sk1-r*T1)))
      opHat_fun = splinefun(kHat, opHat)
      
      #jpeg(paste("O_PLS_emp",k, ".jpeg",sep=""), width = 2, height = 4, units = 'in', res = 600)
      plot(gridd, O_fun(gridd), col = "grey60", xlab = TeX("$x$"), ylab=TeX("Option Function $O(x)$"))
      grid()
      lines(gridd,opHat_fun(gridd), col = "black", lwd = 2)
      #dev.off()
      
      residuals2 = function(x){(((O_fun)(x)-opHat_fun(x))/opHat_fun(x))^2}
      plot(gridd,residuals2(gridd))
  
      estNoiseLevel <-sqrt(mean(residuals2(gridd)))
      
      #estNoise <- sqrt((opHat-snop1+pmax(0.,1-exp(sk1-r*T1)))^2)
      #O_PLS = snop1-pmax(0.,1-exp(sk1-r*T1))
      #plot(opHat)
      #plot(O_PLS)
      #estNoiseLevel <-sqrt(1/length(opHat)*sum(((opHat-O_PLS)/O_PLS)^2))
      #estNoiseLevel <-sqrt(1/length(opHat)*sum(((opHat-snop1+pmax(0.,1-exp(sk1-r*T1)))/opHat)^2))
    }
  }
  # return
  print(estNoiseLevel)
  print(c(U,sigmaHat,gammaHat,lambdaHat))
  estNoiseLevel*opHat
}

generatePaths <- function(r,NoOfSteps){
  if(model == "simulations"){
    #we take the median of all monteCarlo simulations, this is less affected by outlier than the expected value
    gammaHat <- c(median(gammaHat1),median(gammaHat2), median(gammaHat3))
    sigmaHat <- c(median(sigmaHat1),median(sigmaHat2), median(sigmaHat3))
    lambdaHat <- c(median(lambdaHat1),median(lambdaHat2), median(lambdaHat3))
    nuHat <- rbind(nuHat1[1,],nuHat2[1,],nuHat3[1,])
    x <- rbind(x1[1,],x2[1,],x3[1,])
    #How to choose right nuHat from MC simulations?
    T <- c(0,time)
    S <- 100
  }
  if(model == "real_data"){
    T <- c(0,as.numeric(difftime(expiration, quoteDate, units = "weeks"))/52)
  }
  dt <- T[length(T)]/NoOfSteps
  t <- c(0)
  X <- c(0)
  for(k in 1:(length(T)-1)){
    tBegin <- T[k]
    tEnd <- T[k+1]
    deltaT <- tEnd - tBegin
    tNew <- seq(0,deltaT,dt)
    G <- rnorm(NoOfSteps, mean = 0.0, sd = 1.0)
    Z <- (G - mean(G)) / sqrt(var(G)) #Check if indeed N(0,1)
    G <- sigmaHat[k]*sqrt(dt) * Z
    N <- rpois(1,pmax(0,lambdaHat[k])*deltaT)
    U <- runif(N,min = 0, max = deltaT)
    rDist <- new_r(exp(-x[k,])*nuHat[k,]/lambdaHat[k], type="continuous") #From exponential weighted to normal weighted
    Y <- rDist(N)
    XNew <- rep(0,length(tNew))
    for (i in 1:(length(tNew))){
      XNew[i] <- gammaHat[k]*tNew[i]+sum(G[1:i])+sum((U<tNew[i])*Y)
    }
    XNew <- X[length(X)]+XNew
    t <- c(t,tBegin+tNew)
    X <- c(X,XNew)
  }
  St <- S*exp(r*t+X)
  return(list(t = t, St = St))
}


### Main Calculation

#Model specification, choose simulations for simulations and real_data for real data
model <- "simulations"  

#Interpolation Parameters
linear <- TRUE  #Linear = TRUE gives linear interpolation of O, otherwise quadratic B-splines
design <- "deterministic" 

#Confidence Interval Parameters
alpha <- 0.05
upperQ <- qnorm(1-alpha/2)
lowerQ <- qnorm(alpha/2)

### Simulations using Merton or Kou model
if(model == "simulations"){
  # We are going to make a inhomogeneous model of (Merton, Kou, Merton,Kou) with all different parameters
  #Calibration parameters
  mode <- "flat"   #oracle only usable for simulations, other choices: "flat","pls","fix"
  noise <- "real" #Noise = "estimated" gives PLS noise, noise="real" gives real noise (only with simulations)
  # Ufix <- 23
  # UfixNu <-8
  
  #Model Parameter, we are going to make a model with (Merton,Kou,Merton,Kou)
  r <- 0.06
  maturity <- 12/365 #every week a maturity
  T1 <- maturity*1/3
  T2 <- maturity*2/3
  T3 <- maturity*3/3
  #T4 <- maturity*4/4
  time <- c(T1,T2,T3)
  
  
  #Simulation parameters
  noiseLevel <- 0.010 #in paper 0.015<noiseLevel<0.03
  sampleSize <- 100
  monteCarlo <- 100
  N <- 4096
  
  #Undersmoothing 
  undersmoothing <- FALSE #Turn on or off
  zeta <- 1.0 #Undersmoothing parameter, 4/3, 3/2, 1.35
  
  # looping over the different models
  for(k in 1:3){
    #appropriate smoothness for Kou and Merton models
    if(k==1 | k==3){s <- 6}
    if(k==2 | k==4){s <- 6}
    
    
    if(k==1){modelList <- modelParameters(r,time[k],k)
    } else{modelList <- modelParameters(r,time[k]-time[k-1],k)}
    sigma <- modelList$sigma
    lambda <- modelList$lambda
    gamma <- modelList$gamma
    Nu <- modelList$Nu
    mu <- modelList$mu
    delta <- modelList$delta
    
    assign(paste("sigma",k,sep=""),sigma)
    assign(paste("lambda",k,sep=""),lambda)
    assign(paste("gamma",k,sep=""),gamma)
    assign(paste("Nu",k,sep=""),Nu)
    
    Multiply=function(a,b){
      force(a)
      force(b)
      function(x){a(x)*b(x)}
    }
    if(k==1){phi1 <- modelList$phi 
    phi <- phi1
    } else if(k==2){phi2 <- modelList$phi
    phi <- Multiply(phi1,phi2)
    } else if(k==3){phi3 <- modelList$phi
    phi <- Multiply(Multiply(phi1,phi2),phi3)
    } else{phi4 <- modelList$phi
    phi <- Multiply(Multiply(Multiply(phi1,phi2),phi3),phi4)
    }
    testPhi <- modelList$phi 
    
    
    #initialization
    U <- rep(0,monteCarlo)
    UNu <- rep(0,monteCarlo)
    sigmaHat <- rep(0,monteCarlo)
    gammaHat <- rep(0,monteCarlo)
    lambdaHat <- rep(0,monteCarlo)
    x <- matrix(0,monteCarlo,N)
    nuHat <- matrix(0,monteCarlo,N)
    nuError2 <- rep(0,monteCarlo)
    stdNoise <- matrix(0,monteCarlo,sampleSize) #Keep track of std of noise for confidence intervals
    stdNoiseReal <- matrix(0,monteCarlo,sampleSize)
    sk1 <- matrix(0,monteCarlo,sampleSize)
    snop1 <- matrix(0,monteCarlo,sampleSize)
    varSigma2Hat <- rep(0,monteCarlo)
    varGammaHat <- rep(0,monteCarlo)
    varLambdaHat <- rep(0,monteCarlo)
    varNuHat <- matrix(0,monteCarlo,N)
  
    #Simulation of Option Prices
    for(mi in 1:monteCarlo){
      optionList <- simulateOptionPrices(r,time[k], model, linear, sampleSize, noiseLevel, design)
      snop1[mi,] <- optionList$snop
      sk1[mi,] <- optionList$sk
      if(noise == "real"){
        stdNoise[mi,] <- optionList$stdNoise}
      ##Assign values in every Monte-Carlo iteration from calibration
      if(k ==1){
        calibrationList <- calibrationHom(sk1[mi,],snop1[mi,],mode, time[k])
        if(undersmoothing ==TRUE){
          Ufix <- zeta*calibrationList$U
          #UfixNu <- zeta*calibrationList$UNu
          UfixNu <- calibrationList$UNu
          calibrationList <- calibrationHom(sk1[mi,],snop1[mi,],"fix", time[k])
        }
        if(noise == "estimated"){
          stdNoiseReal[mi,] <- optionList$stdNoise
          stdNoise[mi,] <- estNoiseHom(sk1[mi,],snop1[mi,],time[k])
        }
      } else{
        calibrationList <- calibrationInHom(sk1[mi,], snop1[mi,],time[k], sk0[mi,], snop0[mi,], time[k-1], mode)
        v0 <- calibrationList$v0
        phiNum0 <- calibrationList$phiNum0
        phiNu0 <- calibrationList$phiNu0
        if(undersmoothing ==TRUE){
          Ufix <- zeta*calibrationList$U
          #UfixNu <- zeta*calibrationList$UNu
          UfixNu <- calibrationList$UNu
          calibrationList <- calibrationInHom(sk1[mi,], snop1[mi,],time[k], sk0[mi,], snop0[mi,], time[k-1], "fix")
        }
        if(noise == "estimated"){
          stdNoiseReal[mi,] <- optionList$stdNoise
          stdNoise[mi,] <- estNoiseInHom(sk1[mi,], snop1[mi,],time[k], sk0[mi,], snop0[mi,], time[k-1])
        }
      }
      U[mi] <- calibrationList$U
      UNu[mi] <- calibrationList$UNu
      sigmaHat[mi] <- calibrationList$sigmaHat
      gammaHat[mi] <- calibrationList$gammaHat
      lambdaHat[mi] <- calibrationList$lambdaHat
      x[mi,] <- calibrationList$x
      nuHat[mi,] <- calibrationList$nuHat
      nuError2[mi] <- calibrationList$nuError2
      v <- calibrationList$v
      phiNum <- calibrationList$phiNum
      phiNu <- calibrationList$phiNum
      
      #Variance of parameters:
      confInt <- confIntervals(k,mi)
      varSigma2Hat[mi] <- confInt$varSigma2
      varGammaHat[mi] <- confInt$varGamma
      varLambdaHat[mi] <- confInt$varLambda
      varNuHat[mi,] <- confInt$varNu
    }
    assign(paste("U",k,sep=""),U)
    assign(paste("UNu",k,sep=""),UNu)
    assign(paste("sigmaHat",k,sep=""),sigmaHat)
    assign(paste("gammaHat",k,sep=""),gammaHat)
    assign(paste("lambdaHat",k,sep=""),lambdaHat)
    assign(paste("x",k,sep=""), x)
    assign(paste("nuHat",k,sep=""), nuHat)
    assign(paste("nuError2",k,sep=""), nuError2)
    assign(paste("nuTrue",k,sep=""), calibrationList$nuTrue)
    
    #Conf Intervals
    assign(paste("stdNoise",k,sep=""), stdNoise)
    if(noise == "estimated"){
      assign(paste("stdNoiseReal",k,sep=""), stdNoiseReal)
    }
    assign(paste("varSigma2Hat",k,sep=""), varSigma2Hat)
    assign(paste("varGammaHat",k,sep=""), varGammaHat)
    assign(paste("varLambdaHat",k,sep=""), varLambdaHat)
    assign(paste("varNuHat",k,sep=""), varNuHat)
    assign(paste("covProb",k,sep=""), coverageProbability(k,mu,mu+delta))
    
    
    sk0 <- sk1
    snop0 <- snop1
  }
    #Delete many unneeded variables
    rm(mu, delta,U,UNu,sigmaHat,gammaHat,lambdaHat,x,nuHat,nuError2,stdNoise,varSigma2Hat,varGammaHat,varLambdaHat,varNuHat,sigma,gamma,lambda,T1,T2,T3,T4,phiNu,phiNum,v,sk0,sk1,snop0,snop1)
    #Calculation of RMSE for estimators with simulations
    rRMSE <- function(actual,predicted){sqrt(mean((actual-predicted)^2)/mean(predicted^2))}
    sigma1RMSE <- rRMSE(sigma1,sigmaHat1)
    sigma2RMSE <- rRMSE(sigma2,sigmaHat2)
    sigma3RMSE <- rRMSE(sigma3,sigmaHat3)
    #sigma4RMSE <- RMSE(sigma4,sigmaHat4)
    gamma1RMSE <- rRMSE(gamma1,gammaHat1)
    gamma2RMSE <- rRMSE(gamma2,gammaHat2)
    gamma3RMSE <- rRMSE(gamma3,gammaHat3)
    #gamma4RMSE <- RMSE(gamma4,gammaHat4)
    lambda1RMSE <- rRMSE(lambda1,lambdaHat1)
    lambda2RMSE <- rRMSE(lambda2,lambdaHat2)
    lambda3RMSE <- rRMSE(lambda3,lambdaHat3)
    #lambda4RMSE <- RMSE(lambda4,lambdaHat4)
    
    RMSENu <- function(x, nuHat, nu){
      w <- rep(1,length(x))
      w[1] <- 0.5
      w[length(x)]<- 0.5
      norm <- 0
      mi <- dim(nuHat)[1]
      for(i in 1:mi){
        norm <- norm + (x[i,2]-x[i,1])*sum((nu(x1[i,])-nuHat[i,])^2*w)
        print((x[i,2]-x[i,1])*sum((nu(x1[i,])-nuHat[i,])^2*w))
      }
      sqrt(1/mi*norm)
    }
    
    #Plots after Monte Carlo Simulations
    df <- data.frame(sigmaHat1,sigmaHat2,sigmaHat3,gammaHat1,gammaHat2,gammaHat3,lambdaHat1, lambdaHat2,lambdaHat3,varSigma2Hat1,varSigma2Hat2,varSigma2Hat3,varGammaHat1,varGammaHat2,varGammaHat3,varLambdaHat1,varLambdaHat2,varLambdaHat3)
    ggSigma1 <- ggplot(df, aes(x = seq(1,length(sigmaHat1)),y= sigmaHat1)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("volatility $\\sigma_1$")) +xlab('') + 
      geom_hline(yintercept = sigma1, color = 'black', size = 1.0) 
      # geom_hline(yintercept = SigmaHat1-lowerQ*sqrt(varSigmaHat1), color = 'black', linetype = "dashed",size = 1.0) +
      # geom_hline(yintercept = varSigmaHat1+upperQ*sqrt(varSigmaHat1), color = 'black', linetype = "dashed",size = 1.0)
    ggSigma1 <- ggMarginal(ggSigma1, type="boxplot", fill="gray", col = "black", margins="y")
    ggSigma2 <- ggplot(df, aes(x = seq(1,length(sigmaHat2)),y= sigmaHat2)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("volatility $\\sigma_2$")) +xlab('') + 
      geom_hline(yintercept = sigma2, color = 'black', size = 1.0)
      # geom_hline(yintercept = sigmaHat2-lowerQ*sqrt(varSigmaHat2), color = 'black', linetype = "dashed",size = 1.0) +
      # geom_hline(yintercept = sigmaHat2+upperQ*sqrt(varSigmaHat2), color = 'black', linetype = "dashed",size = 1.0)
    ggSigma2 <- ggMarginal(ggSigma2, type="boxplot", fill="gray", col = "black", margins="y")
    ggSigma3 <- ggplot(df, aes(x = seq(1,length(sigmaHat3)),y= sigmaHat3)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("volatility $\\sigma_3$")) +xlab('') +
      geom_hline(yintercept = sigma3, color = 'black', size = 1.0)
      # geom_hline(yintercept = sigmaHat3-lowerQ*sqrt(varSigmaHat3), color = 'black', linetype = "dashed",size = 1.0) +
      # geom_hline(yintercept = sigmaHat3+upperQ*sqrt(varSigmaHat3), color = 'black', linetype = "dashed",size = 1.0)
    ggSigma3 <- ggMarginal(ggSigma3, type="boxplot", fill="gray", col = "black", margins="y")
    # #ggSigma4 <- ggplot(df, aes(x = seq(1,length(sigmaHat4)),y= sigmaHat4)) + geom_point(alpha=1/3,colour = "black") +
    #   theme_bw()+ylab(TeX("volatility $\\sigma_4$")) +xlab('') +
    #   geom_hline(yintercept = sigma4, color = 'black', size = 1.0)
    #   # geom_hline(yintercept = sigmaHat4-lowerQ*sqrt(varSigmaHat4), color = 'black', linetype = "dashed",size = 1.0) +
    #   # geom_hline(yintercept = sigmaHat4+upperQ*sqrt(varSigmaHat4), color = 'black', linetype = "dashed",size = 1.0)
    # ggSigma4 <- ggMarginal(ggSigma4, type="boxplot", fill="gray", col = "black", margins="y")

    
    ggGamma1 <- ggplot(df, aes(x = seq(1,length(gammaHat1)),y= gammaHat1)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("drift $\\gamma_1$")) +xlab('') +
      geom_hline(yintercept = gamma1, color = 'black', size = 1.0)
    ggGamma1 <- ggMarginal(ggGamma1, type="boxplot", fill="gray", col = "black", margins="y")
    ggGamma2 <- ggplot(df, aes(x = seq(1,length(gammaHat2)),y= gammaHat2)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("drift $\\gamma_2$")) +xlab('') +
      geom_hline(yintercept = gamma2, color = 'black', size = 1.0)
    ggGamma2 <- ggMarginal(ggGamma2, type="boxplot", fill="gray", col = "black", margins="y")
    ggGamma3 <- ggplot(df, aes(x = seq(1,length(gammaHat3)),y= gammaHat3)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("drift $\\gamma_3$")) +xlab('') +
      geom_hline(yintercept = gamma3, color = 'black', size = 1.0)
    ggGamma3 <- ggMarginal(ggGamma3, type="boxplot", fill="gray", col = "black", margins="y")
    # ggGamma4 <- ggplot(df, aes(x = seq(1,length(gammaHat4)),y= gammaHat4)) + geom_point(alpha=1/3,colour = "black") +
    #   theme_bw()+ylab(TeX("drift $\\gamma_4$")) +xlab('') +
    #   geom_hline(yintercept = gamma4, color = 'black', size = 1.0)
    # ggGamma4 <- ggMarginal(ggGamma4, type="boxplot", fill="gray", col = "black", margins="y")
    
    ggLambda1 <- ggplot(df, aes(x = seq(1,length(lambdaHat1)),y= lambdaHat1)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("intensity $\\lambda_1$")) +xlab('Monte Carlo Simulations') +
      geom_hline(yintercept = lambda1, color = 'black', size=1.0)
    ggLambda1 <- ggMarginal(ggLambda1, type="boxplot", fill="gray", col = "black", margins="y")
    ggLambda2 <- ggplot(df, aes(x = seq(1,length(lambdaHat2)),y= lambdaHat2)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("intensity $\\lambda_2$")) +xlab('Monte Carlo Simulations') +
      geom_hline(yintercept = lambda2, color = 'black', size=1.0)
    ggLambda2 <- ggMarginal(ggLambda2, type="boxplot", fill="gray", col = "black", margins="y")
    ggLambda3 <- ggplot(df, aes(x = seq(1,length(lambdaHat3)),y= lambdaHat3)) + geom_point(alpha=1/3,colour = "black") +
      theme_bw()+ylab(TeX("intensity $\\lambda_3$")) +xlab('Monte Carlo Simulations') +
      geom_hline(yintercept = lambda3, color = 'black', size=1.0)
    ggLambda3 <- ggMarginal(ggLambda3, type="boxplot", fill="gray", col = "black", margins="y")
    # ggLambda4 <- ggplot(df, aes(x = seq(1,length(lambdaHat4)),y= lambdaHat4)) + geom_point(alpha=1/3,colour = "black") +
    #   theme_bw()+ylab(TeX("intensity $\\lambda_4$")) +xlab('Monte Carlo Simulations') +
    #   geom_hline(yintercept = lambda4, color = 'black', size=1.0)
    # ggLambda4 <- ggMarginal(ggLambda4, type="boxplot", fill="gray", col = "black", margins="y")
    
    grid.arrange(ggSigma1,ggSigma2, ggSigma3, ggGamma1, ggGamma2,ggGamma3,ggLambda1,ggLambda2,ggLambda3, ncol=3)
    
    
    
    ggEstimator2 <- arrangeGrob(ggSigma2,ggGamma2,ggLambda2)
    #ggsave("estimator2.png", plot=ggEstimator2, dpi=500)
  
    ggU3 <- ggplot(data.frame(U3), aes(x=seq(1,length(U3)),y=U3))+geom_point(alpha=1/2,colour="black") + theme_bw()+ylab(TeX("cut-off $\\U_3$")) +xlab('Monte Carlo Simulations')
    ggU3
    #ggsave("cutoff1.png", plot=ggU1,dpi=500)
    ###ErrorBar Plot 
    ggSigma1 <- ggplot(data.frame(sigmaHat1,varSigma2Hat1), aes(x = seq(1,length(sigmaHat1)),y= sigmaHat1^2)) +
      theme_bw()+ylab(TeX("volatility $\\sigma_1^2$")) +xlab('')+ 
      geom_hline(yintercept = sigma1^2, color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=sigmaHat1^2+lowerQ*sqrt(varSigma2Hat1), ymax=sigmaHat1^2+upperQ*sqrt(varSigma2Hat1)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggSigma2 <- ggplot(data.frame(sigmaHat2,varSigma2Hat2), aes(x = seq(1,length(sigmaHat2)),y= sigmaHat2^2)) +
      theme_bw()+ylab(TeX("volatility $\\sigma_2^2$")) +xlab('')+ 
      geom_hline(yintercept = sigma2^2, color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=sigmaHat2^2+lowerQ*sqrt(varSigma2Hat2), ymax=sigmaHat2^2+upperQ*sqrt(varSigma2Hat2)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggSigma3 <- ggplot(data.frame(sigmaHat3,varSigma2Hat3), aes(x = seq(1,length(sigmaHat3)),y= sigmaHat3^2)) +
      theme_bw()+ylab(TeX("volatility $\\sigma_3^2$")) +xlab('')+ 
      geom_hline(yintercept = sigma3^2, color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=sigmaHat3^2+lowerQ*sqrt(varSigma2Hat3), ymax=sigmaHat3^2+upperQ*sqrt(varSigma2Hat3)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    # ggSigma4 <- ggplot(df, aes(x = seq(1,length(sigmaHat4)),y= sigmaHat4^2)) +
    #   theme_bw()+ylab(TeX("volatility $\\sigma_4^2$")) +xlab('')+ 
    #   geom_hline(yintercept = sigma4^2, color = 'black', size = 1.0,alpha =0.5) +
    #   geom_errorbar(aes(ymin=sigmaHat4^2+lowerQ*sqrt(varSigma2Hat4), ymax=sigmaHat4^2+upperQ*sqrt(varSigma2Hat4)),width=0.1, colour="aquamarine4", size=0.8)+
    #   geom_point(colour = "black", size=2)
    ggGamma1 <- ggplot(data.frame(gammaHat1,varGammaHat1), aes(x = seq(1,length(gammaHat1)),y= gammaHat1)) +
      theme_bw()+ylab(TeX("drift $\\gamma_1$")) +xlab('') +
      geom_hline(yintercept = gamma1, color = 'black', size = 1.0, alpha = 0.5)+
      geom_errorbar(aes(ymin=gammaHat1+lowerQ*sqrt(varGammaHat1), ymax=gammaHat1+upperQ*sqrt(varGammaHat1)),width=0.1, colour="aquamarine4", size=0.8)+ 
      geom_point(colour = "black",size=2)
    ggGamma2 <- ggplot(data.frame(gammaHat2,varGammaHat2), aes(x = seq(1,length(gammaHat2)),y= gammaHat2)) +
      theme_bw()+ylab(TeX("drift $\\gamma_2$")) +xlab('') +
      geom_hline(yintercept = gamma2, color = 'black', size = 1.0, alpha = 0.5)+
      geom_errorbar(aes(ymin=gammaHat2+lowerQ*sqrt(varGammaHat2), ymax=gammaHat2+upperQ*sqrt(varGammaHat2)),width=0.1, colour="aquamarine4", size=0.8)+ 
      geom_point(colour = "black",size=2)
    ggGamma3 <- ggplot(data.frame(gammaHat3,varGammaHat3), aes(x = seq(1,length(gammaHat3)),y= gammaHat3)) +
      theme_bw()+ylab(TeX("drift $\\gamma_3$")) +xlab('') +
      geom_hline(yintercept = gamma3, color = 'black', size = 1.0, alpha = 0.5)+
      geom_errorbar(aes(ymin=gammaHat3+lowerQ*sqrt(varGammaHat3), ymax=gammaHat3+upperQ*sqrt(varGammaHat3)),width=0.1, colour="aquamarine4", size=0.8)+ 
      geom_point(colour = "black",size=2)
    # ggGamma4 <- ggplot(df, aes(x = seq(1,length(gammaHat1)),y= gammaHat4)) +
    #   theme_bw()+ylab(TeX("drift $\\gamma_4$")) +xlab('') +
    #   geom_hline(yintercept = gamma4, color = 'black', size = 1.0, alpha = 0.5)+
    #   geom_errorbar(aes(ymin=gammaHat4+lowerQ*sqrt(varGammaHat4), ymax=gammaHat4+upperQ*sqrt(varGammaHat4)),width=0.1, colour="aquamarine4", size=0.8)+ 
    #   geom_point(colour = "black",size=2)
    ggLambda1 <- ggplot(data.frame(lambdaHat1,varLambdaHat1), aes(x = seq(1,length(lambdaHat1)),y= lambdaHat1)) +
      theme_bw()+ylab(TeX("drift $\\lambda_1$")) +xlab('') +
      geom_hline(yintercept = lambda1, color = 'black', size = 1.0, alpha = 0.5)+
      geom_errorbar(aes(ymin=lambdaHat1+lowerQ*sqrt(varLambdaHat1), ymax=lambdaHat1+upperQ*sqrt(varLambdaHat1)),width=0.1, colour="aquamarine4", size=0.8)+ 
      geom_point(colour = "black",size=2)
    ggLambda2 <- ggplot(data.frame(lambdaHat2,varLambdaHat2), aes(x = seq(1,length(lambdaHat2)),y= lambdaHat2)) +
      theme_bw()+ylab(TeX("drift $\\lambda_2$")) +xlab('') +
      geom_hline(yintercept = lambda2, color = 'black', size = 1.0, alpha = 0.5)+
      geom_errorbar(aes(ymin=lambdaHat2+lowerQ*sqrt(varLambdaHat2), ymax=lambdaHat2+upperQ*sqrt(varLambdaHat2)),width=0.1, colour="aquamarine4", size=0.8)+ 
      geom_point(colour = "black",size=2)
    ggLambda3 <- ggplot(data.frame(lambdaHat3,varLambdaHat3), aes(x = seq(1,length(lambdaHat3)),y= lambdaHat3)) +
      theme_bw()+ylab(TeX("drift $\\lambda_3$")) +xlab('') +
      geom_hline(yintercept = lambda3, color = 'black', size = 1.0, alpha = 0.5)+
      geom_errorbar(aes(ymin=lambdaHat3+lowerQ*sqrt(varLambdaHat3), ymax=lambdaHat3+upperQ*sqrt(varLambdaHat3)),width=0.1, colour="aquamarine4", size=0.8)+ 
      geom_point(colour = "black",size=2)
    # ggLambda4 <- ggplot(df, aes(x = seq(1,length(lambdaHat4)),y= lambdaHat4)) +
    #   theme_bw()+ylab(TeX("drift $\\lambda_4$")) +xlab('') +
    #   geom_hline(yintercept = lambda4, color = 'black', size = 1.0, alpha = 0.5)+
    #   geom_errorbar(aes(ymin=lambdaHat4+lowerQ*sqrt(varLambdaHat4), ymax=lambdaHat4+upperQ*sqrt(varLambdaHat4)),width=0.1, colour="aquamarine4", size=0.8)+ 
    #   geom_point(colour = "black",size=2)
    
    ggNuMu1 <- ggplot(data.frame(nuHat0 = nuHat1[,which.min(x1[1,]< -0.1)],varNuHat0 = varNuHat1[,which.min(x1[1,]< -0.1)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("density $\\nu_1 (\\mu_1)$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu1(-0.1), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggNuDelta1 <- ggplot(data.frame(nuHat0 = nuHat1[,which.min(x1[1,]< 0.1)],varNuHat0 = varNuHat1[,which.min(x1[1,]< 0.1)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("density $\\nu_1 (\\mu_1+\\delta_1)$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu1(0.1), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggNuMu2 <- ggplot(data.frame(nuHat0 = nuHat2[,which.min(x2[1,]< -0.2)],varNuHat0 = varNuHat2[,which.min(x2[1,]< -0.2)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("density $\\nu_2 (\\mu_2)$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu2(-0.2), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggNuDelta2 <- ggplot(data.frame(nuHat0 = nuHat2[,which.min(x2[1,]< 0.2)],varNuHat0 = varNuHat2[,which.min(x2[1,]< 0.2)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("density $\\nu_2 (\\mu_2+\\delta_2)$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu2(0.2), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggNuMu3 <- ggplot(data.frame(nuHat0 = nuHat3[,which.min(x3[1,]< -0.1)],varNuHat0 = varNuHat3[,which.min(x3[1,]< -0.1)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("density $ \\nu_3 (\\mu_3)$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu3(-0.1), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggNuDelta3 <- ggplot(data.frame(nuHat0 = nuHat3[,which.min(x3[1,]< 0.2)],varNuHat0 = varNuHat3[,which.min(x3[1,]< 0.2)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("density $\\nu_3 (\\mu_3+\\delta_3)$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu3(0.2), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    
    grid.arrange(ggSigma1,ggSigma2, ggSigma3, ggGamma1, ggGamma2,ggGamma3,ggLambda1,ggLambda2,ggLambda3, ggNuMu1, ggNuMu2, ggNuMu3, ggNuDelta1, ggNuDelta2, ggNuDelta3, ncol=3)
    
    ggConf1 <- arrangeGrob(ggSigma1,ggGamma1,ggLambda1, ggNuMu1, ggNuDelta1)
    plot(ggConf1)
    #ggsave("conf1Zeta.png", plot=ggConf1, dpi=500, width=4, height =4)
    ggConf2 <- arrangeGrob(ggSigma2,ggGamma2,ggLambda2, ggNuMu2, ggNuDelta2)
    plot(ggConf2)
    #ggsave("conf2Zeta.png", plot=ggConf2, dpi=500, width=4, height =4)
    ggConf3 <- arrangeGrob(ggSigma3,ggGamma3,ggLambda3, ggNuMu3, ggNuDelta3)
    plot(ggConf3)
    #ggsave("conf3Zeta.png", plot=ggConf3, dpi=500, width=4, height =4)
    
    # ggU3 <- ggplot(data.frame(U3), aes(x=seq(1,length(U3)),y=U3))+geom_point(alpha=1/2,colour="black") + theme_bw()+ylab(TeX("cut-off $\\U_3$")) +xlab('Monte Carlo Simulations')
    # ggU3
    # ggsave("cutoff3.png", plot=ggU3,dpi=500)
    ###ErrorBar Plot 
    
    ###Levy density plots

    #plot(x2[,],nuHat2[,], type = "l", col = 'grey50', xlab= "grid", ylab=TeX("Levy density $\\nu"))
    #lines(x2[2,],Nu2(x2[2,]), col = "black")
    #jpeg("nu1.jpeg", width = 4, height = 4, units = 'in', res = 900)
    par(mfrow=c(1,1))
    matplot(t(x1[,]),t(nuHat1[,]), type = "l", col = 'grey60', xlab= TeX("$x$"), ylab=TeX("$\\nu_1(x)$"),xlim=c(-1.5,1.5)) 
    grid()
    lines(x1[1,],Nu1(x1[1,]), col = "black", lwd=2)
    #dev.off()
    
    
    # ggMatPlot1<-ggmatplot(t(x1[1:10,]),t(nuHat1[1:10,]), plot_type = "line", linetype = "solid", color="grey60", size=1, xlab= TeX("$x$"), ylab=TeX("$\\nu_1(x)$"), show.legend = FALSE)+theme_bw()+geom_line(data=data.frame(x=x1[1,],y=Nu1(x1[1,])), aes(x=x,y=y), inherit.aes = FALSE, color="black", size=1)
    # ggMatPlot1
    # ggsave("nu1.jpg",ggMatPlot1, width = 6, height=6, dpi=500)
    
    
    #jpeg("nu2.jpeg", width = 4, height = 4, units = 'in', res = 900)
    matplot(t(x2[,]),t(nuHat2[,]), type = "l", col = 'grey60', xlab= "x", ylab=TeX("$\\nu_2(x)$"),xlim=c(-1.5,1.5)) 
    grid()
    lines(x2[1,],Nu2(x2[1,]), col = "black",lwd=2)
    #dev.off()
    
    #jpeg("nu3.jpeg", width = 4, height = 4, units = 'in', res = 900)
    matplot(t(x3[,]),t(nuHat3[,]), type = "l", col = 'grey60', xlab= "x", ylab=TeX("$\\nu_3(x)$"),xlim=c(-1.5,1.5))  
    grid()
    lines(x3[1,],Nu3(x3[1,]), col = "black", lwd=2)
    #dev.off()
  
    #matplot(t(x4[,]),t(nuHat4[,]), type = "l", col = 'grey50', xlab= "grid", ylab=TeX("Levy density $\\nu$")) 
    #lines(x4[1,],Nu4(x4[1,]), col = "black")
    # grid()
    
    #jpeg("nu1.jpeg", width = 4, height = 4, units = 'in', res = 600)
    par(mfrow=c(1,1))
    matplot(t(x1[,]),t(nuHat1[,]), type = "l", col = 'grey60', xlab= TeX("$x$"), ylab=TeX("Levy density $\\nu_1(x)$")) 
    grid()
    lines(x1[1,],Nu1(x1[1,]), col = "black", lwd=2)
    #dev.off()
    

    n <- 98
    par(mfrow=c(1,1))
    #jpeg("nuConf1.jpeg", width = 4, height = 4, units = 'in', res = 900)
    plot(x1[n,],nuHat1[n,]+upperQ*sqrt(varNuHat1[n,]), type="l", lwd=2, col = "black", xlab= TeX("x"), ylab=TeX("Levy density $\\nu_1(x)"))
    grid()
    lines(x1[n,],nuHat1[n,], type = "l",lwd=2, col = 'gray50') 
    lines(x1[n,],pmax(0,nuHat1[n,]+lowerQ*sqrt(varNuHat1[n,])),lwd=2, col="black")
    lines(x1[n,],Nu1(x1[n,]), type="l",lwd=2, col="red")
    #dev.off()
    #jpeg("nuConf2.jpeg", width = 4, height = 4, units = 'in', res = 900)
    plot(x2[n,],nuHat2[n,]+upperQ*sqrt(varNuHat2[n,]), type="l", lwd=2, col = "black", xlab= TeX("x"), ylab=TeX("Levy density $\\nu_2(x)"))
    grid()
    lines(x2[n,],nuHat2[n,], type = "l",lwd=2, col = 'grey50') 
    lines(x2[n,],pmax(0,nuHat2[n,]+lowerQ*sqrt(varNuHat2[n,])),lwd=2, col="black")
    lines(x2[n,],Nu2(x2[n,]), type="l",lwd=2, col="red")
    #dev.off()
    #jpeg("nuConf3.jpeg", width = 4, height = 4, units = 'in', res = 900)
    plot(x3[n,],nuHat3[n,]+upperQ*sqrt(varNuHat3[n,]), type="l", lwd=2, col = "black", xlab= TeX("x"), ylab=TeX("Levy density $\\nu_3(x)"))
    grid()
    lines(x3[n,],nuHat3[n,], type = "l",lwd=2, col = 'grey50') 
    lines(x3[n,],pmax(0,nuHat3[n,]+lowerQ*sqrt(varNuHat3[n,])),lwd=2, col="black")
    lines(x3[n,],Nu3(x3[n,]), type="l",lwd=2, col="red")
    #dev.off()
    # plot(x4[n,],nuHat4[n,]+upperQ*sqrt(varNuHat4[n,]), type="l", col = "black", xlab= "grid", ylab=TeX("Levy density $\\nu"))
    # lines(x4[n,],nuHat4[n,], type = "l", col = 'grey50', xlab= "grid", ylab=TeX("Levy density $\\nu")) 
    # lines(x4[n,],pmax(0,nuHat4[n,]+lowerQ*sqrt(varNuHat4[n,])), col="black")
    # lines(x4[n,],Nu4(x4[n,]), type="l", col="red")
    
    ###boxPlot of values of Levy density at x (same as Sohl paper)
    
    x <- -0.1
    ggNux_2 <- ggplot(data.frame(nuHat0 = nuHat3[,which.min(x3[2,]< x)],varNuHat0 = varNuHat3[,which.min(x3[2,]< x)]), aes(x = seq(1,length(nuHat0)),y= nuHat0)) +
      theme_bw()+ylab(TeX(paste("$\\nu (",x,")$"))) +xlab('Monte Carlo Simulations')+ 
      geom_hline(yintercept = Nu3(x), color = 'black', size = 1.0,alpha =0.5) +
      geom_errorbar(aes(ymin=nuHat0+lowerQ*sqrt(varNuHat0), ymax=nuHat0+upperQ*sqrt(varNuHat0)),width=0.1, colour="aquamarine4", size=0.8)+
      geom_point(colour = "black", size=2)
    ggNux_2 <- ggMarginal(ggNux_2, type="boxplot", fill="gray", col = "black", margins="y")
    ggNux_2
    
    #Plot to show misbehaviour of Nu3
    ggBoxPlot1 <-  ggplot(data.frame(nuHat0 = nuHat3[,which.min(x3[1,]<x)],varNuHat0 = varNuHat3[,which.min(x3[1,]<x)]), aes(x = seq(1,length(nuHat0)),y= nuHat0))+geom_boxplot(fill="gray", col = "black")+theme_bw()+ ylab(TeX("$\\tilde{\\nu}_3(\\mu_3)$"))+xlab("")
    ggBoxPlot2 <-  ggplot(data.frame(nuHat0 = nuHat3[,which.min(x3[1,]<x)],varNuHat0 = varNuHat3[,which.min(x3[1,]<x)]), aes(x = seq(1,length(nuHat0)),y= nuHat0))+geom_boxplot(fill="gray", col = "black")+theme_bw()+ ylab(TeX("$\\tilde{\\nu}_3(\\mu_3)$"))+xlab("")+ylim(3.0,5.2)
    ggBoxPlot<- arrangeGrob(ggBoxPlot1,ggBoxPlot2, nrow=1)
    plot(ggBoxPlot)
    #ggsave("boxplot.jpg",ggBoxPlot, width = 4, height=3)
    
    #Some Examples of the misbehaviour
    z <-98
    ggplot(data.frame(x=x3[z,], y=nuHat3[z,]), aes(x=x,y=y))+geom_line() +theme_bw()+ ylab(TeX("Simulation 77 of $\\tilde{\\nu}_3(x)$"))+xlab(TeX("$$x$$"))
    #56, 16, 52
    ggNu3Mis1 <- ggplot(data.frame(x=x3[85,], y=nuHat3[85,]), aes(x=x,y=y))+geom_line() +theme_bw()+ ylab(TeX("Simulation 85 of $\\tilde{\\nu}_3(x)$"))+xlab(TeX("$$x$$"))
    ggNu3Mis2 <- ggplot(data.frame(x=x3[16,], y=nuHat3[16,]), aes(x=x,y=y))+geom_line() +theme_bw()+ ylab(TeX("Simulation 16 of $\\tilde{\\nu}_3(x)$"))+xlab(TeX("$$x$$"))
    ggNu3Mis3 <- ggplot(data.frame(x=x3[10,], y=nuHat3[10,]), aes(x=x,y=y))+geom_line() +theme_bw()+ ylab(TeX("Simulation 10 of $\\tilde{\\nu}_3(x)$"))+xlab(TeX("$$x$$"))
    ggNu3Mis <- arrangeGrob(ggNu3Mis1, ggNu3Mis2,ggNu3Mis3, nrow=1)
    plot(ggNu3Mis)
    #ggsave("Nu3MisBehavior.jpg",ggNu3Mis, width = 4, height=2)
    
    #Investigation of Estimated Noise
    plot(stdNoiseReal1[1,])
    plot(stdNoise1[1,])
    plot(stdNoiseReal1[1,]-stdNoise1[1,])
    
    #Generated Paths
    NoOfSteps <- 1e4
    paths <- generatePaths(r,NoOfSteps)
    t <- paths$t
    St <- paths$St
    ggPath <- ggplot(data = data.frame(t,St), aes(x=t,y=St)) +geom_line()+theme_bw()+ylab(TeX("stock price $S_t$")) +xlab(TeX('time $t$'))+
      geom_vline(xintercept = c(0,time),linetype="solid", color = "grey30", size=0.5) 
    ggPath
    #ggsave("samplePath.png", plot=ggPath, dpi=500, width=6, height =4)
}

dataType = "historical"
### REAL DATA, choose ticker in function
if(model == "real_data"){
  #Calibration parameters
  #r <- 0.03
  s <- 2          #Smoothness of jump density, is s=2 a realistic assumption?
  mode <- "PLS"    #oracle only usable for simulations, other choices: "flat","PLS","fix"
  noise <- "estimated"
  # Ufix <- 23
  # UfixNu <-8
  
  if(dataType == "live"){
    # First Choose tickr and expiration date, the rest of the script loads all the data.
    tickr <- "SPY" #SPY = ETF S&P 500, AAPL = Apple, QQQ= Nasdaq, GLD = Gold ETF 
    expiration <-c("2022-10-03")
    S <- getQuote(tickr, src = "yahoo")[,2] #Current price of stock
    #Calculate Maturity 
    quoteDate = Sys.Date()
    time <- as.numeric(difftime(expiration, quoteDate, units = "weeks"))/52
  }
  if (dataType == 'historical'){
    #Historical data for SPX
    #histData <- read.csv("/Users/loekkoorevaar/Downloads/Batch_PRO_Sample_PHC/SPX_20180904_to_20180928.csv", header = TRUE)
    #histData <- histData[histData$QuoteDate == '09/28/2018',]
    #S <- histData$UnderlyingPrice[1]
    #quoteDate <- as.POSIXct(histData$QuoteDate[1],  format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    #expiration1 <- unique(histData$Expiration)[1:2]
    #expiration <- as.POSIXct(expiration1, format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    #time <- as.numeric(difftime(expiration,quoteDate,units="weeks"))/52
    
    #histData <- read.csv("/Users/loekkoorevaar/Downloads/Sample_L2_20190815/options_20190815.csv", header = TRUE)
    #histData <- histData[histData$DataDate == "08/15/2019 16:00",]
    #histData <- histData[histData$UnderlyingSymbol == "A",]
    #S  <- histData$UnderlyingPrice[1]
    #quoteDate <- as.POSIXct(histData$DataDate[1],  format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    #expiration1 <- unique(histData$Expiration)[1:2]
    #expiration <- as.POSIXct(expiration1, format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    #time <- as.numeric(difftime(expiration,quoteDate,units="weeks"))/52
    
    #histData <- read.csv("/Users/loekkoorevaar/Downloads/bb_options_20190815.csv", header = TRUE)
    #histData <- histData[histData$DataDate == "08/15/2019",]
    #histData <- histData[histData$UnderlyingSymbol == "SPY",]
    #S  <- histData$UnderlyingPrice[1]
    #quoteDate <- as.POSIXct(histData$DataDate[1],  format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    #expiration1 <- unique(histData$Expiration)[1:2]
    #expiration <- as.POSIXct(expiration1, format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    #time <- as.numeric(difftime(expiration,quoteDate,units="weeks"))/52
    
    histData <- read.csv("/Users/loekkoorevaar/Downloads/Batch_3REYDXXRXE/SPY_2022.csv", header = TRUE)
    histData <- histData[histData$quotedate == "08/02/2022",]
    histData <- histData[histData$volume>=5,]

    S  <- histData$underlying_last[1]
    quoteDate <- as.POSIXct(histData$quotedate[1],  format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    expiration1 <- c() #mae 
    for (exp in unique(histData$expiration)){
      if (nrow(histData[histData$expiration == exp,])>=50){
        expiration1 <- append(expiration1,exp)
      }
    }
    #expiration1 <- expiration1[c(3,6,9)]
    expiration <- as.POSIXct(expiration1, format = '%m/%d/%Y', tz = "Europe/Amsterdam")
    time <- as.numeric(difftime(expiration,quoteDate,units="weeks"))/52
  }
  
  #Initilialization
  
  N <- 4096
  U <- rep(0,length(expiration))
  r_list <- rep(0,length(expiration))
  UNu <- rep(0,length(expiration))
  sigmaHat <- rep(0,length(expiration))
  gammaHat <- rep(0,length(expiration))
  lambdaHat <- rep(0,length(expiration))
  x <- matrix(0,length(expiration),N)
  nuHat <- matrix(0,nrow=length(expiration),ncol=N)
  varSigma2Hat <-  rep(0,length(expiration))
  varGammaHat <-  rep(0,length(expiration))
  varLambdaHat <-  rep(0,length(expiration))
  varNuHat <-  matrix(0,nrow=length(expiration),ncol=N)
  
  for(k in 1:length(expiration)){
    if(k == 1){T0 <- 0} else{T0 <- time[k-1]}
    T1 <- time[k]
    if (dataType == "live"){
      optionData <- getOptionChain(tickr, Exp = expiration[k], src="yahoo") #option data, try orats
      optionData$calls = optionData$calls[optionData$calls[["LastTradeTime"]]>"2022-09-23 09:30:00 EDT",]
      #optionData$calls = optionData$calls[optionData$calls[["Vol"]]>5,]
      optionData$puts = optionData$puts[optionData$puts[["LastTradeTime"]]>"2022-09-23 09:30:00 EDT",]
      #optionData$puts = optionData$puts[optionData$puts[["Vol"]]>5,]
      kCalls <- optionData$calls[,'Strike'] #Strike prices for calls 
      kPuts <- optionData$puts[,'Strike']
      #opCalls <-optionData$calls[,'Last']
      #opPuts <- optionData$puts[,'Last']
      opCalls <- (optionData$calls[,'Bid']+optionData$calls[,'Ask'])/2
      opPuts <- (optionData$puts[,'Bid']+optionData$puts[,'Ask'])/2
    }
    if(dataType == "historical"){
      #optionData <- histData[histData$Expiration == expiration1[k],]
      #kCalls <- optionData[optionData$OptionType=='call',][,'Strike']
      #kPuts <- optionData[optionData$OptionType=='put',][,'Strike']
      #opCalls <- (optionData[optionData$OptionType=='call',][,'Bid']+optionData[optionData$OptionType=='call',][,'Ask'])/2
      #opPuts <- (optionData[optionData$OptionType=='put',][,'Bid']+optionData[optionData$OptionType=='put',][,'Ask'])/2
      
      #optionData <- histData[histData$Expiration == expiration1[k],]
      #kCalls <- optionData[optionData$Type=='call',][,'Strike']
      #kPuts <- optionData[optionData$Type=='put',][,'Strike']
      #opCalls <- (optionData[optionData$Type=='call',][,'Bid']+optionData[optionData$Type=='call',][,'Ask'])/2
      #opPuts <- (optionData[optionData$Type=='put',][,'Bid']+optionData[optionData$Type=='put',][,'Ask'])/2
      
      #optionData <- histData[histData$Expiration == expiration1[k],]
      #kCalls <- optionData[optionData$Type=='call',][,'Strike']
      #kPuts <- optionData[optionData$Type=='put',][,'Strike']
      #opCalls <- (optionData[optionData$Type=='call',][,'Bid']+optionData[optionData$Type=='call',][,'Ask'])/2
      #opPuts <- (optionData[optionData$Type=='put',][,'Bid']+optionData[optionData$Type=='put',][,'Ask'])/2
      
      optionData <- histData[histData$expiration == expiration1[k],]
      kCalls <- optionData[optionData$type=='call',][,'strike']
      kPuts <- optionData[optionData$type=='put',][,'strike']
      opCalls <- (optionData[optionData$type=='call',][,'bid']+optionData[optionData$type=='call',][,'ask'])/2
      opPuts <- (optionData[optionData$type=='put',][,'bid']+optionData[optionData$type=='put',][,'ask'])/2
    }
    
    ### Find interest rate
    # Find both with same strike
    opCallsR <- opCalls[is.element(kCalls,kPuts)]
    opPutsR <- opPuts[is.element(kPuts,kCalls)]
    kR <- kCalls[is.element(kCalls,kPuts)]
    r <- mean(-1/T1*log((S-opCallsR+opPutsR)/kR))
    r_list[k] <- r
    
    xCalls <- log(kCalls/S)-r*T1 # negative log forward moneyness
    #xPuts <- log(kPuts/S) -r*T1
    xPuts <- log(kPuts/S) -r*T1
    
    #Check Call-Put parity
    opCallsFromPuts <- opPuts+S*(1-exp(xPuts))
    opPutsFromCalls <- opCalls - S*(1-exp(xCalls))
    
    #Call Plot
    #jpeg("putCall.jpeg", width = 4, height = 4, units = 'in', res = 900)
    plot(kCalls,opCalls, xlim = c(275,490), ylim = c(0.0,134), ylab="Calls", xlab= "strike")
    grid()
    #dev.off()
    #jpeg("putCall2.jpeg", width=4, height=4, units = 'in', res=900)
    plot(kPuts,opCallsFromPuts,xlim = c(275,490),  ylim = c(0.0,134), ylab = "Calls from Puts", xlab="strike")
    grid()
    #dev.off()
    ggCallPutParity <- ggplot()+
      geom_point(data=data.frame(kPuts,opCallsFromPuts), aes(x=kPuts,y=opCallsFromPuts, color = "darkseagreen4"),  shape = 19, size = 4)+xlab("Strike Prices")+theme_bw()+
      geom_point(data=data.frame(kCalls,opCalls), aes(x=kCalls,y=opCalls, color = 'firebrick'), shape =17, size = 2)+ ylab("Call Option Prices")+
      scale_colour_manual(name = '', values =c('darkseagreen4'='darkseagreen4','firebrick'='firebrick'), labels = c('Call Prices From Puts','Call Prices'))
    ggCallPutParity
    #Put Plot
    plot(kPuts,opPuts)
    grid()
    plot(kCalls,opPutsFromCalls)
    grid()
    ggCallPutParity2 <- ggplot()+
      geom_point(data=data.frame(kPuts,opPuts), aes(x=kPuts,y=opPuts, color = "darkseagreen4"), shape = 16, size = 2)+xlab("Strike Prices")+theme_bw()+
      geom_point(data=data.frame(kCalls,opPutsFromCalls), aes(x=kCalls,y=opPutsFromCalls, color = 'firebrick'),shape = 15, size = 2)+ ylab("Put Option Prices")+
      scale_colour_manual(name = '', values =c('darkseagreen4'='darkseagreen4','firebrick'='firebrick'), labels = c('Put Prices','Put Prices from Calls'))
    ggCallPutParity2
    
    
    if(dataType == "live"){
      #For now only Call Options, how can we also use Put Options? Sometimes you get 2 data points when using parity
      #kBoth <- c(kPuts[kPuts<S],kCalls[kCalls>=S])
      kBoth <- kCalls
      #kBoth <- c(kPuts,kCalls)
      #kBoth <- c(kPuts[kPuts>S],kCalls[kCalls<=S])
      #opBoth <- c(opCallsFromPuts[kPuts<S],opCalls[kCalls>=S])
      #opBoth <- c(opCallsFromPuts, opCalls)
      #opBoth <- c(opCallsFromPuts[kPuts>S],opCalls[kCalls<=S])
      opBoth <- opCalls
    }
    if(dataType == "historical"){
      #kBoth <- c(kCalls,kPuts)#, kPuts)
      #opBoth <- c(opCalls,opCallsFromPuts)
      kBoth1 <- c(kCalls, kPuts)
      kBoth <- kBoth1[order(kBoth1)]
      opBoth1 <- c(opCalls,opCallsFromPuts)
      opBoth <- opBoth1[order(opBoth1, decreasing = TRUE)]
      #Puts gebruiken, Calls gebruiken, vraag Sohl?
      #kBoth <- c(kPuts[kPuts<S],kCalls[kCalls>=S])
      #opBoth <- c(opCallsFromPuts[kPuts<S],opCalls[kCalls>=S])
    }
    plot(kBoth,opBoth)
    snop1 <- opBoth/S
    sk1 <- log(kBoth/S)
    plot(sk1, snop1)
    plot(snop1-pmax(0.,1-exp(sk1-r*T1)))
    
    # Calibration
    if (k==1){
      calibrationList <- calibrationHom(sk1,snop1,mode,T1)
      stdNoise1 <- estNoiseHom(sk1,snop1,T1)
    } else{
      T0 <- time[k-1]
      calibrationList <- calibrationInHom(sk1, snop1,T1, sk0, snop0,T0, mode)
      stdNoise1 <- estNoiseInHom(sk1, snop1,T1, sk0, snop0,T0)
      v0 <- calibrationList$v0
      phiNum0 <- calibrationList$phiNum0
      phiNu0 <- calibrationList$phiNu0
    }
    U[k] <- calibrationList$U
    UNu[k] <- calibrationList$UNu
    sigmaHat[k] <- calibrationList$sigmaHat
    gammaHat[k] <- calibrationList$gammaHat
    lambdaHat[k] <- calibrationList$lambdaHat
    x[k,] <- as.vector(calibrationList$x)
    nuHat[k,] <- as.vector(calibrationList$nuHat)
    v <- calibrationList$v
    phiNum <- calibrationList$phiNum
    phiNu <- calibrationList$phiNum
    
    #Variance of parameters:
    
    #stdNoise1 <- (snop1-pmax(0.,1-exp(sk1-r*T1)))*0.05 #Cont and Tankov, p439, However I want to find a better noise estimate
    confInt <- confIntervals(k,1)
    varSigma2Hat[k] <- confInt$varSigma2
    varGammaHat[k] <- confInt$varGamma
    varLambdaHat[k] <- confInt$varLambda
    varNuHat[k,] <- confInt$varNu
     
    sk0 <- sk1
    snop0 <- snop1
    stdNoise0 <- stdNoise1
    assign(paste("estNoise",k),stdNoise1)
  }
  
  # Plots of Estimators
  expirationDay <- str_sub(expiration, start= -5) #Displays days of expiration, easy for plots
  df <- data.frame(expirationDay, sigmaHat,gammaHat,lambdaHat,varSigma2Hat,varGammaHat,varLambdaHat)
  df$expirationDay<- factor(df$expirationDay, expirationDay)
  
  ggSigma <- ggplot(df, aes(x = expirationDay ,y= sigmaHat)) + geom_point(alpha=1/3,colour = "black") +
    theme_bw()+ylab(TeX("volatility $\\sigma$")) +xlab('Monte Carlo Simulations')
  ggSigma <- ggMarginal(ggSigma, type="boxplot", fill="gray", col = "black", margins="y")
  
  ggGamma <- ggplot(df, aes(x = expirationDay,y= gammaHat)) + geom_point(alpha=1/3,colour = "black") +
    theme_bw()+ylab(TeX("drift $\\gamma$")) +xlab('Monte Carlo Simulations') 
  ggGamma <- ggMarginal(ggGamma, type="boxplot", fill="gray", col = "black", margins="y")
  
  ggLambda <- ggplot(df, aes(x = expirationDay,y= lambdaHat)) + geom_point(alpha=1/3,colour = "black") +
    theme_bw()+ylab(TeX("intensity $\\lambda$")) +xlab('Monte Carlo Simulations') 
  ggLambda <- ggMarginal(ggLambda, type="boxplot", fill="gray", col = "black", margins="y")
  
  grid.arrange(ggSigma, ggGamma,ggLambda, nrow=3)
  
  ###errorbar Plots
  ggSigmaErr <- ggplot(df, aes(x =  expirationDay, y= sigmaHat^2)) +
    theme_bw()+ylab(TeX("volatility $\\tilde{\\sigma}_j^2$")) +xlab('')+
    geom_errorbar(aes(ymin=sigmaHat^2+lowerQ*sqrt(varSigma2Hat), ymax=sigmaHat^2+upperQ*sqrt(varSigma2Hat)),width=0.1, colour="aquamarine4", size=0.8)+
    geom_point(colour = "black", size=2)
  #ggsave("sigmaEmp.png", plot=ggSigmaErr, dpi=500, width=4, height =4)
  ggGammaErr <- ggplot(df, aes(x =  expirationDay, y= gammaHat)) +
    theme_bw()+ylab(TeX("drift $\\tilde{\\gamma}_$")) +xlab('')+
    geom_errorbar(aes(ymin=gammaHat+lowerQ*sqrt(varGammaHat), ymax=gammaHat+upperQ*sqrt(varGammaHat)),width=0.1, colour="aquamarine4", size=0.8)+
    geom_point(colour = "black", size=2)
  #ggsave("gammaEmp.png", plot=ggGammaErr, dpi=500, width=4, height =4)
  ggLambdaErr <- ggplot(df, aes(x =  expirationDay, y= lambdaHat)) +
    theme_bw()+ylab(TeX("intensity $\\tilde{\\lambda}_j$")) +xlab('')+
    geom_errorbar(aes(ymin=lambdaHat+lowerQ*sqrt(varLambdaHat), ymax=lambdaHat+upperQ*sqrt(varLambdaHat)),width=0.1, colour="aquamarine4", size=0.8)+
    geom_point(colour = "black", size=2)
  #ggsave("lambdaEmp.png", plot=ggLambdaErr, dpi=500, width=4, height =4)
  grid.arrange(ggSigmaErr,ggGammaErr,ggLambdaErr)
  
  
  #Plots of Levy densities
  gg <- list()
  for(i in 1:20){
    ggNew <- ggplot(data = data.frame(x=x[i,], y=nuHat[i,],z=varNuHat[i,]), aes(x,y))+geom_line()+ theme_bw()+
      geom_line(aes(x,y+upperQ*sqrt(z)),colour = "gray50", lty=2)+
      geom_line(aes(x,pmax(0,y+lowerQ*sqrt(z))),colour = "gray50", lty=2)+xlim(-2.0,0.5)+
      labs(x = TeX("$x$"), y = TeX("$\\hat{\\nu}(x)$"), title = df$expirationDay[i])
    #ggNew
    #ggsave("NuEmp3.png", plot=ggNew, dpi=500, width=4, height =4)
    gg <- c(gg, list(ggNew))
  }
  ggNu <- do.call("grid.arrange", c(gg, nrow =4 , ncol =5))
  ggNu
  #ggsave("NuTime.png", plot=ggNu, dpi=500, width=12, height =12)
  
  frame()
  par(mfrow=c(1,1))
  matplot(t(x),t(nuHat/df$lambdaHat), type = "l", xlab= "x", ylab=TeX("$\\hat{\\nu}(x)/\\tilde{\\lambda}$"), xlim = c(-2.0,0.7)) 
  grid()
  
  #Generated Paths
  NoOfSteps <- 1e4
  paths <- generatePaths(r,NoOfSteps)
  t <- paths$t
  St <- paths$St
  tDate <- str_sub(as.Date(t*52*365, origin =  Sys.Date()), start = -5)
  ggPath <- ggplot(data = data.frame(t,St), aes(x=t*52,y=St)) +geom_line()+theme_bw()+ylab("S(t)  (dollars)") +xlab('t (weeks)')+
    geom_vline(xintercept = c(0,as.numeric(difftime(expiration, quoteDate, units = "weeks"))/52)*52, 
               linetype="solid", color = "grey30", size=0.5)
  ggPath
  #ggsave("pathTime1.png", plot=ggPath, dpi=500, width=4, height =4)
  
  
  ggSigmaTime <- ggplot(df, aes(x=expirationDay, y=sigmaHat^2, group=1)) + theme_bw()+geom_line() + 
    geom_point()+geom_line(aes(x = expirationDay, y = sigmaHat^2+upperQ*sqrt(varSigma2Hat)), colour = 'gray')+
    geom_line(aes(x = expirationDay, y = pmax(0,sigmaHat^2+lowerQ*sqrt(varSigma2Hat))), colour = 'gray') +
    labs(x=TeX("$t$"),y = TeX("$ \\tilde{\\sigma}_j^2(t)$"))
  ggSigmaTime
  #ggsave("sigmaTime.png", plot=ggSigmaTime, dpi=500, width=10, height =4)
  
  ggGammaTime <- ggplot(df, aes(x=expirationDay, y=gammaHat, group=1)) + theme_bw()+geom_line() + 
    geom_point()+geom_line(aes(x = expirationDay, y = gammaHat+upperQ*sqrt(varGammaHat)), colour = 'gray')+
    geom_line(aes(x = expirationDay, y = gammaHat+lowerQ*sqrt(varGammaHat)), colour = 'gray') +
    labs(x=TeX("$t$"),y = TeX("$ \\tilde{\\gamma}_j(t)$"))
  ggGammaTime
  #ggsave("gammaTime.png", plot=ggGammaTime, dpi=500, width=10, height =4)
  
  ggLambdaTime <- ggplot(df, aes(x=expirationDay, y=lambdaHat, group=1)) + theme_bw()+geom_line() + 
    geom_point()+geom_line(aes(x = expirationDay, y = lambdaHat+upperQ*sqrt(varLambdaHat)), colour = 'gray')+
    geom_line(aes(x = expirationDay, y = lambdaHat+lowerQ*sqrt(varLambdaHat)), colour = 'gray') +
    labs(x=TeX("$t$"),y = TeX("$ \\tilde{\\lambda}_j(t)$"))
  ggLambdaTime
  #ggsave("lambdaTime.png", plot=ggLambdaTime, dpi=500, width=10, height =4)
  
  library(plotly)
  axx <- list(nticks = 4,range = c(-4,1))
  axy <- list(nticks = 4,range = c(0,16))
  axz <- list(nticks = 4,range = c(0,130))
  plotly_figure <- plot_ly(x = x, y = as.numeric(df$expirationDay), z = nuHat, type= 'surface')
  plotly_figure <- plotly_figure  %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  plotly_figure
}
# Export PDF
#pdf(file = "histData.pdf", height = 5, width = 13)
grid.table(histData[250:255,c(1,2,6,7,8,9,10,11,12,13)])
#dev.off()

