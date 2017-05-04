require(coda)
require(doParallel)
require(ggplot2)
require(foreach)
require(tidyverse)
require(plyr)
require(dplyr)
require(tidyr)

###### Plotting Functions

## Plots GR Diagnostic Values
Chains <- Result1
Chains2 <- Result2
GelmanPlot <- function(Chains, Chains2, pargroup = 1, parnames, name1, name2, p){
  
  #It depends on if you have one group or two how the list structure is
  if(pargroup != 1){
    G <- Chains[[3]][[pargroup]]
    G2 <- Chains2[[3]][[pargroup]]
    p <- dim(Chains[[1]][[1]][[1]])[2]
  }
  
  #Collecting Chains GR Info
  if(pargroup == 1){
    G <- Chains[[3]]
    G2 <- Chains2[[3]]
    p <- dim(as.matrix(Chains[[1]][[1]]))[2]
  
  }
  #Bind them into one data frame
  Bout <- as.data.frame(rbind(G, G2))
  
  #Include Indicator for which model it came from
  Bout <- as.data.frame(cbind(Bout, rep(c(name1, name2), each = dim(G)[1])))
  
  #Naming model column appropriately
  colnames(Bout)[length(colnames(Bout))] <- "Method" 

  ## Breaking them up by parameter
  kk <- kkFunc(Bout, p, parnames, name1, name2)
  
  #Plotting
  z <- plotGR(kk, p)
  
  return(z)
  
}


## Plots ACF Values
ACFPlot <- function(Chains, Chains2){
  #Binding together the ACF values
  ACF <- rbind(Chains$ACF, Chains2$ACF)

  P <- ggplot(ACF, aes(x = Lag, y = values, colour = Method)) + ylab("") + geom_line(aes(linetype = Method)) + facet_grid(Chain ~ Beta, labeller = label_parsed) + ggtitle("ACF D2") +  theme_classic() + scale_x_continuous(breaks = c(0, 20, 40)) + scale_linetype_manual(values = c(1,3)) + scale_fill_manual(values=c("red", "dark green")) + theme(strip.text=element_text(size=10), axis.text = element_text(size = 10), plot.title = element_text(size = 14) )

  return(P)
  
}

#################################### MCMC


####ASIS Probit
#Initial is one starting set
#K is the length of chain
#p is the number of parameters
#x and y are the data

AProbit <- function(initial, K, p, x, y){
  
  #Turning data into -1 and 1, as seen later why.
  
  ys <- y
  for(i in 1:length(y)){
    if(y[i] == 0){ys[i] <- -1}
  }
  
  

  n <- length(y)
  eta <- NA
  z <- NA
  phi <- NA
  lambda <- NA
  phi <- NA

  ASISprobit <- array(-999, dim = c(K+1, p))
  
  #Required libraries
  require(tmvtnorm)
  require(MASS)
  require(truncnorm)
  
  #Instigatng Initial Values
  Beta <- as.matrix(initial)
  
  
  ASISprobit[1,] <- Beta
  
  for(k in 1:K){
      #Finding bounds for phi
      bound <- bounds2(y, Beta, x, n)
      
      phi <- rtruncnorm(1, mean = c((x%*%Beta)), sd = 1, a = bound[,1], b = bound[,2])
    
    
    for(j in 1:p){
      
      
      eta <- phi - x%*%Beta
      
      minim <- NA
      
      maxim <- NA
      
      mincount <- 0
      
      maxcount <- 0
      
      #Allows an easy mathematical expression for bounds for Beta
      
      minims <- which(ys*x[,j] > 0)
      maxims <- which(ys*x[,j] < 0)
      
      vals <- -(eta+x[,-j]%*%Beta[-j,])/(x[,j])
      
      minim <- max(vals[minims])
      maxim <- min(vals[maxims])
      
      Beta[j] <- runif(1,minim,maxim)
      
    }
    
    ASISprobit[k+1,] <- Beta
  }
  
  return(ASISprobit)
}

#### sandwich Probit
M <- 3
K <- 60
x <- X
initial <- rnorm(3)
SandProbit <- function(initial, M, K, p, x, y){
  
  #Initializing values
  n <- length(y)
  Sand2Beta <- array(-9999, dim = c(K+1, p))
  Mu <- c(0,0,0)
  n <- length(y)
  mu <- c(0.0,0)
  eta <- NA
  z <- NA
  phi <- NA
  lambda <- NA
  g <- NA
  phi <- NA
  
  #Libraries required
  require(tmvtnorm)
  require(MASS)
  require(truncnorm)
  
  
  Beta <- as.matrix(initial)
  
  Sand2Beta[1,] <- Beta
  
  for(k in 1:K){
    
    bounds <- bounds2(y, Beta, x, n)
    
    z <- rtruncnorm(1, bounds[,1], bounds[,2], x%*%Beta, 1)
    
    b <- 1/2*sum((z[i]-t(x[i,])%*%solve(t(x)%*%x)%*%t(x)*z)^2)
    
    g <- (rgamma(1, n/2, b))^.5
    
    zprime <- g*z
    
    Betazprime <- solve(t(x)%*%x)%*%t(x)%*%as.matrix(zprime)
    
    Beta <- t(as.matrix(rtmvnorm(1, c(Betazprime), solve(t(x)%*%x))))
    
    Sand2Beta[k+1,] <- Beta
  }
  
  return(Sand2Beta)
}


##### ASIS Robit Model Algorithm



ASISRobit <- function(initial, M, K, p, x, y, nu, nunot, c){
  
  ### See Sand for notes.
  ys <- y
  for(i in 1:length(y)){
    if(y[i] == 0){ys[i] <- -1}
  }
  
  n <- length(y)
  Mu <- c(0,0,0)
  Sigmanot <- c*t(x)%*%x
  Betastore <- array(-9999, dim = c(K+1, p))
  phi <- NA
  lambda <- NA
  require(tmvtnorm)
  require(MASS)
  require(truncnorm)
  
  Beta <- as.matrix(initial)
  
  zz <- 1
  
  Betastore[zz,] <- Beta
  zz <- zz + 1 
  
  for(i in 1:K){
    
    
    low <- -Inf
    high <- Inf
    
    bound <- bounds(y, Beta, x, n)
    
    eta <- rtmvt(1, mean = rep(0, n), lower = bound[,1], upper = bound[,2], algorithm = "gibbs")
    
    a <- (nu+1)/2 
    b <- (nu+(eta^2))/2
    
    lambda <- apply(as.matrix(b),2,gammaout, a = a)
    
    
    
    atau <- (nunot + p)/2
    btau <- (nunot + t(Beta)%*%Sigmanot%*%Beta)/2
    
    tau <- rgamma(1, atau, btau)
    
    
    
    for(j in 1:p){
      
      minims <- which(ys*x[,j] > 0)
      maxims <- which(ys*x[,j] < 0)
      
      vals <- -(eta+x[,-j]%*%Beta[-j,])/(x[,j])
      
      minim <- max(vals[minims])
      maxim <- min(vals[maxims])
      
      out <- conditionalmeansig(j,solve(tau*Sigmanot), Mu, Beta)
      
      condmu <- out[[1]]
      
      condsig <- out[[2]]
      
      Beta[j,] <- rtruncnorm(1, a = minim, b = maxim, mean = condmu, sd = condsig)
      
    }
    
    L <- diag(lambda)
    
    z <- eta + x%*%Beta
    
    Beta <- as.matrix(mvrnorm(n = 1, mu = ginv(t(x)%*%L%*%x+tau*Sigmanot)%*%t(x)%*%L%*%z, Sigma = ginv(t(x)%*%L%*%x+tau*Sigmanot)))
    
    Betastore[zz,] <- Beta
    
    zz <- zz + 1
    
  }
  
  
  return((Betastore))
  
}


####### Sandwich Robit Algorithm

SandRobit <- function(initial, M, K, p, x, y, nu, nunot, c){
  z <- NA
  Sigmanot <- c*t(x)%*%x
  
  Betastore <- array(-9999, dim = c(K+1, p))
  phi <- NA
  eta <- NA
  lambda <- NA
  require(tmvtnorm)
  require(MASS)
  require(truncnorm)
  g <- NA
  count <- 0
  n <- length(y)
  
  Beta <- as.matrix(initial)
  
  zz <- 1
  
  Betastore[zz,] <- Beta
  
  zz <- zz + 1
  
  for(i in 1:K){
    
    bound <- bounds2(y, Beta, x, n)
    
    z <- rtmvt(1, mean = c(x%*%Beta), df = 3, lower = bound[,1], upper = bound[,2], algorithm = "gibbs")
    
    a <- (nu+1)/2 
    
    b <- (nu+(z-(x%*%Beta)^2))/2
    
    lambda <- apply(b, 1, gammaout, a = a)
    
    atau <- (nunot + p)/2
    btau <- (nunot + t(Beta)%*%Sigmanot%*%Beta)/2
    
    tau <- rgamma(1, atau, btau)
    
    L <- diag(lambda)
    
    Q <- L^.5%*%x%*%solve(t(x)%*%L%*%x+tau*Sigmanot)%*%t(x)%*%L^.5
    
    I <- diag(rep(1, dim(Q)[1]))
    
    
    g <- (rgamma(1,n/2, t(z)%*%L^.5%*%(I - Q)%*%L^.5%*%z/2))^.5
    
    
    
    
    Beta <- as.matrix(mvrnorm(n = 1, mu = g*solve(t(x)%*%L%*%x+tau*Sigmanot)%*%t(x)%*%L%*%z, Sigma = solve(t(x)%*%L%*%x+tau*Sigmanot)))
    
    Betastore[zz,] <- Beta
    
    zz <- zz + 1
  }
  
  
  return((Betastore))
  
}


### Conditional Mean and Variance for a normal distribution

conditionalmeansig <- function(i, Sigma, mu,Beta){
  Sigma11 <- Sigma[i,i]
  Sigma12 <- Sigma[i,]
  Sigma21 <- Sigma[-i,i]
  Sigma22 <- Sigma[-i,-i]
  Sigma12 <- Sigma[i,-i]
  
  condmu <- mu[i]+Sigma12%*%solve(Sigma22)%*%(Beta[-i,]-as.matrix(mu[-i]))
  
  condsig <- Sigma11 - Sigma12%*%solve(Sigma22)%*%Sigma21
  return(list(condmu, condsig))
}




############ Initial Value Generator


InitialGenerator <- function(Betastore, M, p){
  seeds <- list(1)
  
  sdB1 <- sd(Betastore[,1])
  sdB2 <- sd(Betastore[,2])
  sdB3 <- sd(Betastore[,3])
  
  mnB1 <- mean(Betastore[,1])
  mnB2 <- mean(Betastore[,2])
  mnB3 <- mean(Betastore[,3])
  
  Beta.Var <- var(Betastore)
  
  for(i in 1:M){
    seeds[[i]]  <-  rmvt(1, sigma = Beta.Var, delta = c(mnB1, mnB2,mnB3), df = c(8,8,8))
  }
  
  
  return(seeds) 
}

initfunction <- function(Initials, Model, MCMC, M, K, c, nu, nunot, x, y,p){
  
  if(Model == "ASISProbit"){
    out <- AProbit(Initials, K, p, x, y)
  }
  
  if(Model == "SandProbit"){
    out <- SandProbit(Initials, 1, K, p, x, y)  
  }
  
  if(Model == "ASISRobit"){
    out <- ASISRobit(Initials, 1, K, p, x, y, 3, 3, 1)  
  }
  
  if(Model == "SandRobit"){
    out <- SandRobit(Initials, 1, K, p, x, y, 3,3, 1)  
  }
  
  if(Model == "NA"){
    out <- MCMC(Initials, M = 1, ...)
  }
  
  x <- InitialGenerator(out, M, p)
  return(x)
}



#### Aggregate Parameters 

BetaAggregate <- function(BetaChains, name, p = 3){
  
  x <- sapply(BetaChains, ACFWrapper)
  
  x.2 <- data.frame(x, rep(paste("Par[", 1:p-1, "]", sep = ""), each = length(x[,1])/p))
  
  colnames(x.2) <- c(paste("C", 1:length(BetaChains), sep = ""), "Beta")
  
  x.2 <- x.2 %>% gather(key = Beta, value = values)
  
  colnames(x.2)[2] <- "Chain"
  
  x.2 <- data.frame(x.2, paste(name))
  
  colnames(x.2) <- c(colnames(x.2)[-length(colnames(x.2))], "Method")
  
  x.2 <- data.frame(x.2, 1:(dim(x.2)/(length(BetaChains)*p))[1])
  
  colnames(x.2) <- c(colnames(x.2)[-length(colnames(x.2))], "Lag")
  
  return(x.2)
  
}



### Inner function for ACF

ACF.Gather <- function(ChainList){
  Chains <- sapply(ChainList, ACFWrapper)
  
  Chains.2 <- data.frame(Chains, rep(paste("Beta", 1:p-1), each = length(Suff.ACF[,1])/p))
  
  colnames(Chains.2) <- c(paste("Chain", 1:length(ChainList), sep = "."), "Beta")
  
  return(Chains.2)
}

### Picks up the first k values for the GR Diagnostic

reduce <- function(x,k){
  x <- x[1:k,]
  return(x)
}

### Simple mean calculator
mean.func <- function(x){
  out <- apply(x, 2, mean)
}

### Simple variance calculator
vars <- function(x){
  out <- apply(x, 2, var)
  return(out)
}

### Calculates GR Diagnostic values

Gelman <- function(data, M, k){
  N <- k
  
  data <- lapply(data, reduce, k = k)
  
  Thetas <- sapply(data, mean.func)
  
  Thetahat <- apply(Thetas,1,mean)
  
  B <-  N/(M-1)*apply((Thetas-Thetahat)^2, 1, sum)
  
  varvals <- sapply(data,vars)
  
  W <- apply(varvals,1,mean)
  
  V <- (N-1)/N*W + (M+1)/(M*N)*B
  
  scale <- (V/W)^.5
  
  return(scale)
}


#Changes a list to dataframe

list.to.dataframe <- function(x){
  c <- length(x)
  cbind(list[[]])
}

#An inner, inner, inner ACF calculation. 

ACF.Inner <- function(x, plot){
  out <- acf(x, plot = plot)$acf
  return(out)
}

#An Inner inner ACF calculation

ACFWrapper <- function(x){
  out <- apply(x,2,ACF.Inner, plot = FALSE)
  return(out)
}

#Calculates Gelman values for each cutoff values ks

GelmanWrapper <- function(data, M, ks){
  out <- array(-999, dim = c(3,length(ks)))
  for(i in 1:length(ks)){
    out[,i] <- Gelman(data, M, ks[i])
  }
  
  colnames(out) <- paste(ks)
  return(out)
}




##### Geweke functions

Gewekedataframe2 <- function(x, frac1, frac2, beta, Method, chains){
  z <- 0
  Geweke <- data.frame(0,0,0,0,0,0)
  for(i in 1:1){
    for(j in 1:1){
      for(k in 1:3){
        for(l in 1:2){
          for(m in 1:5){
            z <- z + 1
            Geweke[z,] <- c(x[i,j,k,l,m], frac1[i], frac2[j], beta[k], Method[l], chains[m])
          }
        }
      }
    }
    
  }
  
  colnames(Geweke) <- c("Values", "Init.Fraction", "End.Fraction", "Beta", "Method", "Chain")
  
  return(Geweke)
}



geweke <- function(x, frac1 = frac1, frac2 = frac2){
  
  out <- array(-999, dim = c(1,1,3))
  
  for(i in 1:1){
    for(j in 1:1){
      out[i,j,] <- geweke.diag(x)$z
    }
  }
  
  return(out)  
}


gewekewrapper <- function(x, frac1, frac2){
  out <- array(-999, dim = c(1,1,3,4,5))
  
  for(i in 1:4){
    for(j in 1:5){
      out[,,,i,j] <- geweke(x[[i]][[j]],frac1, frac2)
    }
  }
  return(out)
}


gewekewrapper2 <- function(x, frac1, frac2, M, p){
  out <- array(-999, dim = c(p,M))
  
  for(j in 1:M){
    out[,j] <- geweke(x[[1]][[j]],frac1, frac2)
  }
  
  return(out)
}



##########

#Creates Mcmclist for coda package for Geweke values.

MCMCfunc <- function(out, M){
  mcmclist <- list(1)
  if(M == 3){
    mcmclist[[1]] <- mcmc.list(mcmc(out[[1]]), mcmc(out[[2]]), mcmc(out[[3]]))
  }
  if(M==4){
    mcmclist[[1]] <- mcmc.list(mcmc(out[[1]]), mcmc(out[[2]]), mcmc(out[[3]]), mcmc(out[[4]]))
  }
  if(M == 5){
    mcmclist[[1]] <- mcmc.list(mcmc(out[[1]]), mcmc(out[[2]]), mcmc(out[[3]]), mcmc(out[[4]]), mcmc(out[[5]]))
  }
  if(M == 6){
    mcmclist[[1]] <- mcmc.list(mcmc(out[[1]]), mcmc(out[[2]]), mcmc(out[[3]]), mcmc(out[[4]]), mcmc(out[[5]]), mcmc(out[[6]]))
  }
  return(mcmclist)
}


### Isolates the parameters into separate list components for graphing GR diagnostics


kkFunc <- function(Bout, p, parnames, method1, method2){
  
  
  if(p == 1){
    kk <- data.frame(c(Bout$X1), Bout$Iteration, Bout$Method, rep(paste("Par", 1:p, sep = "."), each = length(Bout$X1)))
    
    colnames(kk) <- c("Value", "Iteration", "Method", "Par")
    
    kk$Iteration <- as.numeric(as.character(kk$Iteration))
    
    kk.B1 <- kk[kk$Par == "Par.1",]
    
    kk.B1$Method <- factor(kk.B1$Method, levels = c(method1, method2))
    
    
    kk <- list(kk.B1)
  }  
  
  
  if(p == 2){
    kk <- data.frame(c(Bout$X1, Bout$X2), Bout$Iteration, Bout$Method, rep(paste("Par", 1:p, sep = "."), each = length(Bout$X1)))
    
    colnames(kk) <- c("Value", "Iteration", "Method", "Par")
    
    kk$Iteration <- as.numeric(as.character(kk$Iteration))
    
    kk.B1 <- kk[kk$Par == "Par.1",]
    
    kk.B2 <- kk[kk$Par == "Par.2",]
    
    
    kk.B1$Method <- factor(kk.B1$Method, levels = c(method1, method2))
    kk.B2$Method <- factor(kk.B2$Method, levels = c(method1, method2))
    kk <- list(kk.B1, kk.B2)
  }  
  
  if(p == 3){
    kk <- data.frame(c(Bout$X1, Bout$X2, Bout$X3), Bout$Iteration, Bout$Method, rep(paste("Par", 1:p, sep = "."), each = length(Bout$X1)))
    
    colnames(kk) <- c("Value", "Iteration", "Method", "Par")
    
    kk$Iteration <- as.numeric(as.character(kk$Iteration))
    
    kk.B1 <- kk[kk$Par == "Par.1",]
    
    kk.B2 <- kk[kk$Par == "Par.2",]
    
    kk.B3 <- kk[kk$Par == "Par.3",]
    
    kk.B1$Method <- factor(kk.B1$Method, levels = c(method1, method2))
    kk.B2$Method <- factor(kk.B2$Method, levels = c(method1, method2))
    kk.B3$Method <- factor(kk.B3$Method, levels = c(method1, method2))
    kk <- list(kk.B1, kk.B2, kk.B3)
  }
  
  
  if(p == 4){
    kk <- data.frame(c(Bout$X1, Bout$X2, Bout$X3, Bout$X4), Bout$Iteration, Bout$Method, rep(paste("Par", 1:p, sep = "."), each = length(Bout$X1)))
    
    colnames(kk) <- c("Value", "Iteration", "Method", "Par")
    
    kk$Iteration <- as.numeric(as.character(kk$Iteration))
    kk.B1 <- kk[kk$Par == "Par.1",]
    
    kk.B2 <- kk[kk$Par == "Par.2",]
    
    kk.B3 <- kk[kk$Par == "Par.3",]
    
    kk.B4 <- kk[kk$Par == "Par.4",]
    
    kk.B1$Method <- factor(kk.B1$Method, levels = c(method1, method2))
    kk.B2$Method <- factor(kk.B2$Method, levels = c(method1, method2))
    kk.B3$Method <- factor(kk.B3$Method, levels = c(method1, method2))
    kk.B4$Method <- factor(kk.B4$Method, levels = c(method1, method2))
    
    kk <- list(kk.B1, kk.B2, kk.B3, kk.B4)
  }  
  
  
  if(p == 5){
    kk <- data.frame(c(Bout$X1, Bout$X2, Bout$X3, Bout$X4, Bout$X5), Bout$Iteration, Bout$Method, rep(paste("Par", 1:p, sep = "."), each = length(Bout$X1)))
    
    colnames(kk) <- c("Value", "Iteration", "Method", "Par")
    
    kk$Iteration <- as.numeric(as.character(kk$Iteration))
    
    kk.B1 <- kk[kk$Par == "Par.1",]
    
    kk.B2 <- kk[kk$Par == "Par.2",]
    
    kk.B3 <- kk[kk$Par == "Par.3",]
    
    kk.B4 <- kk[kk$Par == "Par.4",]
    
    kk.B5 <- kk[kk$Par == "Par.5",]
    
    
    
    kk.B1$Method <- factor(kk.B1$Method, levels = c(method1, method2))
    kk.B2$Method <- factor(kk.B2$Method, levels = c(method1, method2))
    kk.B3$Method <- factor(kk.B3$Method, levels = c(method1, method2))
    kk.B4$Method <- factor(kk.B4$Method, levels = c(method1, method2))
    kk.B5$Method <- factor(kk.B5$Method, levels = c(method1, method2))
    
    
    
    kk <- list(kk.B1, kk.B2, kk.B3, kk.B4, kk.B5)
  }  
  
  if(p == 6){
    kk <- data.frame(c(Bout$X1, Bout$X2, Bout$X3, Bout$X4, Bout$X5, Bout$X6), Bout$Iteration, Bout$Method, rep(paste("Par", 1:p, sep = "."), each = length(Bout$X1)))
    
    colnames(kk) <- c("Value", "Iteration", "Method", "Par")
    
    kk$Iteration <- as.numeric(as.character(kk$Iteration))
    
    kk.B1 <- kk[kk$Par == "Par.1",]
    
    kk.B2 <- kk[kk$Par == "Par.2",]
    
    kk.B3 <- kk[kk$Par == "Par.3",]
    
    kk.B4 <- kk[kk$Par == "Par.4",]
    
    kk.B5 <- kk[kk$Par == "Par.5",]
    
    kk.B5 <- kk[kk$Par == "Par.6",]
    
    kk.B1$Method <- factor(kk.B1$Method, levels = c(method1, method2))
    kk.B2$Method <- factor(kk.B2$Method, levels = c(method1, method2))
    kk.B3$Method <- factor(kk.B3$Method, levels = c(method1, method2))
    kk.B4$Method <- factor(kk.B4$Method, levels = c(method1, method2))
    kk.B5$Method <- factor(kk.B5$Method, levels = c(method1, method2))
    kk.B6$Method <- factor(kk.B6$Method, levels = c(method1, method2))   
    
    
    kk <- list(kk.B1, kk.B2, kk.B3, kk.B4, kk.B5, kk.B6)
  } 
  return(kk)
}




### Plots GR Diagnostic values

plotGR <- function(kk, p){
  
  
  if(p == 1){
    kk.B1 <- kk[[1]]
    
    P1 <- ggplot(kk.B1, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 1"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    return(P1)
  }
  
  if(p == 2){
    kk.B1 <- kk[[1]]
    kk.B2 <- kk[[2]]
    
    
    P1 <- ggplot(kk.B1, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 1"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P2 <- ggplot(kk.B2, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 2"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    
    return(list(P1, P2))  
  }
  
  if(p == 3){
    
    
    kk.B1 <- kk[[1]]
    kk.B2 <- kk[[2]]
    kk.B3 <- kk[[3]]
    
    P1 <-ggplot(kk.B1, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 1"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P2 <- ggplot(kk.B2, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 2"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P3 <- ggplot(kk.B3, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 3"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    
    return(list(P1,P2,P3))
  }
  
  if(p == 4){
    kk.B1 <- kk[[1]]
    kk.B2 <- kk[[2]]
    kk.B3 <- kk[[3]]
    kk.B4 <- kk[[4]]
    
    P1 <-ggplot(kk.B1, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 1"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P2 <- ggplot(kk.B2, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 2"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P3 <- ggplot(kk.B3, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 3"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P4 <- ggplot(kk.B4, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 4"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    return(list(P1, P2, P3, P4))  
    
  }
  
  if(p == 5){
    kk.B1 <- kk[[1]]
    kk.B2 <- kk[[2]]
    kk.B3 <- kk[[3]]
    kk.B4 <- kk[[4]]
    kk.B5 <- kk[[5]]
    
    
    P1 <-ggplot(kk.B1, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 1"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P2 <- ggplot(kk.B2, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 2"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P3 <- ggplot(kk.B3, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 3"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P4 <- ggplot(kk.B4, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 4"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P5 <- ggplot(kk.B5, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 5"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    
    return(list(P1, P2, P3, P4, P5))
  }
  
  if(p == 6){
    kk.B1 <- kk[[1]]
    kk.B2 <- kk[[2]]
    kk.B3 <- kk[[3]]
    kk.B4 <- kk[[4]]
    kk.B5 <- kk[[5]]
    kk.B6 <- kk[[6]]
    
    P1 <-ggplot(kk.B1, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 1"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P2 <- ggplot(kk.B2, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 2"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P3 <- ggplot(kk.B3, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 3"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P4 <- ggplot(kk.B4, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 4"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P5 <- ggplot(kk.B5, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 5"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    P6 <- ggplot(kk.B6, aes(x = Iteration, y = Value, colour = Method)) + geom_line(aes(linetype = Method)) + ggtitle(expression(paste("GR Diag for Par 6"))) + ylab("GR") + theme_classic() + scale_linetype_manual(values = c(1,3,4,5)) + scale_colour_manual(values=c("red", "dark green"))
    
    return(list(P1, P2, P3, P4, P5, P6))
  }
  
  
  
}


### calculates bounds for MCMC algorithm methods

bounds2 <- function(y, Beta, x, n){
  bounds <- array(-9999, dim = c(n,2))
  zeroes <- which(y == 0)
  
  bounds[zeroes,] <-  cbind(rep(-Inf,length(zeroes)),rep(0, length(zeroes)))
  
  
  ones <- which(y == 1)
  
  bounds[ones,] <- cbind(rep(0, length(ones)),rep(Inf, length(ones)))
  
  return(bounds)
}

#Same as above

bounds <- function(y, Beta, x,n){
  bounds <- array(-9999, dim = c(n,2))
  zeroes <- which(y == 0)
  high <- -(x%*%Beta)[zeroes]
  bounds[zeroes,] <-  cbind(rep(-Inf,length(zeroes)),high)
  
  
  ones <- which(y == 1)
  low <- -(x%*%Beta)[ones]
  bounds[ones,] <- cbind(low,rep(Inf, length(ones)))
  
  return(bounds)
}

### Gamma function used within apply to avoid loops

gammaout <- function(b, a){
  lambdatemp <- rgamma(1, a, b)
  return(lambdatemp)
}
