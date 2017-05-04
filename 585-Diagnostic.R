Diagnostic <- function(Initials = NA, MCMC = NA, M = 5, K, x, y, Initialize = 0, Model = NA, nu = 3, nunot = 3, c = 1, pargroups = 1, Methodname, GewInit = .1, GewEnd = .5, export = NA, ...){
  
  n <- length(y)
  p <- dim(x)[2]
  if(Initialize == 1){
    if(is.na(M) == T){
      M <- 5
    }
    initial <- as.matrix(rnorm(dim(x)[2]))
    Initials <- initfunction(initial, Model, MCMC, M, K, c, nu, nunot, x, y, p)
  }
  
  if(Model == "ASISProbit"){
    out <- foreach(i = 1:M, .export = c("AProbit", "bounds2", "bounds")) %dopar% {AProbit(t(as.matrix(Initials[[i]])), K, p, x, y)}
  }
  
  
  if(Model == "SandProbit"){
    out <- foreach(i = 1:M, .export = c("SandProbit", "bounds2", "bounds")) %dopar% SandProbit(t(as.matrix(Initials[[i]])), M, K, p, x, y)
  }
  
  if(Model == "ASISRobit"){
    out <- foreach(i = 1:M, .export = c("ASISRobit", "conditionalmeansig", "bounds2", "bounds", "gammaout")) %dopar% ASISRobit(t(as.matrix(Initials[[i]])), M, K, p, x, y, nu, nunot, c)  
  }
  
  if(Model == "SandRobit"){
    out <- foreach(i = 1:M, .export = c("SandRobit", "conditionalmeansig", "bounds", "bounds2", "gammaout")) %dopar% SandRobit(t(as.matrix(Initials[[i]])), M, K, p, x, y, nu, nunot, c)  
  }
  
  if(Model == "NA"){
    out <- foreach(i=1:M, .export = export) %dopar% MCMC(Initials[[i]], M= M, K = K, x, y, ...)
  }
  
  
  
  for(i in 1:M){
    colnames(out[[i]]) <- c(expression(beta[0]), expression(beta[1]), expression(beta[2]))
  }
  
  if(pargroups == 1){
    
    ACF <- BetaAggregate(out, Methodname)
    
    ratio <- floor(n/40)
    
    ks <- c(ratio*(1:40))
    
    Gel <- GelmanWrapper(out, M, ks)
    
    Gel <- data.frame(t(Gel), colnames(Gel))
    
    colnames(Gel)[length(colnames(Gel))] <- "Iteration"
    
    mcmclist <- MCMCfunc(out, M)
    
    
    frac1 <- c(GewInit)
    frac2 <- c(GewEnd)
    
    gew <- gewekewrapper2(mcmclist, frac1, frac2, M, p)
    rownames(gew) <- paste("Par", 1:p)
    colnames(gew) <- paste("Chain", 1:M)
  }
  
  if(pargroups != 1){
    ACF <- list(1)
    GEl <- list(1)
    
    for(i in 1:pargroups){
      
      
      ACF[[i]] <- BetaAggregate(out[[i]], Methodname)
      
      ratio <- floor(n/40)
      
      ks <- c(ratio*(1:40))
      
      Gel[[i]] <- GelmanWrapper(out[[i]], M, ks)
      
      Gel[[i]] <- data.frame(t(Gel[[i]]), colnames(Gel[[i]]))
      
      colnames(Gel[[i]])[length(colnames(Gel[[i]]))] <- "Iteration"
      
      
      
      mcmcmlist[[i]] <- MCMCfunc(out[[i]], M)
      frac1 <- c(GewInit)
      frac2 <- c(GewEnd)
      
      gew[[i]] <- gewekewrapper2(mcmclist[[i]], frac1, frac2, M, p)
      rownames(gew[[i]]) <- paste("Par", 1:p)
      colnames(gew[[i]]) <- paste("Chain", 1:M)
    }
    
    
    
  }
  
  
  
  finalout <- list(Chains = out,ACF =  ACF,GR = Gel, Geweke = list(Geweke = gew, Fractions = c(GewInit, GewEnd)), Initial.Values = Initials)
  
  return(finalout)  
  
}

