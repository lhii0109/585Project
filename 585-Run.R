

Clusters <- detectCores()-1
cl <- makePSOCKcluster(Clusters)
registerDoParallel(cl)

source("background-functions.R")
  
source("585-Diagnostic.R")
  
Result1 <- Diagnostic(M = 3, K = 80, x = X, y = y, Initialize = 1, Model = "SandProbit", Methodname = "SandwichProbit")
Result2 <- Diagnostic(M = 3, K = 80, x = X, y = y, Initialize = 1, Model = "ASISProbit", Methodname = "ASISProbit")

Result1
Result2

ACFPlot(Result1, Result2)
GelmanPlot(Result1, Result2, parnames = C("Par[0]", "Par[1]", "Par[2]"), name1 = "SandwichProbit", name2 = "ASISProbit", p = 3)

stopCluster(cl)

