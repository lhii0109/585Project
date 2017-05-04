source("background-functions.R")

source("585-Diagnostic.R")

Clusters <- detectCores()-1
cl <- makePSOCKcluster(Clusters)
registerDoParallel(cl)

Diagnostic()

ACFPlot()

GelmanPlot()

stopCluster(cl)

