# 585Project

This is where I am going to store my 585 project! :)


Instructions for use: 

-read in "Backgroundfunctions.R"

-read in "585Diagnostic.R"

-Use 585Run.R to source in the previous functions, create clusters

-Use Diagnostic function to create output
    Inputs. 
    -Initials are the initial values for the chains, if provided. Can be created if you choose Initialize = 1.
    -MCMC is the option to input your own algorithm. It is important that M is referred to as the number of chains desired, and that the output is a list where each components is the chain from the Mth run, and that you pass into export a vector of function names used within the algorithm.
    -K is the length of chain that you want (must be bigger than 40)
    -x is the covariates
    -y is the responses
    -Model is one of four options: "SandRobit", "ASISRobit", "SandProbit" and "ASISProbit" for the methods used for comparison, reproducibility, and running.
    -Nu, nunot, and c are all values set for the Robit model Algorithms that are default options.
    -pargroups. You may wish to have more than one parameter group, eg, beta parameters and variance parameters. You can choose as many groups as you like. Note that this may fail and have some bugs in it still as it is hard to test for all cases. Output for differnet parameter groups must be listed as sublists. I.e, the main list component is the Mth chain, the first component of which is the first parameter group, the second the second parameter group, etc.
    -Methodname is the name of your algorithm. Must be put in quotes.
    -GewInit and GewEnd are the proportions of the chain to use at the beginning and end for the Geweke diagnostic values.

-Use ACFPlot to create ACF Plots to compare two methods Chains and Chains2 are the chains of interest.

-Use Gelman Plot to create Gelman Diagnostic Plots where Chains and Chains2 are the Chains of interest.

-Look at Results$Gewek (output for Diagnostic) for Geweke Results




'Vignette':

load("X.R")
load("y.R")

p <- dim(X)[2]

source("Backgroundfunctions.R")
source("585Diagnostic.R")


Results <- Diagnostic(M = 3, K = 80, x = X, y =  y, Model = "ASISProbit", Methodname = "test", Initialize = 1, nu = 9, nunot = 9, GewInit = .1, GewEnd = .5)

Results2 <- Diagnostic(M = 3, K = 80, x = X, y =  y, Model = "SandProbit", Methodname = "test2", Initialize = 1, nu = 9, nunot = 9, GewInit = .1, GewEnd = .5)


  
ACFPlot(Results, Results2)

GelmanPlot(Results, Results2, pargroup = 1, parnames = C("A", "B", "C"), name1 = "test1", name2 = "test2", p)
