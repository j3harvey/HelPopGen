dat00 <- read.csv(file="ag5am5era4_allGenes_ambig_0.0.csv",head=TRUE,sep=",")
dat04 <- read.csv(file="ag5am5era4_allGenes_ambig_0.04.csv",head=TRUE,sep=",")
dat08 <- read.csv(file="ag5am5era4_allGenes_ambig_0.08.csv",head=TRUE,sep=",")
dat12 <- read.csv(file="ag5am5era4_allGenes_ambig_0.12.csv",head=TRUE,sep=",")
dat16 <- read.csv(file="ag5am5era4_allGenes_ambig_0.16.csv",head=TRUE,sep=",")
dat20 <- read.csv(file="ag5am5era4_allGenes_ambig_0.20.csv",head=TRUE,sep=",")
dat24 <- read.csv(file="ag5am5era4_allGenes_ambig_0.24.csv",head=TRUE,sep=",")
dat28 <- read.csv(file="ag5am5era4_allGenes_ambig_0.28.csv",head=TRUE,sep=",")
dat32 <- read.csv(file="ag5am5era4_allGenes_ambig_0.32.csv",head=TRUE,sep=",")
dat36 <- read.csv(file="ag5am5era4_allGenes_ambig_0.36.csv",head=TRUE,sep=",")
dat40 <- read.csv(file="ag5am5era4_allGenes_ambig_0.4.csv",head=TRUE,sep=",")

alphaBar <- function(data){
  P_N = data$P_N
  P_S = data$P_S
  D_N = data$D_N
  D_S = data$D_S
  return( 1 - mean(D_N)*mean( P_S / (P_N + 1) )/mean(D_S) )
}

alphaBar2 <- function(data){
  P_N = data$P_N
  P_S = data$P_S
  D_N = data$D_N
  D_S = data$D_S
  return( 1 - mean(D_N)*mean(P_S)/(mean(P_N)*mean(D_S)))
}


