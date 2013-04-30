dat00 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_00.csv",head=TRUE,sep=",")
dat02 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_02.csv",head=TRUE,sep=",")
dat04 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_04.csv",head=TRUE,sep=",")
dat06 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_06.csv",head=TRUE,sep=",")
dat08 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_08.csv",head=TRUE,sep=",")
dat10 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_10.csv",head=TRUE,sep=",")
dat12 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_12.csv",head=TRUE,sep=",")
dat16 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_16.csv",head=TRUE,sep=",")
dat20 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_20.csv",head=TRUE,sep=",")
dat24 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_24.csv",head=TRUE,sep=",")
dat28 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_28.csv",head=TRUE,sep=",")
dat32 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_32.csv",head=TRUE,sep=",")
dat36 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_36.csv",head=TRUE,sep=",")
dat40 <- read.csv(file="ag5am5era4_allGenes_ambig_agam_40.csv",head=TRUE,sep=",")

date00 <- read.csv(file="ag5am5era4_allGenes_ambig_0.0.csv",head=TRUE,sep=",")
date04 <- read.csv(file="ag5am5era4_allGenes_ambig_0.04.csv",head=TRUE,sep=",")
date08 <- read.csv(file="ag5am5era4_allGenes_ambig_0.08.csv",head=TRUE,sep=",")
date12 <- read.csv(file="ag5am5era4_allGenes_ambig_0.12.csv",head=TRUE,sep=",")
date16 <- read.csv(file="ag5am5era4_allGenes_ambig_0.16.csv",head=TRUE,sep=",")
date20 <- read.csv(file="ag5am5era4_allGenes_ambig_0.2.csv",head=TRUE,sep=",")
date24 <- read.csv(file="ag5am5era4_allGenes_ambig_0.24.csv",head=TRUE,sep=",")
date28 <- read.csv(file="ag5am5era4_allGenes_ambig_0.28.csv",head=TRUE,sep=",")
date32 <- read.csv(file="ag5am5era4_allGenes_ambig_0.32.csv",head=TRUE,sep=",")
date36 <- read.csv(file="ag5am5era4_allGenes_ambig_0.36.csv",head=TRUE,sep=",")
date40 <- read.csv(file="ag5am5era4_allGenes_ambig_0.4.csv",head=TRUE,sep=",")

date00 <- read.csv(file="third_run_00.csv",head=TRUE,sep=",")
date02 <- read.csv(file="third_run_02.csv",head=TRUE,sep=",")
date04 <- read.csv(file="third_run_04.csv",head=TRUE,sep=",")
date06 <- read.csv(file="third_run_06.csv",head=TRUE,sep=",")
date08 <- read.csv(file="third_run_08.csv",head=TRUE,sep=",")
date10 <- read.csv(file="third_run_10.csv",head=TRUE,sep=",")
date12 <- read.csv(file="third_run_12.csv",head=TRUE,sep=",")
date16 <- read.csv(file="third_run_16.csv",head=TRUE,sep=",")
date20 <- read.csv(file="third_run_20.csv",head=TRUE,sep=",")
date24 <- read.csv(file="third_run_24.csv",head=TRUE,sep=",")
date28 <- read.csv(file="third_run_28.csv",head=TRUE,sep=",")
date32 <- read.csv(file="third_run_32.csv",head=TRUE,sep=",")
date36 <- read.csv(file="third_run_36.csv",head=TRUE,sep=",")

datAm00 <- read.csv(file="amIn_0.0.csv", head=TRUE, sep=",")
datAm01 <- read.csv(file="amIn_0.1.csv", head=TRUE, sep=",")
datAm02 <- read.csv(file="amIn_0.2.csv", head=TRUE, sep=",")
datAm03 <- read.csv(file="amIn_0.3.csv", head=TRUE, sep=",")

datAg00 <- read.csv(file="agIn_0.0.csv", head=TRUE, sep=",")
datAg01 <- read.csv(file="agIn_0.1.csv", head=TRUE, sep=",")
datAg02 <- read.csv(file="agIn_0.2.csv", head=TRUE, sep=",")
datAg03 <- read.csv(file="agIn_0.3.csv", head=TRUE, sep=",")

alphaBar <- function(data){
  P_N = data$P_N
  P_S = data$P_S
  D_N = data$D_N
  D_S = data$D_S
  sl = data$sequence_length
  return( 1 - weighted.mean(D_S,sl)*weighted.mean(P_N,sl)/(weighted.mean(P_S,sl)*weighted.mean(D_N,sl)))
}

alphaBar2 <- function(data){
  P_N = data$P_N
  P_S = data$P_S
  D_N = data$D_N
  D_S = data$D_S
  sl = data$sequence_length
  return( 1 - weighted.mean(D_S,sl)*weighted.mean(P_N/(P_S + 1),sl)/weighted.mean(D_N,sl))
}



