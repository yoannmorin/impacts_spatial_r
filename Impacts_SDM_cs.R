####################################################################
####################################################################
######                                                         #####
######    Impact calculation in Spatial Durbin models in R     #####
######                                                         #####
####################################################################
####################################################################



# Install Packages

install.packages("MASS")
install.packages("maptools")
install.packages("sp")
install.packages("spdep")
install.packages("Ecdat")
install.packages("splm")



# Load packages 

library(MASS)
library(maptools)
library(sp)
library(spdep)
library(Ecdat)
library(splm)



# Working directory

setwd("G:/Dropbox/Boulot recherche_autres/Gitlab/Codes/impacts_spatial_r")






########################
##### Prepare Data #####
########################


# Load database from ecdat package

data("Produc", package = "Ecdat")
print(Produc)



# A binary contiguity spatial weights matrix for the US states

data("usaww")



# Contiguiity matrix to listw

usalw<-mat2listw(usaww, style="W")



# Cross section estimation on the mean of values in the database Produc

attach(Produc)
munnell_cross<-aggregate(Produc, by=list(state), FUN='mean')
munnell_cross$state<-munnell_cross$Group.1
munnell_cross <- subset(munnell_cross, select = -c(year,Group.1) )






######################
##### Estimation #####
######################


# Model to esitmate 
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 

durbin_var <- ~ log(pcap) + log(pc)



# Estimation of SDM using maximum likelihood
munnell_SAR_ML<-lagsarlm(fm, data=munnell_cross, listw=usalw, type = "lag", Durbin=durbin_var)
summary(munnell_SAR_ML)








#######################################
##### Impact Estimation : my code #####
#######################################



### Identification des coefficients ###
COEF1 <- munnell_SAR_ML$rho             ### rho
COEF2 <- munnell_SAR_ML$coefficients[2] ### log(pcap)
COEF3 <- munnell_SAR_ML$coefficients[3] ### log(pc)
COEF4 <- munnell_SAR_ML$coefficients[4] ### log(emp)
COEF5 <- munnell_SAR_ML$coefficients[5] ### unemp
COEF6 <- munnell_SAR_ML$coefficients[6] ### spatial lag log(pcap)
COEF7 <- munnell_SAR_ML$coefficients[7] ### spatial lag log(pc)


# Inverse (I - rho W)^-1
iIrW <- invIrW(usalw, COEF1)




### Construction the (dense) S(W)_r matrices ###

S_pcap <- iIrW %*% ((COEF2*diag(nrow(munnell_cross)))+ (COEF6*listw2mat(usalw)))
S_pc <- iIrW %*% ((COEF3*diag(nrow(munnell_cross)))+ (COEF7*listw2mat(usalw)))
S_emp <- iIrW %*% ((COEF4*diag(nrow(munnell_cross))))
S_unemp <- iIrW %*% ((COEF5*diag(nrow(munnell_cross))))


### Direct impacts ###
dir_pcap <- sum(diag(S_pcap))/nrow(munnell_cross)
dir_pc <- sum(diag(S_pc))/nrow(munnell_cross)
dir_emp  <- sum(diag(S_emp))/nrow(munnell_cross)
dir_unemp <- sum(diag(S_unemp))/nrow(munnell_cross)


### Total impacts ####

tot_pcap <- sum(c(S_pcap))/nrow(munnell_cross)
tot_pc <- sum(c(S_pc))/nrow(munnell_cross)
tot_emp <- sum(c(S_emp))/nrow(munnell_cross)
tot_unemp <- sum(c(S_unemp))/nrow(munnell_cross)


### Indirect impact

indir_pcap <- tot_pcap - dir_pcap
indir_pc <- tot_pc - dir_pc
indir_emp <- tot_emp - dir_emp
indir_unemp <- tot_unemp - dir_unemp





####################################################################
##### Simulations : dispersion measures and significance tests #####
####################################################################


R=10000 # number of simulations


select_beta <- munnell_SAR_ML$coefficients   # select coefficients
select_rho <- munnell_SAR_ML$rho  # select rho
select_sigma <- munnell_SAR_ML$s2   # select sigma
mu <- c(select_sigma, select_rho, select_beta)


Sigma <- munnell_SAR_ML$resvar # select variance-covariance matrix


tol=1e-10
empirical=FALSE
set.seed(12345) # fix seed


# Simulate values for coefficients

samplesA1 <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,empirical=empirical)

# Simulated values for rho have to be in the [-1;1[ range, we remove simulated observations with values for rho outside this range.

samplesA1 <- subset(samplesA1, samplesA1[,2] < 1 )
samplesA1 <- subset(samplesA1, samplesA1[,2] >= -1 )



# Initialize matrix impacts simulations will be stored in

simuA1_pcap <- matrix(nrow=nrow(samplesA1),ncol=3)  
simuA1_pc <- matrix(nrow=nrow(samplesA1),ncol=3)  
simuA1_emp <- matrix(nrow=nrow(samplesA1),ncol=3)  
simuA1_unemp <- matrix(nrow=nrow(samplesA1),ncol=3)  




# Compute impact measures for each simulated values

check <- seq(1,nrow(samplesA1),100) # sequence to print progress in the loop



for ( i in 1:nrow(samplesA1)) { # loop through simulated parameters values.
  
  if (i%in% check){print(i)}  # print progress in the loop
  
  iIrW <- invIrW(usalw, samplesA1[i,2]) # Inverse (I - rho W)^-1
  
  
  # Each bloc compute for one of the variable : the (dense) S(W)_r matrices, direct impact, total impact and indirect impact.
  
  simu_effectA1_pcap <- iIrW %*% ((samplesA1[i,4]*diag(nrow(munnell_cross)))+ (samplesA1[i,8]*listw2mat(usalw))) 
  simuA1_pcap[i,1] <- sum(diag(simu_effectA1_pcap))/nrow(munnell_cross) 
  simuA1_pcap[i,2] <- sum(c(simu_effectA1_pcap))/nrow(munnell_cross)
  simuA1_pcap[i,3] <- simuA1_pcap[i,2] - simuA1_pcap[i,1] 
  
  simu_effectA1_pc <- iIrW %*% ((samplesA1[i,5]*diag(nrow(munnell_cross)))+ (samplesA1[i,9]*listw2mat(usalw))) 
  simuA1_pc[i,1] <- sum(diag(simu_effectA1_pc))/nrow(munnell_cross) 
  simuA1_pc[i,2] <- sum(c(simu_effectA1_pc))/nrow(munnell_cross)
  simuA1_pc[i,3] <- simuA1_pc[i,2] - simuA1_pc[i,1] 
  
  simu_effectA1_emp <- iIrW %*% ((samplesA1[i,6]*diag(nrow(munnell_cross)))) 
  simuA1_emp[i,1] <- sum(diag(simu_effectA1_emp))/nrow(munnell_cross) 
  simuA1_emp[i,2] <- sum(c(simu_effectA1_emp))/nrow(munnell_cross)
  simuA1_emp[i,3] <- simuA1_emp[i,2] - simuA1_emp[i,1] 
  
  simu_effectA1_unemp <- iIrW %*% ((samplesA1[i,7]*diag(nrow(munnell_cross)))) 
  simuA1_unemp[i,1] <- sum(diag(simu_effectA1_unemp))/nrow(munnell_cross) 
  simuA1_unemp[i,2] <- sum(c(simu_effectA1_unemp))/nrow(munnell_cross)
  simuA1_unemp[i,3] <- simuA1_unemp[i,2] - simuA1_unemp[i,1] 
}




### Compute means, standard-deviations, z-stats, et p-values

# Compute stadanrd errors of the sample
sd_dir_pcap <- sd(simuA1_pcap[,1])
sd_tot_pcap <- sd(simuA1_pcap[,2])
sd_indir_pcap <- sd(simuA1_pcap[,3])

sd_dir_pc <- sd(simuA1_pc[,1])
sd_tot_pc <- sd(simuA1_pc[,2])
sd_indir_pc <- sd(simuA1_pc[,3])

sd_dir_emp <- sd(simuA1_emp[,1])
sd_tot_emp <- sd(simuA1_emp[,2])
sd_indir_emp <- sd(simuA1_emp[,3])

sd_dir_unemp <- sd(simuA1_unemp[,1])
sd_tot_unemp <- sd(simuA1_unemp[,2])
sd_indir_unemp <- sd(simuA1_unemp[,3])


# compute means of impact measures
mean_dir_pcap <- mean(simuA1_pcap[,1])
mean_tot_pcap <- mean(simuA1_pcap[,2])
mean_indir_pcap <- mean(simuA1_pcap[,3])

mean_dir_pc <- mean(simuA1_pc[,1])
mean_tot_pc <- mean(simuA1_pc[,2])
mean_indir_pc <- mean(simuA1_pc[,3])

mean_dir_emp <- mean(simuA1_emp[,1])
mean_tot_emp <- mean(simuA1_emp[,2])
mean_indir_emp <- mean(simuA1_emp[,3])

mean_dir_unemp <- mean(simuA1_unemp[,1])
mean_tot_unemp <- mean(simuA1_unemp[,2])
mean_indir_unemp <- mean(simuA1_unemp[,3])


# Compute z-stats
z_dir_pcap <- mean_dir_pcap / (sd_dir_pcap)
z_tot_pcap <- mean_tot_pcap / (sd_tot_pcap)
z_indir_pcap <- mean_indir_pcap / (sd_indir_pcap)

z_dir_pc <- mean_dir_pc / (sd_dir_pc)
z_tot_pc <- mean_tot_pc / (sd_tot_pc)
z_indir_pc <- mean_indir_pc / (sd_indir_pc)

z_dir_emp <- mean_dir_emp / (sd_dir_emp)
z_tot_emp <- mean_tot_emp / (sd_tot_emp)
z_indir_emp <- mean_indir_emp / (sd_indir_emp)

z_dir_unemp <- mean_dir_unemp / (sd_dir_unemp)
z_tot_unemp <- mean_tot_unemp / (sd_tot_unemp)
z_indir_unemp <- mean_indir_unemp / (sd_indir_unemp)


# compute p-values
pvalue_dir_pcap<-2*(1-pnorm(abs(z_dir_pcap),lower.tail=TRUE))
pvalue_tot_pcap<-2*(1-pnorm(abs(z_tot_pcap),lower.tail=TRUE))
pvalue_indir_pcap<-2*(1-pnorm(abs(z_indir_pcap),lower.tail=TRUE))

pvalue_dir_pc<-2*(1-pnorm(abs(z_dir_pc),lower.tail=TRUE))
pvalue_tot_pc<-2*(1-pnorm(abs(z_tot_pc),lower.tail=TRUE))
pvalue_indir_pc<-2*(1-pnorm(abs(z_indir_pc),lower.tail=TRUE))

pvalue_dir_emp<-2*(1-pnorm(abs(z_dir_emp),lower.tail=TRUE))
pvalue_tot_emp<-2*(1-pnorm(abs(z_tot_emp),lower.tail=TRUE))
pvalue_indir_emp<-2*(1-pnorm(abs(z_indir_emp),lower.tail=TRUE))

pvalue_dir_unemp<-2*(1-pnorm(abs(z_dir_unemp),lower.tail=TRUE))
pvalue_tot_unemp<-2*(1-pnorm(abs(z_tot_unemp),lower.tail=TRUE))
pvalue_indir_unemp<-2*(1-pnorm(abs(z_indir_unemp),lower.tail=TRUE))




##### Tables to resume results #####

# True values of impact estimates
impact_estimates =  matrix(  c(dir_pcap, indir_pcap, tot_pcap,    dir_pc, indir_pc, tot_pc,    dir_emp, indir_emp, tot_emp,    dir_unemp, indir_unemp, tot_unemp), nrow=4, ncol=3, byrow = TRUE) 

dimnames(impact_estimates) = list(c("pcap", "pc", "emp", "unemp"), c("direct", "indirect", "total"))


# Mean of impact measures estimated in the simulation
impact_mean = matrix(  c(mean_dir_pcap, mean_indir_pcap, mean_tot_pcap,    mean_dir_pc, mean_indir_pc, mean_tot_pc,    mean_dir_emp, mean_indir_emp, mean_tot_emp,    mean_dir_unemp, mean_indir_unemp, mean_tot_unemp), nrow=4, ncol=3, byrow = TRUE) 

dimnames(impact_mean) = list(c("pcap", "pc", "emp", "unemp"), c("direct", "indirect", "total"))


# Standard deviation of impact measures estimated in the simulation
impact_sd = matrix(  c(sd_dir_pcap, sd_indir_pcap, sd_tot_pcap,    sd_dir_pc, sd_indir_pc, sd_tot_pc,    sd_dir_emp, sd_indir_emp, sd_tot_emp,    sd_dir_unemp, sd_indir_unemp, sd_tot_unemp), nrow=4, ncol=3, byrow = TRUE) 

dimnames(impact_sd) = list(c("pcap", "pc", "emp", "unemp"), c("direct", "indirect", "total"))


# Z-stats of impact measures estimated in the simulation
impact_z = matrix(  c(z_dir_pcap, z_indir_pcap, z_tot_pcap,    z_dir_pc, z_indir_pc, z_tot_pc,    z_dir_emp, z_indir_emp, z_tot_emp,    z_dir_unemp, z_indir_unemp, z_tot_unemp), nrow=4, ncol=3, byrow = TRUE) 

dimnames(impact_z) = list(c("pcap", "pc", "emp", "unemp"), c("direct", "indirect", "total"))


# P-values of impact measures estimated in the simulation
impact_pvalue = matrix(  c(pvalue_dir_pcap, pvalue_indir_pcap, pvalue_tot_pcap,    pvalue_dir_pc, pvalue_indir_pc, pvalue_tot_pc,    pvalue_dir_emp, pvalue_indir_emp, pvalue_tot_emp,    pvalue_dir_unemp, pvalue_indir_unemp, pvalue_tot_unemp), nrow=4, ncol=3, byrow = TRUE) 

dimnames(impact_pvalue) = list(c("pcap", "pc", "emp", "unemp"), c("direct", "indirect", "total"))




# Print tables
impact_estimates
impact_mean
impact_sd
impact_z
impact_pvalue
























