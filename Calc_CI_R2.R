# Shinichi calculates R2 CI 
####################################################

# A. Preparation

####################################################

# Note that data generation appears below the analysis section.

# You can use the simulated data table from the supplementary files to reproduce exactly the same results as presented in the paper.
# Set the work directy that is used for rading/saving data tables

# setwd("/Users/R2")


# load R required packages

# If this is done for the first time, it might need to first download and install the package

# install.packages("arm")

library(arm)

# install.packages("lme4")

# the verson 1.0-5

library(lme4)

# install.packages("MCMCglmm")

library(MCMCglmm)
####################################################

# B. Analysis

####################################################
# 1. Analysis of body size (Gaussian mixed models)

#---------------------------------------------------
# Clear memory

rm(list = ls())
# Read body length data (Gaussian, available for both sexes)

Data <- read.csv("BeetlesBody.csv")
# Fit null model without fixed effects (but including all random effects)

m0 <- lmer(BodyL ~ 1 + (1 | Population) + (1 | Container), data = Data)
# MCMCglmm model
mm0<-MCMCglmm(BodyL ~ 1 , random = ~ Population + Container, data=Data)
# Fit alternative model including fixed and all random effects

mF <- lmer(BodyL ~ Sex + Treatment + Habitat + (1 | Population) + (1 | Container), data = Data)
# MCMCglmm model
mmF<-MCMCglmm(BodyL ~ Sex + Treatment + Habitat , random = ~ Population + Container, data=Data)
# View model fits for both models

summary(m0)

summary(mm0)
summary(mF)

summary(mmF)
# Extraction of fitted value for the alternative model

# fixef() extracts coefficents for fixed effects

# mF@pp$X returns fixed effect design matrix

Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3] + fixef(mF)[4] * mF@pp$X[, 4]


# MCMCglmm (it is probably better to get a posterior distribuiton of R2 rather than getting each varaince component - we do this below as an alternative)

mFixed <- mean(mmF$Sol[,2]) * mmF$X[, 2] + mean(mmF$Sol[, 3]) * mmF$X[, 3] + mean(mmF$Sol[ ,4]) * mmF$X[, 4]


# Calculation of the variance in fitted values

VarF <- var(Fixed)

mVarF<- var(mFixed)


# An alternative way for getting the same result

VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))

mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))
# R2GLMM(m) - marginal R2GLMM

# Equ. 26, 29 and 30

# VarCorr() extracts variance components

# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance

VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + attr(VarCorr(mF), "sc")^2)


# MCMCglmm - marginal
mVarF/(mVarF+sum(apply(mmF$VCV,2,mean)))


# alternative with crebile intervals
vmVarF<-numeric(1000)

for(i in 1:1000){
    
    Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
    
    vmVarF[i]<-Var
    }


R2m<-vmVarF/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])

mean(R2m)

posterior.mode(R2m)
HPDinterval(R2m)

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30

(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + (attr(VarCorr(mF), "sc")^2))


# MCMCglmm - conditional

(mVarF+sum(apply(mmF$VCV,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals

R2c<-(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])

mean(R2c)

posterior.mode(R2c)

HPDinterval(R2c) 

#How can I calculate R2 for an Bayesian MCMC multilevel model?. Available from: https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model [accessed Jul 24, 2017].