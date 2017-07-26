
mcmcR2 <- function(mod, type = "marginal", family = "gaussian"){
    if (type != "marginal") stop("At the moment, there is just the marginal R2")
    if (family != "gaussian") stop("At the moment just gaussian errors are supported")
    # Shinichis answer on Researchgate
    
    mVarF <- var(as.vector(apply(mod$Sol,2,mean) %*% t(mod$X)))
    R2 <- mVarF/(mVarF+sum(apply(mod$VCV,2,mean)))
    
    # alternative with crebile intervals
    mcmc_chain_length <- nrow(mod$VCV)
    vmVarF<-numeric( mcmc_chain_length)
    
    for(i in 1:mcmc_chain_length){
        Var <- var(as.vector(mod$Sol[i,] %*% t(mod$X)))
        vmVarF[i]<-Var
    }
    
    R2m<-vmVarF/(vmVarF+ rowSums(mod$VCV)) # include here all random effects plus errors
    class(R2m) <- "mcmc"
    R2m<-vmVarF/(vmVarF+mod$VCV[,1]+mod$VCV[,2])
    
    mean(R2m)
    posterior.mode(R2m)
    data.frame(HPDinterval(R2m), row.names = NULL)
    
    out <- data.frame("meanR2" = mean(R2m), "modeR2" =  posterior.mode(R2m), data.frame(HPDinterval(R2m)),  row.names = NULL)
    
}



partR2 <- function(mod, partvars = NULL){
    
    if (is.null(partvars)) stop("partvars has to contain the variables for the commonality analysis")
        
    model_formula <- mod$Fixed$formula
    
    # just unique effects
    all_unique_R2 <- c()
    R2_full <- mcmcR2(mod)
    for (i in partvars){
        new_formula <- update.formula(model_formula, paste(". ~ . -", i, sep=""))

        mod_red <- MCMCglmm(new_formula,
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=test ,nitt=100000,burnin=1000,thin=20)

        R2_red <- mcmcR2(mod_red)
        all_unique_R2[i] <- R2_full$modeR2 - R2_red$modeR2
    }
    
    names(all_unique_R2) <- partvars
    all_unique_R2
}

out2 <- partR2(mod, partvars = partvars)



mod0 <- lm(TPM70_ratio ~ obs_het_mean + num_alleles_mean + mean_allele_range + prop_low_afs_mean, data = stats_mod) #prop_low_afs_mean + 
summary(mod0)
calc.yhat(mod0)


mod1 <- lm(TPM70_ratio ~ obs_het_mean + num_alleles_mean + mean_allele_range , data = stats_mod)

sum_mod0 <- summary(mod0)
sum_mod1 <- summary(mod1)

sum_mod0$r.squared - sum_mod1$r.squared
