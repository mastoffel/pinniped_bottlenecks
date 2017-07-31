
mod <- mod1


# function to calculate the R2 for gaussian MCMCglmm models
mcmcR2 <- function(mod, type = "marginal", family = "gaussian"){
    if (type != "marginal") stop("At the moment, there is just the marginal R2")
    if (family != "gaussian") stop("At the moment just gaussian errors are supported")
    # Shinichis answer on Researchgate
    mVarF <- var(as.vector(apply(mod$Sol,2,mean) %*% t(mod$X)))
    R2 <- mVarF/(mVarF+sum(apply(mod$VCV,2,mean)))
    
    # alternative with crebile intervals
    mcmc_chain_length <- nrow(mod$VCV)
    vmVarF<-numeric(mcmc_chain_length)
    
    vmVarF <- vapply(1:mcmc_chain_length, function(x) out <- var(as.vector(mod$Sol[x,] %*% t(mod$X))), 
                    FUN.VALUE = numeric(length(mcmc_chain_length)))
    
    R2m<-vmVarF/(vmVarF+ rowSums(mod$VCV)) # include here all random effects plus errors
    outR2m <- R2m
    class(R2m) <- "mcmc"
    R2m<-vmVarF/(vmVarF+mod$VCV[,1]+mod$VCV[,2])
    
    #mean(R2m)
    #posterior.mode(R2m)
    #data.frame(HPDinterval(R2m), row.names = NULL)
    
    out <- list("partR2" = data.frame("meanR2" = mean(R2m), "modeR2" =  posterior.mode(R2m), data.frame(HPDinterval(R2m)), row.names = NULL),
                "R2_chain" = outR2m)
}

out <- mcmcR2(mod)
out

# partition the R2 into variation unique and common to the predictors
partR2 <- function(mod, partvars = NULL, data = NULL, inv_phylo = NULL, prior = NULL, nitt=10000,burnin=1000, thin=50){
    
    if (is.null(partvars)) stop("partvars has to contain the variables for the commonality analysis")
        
    model_formula <- mod$Fixed$formula
    
    # just unique effects
    all_unique_R2 <- c()
    R2_full <- mcmcR2(mod)
    
    all_comb <- lapply(1:(length(partvars)-1), function(x) combn(partvars, x))
    all_comb2 <- lapply(all_comb, function(x) apply(x, 2, function(x) out <- list(x)))
    all_comb3 <- unlist(unlist(all_comb2, recursive = FALSE), recursive = FALSE)
    
    CI <- 0.95
    calc_CI <- function(x) {
        out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
    }
    
    calc_partR2 <- function(var_to_red, R2_full, data) {
        to_del <- paste(paste("-", var_to_red, sep= ""), collapse = " ")
        new_formula <- update.formula(model_formula, paste(". ~ . ", to_del, sep=""))
        
        mod_red <- MCMCglmm(new_formula,
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=data,nitt=nitt,burnin=burning,thin=thin)
        
        R2_red <- mcmcR2(mod_red)
        # R2 mode
        R2_diff <- R2_full$partR2$modeR2 - R2_red$partR2$modeR2
        # R2 CI
        R2_diff_vec <- calc_CI(R2_full$R2_chain - R2_red$R2_chain) 
        out <- c("R2mode" = R2_diff, "CIlow" = R2_diff_vec[1], "CIhigh" = R2_diff_vec[2])
    }

    R2_out <- do.call(rbind, lapply(all_comb3, calc_partR2, R2_full, data))
    # all_vars <- lapply(all_comb3, function(x) gsub('(a|e|i|o|u)', '', x))
    all_comb_names <- unlist(lapply(all_comb3, function(x) paste(x, collapse = " & ")))
    
    out <- data.frame("combinations" = all_comb_names, R2_out)
    out
}

out2 <- partR2(mod, partvars = partvars)



gsub('(a|e|i|o|u)', '', x = "BreedingType")




mod0 <- lm(TPM70_ratio ~ obs_het_mean + num_alleles_mean + mean_allele_range + prop_low_afs_mean, data = stats_mod) #prop_low_afs_mean + 
summary(mod0)
calc.yhat(mod0)


mod1 <- lm(TPM70_ratio ~ obs_het_mean + num_alleles_mean + mean_allele_range , data = stats_mod)

sum_mod0 <- summary(mod0)
sum_mod1 <- summary(mod1)

sum_mod0$r.squared - sum_mod1$r.squared
