## modified function for the strataG package. Allows some of the (non-maximal sized and non-minimal )

mRatio_new <- function (g, by.strata = TRUE, rpt.size = 8:2, prop_corr_rep_length = 0.5) 
{
    calc.mratio <- function(freqs) {
        freqs <- freqs[, 1]
        if (length(freqs) == 1) {
            warning("only one allele")
            return(NA)
        }
        # if (!all.is.numeric(names(freqs))) {
        #     warning("allele names are non-numeric")
        #     return(NA)
        # }
        if (all(freqs == 0)) {
            warning("all frequencies are 0")
            NA
        }
        else {
            freqs <- freqs[order(as.numeric(names(freqs)))]
            sizes <- as.numeric(names(freqs))
            size.diff <- diff(sizes)
            rpt.found <- FALSE
            for (r in sort(rpt.size, decreasing = TRUE)) {
                # allow some loci not to have the correct repeats in between
                
                # check if smallest and largest allele can be divided by the repeat length
                if (all(size.diff[c(1, length(size.diff))] %% r == 0)) {
                    if ((sum(size.diff %% r == 0) / length(size.diff)) > prop_corr_rep_length) {
                    rpt.found <- TRUE
                    break
                    }
                }
                
            }
            if (!rpt.found) {
                warning("valid repeat length not found")
                return(NA)
            }
            smallest <- sizes[which.min(sizes[freqs > 0])]
            largest <- sizes[which.max(sizes[freqs > 0])]
            n <- (largest - smallest)/r
            sum(freqs > 0)/(n + 1)
        }
    }
    if (nStrata(g) == 1 & by.strata) {
        by.strata <- FALSE
        g <- g[, , strataNames(g)]
    }
    if (by.strata) {
        freqs <- alleleFreqs(g, by.strata = TRUE)
        do.call(rbind, lapply(freqs, function(loc) apply(loc, 
            3, calc.mratio)))
    }
    else {
        freqs <- alleleFreqs(g, by.strata = FALSE)
        sapply(freqs, calc.mratio)
    }
}