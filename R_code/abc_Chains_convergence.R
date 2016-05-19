####Based on a script provided by Vitor Sousa (Sousa et al 2010, Animal Conservation). Modified by Erwan Quéméré and Paz-Vinas Ivan###
 
library(boa)
 
 #Reading the files
	hpars.run01 <- matrix(scan("./run01/hpars.dat"), byrow=T, ncol=11)

	# Just save 
	# LogN0 (present day pop size)
	# LogN1 (past pop size)
	# LogU (average mut rate)
	# LogTF (time since pop change)
	
	run01 <- hpars.run01[,c(4,6,8,10)]

	# Save the variance among loci of LogN0, LogN1, LogU and LogTF
	varRun01 <- hpars.run01[,c(5,7,9,11)]

	# Remove the hpars table - CLEAN MEMORY
	rm(hpars.run01)

	# The same procedure as for run02
	hpars.run02 <- matrix(scan("./run02/hpars.dat"), byrow=T, ncol=11)
	run02 <- hpars.run02[,c(4,6,8,10)]
	varRun02 <- hpars.run02[,c(5,7,9,11)]
	rm(hpars.run02)

	# The same procedure as for run03
	hpars.run03 <- matrix(scan("./run03/hpars.dat"), byrow=T, ncol=11)
	run03 <- hpars.run03[,c(4,6,8,10)]
	varRun03 <- hpars.run03[,c(5,7,9,11)]
	rm(hpars.run03)

	# The same procedure as for run04
	hpars.run04 <- matrix(scan("./run04/hpars.dat"), byrow=T, ncol=11)
	run04 <- hpars.run04[,c(4,6,8,10)]
	varrun04 <- hpars.run04[,c(5,7,9,11)]
	rm(hpars.run04)


# Creating MCMC objects
        library(coda)

        # Length of runs
        l1 <- length(run01[,1])
        l2 <- length(run02[,1])
        l3 <- length(run03[,1])
        l4 <- length(run04[,1])


        # Get minimum length
        minL <- min(l1, l2, l3, l4)

        # Create mcmc objects
        mcmcRun1 <- mcmc(run01[(length(run01[,1])-minL+1):(length(run01[,1])),])
        mcmcRun2 <- mcmc(run02[(length(run02[,1])-minL+1):(length(run02[,1])),])

        mcmcRun3 <- mcmc(run03[(length(run03[,1])-minL+1):(length(run03[,1])),])
        mcmcRun4 <- mcmc(run04[(length(run04[,1])-minL+1):(length(run04[,1])),])

        # Check the length of chains (they all must be the same)

        length(mcmcRun1[,1])
        length(mcmcRun2[,1])
        length(mcmcRun3[,1])
        length(mcmcRun4[,1])

        # Create a list of mcmc runs
        mcmcList <- mcmc.list(mcmcRun1, mcmcRun2, mcmcRun3, mcmcRun4)


        # Gelman diagnosis
        result <- gelman.diag(mcmcList)
        result
        write.table(data.frame(result[[1]]), "GelmanDiag.txt", dec=",",
row.names=FALSE, col.names=TRUE)

 #gelman.plot(mcmcList)


        # CLEAN MEMORY
        rm( mcmcRun1, mcmcRun2, mcmcRun3, result, minL, mcmcList, l1, l2, l3)