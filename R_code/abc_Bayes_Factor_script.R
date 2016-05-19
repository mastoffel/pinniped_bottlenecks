####Based on a script provided by Vitor Sousa (Sousa et al 2010, Animal Conservation). Modified by Erwan Quéméré and Paz-Vinas Ivan###

library(locfit)

	# Reading the files
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


#################
# REMOVE BURNIN #
#################
	#par(mfrow=c(4,1))
	#plot(run04[,1])
	##plot(run04[,2])
	#plot(run04[,3])
	#plot(run04[,4])
	
	dim(run01)
	dim(run02)
	dim(run03)
	dim(run04)

	run01 <- run01[-(1:(length(run01[,1])/10)),]
	run02 <- run02[-(1:(length(run02[,1])/10)),]
	run03 <- run03[-(1:(length(run03[,1])/10)),]
	run04 <- run04[-(1:(length(run04[,1])/10)),]

	dim(run01)
	dim(run02)
	dim(run03)
	dim(run04)


	# N0
	run01N0 <- 10^run01[,1]
	run02N0 <- 10^run02[,1]
	run03N0 <- 10^run03[,1]
	run04N0 <- 10^run04[,1]

	# N1
	run01N1 <- 10^run01[,2]
	run02N1 <- 10^run02[,2]
	run03N1 <- 10^run03[,2]
	run04N1 <- 10^run04[,2]

	# U
	run01U <- 10^run01[,3]
	run02U <- 10^run02[,3]
	run03U <- 10^run03[,3]
	run04U <- 10^run04[,3]


	# TF
	run01TF <- 10^run01[,4]
	run02TF <- 10^run02[,4]
	run03TF <- 10^run03[,4]
	run04TF <- 10^run04[,4]
	
################
# CLEAN MEMORY #
################
rm(run01, run02, run03, run04)

##########################
# PUTS THE RUNS TOGETHER #
##########################

	runN0 <- c( run01N0, run02N0, run03N0, run04N0, run05N0)
	runN1 <- c(run01N1, run02N1, run03N1, run04N1, run05N1)
	runU <- c(run01U, run02U, run03U, run04U, run05U)
	runTF <- c(run01TF, run02TF, run03TF, run04TF, run05TF)
	
	length(runN0)
	length(runN1)
	length(runU)
	length(runTF)
	
	#quantilesNatural <- array(dim=c(4,7))

	#quantilesNatural[1,] <- quantile(runN0, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99) )
	#quantilesNatural[2,] <- quantile(runN1, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99) )
	#quantilesNatural[3,] <- quantile(runU, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99) )
	#quantilesNatural[4,] <- quantile(runTF, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99) )

	#round(quantilesNatural)
	
	#write.table(quantilesNatural, "quantilesNaturalScaleAllTog.txt", dec=",")

####################################################
# "APPROXIMATE" BAYES FACTOR FOR THE TIME TF       #
####################################################


##########################################################################################################################################


	bayes.factor.interval <- function(interval, prior, posterior) {

		hist_breaks = c(0,interval, 1e35)

		post <- hist(posterior, plot=F, breaks=hist_breaks)
		prior <- hist(prior, plot=F, breaks=hist_breaks)

		cpost <- length(post$count)
		cprior <- length(prior$count)

		if(cpost == cprior) {

			sum_post <- sum(post$count[1:3])
			sum_prior <- sum(prior$count[1:3])

			bf <- ( post$count[2]/(sum_post - post$count[2]) ) / (
prior$count[2]/(sum_prior - prior$count[2]) )

		}
		else {
			bf == "NA"
		}

		bf

	}

####################################################################################################################################################


norm5.3 <- 10^rnorm(200000, mean=5, sd=3)

vec = c(seq(from=0,to=100, by=10),seq(from=200, to=10000, by=100))

bfEastClusterMaxgroup.R1 <- vector(length=0)
for (i in 1:(length(vec)-1)) {
	int <- c(vec[i], vec[i+1])
	bfEastClusterMaxgroup.R1[i] <- bayes.factor.interval( int, norm5.3, runTF)
}

bfEastClusterMaxgroup.R1

vec[bfEastClusterMaxgroup.R1==max(bfEastClusterMaxgroup.R1)]
max(bfEastClusterMaxgroup.R1)

pdf(file="Viaur_Phox_Full_BF.pdf")

plot(vec[1:(length(vec)-1)], bfEastClusterMaxgroup.R1,xlim=c(0,10000),ylim=c(0,13),main="Viaur_Phox_Tout BF",xlab="time", ylab="Bayes Factor",type="l", col="green")
abline(v=100, lty=2, col="dark grey")
abline(v=500, lty=2, col="dark grey")
abline(v=1000, lty=2, col="dark grey")
abline(v=2000, lty=2, col="dark grey")
abline(v=3000, lty=2, col="dark grey")
abline(v=4000, lty=2, col="dark grey")
abline(v=5000, lty=2, col="dark grey")
abline(v=6000, lty=2, col="dark grey")
abline(v=7000, lty=2, col="dark grey")
abline(v=8000, lty=2, col="dark grey")
abline(v=9000, lty=2, col="dark grey")
abline(v=10000, lty=2, col="dark grey")

abline(v=vec[bfEastClusterMaxgroup.R1==max(bfEastClusterMaxgroup.R1)])

abline(h=4, lty=2, col="dark grey")
abline(h=7, lty=2, col="dark grey")
abline(h=8, lty=2, col="dark grey")
abline(h=9, lty=2, col="dark grey")
abline(h=10, lty=2, col="dark grey")

vec[bfEastClusterMaxgroup.R1==max(bfEastClusterMaxgroup.R1)]
max(bfEastClusterMaxgroup.R1)

dev.off()

####################################################################################################################################################
