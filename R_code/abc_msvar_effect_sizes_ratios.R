####Based on a script provided by Vitor Sousa (Sousa et al 2010, Animal Conservation). Modified by Erwan Quéméré and Paz-Vinas Ivan###


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
	#plot(run01[,1])
	#plot(run01[,2])
	#plot(run01[,3])
	#plot(run01[,4])

	#par(mfrow=c(4,1))
	#plot(run02[,1])
	#plot(run02[,2])
	#plot(run02[,3])
	#plot(run02[,4])

	#par(mfrow=c(4,1))
	#plot(run03[,1])
	#plot(run03[,2])
	#plot(run03[,3])
	#plot(run03[,4])

	#par(mfrow=c(4,1))
	#plot(run04[,1])
	#plot(run04[,2])
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
	run01N0 <- run01[,1]
	run02N0 <- run02[,1]
	run03N0 <- run03[,1]
	run04N0 <- run04[,1]

	# N1
	run01N1 <- run01[,2]
	run02N1 <- run02[,2]
	run03N1 <- run03[,2]
	run04N1 <- run04[,2]

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

	runN0 <- c(run01N0, run02N0, run03N0, run04N0)
	runN1 <- c(run01N1, run02N1, run03N1, run04N1)
	runTF <- c(run01TF, run02TF, run03TF, run04TF)

######################################
# MEDIAN N0 and N1 (for each run)#####
######################################	

10^median(run01N0)
10^median(run02N0)
10^median(run03N0)
10^median(run04N0)

10^median(run01N1)
10^median(run02N1)
10^median(run03N1)
10^median(run04N1)

#############################
# MEDIAN TF1(for each run)###
#############################	

median(run01TF)
median(run02TF)
median(run03TF)
median(run04TF)

###############################
#########Effect Sizes##########
###############################

Effects=rep(0,4)
I95Inf=rep(0,4)
I95Sup= rep(0,4)

for(i in 1:4){

TAILLENE=c(length(run01N0),length(run02N0),length(run03N0),length(run04N0))
TAILLENC=c(length(run01N1),length(run02N1),length(run03N1),length(run04N1))
VALSE=c(sd(run01N0),sd(run02N0),sd(run03N0),sd(run04N0))
VALXE=c(mean(run01N0),mean(run02N0),mean(run03N0),mean(run04N0))
VALSC=c(sd(run01N1),sd(run02N1),sd(run03N1),sd(run04N1))
VALXC=c(mean(run01N1),mean(run02N1),mean(run03N1),mean(run04N1))

NE=TAILLENE[i]
SE=VALSE[i]
XE=VALXE[i]
NC=TAILLENE[i]
SC=VALSC[i]
XC=VALXC[i]

J=1-(3/(4*(NC+NE-2)-1))
S=sqrt(((NE-1)*(SE)^2+(NC-1)*(SC)^2)/(NE+NC-2))
G=((XE-XC)/S)*J
SE=sqrt(((NC+NE)/(NC+NE))+((G^2)/(2*(NC+NE-2))))

Effects[i]=G
I95Inf[i]=G-1.96*SE
I95Sup[i]=G+1.96*SE

}

Effects
I95Inf
I95Sup
TAILLENE
TAILLENC
VALSE
VALXE
VALSC
VALXC

lesrun0=c(10^median(run01N0),10^median(run02N0),10^median(run03N0),10^median(run04N0))
lesrun1=c(10^median(run01N1),10^median(run02N1),10^median(run03N1),10^median(run04N1))
lesrunt=c(median(run01TF),median(run02TF),median(run03TF),median(run04TF))

###############################
###########Ratios##############
###############################

rat1=na.omit((log(run01N0/run01N1)),T)
rat1=rat1[is.finite(rat1)]

rat2=na.omit((log(run02N0/run02N1)),T)
rat2=rat2[is.finite(rat2)]

rat3=na.omit((log(run03N0/run03N1)),T)
rat3=rat3[is.finite(rat3)]

rat4=na.omit((log(run04N0/run04N1)),T)
rat4=rat4[is.finite(rat4)]

e=c(mean(rat1),mean(rat2),mean(rat3),mean(rat4))
mr=mean(e)
sdmr=sd(e)


ratciINF=mr+(qnorm(0.025))*(sdmr)/(length(e))
ratciSUP=mr+(qnorm(0.975))*(sdmr)/(length(e))

mr
sdmr
ratciINF
ratciSUP
varrat=var(e)

###############################
###############################
###############################

eff=matrix(nrow=4,ncol=18)
colnames(eff)=c("NE","SE","XE","NC","SC","XC","Hedges","Inf","Sup","N0","N1","TF","MedRat","MeanRat","RINF","RSUP","SD","Var")

for(i in 1:4){
eff[i,1]=TAILLENE[i]
eff[i,2]=VALSE[i]
eff[i,3]=VALXE[i]
eff[i,4]=TAILLENC[i]
eff[i,5]=VALSC[i]
eff[i,6]=VALXC[i]
eff[i,7]=Effects[i]
eff[i,8]=I95Inf[i]
eff[i,9]=I95Sup[i]
eff[i,10]=lesrun0[i]
eff[i,11]=lesrun1[i]
eff[i,12]=lesrunt[i]
eff[i,13]=median(e)
eff[i,14]=mr
eff[i,15]=ratciINF
eff[i,16]=ratciSUP
eff[i,17]=sdmr
eff[i,18]=var(e)
}

write.table(eff,"effects.txt",dec=".",sep=" ")
