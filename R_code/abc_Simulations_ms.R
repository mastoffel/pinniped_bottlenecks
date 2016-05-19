vectrap=seq(1,15,0.1) ##Asymetry a
vectmig=seq(0.01,0.1,0.005) ##downstream-directed migration rate
vecty=seq(100,10000,20)

resultat=matrix(nrow=(length(vectrap)*length(vectmig)*length(vecty)),ncol=11)
colnames(resultat)=c("theta","N","mAm","mAv","rapm","m","cor","var","moy","alpha","tstu")
co=0 

length(vectrap)*length(vectmig)*length(vecty)

for(y in vecty){
for(x in 1:length(vectrap)){
for(z in 1:length(vectmig)){

co=co+1

theta=4*0.000556*y
mAv=theta/0.000556*(vectmig[z])
mAm=mAv*(1/vectrap[x])

##passing parameter values to the ms command, then executing ms as an external program

eval(parse(text=paste('commande=paste(Sys.getenv("COMSPEC"),"/c","ms 440 15 -t ",theta," -I 10 44 44 44 44 44 44 44 44 44 44 -n 1 1 -n 2 1 -n 3 1 -n 4 1 -n 5 1 -n 6 1 -n 7 1 -n 8 1 -n 9 1 -n 10 1 -m 1 2 ",mAm," -m 2 3 ",mAm," -m 3 4 ",mAm," -m 4 5 ",mAm," -m 5 6 ",mAm," -m 6 7 ",mAm," -m 7 8 ",mAm," -m 8 9 ",mAm," -m 9 10 ",mAm," -m 10 9 ",mAv," -m 9 8 ",mAv," -m 8 7 ",mAv," -m 7 6 ",mAv," -m 6 5 ",mAv," -m 5 4 ",mAv," -m 4 3 ",mAv," -m 3 2 ",mAv," -m 2 1 ",mAv," | microsat.exe > Test_ms.txt")',sep='')))
system(commande)

resultat[co,1]=theta
resultat[co,2]=(theta/0.000556)/4
resultat[co,3]=mAm
resultat[co,4]=mAv
resultat[co,5]=vectrap[x]
resultat[co,6]=vectmig[z]

sim=read.table("Test_ms.txt",h=F) #reads the ms output

pop1=sim[,1:44]
for(i in 2:10){
eval(parse(text=paste("pop",i,"=sim[,(((",i,"-1)*44)+1):(",i,"*44)]",sep="")))
}

ARperlocus=rep(0,15)
MeanARperpop=rep(0,10)

for(k in 1:10){
for(b in 1:15){
eval(parse(text=paste("ARperlocus[",b,"]=length(unique(as.character(pop",k,"[",b,",])))",sep="")))
eval(parse(text=paste("MeanARperpop[",k,"]=mean(ARperlocus)",sep="")))
}
}
calc=matrix(nrow=2,ncol=10)
calc[1,1:10]=seq(1,10,by=1)
calc[2,1:10]=MeanARperpop

resultat[co,7]=cor(calc[1,],calc[2,])
resultat[co,8]=var(MeanARperpop)
resultat[co,9]=mean(MeanARperpop)
mod=lm(scale(calc[2,])~calc[1,])
resultat[co,10]=mod$coefficients[2]
tstudent=summary(mod)
resultat[co,11]=tstudent$coefficients[2,3]
}
}
}

write.table(resultat,"Simuls.txt", dec=".",sep=" ")

#####################################################################################################
#####################################################################################################


