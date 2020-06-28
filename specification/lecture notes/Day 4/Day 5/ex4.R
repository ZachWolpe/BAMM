#SIR model for measles
library(deSolve)
#Exercise 4 Part 1
####################################################################
#Part 1
t<-seq(0,1000,1)
beta1<-1
beta0<-1
betat<-beta0*(1+beta1*cos(2*pi*t))
par(mfrow=c(1,1))
plot(t, betat, type="l")

beta0range<-seq(0,5,0.25)
beta1range<-seq(0,1,0.025)

for (i in 1:length(beta0range)){
	plot(t,beta0range[i]*(1+beta1*cos(2*pi*t/365)),ylim=c(0,11), type="l")
	Sys.sleep(0.2)
}

for (i in 1:length(beta1range)){
	plot(t,beta0*(1+beta1range[i]*cos(2*pi*t/365)),ylim=c(0,2.5), type="l")
	Sys.sleep(0.2)
}


####################################################################

#SIR model with treatment in two patches
sirsbdtp2s <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    N1=S1+I1+T1+R1
    N2=S2+I2+T2+R2
    beta1t<-beta1*(1+0.5*cos(2*pi*t/365))
    beta2t<-beta2*(1+0.5*cos(2*pi*t/365))
    #Patch 1 equations
    dS1=-beta1t*(I1/N1)*S1-mu*S1+p*R1 + mu*N1+m21*S2 -m12*S1
    dI1=beta1t*I1/N1*S1-(1-pi1)*gamma*I1 - mu*I1 - tau*pi1*I1 - m12*I1 +m21*I2
    dT1= tau*pi1*I1 - mu*T1 -r1*T1 -m12*T1 +m21*T2
    dR1=(1-pi1)*gamma*I1 -p*R1 -mu*R1 +r1*T1 -m12*R1 +m21*R2
    
    #Patch 2 equations
    dS2=-beta2t*(I2/N2)*S2-mu*S2 +p*R2 + mu*N2-m21*S2 +m12*S1
    dI2=beta2t*I2/N2*S2-(1-pi2)*gamma*I2 - mu*I2 - tau*pi2*I2 + m12*I1 -m21*I2
    dT2= tau*pi2*I2 - mu*T2 -r2*T2 +m12*T1 -m21*T2
    dR2=(1-pi2)*gamma*I2 -mu*R2 -p*R2 +r2*T2 +m12*R1 -m21*R2
    
    output <- c(dS1, dI1,dT1,dR1, dS2, dI2, dT2, dR2)
    list(output)
  })
}

#the Initial values
start<-c(S1=1000, I1=0,T1=0, R1=0, S2=999, I2=1, T2=0, R2=0 )

## The parameters 
parms <- c(beta1=0.8,beta2=0.8, gamma=1/10, mu=1/100, tau=1/3, r1=1/2,r2=1/5,p=1/50, pi1=0.3, pi2=0.3, m21=1/20, m12=1/20)

## vector of timesteps
times <- seq(0, 1000, 1)

run_d<-ode(times=times, y=start, func=sirsbdtp2s,parms=parms)

pop1<-run_d[,2]+run_d[,3]+run_d[,4]+run_d[,5]
pop2<-run_d[,6]+run_d[,7]+run_d[,8]+run_d[,9]
  
par(mfrow=c(1,2))

plot(times,run_d[,2], col="red", ylim=c(0,2300), type="l", main="Patch 1: r1=1/2")
lines(times,run_d[,3], col="blue")
lines(times, run_d[,4], col="green")
lines(times,run_d[,5], col="purple")
lines(times,pop1,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))

plot(times,run_d[,6], col="red", ylim=c(0,2300), type="l", main="Patch 2: r2=1/5")
lines(times,run_d[,7], col="blue")
lines(times, run_d[,8], col="green")
lines(times,run_d[,9], col="purple")
lines(times,pop2,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))


matplot(run_d,type="l")


####################################################################
#Part 2
t<-seq(0,1000,1)
beta1<-1
beta0<-1
betat<-beta0*(1+beta1*cos(2*pi*t/365))
plot(t, betat, type="l")

#Import elnino data
elnino_sampledata <- read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/SHEETAL IDM/Mathematical Modelling/2019/Day 4/elnino_sampledata.txt")
plot(elnino_sampledata$nino34, type="l")
length(elnino_sampledata$nino34) #monthly
#132, 11 years
t<-seq(0,132*30,1)
betat<-beta0*(1+beta1*cos(2*pi*t/365))
plot(t, betat, type="l")

tm<-seq(0,131,1)
#Compute new beta
betat_el<-elnino_sampledata$nino34*(beta0*(1+beta1*cos(2*pi*tm/12)))
plot(tm, betat_el, type="l")

#Compare side by side
par(mfrow=c(1,2))
plot(t, betat, type="l")
plot(tm, betat_el, type="l")

#Use standardised el nino figures 
betat_el_std<-elnino_sampledata$stdnino34*(beta0*(1+beta1*cos(2*pi*tm/12)))
plot(t, betat, type="l")
plot(tm, betat_el_std, type="l")

#Compare to el Nino pattern
plot(elnino_sampledata$nino34, type="l")
plot(elnino_sampledata$stdnino34, type="l")

#Smooth el nino
par(mfrow=c(4,4))
plot(tm, runmed(elnino_sampledata$stdnino34, 1), type="l", ylab="", col=1)
for (i in 1:15){
plot(tm, runmed(elnino_sampledata$stdnino34, 1+2*i), ylab="", col=1, type="l")
}

#Run your model 
# parmsM <- c(beta1=30*0.8,beta2=30*0.8, gamma=1/(10/30), mu=1/(100/30), tau=1/(3/30), r1=1/(2/30),r2=1/(5/30),p=1/(50/30), pi1=0.3, pi2=0.3, m21=1/(20/30), m12=1/(20/30))


sirsbdtp2sEN <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    N1=S1+I1+T1+R1
    N2=S2+I2+T2+R2
    eln<-approx(tm,elnino_sampledata$stdnino34,t/30)$y
    beta1t<-eln*(beta1*(1+0.5*cos(2*pi*t/365)))
    beta2t<-eln*(beta2*(1+0.5*cos(2*pi*t/365)))

    #Patch 1 equations
    dS1=-beta1t*(I1/N1)*S1-mu*S1+p*R1 + mu*N1+m21*S2 -m12*S1
    dI1=beta1t*I1/N1*S1-(1-pi1)*gamma*I1 - mu*I1 - tau*pi1*I1 - m12*I1 +m21*I2
    dT1= tau*pi1*I1 - mu*T1 -r1*T1 -m12*T1 +m21*T2
    dR1=(1-pi1)*gamma*I1 -p*R1 -mu*R1 +r1*T1 -m12*R1 +m21*R2
    
    #Patch 2 equations
    dS2=-beta2t*(I2/N2)*S2-mu*S2 +p*R2 + mu*N2-m21*S2 +m12*S1
    dI2=beta2t*I2/N2*S2-(1-pi2)*gamma*I2 - mu*I2 - tau*pi2*I2 + m12*I1 -m21*I2
    dT2= tau*pi2*I2 - mu*T2 -r2*T2 +m12*T1 -m21*T2
    dR2=(1-pi2)*gamma*I2 -mu*R2 -p*R2 +r2*T2 +m12*R1 -m21*R2
    
    output <- c(dS1, dI1,dT1,dR1, dS2, dI2, dT2, dR2)
    list(output)
  })
}
times2 <- seq(0, 131*30, 1)

run_d<-ode(times=times2, y=start, func=sirsbdtp2sEN,parms=parms)

pop1<-run_d[,2]+run_d[,3]+run_d[,4]+run_d[,5]
pop2<-run_d[,6]+run_d[,7]+run_d[,8]+run_d[,9]

par(mfrow=c(1,2))

plot(times2,run_d[,2], col="red", ylim=c(0,2300), type="l", main="Patch 1: r1=1/2")
lines(times2,run_d[,3], col="blue")
lines(times2, run_d[,4], col="green")
lines(times2,run_d[,5], col="purple")
lines(times2,pop1,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))

plot(times2,run_d[,6], col="red", ylim=c(0,2300), type="l", main="Patch 2: r2=1/5")
lines(times2,run_d[,7], col="blue")
lines(times2, run_d[,8], col="green")
lines(times2,run_d[,9], col="purple")
lines(times2,pop2,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))

par(mfrow=c(3,1))
plot(times2,run_d[,3], col="blue", ylim=c(0,400), type="l", ylab="Infections")
plot(elnino_sampledata$stdnino34, type="l", ylab="elnin0std")
plot(betat, type="l")


