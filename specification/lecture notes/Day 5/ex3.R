####################################################################
#exercise 3
#SIR model with treatment in two patches
library(deSolve)
sirsbdtp2 <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    N1=S1+I1+T1+R1
    N2=S2+I2+T2+R2
    #Patch 1 equations
    dS1=-beta1*(I1/N1)*S1-mu*S1 +p*R1 + mu*N1+m21*S2 -m12*S1
    dI1=beta1*I1/N1*S1-(1-pi1)*gamma*I1 - mu*I1 - tau*pi1*I1 - m12*I1 +m21*I2
    dT1= tau*pi1*I1 - mu*T1 -r1*T1 -m12*T1 +m21*T2
    dR1=(1-pi1)*gamma*I1 -p*R1 -mu*R1 +r1*T1 -m12*R1 +m21*R2

    #Patch 2 equations
    dS2=-beta2*(I2/N2)*S2-mu*S2 +p*R2 + mu*N2-m21*S2 +m12*S1
    dI2=beta2*I2/N2*S2-(1-pi2)*gamma*I2 - mu*I2 - tau*pi2*I2 + m12*I1 -m21*I2
    dT2= tau*pi2*I2 - mu*T2 -r2*T2 +m12*T1 -m21*T2
    dR2=(1-pi2)*gamma*I2 -p*R2 -mu*R2 +r2*T2 +m12*R1 -m21*R2

    output <- c(dS1, dI1,dT1,dR1, dS2, dI2, dT2, dR2)
    list(output)
  })
}

#the Initial values
start<-c(S1=900, I1=20,T1=10, R1=70, S2=900, I2=20, T2=10, R2=70 )

## The parameters 
parms <- c(beta1=0.4,beta2=0.4, gamma=1/50, mu=1/100, tau=1/2, r1=1/5,r2=1/5, p=1/50, pi1=0.6, pi2=0.3, m21=1/20, m12=1/20)

## vector of timesteps
times <- seq(0, 1000, 1)

run_d<-ode(times=times, y=start, func=sirsbdtp2,parms=parms)

pop1<-run_d[,2]+run_d[,3]+run_d[,4]+run_d[,5]
pop2<-run_d[,6]+run_d[,7]+run_d[,8]+run_d[,9]

par(mfrow=c(1,2))
#Number plots
plot(times,run_d[,2], col="red", ylim=c(0,2000), type="l", main="Patch 1")
lines(times,run_d[,3], col="blue")
lines(times, run_d[,4], col="green")
lines(times,run_d[,5], col="purple")
lines(times,pop1,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))

plot(times,run_d[,6], col="red", ylim=c(0,2000), type="l", main="Patch 2")
lines(times,run_d[,7], col="blue")
lines(times, run_d[,8], col="green")
lines(times,run_d[,9], col="purple")
lines(times,pop2,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))
# 
#Percentage plots
plot(times,run_d[,2]/pop1, col="red", ylim=c(0,1), type="l", main="Patch 1")
lines(times,run_d[,3]/pop1, col="blue")
lines(times, run_d[,4]/pop1, col="green")
lines(times,run_d[,5]/pop1, col="purple")
lines(times,pop1,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))

plot(times,run_d[,6]/pop2, col="red", ylim=c(0,1), type="l", main="Patch 2")
lines(times,run_d[,7]/pop2, col="blue")
lines(times, run_d[,8]/pop2, col="green")
lines(times,run_d[,9]/pop2, col="purple")
lines(times,pop2,type="l")
#legend("topright",legend=c("S", "I","T","R"), col=c("red", "blue", "green","purple"), lty=c(1,1,1,1))

