#Exercise 6: Interventions
library(deSolve)
# SIR model with birth and death
sirsv <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    betat=beta*(1+beta1*cos(2*pi*t))
    if (t>10) {propvt=propv}
    else (propvt=0)
    dS=mu*1000-betat*I/1000*S*(1-propvt) -mu*S+rho*R -v*propvt*S
    dI=betat*I/1000*S*(1-propvt) -gamma*I -mu*I
    dR=gamma*I -mu*R - rho*R
    dV=v*propvt*S - mu*V
    dVacpop=v*propvt*S
    output <- c(dS, dI, dR, dV, dVacpop)
    list(output)
  })
}

#the Initial values
start<-c(S=999, I=1,R=0, V=0, Vacpop=0)

## The parameters 
parms <- c(beta=300, beta1=0.4,gamma=365/10, mu=365/100, rho=365/30, propv=0.2,v=365/7)

## vector of timesteps
times <- seq(0, 20, 1/365)

run_d<-ode(times=times, y=start, func=sirsv,parms=parms)
matplot(run_d[,2:5], type="l")
 
plot(times, run_d[,6], col="purple", type="l")

#Plot for different levels of vaccination
parms <- c(beta=300, beta1=0.4,gamma=365/10, mu=365/100, rho=365/30, propv=0,v=365/7)
run_d<-ode(times=times, y=start, func=sirsv,parms=parms)
plot(times, run_d[,3], type="l")
for (i in 1:10){
  parms <- c(beta=300, beta1=0.4,gamma=365/10, mu=365/100, rho=365/30, propv=i*0.1,v=365/7)
  run_d<-ode(times=times, y=start, func=sirsv,parms=parms)
  lines(times, run_d[,3], col=(i+1))
}

#Plot for Cumulative counts of vaccinations
parms <- c(beta=300, beta1=0.4,gamma=365/10, mu=365/100, rho=365/30, propv=0,v=365/7)
run_d<-ode(times=times, y=start, func=sirs,parms=parms)
plot(times, run_d[,6], type="l", ylim=c(0,1000))
for (i in 1:10){
  parms <- c(beta=300, beta1=0.4,gamma=365/10, mu=365/100, rho=365/30, propv=i*0.1,v=365/7)
  run_d<-ode(times=times, y=start, func=sirs,parms=parms)
  lines(times, run_d[,6], col=(i+1))
}

#calculate Incidence of vaccines
vinc<-c(run_d[1,6])
for (tt in 2:length(times)){
vinc<-c(vinc,run_d[tt,6] - run_d[tt-1,6])
}

#######################################################################

#Part B

# SIR model with birth and death
sirsv <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    betat=beta*(1+beta1*cos(2*pi*t))
    if (t>10) {propvt=propv}
    else (propvt=0)
    dS=mu*1000-betat*I/1000*S*(1-propvt) -mu*S+rho*R -v*propvt*S
    dI=betat*I/1000*S*(1-propvt)-gamma*I -mu*I
    dR=gamma*I -mu*R - rho*R
    dV=v*propvt*S - mu*V
    dVacpop=v*propvt*S
    output <- c(dS, dI, dR, dV, dVacpop)
    list(output)
  })
}


#the Initial values
start<-c(S=999, I=1,R=0, V=0, Vacpop=0)

## The parameters 
parms <- c(beta=300, beta1=0.4,gamma=365/10, mu=365/100, rho=365/30, propv=0.2,v=365/7)

## vector of timesteps
times <- seq(0, 20, 1/365)

run_d<-ode(times=times, y=start, func=sirsv,parms=parms)

par(mfrow=c(1,1))
plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
lines(run_d[,3], col="blue")
lines(run_d[,4], col="green")
lines(run_d[,5], col="purple")
lines(run_d[,2]+run_d[,3]+run_d[,4]+run_d[,5])
legend("topright",legend=c("S", "I", "R", "V", "N"), col=c("red", "blue", "green", "purple", "black"), lty=c(1,1,1,1,1))


for ( i in 1:10){
  parms <- c(beta=runif(1,250,350), beta1=runif(1,0.3, 0.5),gamma=365/runif(1,8,12), mu=365/runif(1,80,120), rho=365/runif(1,25,35), propv=0.2,v=365/runif(1,5,9))
  run_d<-ode(times=times, y=start, func=sirsv,parms=parms)
  lines(run_d[,2], col="red")
  lines(run_d[,3], col="blue")
  lines(run_d[,4], col="green")
  lines(run_d[,5], col="purple")
  
}

CIdata<-NULL
for ( i in 1:20){
  parms <- c(beta=runif(1,250,350), beta1=runif(1,0.3, 0.5),gamma=365/runif(1,8,12), mu=365/runif(1,80,120), rho=365/runif(1,25,35), propv=0.2,v=365/runif(1,5,9))
  run_d<-ode(times=times, y=start, func=sirsv,parms=parms)
  CIdata<-cbind(CIdata,run_d[,3])
}
CIsd<-CIuci<-CIlci<-NULL
CIMean<-rowMeans(CIdata)
for (i in 1:(dim(CIdata)[1])){
  CIsd[i]<-sd(CIdata[i,])
  CIuci[i]<-CIMean[i]+1.96*CIsd[i]/sqrt(20)
  CIlci[i]<-CIMean[i]-1.96*CIsd[i]/sqrt(20)
}

plot(times, CIMean, type="l")
lines(times, CIlci, col="red")
lines(times, CIuci, col="red")
polygon(c(times,rev(times)), c(CIuci, rev(CIlci)), col="red")

