#Data fitting
#Least squares algorithm

data<-read.csv("/Users/sheetalsilal/Library/Mobile Documents/com~apple~CloudDocs/SHEETAL IDM/Mathematical Modelling/2019/Day 5/Measles.csv", header=T)
plot(data$measles)
plot(data$measles, type="b")

library(deSolve)
library(gtools)
sir <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    beta=exp(logbeta)
    mu=exp(logmu)
    gamma=exp(loggamma)
    rho=exp(logrho)
    amp=inv.logit(logitamp)
    phi=inv.logit(logitphi)
    betat=beta*(1+amp*cos((2*pi*(t/30 - phi))))
  
    dS=mu*pop-betat*I/pop*S +rho*R -mu*S
    dI=betat*I/pop*S-gamma*I -mu*I
    dR=gamma*I-mu*R -rho*R
    dInc=betat*I/pop*S
    output <- c(dS, dI, dR, dInc)
    list(output)
  })
}

pop=15500
#the Initial values
start<-c(S=pop-22, I=22,R=0, Inc=0 )

## The parameters 
parms <- c(logbeta=log(0.7),loggamma=log(1/10), logrho=log(1/12), logmu=log(1/100), logitamp=logit(1), logitphi=logit(0.03))

## vector of timesteps
times <- seq(0, 57, 1)

run_d<-ode(times=times, y=start, func=sir,parms=parms)
plot(data$measles, type="b", ylim=c(0, max(data$measles, run_d[,3])))
lines(times,run_d[,3],col="blue")
#summary(run_d)

inc<-c(run_d[1,5])
for (tt in 2:length(times)){
  inc<-c(inc,run_d[tt,5] - run_d[tt-1,5])
}

rep=0.7
plot(data$measles, type="b", ylim=c(0,max(data$measles, rep*inc)))
lines(times,rep*inc,col="blue")


sir.sse<-function(data,parms, rep){
  model<-ode(times=times, y=start, func=sir,parms=parms)
  inc<-c(model[1,5])
  for (tt in 2:length(times)){
    inc<-c(inc,model[tt,5] - model[tt-1,5])
  }
  error<-rep*inc[-1]-data$measles
  sse<-sum(error^2)
  return(sse)
}

sir.sse(data,parms, rep=0.7)

parms <- c(logbeta=log(0.7),
           loggamma=log(1/6), 
           logrho=log(1/25), 
           logmu=log(1/(50*52)), 
           logitamp=logit(0.9), 
           logitphi=logit(0.03))

run_d<-ode(times=times, y=start, func=sir,parms=parms)
matplot(run_d, type="l")

fit_sse<-optim(parms,sir.sse,data=data, rep=0.7)


parmtr<-function(parms){
  beta=exp(parms[1])
  mu=exp(parms[4])
  gamma=exp(parms[2])
  rho = exp(parms[3])
  amp=inv.logit(parms[5])
  phi=inv.logit(parms[6])
 output<-c(beta=beta, gamma=gamma, rho=rho, mu=mu, amp=amp, phi=phi) 
 return(output)
 }

bestpar<-parmtr(fit_sse$par)
startpar<-parmtr(parms)
cbind(startpar, bestpar)

run_d<-ode(times=times, y=start, func=sir,parms=parms)
inc<-c(run_d[1,5])
for (tt in 2:length(times)){
  inc<-c(inc,run_d[tt,5] - run_d[tt-1,5])
}

run_fit<-ode(times=times, y=start, func=sir,parms=fit_sse$par)
incfit<-c(run_fit[1,5])
for (tt in 2:length(times)){
  incfit<-c(incfit,run_fit[tt,5] - run_fit[tt-1,5])
}

par(mfrow=c(1,2))
plot(data$measles, type="b", ylim=c(0, max(rep*inc, data$measles)), main="Visual")
lines(times,rep*inc,col="blue")

plot(data$measles, type="b", ylim=c(0, max(rep*incfit, data$measles)), main="SSE")
lines(times,rep*incfit,col="blue")

##################################################################
#Poisson Log likelihood

pois_nll<-function(data,parms, rep){
  model<-ode(times=times, y=start, func=sir,parms=parms)
  inc<-c(model[1,5])
  for (tt in 2:length(times)){
    inc<-c(inc,model[tt,5] - model[tt-1,5])}
    ll<-sum(dpois(x=data$measles, lambda=rep*inc[-1],log=T ))
    nll<--1*ll
return(nll)  
}

pois_nll(data,parms, rep=0.7)

fit_ll<-optim(parms,pois_nll,data=data, rep=0.7)


parmtr<-function(parms){
  beta=exp(parms[1])
  mu=exp(parms[4])
  gamma=exp(parms[2])
  rho = exp(parms[3])
  amp=inv.logit(parms[5])
  phi=inv.logit(parms[6])
  output<-c(beta=beta, gamma=gamma, rho=rho, mu=mu, amp=amp, phi=phi) 
  return(output)
}

bestpar_sse<-bestpar#parmtr(fit_sse$par)
bestpar_ll<-parmtr(fit_ll$par)
cbind(bestpar_sse, bestpar_ll)

run_sse<-ode(times=times, y=start, func=sir,parms=fit_sse$par)
inc_sse<-c(run_sse[1,5])
for (tt in 2:length(times)){
  inc_sse<-c(inc,run_sse[tt,5] - run_sse[tt-1,5])
}

run_ll<-ode(times=times, y=start, func=sir,parms=fit_ll$par)
inc_ll<-c(run_ll[1,5])
for (tt in 2:length(times)){
  inc_ll<-c(incfit,run_ll[tt,5] - run_ll[tt-1,5])
}

par(mfrow=c(1,2))
plot(data$measles, type="b", ylim=c(0, max(rep*inc, data$measles)), main="SSE")
lines(times,rep*inc_sse[-1],col="blue")

plot(data$measles, type="b", ylim=c(0, max(rep*incfit, data$measles)), main="PoisLL")
lines(times,rep*inc_ll[-1],col="blue")
