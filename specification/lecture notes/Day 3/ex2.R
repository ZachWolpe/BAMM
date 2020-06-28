# =============================================================================================================================================
# Attempt
library(deSolve)


SIR_bd <- function(t, x, parms) {
  with(list(c(parms, x)), {
    dS <- mu*1000 - beta*I/1000*S +p*R - mu*S
    dI <- beta*I/1000*S - r*I - mu*I
    dR <- r*I - p*R - mu*I
    output <- c(dS, dI, dR)
    list(output)
  })
}






start <- c(S=999, I=1, R=0)
parms <- c(mu=0.01, beta=0.8, r=1/10)
times <- seq(0, 1000, 1)
run_d <- ode(times=times, y=start, func=SIR_bd, parms=parms)
summary(run_d)



matplot(run_d[,2:4], type="l")
plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model", ylab="Population")
lines(run_d[,3], col="blue")
lines(run_d[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))
















# =============================================================================================================================================
# Solution


#Exercise 2: Comparing 3 models for treatment
library(deSolve)

# SIRS model with birth and death
sirsbd <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dS=mu*1000-beta*I/1000*S + p*R -mu*S
    dI=beta*I/1000*S-r*I -mu*I
    dR=r*I -p*R -mu*R
    output <- c(dS, dI, dR)
    list(output)
  })
}

#the Initial values
start<-c(S=999, I=1,R=0 )

## The parameters 
parms <- c(beta=0.8,r=1/25, p=1/50, mu=1/100)

## vector of timesteps
times <- seq(0, 1000, 1)

run_d<-ode(times=times, y=start, func=sirsbd,parms=parms)
#summary(run_d)

plot(run_d[,2], col="red", ylim=c(0,1200), type="l", main="SIR Model with b/d")
lines(run_d[,3], col="blue")
lines(run_d[,4], col="green")
lines(run_d[,2]+run_d[,3]+run_d[,4])
legend("topright",legend=c("S", "I", "R", "N"), col=c("red", "blue", "green", "black"), lty=c(1,1,1,1))

#************************************************************************

#Model 1

sirsbdt1 <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    dS=mu-beta*I*S +p*R -mu*S
    dI=beta*I*S-ppi*r*I-(1-ppi)*a*I-mu*I
    dR=ppi*r*I+(1-ppi)*a*I-p*R -mu*R
    output <- c(dS, dI, dR)
    list(output)
  })
}

#the Initial values
start<-c(S=0.8, I=0.2,R=0 )

## The parameters 
parms <- c(beta=0.8,r=1/(5+2), p=1/50, mu=1/100, ppi=0.1,a=1/25)

## vector of timesteps
times <- seq(0, 1000, 1)

t1<-ode(times=times, y=start, func=sirsbdt1,parms=parms)

plot(t1[,2], col="red", ylim=c(0,1), type="l", main="SIRS bd t1")
lines(t1[,3], col="blue")
lines(t1[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))

#Steady state conditions
beta=0.8;r=1/(5+2); p=1/50; mu=15/100; k=0.5; ppi=0.3; a=1/20
ss<-c(S=(ppi*r+(1-ppi)*a+mu)/beta, RdivI=(ppi*r+(1-ppi)*a)/(p+mu) )
t1[1001,2]
t1[1001,4]/t1[1001,3]
ss
R0<-beta/(ppi*r+(1-ppi)*a+mu)
R0

##################################################
sirsbdt2 <- function(t, x, parms2, input)  {
  with(as.list(c(parms2, x)), {
    dS=mu-beta*(I)*S +p*R -mu*S
    dI=beta*(I)*S-ppi*tau*I-(1-ppi)*a*I-mu*I
    dIt=ppi*tau*I-t*It -mu*It
    dR=t*It+(1-ppi)*a*I-p*R -mu*R
    output <- c(dS, dI, dIt, dR)
    list(output)
  })
}

#the Initial values
start2<-c(S=0.999, I=0.001,It=0, R=0 )

## The parameters 
parms2 <- c(beta=0.8,tau=1/2, t=1/5, p=1/50, mu=1/100, k=0.5, ppi=0,a=1/25)

## vector of timesteps
times <- seq(0, 1000, 1)

t2<-ode(times=times, y=start2, func=sirsbdt2,parms2=parms2)

plot(t2[,2], col="red", ylim=c(0,1), type="l", main="SIRS bd t2")
lines(t2[,3], col="blue")
lines(t2[,4], col="purple")
lines(t2[,5], col="green")
legend("topright",legend=c("S", "I","It", "R"), col=c("red", "blue","purple", "green"), lty=c(1,1,1,1))

#Steady state conditions
beta=0.8;tau=1/2; t=1/5; p=1/50; mu=15/100; k=0.5; ppi=0.3;a=1/20
t2[1001,]
ss2<-c(S=((1-ppi)*a+ppi*tau+mu)*(t+mu)/(ppi*tau*k*beta)/(1+(t+mu)/(k*ppi*tau)), ItdivI=ppi*tau/(t+mu))
ss2
t2[1001,4]/t2[1001,3]
R0<-beta/((1-ppi)*a+ppi*tau+mu)
R0


##########################################################
sirsbdt3 <- function(t, x, parms3, input)  {
  with(as.list(c(parms3, x)), {
    dS=mu-beta*(I)*S +p*R -mu*S
    dI=(1-ppi)*beta*(I)*S-a*I-mu*I
    dIt=ppi*beta*(I)*S-r*It -mu*It
    dR=r*It+a*I-p*R -mu*R
    output <- c(dS, dI, dIt, dR)
    list(output)
  })
}

#the Initial values
start3<-c(S=0.8, I=0.2,It=0, R=0 )

## The parameters 
parms3 <- c(beta=0.8,r=1/7, p=1/50, mu=15/100, k=0.5, ppi=0.3,a=1/20)

## vector of timesteps
times <- seq(0, 1000, 1)

t3<-ode(times=times, y=start3, func=sirsbdt3,parms3=parms3)

plot(t3[,2], col="red", ylim=c(0,1), type="l", main="SIRS bd t3")
lines(t3[,3], col="blue")
lines(t3[,4], col="purple")
lines(t3[,5], col="green")
legend("topright",legend=c("S", "I","It", "R"), col=c("red", "blue","purple", "green"), lty=c(1,1,1,1))

#Steady state conditions
beta=0.8;r=1/7; p=1/50; mu=15/100; k=0.5; ppi=0.3;a=1/20
t3[1001,]
ss3<-c(S=(ppi*beta/(r+mu)*(k+(1-ppi)*(r+mu)/(ppi*(a+mu))))^-1,  ItdivI=(ppi*(a+mu))/((1-ppi)*(r+mu)))
ss3
t3[1001,4]/t3[1001,3]

#R0 comparison
beta=0.8;r_comp=1/(5+2); p=1/50; mu=15/100; ppi=0.3;a=1/20; tau=1/2; r=1/5;

R0_1=beta/(ppi*r_comp+(1-ppi)*a+mu)
R0_2=beta/(ppi*tau+(1-ppi)*a+mu)
R0_3=(1-ppi)*beta/(a+mu)

trt<-seq(0,1,0.02)
R0_1v<-R0_2v<-R0_3v<-c()

for (i in 1:length(trt)){
  R0_1v=c(R0_1v,beta/(trt[i]*r_comp+(1-trt[i])*a+mu))
  R0_2v=c(R0_2v,beta/(trt[i]*tau+(1-trt[i])*a+mu))
  R0_3v=c(R0_3v,(1-trt[i])*beta/(a+mu))
}
plot(trt, R0_1v, lwd=2,type="b", col="red", main="R0 comparison by treatment probability", ylim=c(0,5))
lines(trt, R0_2v,lwd=2, type="b", col="blue")
lines(trt, R0_3v,lwd=2, type="b", col="green")
legend("topright",legend=c("M1", "M2", "M3"), col=c("red", "blue", "green"), lty=c(1,1,1))
