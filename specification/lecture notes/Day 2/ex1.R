# =============================================================================================================================================
# Attempt
library(deSolve)

SIR <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    dS <- -beta*I/1000*S
    dI <- beta*I/1000*S - r*I
    dR <- r*I
    output <- c(dS, dI, dR)
    list(output)
  })
}



start <- c(S=999, I=1, R=0)
parms <- c(beta=0.8, r=1/10)
times <- seq(0, 100, 1)
run_d <- ode(times=times, y=start, func=SIR, parms=parms)
summary(run_d)



matplot(run_d[,2:4], type="l")
plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model", ylab="Population")
lines(run_d[,3], col="blue")
lines(run_d[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))



library(animation)
saveGIF({
  beta <- seq(0, 5, 0.1)
  for (i in 2:length(beta)) {
    parms[1] <- beta[i]
    run_d <- ode(times=times, y=start, func=sir, parms = parms)
    plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
    lines(run_d[,3], col="blue")
    lines(run_d[,4], col="green")
    legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))
    Sys.sleep(0.3)
  }
},movie.name = "beta_sir.gif", internal=0.2)





# ------------------------------------------------ SEIR Model ------------------------------------------------ x
SEIR <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    dS <- -beta*I/1000*S
    dE <- beta*I/1000*S - a*E
    dI <- a*E - r*I
    dR <- r*I
    output <- c(dS, dE, dI, dR)
    list(output)
  })
}


start <- c(S=999, E=0, I=1, R=0)
parms <- c(beta=0.8, a=1/12, r=1/10)
times <- seq(0, 100, 1)
run_d <- ode(times=times, y=start, func=SEIR, parms=parms)

matplot(run_d[,2:4], type="l")
plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model", ylab="Population")
lines(run_d[,3], col="orange")
lines(run_d[,4], col="blue")
lines(run_d[,5], col="green")
legend("topright",legend=c("S", 'E', "I", "R"), col=c("red", 'orange', "blue", "green"), lty=c(1,1,1))


beta<-seq(0,5,0.1)
for (i in 1:length(beta)){
  parms[1]<-beta[i]
  run_d<-ode(times=times, y=start, func=SEIR,parms=parms)
  plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SEIR")
  lines(run_d[,3], col="blue")
  lines(run_d[,4], col="green")
  lines(run_d[,5], col="orange")
  legend("topright",legend=c("S","E", "I", "R"), col=c("red", "blue", "green", "orange"), lty=c(1,1,1,1))
  Sys.sleep(0.2)
}










# =============================================================================================================================================
# Solution



#IDM Exercise 1

#Differential Equations 
#SIR model for measles
library(deSolve)
sir <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dS=-beta*I/1000*S
    dI=beta*I/1000*S-r*I
    dR=r*I
    output <- c(dS, dI, dR)
    list(output)
  })
}

#the Initial values
start<-c(S=999, I=1,R=0 )

## The parameters 
parms <- c(beta=0.8,r=1/10)

## vector of timesteps
times <- seq(0, 100, 1)

run_d<-ode(times=times, y=start, func=sir,parms=parms)
#summary(run_d)

matplot(run_d[,2:4], type="l")
plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model", ylab="Population")
lines(run_d[,3], col="blue")
lines(run_d[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))


#(c)
library(animation)
saveGIF({
beta<-seq(0,5,0.1)
#r<-seq(1/20,1,0.05)
for (i in 1:length(beta)){
  parms[1]<-beta[i]
  # parms[2]<-r[i]
 run_d<-ode(times=times, y=start, func=sir,parms=parms)
 plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
 lines(run_d[,3], col="blue")
 lines(run_d[,4], col="green")
 legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))
 Sys.sleep(0.3)
}
}, movie.name = "beta_sir.gif", internal=0.2)

#SEIR model for measles
library(deSolve)
seir <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dS=-beta*I/1000*S
    dE=beta*I/1000*S-f*E
    dI=f*E-r*I
    dR=r*I
    output <- c(dS, dE, dI, dR)
    list(output)
  })
}

#the Initial values
start<-c(S=999, E=0, I=1,R=0 )

## The parameters 
parms <- c(beta=0.8,f=1/12,r=1/10)

## vector of timesteps
times <- seq(0, 100, 1)

run_d<-ode(times=times, y=start, func=seir,parms=parms)
#plot(run_d, ylim=c(0,1))
#summary(run_d)

plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SEIR")
lines(run_d[,3], col="blue")
lines(run_d[,4], col="green")
lines(run_d[,5], col="orange")
legend("topright",legend=c("S","E", "I", "R"), col=c("red", "blue", "green", "orange"), lty=c(1,1,1,1))

beta<-seq(0,5,0.1)
for (i in 1:length(beta)){
  parms[1]<-beta[i]
  run_d<-ode(times=times, y=start, func=seir,parms=parms)
  plot(run_d[,2], col="red", ylim=c(0,1000), type="l", main="SEIR")
  lines(run_d[,3], col="blue")
  lines(run_d[,4], col="green")
  lines(run_d[,5], col="orange")
  legend("topright",legend=c("S","E", "I", "R"), col=c("red", "blue", "green", "orange"), lty=c(1,1,1,1))
  Sys.sleep(0.2)
}


