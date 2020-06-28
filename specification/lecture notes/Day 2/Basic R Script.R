#Basic Differential Equation Script

#Differential Equations 
#SIR model for measles
library(deSolve)

sir <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dS=-beta*I/1000*S
	.....
    output <- c(dS, ....)
    list(output)
  })
}

#the Initial values
start<-c(S=999, .... )

## The parameters 
parms <- c(beta=0.8, ....)

## vector of timesteps
times <- seq(0, 50, 1/30)

run<-ode(times=times, y=start, func=sir,parms=parms)

plot(run[,2], col="red", ylim=c(0,1000), type="l", main="SIR Model")
lines(run[,3], col="blue")
lines(run[,4], col="green")
legend("topright",legend=c("S", "I", "R"), col=c("red", "blue", "green"), lty=c(1,1,1))
