source('functions.R')

######################################################################
## Problem 1                                                        ##
######################################################################

s <- 10
k <- 10
r <- .04
t <- .5
sd <- .2
dt <- .002
dx <- c(sd*sqrt(dt), sd*sqrt(3*dt), sd*sqrt(4*dt))
type <- 'put'; style <- 'european'

## Black-Scholes-Merton for comparison
bsm.option(s, k, r, t, sd, type)

## Regular implementation.
fde.log(s, k, r, t, sd, type=type, style=style)
fdi.log(s, k, r, t, sd, type=type, style=style)
fdcn.log(s, k, r, t, sd, type=type, style=style)

## Project requirements.
n <- round(t/dt)
m <- round(6*sd*sqrt(t)/dx[1])          # dx = sd*sqrt(dt)
fde.log(s, k, r, t, sd, n, m, type, style)
fdi.log(s, k, r, t, sd, n, m, type, style)
fdcn.log(s, k, r, t, sd, n, m, type, style)

m <- round(6*sd*sqrt(t)/dx[2])          # dx = sd*sqrt(3*dt)
fde.log(s, k, r, t, sd, n, m, type, style)
fdi.log(s, k, r, t, sd, n, m, type, style)
fdcn.log(s, k, r, t, sd, n, m, type, style)

m <- round(6*sd*sqrt(t)/dx[3])          # dx = sd*sqrt(4*dt)
fde.log(s, k, r, t, sd, n, m, type, style)
fdi.log(s, k, r, t, sd, n, m, type, style)
fdcn.log(s, k, r, t, sd, n, m, type, style)
