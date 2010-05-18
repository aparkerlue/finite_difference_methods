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


######################################################################
## Problem 2                                                        ##
######################################################################

s <- 10
k <- 10
r <- .04
t <- .5
sd <- .2
dt <- .002
ds <- c(0.5, 1, 1.5)
type <- 'put'; style <- 'american'

## Binomial or Black-Scholes-Merton for comparison
if (type == 'put') {
  binom.american.put(s, k, r, t, sd, n=1e3)
} else
  bsm.option(s, k, r, t, sd, type)

## Regular implementation.
fde(s, k, r, t, sd, type=type, style=style)
fdi(s, k, r, t, sd, type=type, style=style)

## Project requirements.
n <- round(t/dt)
m <- round(2*s/ds[1])                   # ds = 0.5
fde(s, k, r, t, sd, n, m, type, style)
fdi(s, k, r, t, sd, n, m, type, style)

m <- round(2*s/ds[2])                   # ds = 1
fde(s, k, r, t, sd, n, m, type, style)
fdi(s, k, r, t, sd, n, m, type, style)

m <- round(2*s/ds[3])                   # ds = 1.5
fde(s, k, r, t, sd, n, m, type, style)
fdi(s, k, r, t, sd, n, m, type, style)
