source('functions.R')

k <- 10
r <- 0.04
sd <- 0.20
t <- 0.5
s <- 10
dt <- 0.002
dx <- sd*sqrt(4*dt)

n <- t / dt
m <- round(log(s) / dx)
fde(s, k, r, t, sd, n, m)

##

optionrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='put',
                 style='american')
with(optionrc, {
  p.bsm <- bsm.option(s=s, k=k, r=r, t=t, sd=sd, type=type)
  p.fde <- fde(s=s, k=k, r=r, t=t, sd=sd, n=10, m=20, type=type, style=style)
  cat('bsm = ', p.bsm, '\n', sep="")
  cat('fde = ', p.fde, '\n', sep="")
})
