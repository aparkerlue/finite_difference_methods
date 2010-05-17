source('functions.R')
source('util.R')

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

optionrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='p', style='a')
with(optionrc, {
  p.bsm <- bsm.option(s=s, k=k, r=r, t=t, sd=sd, type=type)
  p.fde <- fde(s=s, k=k, r=r, t=t, sd=sd, n=10, m=20, type=type, style=style)
  cat('bsm = ', p.bsm, '\n', sep="")
  cat('fde =\n')
  print(p.fde)
})
with(optionrc, fde(s, k, r, t, sd, n=10, m=20, type=type, style=style))
with(optionrc, rotate(fde(s,k,r,t,sd,n=10,m=20,type,style,grid=T), digits=2))


optionrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='p', style='a')
## Binomial
with(optionrc, binom.american.put(s, k, r, t, sd, n=1e3))
## Explicit finite difference
with(optionrc, fde(s, k, r, t, sd, n=10, m=20, type=type, style=style))
with(optionrc, rotate(fde(s,k,r,t,sd,n=10,m=20,type,style,grid=T), digits=2))
## Explicit finite difference, log transform
with(optionrc, fde.log(s, k, r, t, sd, n=10, m=20, type, style))
with(optionrc, rotate(fde.log(s,k,r,t,sd,n=10,m=20,type,style,grid=T), 2))
