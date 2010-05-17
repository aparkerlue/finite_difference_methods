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

optrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='p', style='a')
with(optrc, {
  p.bsm <- bsm.option(s=s, k=k, r=r, t=t, sd=sd, type=type)
  p.fde <- fde(s=s, k=k, r=r, t=t, sd=sd, n=10, m=20, type=type, style=style)
  cat('bsm = ', p.bsm, '\n', sep="")
  cat('fde =\n')
  print(p.fde)
})
with(optrc, fde(s, k, r, t, sd, n=10, m=20, type=type, style=style))
with(optrc, rotate(fde(s,k,r,t,sd,n=10,m=20,type,style,grid=T), digits=2))


optrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='p', style='a')
## Binomial
with(optrc, binom.american.put(s, k, r, t, sd, n=1e3))
## Explicit finite difference
with(optrc, fde(s, k, r, t, sd, n=10, m=20, type=type, style=style))
with(optrc, rotate(fde(s,k,r,t,sd,n=10,m=20,type,style,grid=T), digits=2))
## Explicit finite difference, log transform
with(optrc, fde.log(s, k, r, t, sd, n=10, m=20, type, style))
with(optrc, rotate(fde.log(s,k,r,t,sd,n=10,m=20,type,style,grid=T), 2))

binom.american.put(exp(3.145760), 50, .1, 5/12, .4, n=1e3)
binom.american.put(exp(7.694953), 50, .1, 5/12, .4, n=1e3)



