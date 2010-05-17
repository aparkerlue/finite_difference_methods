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

s <- 50
k <- 50
r <- 0.10
sd <- 0.40
t <- 5/12

bsm.option(s, k, r, t, sd, type='call')
fde.s(s, k, r, t, sd, n=10, m=20, type='put', style='american')
