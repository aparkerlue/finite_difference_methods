source('functions.R')

k <- 10
r <- 0.04
sd <- 0.20
t <- 0.5
s <- 10
dt <- 0.002
dx <- sd*sqrt(4*dt)

m <- log(s) / dx
n <- t / dt
fde(k, r, sd, t, n, s, m)
