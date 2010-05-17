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
