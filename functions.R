######################################################################
## Explicit Finite Difference Method                                ##
######################################################################

fde <- function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(n)),
                type = c("call", "put"), style = c("european", "american")) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type)
  style <- match.arg(style)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  s.lim <- c(max=2*s, min=0)
  ds <- unname(s.lim['max'] - s.lim['min']) / m
  s.seq <- s.lim['min'] + 0:m*ds            # vector, m+1 elements

  f <- matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call = pmax(s.seq - k, 0), put = pmax(k - s.seq, 0))

  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    for (j in (m-1):1) {
      a <- (1 + r*dt)^-1 * (-1/2*r*j*dt + 1/2*sd^2*j^2*dt)
      b <- (1 + r*dt)^-1 * (1 - sd^2*j^2*dt)
      c <- (1 + r*dt)^-1 * (1/2*r*j*dt + 1/2*sd^2*j^2*dt)
      f[i,g2m(j)] <- t(c(a,b,c)) %*% f[i+1,g2m((j-1):(j+1))]
    }
    if (style == 'american')
      f[i,g2m((m-1):1)] <-
        pmax(f[i,g2m((m-1):1)], switch(type, call = s.seq[g2m((m-1):1)] - k,
                                       put = k - s.seq[g2m((m-1):1)]))
    f[i,g2m(m)] <- switch(type,         # call: ∂C/∂S ≈ 1; put: ∂C/∂S ≈ 0
                          call = f[i,g2m(m-1)]+s.seq[g2m(m)]-s.seq[g2m(m-1)],
                          put = f[i,g2m(m-1)])
    f[i,g2m(0)] <- switch(type,         # call: ∂C/∂S ≈ 0; put: ∂C/∂S ≈ 1
                          call = f[i,g2m(1)],
                          put = f[i,g2m(1)] - (s.seq[g2m(m)] - s.seq[g2m(m-1)]))
    if (style == 'american') {
      f[i,g2m(m)] <- pmax(f[i,g2m(m)], switch(type, call = s.seq[g2m(m)] - k,
                                              put = k - s.seq[g2m(m)]))
      f[i,g2m(0)] <- pmax(f[i,g2m(0)], switch(type, call = s.seq[g2m(0)] - k,
                                              put = k - s.seq[g2m(0)]))
    }
  }

  f[g2m(0), g2m(m/2)]
}

fde.log <- function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(n)),
                    type = c("call", "put"), optim = TRUE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim <- c(max=log(s) + (r-1/2*sd^2)*t + sd*3*sqrt(t),
             min=log(s) + (r-1/2*sd^2)*t - sd*3*sqrt(t))
  dz <- unname(z.lim['max'] - z.lim['min']) / m
  z <- z.lim['min'] + 0:m*dz            # vector, m+1 elements

  f <- matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z) - k, 0), put=pmax(k - exp(z), 0))

  a <- (1 + r*dt)^-1 * (-dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  b <- (1 + r*dt)^-1 * (1 - dt/dz^2*sd^2)
  c <- (1 + r*dt)^-1 * (dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    if (optim)
      f[i,g2m(1:(m-1))] <- matrix(f[i+1,rep(g2m(0:2),m-1)+rep(0:(m-2),each=3)],
                                  ncol=3, byrow=TRUE) %*% c(a,b,c)
    else
      for (j in g2m((m-1):1))
        f[i,j] <- t(c(a,b,c)) %*% f[i+1,(j-1):(j+1)]
    f[i,g2m(m)] <- switch(type,         # call: ∂C/∂S ≈ 1; put: ∂C/∂S ≈ 0
                          call = f[i,g2m(m-1)]+exp(z[g2m(m)])-exp(z[g2m(m-1)]),
                          put = f[i,g2m(m-1)])
    f[i,g2m(0)] <- switch(type,         # call: ∂C/∂S ≈ 0; put: ∂C/∂S ≈ 1
                          call = f[i,g2m(1)],
                          put = f[i,g2m(1)] - (exp(z[g2m(m)])-exp(z[g2m(m-1)])))
  }

  f[g2m(0), g2m(m/2)]
}

######################################################################
## Black-Scholes-Mertons Option Pricing                             ##
######################################################################

bsm.d1 <- function (s, k, r, t, sd) {
  (log(s/k) + (r + 1/2*sd^2)*t) / (sd*sqrt(t))
}

bsm.d2 <- function(s, k, r, t, sd) {
  bsm.d1(s, k, r, t, sd) - sd*sqrt(t)
}

bsm.option <- function (s, k, r, t, sd, type = c("call", "put")) {
  type <- match.arg(type)
  d1 <- bsm.d1(s, k, r, t, sd)
  d2 <- bsm.d2(s, k, r, t, sd)
  r <- switch(type,
              call = s*pnorm(d1) - k*exp(-r*t)*pnorm(d2),
              put = k*exp(-r*t)*pnorm(d2,lower.tail=F)-s*pnorm(d1,lower.tail=F))
  return(r)
}
