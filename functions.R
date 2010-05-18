source('util.R')

######################################################################
## Implicit Finite Difference Method                                ##
######################################################################

fdi <- function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                type = c("call", "put"), style = c("european", "american"),
                grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type); style <- match.arg(style)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## FIXME: s.max depends on k, t, sd
  s.lim <- c(max=2*s, min=0)
  ds <- unname(s.lim['max'] - s.lim['min']) / m
  s.seq <- s.lim['min'] + 0:m*ds        # vector, m+1 elements

  f <- matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call = pmax(s.seq - k, 0), put = pmax(k - s.seq, 0))
  f[,g2m(m)] = switch(type, call = s.seq[g2m(m)] - k,  put = 0)
  f[,g2m(0)] = switch(type, call = 0,                  put = k)

  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    A <- matrix(0, nrow=m-1, ncol=m-1)
    B <- f[i+1,g2m(1:(m-1))]
    for (j in 1:(m-1)) {
      a <- 1/2*r*j*dt - 1/2*sd^2*j^2*dt
      b <- 1 + sd^2*j^2*dt + r*dt
      c <- -1/2*r*j*dt - 1/2*sd^2*j^2*dt
      if (j <= 1) {                     # j == 1
        A[j,1:2] <- c(b, c)
        B[1] <- B[1] - a*f[i,g2m(0)]
      }
      else if (j < m-1)
        A[j,(j-1):(j+1)] <- c(a, b, c)
      else {                            # j == m-1
        A[j,(m-2):(m-1)] <- c(a, b)
        B[m-1] <- B[m-1] - c*f[i,g2m(m)]
      }
    }
    f[i,g2m(1:(m-1))] <- solve(A, B)
    
    if (type == 'put' && style == 'american')
      f[i,] <- pmax(f[i,], k - s.seq)
  }

  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

fdi.log <- function(s, k, r, t, sd,
                    n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                    type = c("call", "put"), style = c("european", "american"),
                    grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type); style <- match.arg(style)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim <- log(s) + 3*sd*sqrt(t)*c(min=-1, max=1)
  dz <- unname(diff(z.lim)) / m
  z.seq <- z.lim['min'] + 0:m*dz        # vector, m+1 elements

  f <- matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z.seq)-k,0), put=pmax(k-exp(z.seq),0))
  f[,g2m(m)] = switch(type, call=exp(z.seq[g2m(m)])-k, put=0)
  f[,g2m(0)] = switch(type, call=0,                    put=k-exp(z.seq[g2m(0)]))

  a <- dt/(2*dz)*(r - 1/2*sd^2) - dt/(2*dz^2)*sd^2
  b <- 1 + dt/dz^2*sd^2 + r*dt
  c <- -dt/(2*dz)*(r - 1/2*sd^2) - dt/(2*dz^2)*sd^2
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    A <- matrix(0, nrow=m-1, ncol=m-1)
    B <- f[i+1,g2m(1:(m-1))]
    for (j in 1:(m-1)) {
      if (j <= 1) {                     # j == 1
        A[j,1:2] <- c(b, c)
        B[1] <- B[1] - a*f[i,g2m(0)]
      }
      else if (j < m-1)
        A[j,(j-1):(j+1)] <- c(a, b, c)
      else {                            # j == m-1
        A[j,(m-2):(m-1)] <- c(a, b)
        B[m-1] <- B[m-1] - c*f[i,g2m(m)]
      }
    }
    f[i,g2m(1:(m-1))] <- solve(A, B)
    
    if (type == 'put' && style == 'american')
      f[i,] <- pmax(f[i,], k - exp(z.seq))
  }

  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

######################################################################
## Explicit Finite Difference Method                                ##
######################################################################

fde <- function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                type = c("call", "put"), style = c("european", "american"),
                grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type); style <- match.arg(style)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## FIXME: s.max depends on k, t, sd
  s.lim <- c(max=2*s, min=0)
  ds <- unname(s.lim['max'] - s.lim['min']) / m
  s.seq <- s.lim['min'] + 0:m*ds        # vector, m+1 elements

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

    if (type == 'call') {               # m: ∂C/∂S ≈ 1
      f[i,g2m(m)] <- f[i,g2m(m-1)] + s.seq[g2m(m)] - s.seq[g2m(m-1)]
      f[i,g2m(0)] <- f[i,g2m(1)]        # 0: ∂C/∂S ≈ 0
    }
    else if (type == 'put') {           # m: ∂C/∂S ≈ 0
      f[i,g2m(m)] <- f[i,g2m(m-1)]      # 0: ∂C/∂S ≈ 1
      f[i,g2m(0)] <- f[i,g2m(1)] - (s.seq[g2m(1)] - s.seq[g2m(0)])
      if (style == 'american')
        f[i,] <- pmax(f[i,], k - s.seq)
    }
  }

  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

fde.log <- function(s, k, r, t, sd,
                    n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                    type = c("call", "put"), style = c("european", "american"),
                    grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type); style <- match.arg(style)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim <- log(s) + 3*sd*sqrt(t)*c(min=-1, max=1)
  dz <- unname(diff(z.lim)) / m
  z.seq <- z.lim['min'] + 0:m*dz        # vector, m+1 elements

  f <- matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z.seq)-k,0), put=pmax(k-exp(z.seq),0))

  a <- (1 + r*dt)^-1 * (-dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  b <- (1 + r*dt)^-1 * (1 - dt/dz^2*sd^2)
  c <- (1 + r*dt)^-1 * (dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    j.seq <- rep(g2m(0:2), times=m-1) + rep(0:(m-2), each=3)
    f[i,g2m(1:(m-1))] <- matrix(f[i+1,j.seq], ncol=3, byrow=TRUE) %*% c(a,b,c)

    if (type == 'call') {               # m: ∂C/∂S ≈ 1
      f[i,g2m(m)] <- f[i,g2m(m-1)] + exp(z.seq[g2m(m)]) - exp(z.seq[g2m(m-1)])
      f[i,g2m(0)] <- f[i,g2m(1)]        # 0: ∂C/∂S ≈ 0
    }
    else if (type == 'put') {           # m: ∂C/∂S ≈ 0
      f[i,g2m(m)] <- f[i,g2m(m-1)]      # 0: ∂C/∂S ≈ 1
      f[i,g2m(0)] <- f[i,g2m(1)] - (exp(z.seq[g2m(m)]) - exp(z.seq[g2m(m-1)]))
      if (style == 'american')
        f[i,] <- pmax(f[i,], k - exp(z.seq))
    }
  }

  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

######################################################################
## Crank-Nicolson Scheme                                            ##
######################################################################

fdcn.log <- function(s, k, r, t, sd,
                    n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                    type = c("call", "put"), style = c("european", "american"),
                    grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  if (!is.wholenumber(n) || n <= 0) stop("n = ",n," is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ",m," is not a positive integer!")
  type <- match.arg(type); style <- match.arg(style)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim <- log(s) + 3*sd*sqrt(t)*c(min=-1, max=1)
  dz <- unname(diff(z.lim)) / m
  z.seq <- z.lim['min'] + 0:m*dz        # vector, m+1 elements

  f <- matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z.seq)-k,0), put=pmax(k-exp(z.seq),0))
  f[,g2m(m)] = switch(type, call=exp(z.seq[g2m(m)])-k, put=0)
  f[,g2m(0)] = switch(type, call=0,                    put=k-exp(z.seq[g2m(0)]))

  a <- dt/dz^2
  e <- diag(sd^2*a/2, m-1) - rbind(cbind(0, diag(sd^2*a/4, m-2)), 0) -
    rbind(0, cbind(diag(sd^2*a/4, m-2), 0))
  c <- diag(1, m-1) + e
  d <- diag(1, m-1) - e
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    b <- c(sum(f[i:(i+1),g2m(0)])/2, rep(0, m-3), sum(f[i:(i+1),g2m(m)])/2)
    rhs <- d %*% f[i+1,g2m(1:(m-1))] + sd^2*a/2*b
    f[i,g2m(1:(m-1))] <- solve(c, rhs)
    
    if (type == 'put' && style == 'american')
      f[i,] <- pmax(f[i,], k - exp(z.seq))
  }

  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

######################################################################
## Black-Scholes-Mertons Option Pricing                             ##
######################################################################

bsm.d1 <- function(s, k, r, t, sd) {
  (log(s/k) + (r + 1/2*sd^2)*t) / (sd*sqrt(t))
}

bsm.d2 <- function(s, k, r, t, sd) {
  bsm.d1(s, k, r, t, sd) - sd*sqrt(t)
}

bsm.option <- function(s, k, r, t, sd, type = c("call", "put")) {
  type <- match.arg(type)
  d1 <- bsm.d1(s, k, r, t, sd)
  d2 <- bsm.d2(s, k, r, t, sd)
  r <- switch(type,
              call = s*pnorm(d1) - k*exp(-r*t)*pnorm(d2),
              put = k*exp(-r*t)*pnorm(d2,lower.tail=F)-s*pnorm(d1,lower.tail=F))
  return(r)
}

######################################################################
## Binomial Option Pricing                                          ##
######################################################################

binom.american.put <- function(s, k, r, t, sd, n,
                               u = exp(sd*sqrt(t/n)), d = 1/u,
                               p = (exp(r*t/n) - d) / (u - d)) {
  n.u <- n:0
  s.t <- s * u^n.u * d^(n - n.u)
  v <- pmax(k - s.t, 0)
  for (i in (n-1):0) {
    w <- matrix(rep(v, c(1, rep(2, times=length(v)-2), 1)), ncol=2, byrow=TRUE)
    cv <- exp(-r*t/n) * (p*w[,1] + (1-p)*w[,2])

    n.u <- i:0
    s.i <- s * u^n.u * d^(i - n.u)
    ev <- pmax(k - s.i, 0)

    v <- pmax(cv, ev)
  }
  
  return(v)
}
