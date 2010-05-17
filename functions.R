######################################################################
## Explicit Finite Difference Method                                ##
######################################################################

## All arguments are scalar.
fde <- function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(n)),
                type = c("call", "put"), optim = TRUE) {
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(n) || n <= 0) stop("n = ", n, " is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ", m, " is not a positive integer!")
  type <- match.arg(type)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim <- c(max=log(s) + (r-1/2*sd^2)*t + sd*3*sqrt(t),
             min=log(s) + (r-1/2*sd^2)*t - sd*3*sqrt(t))
  dz <- unname(z.lim['max'] - z.lim['min']) / m
  z <- z.lim['min'] + 0:m*dz          # vector, m+1 elements

  f <- matrix(rep(NA,(n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index to matrix index
  if (type == 'call') {
    f[g2m(n),] <- pmax(exp(z) - k, 0)
    ##f[,g2m(m)] <- pmax(exp(z.lim['max']) - k, 0)
    ##f[,g2m(0)] <- pmax(exp(z.lim['min']) - k, 0)
  }
  else if (type == 'put') {
    f[g2m(n),] <- pmax(k - exp(z), 0)
    f[,g2m(m)] <- pmax(k - exp(z.lim['max']), 0)
    f[,g2m(0)] <- pmax(k - exp(z.lim['min']), 0)
  }
  else
    stop("Unknown option type '", type, "'")

  a <- (1 + r*dt)^-1 * (-dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  b <- (1 + r*dt)^-1 * (1 - dt/dz^2*sd^2)
  c <- (1 + r*dt)^-1 * (dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  for (i in g2m((n-1):0)) {
    if (optim) {
      f.next <- matrix(f[i+1,rep(g2m(0:2),m-1)+rep(0:(m-2),each=3)],ncol=3,byrow=T)
      f[i,g2m(1:(m-1))] <- f.next %*% c(a,b,c)
    }
    else
      for (j in g2m((m-1):1))
        f[i,j] <- t(c(a,b,c)) %*% f[i+1,(j-1):(j+1)]
    f[i,g2m(m)] = f[i,g2m(m-1)] + exp(z[g2m(m)]) - exp(z[g2m(m-1)])
    f[i,g2m(0)] = f[i,g2m(1)]           # ∂C/∂S ≈ 0
  }

  f[g2m(0),g2m(m/2)]
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
              put = k*exp(-r*t)*pnorm(d2,lower.tail=F) - s*pnorm(d1,lower.tail=F))
  return(r)
}
