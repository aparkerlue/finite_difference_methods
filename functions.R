## All arguments are scalar.
fde <- function(s, k, r, t, sd, n = ceiling(1e3*t), m = ceiling(2*sqrt(n)),
                type = c("call", "put")) {
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(n) || n <= 0) stop("n = ", n, " is not a positive integer!")
  if (!is.wholenumber(m) || m <= 0) stop("m = ", m, " is not a positive integer!")
  type <- match.arg(type)

  dt <- t / n
  m <- m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim <- c(max=log(s) + (r-1/2*sd^2)*t + sd*sqrt(t)*3,
             min=log(s) + (r-1/2*sd^2)*t - sd*sqrt(t)*3)
  dz <- unname(z.lim['max'] - z.lim['min']) / m
  z.t <- z.lim['min'] + 0:m*dz          # vector, m+1 elements

  f <- matrix(rep(NA,(n+1)*(m+1)), nrow=n+1)
  g2m <- function(i)  i + 1             # grid index -> matrix index
  if (type == 'call') {
    f[g2m(n),] <- pmax(exp(z.t) - k, 0)
    f[,g2m(m)] <- pmax(exp(z.lim['max']) - k, 0)
    ##f[,g2m(0)] <- pmax(exp(z.lim['min']) - k, 0)
  }
  else if (type == 'put') {
    f[g2m(n),] <- pmax(k - exp(z.t), 0)
    f[,g2m(m)] <- pmax(k - exp(z.lim['max']), 0)
    f[,g2m(0)] <- pmax(k - exp(z.lim['min']), 0)
  }
  else
    stop("Unknown option type '", type, "'")

  a <- (1 + r*dt)^-1 * (-dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  b <- (1 + r*dt)^-1 * (1 - dt/dz^2*sd^2)
  c <- (1 + r*dt)^-1 * (dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  for (i in g2m((n-1):0)) {
    for (j in g2m((m-1):1))
      f[i,j] <- t(c(a,b,c)) %*% f[i+1,(j-1):(j+1)]
    f[,g2m(0)] = f[,g2m(1)]             # ∂C/∂S ≈ 0
  }

  return(f[g2m(0),g2m(m/2)])
}
