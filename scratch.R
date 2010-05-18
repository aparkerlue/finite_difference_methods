######################################################################
## Tests                                                            ##
######################################################################

source('functions.R')

optrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='p', style='a')

## Binomial model.
##
## For n = 1000, the price of an American put is $4.28. 
with(optrc, binom.american.put(s, k, r, t, sd, n=1e3))

## Implicit finite difference method.
## 
## American put yields $4.07; European put yields $3.91.
with(optrc, {
  grid.a <- rotate(fdi(s, k, r, t, sd, n=10, m=20, type, style, grid=T), 2)
  price.a <- fdi(s, k, r, t, sd, n=10, m=20, type, style)
  price.e <- fdi(s, k, r, t, sd, n=10, m=20, type, style='e')
  print(grid.a)
  cat('American put:', round(price.a, 2), '\n')
  cat('European put:', round(price.e, 2), '\n')
})

## Implicit finite difference method, using log transform.
##
## American put yields $4.11, but it should be around $4.28.
with(optrc, {
  grid.a <- rotate(fdi.log(s, k, r, t, sd, n=10, m=20, type, style, grid=T), 2)
  price.a <- fdi.log(s, k, r, t, sd, n=10, m=20, type, style)
  print(grid.a)
  cat('American put:', round(price.a, 2), '\n')
})
## Using programmed default values for n and m, computation converges
## to correct value of $4.28.
with(optrc, fdi.log(s, k, r, t, sd, type=type, style=style))

## Explicit finite difference method.
##
## American put yields $4.26.
with(optrc, {
  grid.a <- rotate(fde(s, k, r, t, sd, n=10, m=20, type, style, grid=T), 2)
  price.a <- fde(s, k, r, t, sd, n=10, m=20, type, style)
  print(grid.a)
  cat('American put:', round(price.a, 2), '\n')
})

## Explicit finite difference method, using log transform.
##
## American put yields $3.62, but it should be around $4.28.
with(optrc, {
  grid.a <- rotate(fde.log(s, k, r, t, sd, n=10, m=20, type, style, grid=T), 2)
  price.a <- fde.log(s, k, r, t, sd, n=10, m=20, type, style)
  print(grid.a)
  cat('American put:', round(price.a, 2), '\n')
})
## Using programmed default values for n and m, computation converges
## to correct value of $4.28.
with(optrc, fde.log(s, k, r, t, sd, type=type, style=style))

## Crank-Nicolson scheme, using log transform.
##
## American put yields $4.21, but it should be around $4.28.
with(optrc, {
  grid.a <- rotate(fdcn.log(s, k, r, t, sd, n=10, m=20, type, style, grid=T), 2)
  price.a <- fdcn.log(s, k, r, t, sd, n=10, m=20, type, style)
  print(grid.a)
  cat('American put:', round(price.a, 2), '\n')
})
## Using programmed default values for n and m, computation converges
## to correct value of $4.28.
with(optrc, fdcn.log(s, k, r, t, sd, type=type, style=style))

######################################################################
## Scratchwork                                                      ##
######################################################################

testfn.stop <- function(type = c("call", "put")) {
  type <- match.arg(type)
  if (type == "call" || type == "put")
    return(type)
  else
    stop("Unknown option type '", type, "'.")
}

testfn.innerfn <- function() {
  fn <- function(i) i+1
  s <- 3:5
  return(fn(s))
}
