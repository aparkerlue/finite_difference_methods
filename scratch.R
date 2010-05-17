######################################################################
## Tests                                                            ##
######################################################################

source('functions.R')

optrc <- list(s=50, k=50, r=0.10, t=5/12, sd=0.40, type='p', style='a')

## Implicit finite difference method.
## 
## American put yields $4.07; European put yields $3.91.
with(optrc, {
  grid.a <- rotate(fdi(s, k, r, t, sd, n=10, m=20, type, style, grid=T), 2)
  price.a <- fdi(s, k, r, t, sd, n=10, m=20, type, style)
  price.e <- fdi(s, k, r, t, sd, n=10, m=20, type, 'e')
  print(grid.a)
  cat('American put:', round(price.a, 2), '\n')
  cat('European put:', round(price.e, 2), '\n')
})

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
