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
